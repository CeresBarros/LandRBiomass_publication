## ------------------------------------------------------
## LandR Biomass PUBLICATION SIMULATIONS
##
## Ceres: June 2018
## ------------------------------------------------------

## /!\ PLEASE MAKE SURE YOU ARE USING R v4.0.5 /!\

## set CRAN repos; use binary linux packages if on Ubuntu
local({
  options("repos" = c(CRAN = "https://cran.rstudio.com"))

  if (Sys.info()["sysname"] == "Linux" && grepl("Ubuntu", utils::osVersion)) {
    .os.version <- strsplit(system("lsb_release -c", intern = TRUE), ":\t")[[1]][[2]]
    .user.agent <- paste0(
      "R/", getRversion(), " R (",
      paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"]),
      ")"
    )

    .os.version <- system("lsb_release -cs", intern = TRUE)
    options(
      repos = c(
        CRAN = if (!grepl("R Under development", R.version.string) && getRversion() >= "4.1") {
          paste0("https://packagemanager.rstudio.com/all/__linux__/", .os.version, "/latest")
        } else {
          "https://cloud.r-project.org"
        }
      )
    )

    options(HTTPUserAgent = .user.agent)
  }
})

## package installation location
pkgDir <- file.path("packages", version$platform,
                    paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1]))
dir.create(pkgDir, recursive = TRUE)
.libPaths(pkgDir, include.site = FALSE)

if (!"remotes" %in% installed.packages()) {
  install.packages("remotes")
}
remotes::install_github("PredictiveEcology/Require@archivedPkg", upgrade = FALSE)
remotes::install_github("PredictiveEcology/SpaDES.install@development", upgrade = FALSE)
options("Require.RPackageCache" = "~/.cache/RPackages/")

# options("Require.standAlone" = TRUE,
#         "Require.unloadNamespaces" = FALSE)
modulePath <- "R/SpaDES/m"
SpaDES.install::getModule(modulePath = modulePath,
                          c("CeresBarros/Biomass_speciesData@master",
                            "CeresBarros/Biomass_borealDataPrep@master",
                            "CeresBarros/Biomass_core@master",
                            "CeresBarros/Biomass_validationKNN@master",
                            "CeresBarros/Biomass_speciesParameters@temp"))

outs <- SpaDES.install::packagesInModules(modulePath = modulePath)
Require::Require(c(unname(unlist(outs)), "PredictiveEcology/SpaDES.experiment@development"),
                 require = FALSE, standAlone = TRUE)

## load packages and make sure minimum versions are installed.
Require::Require(c("raster", "dplyr", "data.table", "future",
                   "SpaDES.core", "SpaDES.experiment",
                   "LandR",
                   "reproducible"), upgrade = FALSE, install = FALSE)

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------

## Set up modelling parameters  ---------------------------
options("reproducible.useNewDigestAlgorithm" = 2)
options("spades.moduleCodeChecks" = FALSE)
options("reproducible.useCache" = TRUE)
options("reproducible.inputPaths" = file.path("R/SpaDES/inputs"))  ## store everything in inputs/ so that there are no duplicated files across modules
options("reproducible.destinationPath" = file.path("R/SpaDES/inputs"))
options("reproducible.useGDAL" = FALSE)
options("spades.useRequire" = FALSE)
options("reproducible.useTerra" = FALSE)

eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## TEST STUDY AREAS
## use studyAreaS or L to parameterise AND run simualations in the SAME area
## one study area from set A
# runName <- "studyAreaS"
# runName <- "studyAreaL"

## SIMULATIONS IN BARROS ET AL 2022 MEE
## use one of the following to parameterise the model in a larger study area than the simulation area
## baseCase uses set A of study areas, parameterises using Biomass_borealDataPrep and Biomass_speciesParameters
## studyAreaChange uses set B of study areas, parameterises using Biomass_borealDataPrep and Biomass_speciesParameters
## altParameterisation uses set A of study areas, parameterises using Biomass_borealDataPrep only
runName <- "baseCase"
# runName <- "studyAreaChange"
# runName <- "altParameters"

## paths
simDirName <- "mar2022Runs"
simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName)
                 , modulePath = modulePath
                 , inputPath = file.path("R/SpaDES/inputs")
                 , outputPath = file.path("R/SpaDES/outputs", simDirName, runName))

figDir <- "R/SpaDES/outputs/GeneralFigs"
dir.create(figDir, recursive = TRUE)

## Get necessary objects like the study area.
source("R/SpaDES/1_simObjects.R")

## Run Biomass_speciesData to get species layers
## running this separately from other modules makes switching
## between using a large and a smaller study area easier when the smaller SA is within the large one,
## as it keeps the data in separate folders that can be used across simulations/scenarios
source("R/SpaDES/2_speciesLayers.R")

## check species layers:
# plot(simOutSpeciesLayers$speciesLayers)
## Populus grandidentata shouldn't be in Saskatchewan (and has only v. few pixels in the layer) and will be excluded
toRm <- which(names(simOutSpeciesLayers$speciesLayers) %in% c("Popu_Gra"))
simOutSpeciesLayers$speciesLayers <- dropLayer(simOutSpeciesLayers$speciesLayers, i = toRm)
rm(toRm)

## subset sppEquivalencies
sppEquivalencies_CA <- sppEquivalencies_CA[Boreal %in% names(simOutSpeciesLayers$speciesLayers)]

## Get land-cover raster now that we have a rasterToMatchLarge
rstLCC2005 <- Cache(prepInputs,
                    targetFile = "NA_LandCover_2005_V3_25haMMU.tif",
                    url = "http://www.cec.org/wp-content/uploads/wpallimport/files/Atlas/Files/Land_Cover_2005/Land_Cover_2005v3_TIFF.zip",
                    destinationPath = simPaths$inputPath,
                    studyArea = simOutSpeciesLayers$studyAreaLarge,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
                    rasterToMatch = simOutSpeciesLayers$rasterToMatchLarge,
                    filename2 = .suffix("rstLCC.tif", paste0("_", SAname)),
                    overwrite = TRUE,
                    cacheRepo = simPaths$cachePath,
                    userTags = c("rstLCC", SAname),
                    omitArgs = c("userTags"))

## simulation params
simTimes <- list(start = 2001, end = 2031)
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
successionTimestep <- 10L  # for dispersal and age reclass.

## list of values of species shade tolerance traits
speciesParams <- list(
  "shadetolerance" = list(
    Betu_Pap = 1,
    Lari_Lar = 1,
    Pice_Gla = 2,
    Pice_Mar = 3,
    Pinu_Ban = 1.5,
    Popu_Spp = 1
  )
)

if (runName %in% c("baseCase", "studyAreaChange", "studyAreaS", "studyAreaL")) {
  simModules <- list("Biomass_borealDataPrep"
                     , "Biomass_speciesParameters"
                     , "Biomass_core"
  )
} else if (runName == "altParameters") {
  ## no need to change parameters list, superfluous params. will simply not be used
  ## same for objects
  simModules <- list("Biomass_borealDataPrep"
                     , "Biomass_core"
  )
} else {
  stop("runName must be one of 'baseCase', 'studyAreaChange', 'altParameters', 'studyAreaS' or 'studyAreaL'")
}

simParams <- list(
  .globals = list("dataYear" = 2001L    ## will not be used as the layers have been pre-preped, but just in case...
                  , "initialB" = NA     ## use LANDIS approach to estimate initial cohort B
                  , "sppEquivCol" = sppEquivCol
                  , "vegLeadingProportion" = vegLeadingProportion
                  , ".sslVerify" = 0L
                  , ".useCache" = eventCaching),
  Biomass_borealDataPrep = list(
    "fitDeciduousCoverDiscount" = TRUE
    , "subsetDataAgeModel" = FALSE
    , "subsetDataBiomassModel" = FALSE
    , "exportModels" = "all"
    , "fixModelBiomass" = TRUE
    , "speciesTableAreas" = c("BSW", "BP")
    , "speciesUpdateFunction" = list(
      quote(LandR::speciesTableUpdate(sim$species, sim$speciesTable, sim$sppEquiv, P(sim)$sppEquivCol)),
      quote(LandR::updateSpeciesTable(speciesTable = sim$species, params = sim$speciesParams)))
    # next two are used when assigning pixelGroup membership; what resolution for
    #   age and biomass
    , "pixelGroupAgeClass" = successionTimestep
    , "pixelGroupBiomassClass" = 100
    , "useCloudCacheForStats" = FALSE
    , "cloudFolderID" = NA
    , ".plots" = c("object", "raw")
  )
  , Biomass_speciesParameters = list(
    "quantileAgeSubset" = list(Betu_Pap = 95, Lari_Lar = 95, Pice_Gla = 95, Pice_Mar = 95, Pinu_Ban = 99, Popu_Spp = 99)
  )
  , Biomass_core = list(
    "calcSummaryBGM" = c("start")
    , "initialBiomassSource" = "cohortData"
    , "plotOverstory" = TRUE
    , "seedingAlgorithm" = "wardDispersal"
    , "successionTimestep" = successionTimestep
    , ".plotInitialTime" = simTimes$start
    , ".plotInterval" = 1L
    , ".plots" = c("object", "raw")
    , ".plotMaps" = FALSE
    , ".saveInitialTime" = NA
    , ".useCache" = eventCaching[1] # experiment doesn't like when init is cached
    , ".useParallel" = useParallel
  )
)

simObjects <- list(
  "rstLCC" = rstLCC2005
  , "sppEquiv" = sppEquivalencies_CA
  , "sppColorVect" = sppColorVect
  , "speciesLayers" = simOutSpeciesLayers$speciesLayers
  , "speciesParams" = speciesParams
  , "treed" = simOutSpeciesLayers$treed
  , "numTreed" = simOutSpeciesLayers$numTreed
  , "nonZeroCover" = simOutSpeciesLayers$nonZeroCover
)

if (grepl("studyArea(S|L)$", runName)) {
  simObjects$studyArea <- get(runName)
} else {
  simObjects$studyArea <- studyAreaS
  simObjects$studyAreaLarge <- studyAreaL
}

## objects will be saved at the start of the simulation (so they reflect the previous year)
simOutputs <- data.frame(expand.grid(objectName = "cohortData",
                                     saveTime = unique(seq(simTimes$start, simTimes$end, by = 1)),
                                     eventPriority = 1,
                                     stringsAsFactors = FALSE))
simOutputs <- rbind(simOutputs, data.frame(objectName = "pixelGroupMap",
                                           saveTime = unique(seq(simTimes$start, simTimes$end, by = 1)),
                                           eventPriority = 1))
simOutputs <- rbind(simOutputs, data.frame(objectName = "biomassMap",
                                           saveTime = simTimes$start,
                                           eventPriority = 1))

## in the first year, eventPriorities need to be set to AFTER the init event (which has priority 1)
simOutputs$eventPriority[simOutputs$saveTime == simTimes$start] <- 1.5

## make a initialisation simList and run init events too
LandRBiomass_simInit <- Cache(simInitAndSpades
                              , times = simTimes
                              , params = simParams
                              , modules = simModules
                              , loadOrder = unlist(simModules)
                              , objects = simObjects
                              , paths = simPaths
                              , outputs = simOutputs
                              , events = "init"
                              # , .plots = "screen"
                              , .studyAreaName = SAname
                              , userTags = c("simInitAndInits", SAname)
                              , cacheRepo = simPaths$cachePath
                              , omitArgs = c("userTags", ".plotInitialTime"))

saveSimList(LandRBiomass_simInit, file.path(simPaths$outputPath, paste0("simInit", runName, ".qs")))   ## only save in first runs

amc::.gc()  ## clean ws
if (Sys.info()["sysname"] == "Windows") {
  plan("multisession", workers = 2)   ## each worker consuming roughly 6Gb
} else {
  plan("multicore", workers = 2)
}
LandRBiomass_sim <- experiment2(
  sim1 = LandRBiomass_simInit,
  clearSimEnv = TRUE,
  replicates = 10)
future:::ClusterRegistry("stop")

## save simLists object.
qs::qsave(LandRBiomass_sim, file.path(simPaths$outputPath, paste0("simList_LandRBiomass_sim_", runName, ".qs")))

## VALIDATION
## get the  land-cover change map (need to have a rasterToMatch first, so get it from the simInitList)
## /!\ it is assumed that the filename of the raster in the simList corresponds to the raster found in disk.
## this may not be the case if the simulations were run in another machine and saved rasters were not imported.

validationPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName)
                        , modulePath = file.path("R/SpaDES/m")
                        , inputPath = file.path("R/SpaDES/inputs")
                        , outputPath = file.path("R/SpaDES/validation", simDirName, runName))

source("R/SpaDES/3_simObjects4Valid.R")

## PARAMETERS FOR VALIDATION MODULE
## in this case all reps have the same parameters, so we can use the first rep to get the values
validationTimes <- list(start = 1, end = 1)
validationParams <- list(
  Biomass_validationKNN = list(
    "minCoverThreshold" = params(LandRBiomass_simInit)$Biomass_borealDataPrep$minCoverThreshold
    , "pixelGroupBiomassClass" = params(LandRBiomass_simInit)$Biomass_borealDataPrep$pixelGroupBiomassClass
    , "deciduousCoverDiscount" = params(LandRBiomass_simInit)$Biomass_borealDataPrep$deciduousCoverDiscount
    , "sppEquivCol" = params(LandRBiomass_simInit)$Biomass_borealDataPrep$sppEquivCol
    , "validationReps" = as.integer(1:10)  ## or length of simLists
    , "validationYears" = as.integer(c(2001, 2011))
    , ".plots" = c("object", "raw", "png")
    , ".useCache" = eventCaching
  )
)

validationObjects <- list(
  "biomassMap" = biomassMap
  , "rasterToMatch" = rasterToMatch
  , "rawBiomassMapStart" = rawBiomassMap
  , "simulationOutputs" = simulationOutputs
  , "speciesLayersStart" = speciesLayers
  , "sppColorVect" = LandRBiomass_simInit$sppColorVect
  , "sppEquiv" = LandRBiomass_simInit$sppEquiv
  , "standAgeMapStart" = standAgeMap
  , "studyArea" = LandRBiomass_simInit$studyArea
)

## the following objects are only saved once at the end of year 0/beggining of year 1 (they don't change)
validationOutputs <- data.frame(expand.grid(objectName = c("rawBiomassMapStart"),
                                            saveTime = c(validationTimes$start),
                                            eventPriority = 1),
                                stringsAsFactors = FALSE)
validationOutputs <- rbind(validationOutputs, data.frame(objectName = "rawBiomassMapEnd",
                                                         saveTime = c(validationTimes$start),
                                                         eventPriority = 1))
validationOutputs <- rbind(validationOutputs, data.frame(objectName = "standAgeMapStart",
                                                         saveTime = c(validationTimes$start),
                                                         eventPriority = 1))
validationOutputs <- rbind(validationOutputs, data.frame(objectName = "standAgeMapEnd",
                                                         saveTime = c(validationTimes$start),
                                                         eventPriority = 1))
validationOutputs <- rbind(validationOutputs, data.frame(objectName = "speciesLayersStart",
                                                         saveTime = c(validationTimes$start),
                                                         eventPriority = 1))
validationOutputs <- rbind(validationOutputs, data.frame(objectName = "speciesLayersEnd",
                                                         saveTime = c(validationTimes$start),
                                                         eventPriority = 1))
SAname <- studyAreaName(validationObjects$studyArea)
LandRBiomass_validation <- simInitAndSpades(times = validationTimes
                                            , params = validationParams
                                            , modules = "Biomass_validationKNN"
                                            , objects = validationObjects
                                            , outputs = validationOutputs
                                            , paths = validationPaths
                                            , .studyAreaName = SAname)

saveSimList(LandRBiomass_validation, file.path(validationPaths$outputPath, paste0("simValid", runName, ".qs")))   ## only save in first runs

## -----------------------------------------------
## POST-HOC ANALYSIS - figures
## -----------------------------------------------
source("R/SpaDES/pubFigures.R")

q("no")

