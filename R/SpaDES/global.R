## ------------------------------------------------------
## LandR Biomass PUBLICATION SIMULATIONS
##
## Ceres Barros: last updated September 2022
## ------------------------------------------------------

## /!\    PLEASE MAKE SURE YOU ARE USING R v4.2   /!\
## /!\ PLEASE MAKE SURE YOU HAVE A GOOGLE ACCOUNT /!\

## set CRAN repo
options(repos = c(CRAN = "https://cloud.r-project.org"))

## package installation location
pkgPath <- file.path("packages", version$platform,
                     paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1]))
dir.create(pkgPath, recursive = TRUE)
.libPaths(pkgPath, include.site = FALSE)

if (!"remotes" %in% installed.packages()) {
  install.packages("remotes")
}

if (!"Require" %in% installed.packages(lib.loc = pkgPath) ||
    packageVersion("Require", lib.loc = pkgPath) < "0.1.6.9014") {
  remotes::install_github("PredictiveEcology/Require@5c44205bf407f613f53546be652a438ef1248147",
                          upgrade = FALSE, force = TRUE)
}
## use binary linux packages if on Ubuntu
Require::setLinuxBinaryRepo()

Require::Require("PredictiveEcology/SpaDES.project@6d7de6ee12fc967c7c60de44f1aa3b04e6eeb5db",
                 require = FALSE, upgrade = FALSE, standAlone = TRUE)

modulePath <- "R/SpaDES/m"
SpaDES.project::getModule(modulePath = modulePath,
                          c("PredictiveEcology/Biomass_speciesData@505fa065399da93e817424373ee1160e46703ce3",
                            "PredictiveEcology/Biomass_borealDataPrep@6cd0c1242cdb95f7432d37becfb2f8dd03642e76",
                            "PredictiveEcology/Biomass_core@5f7691af755a651579408f37bf292a7274b7678f",
                            "PredictiveEcology/Biomass_validationKNN@ec9b5362aa0d07eee844fcdb74e33abeeda89b4b",
                            "CeresBarros/Biomass_speciesParameters@31b66f22f915cd8b5774e5b56359dc2406691895"))

pkgSnapshotFile <- file.path("packages",
                             paste0("pkgSnapshot_",
                                    paste0(version$major, "_", strsplit(version$minor, "[.]")[[1]][1]),
                                    ".txt"))

if (file.exists(pkgSnapshotFile)) {
  Require::Require(packageVersionFile = pkgSnapshotFile,
                   require = FALSE, standAlone = TRUE)
} else {
  outs <- SpaDES.project::packagesInModules(modulePath = modulePath)
  Require::Require(c(unname(unlist(outs)),
                     "PredictiveEcology/SpaDES.experiment@91bfad98d67ea2b7fcee3ea0115f8746e47534ad",
                     "devtools", "ggspatial", "ggpubr", "cowplot"),
                   require = FALSE, standAlone = TRUE)

  ## the next line is commented to prevent accidental execution, but can be run to generate a missing pkg snapshot
  Require::pkgSnapshot(pkgSnapshotFile, libPaths = pkgPath,
                       standAlone = TRUE, exact = TRUE, purge = TRUE,
                       includeBase = FALSE)
}

## load packages
Require::Require(c("raster", "terra", "dplyr", "data.table", "future",
                   "SpaDES.core", "LandR", "reproducible",
                   "ggspatial", "ggpubr", "cowplot"),
                 upgrade = FALSE, install = FALSE)

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------

## Set up modelling parameters  ---------------------------
options("reproducible.useNewDigestAlgorithm" = 2)
options("spades.moduleCodeChecks" = FALSE)
options("spades.inputPath" = Require::normPath(file.path("R/SpaDES/inputs")))  ## store everything in inputs/ so that there are no duplicated files across modules
options("reproducible.useCache" = TRUE)
options("reproducible.destinationPath" = Require::normPath(file.path("R/SpaDES/inputs")))
options("reproducible.useGDAL" = FALSE)
options("spades.useRequire" = FALSE)
options("Require.unloadNamespaces" = FALSE)
options("reproducible.useTerra" = FALSE)

eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## TEST STUDY AREAS
## use studyAreaS or L to parameterise AND run simulations in the SAME area
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
simPathName <- "oct2022Runs"
simPaths <- list(cachePath = Require::normPath(file.path("R/SpaDES/cache", simPathName))
                 , modulePath = Require::normPath(modulePath)
                 , inputPath = Require::normPath(file.path("R/SpaDES/inputs"))
                 , outputPath = Require::normPath(file.path("R/SpaDES/outputs", simPathName, runName)))

figPath <- "R/SpaDES/outputs/GeneralFigs"
dir.create(figPath, recursive = TRUE)

## prompt googledrive authorization - you can replace NULL with your Google Account email, to run non-interactively
GAemail <- NULL
googledrive::drive_auth(email = GAemail)

## Get necessary objects like the study area.
x11() ## open a new plotting window - avoids errors if the current one is too small.
devtools::source_url(paste0("https://raw.githubusercontent.com/CeresBarros/",
                            "LandRBiomass_publication/repPkgInstall/R/SpaDES/",
                            "1_simObjects.R?raw=TRUE"))
## Run Biomass_speciesData to get species layers
## running this separately from other modules makes switching
## between using a large and a smaller study area easier when the smaller SA is within the large one,
## as it keeps the data in separate folders that can be used across simulations/scenarios
devtools::source_url(paste0("https://raw.githubusercontent.com/CeresBarros/",
                            "LandRBiomass_publication/repPkgInstall/R/SpaDES/",
                            "2_speciesLayers.R?raw=TRUE"))

## check species layers:
plot(simOutSpeciesLayers$speciesLayers)

## subset sppEquivalencies and colorVector
sppEquivalencies_CA <- sppEquivalencies_CA[Boreal %in% names(simOutSpeciesLayers$speciesLayers)]
sppColorVect <- sppColorVect[c(names(simOutSpeciesLayers$speciesLayers), "Mixed")]

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
    , "fireURL" = "https://drive.google.com/file/d/1YIc_BSkPKqW60SmfpR2vDpeRGrwOFKso/view?usp=sharing"  ## use a frozen version of fire perimeter data
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
  , "rasterToMatchLarge" = simOutSpeciesLayers$rasterToMatchLarge
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

## Make a rasterToMatch now that we have a rasterToMatchLarge and a studyArea -- use terra here.
RTM <- try(project(rast(simObjects$rasterToMatchLarge), y = crs(vect(simObjects$studyArea))))
RTM <- crop(RTM, vect(simObjects$studyArea), mask = TRUE)
RTM <- raster(RTM)
RTM[!is.na(RTM[])] <- 1L

simObjects$rasterToMatch <- RTM

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
  # plan("multisession", workers = 2) ## not working, options not passed to futures
  plan("sequential")
} else {
  plan("multicore", workers = 2)
}

LandRBiomass_sim <- SpaDES.experiment::experiment2(
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

validationPaths <- list(cachePath = Require::normPath(file.path("R/SpaDES/cache", simPathName))
                        , modulePath = Require::normPath(file.path("R/SpaDES/m"))
                        , inputPath = Require::normPath(file.path("R/SpaDES/inputs"))
                        , outputPath = Require::normPath(file.path("R/SpaDES/validation", simPathName, runName)))

devtools::source_url(paste0("https://raw.githubusercontent.com/CeresBarros/",
                            "LandRBiomass_publication/repPkgInstall/R/SpaDES/",
                            "3_simObjects4Valid.R?raw=TRUE"))

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
LandRBiomass_validation <- Cache(simInitAndSpades
                                 , times = validationTimes
                                 , params = validationParams
                                 , modules = "Biomass_validationKNN"
                                 , objects = validationObjects
                                 , outputs = validationOutputs
                                 , paths = validationPaths
                                 , .studyAreaName = SAname
                                 , userTags = c("validation", SAname)
                                 , cacheRepo = validationPaths$cachePath
                                 , omitArgs = c("userTags"))

saveSimList(LandRBiomass_validation, file.path(validationPaths$outputPath, paste0("simValid", runName, ".qs")))   ## only save in first runs

## -----------------------------------------------
## POST-HOC ANALYSIS - figures
## -----------------------------------------------

## can only be run after all simulations are complete
if (all(c("baseCase", "studyAreaChange", "altParameters") %in% dir(file.path(validationPaths$outputPath, "..")))) {
  devtools::source_url(paste0("https://raw.githubusercontent.com/CeresBarros/",
                              "LandRBiomass_publication/repPkgInstall/R/SpaDES/",
                              "pubFigures.R?raw=TRUE"))
} else {
  warning("Not all simulations were validated.")
}

q("no")

