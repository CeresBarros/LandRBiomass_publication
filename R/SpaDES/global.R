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
.libPaths(pkgDir)

## maybe restart with command here?

if (!require("remotes")) {   ## had to remove libloc - that was what was causing the restart because remotes is loaded thanks to the ::
  install.packages("remotes", lib = pkgDir)
  library(remotes)
}

remotes::install_github("PredictiveEcology/SpaDES.install@a120e301eed8e427dc10b1ac2cdb34156b50af28", lib = pkgDir, upgrade = FALSE) ## install this first, it brings Require along
remotes::install_github("PredictiveEcology/Require@6f2276ca64eee8365045363ccb8fc6f7935225a8", lib = pkgDir, upgrade = FALSE) ## change to sha

options("Require.unloadNamespaces" = FALSE)
Require::Require("checkpoint")
checkpoint("2022-06-01", checkpoint_location = pkgDir, r_version = "4.2.1",
           scan_now = FALSE, scan_rprofile = FALSE)
.libPaths(c(.libPaths(), pkgDir))    ## we need SpaDES.install, and not a bad idea to take advantage of other installed deps. This also seems to prevent issus with Rcpp/rgdal/magrittr install

## SpaDES/LandR pkg installation:
install.packages(c("magrittr", "pryr", "RCurl", "XML"))   ## need to be separate, too finicky and fail to install in subsequent calls
Require::Require(c("SpaDES"), require = FALSE)

SpaDES.install::makeSureAllPackagesInstalled(modulePath = "R/SpaDES/m")

## some packages need specific versions -- a restart may be needed, if so rerun lines 1-58
## make sure the appropriate RTools is installed and the path added to your .Renviron (see https://stackoverflow.com/a/71751606/11969696)
remotes::install_github("PredictiveEcology/reproducible@ed25a8024d71b4ffa931a2e58e7a919257eec7e1")
remotes::install_github("PredictiveEcology/LandR@65f3b8e8da8df22f3cf0168e2d505d5da49067c8")
remotes::install_github("ianmseddy/LandR.CS@02b5610366f1cac53011a72a43805a906c2437e0")
remotes::install_github("ianmseddy/PSPclean@3c3f0e7082e14c111a607c3ba803abf0396343e6")
remotes::install_github("PredictiveEcology/SpaDES.experiment@5a23c40f8aa9a9efc6dc16e040f8771561059152")

## after installing everything, please restart R again.
## when you do so rerun lines 1-58

# install.packages("SpaDES", dependencies = TRUE, lib = pkgDir) # yes restart if need be, yes install from soruce
# remotes::install_github("PredictiveEcology/SpaDES.core@8a7886a6afd7f3b90df10ea6b87caae8661f8709", lib = pkgDir)
# remotes::install_github("PredictiveEcology/LandR@093c39898912a6e89ac9b6e862733052a0fae407", lib = pkgDir)
# remotes::install_github("ianmseddy/LandR.CS@b39c8c72d20189fa6b6aeb057cdc751a631e0efa", lib = pkgDir)
# remotes::install_github("PredictiveEcology/reproducible@aedea49637a6ebd0db6897f1d33f53959f41bee2", lib = pkgDir)
# remotes::install_github("PredictiveEcology/SpaDES.install@80c43dcb94d897d25545105a7b83111cf634a556", lib = pkgDir)
# remotes::install_github("PredictiveEcology/SpaDES.experiment@5a23c40f8aa9a9efc6dc16e040f8771561059152", lib = pkgDir)

# Require::pkgSnapshot("packages/pkgSnapshot.txt", libPaths = "packages/x86_64-w64-mingw32/4.0/")
# Much later on a different or same machine:
# options("Require.unloadNamespaces" = FALSE)
# Require::Require(packageVersionFile = "packages/pkgSnapshot.txt", libPaths = pkgDir)
# Require::Require(packageVersionFile = "packages/pkgSnapshot.txt", libPaths = pkgDir)  ## run a second time -- some pkgs may fail to install the first time around
# sink(type = "message")
# sink()

# isRstudio <- Sys.getenv("RSTUDIO") == 1 ||
#   .Platform$GUI == "RStudio" ||
#   if (suppressWarnings(requireNamespace("rstudioapi", quietly = TRUE))) {
#     rstudioapi::isAvailable()
#   } else {
#     FALSE
#   }
# if (isRstudio)
#   rstudioapi::restartSession(command = paste0(".libPaths('", pkgDir, "')"))   ## this is not fixing the pkg loading at start-up issue

## the previous line may result in a few errors on Windows when installing packages from scratch
## PLEASE RESTART R HERE, and re-run all code up to this point again. All errors should be resolved and
## any packages missing should be installed the second time around.

## Try running sim without this:
# out <- SpaDES.install::makeSureAllPackagesInstalled(modulePath = "R/SpaDES/m")
## the previous line may result in a few errors on Windows when installing packages from scratch
## PLEASE RESTART R HERE, and re-run all code up to this point again. All errors should be resolved and
## any packages missing should be installed the second time around.

## load packages and make sure minimum versions are installed.
Require::Require(c("SpaDES",
                   "raster", "dplyr", "data.table", "future",
                   "SpaDES.experiment",
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
                 , modulePath = file.path("R/SpaDES/m")
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
# httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))
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
  Biomass_borealDataPrep = list(
    "sppEquivCol" = sppEquivCol
    , "fitDeciduousCoverDiscount" = TRUE
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
    , ".useCache" = eventCaching
  )
  , Biomass_speciesParameters = list(
    "quantileAgeSubset" = list(Betu_Pap = 95, Lari_Lar = 95, Pice_Gla = 95, Pice_Mar = 95, Pinu_Ban = 99, Popu_Spp = 99)
    , "sppEquivCol" = sppEquivCol
    , ".useCache" = eventCaching
  )
  , Biomass_core = list(
    "calcSummaryBGM" = c("start")
    , "initialBiomassSource" = "cohortData"
    , "plotOverstory" = TRUE
    , "seedingAlgorithm" = "wardDispersal"
    , "sppEquivCol" = sppEquivCol
    , "successionTimestep" = successionTimestep
    , "vegLeadingProportion" = vegLeadingProportion
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

