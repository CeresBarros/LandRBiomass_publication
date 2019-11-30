## ------------------------------------------------------
## LandR Biomass PUBLICATION SIMULATIONS
##
## Ceres: June 2018
## ------------------------------------------------------

## clean workspace
rm(list=ls()); amc::.gc()
stopifnot(utils::packageVersion("googledrive") == "1.0.0")

## requires as of June 10th 2019
# loading reproducible     0.2.8.9001
# loading quickPlot        0.1.6.9000
# loading SpaDES.core      0.2.5.9004
# loading SpaDES.tools     0.3.2.9000
# loading SpaDES.addins    0.1.2

library(SpaDES)
library(SpaDES.experiment)
library(LandR)
library(raster)

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------

## Get necessary objects -----------------------
source("R/SpaDES/1_simObjects.R")

## Set up modelling parameters  ---------------------------
options('reproducible.useNewDigestAlgorithm' = TRUE)
# runName <- "studyAreaS"
# runName <- "studyAreaL"
runName <- "parametriseSALarge"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths
simPaths <- list(cachePath = file.path("R/SpaDES/cache", runName),
                 modulePath = file.path("R/SpaDES/m"),
                 inputPath = file.path("R/SpaDES/inputs"),
                 outputPath = file.path("R/SpaDES/outputs", runName))
## simulation params
simTimes <- list(start = 0, end = 30)
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
successionTimestep <- 10L  # for dispersal and age reclass.

simModules <- list("Biomass_borealDataPrep"
                   , "Biomass_core"
                   , "Biomass_validationKNN"
)

simParams <- list(
  Biomass_borealDataPrep = list(
    "sppEquivCol" = sppEquivCol
    , "forestedLCCClasses" = c(1:15, 34:36)
    # next two are used when assigning pixelGroup membership; what resolution for
    #   age and biomass
    , "pixelGroupAgeClass" = successionTimestep
    , "pixelGroupBiomassClass" = 100
    , "useCloudCacheForStats" = FALSE
    , "cloudFolderID" = NA
    , ".useCache" = eventCaching
  )
  , Biomass_core = list(
    "calcSummaryBGM" = c("start")
    , "initialBiomassSource" = "cohortData" # can be 'biomassMap' or "spinup" too
    # , ".plotInitialTime" = simTimes$start
    , ".plotInitialTime" = NA ## for experiment
    , "plotOverstory" = TRUE
    , "seedingAlgorithm" = "wardDispersal"
    , "sppEquivCol" = sppEquivCol
    , "successionTimestep" = successionTimestep
    , "vegLeadingProportion" = vegLeadingProportion
    , ".plotInterval" = 1
    , ".plotMaps" = TRUE
    , ".saveInitialTime" = NA
    # , ".useCache" = eventCaching
    , ".useCache" = eventCaching[1] # experiment doesn't like when init is cached
    , ".useParallel" = useParallel
  )
  , Biomass_validationKNN = list(
    "sppEquivCol" = sppEquivCol
    , ".useCache" = eventCaching
  )
)

## Run Biomass_speciesData to get species layers
source("R/SpaDES/2_speciesLayers.R")

## check species layers:
plot(sim$speciesLayers)
## Pinus contorta and Populus grandidentata should not be in SK, and will be excluded
## Popu_sp has v. few (one?) pixels and should also be dropeed
toRm <- which(names(simOutSpeciesLayers$speciesLayers) %in% c("Pinu_Con_Lat", "Popu_Gra", "Popu_Spp"))
simOutSpeciesLayers$speciesLayers <- dropLayer(simOutSpeciesLayers$speciesLayers, i = toRm)
rm(toRm)

## subset sppEquivalencies
sppEquivalencies_CA <- sppEquivalencies_CA[Boreal %in% names(simOutSpeciesLayers$speciesLayers)]

simOutputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                  saveTime = c(0:15, seq(20, 100, 5)),
                                  eventPriority = 10,
                                  stringsAsFactors = FALSE))
simOutputs[1, "eventPriority"] <- 5.5  ## in the beggining, after init events, before mortalityAndGrowth
simOutputs <- rbind(simOutputs, data.frame(objectName = "vegTypeMap",
                                     saveTime = c(0:15, seq(20, 100, 5)),
                                     eventPriority = 10))
simOutputs <- rbind(simOutputs, data.frame(objectName = "pixelGroupMap",
                                     saveTime = c(0:15, seq(20, 100, 5)),
                                     eventPriority = 10))
simOutputs <- rbind(simOutputs, data.frame(objectName = "rstDisturbedPix",
                                           saveTime = c(0),
                                           eventPriority = 10))
simOutputs <- rbind(simOutputs, data.frame(objectName = "rawBiomassMapValidation",
                                           saveTime = c(0),
                                           eventPriority = 10))
simOutputs <- rbind(simOutputs, data.frame(objectName = "standAgeMapValidation",
                                           saveTime = c(0),
                                           eventPriority = 10))
simOutputs <- rbind(simOutputs, data.frame(objectName = "speciesLayersValidation",
                                           saveTime = c(0),
                                           eventPriority = 10))

if (runName == "parametriseSALarge") {
  simObjects <- list("studyArea" = studyAreaS
                     , "studyAreaLarge" = studyAreaL
                     , "sppEquiv" = sppEquivalencies_CA
                     , "sppColorVect" = sppColorVect
                     , "speciesLayers" = simOutSpeciesLayers$speciesLayers
                     , "treed" = simOutSpeciesLayers$treed
                     , "numTreed" = simOutSpeciesLayers$numTreed
                     , "nonZeroCover" = simOutSpeciesLayers$nonZeroCover
  )
} else {
  simObjects <- list("studyArea" = if (grepl("study", runName)) get(runName) else studyAreaS
                     , "sppEquiv" = sppEquivalencies_CA
                     , "sppColorVect" = sppColorVect
                     , "speciesLayers" = simOutSpeciesLayers$speciesLayers
                     , "treed" = simOutSpeciesLayers$treed
                     , "numTreed" = simOutSpeciesLayers$numTreed
                     , "nonZeroCover" = simOutSpeciesLayers$nonZeroCover
  )
}

startTime <- date()
options(spades.moduleCodeChecks = FALSE)
options("reproducible.useCache" = TRUE)
# reproducible::clearCache(simPaths$cachePath)
Biomass_core_testSim <- simInit(times = simTimes
                                , params = simParams
                                , modules = simModules
                                , objects = simObjects
                                , outputs = simOutputs
                                , paths = simPaths)

## to avoid synonim bug run the spades call once for 1 year.
end(Biomass_core_testSim) <- 0
spades(Biomass_core_testSim
       , debug = TRUE
       , .plotInitialTime = NA)
end(Biomass_core_testSim) <- 30   ## now change back for experiment.
# saveRDS(Biomass_core_testSim, file.path(simPaths$outputPath, paste0("simList_", runName, ".rds")))
# unlink(file.path(simPaths$outputPath, "figures"), recursive = TRUE) ## remove unnecessary figures


library(future)
plan("multiprocess", workers = 6)
factorialSimulations <- experiment2(
  sim1 = Biomass_core_testSim,
  clearSimEnv = TRUE,
  replicates = 10)

saveRDS(factorialSimulations, file.path(simPaths$outputPath, paste0("simList_factorialSimulations", runName, ".rds")))
