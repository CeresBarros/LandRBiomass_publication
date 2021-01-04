## ------------------------------------------------------
## LandR Biomass PUBLICATION SIMULATIONS
##
## Ceres: June 2018
## ------------------------------------------------------

## clean workspace
rm(list=ls()); amc::.gc()

## requires as of Dec 22 2020
# reproducible 1.2.3
# quickPlot_0.1.7.9001
# SpaDES.core 1.0.4
# SpaDES.tools_0.3.7
# SpaDES.addins_0.1.2
# SpaDES.experiment 0.0.2.9000
# LandR_0.0.10.9001 ## dispersalRcpp branch
library(SpaDES)
library(SpaDES.experiment)
library(raster)
library(LandR)

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------

## Set up modelling parameters  ---------------------------
options("reproducible.useNewDigestAlgorithm" = TRUE)
options("spades.moduleCodeChecks" = FALSE)
options("reproducible.useCache" = TRUE)
options("reproducible.inputPaths" = file.path("R/SpaDES/inputs"))  ## store everything in data/ so that there are no duplicated files across modules
options("reproducible.destinationPath" = file.path("R/SpaDES/inputs"))
options("reproducible.useGDAL" = FALSE)
# runName <- "studyAreaS"
# runName <- "studyAreaL"
runName <- "parametriseSALarge"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths
simDirName <- "dec2020Runs"
simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName)
                 , modulePath = file.path("R/SpaDES/m")
                 , inputPath = file.path("R/SpaDES/inputs")
                 , outputPath = file.path("R/SpaDES/outputs", simDirName))

## Get necessary objects -----------------------
source("R/SpaDES/1_simObjects.R")

## simulation params
simTimes <- list(start = 2001, end = 2031)
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
successionTimestep <- 10L  # for dispersal and age reclass.

simModules <- list("Biomass_borealDataPrep"
                   , "Biomass_core"
                   , "Biomass_validationKNN"
                   , "Biomass_speciesParameters"
)

simParams <- list(
  Biomass_speciesParameters = list(
    "sppEquivCol" = sppEquivCol
    , ".useCache" = eventCaching
  ),
  Biomass_borealDataPrep = list(
    "sppEquivCol" = sppEquivCol
    , "forestedLCCClasses" = c(1:15, 34:35)
    , "LCCClassesToReplaceNN" = c(34:35)
    , "fitDeciduousCoverDiscount" = TRUE
    , "exportModels" = "all"
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
# plot(simOutSpeciesLayers$speciesLayers)
## Populus grandidentata shouldn't be inin SK (and has only v. few pixels in the layer) and will be excluded
toRm <- which(names(simOutSpeciesLayers$speciesLayers) %in% c("Popu_Gra"))
simOutSpeciesLayers$speciesLayers <- dropLayer(simOutSpeciesLayers$speciesLayers, i = toRm)
rm(toRm)

## subset sppEquivalencies
sppEquivalencies_CA <- sppEquivalencies_CA[Boreal %in% names(simOutSpeciesLayers$speciesLayers)]

## objects will be saved at the start of the simulation (so they reflect the previous year)
simOutputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                     saveTime = unique(seq(simTimes$start, simTimes$end, by = 5)),
                                     eventPriority = 10,
                                     stringsAsFactors = FALSE))
simOutputs <- rbind(simOutputs, data.frame(objectName = "pixelGroupMap",
                                           saveTime = unique(seq(simTimes$start, simTimes$end, by = 5)),
                                           eventPriority = 10))
simOutputs <- rbind(simOutputs, data.frame(objectName = "rstDisturbedPix",
                                           saveTime = c(simTimes$start),
                                           eventPriority = 1))
simOutputs <- rbind(simOutputs, data.frame(objectName = "vegTypeMap",
                                           saveTime = unique(seq(simTimes$start, simTimes$end, by = 5)),
                                           eventPriority = 10))
## the following objects are only saved once at the end of year 0/beggining of year 1 (they don't change)
simOutputs <- rbind(simOutputs, data.frame(objectName = "rawBiomassMapValidation",
                                           saveTime = c(simTimes$start),
                                           eventPriority = 1))
simOutputs <- rbind(simOutputs, data.frame(objectName = "standAgeMapValidation",
                                           saveTime = c(simTimes$start),
                                           eventPriority = 1))
simOutputs <- rbind(simOutputs, data.frame(objectName = "speciesLayersValidation",
                                           saveTime = c(simTimes$start),
                                           eventPriority = 1))

## in the first year, objects have to be saved after init events, before mortalityAndGrowth
## note that vegTypeMap won't be saved at yr 0, because it doens't exist at priority 5.5
## and save is not scheduled twice (even if we changed the priority)
## (not doing this now and simply acknowledging that cohortData will show variation across reps in year 1)
# simOutputs[simOutputs$saveTime == 0, "eventPriority"] <- 5.5

## get the  land-cover change map (needed to have an RTM first.)
source("R/SpaDES/3_simObjects4Valid.R")

simObjects <- list(
  # "sppNamesVect" = names(simOutSpeciesLayers$speciesLayers)
  "sppEquiv" = sppEquivalencies_CA
  , "sppColorVect" = sppColorVect
  , "speciesLayers" = simOutSpeciesLayers$speciesLayers
  , "treed" = simOutSpeciesLayers$treed
  , "numTreed" = simOutSpeciesLayers$numTreed
  , "nonZeroCover" = simOutSpeciesLayers$nonZeroCover
  , "PSPgis" = PSPgis
  , "PSPmeasure" = PSPmeasure
  , "PSPplot" = PSPplot
  , "rstLCChange" = rstLCChangeAllbin
)

if (runName == "parametriseSALarge") {
  simObjects$studyArea <- studyAreaS
  simObjects$studyAreaLarge <- studyAreaL
} else {
  simObjects$studyArea <- get(runName)
}

startTime <- date()

# reproducible::clearCache(simPaths$cachePath)
LandRBiomass_sim <- simInit(times = simTimes
                            , params = simParams
                            , modules = simModules
                            , objects = simObjects
                            , outputs = simOutputs
                            , paths = simPaths)

saveRDS(LandRBiomass_sim,
        file.path(simPaths$outputPath, paste0("simInitList_", runName, ".rds")))

library(future)
plan("multiprocess", workers = 1)   ## each worker consumming roughly 15Gb.
factorialSimulations <- experiment2(
  sim1 = LandRBiomass_sim,
  clearSimEnv = TRUE,
  replicates = 10)

saveRDS(factorialSimulations, file.path(simPaths$outputPath, paste0("simList_factorialSimulations_", runName, ".rds")))

## save outputfiles table. use lapply for backwards compatibility.
outputFiles <- lapply(factorialSimulations, outputs)
saveRDS(outputFiles, file.path(simPaths$outputPath, paste0("outputFiles_", runName, ".rds")))

q("no")
