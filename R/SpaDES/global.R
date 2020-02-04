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
library(raster)
library(LandR)

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------

## Set up modelling parameters  ---------------------------
options("reproducible.useNewDigestAlgorithm" = TRUE)
options("spades.moduleCodeChecks" = FALSE)
options("reproducible.useCache" = TRUE)
options("reproducible.inputPaths" = normPath("data"))  ## store everything in data/ so that there are no duplicated files across modules
options("reproducible.destinationPath" = normPath("data")) 
# runName <- "studyAreaS"
# runName <- "studyAreaL"
runName <- "parametriseSALarge"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths
simPaths <- list(cachePath = normPath(file.path("R/SpaDES/cache", runName)),
                 modulePath = file.path("R/SpaDES/m"),
                 inputPath = normPath("data"),
                 outputPath = normPath(file.path("R/SpaDES/outputs", runName)))

## Get necessary objects -----------------------
source("R/SpaDES/1_simObjects.R")

## simulation params
simTimes <- list(start = 0, end = 30)
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
successionTimestep <- 10L  # for dispersal and age reclass.

simModules <- list("Biomass_borealDataPrep"
                   , "Biomass_core"
                   , "Biomass_validationKNN"
                   , "LandR_speciesParameters"
)

simParams <- list(
  LandR_speciesParameters = list(
    "sppEquivCol" = sppEquivCol
    , ".useCache" = eventCaching
  ),
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
# plot(simOutSpeciesLayers$speciesLayers)
## Populus grandidentata shouldn't be inin SK (and has only v. few pixels in the layer) and will be excluded
toRm <- which(names(simOutSpeciesLayers$speciesLayers) %in% c("Popu_Gra"))
simOutSpeciesLayers$speciesLayers <- dropLayer(simOutSpeciesLayers$speciesLayers, i = toRm)
rm(toRm)

## subset sppEquivalencies
sppEquivalencies_CA <- sppEquivalencies_CA[Boreal %in% names(simOutSpeciesLayers$speciesLayers)]

## objects will be saved at the start of the simulation (so they reflect the previous year)
simOutputs <- data.frame(expand.grid(objectName = c("cohortData"),
                                     saveTime = c(0:15, seq(20, 100, 5)),
                                     eventPriority = 1,
                                     stringsAsFactors = FALSE))
simOutputs <- rbind(simOutputs, data.frame(objectName = "pixelGroupMap",
                                           saveTime = c(0:15, seq(20, 100, 5)),
                                           eventPriority = 1))
simOutputs <- rbind(simOutputs, data.frame(objectName = "rstDisturbedPix",
                                           saveTime = c(0),
                                           eventPriority = 1))
simOutputs <- rbind(simOutputs, data.frame(objectName = "vegTypeMap",
                                           saveTime = c(0:15, seq(20, 100, 5)),
                                           eventPriority = 1))
## the following objects are only saved once at the end of year 0/beggining of year 1 (they don't change)
simOutputs <- rbind(simOutputs, data.frame(objectName = "rawBiomassMapValidation",
                                           saveTime = c(1),
                                           eventPriority = 1))
simOutputs <- rbind(simOutputs, data.frame(objectName = "standAgeMapValidation",
                                           saveTime = c(1),
                                           eventPriority = 1))
simOutputs <- rbind(simOutputs, data.frame(objectName = "speciesLayersValidation",
                                           saveTime = c(1),
                                           eventPriority = 1))

## in the first year, objects have to be saved after init events, before mortalityAndGrowth
## note that vegTypeMap won't be saved at yr 0, because it doens't exist at priority 5.5
## and save is not scheduled twice (even if we changed the priority)
## (not doing this now and simply acknowledging that cohortData will show variation across reps in year 1)
# simOutputs[simOutputs$saveTime == 0, "eventPriority"] <- 5.5

## get the  land-cover change map (needed to have an RTM first.)
source("R/SpaDES/3_simObjects4Valid.R")

simObjects <- list("sppEquiv" = sppEquivalencies_CA
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

saveRDS(Biomass_core_testSim,
        file.path(simPaths$outputPath, paste0("simList_", runName, ".rds")))
end(Biomass_core_testSim) <- 30   ## now change back for experiment.
 unlink(file.path(simPaths$outputPath, "figures"), recursive = TRUE) ## remove unnecessary figures

library(future)
plan("multiprocess", workers = 3)
factorialSimulations <- experiment2(
  sim1 = Biomass_core_testSim,
  clearSimEnv = TRUE,
  replicates = 10)

saveRDS(factorialSimulations, file.path(simPaths$outputPath, paste0("simList_factorialSimulations", runName, ".rds")))

q("no")
