## ------------------------------------------------------
## LandR Biomass PUBLICATION SIMULATIONS
##
## Ceres: June 2018
## ------------------------------------------------------

## clean workspace
rm(list=ls()); amc::.gc()

## requires as of Jan 22 2020
# reproducible 1.2.5.9000
# quickPlot_0.1.7.9001
# SpaDES.core 1.0.5
# SpaDES.tools_0.3.7
# SpaDES.addins_0.1.2
# SpaDES.experiment 0.0.2.9000
# LandR_0.0.11.9008 ## development branch
library(SpaDES)
library(SpaDES.experiment)
library(raster)
library(LandR)
library(dplyr)
library(data.table)

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

speciesParams <- list(
  "shadetolerance" = list(
    Betu_Pap = 1.5,
    Lari_Lar = 1.5,
    Pice_Gla = 2,
    Pice_Mar = 3,
    Pinu_Ban = 1,
    Popu_Spp = 1
  )
)

simModules <- list("Biomass_borealDataPrep"
                   , "Biomass_core"
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
    ,"speciesUpdateFunction" = list(
      quote(LandR::speciesTableUpdate(sim$species, sim$speciesTable, sim$sppEquiv, P(sim)$sppEquivCol)),
      quote(LandR::updateSpeciesTable(speciesTable = sim$species, params = sim$speciesParams)))
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
)

## Run Biomass_speciesData to get species layers
## running this separately from other modules makes switching
## between using a large and a smaller study area easier when the smaller SA is within the large one,
## as it keeps the data in separate folders that can be used across simulations/scenarios
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
simOutputs <- rbind(simOutputs, data.frame(objectName = "vegTypeMap",
                                           saveTime = unique(seq(simTimes$start, simTimes$end, by = 5)),
                                           eventPriority = 10))
## other outputs needed to run validation:
simOutputs <- rbind(simOutputs, data.frame(objectName = "biomassMap",
                                           saveTime = simTimes$start,
                                           eventPriority = 1))

## in the first year, objects have to be saved after init events, before mortalityAndGrowth
## note that vegTypeMap won't be saved at yr 0, because it doens't exist at priority 5.5
## and save is not scheduled twice (even if we changed the priority)
## (not doing this now and simply acknowledging that cohortData will show variation across reps in year 1)
# simOutputs[simOutputs$saveTime == 0, "eventPriority"] <- 5.5

simObjects <- list(
  "sppEquiv" = sppEquivalencies_CA
  , "sppColorVect" = sppColorVect
  , "speciesLayers" = simOutSpeciesLayers$speciesLayers
  , "speciesParams" = speciesParams
  , "treed" = simOutSpeciesLayers$treed
  , "numTreed" = simOutSpeciesLayers$numTreed
  , "nonZeroCover" = simOutSpeciesLayers$nonZeroCover
  , "PSPgis" = PSPgis
  , "PSPmeasure" = PSPmeasure
  , "PSPplot" = PSPplot
)

if (runName == "parametriseSALarge") {
  simObjects$studyArea <- studyAreaS
  simObjects$studyAreaLarge <- studyAreaL
} else {
  simObjects$studyArea <- get(runName)
}

# reproducible::clearCache(simPaths$cachePath)
LandRBiomass_sim <- simInit(times = simTimes
                            , params = simParams
                            , modules = simModules
                            , objects = simObjects
                            , outputs = simOutputs
                            , paths = simPaths)

saveRDS(LandRBiomass_sim,
        file.path(simPaths$outputPath, paste0("simInitList_", runName, ".rds")))

## just one rep - test sim
# LandRBiomass_simTest <- spades(LandRBiomass_sim, .plotInitialTime = simTimes$start)

## SIMULATION WITH 10 REPS
# options("reproducible.useCache" = "overwrite")
library(future)
plan("multiprocess", workers = 1)   ## each worker consumming roughly 15Gb.
factorialSimulations <- experiment2(
  sim1 = LandRBiomass_sim,
  clearSimEnv = TRUE,
  replicates = 10)

saveRDS(factorialSimulations, file.path(simPaths$outputPath, paste0("simList_factorialSimulations_", runName, ".rds")))
## VALIDATION
## get the  land-cover change map (needed to have an RTM first, so get it from the simInitList)
## /!\ it is assumed that the filename of the raster in the simList corresponds to the raster found in disk.
## this may not be the case if the simulations were run in another machine and saved rasters were not imported.

## make objects again in case only this part of the script is being run:
if (!exists("simDirName"))
  simDirName <- "dec2020Runs"

if (!exists("simPaths"))
  simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName)
                   , modulePath = file.path("R/SpaDES/m")
                   , inputPath = file.path("R/SpaDES/inputs")
                   , outputPath = file.path("R/SpaDES/outputs", simDirName))

validationPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName)
                        , modulePath = file.path("R/SpaDES/m")
                        , inputPath = file.path("R/SpaDES/inputs")
                        , outputPath = file.path("R/SpaDES/validation", simDirName))

source("R/SpaDES/3_simObjects4Valid.R")

## PARAMETERS FOR VALIDATION MODULE
## in this case all reps have the same parameters, so we can use the first rep to get the values
validationTimes <- list(start = 1, end = 1)
validationParams <- list(
  Biomass_validationKNN = list(
    "minCoverThreshold" = params(factorialSimulations$sim1_rep01)$Biomass_borealDataPrep$minCoverThreshold
    , "pixelGroupBiomassClass" = params(factorialSimulations$sim1_rep01)$Biomass_borealDataPrep$pixelGroupBiomassClass
    , "deciduousCoverDiscount" = params(factorialSimulations$sim1_rep01)$Biomass_borealDataPrep$deciduousCoverDiscount
    , "sppEquivCol" = params(factorialSimulations$sim1_rep01)$Biomass_borealDataPrep$sppEquivCol
    , "validationReps" = as.integer(1:10)  ## or length of simLists
    , "validationYears" = as.integer(c(simTimes$start, 2011))
    , ".useCache" = eventCaching
  )
)

validationObjects <- list(
  "biomassMap" = rawBiomassMap  ## to change when outputs are ready
  , "rasterToMatch" = rasterToMatch  ## it's retrieved from LandRBiomass_sim in the sourced script
  , "rawBiomassMapStart" = rawBiomassMap
  , "rstLCChange" = rstLCChangeAllbin
  , "simulationOutputs" = simulationOutputs
  , "speciesLayersStart" = speciesLayers
  , "sppColorVect" = LandRBiomass_sim$sppColorVect
  , "sppEquiv" = LandRBiomass_sim$sppEquiv
  , "standAgeMapStart" = standAgeMap
  , "studyArea" = LandRBiomass_sim$studyArea
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
devtools::load_all("../LandR/")
LandRBiomass_sim <- simInitAndSpades(times = validationTimes
                                     , params = validationParams
                                     , modules = "Biomass_validationKNN"
                                     , objects = validationObjects
                                     , outputs = validationOutputs
                                     , paths = )   ## note that simPaths must respect the simulation outputPath

q("no")

