---
title: "LandR Publication Simulations"
author: "Ceres Barros"
date: "22/01/2021"
output:
  html_document:
    keep_md: yes
editor_options:
  chunk_output_type: console
---




```r
library(Require)
Require(c(
  "data.table", "dplyr", "ggplot2",
  "quickPlot", "raster", "sp", "SpaDES",
  "PredictiveEcology/SpaDES.experiment",
  "PredictiveEcology/LandR",
))

options("reproducible.useNewDigestAlgorith" = 2)
options("spades.moduleCodeChecks" = FALSE)
options("reproducible.useCache" = TRUE)
options("reproducible.inputPaths" = file.path("R/SpaDES/inputs"))
options("reproducible.destinationPath" = file.path("R/SpaDES/inputs"))
options("reproducible.useGDAL" = FALSE)

## set run name and paths
runName <- "parametriseSALarge"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths
simDirName <- "test"
simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName)
                 , modulePath = file.path("R/SpaDES/m")
                 , inputPath = file.path("R/SpaDES/inputs")
                 , outputPath = file.path("R/SpaDES/outputs", simDirName))
```

## Simulation setup - part 1 

### Get study area and other necessary objects 


```r
## Get necessary objects -----------------------
source("R/SpaDES/1_simObjects.R")
```

## Simulation setup - part 2 - species layers


```r
## Run Biomass_speciesData to get species layers
## running this separately from other modules makes switching
## between using a large and a smaller study area easier when the smaller SA is within the large one,
## as it keeps the data in separate folders thatn can be used across simulations/scenarios
source("R/SpaDES/2_speciesLayers.R")

## check species layers:
# plot(simOutSpeciesLayers$speciesLayers)
## Populus grandidentata shouldn't be in SK (and has v. few pixels in the layer) and will be excluded "manually"
toRm <- which(names(simOutSpeciesLayers$speciesLayers) %in% c("Popu_Gra"))
simOutSpeciesLayers$speciesLayers <- dropLayer(simOutSpeciesLayers$speciesLayers, i = toRm)
rm(toRm)
```

## Simulation setup - part 3 - module parameters and outputs

* `vegLeadingProportion` indicates what proportion the stand must be in one species group for it to be leading.
  If 0, then there is never a notion of "mixed" vegetation types and a species is a leading species if it has the highest relative biomass in the pixel.
* `successionTimestep` defines the frequency at which dispersal and age reclassification occurs - every 10 years is the default LANDIS behaviour. 


```r
## simulation params
simTimes <- list(start = 2001, end = 2031)
vegLeadingProportion <- 0 
successionTimestep <- 10L  

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

simModules <- list("Biomass_borealDataPrep"
                   , "Biomass_speciesParameters"
                   , "Biomass_core"
)

simParams <- list(
  Biomass_borealDataPrep = list(
    "sppEquivCol" = sppEquivCol
    , "forestedLCCClasses" = c(1:15, 34:35)
    , "LCCClassesToReplaceNN" = c(34:35)
    , "fitDeciduousCoverDiscount" = TRUE
    , "subsetDataAgeModel" = FALSE
    , "subsetDataBiomassModel" = FALSE
    , "biomassModel" =  quote(lme4::lmer(B ~ logAge * speciesCode + cover * speciesCode +
                                           (logAge + cover | ecoregionGroup)))
    , "exportModels" = "all"
    , "fixModelBiomass" = TRUE
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
  , Biomass_speciesParameters = list(
    "sppEquivCol" = sppEquivCol
    , ".useCache" = eventCaching
  )
  , Biomass_core = list(
    "calcSummaryBGM" = c("start")
    , "initialBiomassSource" = "cohortData" 
    , ".plotInitialTime" = simTimes$start
    , "plotOverstory" = TRUE
    , "seedingAlgorithm" = "wardDispersal"
    , "sppEquivCol" = sppEquivCol
    , "successionTimestep" = successionTimestep
    , "vegLeadingProportion" = vegLeadingProportion
    , ".plotInterval" = 1
    , ".plotMaps" = TRUE
    , ".saveInitialTime" = NA
    , ".useCache" = eventCaching[1]
    , ".useParallel" = useParallel
  )
)

## subset sppEquivalencies
sppEquivalencies_CA <- sppEquivalencies_CA[Boreal %in% names(simOutSpeciesLayers$speciesLayers)]

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
simOutputs <- rbind(simOutputs, data.frame(objectName = "speciesEcoregion",
                                           saveTime = simTimes$start,
                                           eventPriority = 1))
simOutputs <- rbind(simOutputs, data.frame(objectName = "species",
                                           saveTime = simTimes$start,
                                           eventPriority = 1))

## in the first year, eventPriorities need to be set to AFTER the init event (which has priority 1)
simOutputs$eventPriority[simOutputs$saveTime == simTimes$start] <- 1.5

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

simObjects$studyAreaLarge <- studyAreaL

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
                              , .plotInitialTime = NA
                              , userTags = "simInitAndInits"
                              , omitArgs = c("userTags", ".plotInitialTime"))

## save the simList
saveSimList(LandRBiomass_simInit, file.path(simPaths$outputPath, paste0("simInit", runName)))
```

## Run simulation

Here we run just one repetition


```r
LandRBiomass_sim <- spades(LandRBiomass_simInit, .plotInitialTime = simTimes$start)
```

If we were to run several repetitions, this would be how:


```r
library(future)
plan("multiprocess", workers = 1)   ## each worker consumes roughly 10Gb
LandRBiomass_sim <- experiment2(
  sim1 = LandRBiomass_sim,
  clearSimEnv = TRUE,
  replicates = 10)
```

## Inspect simulation objects

We can use the `simList` objects to plot simulation objects, such as the input layers used for parameterisation.


```r
## study area within Saskatchewan province
## get Canadian provinces and subset to SK
can1 <- raster::getData('GADM', country = "CAN", level = 1, path = tempdir())
can1 <- can1[can1$NAME_1 == "Saskatchewan",]
can1 <- spTransform(can1, CRSobj = crs(LandRBiomass_simInit$studyArea))

plotStudyAreas <- ggplot() +
  layer_spatial(can1, fill = "grey90") + 
  layer_spatial(LandRBiomass_simInit$studyAreaLarge, fill = "darkgreen") +
  layer_spatial(LandRBiomass_simInit$studyArea, fill = "green") +
  labs(fill = "") +
  theme_void()

clearPlot(force = TRUE)
Plot(plotStudyAreas, title = " ")

## input stand biomass and age, ecological zonation (ecodistricts) and land-cover (LCC 2005)
plotEcodist <- ggplot() +
  layer_spatial(LandRBiomass_simInit$ecoregionLayer, aes(fill = as.factor(ECODISTRIC))) + 
  labs(fill = "") +
  theme_void()

clearPlot(force = TRUE)
Plot(LandRBiomass_simInit$rawBiomassMap,
     LandRBiomass_simInit$standAgeMap, 
     LandRBiomass_simInit$rstLCC,
     plotEcodist,
     title = c("kNN stand biomass", "kNN stand age", "land-cover", "ecodistricts"))

## species % cover
clearPlot(force = TRUE)
Plot(simOutSpeciesLayers$speciesLayers,
     title = names(simOutSpeciesLayers$speciesLayers))
```

Similarly, we can have a look at the species traits values used in the simulation directly from the `simList` object (although we also chose to save them).


```r
LandRBiomass_simInit$species
LandRBiomass_simInit$minRelativeB ## can be spatially varying, but identical across ecolocation (AKA ecoregion)
LandRBiomass_simInit$sufficientLight
LandRBiomass_simInit$speciesEcoregion
```

For example, these were the (spatially) invariant species traits used in the simulation:



# Validation

Here we run the validation on the outputs of just one repetition.

We begin by preparing all the inputs necessary for the `Biomass_validationKNN` module.


```r
## get the  land-cover change map (needed to have an RTM first, so get it from the simInitList)
## /!\ it is assumed that the filename of the raster in the simList corresponds to the raster found in disk.
## this may not be the case if the simulations were run in another machine and saved rasters were not imported.

## make objects again in case only this part of the script is being run:
if (!exists("simDirName"))
  simDirName <- "test"

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
    "minCoverThreshold" = params(LandRBiomass_simInit)$Biomass_borealDataPrep$minCoverThreshold
    , "pixelGroupBiomassClass" = params(LandRBiomass_simInit)$Biomass_borealDataPrep$pixelGroupBiomassClass
    , "deciduousCoverDiscount" = params(LandRBiomass_simInit)$Biomass_borealDataPrep$deciduousCoverDiscount
    , "sppEquivCol" = params(LandRBiomass_simInit)$Biomass_borealDataPrep$sppEquivCol
    , "validationReps" = as.integer(1:10)  ## or length of simLists
    , "validationYears" = as.integer(c(2001, 2011))
    , ".useCache" = eventCaching
  )
)

validationObjects <- list(
  "biomassMap" = biomassMap
  , "rasterToMatch" = rasterToMatch
  , "rawBiomassMapStart" = rawBiomassMap
  , "rstLCChange" = rstLCChangeAllbin
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
```

We can now run the validation:


```r
LandRBiomass_validation <- simInitAndSpades(times = validationTimes
                                            , params = validationParams
                                            , modules = "Biomass_validationKNN"
                                            , objects = validationObjects
                                            , outputs = validationOutputs
                                            , paths = validationPaths
                                            , .plotInitialTime = validationTimes$start) 
```