## ------------------------------------------------------
## LBMR PUBLICATION SIMULATIONS
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
library(LandR)
library(raster)

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------

## Get necessary objects -----------------------
source("1_simObjects.R")

## Set up modelling parameters  ---------------------------
options('reproducible.useNewDigestAlgorithm' = TRUE)
runName <- "studyAreaS"
# runName <- "studyAreaL"
# runName <- "parametriseSALarge"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths
simPaths <- list(cachePath = file.path("R/SpaDES/cache", runName),
                 modulePath = file.path("R/SpaDES/m"),
                 inputPath = file.path("R/SpaDES/inputs"),
                 outputPath = file.path("R/SpaDES/outputs", runName))
## simulation params
simTimes <- list(start = 0, end = 100)
vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
successionTimestep <- 10L  # for dispersal and age reclass.

simModules <- list("Boreal_LBMRDataPrep"
                   , "LBMR"
)

simParams <- list(
  Boreal_LBMRDataPrep = list(
    "sppEquivCol" = sppEquivCol
    , "forestedLCCClasses" = c(1:15, 34:36)
    # next two are used when assigning pixelGroup membership; what resolution for
    #   age and biomass
    , "pixelGroupAgeClass" = successionTimestep
    , "pixelGroupBiomassClass" = 100
    , "runName" = runName
    , "useCloudCacheForStats" = FALSE
    , "cloudFolderID" = NA
    , ".useCache" = eventCaching
  )
  , LBMR = list(
    "calcSummaryBGM" = c("start")
    , "initialBiomassSource" = "cohortData" # can be 'biomassMap' or "spinup" too
    , ".plotInitialTime" = simTimes$start
    , "seedingAlgorithm" = "wardDispersal"
    , "sppEquivCol" = sppEquivCol
    , "successionTimestep" = successionTimestep
    , "vegLeadingProportion" = vegLeadingProportion
    , ".plotInterval" = 1
    , ".plotMaps" = TRUE
    , ".saveInitialTime" = NA
    , ".useCache" = eventCaching # seems slower to use Cache for both
    , ".useParallel" = useParallel
  )
)

## Run Biomass_speciesData to get species layers
source("2_speciesLayers.R")

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
LBMR_testSim <- simInitAndSpades(times = simTimes
                                 , params = simParams
                                 , modules = simModules
                                 , objects = simObjects
                                 , paths = simPaths
                                 , debug = TRUE
                                 , .plotInitialTime = NA
)
endTime <- date()
cat(paste0("Took: ", endTime - startTime))
# End time: Wed Aug 28 17:16:08 2019

unlink(file.path(Paths$outputPath, "figures"), recursive = TRUE) ## remove unnecessary figures
