## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list=ls()); amc::.gc()

## requires as of June 10th 2019
# loading reproducible     0.2.8.9001
# loading quickPlot        0.1.6.9000
# loading SpaDES.core      0.2.5.9004
# loading SpaDES.tools     0.3.2.9000
# loading SpaDES.addins    0.1.2

library(SpaDES)
library(LandR)
library(raster)

## define paths
setPaths(modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------
## get larger study area and create a smaller one (half extent)
studyAreaL <- shapefile("C:/Ceres/SKStudyAreaL")

largeExtent <- extent(studyAreaL)
smallExtent <- extent(studyAreaL)
smallExtent@xmax <- largeExtent@xmin + (largeExtent@xmax - largeExtent@xmin)/2
smallExtent@ymax <- largeExtent@ymin + (largeExtent@ymax - largeExtent@ymin)/2

studyAreaS <- as(smallExtent, "SpatialPolygons")
studyAreaS <-  SpatialPolygonsDataFrame(studyAreaS, data.frame(id = 1:length(studyAreaS)))
crs(studyAreaS) <- crs(studyAreaL)

studyAreaL <- spTransform(studyAreaL,
                         "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
studyAreaS <- spTransform(studyAreaS,
                          "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")


plot(studyAreaL); plot(studyAreaS, add = TRUE, col = "red")

## Set up sppEquiv  ------------------------------
data("sppEquivalencies_CA", package = "LandR")
sppEquivalencies_CA[grep("Pin", LandR), `:=`(EN_generic_short = "Pine",
                                             EN_generic_full = "Pine",
                                             Leading = "Pine leading")]

## define spp column to use for model
sppEquivCol <- "Boreal"
sppEquivalencies_CA <- na.omit(sppEquivalencies_CA, cols = sppEquivCol)

## create color palette for species used in model
sppColorVect <- sppColors(sppEquivalencies_CA, sppEquivCol,
                          newVals = "Mixed", palette = "Accent")


## Set up modelling parameters  ---------------------------
options('reproducible.useNewDigestAlgorithm' = TRUE)
runName <- "studyAreaL"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths
pathsSim <- getPaths()
pathsSim$outputPath <- file.path(pathsSim$outputPath, runName)
pathsSim$cachePath <- file.path("R/SpaDES/cache", runName)

## simulation params
timesSim <- list(start = 0, end = 10)
eventCaching <- c(".inputObjects", "init")


vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
successionTimestep <- 10L  # for dispersal and age reclass.

modulesSim <- list("BiomassSpeciesData"
                   , "Boreal_LBMRDataPrep"
                   , "LBMR"
)

objectsSim <- list("studyArea" = get(runName)
                   , "sppEquiv" = sppEquivalencies_CA
                   , "sppColorVect" = sppColorVect
)

paramsSim <- list(
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
    , ".plotInitialTime" = timesSim$start
    , "seedingAlgorithm" = "wardDispersal"
    , "sppEquivCol" = sppEquivCol
    , "successionTimestep" = successionTimestep
    , "vegLeadingProportion" = vegLeadingProportion
    , ".plotInterval" = 1
    , ".plotMaps" = FALSE
    , ".saveInitialTime" = NA
    , ".useCache" = eventCaching[eventCaching] # seems slower to use Cache for both
    , ".useParallel" = useParallel
  )
  , BiomassSpeciesData = list(
    "types" = c("KNN")
    , "sppEquivCol" = sppEquivCol
    , ".useCache" = TRUE
  )
)

options(spades.moduleCodeChecks = FALSE)
graphics.off()
LBMR_testSim <- simInitAndSpades(times = timesSim
                                 , params = paramsSim
                                 , modules = modulesSim
                                 , objects = objectsSim
                                 , paths = pathsSim
                                 , debug = TRUE
                                 , .plotInitialTime = NA
)


