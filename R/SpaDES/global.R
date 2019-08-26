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

## define paths
setPaths(modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------
## create a larger study area and create a smaller one (half extent)
## note that projection of the orignal CRS is always necessary
originalcrs <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
largeExtent <- extent(-104.757, -104.2197, 55.68663, 56.20319)
smallExtent <- largeExtent

smallExtent@xmax <- largeExtent@xmin + (largeExtent@xmax - largeExtent@xmin)*0.5   ## this is minimum size (10 000 pix at 250m res) -- need to increase large are to double
smallExtent@ymax <- largeExtent@ymin + (largeExtent@ymax - largeExtent@ymin)*0.5

studyAreaL <- as(largeExtent, "SpatialPolygons")
studyAreaL <-  SpatialPolygonsDataFrame(studyAreaL, data.frame(id = 1:length(studyAreaL)))
crs(studyAreaL) <- originalcrs
studyAreaL <- spTransform(studyAreaL, originalcrs)

studyAreaS <- as(smallExtent, "SpatialPolygons")
studyAreaS <-  SpatialPolygonsDataFrame(studyAreaS, data.frame(id = 1:length(studyAreaS)))
crs(studyAreaS) <- originalcrs
studyAreaS <- spTransform(studyAreaS, originalcrs)

## now reproject to LBMR standard
LBMRcrs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
studyAreaL <- spTransform(studyAreaL, LBMRcrs)
studyAreaS <- spTransform(studyAreaS, LBMRcrs)

plot(studyAreaL); plot(studyAreaS, add = TRUE, col = "red")

## Set up sppEquiv  ------------------------------
data("sppEquivalencies_CA", package = "LandR")
sppEquivalencies_CA[grep("Pin", LandR), `:=`(EN_generic_short = "Pine",
                                             EN_generic_full = "Pine",
                                             Leading = "Pine leading")]
# ## these species will be retained as separate
sppEquivalencies_CA[grep("Betu_pap", LandR), `:=`(EN_generic_short = "Birch",
                                             EN_generic_full = "Birch",
                                             Leading = "Birch leading")]
sppEquivalencies_CA[grep("Popu_bal", LandR), `:=`(EN_generic_short = "Poplar",
                                                  EN_generic_full = "Poplar",
                                                  Leading = "Poplar leading")]


## define spp column to use for model
sppEquivCol <- "Boreal"
sppEquivalencies_CA <- sppEquivalencies_CA[!"", on = sppEquivCol]
sppEquivalencies_CA <- na.omit(sppEquivalencies_CA, sppEquivCol)

## create color palette for species used in model
sppColorVect <- sppColors(sppEquivalencies_CA, sppEquivCol,
                          newVals = "Mixed", palette = "Accent")


## Set up modelling parameters  ---------------------------
options('reproducible.useNewDigestAlgorithm' = TRUE)
# runName <- "testAllModules"
# runName <- "testLBMRonly"
runName <- "studyAreaS"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths
pathsSim <- getPaths()
pathsSim$outputPath <- file.path(pathsSim$outputPath, runName)
pathsSim$cachePath <- file.path("R/SpaDES/cache", runName)

## simulation params
timesSim <- list(start = 0, end = 500)
eventCaching <- c(".inputObjects", "init")

vegLeadingProportion <- 0 # indicates what proportion the stand must be in one species group for it to be leading.
successionTimestep <- 10L  # for dispersal and age reclass.

modulesSim <- list("BiomassSpeciesData"
                   , "Boreal_LBMRDataPrep"
                   , "LBMR"
)

objectsSim <- list("studyArea" = if (grepl("study", runName)) get(runName) else studyAreaS
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
    , ".plotMaps" = TRUE
    , ".saveInitialTime" = NA
    , ".useCache" = eventCaching # seems slower to use Cache for both
    , ".useParallel" = useParallel
  )
  , BiomassSpeciesData = list(
    "types" = c("KNN")
    , "sppEquivCol" = sppEquivCol
    , ".useCache" = eventCaching
  )
)

options(spades.moduleCodeChecks = FALSE)
options("reproducible.useCache" = TRUE)
LBMR_testSim <- simInitAndSpades(times = timesSim
                                 , params = paramsSim
                                 , modules = modulesSim
                                 , objects = objectsSim
                                 , paths = pathsSim
                                 , debug = TRUE
                                 # , .plotInitialTime = NA
)
unlink(file.path(Paths$outputPath, "figures"), recursive = TRUE) ## remove unnecessary figures
