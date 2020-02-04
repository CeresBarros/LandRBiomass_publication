library(LandR)
library(raster)
library(quickPlot)

resDir <- "R/SpaDES/outputs/parametriseSALarge/"
outDir <- "C:/Ceres/temp"
dir.create(outDir)

pixelGroupMapFiles <- list.files(resDir, pattern = "pixelGroupMap")
cohortDataFiles <- list.files(resDir, pattern = "cohortData")

originalcrs <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
largeExtent <- extent(-104.757, -104.2197, 55.68663, 56.20319)
smallExtent <- largeExtent
smallExtent@xmax <- largeExtent@xmin + (largeExtent@xmax - largeExtent@xmin)*0.5   ## this is minimum size (10 000 pix at 250m res) -- need to increase large are to double
smallExtent@ymax <- largeExtent@ymin + (largeExtent@ymax - largeExtent@ymin)*0.5
studyAreaS <- as(smallExtent, "SpatialPolygons")
crs(studyAreaS) <- originalcrs
studyAreaS <- spTransform(studyAreaS, originalcrs)
Biomass_corecrs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
studyAreaS <- spTransform(studyAreaS, Biomass_corecrs)

years <- sub(".rds", "", sub(".*year", "year", pixelGroupMapFiles))

data("sppEquivalencies_CA", package = "LandR")
sppEquivalencies_CA[grep("Pin", LandR), `:=`(EN_generic_short = "Pine",
                                             EN_generic_full = "Pine",
                                             Leading = "Pine leading")]
sppEquivalencies_CA[grep("Betu_pap", LandR), `:=`(EN_generic_short = "Birch",
                                                  EN_generic_full = "Birch",
                                                  Leading = "Birch leading")]
## all Popu will be merged
sppEquivalencies_CA[grep("Popu_", LandR), `:=`(EN_generic_short = "Poplar",
                                               EN_generic_full = "Poplar",
                                               Leading = "Poplar leading")]

sppEquivalencies_CA[grep("Popu_", LandR), Boreal := "Popu_Spp"]

## define spp column to use for model
sppEquivCol <- "Boreal"
sppEquivalencies_CA <- sppEquivalencies_CA[!"", on = sppEquivCol]
sppEquivalencies_CA <- na.omit(sppEquivalencies_CA, sppEquivCol)

## create color palette for species used in model
sppColorVect <- sppColors(sppEquivalencies_CA, sppEquivCol,
                          newVals = "Mixed", palette = "Accent")

for(yr in years) {
  pixelGroupMap <- readRDS(file.path(resDir, paste0("pixelGroupMap_", yr, ".rds")))
  pixelGroupMap <- crop(pixelGroupMap, studyAreaS)
  cohortData <- readRDS(file.path(resDir, paste0("cohortData_", yr, ".rds")))
  # all(pixelGroupMap[!is.na(pixelGroupMap[])] %in% cohortData$pixelGroup)

  vegTypeMap <- vegTypeMapGenerator(cohortData, pixelGroupMap,
                                    0.6, mixedType = 2,
                                    sppEquiv = sppEquivalencies_CA, sppEquivCol = "Boreal",
                                    colors = sppColorVect,
                                    doAssertion = getOption("LandR.assertions", TRUE))
  levs <- raster::levels(vegTypeMap)[[1]]
  levelsName <- names(levs)[2]
  levsLeading <- equivalentName(levs[[levelsName]], sppEquivalencies_CA, "Leading")
  shortNames <- equivalentName(levsLeading, sppEquivalencies_CA, "EN_generic_short")
  whMixedLevs <- which(rownames(levs) == "Mixed")
  whMixedSppColors <- which(names(sppColorVect) == "Mixed")
  shortNames[whMixedLevs] <- "Mixed"
  levs[[levelsName]] <- shortNames
  levels(vegTypeMap) <- levs

  png(filename = file.path(outDir, paste0("vegTypeMap_", yr, ".png")))
  Plot(vegTypeMap)
  dev.off()
}
