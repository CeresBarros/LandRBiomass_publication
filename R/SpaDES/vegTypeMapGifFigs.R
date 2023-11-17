pkgPath <- file.path("packages", version$platform,
                     paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1]))
.libPaths(pkgPath, include.site = FALSE)

out <- Require::Require(c("LandR", "raster", "quickPlot", "purrr", "magick", "ggplot2"),
                 upgrade = FALSE, install = FALSE, standAlone = TRUE)

resDir <- "R/SpaDES/outputs/oct2022Runs/baseCase/sim1_rep01/"
outDir <- "~/temp"
dir.create(outDir)

pixelGroupMapFiles <- list.files(resDir, pattern = "pixelGroupMap")

originalcrs <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
largeExtent <- raster::extent(-104.757, -104.2197, 55.68663, 56.20319)
smallExtent <- largeExtent
smallExtent@xmax <- largeExtent@xmin + (largeExtent@xmax - largeExtent@xmin)*0.5   ## this is minimum size (10 000 pix at 250m res) -- need to increase large are to double
smallExtent@ymax <- largeExtent@ymin + (largeExtent@ymax - largeExtent@ymin)*0.5
studyAreaS <- as(smallExtent, "SpatialPolygons")
crs(studyAreaS) <- originalcrs
studyAreaS <- spTransform(studyAreaS, originalcrs)
Biomass_corecrs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
studyAreaS <- spTransform(studyAreaS, Biomass_corecrs)

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

## veg type map
for (pgmFile in pixelGroupMapFiles) {
  yr <- sub(".rds", "", sub(".*year", "year", pgmFile))
  pixelGroupMap <- readRDS(file.path(resDir, pgmFile))
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

  plot1 <- ggplot(as.data.frame(vegTypeMap, xy = TRUE), aes(x = x, y = y, fill = layer_VALUE)) +
    geom_raster() +
    scale_fill_manual(values = na.omit(unique(as.data.frame(vegTypeMap, xy = TRUE)$layer_colors)),
                      na.value = "transparent",
                      breaks = na.omit(unique(as.data.frame(vegTypeMap, xy = TRUE)$layer_VALUE))) +
    theme_bw() +
    coord_equal() +
    labs(title = paste("Year", sub("year", "", yr)), x = "Longitude", y = "Latitude", fill = "Dominant species")
  ggsave(plot = plot1, filename = file.path(outDir, paste0("vegTypeMap_", yr, ".png")), dpi = 300)
}

## stand biomassMap
for (pgmFile in pixelGroupMapFiles) {
  yr <- sub(".rds", "", sub(".*year", "year", pgmFile))
  pixelGroupMap <- readRDS(file.path(resDir, pgmFile))
  pixelGroupMap <- crop(pixelGroupMap, studyAreaS)
  cohortData <- readRDS(file.path(resDir, paste0("cohortData_", yr, ".rds")))
  # all(pixelGroupMap[!is.na(pixelGroupMap[])] %in% cohortData$pixelGroup)

  pcd <- addPixels2CohortData(cohortData, pixelGroupMap)
  pcd[, totalBiomass := sum(B), by = pixelIndex]
  biomassMap <- LandR::makeBiomassMap(pcd, pixelGroupMap)

  plot1 <- ggplot(as.data.frame(biomassMap, xy = TRUE), aes(x = x, y = y, fill = layer)) +
    geom_raster() +
    scale_fill_distiller(palette = "Greens", direction = 1, na.value = "transparent", limits = c(0, 13000)) +
    theme_bw() +
    coord_equal() +
    labs(title = paste("Year", sub("year", "", yr)), x = "Longitude", y = "Latitude", fill = "g/m2")
  ggsave(filename = file.path(outDir, paste0("standBiomassMap_", yr, ".png")), dpi = 300)
}



makeGIF <- function(gif.dir, gifPrefix, pngPrefix, ...) {
  ## get file list and the file numbers (to sort numerically, rather than alphabetically)
  PNGlist <- list.files(path = gif.dir, pattern = pngPrefix)
  fileNos <- as.numeric(sub("\\.png", "", sub("^\\D*(\\d)", "\\1", PNGlist)))

  ## make GIF
  file.path(gif.dir, PNGlist[order(fileNos)]) |>
    map(image_read) |> # reads each path file
    image_join() |> # joins image
    image_animate(...) |> # animates, can opt for number of loops
    image_write(file.path(gif.dir, paste0(gifPrefix, ".gif")), quality = 100)
}

makeGIF(gif.dir = outDir,
        gifPrefix = "vegTypeMapChange",
        pngPrefix = "vegTypeMap_year",
        fps = 5)

makeGIF(gif.dir = outDir,
        gifPrefix = "standBMapChange",
        pngPrefix = "standBiomassMap_year",
        fps = 10)

