## clean workspace
rm(list=ls()); amc::.gc()

library(map)
library(data.table)
library(ggplot2)
library(raster)
library(quickPlot)
library(ggpubr)
library(gridExtra)
library(SpaDES)
library(LandR)
library(SpaDES.experiment)
library(purrr)
library(magick)


## Set up modelling parameters  ---------------------------
options('reproducible.useNewDigestAlgorithm' = TRUE)
# runName <- "studyAreaS"
# runName <- "studyAreaL"
runName <- "parametriseSALarge"

## paths
simPaths <- list(cachePath = file.path("R/SpaDES/cache", runName),
                 modulePath = file.path("R/SpaDES/m"),
                 inputPath = file.path("R/SpaDES/inputs"),
                 outputPath = file.path("R/SpaDES/outputs", runName))

## simLists
factorialSimulations <- readRDS(list.files(simPaths$outputPath, "simList_factorialSimulations",
                                           full.names = TRUE))
simListInit <- readRDS(list.files(simPaths$outputPath, paste0("simList_", runName), full.names = TRUE))

## get files of validation year per rep
outputFiles <- lapply(factorialSimulations, outputs)
outputFiles <- rbindlist(lapply(seq_along(outputFiles), FUN = function(x) {
  DT <- as.data.table(outputFiles[[x]])
  DT <- DT[saveTime %in% c(0, 10)]   ## 2001 is the first year (0), 2011, is year 10
  DT[, rep := x]
}), use.names = TRUE)

## load simulation outputs
allCohortData <- rbindlist(fill = TRUE, use.names = TRUE,
                           l = apply(outputFiles[objectName == "cohortData"], MARGIN = 1, FUN = function(x) {
                             cohortData <- readRDS(x["file"])
                             cohortData[, year := as.numeric(x["saveTime"])]
                             cohortData[, rep := as.numeric(x["rep"])]
                             return(cohortData)
                           }))

pixelGroupMapStk <- stack(apply(outputFiles[objectName == "pixelGroupMap"], MARGIN = 1, FUN = function(x) {
  pixelGroupMap <- readRDS(x["file"])
  names(pixelGroupMap) <- paste0("year", as.numeric(x["saveTime"]), "_rep", as.numeric(x["rep"]))
  pixelGroupMap
}))

vegTypeMapStk <- apply(outputFiles[objectName == "vegTypeMap"], MARGIN = 1, FUN = function(x) {
  vegTypeMap <- try(readRDS(x["file"]), silent = TRUE)
  if (class(vegTypeMap) != "try-error") {
    names(vegTypeMap) <- paste0("year", as.numeric(x["saveTime"]), "_rep", as.numeric(x["rep"]))
    vegTypeMap
  }
})

## cheat - TODO: need to figure out why some reps have no vegMap for year 0
vegTypeMapStk[sapply(vegTypeMapStk, is.null)] <- vegTypeMapStk[1]
vegTypeMapStk <- stack(vegTypeMapStk)
names(vegTypeMapStk)[names(vegTypeMapStk) %in% c("year0_rep1.1", "year0_rep1.2", "year0_rep1.3", "year0_rep1.4")] <- c("year0_rep1", "year0_rep3", "year0_rep6", "year0_rep7")

## load validation layers
rstDisturbedPix <- readRDS(list.files(simPaths$outputPath, "rstDisturbed", full.names = TRUE))
speciesLayersValidation <- readRDS(list.files(simPaths$outputPath, "speciesLayersValidation", full.names = TRUE))
rawBiomassMapValidation <- readRDS(list.files(simPaths$outputPath, "rawBiomassMapValidation", full.names = TRUE))
standAgeMapValidation <- readRDS(list.files(simPaths$outputPath, "standAgeMapValidation", full.names = TRUE))

speciesLayersInit <- simListInit$speciesLayers
rawBiomassMapInit <- simListInit$rawBiomassMap
standAgeMapInit <- simListInit$standAgeMap
## make a validation data tables, corrected for mismatches and with covers rescaled
## use makePixelTable and .createCohortData to rescale covers (to sum to 100)
## and filter bad data.

## year 2001
pixelTable <- makePixelTable(speciesLayers = speciesLayersInit,
                             biomassMap = rawBiomassMapInit,
                             standAgeMap = standAgeMapInit,
                             rasterToMatch = raster(file.path(simPaths$cachePath, "rasterToMatch.tif")))
pixelTable[, initialEcoregionCode := NULL]
pixelTable <- unique(pixelTable)
Bclass <- P(simListInit)$Biomass_borealDataPrep$pixelGroupBiomassClass
validationDataInit <- LandR:::.createCohortData(pixelTable, rescale = TRUE,
                                                pixelGroupBiomassClass = Bclass,
                                                doAssertion = FALSE)
validationDataInit[cover > 0 & age == 0, B := 0L]
validationDataInit[cover == 0 & age > 0, B := 0L]
validationDataInit[, totalBiomass := asInteger(sum(B)), by = "pixelIndex"]
validationDataInit[, relativeAbundValid := B/totalBiomass,
               by = .(pixelIndex, speciesCode)]
validationDataInit[, logAge := NULL]

## year 2011
pixelTable <- makePixelTable(speciesLayers = speciesLayersValidation,
                             biomassMap = rawBiomassMapValidation,
                             standAgeMap = standAgeMapValidation,
                             rasterToMatch = raster(file.path(simPaths$cachePath, "rasterToMatch.tif")))
pixelTable[, initialEcoregionCode := NULL]
pixelTable <- unique(pixelTable)
Bclass <- P(factorialSimulations$sim1_rep01)$Biomass_borealDataPrep$pixelGroupBiomassClass
validationData <- LandR:::.createCohortData(pixelTable, rescale = TRUE,
                                            pixelGroupBiomassClass = Bclass,
                                            doAssertion = FALSE)
validationData[cover > 0 & age == 0, B := 0L]
validationData[cover == 0 & age > 0, B := 0L]
validationData[, totalBiomass := asInteger(sum(B)), by = "pixelIndex"]
validationData[, relativeAbundValid := B/totalBiomass,
               by = .(pixelIndex, speciesCode)]
validationData[, logAge := NULL]

## expand cohortData
allCohortData <- rbindlist(fill = TRUE, use.names = TRUE,
                           l = lapply(as.list(pixelGroupMapStk), FUN = function(pixelGroupMap, allCohortData) {
                             yr <- as.numeric(sub("year", "", sub("_rep.*", "", names(pixelGroupMap))))
                             rp <- as.numeric(sub(".*_rep", "", names(pixelGroupMap)))
                             pixelCohortData <- allCohortData[year == yr & rep == rp]
                             pixelCohortData <- addPixels2CohortData(pixelCohortData, pixelGroupMap)
                             pixelCohortData
                           }, allCohortData = allCohortData))

## JOIN VALIDATION AND SIMULATED DATA
## before joining, summarize cohortData to stand totalB per species and
## and biomass-averaged stand age.
standCohortData <- allCohortData[, .(B, sum(B), age),
                                 by = .(pixelIndex, rep, year, speciesCode)]
standCohortData <- standCohortData[, standAge := sum(age * B, na.rm = TRUE) / sum(B, na.rm = TRUE),
                                   by = .(pixelIndex, rep, year)]

## drop unnecessary columns and remove separate cohorts
standCohortData[, B := V2] ## overwrite
standCohortData[, `:=`(V2 = NULL, age = NULL)]
standCohortData <- unique(standCohortData)


## add validation data to standCohortData
## note that some pixelIndex X spp combinations are lacking
## because the validation data has spp in some pixels that are not found in the simulation
## data, and vice-versa. To make sure that the validation pixel X spp combinations are added to each
## rep/year the validation dataset needs to be extended - other wise some
## times the validation data is only joined to some reps/years, making the validation
## averages "vary" across reps/years
combinationsInit <- as.data.table(expand.grid(list(speciesCode = unique(standCohortData$speciesCode),
                                                   pixelIndex = unique(validationDataInit$pixelIndex),
                                                   rep = 1:10,
                                                   year = 0)))

combinationsValid <- as.data.table(expand.grid(list(speciesCode = unique(standCohortData$speciesCode),
                                                    pixelIndex = unique(validationData$pixelIndex),
                                                    rep = 1:10,
                                                    year = 10)))
validationDataInit <- validationDataInit[combinationsInit,
                                         on = c("pixelIndex", "speciesCode")]
validationData <- validationData[combinationsValid,
                                 on = c("pixelIndex", "speciesCode")]
validationData <- rbindlist(list(validationDataInit, validationData),
                            use.names = TRUE)

## TODO: exclude pixels that are not simulated (?)
validationData <- validationData[pixelIndex %in% standCohortData$pixelIndex]

## clean up
rm(combinationsInit, combinationsValid, validationDataInit)

cols <- c("rep", "year", "pixelIndex", "speciesCode",
          "cover", "B", "totalBiomass", "relativeAbundValid")
standCohortData <- validationData[, ..cols][standCohortData,
                                   on = c("rep", "year", "pixelIndex", "speciesCode")]

setnames(standCohortData, c("cover", "B", "totalBiomass", "i.B"),
         c("coverValid", "BValid", "totalBValid", "B"))


## add vegType
vegTypeMapStk <- vegTypeMapStk[[names(pixelGroupMapStk)]]  ## make sure order is the same
vegTypeTable <- rbindlist(fill = TRUE, use.names = TRUE,
                          l = mapply(FUN = function(pixelGroupMap, vegTypeMap) {
                            yr <- as.numeric(sub("year", "", sub("_rep.*", "", names(pixelGroupMap))))
                            rp <- as.numeric(sub(".*_rep", "", names(pixelGroupMap)))
                            data.table(pixelIndex = 1:ncell(pixelGroupMap),
                                       vegType = getValues(vegTypeMap),
                                       year = yr, rep = rp)
                          }, pixelGroupMap = as.list(pixelGroupMapStk),
                          vegTypeMap = as.list(vegTypeMapStk), SIMPLIFY = FALSE))
vegTypeTable <- na.omit(vegTypeTable)
standCohortData <- vegTypeTable[standCohortData, on = c("rep", "year", "pixelIndex")]

vegTypeLabels <- as.data.table(levels(vegTypeMapStk[[1]]))
standCohortData <- vegTypeLabels[, .(ID, VALUE)][standCohortData, on = "ID==vegType"]
setnames(standCohortData, old = c("ID", "VALUE"),
         new = c("vegTypeID", "vegType"))

## remove disturbed pixels
disturbedPix <- which(!is.na(getValues(rstDisturbedPix)))
standCohortData <- standCohortData[!pixelIndex %in% disturbedPix]

## calculate some summary metrics
standCohortData[, `:=`(noSppPix = as.numeric(length(unique(speciesCode))),
                       standB = sum(B, na.rm = TRUE)),
              by = .(rep, year, pixelIndex)]
standCohortData[, relativeAbund := B/standB,
              by = .(rep, year, pixelIndex, speciesCode)]
standCohortData[, `:=`(landscapeB = sum(B, na.rm = TRUE),
                       landscapeBValid = sum(BValid, na.rm = TRUE)),
                by = .(rep, year)]


## plot labels
speciesLabels <- c("Abie_Bal" = "Fir", "Lari_Lar" = "Larch",
                   "Betu_Pap" = "Birch", "Pice_Gla" = "Wh. spruce",
                   "Pice_Mar" = "Bl. spruce", "Pinu_Ban" = "Jack pine",
                   "Popu_Bal" = "Poplar")
speciesColours <- levels(vegTypeMapStk[[1]])[[1]]$colors
names(speciesColours) <- levels(vegTypeMapStk[[1]])[[1]]$VALUE

## TODO: the number of pixels with B for each species varies
## between reps and years, and combinations os spp, pixel, are not
## the same between sim and valid
## also, I don't think these landscape-wide averages are working.
## because relative abundances won't sum to one, if calculated per pixel.
# standCohortData[, mean(relativeAbundValid, na.rm = TRUE), by = .(rep, year, speciesCode)] %>%
#   +     .[, mean(V1), by = speciesCode] %>%
#   +     .[,sum(V1)]

## calculate everything at landscape level

plot1 <- ggplot(data = plotData,
                aes(x = speciesCode, y = relativeAbund)) +
  stat_summary(fun.y = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(aes(y = relativeAbundValid),
               fun.y = "mean", geom = "point") +
  scale_x_discrete(labels = speciesLabels) +
  theme_pubr() +
  labs(title = "Species mean relative abundance", y = "species B / stand B")

plot2 <- ggplot(data = summaryAllCohortData,
                aes(x = Year, y = log(MortalityBySpecies), colour = speciesCode)) +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  geom_line(size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours,
                      labels = speciesLabels) +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot3 <- ggplot(data = pixelBurnCohortData,
                aes(x = Year, fill = vegType)) +
  geom_area(stat = "count", position = "fill") +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = speciesColours,
                    labels = speciesLabels) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  labs(title = "No. pixels per vegetation type", y = "g/m^2") +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot4 <- ggplot(data = summaryAllCohortData,
                aes(x = Year, y = AgeBySppWeighted, colour = speciesCode)) +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  geom_line(size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours,
                      labels = speciesLabels) +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot5 <- ggplot(data = summaryAllCohortData,
                aes(x = Year, y = noCohorts, colour = speciesCode)) +
  geom_vline(xintercept = 5, size = 1.5, linetype = "dashed", colour = "red") +
  geom_line(size = 1.5) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  scale_colour_manual(values = speciesColours,
                      labels = speciesLabels) +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

ggpubr::ggarrange(plot1, plot2, plot3, nrow = 3,
                  legend = "right", common.legend = TRUE)
ggsave(filename = "R/SpaDES/outputs/blogSep2019_noPM_PM_BMortVegType.tiff",
       width = 10, height = 15)

## GIFS ------------------------------------------------------
## make GIFs of vegetation maps
## individual pics first
speciesLabels <- c("Fir", "Larch", "En. spruce", "Wh. spruce",
                   "Bl. spruce", "Lo. pine", "Aspen","Douglas-fir")
names(speciesLabels) <- levels(vegTypeMapStk_noPM[[1]])[[1]]$ID
speciesColours <- levels(vegTypeMapStk_noPM[[1]])[[1]]$colors
names(speciesColours) <-  levels(vegTypeMapStk_noPM[[1]])[[1]]$ID

foothillsMask <- simList_noPM$rawBiomassMap
foothillsMask[!is.na(foothillsMask)] <- 1
foothillsMaskDF <- as.data.frame(as(foothillsMask, "SpatialPixelsDataFrame"))
names(foothillsMaskDF) <- c("value", "x", "y")

makePNGs <- function(id, rasterStack, filePrefix, gif.dir) {
  suppressWarnings(dir.create(gif.dir, recursive = TRUE))
  cat("Making ", id, "\n")
  rasterVis::gplot(rasterStack[[id]],
                   maxpixels = ncell(rasterStack[[id]])) +
    geom_tile(data = foothillsMaskDF,
              aes(x = x, y = y), fill = "grey95") +
    geom_tile(aes(fill = as.factor(value))) +
    scale_fill_manual(values = speciesColours,
                      labels = speciesLabels) +
    theme_void() + theme(legend.position = "none", text = element_text(size = 20)) +
    labs(title = sub("year", "Year ", names(rasterStack)[id])) +
    coord_equal()
  ggsave(file.path(gif.dir, paste0(filePrefix, id, ".png")),
         device = "png", width=5, height=10, dpi = 300, units = "in")
}

map_df(.x = 1:nlayers(vegTypeMapStk_noPM), .f = makePNGs,
       rasterStack = vegTypeMapStk_noPM,
       filePrefix = "vegTypeMapStk_noPM",
       gif.dir = "R/SpaDES/outputs/blogSep2019_noPM_oneFire/gif")
map_df(.x = 1:nlayers(vegTypeMapStk_PM), .f = makePNGs,
       rasterStack = vegTypeMapStk_PM,
       filePrefix = "vegTypeMapStk_PM",
       gif.dir = "R/SpaDES/outputs/blogSep2019_PM_oneFire/gif")

makeGIF <- function(gif.dir, gifPrefix, ...) {
  ## get file list and the file numbers (to sort numerically, rather than alphabetically)
  PNGlist <- list.files(path = gif.dir, pattern = "*.png")
  fileNos <- as.numeric(sub("\\.png", "", sub("^\\D*(\\d)", "\\1", PNGlist)))

  ## make GIF
  file.path(gif.dir, PNGlist[order(fileNos)]) %>%
    map(image_read) %>% # reads each path file
    image_join() %>% # joins image
    image_animate(...) %>% # animates, can opt for number of loops
    image_write(file.path(gif.dir, paste0(gifPrefix, ".gif")))
}

makeGIF(gif.dir = "R/SpaDES/outputs/blogSep2019_noPM_oneFire/gif",
        gifPrefix = "vegTypeMapStk_noPM",
        fps = 2)
makeGIF(gif.dir = "R/SpaDES/outputs/blogSep2019_PM_oneFire/gif",
        gifPrefix = "vegTypeMapStk_PM",
        fps = 2)

## OTHER PLOTS --------------------------------------------------------
## topo and climate examples
slopeRas <- projectRaster(simList_noPM$slopeRas, foothillsMask)
topoPal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "BrBG"))
plot(foothillsMask, col = "grey90", axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)
plot(slopeRas, col = topoPal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE, add = TRUE)

temperatureRas <- projectRaster(simList_noPM$temperatureRas, foothillsMask)
tempPal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "RdBu"))
plot(foothillsMask, col = "grey90", axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)
plot(temperatureRas, col = tempPal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE, add = TRUE)

precipitationRas <- projectRaster(simList_noPM$precipitationRas, foothillsMask)
plot(foothillsMask, col = "grey90", axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)
plot(precipitationRas, col = tempPal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE, add = TRUE)

## spp cover, age and biomass examples
sppPal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Blues"))
plot(foothillsMask, col = "grey90", axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)
plot(simList_noPM$speciesLayers[[1]], col = sppPal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE, add = TRUE)

plot(foothillsMask, col = "grey90", axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)
plot(simList_noPM$rawBiomassMap, axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE, add = TRUE)

plot(foothillsMask, col = "grey90", axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)
agePal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Greens"))
plot(simList_noPM$standAgeMap, col = agePal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE, add = TRUE)

## ecodistricts
## https://www.statcan.gc.ca/eng/subjects/standard/environment/elc/12-607-x2018001-eng.pdf

ecoDistSF <- sf::st_as_sf(simList_noPM$ecoDistrict)
ecoDistSF$ECODISTRIC <- factor(ecoDistSF$ECODISTRIC,
                               levels = c(798, 800, 801, 799, 793, 750, 631, 1018, 1017, 1019))
canada <- sf::st_as_sf(shapefile("data/CA_admin/gpr_000a11a_e.shp"))
canada <- sf::st_transform(canada, crs = crs(ecoDistSF))
alberta <- canada[canada$PRENAME %in% "Alberta",]

ggplot(ecoDistSF) +
  geom_sf(data = alberta) +
  geom_sf(data = ecoDistSF, aes(fill = as.factor(ECODISTRIC))) +
  scale_fill_brewer(palette = "Paired",
                    labels = c("631" = "W AB upland - Foothills",
                               "750" = "Aspen parkland - Upland",
                               "793" = "Moist mixed grassland - Plain",
                               "798" = "Fescue grassland - Plain",
                               "799" = "Fescue grassland - Upland",
                               "800" = "Fescue grassland - Plain",
                               "801" = "Fescue grassland - Foothills",
                               "1017" = "N Cont. Divide - Mountains",
                               "1018" = "N Cont. Divide - Foothills",
                               "1019" = "N Cont. Divide - Mountains")) +
  theme_void() +
  theme(text = element_text(colour = "white"),
        plot.background = element_rect(fill = "black")) +
  labs(fill = "Ecoregion - ecodistrict") +
  coord_sf()

plot(sf::st_as_sf(simList_noPM$ecoDistrict["ECODISTRIC"]))


save(list = grep("model", ls(), value = TRUE), file = "E:/GitHub/LandscapesInMotion/analyses/modelsGAMLSS_0-3Days_goodSample_Nov1_v2.RData")























## list outputs of one repetition
outputs(factorialSimulations$sim1_rep01)
# Try loading one object from output(sim)
# b <- readRDS(outputs(factorialSimulations$sim1_rep01)$file[1])

# Now run the code snippet to make the figure. This was taken directly from ?experiment2
fn <- quote({
  pixelCohortData <- addNoPixel2CohortData(cohortData, pixelGroupMap)
  pixelCohortData[, list(BiomassBySpecies = as.numeric(sum(B * noPixels, na.rm = TRUE)),
                         MortalityBySpecies = as.numeric(sum(mortality * noPixels, na.rm = TRUE)),
                         aNPPBySpecies = as.numeric(sum(aNPPAct * noPixels, na.rm = TRUE)),
                         AgeBySppWeighted = as.numeric(sum(age * B * noPixels, na.rm = TRUE) /
                                                         sum(B * noPixels, na.rm = TRUE)),
                         noCohorts = as.numeric(length(unique(age)))),
                  by = speciesCode] %>%
    as.list(.)
})
b <- as.data.table(factorialSimulations,
                   vals = fn,
                   objectsFromOutputs = lapply(1:5, function(x) c("cohortData", "pixelGroupMap")))
p <- ggplot(b, aes(x=simList, y=NBurnedPixels, group=simList, color=simList)) +
  stat_summary(geom = "point", fun.y = mean) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.2)
print(p)
