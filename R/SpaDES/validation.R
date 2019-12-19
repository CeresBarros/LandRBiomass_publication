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
                                       year = yr,
                                       rep = rp)
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

## convert NAs to zeros in biomasses
standCohortData[is.na(B), B := 0]
standCohortData[is.na(BValid), BValid := 0]

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

## LANDSCAPE-WIDE COMPARISONS IN A GIVEN YEAR --------------------
## relative abundance (biomass) per species
plotData <- standCohortData[year == 10,
                            list(landRelativeAbund = sum(B, na.rm = TRUE)/unique(landscapeB),
                                 landRelativeAbundValid = sum(BValid, na.rm = TRUE)/unique(landscapeBValid)),
                            by = .(rep, speciesCode)]
plot1 <- ggplot(data = plotData,
                aes(x = speciesCode, y = landRelativeAbund)) +
  stat_summary(fun.y = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(aes(y = landRelativeAbundValid),
               fun.y = "mean", geom = "point") +
  scale_x_discrete(labels = speciesLabels) +
  theme_pubr() +
  labs(title = "Landscape-wide species relative abundances",
       x = "", y = expression(over("species B", "total B")))

## no. pixels with a species
plotData <- standCohortData[year == 10, list(count = sum(B > 0),
                                             countValid = sum(BValid > 0)),
                            by = .(rep, speciesCode)]
plotData <- melt(plotData, measure.vars = c("count", "countValid"),
                 variable.name = "dataType", value.name = "count")

## HERE: TODO: THE NUMBER OF PIXELS IN VALID IS VARYING BETWEENS REPS FOR THE SAME SPECIES.
plot2 <- ggplot(data = plotData[dataType == "count"],
                aes(x = speciesCode, y = count)) +
  stat_summary(fun.y = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = plotData[dataType == "countValid"],
               aes(x = speciesCode, y = count),
               fun.y = "mean", geom = "point") +
  scale_x_discrete(labels = speciesLabels) +
  theme_pubr() +
  labs(title = "Species presences",
       x = "", y = "no. pixels",
       fill = "")

plot3 <- ggplot(data = plotData, aes(x = speciesCode, y = count, fill = dataType)) +
 geom_boxplot() +
  scale_x_discrete(labels = speciesLabels) +
  scale_fill_discrete(labels = c("count" = "simulated",
                                 "countValid" = "observed")) +
  theme_pubr() +
  labs(title = "Species presences",
       x = "", y = "no. pixels",
       fill = "")

## relative abundance per dominant species
## get dominant species - these match with model outputs, I checked
plotData <- unique(standCohortData[year == 10,
                            list(vegTypeID2 = speciesCode[which.max(B)],
                                 B = max(B, na.rm = TRUE),
                                 landscapeB = unique(landscapeB),
                                 vegTypeID2Valid = speciesCode[which.max(BValid)],
                                 BValid = max(BValid, na.rm = TRUE),
                                 landscapeBValid = unique(landscapeBValid)),
                            by = .(rep, pixelIndex)])

## calculate relative abundances and melt data.
plotData[, landRelativeAbund := sum(B, na.rm = TRUE)/unique(landscapeB),
         by = .(rep, vegTypeID2)]
plotData[, landRelativeAbundValid := sum(BValid, na.rm = TRUE)/unique(landscapeBValid),
         by = .(rep, vegTypeID2Valid)]
plotData1 <- plotData[, .(rep, pixelIndex, vegTypeID2, B, landscapeB, landRelativeAbund)]
plotData1[, dataType := "simulated"]
plotData2 <- plotData[, .(rep, pixelIndex, vegTypeID2Valid, BValid,
                          landscapeBValid, landRelativeAbundValid)]
plotData2[, dataType := "observed"]
setnames(plotData2, grep("Valid", names(plotData2), value = TRUE),
         sub("Valid", "", grep("Valid", names(plotData2), value = TRUE)))

plotData <- rbind(plotData1, plotData2)
rm(plotData1, plotData2)

plot2 <- ggplot(data = plotData[dataType == "simulated"],
                aes(x = vegTypeID2, y = landRelativeAbund)) +
  stat_summary(fun.y = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = plotData[dataType != "simulated"],
               aes(x = vegTypeID2, y = landRelativeAbund),
               fun.y = "mean", geom = "point") +
  scale_x_discrete(labels = speciesLabels) +
  theme_pubr() +
  labs(title = "Landscape-wide relative abundances of dominant species",
       x = "", y = expression(over("species B", "total B")))

## no. pixels with a certain dominant species
plot3 <- ggplot(data = plotData, aes(x = vegTypeID2, fill = dataType)) +
  geom_histogram(stat = "count", position = "dodge") +
  scale_x_discrete(labels = speciesLabels) +
  theme_pubr() +
  labs(title = "Landscape-wide species dominance",
       x = "", y = "no. pixels",
       fill = "")

## as previous, but in relative terms
plot4 <- ggplot(data = plotData, aes(x = dataType, fill = vegTypeID2)) +
  geom_histogram(stat = "count", position = "fill") +
  scale_fill_brewer(palette = "Dark2", labels = speciesLabels) +
  theme_pubr() +
  labs(title = "Landscape-wide species dominance", x = "",
       y = "proportion of pixels", fill = "dominant species")

grid.arrange(plot1, plot2, plot3, plot4,
             nrow = 2, ncol = 2)

## COMPARISONS OF DELTA PER PIXEL-------------------
## per species
plotData <- standCohortData[year %in% c(0,10),]
plotData <- plotData[, list(deltaB = B[which(year == 10)] - B[which(year == 0)],
                            deltaBValid = BValid[which(year == 10)] - BValid[which(year == 0)]),
                     , by = .(rep, pixelIndex, speciesCode)]
plotData <- melt(plotData, measure.vars = c("deltaB", "deltaBValid"),
                 variable.name = "dataType", value.name = "deltaB")

plot5 <- ggplot(data = plotData,
                aes(x = speciesCode, y = deltaB, fill = dataType)) +
  geom_boxplot() +
  scale_x_discrete(labels = speciesLabels) +
  scale_fill_discrete(labels = c("deltaB" = "simulated",
                                 "deltaBValid" = "observed")) +
  theme_pubr() +
  labs(title = "Pixel changes in abundances", fill = "",
       x = "", y = expression(paste(Delta, "B")))


plotData <- plotData[, mean(deltaB) , by = .(rep, speciesCode, dataType)]

plot6 <- ggplot(data = plotData[dataType == "deltaB"],
                aes(x = speciesCode, y = V1)) +
  stat_summary(fun.y = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = plotData[dataType == "deltaBValid"],
               aes(x = speciesCode, y = V1),
               fun.y = "mean", geom = "point") +
  scale_x_discrete(labels = speciesLabels) +
  theme_pubr() +
  labs(title = "Landscape-wide changes in abundances",
       x = "", y = expression(paste(Delta, "B")))

grid.arrange(plot5, plot6,
             ncol = 2)
