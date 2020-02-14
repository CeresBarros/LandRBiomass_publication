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
                 inputPath = file.path("data/"),
                 outputPath = file.path("R/SpaDES/outputs", runName))

## define local and remote repo paths for file name subs:
localProjPath <- getwd()
remoteProjPath <- "/mnt/storage/cbarros/LandRBiomass_publication"

## simLists
factorialSimulations <- readRDS(list.files(simPaths$outputPath, "simList_factorialSimulations",
                                           full.names = TRUE))
simListInit <- readRDS(list.files(simPaths$outputPath, paste0("simList_", runName), full.names = TRUE))

## get files of validation year per rep
## 2001 is the first year (0) and 2011, is year 10,
## but choose years 1 and 11 because objects were saved at the beginning of the year, they reflect the end of the previous year
## cohortData is the only that will vary across reps for year 0, because the beggining of year 1 already includes growth/mortality from year 0.
outputFiles <- lapply(factorialSimulations, outputs)
outputFiles <- rbindlist(lapply(seq_along(outputFiles), FUN = function(x) {
  DT <- as.data.table(outputFiles[[x]])
  DT <- DT[saveTime %in% c(1, 11)]
  DT[, rep := x]
}), use.names = TRUE)

## change directory if not on proj directory
if (all(!grepl(getwd(), outputFiles$file)))
  outputFiles$file <- sub("/mnt/storage/cbarros/LandRBiomass_publication",
                          getwd(), outputFiles$file)

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

vegTypeMapStk <- stack(apply(outputFiles[objectName == "vegTypeMap"], MARGIN = 1, FUN = function(x) {
  vegTypeMap <- try(readRDS(x["file"]), silent = TRUE)
  if (class(vegTypeMap) != "try-error") {
    names(vegTypeMap) <- paste0("year", as.numeric(x["saveTime"]), "_rep", as.numeric(x["rep"]))
    vegTypeMap
  }
}))

## load validation layers - they are saved in the outputs once for each rep.
## here the rstDisturbed only has 1s
speciesLayersInit <- if (inMemory(simListInit$speciesLayers))
  simListInit$speciesLayers else {
    sppCodes <- as.character(unique(allCohortData$speciesCode))
    sppCodesStr <- paste(sppCodes, collapse = "|")

    speciesLayersInit <- lapply(simListInit$speciesLayers@layers, function(x) {
      localFileName <- sub(remoteProjPath, localProjPath, x@file@name)
      sppLayer <- raster(localFileName)

      names(sppLayer) <- sub(paste0("(.*)(", sppCodesStr, ")(.*)"), "\\2",
                             basename(localFileName))
      return(sppLayer)
    }
    )
    rm(sppCodes, sppCodesStr)
    stack(speciesLayersInit)
  }

standAgeMapInit <- if (inMemory(simListInit$standAgeMap))
  simListInit$standAgeMap else {
    localFileName <- sub(remoteProjPath, localProjPath,
                         simListInit$standAgeMap@file@name)
    raster(localFileName)
  }

rawBiomassMapInit <- if (inMemory(simListInit$rawBiomassMap))
  simListInit$rawBiomassMap else {
    localFileName <- sub(remoteProjPath, localProjPath,
                         simListInit$rawBiomassMap@file@name)
    raster(localFileName)
  }

rstDisturbedPix <- if (inMemory(simListInit$rstLCChange))
  simListInit$rstLCChange else {
    localFileName <- sub(remoteProjPath, localProjPath,
                         simListInit$rstLCChange@file@name)
    raster(localFileName)
  }

speciesLayersValidation <- if (inMemory(simListInit$speciesLayersValidation))
  simListInit$speciesLayersValidation else {
    sppCodes <- as.character(unique(allCohortData$speciesCode))
    sppCodesStr <- paste(sppCodes, collapse = "|")

    speciesLayersValidation <- lapply(simListInit$speciesLayersValidation@layers, function(x) {
      localFileName <- sub(remoteProjPath, localProjPath, x@file@name)
      sppLayer <- raster(localFileName)

      names(sppLayer) <- sub(paste0("(.*)(", sppCodesStr, ")(.*)"), "\\2",
                             basename(localFileName))
      return(sppLayer)
    }
    )
    rm(sppCodes, sppCodesStr)
    stack(speciesLayersValidation)
  }

standAgeMapValidation <- if (inMemory(simListInit$standAgeMapValidation))
  simListInit$standAgeMapValidation else {
    localFileName <- sub(remoteProjPath, localProjPath,
                         simListInit$standAgeMapValidation@file@name)
    raster(localFileName)
  }

rawBiomassMapValidation <- if (inMemory(simListInit$rawBiomassMapValidation))
  simListInit$rawBiomassMapValidation else {
    localFileName <- sub(remoteProjPath, localProjPath,
                         simListInit$rawBiomassMapValidation@file@name)
    raster(localFileName)
  }

## check names of species layers match
if (!all(names(speciesLayersInit) %in% names(speciesLayersValidation)))
  stop("validation species layers names do not match!")


## year 2001
pixelTable <- makePixelTable(speciesLayers = speciesLayersInit,
                             biomassMap = rawBiomassMapInit,
                             standAgeMap = standAgeMapInit,
                             rasterToMatch = simListInit$rasterToMatch)
pixelTable[, initialEcoregionCode := NULL]
pixelTable <- unique(pixelTable)
Bclass <- P(simListInit)$Biomass_borealDataPrep$pixelGroupBiomassClass
validationDataInit <- LandR:::.createCohortData(pixelTable, rescale = TRUE,
                                                pixelGroupBiomassClass = Bclass,
                                                doAssertion = FALSE)
validationDataInit[cover > 0 & age == 0, B := 0L]
validationDataInit[, totalBiomass := asInteger(sum(B)), by = "pixelIndex"]
validationDataInit[, relativeAbundValid := B/totalBiomass,
                   by = .(pixelIndex, speciesCode)]
validationDataInit[, logAge := NULL]

## year 2011
pixelTable <- makePixelTable(speciesLayers = speciesLayersValidation,
                             biomassMap = rawBiomassMapValidation,
                             standAgeMap = standAgeMapValidation,
                             rasterToMatch = simListInit$rasterToMatch)
pixelTable[, initialEcoregionCode := NULL]
pixelTable <- unique(pixelTable)
Bclass <- P(factorialSimulations$sim1_rep1)$Biomass_borealDataPrep$pixelGroupBiomassClass
validationData <- LandR:::.createCohortData(pixelTable, rescale = TRUE,
                                            pixelGroupBiomassClass = Bclass,
                                            doAssertion = FALSE)
validationData[cover > 0 & age == 0, B := 0L]
validationData[, totalBiomass := asInteger(sum(B)), by = "pixelIndex"]
validationData[, relativeAbundValid := B/totalBiomass,
               by = .(pixelIndex, speciesCode)]
validationData[is.nan(relativeAbundValid), relativeAbundValid := 0]
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
## rep/year the validation dataset needs to be extended - otherwise some
## times the validation data is only joined to some reps/years, making the validation
## averages "vary" across reps/years
combinationsInit <- as.data.table(expand.grid(list(speciesCode = unique(standCohortData$speciesCode),
                                                   pixelIndex = unique(c(validationDataInit$pixelIndex, standCohortData$pixelIndex)),
                                                   rep = 1:10,
                                                   year = 1)))

combinationsValid <- as.data.table(expand.grid(list(speciesCode = unique(standCohortData$speciesCode),
                                                    pixelIndex = unique(c(validationData$pixelIndex, standCohortData$pixelIndex)),
                                                    rep = 1:10,
                                                    year = 11)))
validationDataInit <- validationDataInit[combinationsInit,
                                         on = c("pixelIndex", "speciesCode")]
validationData <- validationData[combinationsValid,
                                 on = c("pixelIndex", "speciesCode")]
validationData <- rbindlist(list(validationDataInit, validationData),
                            use.names = TRUE)

# ## exclude pixels that are not simulated
validationData <- validationData[pixelIndex %in% standCohortData$pixelIndex]

## clean up
rm(combinationsInit, combinationsValid, validationDataInit)

cols <- c("rep", "year", "pixelIndex", "speciesCode",
          "cover", "B", "totalBiomass", "relativeAbundValid")
standCohortData <- standCohortData[validationData[, ..cols],
                                   on = c("rep", "year", "pixelIndex", "speciesCode")]
setnames(standCohortData, c("cover", "i.B", "totalBiomass"),
         c("coverValid", "BValid", "totalBValid"))

## some NAs added on simutated data for year 2011
standCohortData[is.na(B), B := 0]
standCohortData[is.na(BValid), BValid := 0]  ## just to be sure
standCohortData[, standAge := max(standAge, na.rm = TRUE), by = .(year, rep, pixelIndex)]

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
if (!compareRaster(rstDisturbedPix, pixelGroupMapStk[[1]])) ## just check rasters match first
  stop("rstDisturbedPix and pixelGroupMap differ!")
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

## TODO: the number of pixels per dominant spp differs between simulated and observed data for year 1..!
## I don't think these landscape-wide averages are working.
## because relative abundances won't sum to one, if calculated per pixel.
# standCohortData[, mean(relativeAbundValid, na.rm = TRUE), by = .(rep, year, speciesCode)] %>%
#   .[, mean(V1), by = speciesCode] %>%
#   .[,sum(V1)]

## calculate everything at landscape level

## MAPS OF OBSERVED CHANGES IN STAND BIOMASS AND STAND AGE -------
standDeltaBValid <- (rawBiomassMapValidation - rawBiomassMapInit) * 100
standDeltaBValid[is.na(pixelGroupMapStk[[1]])] <- NA
standDeltaBValid[disturbedPix] <- NA
names(standDeltaBValid) <- "stand delta-B"

standDeltaAgeValid <- standAgeMapValidation - standAgeMapInit
standDeltaAgeValid[is.na(pixelGroupMapStk[[1]])] <- NA
standDeltaAgeValid[disturbedPix] <- NA
names(standDeltaAgeValid) <- "stand delta-Age"

Plot(standDeltaBValid, standDeltaAgeValid)

## what is the relationship between the two?
standDeltaValidData <- na.omit(data.table(pixelIndex = 1:ncell(standDeltaBValid),
                                          standDeltaBValid = getValues(standDeltaBValid),
                                          standDeltaAgeValid = getValues(standDeltaAgeValid)))

## try excluding all disturbances before 2001:
# distYears <- raster("R/SpaDES/m/Biomass_validationKNN/data/SmallC2C_change_year_1985_2011.tif")
# yrs <- seq(1980, 2011) - 1900
# pixKeep <- which(getValues(distYears) %in% yrs)
# standDeltaValidData <- standDeltaValidData[pixelIndex %in% pixKeep]

ggplot(standDeltaValidData,
       aes(x = standDeltaAgeValid, y = standDeltaBValid)) +
  geom_point() +
  stat_smooth(method = "lm")

summary(lm(standDeltaBValid ~ standDeltaAgeValid, data = standDeltaValidData))

ggplot(standDeltaValidData,
       aes(x = standDeltaBValid)) +
  geom_histogram()
ggplot(standDeltaValidData,
       aes(y = standDeltaBValid)) +
  geom_boxplot()

ggplot(standDeltaValidData,
       aes(x = standDeltaAgeValid)) +
  geom_histogram()
ggplot(standDeltaValidData,
       aes(y = standDeltaAgeValid)) +
  geom_boxplot()
ggplot(standDeltaValidData[standDeltaAgeValid == 10],
       aes(x = standDeltaBValid)) +
  geom_histogram()

## remove pixels were the aging was too extreme (lower/higher than 25% and 75% quantiles)
# quants <- quantile(standDeltaValidData$standDeltaAgeValid, probs = c(0.25, 0.75))
# extremeDeltaAgePix <- standDeltaValidData[standDeltaValidData$standDeltaAgeValid <= quants["25%"] |
#                                             standDeltaValidData$standDeltaAgeValid >= quants["75%"],
#                                           pixelIndex]
# extremeDeltaAgePix <- standDeltaValidData[standDeltaValidData$standDeltaAgeValid != 10,
#                                           pixelIndex]
#
# standCohortData <- standCohortData[!pixelIndex %in% extremeDeltaAgePix]

standCohortData <- standCohortData[pixelIndex %in% standDeltaValidData[standDeltaBValid > 0 & standDeltaAgeValid > 0,
                                                                       pixelIndex]]

## LANDSCAPE-WIDE COMPARISONS IN A GIVEN YEAR --------------------
## relative abundance (biomass) per species
plotData <- standCohortData[year %in% c(1, 11),
                            list(landRelativeAbund = sum(B, na.rm = TRUE)/unique(landscapeB),
                                 landRelativeAbundValid = sum(BValid, na.rm = TRUE)/unique(landscapeBValid)),
                            by = .(rep, year, speciesCode)]
plot1 <- ggplot(data = plotData,
                aes(x = speciesCode, y = landRelativeAbund)) +
  stat_summary(fun.y = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(aes(y = landRelativeAbundValid),
               fun.y = "mean", geom = "point", size = 2) +
  scale_x_discrete(labels = speciesLabels) +
  theme_pubr() +
  facet_wrap(~ year, labeller = labeller(year = c("1" = "2001", "11" = "2011"))) +
  labs(title = "Landscape-wide species relative abundances",
       x = "", y = expression(over("species B", "total B")))

## pixel-level species'realtive abundances
plotData <- standCohortData[year %in% c(1, 11),
                            list(landRelativeAbund = sum(B, na.rm = TRUE)/unique(landscapeB),
                                 landRelativeAbundValid = sum(BValid, na.rm = TRUE)/unique(landscapeBValid)),
                            by = .(rep, year, pixelIndex, speciesCode)]
plotData <- melt(plotData,
                 id.vars = c("rep", "year", "pixelIndex", "speciesCode"))
plot1.2 <- ggplot(data = plotData,
                aes(x = speciesCode, y = value, fill = variable)) +
  geom_boxplot() +
  scale_x_discrete(labels = speciesLabels) +
  scale_fill_discrete(labels = c("landRelativeAbund" = "simulated",
                                 "landRelativeAbundValid" = "observed")) +
  theme_pubr() +
  facet_wrap(~ year, labeller = labeller(year = c("1" = "2001", "11" = "2011"))) +
  labs(title = "Pixel-level species relative abundances",
       x = "", y = expression(over("species B", "total B")))

## no. pixels with a species
plotData <- standCohortData[year %in% c(1, 11),
                            list(count = sum(B > 0),
                                 countValid = sum(BValid > 0)),
                            by = .(rep, year, speciesCode)]
plotData <- melt(plotData, measure.vars = c("count", "countValid"),
                 variable.name = "dataType", value.name = "count")

plot2 <- ggplot(data = plotData[dataType == "count"],
                aes(x = speciesCode, y = count)) +
  stat_summary(fun.y = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = plotData[dataType == "countValid"],
               aes(x = speciesCode, y = count),
               fun.y = "mean", geom = "point", size = 2) +
  scale_x_discrete(labels = speciesLabels) +
  theme_pubr() +
  facet_wrap(~ year, labeller = labeller(year = c("1" = "2001", "11" = "2011"))) +
  labs(title = "Species presences",
       x = "", y = "no. pixels",
       fill = "")

## relative abundance per dominant species
## get dominant species (those with maxB) - these match with model outputs, I checked
plotData <- unique(standCohortData[year %in% c(1, 11),
                                   list(vegTypeID2 = speciesCode[which.max(B)],
                                        B = max(B, na.rm = TRUE),
                                        landscapeB = unique(landscapeB),
                                        vegTypeID2Valid = speciesCode[which.max(BValid)],
                                        BValid = max(BValid, na.rm = TRUE),
                                        landscapeBValid = unique(landscapeBValid)),
                                   by = .(year, rep, pixelIndex)])

## calculate relative abundances of dominant species and melt data.
plotData[, landRelativeAbund := sum(B, na.rm = TRUE)/unique(landscapeB),
         by = .(year, rep, vegTypeID2)]
plotData[, landRelativeAbundValid := sum(BValid, na.rm = TRUE)/unique(landscapeBValid),
         by = .(year, rep, vegTypeID2Valid)]
plotData1 <- plotData[, .(year, rep, pixelIndex, vegTypeID2, B, landscapeB, landRelativeAbund)]
plotData1[, dataType := "simulated"]
plotData2 <- plotData[, .(year, rep, pixelIndex, vegTypeID2Valid, BValid,
                          landscapeBValid, landRelativeAbundValid)]
plotData2[, dataType := "observed"]
setnames(plotData2, grep("Valid", names(plotData2), value = TRUE),
         sub("Valid", "", grep("Valid", names(plotData2), value = TRUE)))

plotData <- rbind(plotData1, plotData2)
rm(plotData1, plotData2)

plot3 <- ggplot(data = plotData[dataType == "simulated"],
                aes(x = vegTypeID2, y = landRelativeAbund)) +
  stat_summary(fun.y = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = plotData[dataType != "simulated"],
               aes(x = vegTypeID2, y = landRelativeAbund),
               fun.y = "mean", geom = "point", size = 2) +
  scale_x_discrete(labels = speciesLabels) +
  theme_pubr() +
  facet_wrap(~ year, labeller = labeller(year = c("1" = "2001", "11" = "2011"))) +
  labs(title = "Landscape-wide dominant species",
       x = "", y = expression(atop("relative abundance", over("species B", "total B"))))

## no. pixels with a certain dominant species
## HERE: simualted and observed differ in no. of pixels in year 1... WHY?!
## possibly due to Biomass being adjusted?
plot4 <- ggplot(data = plotData[dataType == "simulated"],
                aes(x = vegTypeID2)) +
  geom_bar() +
  geom_point(data = plotData[dataType != "simulated"],
             aes(x = vegTypeID2), stat = "count", size = 2) +
  scale_x_discrete(labels = speciesLabels) +
  theme_pubr() +
  facet_wrap(~ year, labeller = labeller(year = c("1" = "2001", "11" = "2011"))) +
  labs(title = "Landscape-wide dominant species",
       x = "", y = "no. pixels",
       fill = "")

## as previous, but in relative terms
plot5 <- ggplot(data = plotData, aes(x = dataType, fill = vegTypeID2)) +
  geom_histogram(stat = "count", position = "fill") +
  scale_fill_brewer(palette = "Accent", labels = speciesLabels) +
  theme_pubr() +
  facet_wrap(~ year, labeller = labeller(year = c("1" = "2001", "11" = "2011"))) +
  labs(title = "Landscape-wide species dominance", x = "",
       y = "proportion of pixels", fill = "dominant species")

grid.arrange(plot1, plot2, plot3, plot4,# plot5,
             nrow = 2, ncol = 2)

## COMPARISONS OF DELTA PER PIXEL-------------------
## per species
plotData <- standCohortData[year %in% c(1,11),]
plotData <- plotData[, list(deltaB = B[which(year == 11)] - B[which(year == 1)],
                            deltaBValid = BValid[which(year == 11)] - BValid[which(year == 1)]),
                     , by = .(rep, pixelIndex, speciesCode)]
plotData <- melt(plotData, measure.vars = c("deltaB", "deltaBValid"),
                 variable.name = "dataType", value.name = "deltaB")

plot6 <- ggplot(data = plotData,
                aes(x = speciesCode, y = deltaB, fill = dataType)) +
  geom_boxplot() +
  scale_x_discrete(labels = speciesLabels) +
  scale_fill_discrete(labels = c("deltaB" = "simulated",
                                 "deltaBValid" = "observed")) +
  theme_pubr() +
  labs(title = "Pixel-level changes in abundance", fill = "",
       x = "", y = expression(paste(Delta, "B")))

plot6.2 <- ggplot(data = plotData[dataType == "deltaB"],
                  aes(x = speciesCode, y = deltaB)) +
  geom_boxplot() +
  stat_summary(data = plotData[dataType == "deltaBValid"],
               aes(x = speciesCode, y = deltaB, group = speciesCode),
               fun.y = "mean", geom = "point",
               size = 2.5, shape = "square", col = "red") +
  scale_x_discrete(labels = speciesLabels) +
  theme_pubr() +
  labs(title = "Pixel-level changes in abundance", fill = "",
       x = "", y = expression(paste(Delta, "B")))


plotData <- plotData[, mean(deltaB) , by = .(rep, speciesCode, dataType)]

plot7 <- ggplot(data = plotData[dataType == "deltaB"],
                aes(x = speciesCode, y = V1)) +
  stat_summary(fun.y = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = plotData[dataType == "deltaBValid"],
               aes(x = speciesCode, y = V1),
               fun.y = "mean", geom = "point", size = 2) +
  scale_x_discrete(labels = speciesLabels) +
  theme_pubr() +
  labs(title = "Landscape-averaged changes in abundance",
       x = "", y = expression(paste(Delta, "B")))

grid.arrange(plot6, plot7,
             ncol = 2)

## map delta b
plotData <- standCohortData[year %in% c(1,11),]
plotData <- plotData[, list(deltaB = B[which(year == 11)] - B[which(year == 1)],
                            deltaBValid = BValid[which(year == 11)] - BValid[which(year == 1)]),
                     , by = .(rep, pixelIndex, speciesCode)]

sppDeltaBValid <- lapply(unique(plotData$speciesCode), FUN = function(sp) {
  ## subset data
  tempData <- unique(plotData[speciesCode == sp, .(pixelIndex, deltaBValid)])

  ## make raster
  deltaBValid <- pixelGroupMapStk[[1]]
  deltaBValid[] <- NA

  deltaBValid[tempData$pixelIndex] <- tempData$deltaBValid

  names(deltaBValid) <- sp
  return(deltaBValid)
})

sppDeltaBValid <- stack(sppDeltaBValid)
plot(sppDeltaBValid)
