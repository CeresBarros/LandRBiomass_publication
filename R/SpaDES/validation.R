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
options("reproducible.useGDAL" = FALSE)
# runName <- "studyAreaS"
# runName <- "studyAreaL"
runName <- "parametriseSALarge"

## select simulation reps and years
reps <- 1:10
year1 <- 2001
year2 <- 2011

## paths
simDirName <- "dec2020Runs"
simPaths <- list(cachePath = file.path("R/SpaDES/cache", simDirName)
                 , modulePath = file.path("R/SpaDES/m")
                 , inputPath = file.path("R/SpaDES/inputs")
                 , outputPath = file.path("R/SpaDES/outputs", simDirName))

## define local and remote repo paths for file name subs:
localProjPath <- getwd()
remoteProjPath <- "/home/cbarros/GitHub/LandRBiomass_publication"

## init simList
simListInit <- readRDS(file.path(simPaths$outputPath, paste0("simInitList_", runName, ".rds")))

## get files of validation year per rep
## 2001 is the first year (0) and 2011, is year 10,
## but choose years 2001 and 2011 because objects were saved at the beginning of the year, they reflect the end of the previous year
## cohortData/pixelGroupMap, however, will vary across reps for year 1, because the beginning of year 1 already includes growth/mortality from year 0.
## so, using year 2001 and year 2011 (equivalent to 2010)
outputFiles <- readRDS(file.path(simPaths$outputPath, paste0("outputFiles_", runName, ".rds")))
outputFiles <- outputFiles[saveTime %in% c(year1, year2)]

## change paths to local
if (all(!grepl(getwd(), outputFiles$file)))
  outputFiles$file <- sub(paste0(".*", basename(getwd())),
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

## LOAD VALIDATION LAYERS - they are saved in the outputs once for each rep.
## here the rstDisturbed only has 1s
speciesLayersInit <- if (inMemory(simListInit$speciesLayers))
  simListInit$speciesLayers else {
    sppCodes <- as.character(unique(allCohortData$speciesCode))
    sppCodesStr <- paste(sppCodes, collapse = "|")

    speciesLayersInit <- lapply(simListInit$speciesLayers@layers, function(x) {
      localFileName <- list.files(simPaths$inputPath, basename(x@file@name),
                                  recursive = TRUE, full.names = TRUE)
      if (!length(localFileName))
        localFileName <- list.files(simPaths$outputPath, basename(x@file@name),
                                    recursive = TRUE, full.names = TRUE)
      sppLayer <- raster(localFileName)

      names(sppLayer) <- sub(paste0("(.*)(", sppCodesStr, ")(.*)"), "\\2",
                             basename(localFileName))
      return(sppLayer)
    }
    )
    rm(sppCodes, sppCodesStr)
    stack(speciesLayersInit)
  }

if (!compareRaster(speciesLayersInit, simListInit$rasterToMatch, stopiffalse = FALSE)) {
  speciesLayersInit <- postProcess(speciesLayersInit, rasterToMatch = simListInit$rasterToMatch)
}

standAgeMapInit <- if (inMemory(simListInit$standAgeMap))
  simListInit$standAgeMap else {
    localFileName <- list.files(simPaths$inputPath, basename(simListInit$standAgeMap@file@name),
                                recursive = TRUE, full.names = TRUE)
    raster(localFileName)
  }

if (!compareRaster(standAgeMapInit, simListInit$rasterToMatch, stopiffalse = FALSE)) {
  standAgeMapInit <- postProcess(standAgeMapInit, rasterToMatch = simListInit$rasterToMatch)
}

rawBiomassMapInit <- if (inMemory(simListInit$rawBiomassMap))
  simListInit$rawBiomassMap else {
    localFileName <- list.files(simPaths$inputPath, basename(simListInit$rawBiomassMap@file@name),
                                recursive = TRUE, full.names = TRUE)
    raster(localFileName)
  }

if (!compareRaster(rawBiomassMapInit, simListInit$rasterToMatch, stopiffalse = FALSE)) {
  rawBiomassMapInit <- postProcess(rawBiomassMapInit, rasterToMatch = simListInit$rasterToMatch)
}

rstDisturbedPix <- if (inMemory(simListInit$rstLCChange))
  simListInit$rstLCChange else {
    localFileName <- list.files(simPaths$inputPath, basename(simListInit$rstLCChange@file@name),
                                recursive = TRUE, full.names = TRUE)
    raster(localFileName)
  }

if (!compareRaster(rstDisturbedPix, simListInit$rasterToMatch, stopiffalse = FALSE)) {
  rstDisturbedPix <- postProcess(rstDisturbedPix, rasterToMatch = simListInit$rasterToMatch)
}

speciesLayersValidation <- if (inMemory(simListInit$speciesLayersValidation))
  simListInit$speciesLayersValidation else {
    sppCodes <- as.character(unique(allCohortData$speciesCode))
    sppCodesStr <- paste(sppCodes, collapse = "|")

    speciesLayersValidation <- lapply(simListInit$speciesLayersValidation@layers, function(x) {
      localFileName <- list.files(simPaths$inputPath, basename(x@file@name),
                                  recursive = TRUE, full.names = TRUE)
      if (!length(localFileName))
        localFileName <- list.files(simPaths$outputPath, basename(x@file@name),
                                    recursive = TRUE, full.names = TRUE)
      sppLayer <- raster(localFileName)

      names(sppLayer) <- sub(paste0("(.*)(", sppCodesStr, ")(.*)"), "\\2",
                             basename(localFileName))
      return(sppLayer)
    }
    )
    rm(sppCodes, sppCodesStr)
    stack(speciesLayersValidation)
  }
if (!compareRaster(speciesLayersValidation, simListInit$rasterToMatch, stopiffalse = FALSE)) {
  speciesLayersValidation <- postProcess(speciesLayersValidation, rasterToMatch = simListInit$rasterToMatch)
}

standAgeMapValidation <- if (inMemory(simListInit$standAgeMapValidation))
  simListInit$standAgeMapValidation else {
    localFileName <- list.files(simPaths$inputPath, basename(simListInit$standAgeMapValidation@file@name),
                                recursive = TRUE, full.names = TRUE)
    raster(localFileName)
  }
if (!compareRaster(standAgeMapValidation, simListInit$rasterToMatch, stopiffalse = FALSE)) {
  standAgeMapValidation <- postProcess(standAgeMapValidation, rasterToMatch = simListInit$rasterToMatch)
}

rawBiomassMapValidation <- if (inMemory(simListInit$rawBiomassMapValidation))
  simListInit$rawBiomassMapValidation else {
    localFileName <- list.files(simPaths$inputPath, basename(simListInit$rawBiomassMapValidation@file@name),
                                recursive = TRUE, full.names = TRUE)
    raster(localFileName)
  }
if (!compareRaster(rawBiomassMapValidation, simListInit$rasterToMatch, stopiffalse = FALSE)) {
  rawBiomassMapValidation <- postProcess(rawBiomassMapValidation, rasterToMatch = simListInit$rasterToMatch)
}

## check names of species layers match
if (!all(names(speciesLayersInit) %in% names(speciesLayersValidation)))
  stop("validation species layers names do not match!")


## MAKE VALIDATION PIXEL TABLES
## need to reproduce some steps in BBDP
## year 2001
pixelTable <- makePixelTable(speciesLayers = speciesLayersInit,
                             biomassMap = rawBiomassMapInit,
                             standAgeMap = standAgeMapInit,
                             rasterToMatch = simListInit$rasterToMatch)
pixelTable[, initialEcoregionCode := NULL]
pixelTable <- unique(pixelTable)

minCoverThreshold <- P(simListInit)$Biomass_borealDataPrep$minCoverThreshold
deciduousCoverDiscount <- P(simListInit)$Biomass_borealDataPrep$deciduousCoverDiscount
pixelGroupBiomassClass <- P(simListInit)$Biomass_borealDataPrep$pixelGroupBiomassClass
validationDataInit <- LandR:::.createCohortData(pixelTable, rescale = TRUE,
                                                minCoverThreshold = minCoverThreshold)
validationDataInit <- simListInit$.mods$Biomass_borealDataPrep$partitionBiomass(x = deciduousCoverDiscount, validationDataInit)
set(validationDataInit, NULL, "B", asInteger(validationDataInit$B/pixelGroupBiomassClass) * pixelGroupBiomassClass)
set(validationDataInit, NULL, c("decid", "cover2", "logAge"), NULL)
set(validationDataInit, NULL, "cover", asInteger(validationDataInit$cover))

## remove pixels with missing data
validationDataInit <- validationDataInit[!is.na(B)]

## calculate relative B
validationDataInit[, relativeAbundValid := B/totalBiomass,
                   by = .(pixelIndex, speciesCode)]
validationDataInit[totalBiomass == 0, relativeAbundValid := 0]

rm(pixelTable); amc::.gc()

## year 2011
pixelTable <- makePixelTable(speciesLayers = speciesLayersValidation,
                             biomassMap = rawBiomassMapValidation,
                             standAgeMap = standAgeMapValidation,
                             rasterToMatch = simListInit$rasterToMatch)

pixelTable[, initialEcoregionCode := NULL]
pixelTable <- unique(pixelTable)

validationData <- LandR:::.createCohortData(pixelTable, rescale = TRUE,
                                            minCoverThreshold = minCoverThreshold)
validationData <- simListInit$.mods$Biomass_borealDataPrep$partitionBiomass(x = deciduousCoverDiscount, validationData)
set(validationData, NULL, "B", asInteger(validationData$B/pixelGroupBiomassClass) * pixelGroupBiomassClass)
set(validationData, NULL, c("decid", "cover2", "logAge"), NULL)
set(validationData, NULL, "cover", asInteger(validationData$cover))

## remove pixels with missing data
validationData <- validationData[!is.na(B)]

## calculate relative B
validationData[, relativeAbundValid := B/totalBiomass,
                   by = .(pixelIndex, speciesCode)]
validationData[totalBiomass == 0, relativeAbundValid := 0]

rm(pixelTable); amc::.gc()

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
                                   by = .(rep, year, pixelIndex)]
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
                                                   rep = reps,
                                                   year = year1)))

combinationsValid <- as.data.table(expand.grid(list(speciesCode = unique(standCohortData$speciesCode),
                                                    pixelIndex = unique(c(validationData$pixelIndex, standCohortData$pixelIndex)),
                                                    rep = reps,
                                                    year = year2)))
validationDataInit <- validationDataInit[combinationsInit,
                                         on = c("pixelIndex", "speciesCode")]
validationData <- validationData[combinationsValid,
                                 on = c("pixelIndex", "speciesCode")]
validationData <- rbindlist(list(validationDataInit, validationData),
                            use.names = TRUE)

## exclude pixels that are not simulated
validationData <- validationData[pixelIndex %in% standCohortData$pixelIndex]

## clean up
rm(combinationsInit, combinationsValid, validationDataInit)

## change names before joining.
## exclude simulated data that has no match on validation data
setnames(validationData, c("cover", "age", "B", "totalBiomass"),
         c("coverValid", "standAgeValid", "BValid", "standBValid"))
cols <- c("rep", "year", "pixelIndex", "speciesCode", "coverValid",
          "standAgeValid", "BValid", "standBValid", "relativeAbundValid")
standCohortData <- standCohortData[validationData[, ..cols],
                                   on = c("rep", "year", "pixelIndex", "speciesCode"),
                                   nomatch = 0]

## some NAs added on simulated data for year 2011
# standCohortData[is.na(BValid), BValid := 0]  ## just to be sure
standCohortData[, standAge := max(standAge, na.rm = TRUE), by = .(year, rep, pixelIndex)]

## calculate simulated standAge, standB and relative B now that NAs have been taken care of
standCohortData[, standB := asInteger(sum(B, na.rm = TRUE)),
                .(rep, year, pixelIndex)]
standCohortData[, relativeAbund := B/standB,
                by = .(rep, year, pixelIndex, speciesCode)]

## reorder column names
cols <-  c(grep("Valid", names(standCohortData), value = TRUE, invert = TRUE),
           grep("Valid", names(standCohortData), value = TRUE))
standCohortData <- standCohortData[, ..cols]

## remove disturbed pixels
if (!compareRaster(rstDisturbedPix, pixelGroupMapStk[[1]])) ## just check rasters match first
  stop("rstDisturbedPix and pixelGroupMap differ!")
disturbedPix <- which(!is.na(getValues(rstDisturbedPix)))
standCohortData <- standCohortData[!pixelIndex %in% disturbedPix]

## calculate some landscape metrics
standCohortData[, `:=`(landscapeB = sum(B, na.rm = TRUE),
                       landscapeBValid = sum(BValid, na.rm = TRUE)),
                by = .(rep, year)]

## plot labels
speciesLabels <- LandR::equivalentName(unique(standCohortData$speciesCode), simListInit$sppEquiv,
                                       column = "EN_generic_short")
names(speciesLabels) <- unique(standCohortData$speciesCode)

## TESTS
## species relative abundances do not sum to 1 across the landscape for year 11
standCohortData[, mean(relativeAbundValid, na.rm = TRUE), by = .(rep, year, speciesCode)] %>%
  .[, sum(V1), by = .(rep, year)]

## nor for each stand (but they are close) - probably due to deciduous/conifer cover adjustments
standCohortData[, sum(relativeAbundValid, na.rm = TRUE), by = .(rep, year, pixelIndex)] %>%
  lapply(., unique) %>% .["V1"]

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

plot1 <- ggplot(standDeltaValidData,
       aes(x = standDeltaBValid)) +
  geom_histogram() +
  labs(title = expression(paste("observed" ~~ Delta, "B"))) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(axis.title.y = element_blank())

plot2 <- ggplot(standDeltaValidData,
       aes(y = standDeltaBValid)) +
  geom_boxplot() +
  coord_flip() +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(axis.title.y = element_blank())
ggarrange(plot1, plot2, ncol = 1, heights = c(1, 0.5))


plot1 <- ggplot(standDeltaValidData,
                aes(x = standDeltaAgeValid)) +
  geom_histogram() +
  labs(title = expression(paste("observed" ~~ Delta, "age"))) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(axis.title.y = element_blank())

plot2 <- ggplot(standDeltaValidData,
                aes(y = standDeltaAgeValid)) +
  geom_boxplot() +
  coord_flip() +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(axis.title.y = element_blank())
ggarrange(plot1, plot2, ncol = 1, heights = c(1, 0.5))

## delta biomass for the "supposed" age increment - all over the place
plot1 <- ggplot(standDeltaValidData[standDeltaAgeValid == 10],
                aes(x = standDeltaBValid)) +
  geom_histogram() +
  labs(title = expression(atop(paste("observed" ~~ Delta, "B"),
                          paste("for" ~~ Delta, "age" == 10)))) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(axis.title.y = element_blank())

plot2 <- ggplot(standDeltaValidData[standDeltaAgeValid == 10],
                aes(y = standDeltaBValid)) +
  geom_boxplot() +
  coord_flip() +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(axis.title.y = element_blank())
ggarrange(plot1, plot2, ncol = 1, heights = c(1, 0.5))


## MAPS OF OBSERVED CHANGES IN SPECIES BIOMASS AFTER ADJUSTMENTS -------
## by adjustments we mean, the data cleanup by BBDP
plotData <- standCohortData[year %in% c(year1, year2),]
plotData <- plotData[, list(standDeltaB = unique(standB[which(year == year2)]) - unique(standB[which(year == year1)]),
                            standDeltaBValid = unique(standBValid[which(year == year2)]) - unique(standBValid[which(year == year1)]),
                            standDeltaAge = unique(standAge[which(year == year2)]) - unique(standAge[which(year == year1)]),
                            standDeltaAgeValid = unique(standAgeValid[which(year == year2)]) - unique(standAgeValid[which(year == year1)])),
                     , by = .(rep, pixelIndex)]

standDeltaBValidAdj <- pixelGroupMapStk[[1]]
standDeltaBValidAdj[] <- NA
standDeltaBValidAdj[unique(plotData[,.(pixelIndex, standDeltaBValid)])[, pixelIndex]] <- unique(plotData[,.(pixelIndex, standDeltaBValid)])[, standDeltaBValid]

standDeltaAgeValidAdj <- pixelGroupMapStk[[1]]
standDeltaAgeValidAdj[] <- NA
standDeltaAgeValidAdj[unique(plotData[,.(pixelIndex, standDeltaAgeValid)])[, pixelIndex]] <- unique(plotData[,.(pixelIndex, standDeltaAgeValid)])[, standDeltaAgeValid]

Plot(standDeltaBValidAdj, standDeltaAgeValidAdj)

## what is the relationship between the two?
## a little bit better after cleanup/corrections
ggplot(plotData,
       aes(x = standDeltaAgeValid, y = standDeltaBValid)) +
  geom_point() +
  stat_smooth(method = "lm")
summary(lm(standDeltaBValid ~ standDeltaAgeValid, data = plotData))

plot1 <- ggplot(plotData,
                aes(x = standDeltaBValid)) +
  geom_histogram() +
  labs(title = expression(paste("observed" ~~ Delta, "B" ~~ "(adjusted)"))) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(axis.title.y = element_blank())

plot2 <- ggplot(plotData,
                aes(y = standDeltaBValid)) +
  geom_boxplot() +
  coord_flip() +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(axis.title.y = element_blank())
ggarrange(plot1, plot2, ncol = 1, heights = c(1, 0.5))

plot1 <- ggplot(plotData,
                aes(x = standDeltaAgeValid)) +
  geom_histogram() +
  labs(title = expression(paste("observed" ~~ Delta, "age" ~~ "(adjusted)"))) +

  theme_pubr(base_size = 12, margin = FALSE) +
  theme(axis.title.y = element_blank())

plot2 <- ggplot(plotData,
                aes(y = standDeltaAgeValid)) +
  geom_boxplot() +
  coord_flip() +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(axis.title.y = element_blank())
ggarrange(plot1, plot2, ncol = 1, heights = c(1, 0.5))

## delta biomass for the "supposed" age increment - all over the place
plot1 <- ggplot(plotData[standDeltaAgeValid == 10],
                aes(x = standDeltaBValid)) +
  geom_histogram() +
  labs(title = expression(atop(paste("observed" ~~ Delta, "B" ~~ "(adjusted)"),
                               paste("for" ~~ Delta, "age" == 10)))) +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(axis.title.y = element_blank())

plot2 <- ggplot(plotData[standDeltaAgeValid == 10],
                aes(y = standDeltaBValid)) +
  geom_boxplot() +
  coord_flip() +
  theme_pubr(base_size = 12, margin = FALSE) +
  theme(axis.title.y = element_blank())
ggarrange(plot1, plot2, ncol = 1, heights = c(1, 0.5))


## MAPS OF OBSERVED CHANGES IN SPECIES BIOMASS AFTER ADJUSTMENTS -------
plotData <- standCohortData[year %in% c(year1, year2),]
plotData <- plotData[, list(deltaB = B[which(year == year2)] - B[which(year == year1)],
                            deltaBValid = BValid[which(year == year2)] - BValid[which(year == year1)]),
                     , by = .(rep, pixelIndex, speciesCode)]

sppDeltaB <- lapply(unique(plotData$speciesCode), FUN = function(sp) {
  ## subset data
  tempData <- unique(plotData[speciesCode == sp, .(pixelIndex, deltaB)])

  ## make raster
  deltaB <- pixelGroupMapStk[[1]]
  deltaB[] <- NA

  deltaB[tempData$pixelIndex] <- tempData$deltaB

  names(deltaB) <- sp
  return(deltaB)
}) %>%
  stack(.)

sppDeltaBValid <- lapply(unique(plotData$speciesCode), FUN = function(sp) {
  ## subset data
  tempData <- unique(plotData[speciesCode == sp, .(pixelIndex, deltaBValid)])

  ## make raster
  deltaBValid <- pixelGroupMapStk[[1]]
  deltaBValid[] <- NA

  deltaBValid[tempData$pixelIndex] <- tempData$deltaBValid

  names(deltaBValid) <- sp
  return(deltaBValid)
}) %>%
  stack(.)

Plot(sppDeltaB, sppDeltaBValid)


## TESTS:
## remove pixels were the aging was too extreme (lower/higher than 25% and 75% quantiles)
# (this didn't work, no changes)
# tempDT <- unique(standCohortData[, .(year, rep, pixelIndex, standAgeValid)])
# tempDT <- tempDT[, .(standDeltaAgeValid = standAgeValid[year == year2] - standAgeValid[year == year1]),
#                  by = .(rep, pixelIndex)]
# quants <- quantile(tempDT$standDeltaAgeValid, probs = c(0.25, 0.75))
# pixToKeep <- unique(tempDT[standDeltaAgeValid > quants["25%"] &
#                       standDeltaAgeValid < quants["75%"],
#                       pixelIndex])
## remove pixels were observed age changes were not 10yrs.
## this yields only 63 pixels...
# tempDT <- unique(standCohortData[, .(year, rep, pixelIndex, standAgeValid)])
# tempDT <- tempDT[, .(standDeltaAgeValid = standAgeValid[year == year2] - standAgeValid[year == year1]),
#                  by = .(rep, pixelIndex)]
# pixToKeep <- unique(tempDT[standDeltaAgeValid == 10,
#                            pixelIndex])

## Remove pixels were stand age or B decreased, as we cannot account for disturbances that may not have been captured from sat data
## or measurement errors that yielded too high B in 2001
tempDT <- unique(standCohortData[, .(year, rep, pixelIndex, standAgeValid, standBValid)])
tempDT <- tempDT[, .(standDeltaAgeValid = standAgeValid[year == year2] - standAgeValid[year == year1],
                     standDeltaBValid = standBValid[year == year2] - standBValid[year == year1]),
                 by = .(rep, pixelIndex)]
pixToKeep <- unique(tempDT[standDeltaBValid > 0 & standDeltaAgeValid > 0,
                    pixelIndex])

if (!exists("pixToKeep"))
  pixToKeep <- unique(standCohortData$pixelIndex)

## Make a general labels for years
yearLabels <- c("2001", "2011")
names(yearLabels) <- as.character(c(year1, year2))

## LANDSCAPE-WIDE COMPARISONS IN A GIVEN YEAR --------------------
## relative abundance (biomass) per species
## note that error bars are absent because there is literally no variation among reps at year 2011:
if (FALSE) {
  test1 <- readRDS("R/SpaDES/outputs/dec2020Runs/sim1_rep01/cohortData_year2011.rds")
  test2 <- readRDS("R/SpaDES/outputs/dec2020Runs/sim1_rep04/cohortData_year2011.rds")
  identical(test1, test2)

  test3 <- readRDS("R/SpaDES/outputs/dec2020Runs/sim1_rep01/cohortData_year2016.rds")
  test4 <- readRDS("R/SpaDES/outputs/dec2020Runs/sim1_rep04/cohortData_year2016.rds")
  identical(test3, test4)

  ## still no variation even at the end of the 30 yrs of simulation
  test5 <- readRDS("R/SpaDES/outputs/dec2020Runs/sim1_rep01/cohortData_year2031.rds")
  test6 <- readRDS("R/SpaDES/outputs/dec2020Runs/sim1_rep04/cohortData_year2031.rds")
  identical(test5, test6)

  rm(test1, test2, test3, test4, test5, test6)
}

plotData <- standCohortData[pixelIndex %in% pixToKeep]
plotData <- plotData[year %in% c(year1, year2),
                     list(landRelativeAbund = sum(B, na.rm = TRUE)/unique(landscapeB),
                          landRelativeAbundValid = sum(BValid, na.rm = TRUE)/unique(landscapeBValid)),
                     by = .(rep, year, speciesCode)]
plot1 <- ggplot(data = plotData,
                aes(x = speciesCode, y = landRelativeAbund)) +
  stat_summary(fun = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(aes(y = landRelativeAbundValid, colour = "observed"),
               fun = "mean", geom = "point", size = 2) +
  scale_x_discrete(labels = speciesLabels, drop = FALSE) +
  scale_color_manual(values = c("observed" = "black")) +
  theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
  theme(legend.position = "right") +
  facet_wrap(~ year, labeller = labeller(year = yearLabels)) +
  labs(title = "Species relative abundances",
       x = "", y = expression(over("species B", "total B")),
       colour = "")

## no. pixels with a species
plotData <- standCohortData[pixelIndex %in% pixToKeep]
plotData <- plotData[year %in% c(year1, year2),
                     list(count = sum(B > 0, na.rm = TRUE),
                          countValid = sum(BValid > 0, na.rm = TRUE)),
                     by = .(rep, year, speciesCode)]
plotData <- melt(plotData, measure.vars = c("count", "countValid"),
                 variable.name = "dataType", value.name = "count")

plot2 <- ggplot(data = plotData[dataType == "count"],
                aes(x = speciesCode, y = count)) +
  stat_summary(fun = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = plotData[dataType == "countValid"],
               aes(x = speciesCode, y = count, colour = "observed"),
               fun = "mean", geom = "point", size = 2) +
  scale_x_discrete(labels = speciesLabels, drop = FALSE) +
  scale_color_manual(values = c("observed" = "black")) +
  theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
  theme(legend.position = "right") +
  facet_wrap(~ year, labeller = labeller(year = yearLabels)) +
  labs(title = "Species presences",
       x = "", y = "no. pixels",
       colour = "", fill = "")

## landscape-wide relative abundance per dominant species
## get dominant species - these match with model outputs, I checked
## note: don't use melt, because dominant spp differ between valid and simul data.
## note2: mixed pixels get a "mixed" type
## calculate the no. of "dominant" species so that pixels with more than one dominant spp are considered "mixed"
tempDT <- standCohortData[pixelIndex %in% pixToKeep]
tempDT <- tempDT[year %in% c(year1, year2),
                 list(noDoms = sum(relativeAbund == max(relativeAbund, na.rm = TRUE), na.rm = TRUE),
                      noDomsValid = sum(relativeAbundValid == max(relativeAbundValid, na.rm = TRUE), na.rm = TRUE)),
                 by = .(year, rep, pixelIndex)]

plotData <- standCohortData[pixelIndex %in% pixToKeep]
plotData <- unique(plotData[year %in% c(year1, year2),
                            list(vegType = speciesCode[which.max(relativeAbund)],
                                 relativeAbund = max(relativeAbund, na.rm = TRUE),
                                 vegTypeValid = speciesCode[which.max(relativeAbundValid)],
                                 relativeAbundValid = max(relativeAbundValid, na.rm = TRUE)),
                            by = .(year, rep, pixelIndex)])
plotData <- tempDT[plotData, on = .(year, rep, pixelIndex)]

plotData1 <- unique(plotData[, .(year, rep, pixelIndex, vegType, relativeAbund, noDoms)])
plotData1[, dataType := "relativeAbund"]
plotData2 <- unique(plotData[, .(year, rep, pixelIndex, vegTypeValid, relativeAbundValid, noDomsValid)])
plotData2[, dataType := "relativeAbundValid"]
setnames(plotData2, c("vegTypeValid", "relativeAbundValid", "noDomsValid"),
         c("vegType", "relativeAbund", "noDoms"))

plotData <- rbind(plotData1, plotData2)
rm(plotData1, plotData2)

plotData[noDoms > 1, vegType := "Mixed"]

## calculate landscape average
plotData1 <- plotData[, relativeAbund := mean(relativeAbund),
                      by = .(rep, year, dataType, vegType)]
plot3 <- ggplot(data = plotData1[!grepl("Valid", dataType)],
                aes(x = vegType, y = relativeAbund)) +
  stat_summary(fun = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = plotData1[grepl("Valid", dataType)],
               aes(x = vegType, y = relativeAbund, colour = "observed"),
               fun = "mean", geom = "point", size = 2) +
  scale_x_discrete(labels = speciesLabels, drop = FALSE) +
  scale_color_manual(values = c("observed" = "black")) +
  theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
  theme(legend.position = "right") +
  facet_wrap(~ year, labeller = labeller(year = yearLabels)) +
  labs(title = "Dominant species' relative abundances",
       x = "", y = expression(over("species B", "total B")),
       colour = "")

## no. pixels with a certain dominant species
## HERE: simulated and observed differ in no. of pixels in year 1...
## this is because B is adjusted using a statistical model
plotData1 <- plotData[, list(count = length(pixelIndex)),
                      by = .(rep, year, dataType, vegType)]
plot4 <- ggplot(data = plotData1[!grepl("Valid", dataType)],
                aes(x = vegType, y = count)) +
  stat_summary(fun = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  geom_point(data = plotData1[grepl("Valid", dataType)],
             aes(x = vegType, y = count, colour = "observed"), size = 2) +
  scale_x_discrete(labels = speciesLabels, drop = FALSE) +
  scale_color_manual(values = c("observed" = "black")) +
  theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
  theme(legend.position = "right") +
  facet_wrap(~ year, labeller = labeller(year = yearLabels)) +
  labs(title = "Dominant species' presences",
       x = "", y = "no. pixels", fill = "", colour = "")

## as previous, but in relative terms
plot5 <- ggplot(data = plotData1, aes(x = dataType, y = count, fill = vegType)) +
  stat_summary(fun = "mean", geom = "bar", position = "fill") +
  scale_fill_brewer(palette = "Accent", labels = speciesLabels) +
  scale_x_discrete(labels = c("relativeAbund" = "simulated", "relativeAbundValid" = "observed")) +
  theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
  theme(legend.position = "right") +
  facet_wrap(~ year, labeller = labeller(year = yearLabels)) +
  labs(title = "Dominant species presences", x = "",
       y = "proportion of pixels", fill = "")

plot1_4 <- ggarrange(plot1 + scale_y_continuous(limits = c(0,1)),
                     plot2 + scale_y_continuous(limits = c(0,2000)),
                     plot3 + scale_y_continuous(limits = c(0,1)),
                     plot4 + scale_y_continuous(limits = c(0,2000)),
                     common.legend = TRUE, legend = "bottom",
                     nrow = 2, ncol = 2)
annotate_figure(plot1_4,
                top = text_grob("Landscape-averaged comparisons", size = 16))

## PIXEL-LEVEL COMPARISONS IN A GIVEN YEAR --------------------
## pixel-level relative abundances per species
plotData <- standCohortData[pixelIndex %in% pixToKeep]
plotData <- plotData[year %in% c(year1, year2),
                     .(rep, year, pixelIndex, speciesCode,
                       relativeAbund, relativeAbundValid)]
plotData <- melt(plotData,
                 id.vars = c("rep", "year", "pixelIndex", "speciesCode"))
plot6 <- ggplot(data = plotData,
                aes(x = speciesCode, y = value, fill = variable)) +
  geom_boxplot() +
  scale_x_discrete(labels = speciesLabels, drop = FALSE) +
  scale_fill_discrete(labels = c("relativeAbund" = "simulated",
                                 "relativeAbundValid" = "observed")) +
  theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
  facet_wrap(~ year, labeller = labeller(year = yearLabels)) +
  labs(title = "Species relative abundances", fill = "",
       x = "", y = expression(over("species B", "stand B")))

## pixel-level relative abundance per dominant species
## get dominant species (those with maxB) - these match with model outputs, I checked
tempDT <- standCohortData[pixelIndex %in% pixToKeep]
tempDT <- tempDT[year %in% c(year1, year2),
                 list(noDoms = sum(relativeAbund == max(relativeAbund, na.rm = TRUE), na.rm = TRUE),
                      noDomsValid = sum(relativeAbundValid == max(relativeAbundValid, na.rm = TRUE), na.rm = TRUE)),
                 by = .(year, rep, pixelIndex)]

plotData <- standCohortData[pixelIndex %in% pixToKeep]
plotData <- unique(plotData[year %in% c(year1, year2),
                            list(vegType = speciesCode[which.max(relativeAbund)],
                                 relativeAbund = max(relativeAbund, na.rm = TRUE),
                                 vegTypeValid = speciesCode[which.max(relativeAbundValid)],
                                 relativeAbundValid = max(relativeAbundValid, na.rm = TRUE)),
                            by = .(year, rep, pixelIndex)])
plotData <- tempDT[plotData, on = .(year, rep, pixelIndex)]

plotData1 <- unique(plotData[, .(year, rep, pixelIndex, vegType, relativeAbund, noDoms)])
plotData1[, dataType := "relativeAbund"]
plotData2 <- unique(plotData[, .(year, rep, pixelIndex, vegTypeValid, relativeAbundValid, noDomsValid)])
plotData2[, dataType := "relativeAbundValid"]
setnames(plotData2, c("vegTypeValid", "relativeAbundValid", "noDomsValid"),
         c("vegType", "relativeAbund", "noDoms"))

plotData <- rbind(plotData1, plotData2)
rm(plotData1, plotData2)

plotData[noDoms > 1, vegType := "Mixed"]

plot7 <- ggplot(data = plotData,
                aes(x = vegType, y = relativeAbund, fill = dataType)) +
  geom_boxplot() +
  scale_fill_discrete(labels = c("relativeAbund" = "simulated",
                                 "relativeAbundValid" = "observed")) +
  scale_x_discrete(labels = speciesLabels, drop = FALSE) +
  theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
  facet_wrap(~ year, labeller = labeller(year = yearLabels)) +
  labs(title = "Dominant species' relative abundances",
       x = "", y = expression(over("species B", "total B")),
       fill = "")

plot6_7 <- ggarrange(plot6,
                     plot7 + labs(y = " \n "),
                     common.legend = TRUE, legend = "bottom",
                     ncol = 2)
annotate_figure(plot6_7,
                top = text_grob("Stand-level comparisons", size = 16))


## COMPARISONS OF DELTA PER PIXEL-------------------
## per species
plotData <- standCohortData[pixelIndex %in% pixToKeep]
plotData <- plotData[year %in% c(year1, year2)]
plotData <- plotData[, list(deltaB = B[which(year == year2)] - B[which(year == year1)],
                            deltaBValid = BValid[which(year == year2)] - BValid[which(year == year1)],
                            standDeltaB = standB[which(year == year2)] - standB[which(year == year1)],
                            standDeltaBValid = standBValid[which(year == year2)] - standBValid[which(year == year1)]),
                     , by = .(rep, pixelIndex, speciesCode)]

## melt spp and stand delta separately and rbind
cols <- grep("standDelta", names(plotData), invert = TRUE, value = TRUE)
plotData1 <- melt(plotData[, ..cols], measure.vars = c("deltaB", "deltaBValid"),
                  variable.name = "dataType", value.name = "deltaB")
cols <- grep("deltaB", names(plotData), invert = TRUE, value = TRUE)
plotData2 <- melt(plotData[, ..cols], measure.vars = c("standDeltaB", "standDeltaBValid"),
                  variable.name = "dataType", value.name = "deltaB")
plotData2[dataType == "standDeltaB", dataType := "deltaB"]
plotData2[dataType == "standDeltaBValid", dataType := "deltaBValid"]
plotData2 <- unique(plotData2[, speciesCode := "stand"])

plotData <- rbind(plotData1, plotData2, use.names = TRUE)

plot8 <-  ggplot(data = plotData[!grepl("Valid", dataType)],
                 aes(x = speciesCode, y = deltaB, group = rep)) +
  stat_summary(fun = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = plotData[grepl("Valid", dataType)],
               aes(x = speciesCode, y = deltaB, group = rep),
               fun = "mean", geom = "point", size = 2) +
  scale_x_discrete(labels = speciesLabels, drop = FALSE) +
  theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
  labs(title = "Landscape-averaged",
       x = "", y = expression(paste(Delta, "B")))

plot9 <- ggplot(data = plotData,
                aes(x = speciesCode, y = deltaB, fill = dataType)) +
  geom_boxplot() +
  scale_x_discrete(labels = speciesLabels) +
  scale_fill_discrete(labels = c("deltaB" = "simulated",
                                 "deltaBValid" = "observed")) +
  theme_pubr(base_size = 12, margin = FALSE, x.text.angle = 45) +
  labs(title = "Pixel-level", fill = "",
       x = "", y = expression(paste(Delta, "B")))

plot8_9 <- ggarrange(plot8, plot9 + labs(y = " \n "),
                     common.legend = TRUE, legend = "bottom",
                     ncol = 2)
png(file.path(figDir, "LandscapeStandComparisons_deltaB.png"), width = 10, height = 6,
    units = "in", res = 300)
annotate_figure(plot8_9,
                top = text_grob("Changes in abundance", size = 16))
graphics.off()
