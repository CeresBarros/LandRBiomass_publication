## ------------------------------------------------------
## LandR Biomass Publication - additional figures
##
## Ceres: July 2021
## ------------------------------------------------------

## this script should be sourced

## clean workspace
Require(c("raster", "sf", "data.table", "ggplot2",
          "ggspatial", "ggpubr"),
        upgrade = FALSE)

## GET ALL SIMLISTS ------------------------------------
runNames <- c("baseCase", "studyAreaChange", "altParameters")

if (!exists(simOutSpeciesLayers)) {
  simOutSpeciesLayers <- loadFromCache(speciesPaths$cachePath, cacheId = "f13730c0a3d0f3c8")
}

## pre simulation simLists
for (runName in runNames) {
  eval(parse(text = paste0("preSimList", runName, " <- loadSimList(file.path('R/SpaDES/outputs', simDirName, runName, paste0('simInit', runName)))")))
}

## post-simulation simLists
for (runName in runNames) {
  eval(parse(text = paste0("simList", runName, " <- qs::qread(file.path('R/SpaDES/outputs', simDirName, runName, paste0('simList_LandRBiomass_sim_', runName)))")))
}

## post-validation simLists
for (runName in runNames) {
  eval(parse(text = paste0("validSimList", runName, " <- qs::qread(file.path('R/SpaDES/outputs', simDirName, runName, paste0('simValid', runName)))")))
}

## ELAPSED TIMES ----------------------------------
elapsedTime(simOutSpeciesLayers)
elapsedTime(preSimListbaseCase)
elapsedTime(preSimListstudyAreaChange)
elapsedTime(preSimListaltParameters)

## for biomass_core use clocktimes as elapsed time seems to underestimate
timesbaseCase <- rbindlist(lapply(simListbaseCase, completed), idcol = "rep")
timesbaseCase <- timesbaseCase[eventType != "init", max(clockTime) - min(clockTime), by = .(rep, moduleName)]
timesbaseCase[, sum(as.integer(V1)), by = rep]
timesbaseCase[, mean(as.integer(V1))]

timesstudyAreaChange <- rbindlist(lapply(simListstudyAreaChange, completed), idcol = "rep")
timesstudyAreaChange <- timesstudyAreaChange[eventType != "init", max(clockTime) - min(clockTime), by = .(rep, moduleName)]
timesstudyAreaChange[, sum(as.integer(V1)), by = rep]
timesstudyAreaChange[, mean(as.integer(V1))]


timesaltParameters <- rbindlist(lapply(simListaltParameters, completed), idcol = "rep")
timesaltParameters <- timesaltParameters[eventType != "init", max(clockTime) - min(clockTime), by = .(rep, moduleName)]
timesaltParameters[, sum(as.integer(V1)), by = rep]
timesaltParameters[, mean(as.integer(V1))]

## OBJECT DIAGRAM
objectDiagram(preSimListbaseCase, width = 1000, height = 2500)
webshot("http://localhost:20581/session/viewhtml6aec1852565c/index.html",
        file = file.path(figDir, "objectDiagram_baseCase.png"))



## "ENSEMBLE" OF TWO INITIAL CONDITIONS ----------------------------------
baseCaseCohortData <- validSimListbaseCase$standCohortData[, .(year, rep, pixelIndex, speciesCode, B,
                                                               landscapeB, standB, relativeAbund)]
altParametersCohortData <- validSimListaltParameters$standCohortData[, .(year, rep, pixelIndex, speciesCode, B,
                                                               landscapeB, standB, relativeAbund)]
allCohortData <- rbind(baseCaseCohortData, altParametersCohortData, idcol = "scenario")
allCohortData[, scenario := ifelse(scenario == 1, "base case", "Example 2")]

plotData <- allCohortData[, list(B = sum(B)), by = .(speciesCode, year, rep, scenario)]
plotData <- plotData[, list(B = mean(B)), by = .(speciesCode, year, scenario)]

sppColoursVect <- preSimListbaseCase$sppColorVect
sppLabels <- equivalentName(names(sppColoursVect), preSimListbaseCase$sppEquiv, column = "EN_generic_short")
names(sppLabels) <- names(sppColoursVect)
sppLabels <- sppLabels[!is.na(sppLabels)]

ggplot(plotData, aes(x = year, y = B, colour = speciesCode, fill = speciesCode)) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon",
               alpha = 0.5, colour = "transparent") +
  stat_summary(fun = mean, geom = "line", size = 1) +
  scale_colour_manual(values = sppColoursVect,
                      labels = sppLabels) +
  scale_fill_manual(values = sppColoursVect,
                    labels = sppLabels) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_pubr(base_size = 12, legend = "right") +
  labs(colour = "", fill = "")
ggsave(file.path(figDir, "baseCase_altParameters_ensemble.png"), width = 7, height = 5, units = "in", dpi = 300)


## FIGURES ----------------------------------------
## Study areas
## get Canadian provinces and subset to SK
can1 <- raster::getData('GADM', country = "CAN", level = 1, path = tempdir())
can1 <- can1[can1$NAME_1 == "Saskatchewan",]
can1 <- spTransform(can1, CRSobj = crs(preSimListbaseCase$studyArea))

plotStudyAreas <- ggplot() +
  layer_spatial(can1, fill = "grey90") +
  layer_spatial(preSimListbaseCase$studyAreaLarge, fill = "darkgreen") +
  layer_spatial(preSimListbaseCase$studyArea, fill = "green") +
  layer_spatial(preSimListstudyAreaChange$studyAreaLarge, fill = "darkblue") +
  layer_spatial(preSimListstudyAreaChange$studyArea, fill = "lightblue") +
  labs(fill = "") +
  theme_pubr(margin = FALSE)
ggsave(file.path(figDir, "studyAreas_setsAB.png"), width = 7, height = 10, units = "in", dpi = 300)


## input stand biomass and age, ecological zonation (ecodistricts) and land-cover (LCC 2005)
rawBiomassMapSetA <- if (!inMemory(preSimListbaseCase$rawBiomassMap)) {
  list.files("R/SpaDES/inputs/",
             pattern = basename(filename(preSimListbaseCase$rawBiomassMap)),
             recursive = TRUE, full.names = TRUE) %>%
    raster(.)
} else {
  preSimListbaseCase$rawBiomassMap
}

standAgeMapSetA <- if (!inMemory(preSimListbaseCase$standAgeMap)) {
  list.files("R/SpaDES/inputs/",
             pattern = basename(filename(preSimListbaseCase$standAgeMap)),
             recursive = TRUE, full.names = TRUE) %>%
    raster(.)
} else {
  postProcess(preSimListbaseCase$standAgeMap,
              studyArea = preSimListbaseCase$studyAreaLarge,
              useCache = FALSE)  ## can have ages beyond study area from fire polys
}

rstLCCSetA <- if (!inMemory(preSimListbaseCase$rstLCC)) {
  list.files("R/SpaDES/inputs/",
             pattern = basename(filename(preSimListbaseCase$rstLCC)),
             recursive = TRUE, full.names = TRUE) %>%
    raster(.)
} else {
  preSimListbaseCase$rstLCC
}


## input stand biomass and age, ecological zonation (ecodistricts) and land-cover (LCC 2005)
rawBiomassMapSetB <- if (!inMemory(preSimListstudyAreaChange$rawBiomassMap)) {
  list.files("R/SpaDES/inputs/",
             pattern = basename(filename(preSimListstudyAreaChange$rawBiomassMap)),
             recursive = TRUE, full.names = TRUE) %>%
    raster(.)
} else {
  preSimListstudyAreaChange$rawBiomassMap
}

standAgeMapSetB <- if (!inMemory(preSimListstudyAreaChange$standAgeMap)) {
  list.files("R/SpaDES/inputs/",
             pattern = basename(filename(preSimListstudyAreaChange$standAgeMap)),
             recursive = TRUE, full.names = TRUE) %>%
    raster(.)
} else {
  postProcess(preSimListstudyAreaChange$standAgeMap,
              studyArea = preSimListstudyAreaChange$studyAreaLarge,
              useCache = FALSE)  ## can have ages beyond study area from fire polys
}

rstLCCSetB <- if (!inMemory(preSimListstudyAreaChange$rstLCC)) {
  list.files("R/SpaDES/inputs/",
             pattern = basename(filename(preSimListstudyAreaChange$rstLCC)),
             recursive = TRUE, full.names = TRUE) %>%
    raster(.)
} else {
  preSimListstudyAreaChange$rstLCC
}

ecoregionLayerSetA <- preSimListbaseCase$ecoregionLayer
ecoregionLayerSetB <- preSimListstudyAreaChange$ecoregionLayer

## all all possible ecodistricts so that legend agrees between plots
allVals <- unique(c(ecoregionLayerSetA$ECODISTRIC, ecoregionLayerSetB$ECODISTRIC))
ecoregionLayerSetA$ECODISTRIC <- factor(ecoregionLayerSetA$ECODISTRIC, allVals)
ecoregionLayerSetB$ECODISTRIC <- factor(ecoregionLayerSetB$ECODISTRIC, allVals)

plotEcodistSetA <- ggplot() +
  layer_spatial(ecoregionLayerSetA, aes(fill = ECODISTRIC)) +
  scale_fill_brewer(type = "qual", drop = FALSE) +
  labs(fill = "", title = "set A") +
  theme_pubr(margin = FALSE)

lims <- range(getValues(rawBiomassMapSetA), getValues(rawBiomassMapSetB), na.rm = TRUE)
plotrawBiomassMapSetA <- ggplot() +
  layer_spatial(rawBiomassMapSetA) +
  scale_fill_distiller(palette = "YlGn", na.value = "transparent", direction = 1,
                       limits = lims) +
  labs(fill = "", title = "set A") +
  theme_pubr(margin = FALSE)

lims <- range(getValues(standAgeMapSetA), getValues(standAgeMapSetB), na.rm = TRUE)
plotstandAgeMapSetA <- ggplot() +
  layer_spatial(standAgeMapSetA) +
  scale_fill_distiller(palette = "Blues", na.value = "transparent", direction = 1,
                       limits = lims) +
  labs(fill = "", title = "set A") +
  theme_pubr(margin = FALSE)

vals <- c(getValues(rstLCCSetA), getValues(rstLCCSetB), na.rm = TRUE)
plotrstLCCSetA <- ggplot() +
  layer_spatial(rstLCCSetA) +
  scale_fill_distiller(palette = "Paired", na.value = "transparent",
                       guide = "legend", breaks = sort(unique(vals)), limits = range(vals)) +
  labs(fill = "", title = "set A") +
  theme_pubr(margin = FALSE)

plotEcodistSetB <- ggplot() +
  layer_spatial(ecoregionLayerSetB, aes(fill = ECODISTRIC)) +
  scale_fill_brewer(type = "qual", drop = FALSE) +
  labs(fill = "", title = "set B") +
  theme_pubr(margin = FALSE)

lims <- range(getValues(rawBiomassMapSetA), getValues(rawBiomassMapSetB), na.rm = TRUE)
plotrawBiomassMapSetB <- ggplot() +
  layer_spatial(rawBiomassMapSetB) +
  scale_fill_distiller(palette = "YlGn", na.value = "transparent", direction = 1,
                       limits = lims) +
  labs(fill = "", title = "set B") +
  theme_pubr(margin = FALSE)

lims <- range(getValues(standAgeMapSetA), getValues(standAgeMapSetB), na.rm = TRUE)
plotstandAgeMapSetB <- ggplot() +
  layer_spatial(standAgeMapSetB) +
  scale_fill_distiller(palette = "Blues", na.value = "transparent", direction = 1,
                       limits = lims) +
  labs(fill = "", title = "set B") +
  theme_pubr(margin = FALSE)

vals <- c(getValues(rstLCCSetA), getValues(rstLCCSetB), na.rm = TRUE)
plotrstLCCSetB <- ggplot() +
  layer_spatial(rstLCCSetB) +
  scale_fill_distiller(palette = "Paired", na.value = "transparent",
                       guide = "legend",
                       breaks = sort(unique(vals)), limits = range(vals)) +
  labs(fill = "", title = "set B") +
  theme_pubr(margin = FALSE)

biomassPlots <- ggarrange(plotrawBiomassMapSetA,
                          plotrawBiomassMapSetB,
                          ncol = 2, nrow = 1,
                          # labels = c("NFI kNN stand biomass", ""),
                          legend = "right",
                          common.legend = TRUE)
agePlots <- ggarrange(plotstandAgeMapSetA,
                      plotstandAgeMapSetB,
                      ncol = 2, nrow = 1,
                      # labels = c("NFI kNN stand age", ""),
                      legend = "right",
                      common.legend = TRUE)
lccPlots <- ggarrange(plotrstLCCSetA,
                      plotrstLCCSetB,
                      ncol = 2, nrow = 1,
                      # labels = c("land-cover", ""),
                      legend = "right",
                      common.legend = TRUE)
ecodistPlots <- ggarrange(plotEcodistSetA,
                          plotEcodistSetB,
                          ncol = 2, nrow = 1,
                          # labels = c("ecodistricts", ""),
                          legend = "right",
                          common.legend = TRUE)
allPlots <- ggarrange(biomassPlots,
                      agePlots,
                      lccPlots,
                      ecodistPlots, hjust = 0, align = "h",
                      labels = c("NFI kNN stand biomass",
                                 "NFI kNN stand age",
                                 "land-cover",
                                 "ecodistricts"))
ggsave(plot = allPlots, filename = file.path(figDir, "inputMaps.png"),
       width = 14, height = 10, units = "in", dpi = 300)


## species % cover
speciesLayersSetA <- if (!inMemory(preSimListbaseCase$speciesLayers)) {
    stack(file.path("R/SpaDES/outputs/jun2021Runs/baseCase",
                    basename(filename(preSimListbaseCase$speciesLayers[[1]]))))
} else {
  preSimListbaseCase$speciesLayers
}

speciesLayersSetB <- if (!inMemory(preSimListstudyAreaChange$speciesLayers)) {
  stack(file.path("R/SpaDES/outputs/jun2021Runs/studyAreaChange",
                  basename(filename(preSimListstudyAreaChange$speciesLayers[[1]]))))
} else {
  preSimListstudyAreaChange$speciesLayers
}

x11(width = 10, height = 15, pointsize = 15)
clearPlot()
Plot(speciesLayersSetA,
     speciesLayersSetB,
     na.color = "white", zero.color = "white",
     title = c(paste(names(speciesLayersSetA), "- set A"),
               paste(names(speciesLayersSetB), "- set B")),
     new = TRUE)
savePlot(file.path(figDir, "speciesLayers.tiff"), type =  "tiff")
dev.off()



## running times - appendix 1
runTimeMapSize <- prepInputs(url = "https://drive.google.com/file/d/1QMDrpURYUhVsTzoDn0K7RlIpxdcyjvnS/view?usp=sharing",
                             fun = "read.csv",
                             useCache = FALSE)
runTimeMapSize <- as.data.table(runTimeMapSize)

plotData <- melt(runTimeMapSize, measure.vars = c("Time_Landis_mean", "Time_SpaDES_mean", "Faster_Mag"),
                 variable.name = "model", value.name = "time" )

coef <- max(max(runTimeMapSize$Time_Landis_mean / runTimeMapSize$Faster_Mag),
            max(runTimeMapSize$Time_SpaDES_mean / runTimeMapSize$Faster_Mag))

runTimePlot1 <- ggplot(plotData[model != "Faster_Mag"], aes(x = pixels, y = time, col = model)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_line(data = plotData[model == "Faster_Mag"], aes(y = time * coef), size = 1) +
  geom_point(data = plotData[model == "Faster_Mag"], aes(y = time * coef), size = 2) +
  scale_y_continuous(sec.axis = sec_axis(~ . / coef, name = "Running time ratio\n(LBSE/Biomass_core)"),
                     limits = c(0, max(plotData$time) + 5000)) +
  scale_colour_manual(values = c("Time_Landis_mean" = "green",
                                 "Time_SpaDES_mean" = "red",
                                 "Faster_Mag" = "blue"),
                      labels = c("Time_Landis_mean" = "LBSE",
                                 "Time_SpaDES_mean" = "LandR Biomass_core",
                                 "Faster_Mag" = "ratio"),
                      breaks = c("Time_SpaDES_mean", "Time_Landis_mean", "Faster_Mag")) +
  theme_pubr(border = TRUE, base_size = 12, legend = c(0.25, 0.90)) +
  theme(legend.background = element_blank(), axis.title.y.right = element_text(angle = 90)) +
  labs(x = "No. of pixels", y = "Running time\n(seconds)",
       col = "")

runTimePlot2 <- ggplot(plotData[model != "Faster_Mag"],
                       aes(x = pixels, y = time/pixels * 1000, col = model)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("Time_Landis_mean" = "green",
                                 "Time_SpaDES_mean" = "red",
                                 "Faster_Mag" = "blue"),
                      labels = c("Time_Landis_mean" = "LBSE",
                                 "Time_SpaDES_mean" = "LandR Biomass_core",
                                 "Faster_Mag" = "ratio"),
                      breaks = c("Time_SpaDES_mean", "Time_Landis_mean", "Faster_Mag")) +
  theme_pubr(border = TRUE, base_size = 12, legend = c(0.2, 0.95)) +
  theme(legend.background = element_blank()) +
  labs(x = "No. of pixels", y = "Running time\nper 1,000 pixels (seconds)",
       col = "")

savePlot <- ggarrange(runTimePlot1, runTimePlot2 + theme(legend.position = "none"),
                      ncol = 1, align = "v", labels = "auto")
ggsave(plot = savePlot, filename = file.path(figDir, "runTimePlots.png"), device = "png",
       width = 6, height = 8, dpi = 300)
