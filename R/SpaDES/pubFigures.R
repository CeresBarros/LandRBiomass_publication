## ------------------------------------------------------
## LandR Biomass Publication - additional figures
##
## Ceres: July 2021
## ------------------------------------------------------

## this script should be sourced
Require::Require(c("SpaDES", "raster", "sf", "data.table", "ggplot2",
                   "ggspatial", "ggpubr", "reproducible", "cowplot", "qs"),
                 upgrade = FALSE)

## GET ALL SIMLISTS ------------------------------------
runNames <- c("baseCase", "studyAreaChange", "altParameters")

## species layers simLists
for (runName in runNames) {
  eval(parse(text = paste0("speciesLayersSimList", runName, " <- loadSimList(file.path('R/SpaDES/outputs',",
                           "simDirName, runName, paste0('simList_speciesLayers', runName)))")))
}

## pre simulation simLists
for (runName in runNames) {
  eval(parse(text = paste0("preSimList", runName, " <- loadSimList(file.path('R/SpaDES/outputs',",
                           "simDirName, runName, paste0('simInit', runName)))")))
}

## post-simulation simLists
for (runName in runNames) {
  eval(parse(text = paste0("simList", runName, " <- qs::qread(file.path('R/SpaDES/outputs',",
                           "simDirName, runName, paste0('simList_LandRBiomass_sim_', runName)))")))
}

## post-validation simLists
for (runName in runNames) {
  eval(parse(text = paste0("validSimList", runName, " <- qs::qread(file.path('R/SpaDES/validation', simDirName, runName, paste0('simValid', runName)))")))
}

## ELAPSED TIMES ----------------------------------
elapsedTime(speciesLayersSimListbaseCase)
elapsedTime(speciesLayersSimListstudyAreaChange)
elapsedTime(speciesLayersSimListaltParameters)

elapsedTime(preSimListbaseCase)
elapsedTime(preSimListstudyAreaChange)
elapsedTime(preSimListaltParameters)

elapsedTime(validSimListbaseCase)[, sum(as.numeric(elapsedTime))]
elapsedTime(validSimListstudyAreaChange)[, sum(as.numeric(elapsedTime))]
elapsedTime(validSimListaltParameters)[, sum(as.numeric(elapsedTime))]

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

## INSPECT MODEL BIOMASS AND MODEL COVER ----------------------------------
preSimListbaseCase$modelBiomass$rsq
preSimListbaseCase$modelCover$rsq

preSimListstudyAreaChange$modelBiomass$rsq
preSimListstudyAreaChange$modelCover$rsq

preSimListaltParameters$modelBiomass$rsq
preSimListaltParameters$modelCover$rsq

## TABLE OF SNLL DIFFERENCES ----------------------------------
## alternative parameterisation minus base case
## positive values mean alternative parameterisation is worse.
SNLL <- rbind(validSimListbaseCase$logLikelihood, validSimListaltParameters$logLikelihood,
              use.names = TRUE, idcol = "simulation")
SNLL[, simulation := ifelse(simulation == 1, "basecase", "altParameters")]
SNLL[, variable2 := paste(SNLL[[2]], SNLL[[3]], sep = "_")]
SNLL[, variable :=  NULL]; SNLL[, variable :=  NULL]

setnames(SNLL, old = c("2001", "2011", "2011-2001"), new = c("y2001", "y2011", "y2011_2001"))
SNLLdif <- SNLL[, list("delta2001" = y2001[simulation == "altParameters"] - y2001[simulation == "basecase"],
                       "delta2011" = y2011[simulation == "altParameters"] - y2011[simulation == "basecase"],
                       "delta2011-2001" = y2011_2001[simulation == "altParameters"] - y2011_2001[simulation == "basecase"]),
                by = variable2]

## "ENSEMBLE" OF TWO INITIAL CONDITIONS ----------------------------------
## get data:
cohortDataOutputs <- rbindlist(lapply(simListbaseCase, FUN = function(sim) as.data.table(outputs(sim))[objectName == "cohortData"]),
                               idcol = "rep")
cohortDataOutputs[, rep :=  sub("sim1_rep", "", rep)]
baseCase_allCohortData <- rbindlist(fill = TRUE, use.names = TRUE,
                                    l = apply(cohortDataOutputs, MARGIN = 1, FUN = function(x) {
                                      cohortData <- readRDS(x["file"])
                                      cohortData[, year := as.numeric(x["saveTime"])]
                                      cohortData[, rep := as.numeric(x["rep"])]
                                      return(cohortData)
                                    }))

cohortDataOutputs <- rbindlist(lapply(simListaltParameters, FUN = function(sim) as.data.table(outputs(sim))[objectName == "cohortData"]),
                               idcol = "rep")
cohortDataOutputs[, rep :=  sub("sim1_rep", "", rep)]
altParams_allCohortData <- rbindlist(fill = TRUE, use.names = TRUE,
                                     l = apply(cohortDataOutputs, MARGIN = 1, FUN = function(x) {
                                       cohortData <- readRDS(x["file"])
                                       cohortData[, year := as.numeric(x["saveTime"])]
                                       cohortData[, rep := as.numeric(x["rep"])]
                                       return(cohortData)
                                     }))
allCohortData <- rbindlist(list("base case" = baseCase_allCohortData,
                                "Example 2" = altParams_allCohortData), idcol = "scenario")


plotData <- allCohortData[, list(B = sum(B)), by = .(speciesCode, year, rep, scenario)]
plotData <- plotData[, list(B = mean(B)), by = .(speciesCode, year, scenario)]

sppColoursVect <- preSimListbaseCase$sppColorVect
sppLabels <- equivalentName(names(sppColoursVect), preSimListbaseCase$sppEquiv, column = "EN_generic_short")
names(sppLabels) <- names(sppColoursVect)
sppLabels <- sppLabels[!is.na(sppLabels)]

sppColoursVect <- sppColoursVect[names(sppLabels)]

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
  labs(y = expression(Total~~biomass~~(g/m^2)), x = "Year", colour = "", fill = "")
ggsave(file.path(figDir, "baseCase_altParameters_ensemble.png"), width = 7, height = 5, units = "in", dpi = 300)


## GENERAL FIGURES ----------------------------------------
## Study areas
## get Canadian provinces and subset to SK
can1 <- raster::getData('GADM', country = "CAN", level = 1, path = tempdir())
SK <- can1[can1$NAME_1 == "Saskatchewan",]
can1 <- spTransform(can1, CRSobj = crs(preSimListbaseCase$studyArea))
SK <- spTransform(SK, CRSobj = crs(preSimListbaseCase$studyArea))
can1 <- st_as_sf(can1)
SK <- st_as_sf(SK)
insetBox <- st_as_sfc(st_bbox(SK))

plotStudyAreasCanada <- ggplot() +
  layer_spatial(can1, fill = "grey98") +
  layer_spatial(insetBox, color = "red", fill = NA) +
  labs(fill = "") +
  theme_void()

plotStudyAreasSK <- ggplot() +
  layer_spatial(can1, fill = "grey90") +
  layer_spatial(preSimListbaseCase$studyAreaLarge, fill = "darkgreen") +
  layer_spatial(preSimListbaseCase$studyArea, fill = "lightgreen") +
  layer_spatial(preSimListstudyAreaChange$studyAreaLarge, fill = "darkblue") +
  layer_spatial(preSimListstudyAreaChange$studyArea, fill = "lightblue") +
  coord_sf(xlim = st_bbox(SK)[c(1, 3)],
           ylim = st_bbox(SK)[c(2, 4)]) +
  labs(fill = "", x = "Longitude", y = "Latitude") +
  theme_pubr()

insetMap <- ggdraw() +
  draw_plot(plotStudyAreasSK) +
  draw_plot(plotStudyAreasCanada, x = 0.60, y = 0.65, width = 0.4, height = 0.4)
ggsave(plot = insetMap, filename = file.path(figDir, "studyAreas_setsAB.png"), bg = "white",
       width = 6, height = 7, units = "in", dpi = 300)


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
  theme_pubr(margin = FALSE)

lims <- range(getValues(rawBiomassMapSetA), getValues(rawBiomassMapSetB), na.rm = TRUE)
plotrawBiomassMapSetA <- ggplot() +
  layer_spatial(rawBiomassMapSetA) +
  scale_fill_distiller(palette = "YlGn", na.value = "transparent", direction = 1,
                       limits = lims) +
  theme_pubr(margin = FALSE)

lims <- range(getValues(standAgeMapSetA), getValues(standAgeMapSetB), na.rm = TRUE)
plotstandAgeMapSetA <- ggplot() +
  layer_spatial(standAgeMapSetA) +
  scale_fill_distiller(palette = "Blues", na.value = "transparent", direction = 1,
                       limits = lims) +
  theme_pubr(margin = FALSE)

vals <- c(getValues(rstLCCSetA), getValues(rstLCCSetB), na.rm = TRUE)
plotrstLCCSetA <- ggplot() +
  layer_spatial(rstLCCSetA) +
  scale_fill_distiller(palette = "Paired", na.value = "transparent",
                       guide = "legend", breaks = sort(unique(vals)), limits = range(vals)) +
  theme_pubr(margin = FALSE)

plotEcodistSetB <- ggplot() +
  layer_spatial(ecoregionLayerSetB, aes(fill = ECODISTRIC)) +
  scale_fill_brewer(type = "qual", drop = FALSE) +
  theme_pubr(margin = FALSE)

lims <- range(getValues(rawBiomassMapSetA), getValues(rawBiomassMapSetB), na.rm = TRUE)
plotrawBiomassMapSetB <- ggplot() +
  layer_spatial(rawBiomassMapSetB) +
  scale_fill_distiller(palette = "YlGn", na.value = "transparent", direction = 1,
                       limits = lims) +
  theme_pubr(margin = FALSE)

lims <- range(getValues(standAgeMapSetA), getValues(standAgeMapSetB), na.rm = TRUE)
plotstandAgeMapSetB <- ggplot() +
  layer_spatial(standAgeMapSetB) +
  scale_fill_distiller(palette = "Blues", na.value = "transparent", direction = 1,
                       limits = lims) +
  theme_pubr(margin = FALSE)

vals <- c(getValues(rstLCCSetA), getValues(rstLCCSetB), na.rm = TRUE)
plotrstLCCSetB <- ggplot() +
  layer_spatial(rstLCCSetB) +
  scale_fill_distiller(palette = "Paired", na.value = "transparent",
                       guide = "legend",
                       breaks = sort(unique(vals)), limits = range(vals)) +
  theme_pubr(margin = FALSE)

biomassPlots <- ggarrange(plotrawBiomassMapSetA + labs(y = "Latitude", x = "", fill = "", title = "set A"),
                          plotrawBiomassMapSetB + labs(y = "", x = "", fill = "", title = "set B"),
                          ncol = 2, nrow = 1,
                          # labels = c("NFI kNN stand biomass", ""),
                          legend = "right",
                          common.legend = TRUE)
agePlots <- ggarrange(plotstandAgeMapSetA + labs(y = "", x = "", fill = "", title = "set A"),
                      plotstandAgeMapSetB + labs(y = "", x = "", fill = "", title = "set B"),
                      ncol = 2, nrow = 1,
                      # labels = c("NFI kNN stand age", ""),
                      legend = "right",
                      common.legend = TRUE)
lccPlots <- ggarrange(plotrstLCCSetA + labs(y = "Latitude", x = "Longitude", fill = "", title = "set A"),
                      plotrstLCCSetB + labs(y = "", x = "Longitude", fill = "", title = "set B"),
                      ncol = 2, nrow = 1,
                      # labels = c("land-cover", ""),
                      legend = "right",
                      common.legend = TRUE)
ecodistPlots <- ggarrange(plotEcodistSetA + labs(y = "", x = "Longitude", fill = "", title = "set A"),
                          plotEcodistSetB + labs(y = "", x = "Longitude", fill = "", title = "set B"),
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
       width = 14, height = 10, units = "in", dpi = 300, bg = "white")


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


## RUN TIMES ----------------------------------------
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

## SIMULATION FIGURES ----------------------------------------
## Full suite of simulation plots from study area change
speciesB <- qread(file.path(outputPath(simListstudyAreaChange$sim1_rep01), "figures", "biomass_by_species_gg.qs"))
speciesBoverstory <- qread(file.path(outputPath(simListstudyAreaChange$sim1_rep01), "figures", "overstory_biomass_gg.qs"))
speciesAge <- qread(file.path(outputPath(simListstudyAreaChange$sim1_rep01), "figures", "biomass-weighted_species_age_gg.qs"))
speciesDom <- qread(file.path(outputPath(simListstudyAreaChange$sim1_rep01), "figures", "N_pixels_leading_gg.qs"))
speciesANPP <- qread(file.path(outputPath(simListstudyAreaChange$sim1_rep01), "figures", "total_aNPP_by_species_gg.qs"))
landscapeVars <- qread(file.path(outputPath(simListstudyAreaChange$sim1_rep01), "figures", "landscape_biomass_aNPP_max_age_gg.qs"))

speciesB <- speciesB + labs(title = "Total biomass", y = expression(g/m^2), x = "Year", fill = "")
speciesBoverstory <- speciesBoverstory + labs(title = "Overstory biomass", y = expression(g/m^2), x = "Year", fill = "")
speciesAge <- speciesAge + labs(title = "Species age", y = "Years", x = "Year", colour = "")
speciesDom <- speciesDom + labs(title = "Species dominance", y = "No. of pixels", x = "Year", fill = "")
speciesANPP <- speciesANPP + labs(title = "Total aNPP", y = expression(g/m^2), x = "Year", colour = "")


landscapePlots <- align_plots(speciesB + theme(legend.position = "none"),
                              speciesDom + theme(legend.position = "none"),
                              landscapeVars,
                              align = 'v', axis = 'l')
speciesPlots <- plot_grid(landscapePlots[[1]],
                          speciesBoverstory + theme(legend.position = "none"),
                          speciesAge + theme(legend.position = "none"),
                          landscapePlots[[2]],
                          speciesANPP + theme(legend.position = "none"),
                          get_legend(speciesB + theme(legend.direction = "horizontal")))
allPlotsSim <- plot_grid(speciesPlots, landscapePlots[[3]], ncol = 1,
                         labels = c("a)", "b)"),
                         rel_heights = c(0.7, 0.3))
ggsave(plot = allPlotsSim, filename = file.path(figDir, "simulationPlots_studyAreaChange.png"),
       bg = "white", width = 14, height = 10, dpi = 300)


## VALIDATION FIGURES ----------------------------------------
## Full suite of validation plots from study area change
## in this case we need the "data" with is the ggplot object
landscapeRelB <- qread(file.path(outputPath(validSimListstudyAreaChange), "figures", "LandscapeComparisons_relB_data.qs"))
landscapeDom <- qread(file.path(outputPath(validSimListstudyAreaChange), "figures", "LandscapeComparisons_Dom_data.qs"))
landscapePresAbs <- qread(file.path(outputPath(validSimListstudyAreaChange), "figures", "LandscapeComparisons_PresAbs_data.qs"))
landscapeDeltaB <- qread(file.path(outputPath(validSimListstudyAreaChange), "figures", "LandscapeComparisons_deltaB_data.qs"))

pixelRelB <- qread(file.path(outputPath(validSimListstudyAreaChange), "figures", "PixelComparisons_relB_data.qs"))
pixelDeltaB <- qread(file.path(outputPath(validSimListstudyAreaChange), "figures", "PixelComparisons_deltaB_data.qs"))

## for MAD its the GG output
landscapeMAD <- qread(file.path(outputPath(validSimListstudyAreaChange), "figures", "landscapeMAD_gg.qs"))
pixelMAD <- qread(file.path(outputPath(validSimListstudyAreaChange), "figures", "pixelMAD_gg.qs"))


landscapeRelB <- landscapeRelB + labs(title = "Relative abundance")
landscapePresAbs <- landscapePresAbs + labs(title = "Presences", y = "No. of pixels")
landscapeDom <- landscapeDom + labs(title = "Dominance", y = "No. of pixels")
landscapeDeltaB <- landscapeDeltaB + labs(title = "Changes in biomass")

pixelRelB <- pixelRelB + labs(title = "Relative abundance")
pixelDeltaB <- pixelDeltaB + labs(title = "Changes in biomass")

landscapeMAD <- landscapeMAD + labs(title = "Landscape-level")
pixelMAD <- pixelMAD + labs(title = "Pixel-level")

allPlotsValid <- plot_grid(
  plot_grid(
    plot_grid(landscapeRelB + theme(legend.position = "none"),
              landscapePresAbs + theme(legend.position = "none"),
              landscapeDom + theme(legend.position = "none"),
              landscapeDeltaB + theme(legend.position = "none"),
              align = "v", axis = "l", ncol = 2),
    get_legend(landscapeRelB), rel_heights = c(1, 0.1), ncol = 1),
  plot_grid(pixelRelB + theme(legend.position = "none"),
            pixelDeltaB + theme(legend.position = "none"),
            get_legend(pixelDeltaB), ncol = 1, rel_heights = c(1, 1, 0.1),
            align = "v", axis = "l"),
  plot_grid(landscapeMAD + theme(legend.position = "none"),
            pixelMAD + theme(legend.position = "none"),
            get_legend(landscapeMAD), ncol = 1, rel_heights = c(1, 1, 0.1),
            align = "v", axis = "l"),
  ncol = 3, rel_widths = c(1.2, 0.6, 0.6), labels = c("a)", "b)", "c)"))

ggsave(plot = allPlotsValid, filename = file.path(figDir, "validationPlots_studyAreaChange.png"),
       bg = "white", width = 14, height = 8, dpi = 300)


## plots confronting Example 2 simulations
landscapeRelB_BC <- qread(file.path(outputPath(validSimListbaseCase), "figures", "LandscapeComparisons_relB_data.qs"))
landscapeDom_BC <- qread(file.path(outputPath(validSimListbaseCase), "figures", "LandscapeComparisons_Dom_data.qs"))
landscapeDeltaB_BC <- qread(file.path(outputPath(validSimListbaseCase), "figures", "LandscapeComparisons_deltaB_data.qs"))

pixelRelB_BC <- qread(file.path(outputPath(validSimListbaseCase), "figures", "PixelComparisons_relB_data.qs"))
pixelDeltaB_BC <- qread(file.path(outputPath(validSimListbaseCase), "figures", "PixelComparisons_deltaB_data.qs"))

## for MAD its the GG output
landscapeMAD_BC <- qread(file.path(outputPath(validSimListbaseCase), "figures", "landscapeMAD_gg.qs"))
pixelMAD_BC <- qread(file.path(outputPath(validSimListbaseCase), "figures", "pixelMAD_gg.qs"))

## plots confronting Example 2 simulations
landscapeRelB_AP <- qread(file.path(outputPath(validSimListaltParameters), "figures", "LandscapeComparisons_relB_data.qs"))
landscapeDom_AP <- qread(file.path(outputPath(validSimListaltParameters), "figures", "LandscapeComparisons_Dom_data.qs"))
landscapeDeltaB_AP <- qread(file.path(outputPath(validSimListaltParameters), "figures", "LandscapeComparisons_deltaB_data.qs"))

pixelRelB_AP <- qread(file.path(outputPath(validSimListaltParameters), "figures", "PixelComparisons_relB_data.qs"))
pixelDeltaB_AP <- qread(file.path(outputPath(validSimListaltParameters), "figures", "PixelComparisons_deltaB_data.qs"))

## for MAD its the GG output
landscapeMAD_AP <- qread(file.path(outputPath(validSimListaltParameters), "figures", "landscapeMAD_gg.qs"))
pixelMAD_AP <- qread(file.path(outputPath(validSimListaltParameters), "figures", "pixelMAD_gg.qs"))




