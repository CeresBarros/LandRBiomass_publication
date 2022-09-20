## ------------------------------------------------------
## LandR Biomass Publication - additional figures
##
## Ceres: July 2021
## ------------------------------------------------------

## this script should be sourced

## GET ALL SIMLISTS ------------------------------------
runNames <- c("baseCase", "studyAreaChange", "altParameters")

## species layers simLists
for (runName in runNames) {
  eval(parse(text = paste0("speciesLayersSimList", runName, " <- loadSimList(file.path('R/SpaDES/outputs',",
                           "simDirName, runName, paste0('simList_speciesLayers', runName, '.qs')))")))
}

## pre simulation simLists
for (runName in runNames) {
  eval(parse(text = paste0("preSimList", runName, " <- loadSimList(file.path('R/SpaDES/outputs',",
                           "simDirName, runName, paste0('simInit', runName, '.qs')))")))
}

## post-simulation simLists
for (runName in runNames) {
  eval(parse(text = paste0("simList", runName, " <- qs::qread(file.path('R/SpaDES/outputs',",
                           "simDirName, runName, paste0('simList_LandRBiomass_sim_', runName, '.qs')))")))
}

## post-validation simLists
for (runName in runNames) {
  eval(parse(text = paste0("validSimList", runName, " <- qs::qread(file.path('R/SpaDES/validation', simDirName, runName, paste0('simValid', runName, '.qs')))")))
}

plotTheme <- function(majorYlines = TRUE, ...) {
  .theme <- theme_pubr(...)
  if (majorYlines) {
    .theme <- .theme +
      theme(panel.grid.major.y = element_line(colour = "grey", linetype = "dotted", size = 0.5))
  }
  .theme
}

## ELAPSED TIMES ----------------------------------
## in minutes
## Biomass_speciesData
sum(as.numeric(elapsedTime(speciesLayersSimListbaseCase)$elapsedTime))/60
sum(as.numeric(elapsedTime(speciesLayersSimListstudyAreaChange)$elapsedTime))/60
sum(as.numeric(elapsedTime(speciesLayersSimListaltParameters)$elapsedTime))/60

## Biomass_borealDataPrep + speciesParams
elapsedTime(preSimListbaseCase)
elapsedTime(preSimListstudyAreaChange)
elapsedTime(preSimListaltParameters)

## Biomass_validationKNN
elapsedTime(validSimListbaseCase)[, sum(as.numeric(elapsedTime))]/60
elapsedTime(validSimListstudyAreaChange)[, sum(as.numeric(elapsedTime))]/60
elapsedTime(validSimListaltParameters)[, sum(as.numeric(elapsedTime))]/60

## for biomass_core use clocktimes as elapsed time seems to underestimate
timesbaseCase <- rbindlist(lapply(simListbaseCase, completed), idcol = "rep")
timesbaseCase <- timesbaseCase[eventTime %in% seq(start(preSimListbaseCase), end(preSimListbaseCase))] ## exclude the spinup (will have a very differet eventTime in this case)
timesbaseCase <- timesbaseCase[eventType != "init", max(clockTime) - min(clockTime), by = .(rep, moduleName)]
timesbaseCase[, sum(as.integer(V1)), by = .(rep, moduleName)]
timesbaseCase[, mean(as.integer(V1))/60, by = moduleName]

timesstudyAreaChange <- rbindlist(lapply(simListstudyAreaChange, completed), idcol = "rep")
timesstudyAreaChange <- timesstudyAreaChange[eventTime %in% seq(start(preSimListstudyAreaChange), end(preSimListstudyAreaChange))] ## exclude the spinup (will have a very differet eventTime in this case)
timesstudyAreaChange <- timesstudyAreaChange[eventType != "init", max(clockTime) - min(clockTime), by = .(rep, moduleName)]
timesstudyAreaChange[, sum(as.integer(V1)), by = .(rep, moduleName)]
timesstudyAreaChange[, mean(as.integer(V1))/60, by = moduleName]

timesaltParameters <- rbindlist(lapply(simListaltParameters, completed), idcol = "rep")
timesaltParameters <- timesaltParameters[eventTime %in% seq(start(preSimListaltParameters), end(preSimListaltParameters))] ## exclude the spinup (will have a very differet eventTime in this case)
timesaltParameters <- timesaltParameters[eventType != "init", max(clockTime) - min(clockTime), by = .(rep, moduleName)]
timesaltParameters[, sum(as.integer(V1)), by = .(rep, moduleName)]
timesaltParameters[, mean(as.integer(V1))/60, by = moduleName]

## OBJECT DIAGRAM
objectDiagram(preSimListbaseCase, width = 1000, height = 2500)
webshot("http://localhost:20581/session/viewhtml6aec1852565c/index.html",
        file = file.path(figDir, "objectDiagram_baseCase.png"))

## INSPECT MODEL BIOMASS AND MODEL COVER ----------------------------------
## two statistical models used to derive spatially-varying species traits (maxB, maxANPP and SEP)
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

## assume a penalisation of 80 for base case
AIC <- SNLL[simulation == "basecase", list(simulation = simulation,
                                           y2001 = 80 + 2*y2001,
                                           y2011 = 80 + 2*y2011,
                                           y2011_2001 = 80 + 2*y2011_2001,
                                           variable2 = variable2)]
AIC <- rbind(AIC, SNLL[simulation == "altParameters", list(simulation = simulation,
                                                           y2001 = 2*y2001,
                                                           y2011 = 2*y2011,
                                                           y2011_2001 = 2*y2011_2001,
                                                           variable2 = variable2)])
AICdiff <- AIC[, list("delta2001" = y2001[simulation == "altParameters"] - y2001[simulation == "basecase"],
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
# plotData <- plotData[, list(B = mean(B)), by = .(speciesCode, year, scenario)]

sppColoursVect <- preSimListbaseCase$sppColorVect
sppLabels <- equivalentName(names(sppColoursVect), preSimListbaseCase$sppEquiv, column = "EN_generic_short")
names(sppLabels) <- names(sppColoursVect)
sppLabels <- sppLabels[!is.na(sppLabels)]

sppColoursVect <- sppColoursVect[names(sppLabels)]

ggplot() +
  ## mean across reps per scenario
  stat_summary(data = plotData[scenario == "Example 2"],
               aes(x = year, y = B, colour = speciesCode, fill = speciesCode),
               fun.data = mean_sdl, geom = "ribbon",
               alpha = 0.5, colour = "transparent") +
  stat_summary(data = plotData[scenario == "Example 2"],
               aes(x = year, y = B, colour = speciesCode,
                   linetype = "alternative param."),
               fun = mean, geom = "line", size = 1, alpha = 0.8) +
  stat_summary(data = plotData[scenario == "base case"],
               aes(x = year, y = B, colour = speciesCode, fill = speciesCode),
               fun.data = mean_sdl, geom = "ribbon",
               alpha = 0.5, colour = "transparent") +
  stat_summary(data = plotData[scenario == "base case"],
               aes(x = year, y = B, colour = speciesCode,
                   linetype = "base case"),
               fun = mean, geom = "line", size = 1, alpha = 0.8) +
  ## mean across scenarios and reps (ensemble)
  # stat_summary(data = plotData,
  #              aes(x = year, y = B, colour = speciesCode, fill = speciesCode),
  #              fun.data = mean_sdl, geom = "ribbon",
  #              alpha = 0.5, colour = "transparent") +
  stat_summary(data = plotData,
               aes(x = year, y = B, colour = speciesCode,
                   linetype = "ensemble"),
               fun = mean, geom = "line", size = 1) +
  scale_colour_manual(values = sppColoursVect,
                      labels = sppLabels) +
  scale_fill_manual(values = sppColoursVect,
                    labels = sppLabels) +
  scale_linetype_manual(values = c("ensemble" = 1, "base case" = 2, "alternative param." = 3)) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_pubr(base_size = 12, legend = "right") +
  labs(y = expression(Average~~total~~biomass~~(g/m^2)), x = "Year", colour = "",
       fill = "", linetype = "")
ggsave(file.path(figDir, "baseCase_altParameters_ensemble.png"),
       width = 8, height = 4.5, units = "in", dpi = 300)


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

## Simulation plots for example 2
## base case
speciesB_BC <- qread(file.path(outputPath(simListbaseCase$sim1_rep01), "figures", "biomass_by_species_gg.qs"))
speciesBoverstory_BC <- qread(file.path(outputPath(simListbaseCase$sim1_rep01), "figures", "overstory_biomass_gg.qs"))
speciesAge_BC <- qread(file.path(outputPath(simListbaseCase$sim1_rep01), "figures", "biomass-weighted_species_age_gg.qs"))
speciesDom_BC <- qread(file.path(outputPath(simListbaseCase$sim1_rep01), "figures", "N_pixels_leading_gg.qs"))
speciesANPP_BC <- qread(file.path(outputPath(simListbaseCase$sim1_rep01), "figures", "total_aNPP_by_species_gg.qs"))

speciesB_BC <- speciesB_BC + labs(title = "Total biomass", y = expression(g/m^2), x = "Year", fill = "")
speciesBoverstory_BC <- speciesBoverstory_BC + labs(title = "Overstory biomass", y = expression(g/m^2), x = "Year", fill = "")
speciesAge_BC <- speciesAge_BC + labs(title = "Species age", y = "Years", x = "Year", colour = "")
speciesDom_BC <- speciesDom_BC + labs(title = "Species dominance", y = "No. of pixels", x = "Year", fill = "")
speciesANPP_BC <- speciesANPP_BC + labs(title = "Total aNPP", y = expression(g/m^2), x = "Year", colour = "")

## alternative params.
speciesB_AP <- qread(file.path(outputPath(simListaltParameters$sim1_rep01), "figures", "biomass_by_species_gg.qs"))
speciesBoverstory_AP <- qread(file.path(outputPath(simListaltParameters$sim1_rep01), "figures", "overstory_biomass_gg.qs"))
speciesAge_AP <- qread(file.path(outputPath(simListaltParameters$sim1_rep01), "figures", "biomass-weighted_species_age_gg.qs"))
speciesDom_AP <- qread(file.path(outputPath(simListaltParameters$sim1_rep01), "figures", "N_pixels_leading_gg.qs"))
speciesANPP_AP <- qread(file.path(outputPath(simListaltParameters$sim1_rep01), "figures", "total_aNPP_by_species_gg.qs"))

speciesB_AP <- speciesB_AP + labs(title = "Total biomass", y = expression(g/m^2), x = "Year", fill = "")
speciesBoverstory_AP <- speciesBoverstory_AP + labs(title = "Overstory biomass", y = expression(g/m^2), x = "Year", fill = "")
speciesAge_AP <- speciesAge_AP + labs(title = "Species age", y = "Years", x = "Year", colour = "")
speciesDom_AP <- speciesDom_AP + labs(title = "Species dominance", y = "No. of pixels", x = "Year", fill = "")
speciesANPP_AP <- speciesANPP_AP + labs(title = "Total aNPP", y = expression(g/m^2), x = "Year", colour = "")


allPlotsSim_BC_AP <- plot_grid(
  plot_grid(
    plot_grid(speciesBoverstory_BC + theme(legend.position = "none", text = element_text(size = 10)),
              speciesDom_BC + theme(legend.position = "none", text = element_text(size = 10)),
              speciesAge_BC + theme(legend.position = "none", text = element_text(size = 10)),
              speciesANPP_BC + theme(legend.position = "none", text = element_text(size = 10)),
              align = 'v', axis = 'l', ncol = 1),
    plot_grid(speciesBoverstory_AP + theme(legend.position = "none", text = element_text(size = 10)),
              speciesDom_AP + theme(legend.position = "none", text = element_text(size = 10)),
              speciesAge_AP + theme(legend.position = "none", text = element_text(size = 10)),
              speciesANPP_AP + theme(legend.position = "none", text = element_text(size = 10)),
              align = 'v', axis = 'l', ncol = 1),
    ncol = 2, labels = c("a)", "b)")),
  NULL,
  get_legend(speciesBoverstory_BC + theme(legend.direction = "horizontal", text = element_text(size = 10))),
  ncol = 1, rel_heights = c(1, 0, 0.1))
ggsave(plot = allPlotsSim_BC_AP, filename = file.path(figDir, "simulationPlots_example2.png"),
       bg = "white", height = 290, width = 180, dpi = 300, units = "mm")

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
landscapePresAbs <- landscapePresAbs + labs(title = "Presences")
landscapeDom <- landscapeDom + labs(title = "Dominance")
landscapeDeltaB <- landscapeDeltaB + labs(title = "Changes in biomass", y = expression(g/m^2))

pixelRelB <- pixelRelB + labs(title = "Relative abundance")
pixelDeltaB <- pixelDeltaB + labs(title = "Changes in biomass", y = expression(g/m^2))

## remove deltaB
landscapeMAD$data <- landscapeMAD$data[variable != "g/m^2"]
pixelMAD$data <- pixelMAD$data[variable != "g/m^2"]

landscapeMAD <- landscapeMAD + labs(title = "Mean abs. deviation", y = "") +
  scale_color_brewer(palette = "Dark2",
                     labels = landscapeMAD$plot_env$collabs,
                     drop = TRUE)
pixelMAD <- pixelMAD + labs(title = "Mean abs. deviation", y = "") +
  scale_color_brewer(palette = "Dark2",
                     labels = pixelMAD$plot_env$collabs,
                     drop = TRUE)



allPlotsValid <- plot_grid(
  plot_grid(landscapeRelB + theme(legend.position = "none", text = element_text(size = 10)),
            landscapeDom + theme(legend.position = "none", text = element_text(size = 10)),
            align = "v", axis = "l", ncol = 1),
  plot_grid(landscapePresAbs + theme(legend.position = "none", text = element_text(size = 10)),
            landscapeMAD + theme(legend.position = "none", text = element_text(size = 10)),
            align = "v", axis = "l", ncol = 1),
  NULL, NULL,
  get_legend(landscapeRelB), NULL,
  NULL, NULL,
  pixelRelB + theme(legend.position = "none", text = element_text(size = 10)),
  pixelMAD + theme(legend.position = "none", text = element_text(size = 10)),
  NULL, NULL,
  get_legend(pixelRelB),
  get_legend(landscapeMAD + guides(colour = guide_legend(ncol = 2))),
  labels = c("a)", "", "", "", "", "", "", "", "b)", "", "", "", "",""),
  rel_heights = c(1, -0.01, 0.05, 0.05, 0.5, -0.01, 0.1), ncol = 2,
)

ggsave(plot = allPlotsValid, filename = file.path(figDir, "validationPlots_studyAreaChange.png"),
       bg = "white", height = 250, width = 210, dpi = 300, units = "mm")


## plots confronting Example 2 simulations
## here we extract the plot data to make nicer plots
## Base case
landscapeVars_BC_AP <- rbindlist(list("aBC" = validSimListbaseCase$landscapeVars,
                                      "bAP" = validSimListaltParameters$landscapeVars),
                                 idcol = "simulation")
pixelVars_BC_AP <- rbindlist(list("aBC" = validSimListbaseCase$pixelVars,
                                  "bAP" = validSimListaltParameters$pixelVars),
                             idcol = "simulation")
pixelVars_BC_AP[, dataType2 := ifelse(dataType == "simulated", simulation, dataType)]

landscapeMAD_BC_AP <- rbindlist(list("aBC" = validSimListbaseCase$landscapeMAD,
                                     "bAP" = validSimListaltParameters$landscapeMAD),
                                idcol = "simulation")
landscapeMAD_BC_AP <- melt(landscapeMAD_BC_AP,
                           measure.vars = c("meanAbsDevRelAbund", "meanAbsDevCount", "meanAbsDevCountDom", "meanAbsDevDeltaB"),
                           value.name = "MAD")

pixelMAD_BC_AP <- rbindlist(list("aBC" = validSimListbaseCase$pixelMAD,
                                 "bAP" = validSimListaltParameters$pixelMAD),
                            idcol = "simulation")

pixelMAD_BC_AP <- melt(pixelMAD_BC_AP,
                       measure.vars = c("meanAbsDevRelAbund", "meanAbsDevDeltaB"),
                       value.name = "MAD")

## bind to keep all factor levels for colours
landPixelMAD_BC_AP <- rbindlist(list("landscape" = landscapeMAD_BC_AP,
                                     "pixel" = pixelMAD_BC_AP), fill = TRUE,
                                idcol = "level")

landPixelMAD_BC_AP[, variableFacet := factor(variable,
                                             levels = unique(variable),
                                             labels = c("frac('species B', 'total/pixel B')",
                                                        "paste('No. of pixels')",
                                                        "paste('No. of pixels ')",
                                                        "g/m^2"))]

speciesLabels <- validSimListbaseCase@.envir$.mods$Biomass_validationKNN$.objects$speciesLabels
MADLabels <- list("meanAbsDevRelAbund" = "species rel. abundance",
                  "meanAbsDevCount" = "species presences",
                  "meanAbsDevCountDom" = "species dominance",
                  "meanAbsDevDeltaB" = bquote(paste("species", Delta, B)))

cols <- c("speciesCode", "relAbund", "dataType", "year", "simulation")
plotLandRelB_BC_AP <- ggplot(na.omit(landscapeVars_BC_AP[dataType == "simulated", ..cols]),
                             aes(x = speciesCode, y = relAbund, fill = simulation)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = na.omit(landscapeVars_BC_AP[dataType == "observed", ..cols]),
               aes(x = speciesCode, y = relAbund, colour = "observed"),
               fun = "mean", geom = "point", size = 2) +
  scale_fill_brewer(palette = "Paired", labels = c("aBC" = "base case", "bAP" = "alternative param.")) +
  scale_x_discrete(labels = speciesLabels) +
  scale_color_manual(values = c("observed" = "red3")) +
  scale_y_continuous(limits = c(0,1)) +
  plotTheme(base_size = 12, margin = FALSE, x.text.angle = 45, legend = "bottom") +
  facet_grid(~ year) +
  labs(title = "Species relative abundances",
       x = "", y = expression(over("species B", "total B")),
       colour = "", fill = "") +
  guides(fill = guide_legend(override.aes = list(colour = "transparent")))

cols <- c("speciesCode", "count", "dataType", "year", "simulation")
plotLandPres_BC_AP <- ggplot(na.omit(landscapeVars_BC_AP[dataType == "simulated", ..cols]),
                             aes(x = speciesCode, y = count, fill = simulation)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = na.omit(landscapeVars_BC_AP[dataType == "observed", ..cols]),
               aes(x = speciesCode, y = count, colour = "observed"),
               fun = "mean", geom = "point", size = 2) +
  scale_fill_brewer(palette = "Paired", labels = c("aBC" = "base case", "bAP" = "alternative param.")) +
  scale_x_discrete(labels = speciesLabels) +
  scale_color_manual(values = c("observed" = "red3")) +
  plotTheme(base_size = 12, margin = FALSE, x.text.angle = 45, legend = "bottom") +
  facet_grid(~ year) +
  labs(title = "Species' presences",
       x = "", y = "No. of pixels",
       colour = "", fill = "") +
  guides(fill = guide_legend(override.aes = list(colour = "transparent")))

cols <- c("speciesCode", "countDom", "dataType", "year", "simulation")
plotLandDom_BC_AP <- ggplot(na.omit(landscapeVars_BC_AP[dataType == "simulated", ..cols]),
                            aes(x = speciesCode, y = countDom, fill = simulation)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = na.omit(landscapeVars_BC_AP[dataType == "observed", ..cols]),
               aes(x = speciesCode, y = countDom, colour = "observed"),
               fun = "mean", geom = "point", size = 2) +
  scale_fill_brewer(palette = "Paired", labels = c("aBC" = "base case", "bAP" = "alternative param.")) +
  scale_x_discrete(labels = speciesLabels) +
  scale_color_manual(values = c("observed" = "red3")) +
  plotTheme(base_size = 12, margin = FALSE, x.text.angle = 45, legend = "bottom") +
  facet_grid(~ year) +
  labs(title = "Species' dominance",
       x = "", y = "No. of pixels",
       colour = "", fill = "") +
  guides(fill = guide_legend(override.aes = list(colour = "transparent")))

cols <- c("speciesCode", "deltaB", "dataType", "simulation", "rep")
plotLandDeltaB_BC_AP <- ggplot(na.omit(landscapeVars_BC_AP[dataType == "simulated", ..cols]),
                               aes(x = speciesCode, y = deltaB, fill = simulation)) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge") +
  stat_summary(fun.data = "mean_sd", geom = "linerange", size = 1) +
  stat_summary(data = na.omit(landscapeVars_BC_AP[dataType == "observed", ..cols]),
               aes(x = speciesCode, y = deltaB, colour = "observed"),
               fun = "mean", geom = "point", size = 2) +
  scale_fill_brewer(palette = "Paired", labels = c("aBC" = "base case", "bAP" = "alternative param.")) +
  scale_x_discrete(labels = speciesLabels, drop = TRUE) +
  scale_color_manual(values = c("observed" = "red3")) +
  plotTheme(base_size = 12, margin = FALSE, x.text.angle = 45, legend = "bottom") +
  labs(title = expression(paste("Changes in biomass (", Delta, "B)")),
       x = "", y = expression(g/m^2),
       colour = "", fill = "") +
  guides(fill = guide_legend(override.aes = list(colour = "transparent")))

## exclude deltaB (and landscape level from x axis given that relative B, presence and dominance is always at the species level)
# plotLandMAD_BC_AP <- ggplot(landPixelMAD_BC_AP[level == "landscape"],
plotLandMAD_BC_AP <- ggplot(landPixelMAD_BC_AP[level == "landscape" & variable != "meanAbsDevDeltaB" & speciesCode != "landscape"],
                            aes(x = speciesCode, y = MAD, colour = variable, alpha = simulation)) +
  stat_summary(fun.data = "mean_sdl", position = position_dodge(width = 0.5),
               size = 0.4) +
  scale_x_discrete(labels = speciesLabels) +
  scale_y_continuous(breaks = scales::extended_breaks(n = 5)) +
  scale_color_brewer(palette = "Dark2", labels = MADLabels, drop = TRUE) +
  scale_alpha_manual(values = c("aBC" = 0.6, "bAP" = 1.0),
                     labels = c("aBC" = "base case", "bAP" = "alternative param.")) +
  plotTheme(base_size = 12, legend = "top", x.text.angle = 45) +
  labs(title = "Mean abs. deviation",
       colour = "", x = "", y = "", alpha = "") +
  theme(strip.placement = "outside", strip.switch.pad.wrap = unit(0, "cm"),
        strip.background = element_blank(), legend.position = "bottom") +
  facet_wrap(~ variableFacet, ncol = 1, scales = "free_y",
             labeller = label_parsed,
             strip.position = "left")

cols <- c("speciesCode", "relAbund", "year", "dataType2", "simulation", "rep")
plotPixelRelB_BC_AP <- ggplot() +
  geom_boxplot(data = na.omit(pixelVars_BC_AP[, ..cols]),
               aes(x = speciesCode, y = relAbund, fill = dataType2,
                   alpha = speciesCode == "pixel")) +
  scale_x_discrete(labels = c(speciesLabels, "pixel" = "pixel")) +
  scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 1.0), guide = "none") +
  scale_fill_brewer(palette = "Paired",
                    labels = c("aBC" = "base case", "bAP" = "alternative param.", "observed" = "observed")) +
  plotTheme(base_size = 12, margin = FALSE, x.text.angle = 45, legend = "bottom") +
  facet_grid(~ year) +
  labs(title = "Species relative abundances",
       x = "", y = expression(over("species B", "pixel B")),
       colour = "", fill = "")

cols <- c("speciesCode", "deltaB", "dataType2", "simulation", "rep")
plotPixelDeltaB_BC_AP <- ggplot() +
  geom_boxplot(data = na.omit(pixelVars_BC_AP[, ..cols]),
               aes(x = speciesCode, y = deltaB, fill = dataType2,
                   alpha = speciesCode == "pixel")) +
  scale_x_discrete(labels = c(speciesLabels, "pixel" = "pixel")) +
  scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 1.0), guide = "none") +
  scale_fill_brewer(palette = "Paired",
                    labels = c("aBC" = "base case", "bAP" = "alternative param.", "observed" = "observed")) +
  plotTheme(base_size = 12, margin = FALSE, x.text.angle = 45, legend = "bottom") +
  labs(title = expression(paste("Changes in biomass (", Delta, "B)")),
       x = "", y = expression(g/m^2),
       colour = "", fill = "")

## exclude deltaB (and pixel level from x axis given that relative B is always at the species level)
# plotPixelMAD_BC_AP <- ggplot(landPixelMAD_BC_AP[level == "pixel"],
plotPixelMAD_BC_AP <- ggplot(landPixelMAD_BC_AP[level == "pixel" & variable != "meanAbsDevDeltaB" & speciesCode != "pixel"],
                             aes(x = speciesCode, y = MAD, colour = variable, alpha = simulation)) +
  stat_summary(fun.data = "mean_sdl", position = position_dodge(width = 0.5),
               size = 0.4) +
  scale_x_discrete(labels = speciesLabels) +
  scale_y_continuous(breaks = scales::extended_breaks(n = 5)) +
  scale_color_brewer(palette = "Dark2", labels = MADLabels, drop = FALSE) +
  scale_alpha_manual(values = c("aBC" = 0.6, "bAP" = 1.0),
                     labels = c("aBC" = "base case", "bAP" = "alternative param.")) +
  plotTheme(base_size = 12, legend = "top", x.text.angle = 45) +
  labs(title = "Mean abs. deviation",
       colour = "", y = "", x = "", alpha = "") +
  theme(strip.placement = "outside", strip.switch.pad.wrap = unit(0, "cm"),
        strip.background = element_blank(), legend.position = "bottom") +
  facet_wrap(~ variableFacet, ncol = 1, scales = "free_y",
             labeller = label_parsed,
             strip.position = "left")

## make a plot with all plots
allPlotsValid_BC_AP <- plot_grid(
  plotLandRelB_BC_AP + theme(legend.position = "none", text = element_text(size = 10)),
  plotLandDom_BC_AP + theme(legend.position = "none", text = element_text(size = 10)),
  NULL, NULL,
  plotLandPres_BC_AP + theme(legend.position = "none", text = element_text(size = 10)),
  plotLandMAD_BC_AP + theme(legend.position = "none", text = element_text(size = 10)),
  NULL, NULL,
  get_legend(plotLandRelB_BC_AP), NULL,
  NULL, NULL,
  plotPixelRelB_BC_AP + theme(legend.position = "none", text = element_text(size = 10)),
  plotPixelMAD_BC_AP + theme(legend.position = "none", text = element_text(size = 10)),
  NULL, NULL,
  get_legend(plotPixelRelB_BC_AP), get_legend(plotLandMAD_BC_AP + theme(legend.box = "vertical", legend.margin = margin(0,0,0,0))),
  rel_heights = c(0.9, -0.05, 0.9, -0.05, 0.1, 0.05, 0.9, -0.01, 0.2),
  ncol = 2, byrow = TRUE,
  labels = c("a)", "", "","", "", "", "", "", "", "", "", "", "b)"),
  align = "hv", axis = "bl"
)

ggsave(plot = allPlotsValid_BC_AP, filename = file.path(figDir, "validationPlots_example2.png"),
       bg = "white", height = 250, width = 230, dpi = 300, units = "mm")

## SIMULATION & VALIDATION FIGURES
## only total biomass and MAD
yMax <- max(speciesB_BC$data[, sum(BiomassBySpecies), by = year]$V1, speciesB_AP$data[, sum(BiomassBySpecies), by = year]$V1)

alignedPlots <- align_plots(speciesB_BC + scale_y_continuous(limits = c(0, yMax)) +
                              labs(title = "Base case simulation", y = expression(atop("Total biomass", g/m^2))) +
                              theme(plot.margin = margin(b = 0), legend.position = "none", text = element_text(size = 10)),
                            speciesB_AP + scale_y_continuous(limits = c(0, yMax)) +
                              labs(title = "Alternative param. simulation", y = "") +
                              theme(plot.margin = margin(b = 0), legend.position = "none", text = element_text(size = 10)),
                            plotLandMAD_BC_AP + labs(title = "Landscape-level MAD") +
                              theme(plot.margin = margin(t = 0), legend.position = "none", text = element_text(size = 10)),
                            plotPixelMAD_BC_AP + labs(title = "Pixel-level MAD") +
                              theme(plot.margin = margin(t = 0), legend.position = "none", text = element_text(size = 10)),
                            align = "hv", axis = "bl")
B_MAD_PlotsSim_BC_AP <- ggarrange(
  ggarrange(alignedPlots[[1]],
            alignedPlots[[2]],
            ncol = 2, labels = c("a)", "b)")),
  get_legend(speciesB_BC +
               theme(legend.direction = "vertical", text = element_text(size = 10))),
  ggarrange(alignedPlots[[3]],
            alignedPlots[[4]],
            ncol = 2, labels = c("c)", "d)")),
  get_legend(plotLandMAD_BC_AP +
               theme(legend.direction = "vertical", legend.box.just = "left",
                     legend.spacing = unit(0, units = "mm"), legend.box = "vertical", text = element_text(size = 10))),
  nrow = 2, ncol = 2, heights = c(1, 1), widths = c(1, 0.2)
)
ggsave(plot = B_MAD_PlotsSim_BC_AP, filename = file.path(figDir, "simValidPlots_B_MAD_example2.png"),
       bg = "white", height = 160, width = 230, dpi = 300, units = "mm")

