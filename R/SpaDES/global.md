---
title: "LandR Publication Simulations"
author: "Ceres Barros"
date: "Last updated: 2022-10-24"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    toc_depth: 4
    theme: sandstone
    number_sections: false
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
  markdown: 
    wrap: 80
always_allow_html: true
---

This provides a simplified example of publication runs. For full code used in 
publication see `global.R`



# Package installation and setup

```pkginstall
## set CRAN repo
options(repos = c(CRAN = "https://cloud.r-project.org"))

## package installation location
pkgPath <- file.path("packages", version$platform,
                     paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1]))
dir.create(pkgPath, recursive = TRUE)
.libPaths(pkgPath, include.site = FALSE)

if (!"remotes" %in% installed.packages()) {
  install.packages("remotes")
}

if (!"Require" %in% installed.packages(lib.loc = pkgPath) ||
    packageVersion("Require", lib.loc = pkgPath) < "0.1.6.9005") {
  remotes::install_github("PredictiveEcology/Require@8ac756d66bc145f226f39e703e3787d3aed159f2",
                          upgrade = FALSE, force = TRUE)
}

## use binary linux packages if on Ubuntu
Require::setLinuxBinaryRepo()

Require::Require("PredictiveEcology/SpaDES.project@6d7de6ee12fc967c7c60de44f1aa3b04e6eeb5db",
                 require = FALSE, upgrade = FALSE, standAlone = TRUE)

SpaDES.project::getModule(modulePath = modulePath,
                          c("PredictiveEcology/Biomass_speciesData@master",
                            "PredictiveEcology/Biomass_borealDataPrep@master",
                            "PredictiveEcology/Biomass_core@master",
                            "CeresBarros/Biomass_validationKNN@master",
                            "CeresBarros/Biomass_speciesParameters@LandRPub"))

pkgSnapshotFile <- file.path("packages",
                             paste0("pkgSnapshot_",
                                    paste0(version$major, "_", strsplit(version$minor, "[.]")[[1]][1]),
                                    ".txt"))

Require::Require(packageVersionFile = pkgSnapshotFile,
                 require = FALSE, standAlone = TRUE)
```

## Package loading and options, folder directories setup

```r
Require::Require(c("raster", "terra", "dplyr", "data.table", "future",
                   "SpaDES.core", "LandR", "reproducible",
                   "ggspatial", "ggpubr", "cowplot"),
                 upgrade = FALSE, install = FALSE)
```

```
## Loading required package: terra
```

```
## Warning: package 'terra' was built under R version 4.2.1
```

```
## terra 1.6.17
```

```
## 
## Attaching package: 'terra'
```

```
## The following object is masked from 'package:LandR':
## 
##     rescale
```

```
## The following object is masked from 'package:SpaDES.tools':
## 
##     wrap
```

```
## The following object is masked from 'package:SpaDES.core':
## 
##     time<-
```

```
## The following object is masked from 'package:data.table':
## 
##     shift
```

```
## Loading required package: future
```

```
## Warning: package 'future' was built under R version 4.2.1
```

```
## 
## Attaching package: 'future'
```

```
## The following object is masked from 'package:terra':
## 
##     values
```

```
## The following object is masked from 'package:raster':
## 
##     values
```

```
## Loading required package: ggpubr
```

```
## Warning: package 'ggpubr' was built under R version 4.2.1
```

```
## 
## Attaching package: 'ggpubr'
```

```
## The following object is masked from 'package:terra':
## 
##     rotate
```

```
## The following object is masked from 'package:raster':
## 
##     rotate
```

```
## Loading required package: cowplot
```

```
## Warning: package 'cowplot' was built under R version 4.2.1
```

```
## 
## Attaching package: 'cowplot'
```

```
## The following object is masked from 'package:ggpubr':
## 
##     get_legend
```

```r
options("reproducible.useNewDigestAlgorithm" = 2)
options("spades.moduleCodeChecks" = FALSE)
options("spades.inputPath" = Require::normPath(file.path("R/SpaDES/inputs")))  
options("reproducible.useCache" = TRUE)
options("reproducible.destinationPath" = Require::normPath(file.path("R/SpaDES/inputs")))
options("reproducible.useGDAL" = FALSE)
options("spades.useRequire" = FALSE)
options("Require.unloadNamespaces" = FALSE)
options("reproducible.useTerra" = FALSE)

## set run name and paths
runName <- "baseCase"
eventCaching <- c(".inputObjects", "init")
useParallel <- FALSE

## paths
simPathName <- "oct2022Runs"
simPaths <- list(cachePath = Require::normPath(file.path("R/SpaDES/cache", simPathName))
                 , modulePath = Require::normPath(modulePath)
                 , inputPath = Require::normPath(file.path("R/SpaDES/inputs"))
                 , outputPath = Require::normPath(file.path("R/SpaDES/outputs", simPathName, runName)))
```

# Simulation setup

## Part 1 - Get study area and other necessary objects 


```r
## Get necessary objects -----------------------
devtools::source_url(paste0("https://raw.githubusercontent.com/CeresBarros/",
                            "LandRBiomass_publication/repPkgInstall/R/SpaDES/",
                            "1_simObjects.R?raw=TRUE"))
```

## Part 2 - Species layers


```r
## Run Biomass_speciesData to get species layers
## running this separately from other modules makes switching
## between using a large and a smaller study area easier when the smaller SA is within the large one,
## as it keeps the data in separate folders thatn can be used across simulations/scenarios
devtools::source_url(paste0("https://raw.githubusercontent.com/CeresBarros/",
                            "LandRBiomass_publication/repPkgInstall/R/SpaDES/",
                            "2_speciesLayers.R?raw=TRUE"))

## check species layers:
plot(simOutSpeciesLayers$speciesLayers)

## subset sppEquivalencies and colorVector
sppEquivalencies_CA <- sppEquivalencies_CA[Boreal %in% names(simOutSpeciesLayers$speciesLayers)]
sppColorVect <- sppColorVect[c(names(simOutSpeciesLayers$speciesLayers), "Mixed")]

## Get land-cover raster now that we have a rasterToMatchLarge
rstLCC2005 <- Cache(prepInputs,
                    targetFile = "NA_LandCover_2005_V3_25haMMU.tif",
                    url = paste0("http://www.cec.org/wp-content/uploads/",
                                 "wpallimport/files/Atlas/Files/Land_Cover_2005/",
                                 "Land_Cover_2005v3_TIFF.zip"),
                    destinationPath = simPaths$inputPath,
                    studyArea = simOutSpeciesLayers$studyAreaLarge,
                    rasterToMatch = simOutSpeciesLayers$rasterToMatchLarge,
                    filename2 = .suffix("rstLCC.tif", paste0("_", SAname)),
                    overwrite = TRUE,
                    cacheRepo = simPaths$cachePath,
                    userTags = c("rstLCC", SAname),
                    omitArgs = c("userTags"))
```



## Part 3 - Module parameters and outputs

* `vegLeadingProportion` indicates what proportion the stand must be in one species group for it to be leading.
If 0, then there is never a notion of "mixed" vegetation types and a species is a leading species if it has the highest relative biomass in the pixel.
* `successionTimestep` defines the frequency at which dispersal and age reclassification occurs - every 10 years is the default LANDIS behaviour. 


```r
## simulation params
simTimes <- list(start = 2001, end = 2031)
vegLeadingProportion <- 0 
successionTimestep <- 10L  

speciesParams <- list(
  "shadetolerance" = list(
    Betu_Pap = 1,
    Lari_Lar = 1,
    Pice_Gla = 2,
    Pice_Mar = 3,
    Pinu_Ban = 1.5,
    Popu_Spp = 1
  )
)

simModules <- list("Biomass_borealDataPrep"
                   , "Biomass_speciesParameters"
                   , "Biomass_core"
)

simParams <- list(
  .globals = list("dataYear" = 2001L
                  , "sppEquivCol" = sppEquivCol
                  , "vegLeadingProportion" = vegLeadingProportion
                  , ".sslVerify" = 0L
                  , ".useCache" = eventCaching),
  Biomass_borealDataPrep = list(
    "fitDeciduousCoverDiscount" = TRUE
    , "subsetDataAgeModel" = FALSE
    , "subsetDataBiomassModel" = FALSE
    , "exportModels" = "all"
    , "fixModelBiomass" = TRUE
    , "fireURL" = paste0("https://drive.google.com/file/d/",
                         "1YIc_BSkPKqW60SmfpR2vDpeRGrwOFKso/view?usp=sharing")  
    , "speciesTableAreas" = c("BSW", "BP")
    , "speciesUpdateFunction" = list(
      quote(LandR::speciesTableUpdate(sim$species, 
                                      sim$speciesTable, 
                                      sim$sppEquiv, 
                                      P(sim)$sppEquivCol)),
      quote(LandR::updateSpeciesTable(speciesTable = sim$species, 
                                      params = sim$speciesParams)))
    , "pixelGroupAgeClass" = successionTimestep
    , "pixelGroupBiomassClass" = 100
    , "useCloudCacheForStats" = FALSE
    , "cloudFolderID" = NA
    , ".plots" = c("object", "raw")
  )
  , Biomass_speciesParameters = list(
    "quantileAgeSubset" = list(Betu_Pap = 95, 
                               Lari_Lar = 95, 
                               Pice_Gla = 95, 
                               Pice_Mar = 95, 
                               Pinu_Ban = 99, 
                               Popu_Spp = 99)
  )
  , Biomass_core = list(
    "calcSummaryBGM" = c("start")
    , "initialBiomassSource" = "cohortData"
    , "plotOverstory" = TRUE
    , "seedingAlgorithm" = "wardDispersal"
    , "successionTimestep" = successionTimestep
    , ".plotInitialTime" = simTimes$start
    , ".plotInterval" = 1L
    , ".plots" = c("object", "raw")
    , ".plotMaps" = FALSE
    , ".saveInitialTime" = NA
    , ".useCache" = eventCaching[1]
    , ".useParallel" = useParallel
  )
)

simObjects <- list(
  "rstLCC" = rstLCC2005
  , "sppEquiv" = sppEquivalencies_CA
  , "sppColorVect" = sppColorVect
  , "speciesLayers" = simOutSpeciesLayers$speciesLayers
  , "speciesParams" = speciesParams
  , "rasterToMatchLarge" = simOutSpeciesLayers$rasterToMatchLarge
  , "studyArea" = studyAreaS
  , "studyAreaLarge" = studyAreaL
  , "treed" = simOutSpeciesLayers$treed
  , "numTreed" = simOutSpeciesLayers$numTreed
  , "nonZeroCover" = simOutSpeciesLayers$nonZeroCover
)

## Make a rasterToMatch now that we have a rasterToMatchLarge and a studyArea -- use terra here.
RTM <- try(project(rast(simObjects$rasterToMatchLarge), y = crs(vect(simObjects$studyArea))))
RTM <- crop(RTM, vect(simObjects$studyArea), mask = TRUE)
RTM <- raster(RTM)
RTM[!is.na(RTM[])] <- 1L

simObjects$rasterToMatch <- RTM

## objects will be saved at the start of the simulation (so they reflect the previous year)
simOutputs <- data.frame(
  expand.grid(objectName = "cohortData",
              saveTime = unique(seq(simTimes$start, simTimes$end, by = 1)),
              eventPriority = 1,
              stringsAsFactors = FALSE)
  )
simOutputs <- rbind(simOutputs, 
                    data.frame(objectName = "pixelGroupMap",
                               saveTime = unique(seq(simTimes$start, simTimes$end, by = 1)),
                               eventPriority = 1))
simOutputs <- rbind(simOutputs, 
                    data.frame(objectName = "biomassMap",
                               saveTime = simTimes$start,
                               eventPriority = 1))

## in the first year, eventPriorities need to be set to AFTER the init event (which has priority 1)
simOutputs$eventPriority[simOutputs$saveTime == simTimes$start] <- 1.5

## make a initialisation simList and run init events too
LandRBiomass_simInit <- Cache(simInitAndSpades
                              , times = simTimes
                              , params = simParams
                              , modules = simModules
                              , loadOrder = unlist(simModules)
                              , objects = simObjects
                              , paths = simPaths
                              , outputs = simOutputs
                              , events = "init"
                              , .plotInitialTime = NA
                              , .studyAreaName = SAname
                              , userTags = "simInitAndInits"
                              , omitArgs = c("userTags", ".plotInitialTime"))
## save the simList
saveSimList(LandRBiomass_simInit, file.path(simPaths$outputPath, paste0("simInit", runName, ".qs")))
```




```r
# look at the model structure and object exchange
objectDiagram(LandRBiomass_simInit, width = 1000, height = 2500)
```

```{=html}
<div id="htmlwidget-fd07f7ea9bb154448a16" style="width:1000px;height:2500px;" class="DiagrammeR html-widget"></div>
<script type="application/json" data-for="htmlwidget-fd07f7ea9bb154448a16">{"x":{"diagram":"sequenceDiagram\nBiomass_borealDataPrep ->> Biomass_borealDataPrep : rawBiomassMap\nBiomass_borealDataPrep ->> Biomass_borealDataPrep : studyArea\nBiomass_borealDataPrep ->> Biomass_core : biomassMap\nBiomass_borealDataPrep ->> Biomass_core : cohortData\nBiomass_borealDataPrep ->> Biomass_core : ecoregion\nBiomass_borealDataPrep ->> Biomass_core : ecoregionMap\nBiomass_borealDataPrep ->> Biomass_core : minRelativeB\nBiomass_borealDataPrep ->> Biomass_core : pixelGroupMap\nBiomass_borealDataPrep ->> Biomass_core : species\nBiomass_borealDataPrep ->> Biomass_core : speciesEcoregion\nBiomass_borealDataPrep ->> Biomass_core : studyArea\nBiomass_borealDataPrep ->> Biomass_core : sufficientLight\nBiomass_borealDataPrep ->> Biomass_speciesParameters : species\nBiomass_borealDataPrep ->> Biomass_speciesParameters : speciesEcoregion\nBiomass_core ->> Biomass_borealDataPrep : speciesLayers\nBiomass_core ->> Biomass_core : cohortData\nBiomass_core ->> Biomass_core : ecoregionMap\nBiomass_core ->> Biomass_core : lastReg\nBiomass_core ->> Biomass_core : minRelativeB\nBiomass_core ->> Biomass_core : pixelGroupMap\nBiomass_core ->> Biomass_core : species\nBiomass_core ->> Biomass_core : speciesEcoregion\nBiomass_core ->> Biomass_core : speciesLayers\nBiomass_core ->> Biomass_core : treedFirePixelTableSinceLastDisp\nBiomass_core ->> Biomass_speciesParameters : species\nBiomass_core ->> Biomass_speciesParameters : speciesEcoregion\nBiomass_speciesParameters ->> Biomass_core : species\nBiomass_speciesParameters ->> Biomass_core : speciesEcoregion\nBiomass_speciesParameters ->> Biomass_speciesParameters : species\nBiomass_speciesParameters ->> Biomass_speciesParameters : speciesEcoregion\n_INPUT_ ->> Biomass_borealDataPrep : cloudFolderID\n_INPUT_ ->> Biomass_borealDataPrep : columnsForPixelGroups\n_INPUT_ ->> Biomass_borealDataPrep : ecoregionLayer\n_INPUT_ ->> Biomass_borealDataPrep : ecoregionRst\n_INPUT_ ->> Biomass_borealDataPrep : rasterToMatch\n_INPUT_ ->> Biomass_borealDataPrep : rasterToMatchLarge\n_INPUT_ ->> Biomass_borealDataPrep : rstLCC\n_INPUT_ ->> Biomass_borealDataPrep : speciesTable\n_INPUT_ ->> Biomass_borealDataPrep : sppColorVect\n_INPUT_ ->> Biomass_borealDataPrep : sppEquiv\n_INPUT_ ->> Biomass_borealDataPrep : sppNameVector\n_INPUT_ ->> Biomass_borealDataPrep : standAgeMap\n_INPUT_ ->> Biomass_borealDataPrep : studyAreaLarge\n_INPUT_ ->> Biomass_core : cceArgs\n_INPUT_ ->> Biomass_core : rasterToMatch\n_INPUT_ ->> Biomass_core : sppColorVect\n_INPUT_ ->> Biomass_core : sppEquiv\n_INPUT_ ->> Biomass_core : studyAreaReporting\n_INPUT_ ->> Biomass_speciesParameters : PSPgis\n_INPUT_ ->> Biomass_speciesParameters : PSPmeasure\n_INPUT_ ->> Biomass_speciesParameters : PSPplot\n_INPUT_ ->> Biomass_speciesParameters : factorialSpeciesTable\n_INPUT_ ->> Biomass_speciesParameters : reducedFactorialCohortData\n_INPUT_ ->> Biomass_speciesParameters : sppEquiv\n_INPUT_ ->> Biomass_speciesParameters : studyAreaANPP\n"},"evals":[],"jsHooks":[]}</script>
```

# Run simulation

Here we run just one repetition


```r
LandRBiomass_sim <- spades(LandRBiomass_simInit, .plotInitialTime = simTimes$start)
## save
saveSimList(LandRBiomass_sim, file.path(simPaths$outputPath, paste0("simList_LandRBiomass_sim_1rep_", runName, ".qs")))
```



If we were to run several repetitions, this would be how:


```r
plan("sequential")    ## don't parallelize -- future is not passing `options` to the nodes
LandRBiomass_sim <- experiment2(
  sim1 = LandRBiomass_sim,
  clearSimEnv = TRUE,
  replicates = 10)

## save simLists object.
qs::qsave(LandRBiomass_sim, file.path(simPaths$outputPath, paste0("simList_LandRBiomass_sim_", runName, ".qs")))
```

## Inspect simulation objects

We can use the `simList` objects to plot simulation objects, such as the input layers used for parameterisation.


```r
## study area within Saskatchewan province
## get Canadian provinces and subset to SK
can1 <- raster::getData('GADM', country = "CAN", level = 1, path = tempdir())
```

```
## Warning in raster::getData("GADM", country = "CAN", level = 1, path = tempdir()): getData will be removed in a future version of raster
## . Please use the geodata package instead
```

```r
can1 <- can1[can1$NAME_1 == "Saskatchewan",]
can1 <- spTransform(can1, CRSobj = crs(LandRBiomass_simInit$studyArea))

plotStudyAreas <- ggplot() +
  layer_spatial(can1, fill = "grey90") + 
  layer_spatial(LandRBiomass_simInit$studyAreaLarge, fill = "darkgreen") +
  layer_spatial(LandRBiomass_simInit$studyArea, fill = "green") +
  labs(fill = "") +
  theme_void()


plotEcodist <- ggplot() +
  layer_spatial(LandRBiomass_simInit$ecoregionLayer, aes(fill = as.factor(ECODISTRIC))) + 
  labs(fill = "") +
  theme_void()
```


```r
## study area within Saskatchewan province
clearPlot(force = TRUE)
Plot(plotStudyAreas, title = " ")

## input stand biomass and age, ecological zonation (ecodistricts) and land-cover (LCC 2005)
clearPlot(force = TRUE)
Plot(LandRBiomass_simInit$rawBiomassMap,
     LandRBiomass_simInit$standAgeMap, 
     LandRBiomass_simInit$rstLCC,
     plotEcodist,
     title = c("kNN stand biomass", "kNN stand age", "land-cover", "ecodistricts"))

## species % cover
clearPlot(force = TRUE)
Plot(simOutSpeciesLayers$speciesLayers,
     title = paste(names(simOutSpeciesLayers$speciesLayers), collapse = ", "))
```

![](global_files/figure-html/plots2-1.png)<!-- -->![](global_files/figure-html/plots2-2.png)<!-- -->![](global_files/figure-html/plots2-3.png)<!-- -->![](global_files/figure-html/plots2-4.png)<!-- -->


Similarly, we can have a look at the species traits values used in the simulation directly from the `simList` object (although we also chose to save them).


```r
LandRBiomass_simInit$species
LandRBiomass_simInit$minRelativeB ## can be spatially varying, but identical across ecolocation (AKA ecoregion)
LandRBiomass_simInit$sufficientLight
LandRBiomass_simInit$speciesEcoregion
```

For example, these were the (spatially) invariant species traits used in the simulation:


Table: Invariant species traits

|species  |Area | longevity| sexualmature| shadetolerance| firetolerance| seeddistance_eff| seeddistance_max| resproutprob| resproutage_min| resproutage_max|postfireregen | leaflongevity| wooddecayrate| mortalityshape| growthcurve| leafLignin|hardsoft |speciesCode | mANPPproportion| inflationFactor|
|:--------|:----|---------:|------------:|--------------:|-------------:|----------------:|----------------:|------------:|---------------:|---------------:|:-------------|-------------:|-------------:|--------------:|-----------:|----------:|:--------|:-----------|---------------:|---------------:|
|Betu_Pap |BP   |       140|           20|            1.0|             1|              200|             5000|          0.5|              10|              70|resprout      |             1|          0.07|             15|         0.5|        0.1|hard     |Betu_Pap    |            5.00|        1.115822|
|Lari_Lar |BP   |       350|           40|            1.0|             1|               50|              200|          0.0|               0|               0|none          |             1|          0.02|             25|         0.0|        0.2|soft     |Lari_Lar    |            3.33|        1.000000|
|Pice_Gla |BP   |       400|           30|            2.0|             2|              100|              303|          0.0|               0|               0|none          |             3|          0.02|             15|         0.5|        0.2|soft     |Pice_Gla    |            3.25|        1.042535|
|Pice_Mar |BP   |       250|           30|            3.0|             2|               80|              200|          0.0|               0|               0|serotiny      |             3|          0.02|             15|         0.5|        0.2|soft     |Pice_Mar    |            4.25|        1.060445|
|Pinu_Ban |BP   |       150|           20|            1.5|             2|               30|              100|          0.0|               0|               0|serotiny      |             2|          0.01|             15|         0.5|        0.2|soft     |Pinu_Ban    |            2.75|        1.355014|
|Popu_Spp |BP   |       200|           20|            1.0|             2|             1000|             5000|          0.9|              10|             130|resprout      |             1|          0.07|             25|         0.5|        0.1|hard     |Popu_Spp    |            2.00|        1.261034|

# Validation

Here we run the validation on the outputs of just one repetition.

## Validation setup

We begin by preparing all the inputs necessary for the `Biomass_validationKNN` module.


```r
## get the  land-cover change map (needed to have an RTM first, so get it from the simInitList)
## /!\ it is assumed that the filename of the raster in the simList corresponds to the raster found in disk.
## this may not be the case if the simulations were run in another machine and saved rasters were not imported.

validationPaths <- list(cachePath = file.path("R/SpaDES/cache", simPathName)
                        , modulePath = file.path("R/SpaDES/m")
                        , inputPath = file.path("R/SpaDES/inputs")
                        , outputPath = file.path("R/SpaDES/validation", simPathName))

devtools::source_url(paste0("https://raw.githubusercontent.com/CeresBarros/",
                            "LandRBiomass_publication/repPkgInstall/R/SpaDES/",
                            "3_simObjects4Valid.R?raw=TRUE"))

## PARAMETERS FOR VALIDATION MODULE
## in this case all reps have the same parameters, so we can use the first rep to get the values
validationTimes <- list(start = 1, end = 1)
validationParams <- list(
  Biomass_validationKNN = list(
    "minCoverThreshold" = params(LandRBiomass_simInit)$Biomass_borealDataPrep$minCoverThreshold
    , "pixelGroupBiomassClass" = params(LandRBiomass_simInit)$Biomass_borealDataPrep$pixelGroupBiomassClass
    , "deciduousCoverDiscount" = params(LandRBiomass_simInit)$Biomass_borealDataPrep$deciduousCoverDiscount
    , "sppEquivCol" = params(LandRBiomass_simInit)$Biomass_borealDataPrep$sppEquivCol
    , "validationReps" = as.integer(1:10)  ## or length of simLists
    , "validationYears" = as.integer(c(2001, 2011))
    , ".plots" = c("png")
    , ".useCache" = eventCaching
  )
)

validationObjects <- list(
  "biomassMap" = biomassMap
  , "rasterToMatch" = rasterToMatch
  , "rawBiomassMapStart" = rawBiomassMap
  , "simulationOutputs" = simulationOutputs
  , "speciesLayersStart" = speciesLayers
  , "sppColorVect" = LandRBiomass_simInit$sppColorVect
  , "sppEquiv" = LandRBiomass_simInit$sppEquiv
  , "standAgeMapStart" = standAgeMap
  , "studyArea" = LandRBiomass_simInit$studyArea
)

## the following objects are only saved once at the end of year 0/beggining of year 1 (they don't change)
validationOutputs <- data.frame(expand.grid(objectName = c("rawBiomassMapStart"),
                                            saveTime = c(validationTimes$start),
                                            eventPriority = 1),
                                stringsAsFactors = FALSE)
validationOutputs <- rbind(validationOutputs, data.frame(objectName = "rawBiomassMapEnd",
                                                         saveTime = c(validationTimes$start),
                                                         eventPriority = 1))
validationOutputs <- rbind(validationOutputs, data.frame(objectName = "standAgeMapStart",
                                                         saveTime = c(validationTimes$start),
                                                         eventPriority = 1))
validationOutputs <- rbind(validationOutputs, data.frame(objectName = "standAgeMapEnd",
                                                         saveTime = c(validationTimes$start),
                                                         eventPriority = 1))
validationOutputs <- rbind(validationOutputs, data.frame(objectName = "speciesLayersStart",
                                                         saveTime = c(validationTimes$start),
                                                         eventPriority = 1))
validationOutputs <- rbind(validationOutputs, data.frame(objectName = "speciesLayersEnd",
                                                         saveTime = c(validationTimes$start),
                                                         eventPriority = 1))
```

## Validation runs

We can now run the validation:


```r
LandRBiomass_validation <- Cache(simInitAndSpades
                                 , times = validationTimes
                                 , params = validationParams
                                 , modules = "Biomass_validationKNN"
                                 , objects = validationObjects
                                 , outputs = validationOutputs
                                 , paths = validationPaths
                                 , .studyAreaName = SAname
                                 , userTags = c("validation", SAname)
                                 , cacheRepo = validationPaths$cachePath
                                 , omitArgs = c("userTags"))
saveSimList(LandRBiomass_validation, file.path(validationPaths$outputPath, paste0("simValid", runName, ".qs")))
```
