## load simList from simInit()
if (!exists("LandRBiomass_preSim"))
  LandRBiomass_preSim <- loadSimList(file.path(simPaths$outputPath, paste0("preSimList_", runName)))

## load the experiment simLists object
if (!exists("LandRBiomass_sim"))
  LandRBiomass_sim <- qs::qread(file.path(simPaths$outputPath, paste0("simList_LandRBiomass_sim_", runName)))

## get necessary rasters from simList (only rasters aren't kept in memory)
if (!inMemory(LandRBiomass_preSim$rawBiomassMap)) {
  rasFilename <- raster::filename(LandRBiomass_preSim$rawBiomassMap)
  rasFilename <- sub(inputPath(LandRBiomass_preSim), simPaths$inputPath, rasFilename)

  rawBiomassMap <- raster(rasFilename)
} else {
  rawBiomassMap <- LandRBiomass_preSim$rawBiomassMap
}

if (!inMemory(LandRBiomass_preSim$standAgeMap)) {
  rasFilename <- raster::filename(LandRBiomass_preSim$standAgeMap)
  rasFilename <- sub(inputPath(LandRBiomass_preSim), simPaths$inputPath, rasFilename)

  standAgeMap <- raster(rasFilename)
} else {
  standAgeMap <- LandRBiomass_preSim$standAgeMap
}

if (!inMemory(LandRBiomass_preSim$speciesLayers)) {
  speciesLayers <- lapply(unstack(LandRBiomass_preSim$speciesLayers), function(x) {
    rasFilename <- raster::filename(x)
    rasFilename <- sub(outputPath(LandRBiomass_preSim), simPaths$outputPath, rasFilename)
    layerName <- names(x)
    sppLayer <- raster(rasFilename)

    names(sppLayer) <- layerName
    return(sppLayer)
  })
  speciesLayers <- stack(speciesLayers)
} else {
  speciesLayers <- LandRBiomass_preSim$speciesLayers
}


if (!inMemory(LandRBiomass_preSim$rasterToMatch)) {
  rasFilename <- raster::filename(LandRBiomass_preSim$rasterToMatch)
  rasFilename <- sub(inputPath(LandRBiomass_preSim), simPaths$inputPath, rasFilename)

  standAgeMap <- raster(rasterToMatch)
} else {
  rasterToMatch <- LandRBiomass_preSim$rasterToMatch
}

rasterToMatch[!is.na(rasterToMatch)] <- 1

rstLCChangeAllbin <- Cache(prepInputs,
                           targetFile = "change_2001_2011.dat",
                           archive = "LCchange_AlldisturbancesBinary.zip",
                           url = "https://drive.google.com/file/d/1xd46zkJRdZVVNK93AUcfqUvtL-OFb4h4/view?usp=sharing",
                           alsoExtract = "similar",
                           fun = "raster::raster",
                           destinationPath = options()$reproducible.inputPaths,
                           studyArea = LandRBiomass_preSim$studyArea,
                           rasterToMatch = rasterToMatch,
                           useSAcrs = FALSE,
                           maskWithRTM = TRUE,
                           method = "ngb",
                           filename2 = file.path(options()$reproducible.inputPaths, "change_2001_2011_postProcess.tiff"),
                           overwrite = TRUE,
                           cacheRepo = simPaths$cachePath,
                           useCache = TRUE)

## make sure the extent matches the study area (it won't if using a smaller study area than RTMLarge)
rstLCChangeAllbin <- Cache(postProcess,
                           x = rstLCChangeAllbin,
                           studyArea = LandRBiomass_preSim$studyArea,
                           useSAcrs = FALSE,
                           filename2 = file.path(options()$reproducible.inputPaths, "change_2001_2011_postProcess.tiff"),
                           overwrite = TRUE,
                           cacheRepo = simPaths$cachePath,
                           useCache = TRUE)

## convert to mask
rstLCChangeAllbin[getValues(rstLCChangeAllbin) != 1] <- NA

## compile all simulation output tables and replace output paths
if (class(LandRBiomass_sim) == "simLists") {
simulationOutputs <- lapply(LandRBiomass_sim, FUN = function(x, localSimPaths) {
  oldPath <- dirname(outputPath(x)) ## exclude sim*_rep* folder
  DT <- as.data.table(outputs(x))
  DT[, file := sub(oldPath, localSimPaths$outputPath, file)]
  DT
  }, localSimPaths = simPaths)
simulationOutputs <- rbindlist(simulationOutputs)
} else {
  oldPath <- outputPath(LandRBiomass_sim) ## exclude sim*_rep* folder
  simulationOutputs <- as.data.table(outputs(LandRBiomass_sim))
  simulationOutputs[, file := sub(oldPath, simPaths$outputPath, file)]
  rm(oldPath)
}

## get biomassMap (should be the same across all reps)
biomassMap <- readRDS(simulationOutputs[objectName == "biomassMap", file][1])

