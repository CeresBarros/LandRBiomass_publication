## load simList from simInit()
LandRBiomass_simInit <- loadSimList(file.path(simPaths$outputPath, paste0("simInit", runName)))

## load the experiment simLists object
LandRBiomass_sim <- qs::qread(file.path(simPaths$outputPath, paste0("simList_LandRBiomass_sim_", runName)))

## get necessary rasters from simList (only rasters aren't kept in memory)
if (!inMemory(LandRBiomass_simInit$rawBiomassMap)) {
  rasFilename <- raster::filename(LandRBiomass_simInit$rawBiomassMap)
  rasFilename <- sub(inputPath(LandRBiomass_simInit), simPaths$inputPath, rasFilename)

  rawBiomassMap <- raster(rasFilename)
} else {
  rawBiomassMap <- LandRBiomass_simInit$rawBiomassMap
}

## get biomassMap
if (!inMemory(LandRBiomass_simInit$biomassMap)) {
  rasFilename <- raster::filename(LandRBiomass_simInit$biomassMap)
  rasFilename <- sub(ouputPath(LandRBiomass_simInit), simPaths$outputPath, rasFilename)

  biomassMap <- raster(rasFilename)
} else {
  biomassMap <- LandRBiomass_simInit$biomassMap
}

if (!inMemory(LandRBiomass_simInit$standAgeMap)) {
  rasFilename <- raster::filename(LandRBiomass_simInit$standAgeMap)
  rasFilename <- sub(inputPath(LandRBiomass_simInit), simPaths$inputPath, rasFilename)

  standAgeMap <- raster(rasFilename)
} else {
  standAgeMap <- LandRBiomass_simInit$standAgeMap
}

## note that species layers stack can be saved  as a multilayer raster obj (so stack("filename.grd") would work)
## HOWEVER because this may NOT be always the case, we're loading them separately anyways
if (!inMemory(LandRBiomass_simInit$speciesLayers)) {
  speciesLayers <- lapply(unstack(LandRBiomass_simInit$speciesLayers), function(x) {
    rasFilename <- raster::filename(x)
    rasFilename <- sub(outputPath(LandRBiomass_simInit), simPaths$outputPath, rasFilename)
    layerName <- names(x)
    sppLayer <- stack(rasFilename)[[layerName]]

    names(sppLayer) <- layerName
    return(sppLayer)
  })
  speciesLayers <- stack(speciesLayers)
} else {
  speciesLayers <- LandRBiomass_simInit$speciesLayers
}


if (!inMemory(LandRBiomass_simInit$rasterToMatch)) {
  rasFilename <- raster::filename(LandRBiomass_simInit$rasterToMatch)
  rasFilename <- sub(inputPath(LandRBiomass_simInit), simPaths$inputPath, rasFilename)

  standAgeMap <- raster(rasterToMatch)
} else {
  rasterToMatch <- LandRBiomass_simInit$rasterToMatch
}

rasterToMatch[!is.na(rasterToMatch)] <- 1

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


