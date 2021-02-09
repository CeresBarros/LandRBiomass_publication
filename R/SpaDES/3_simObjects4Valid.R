## load simList from simInit()
if (!exists("LandRBiomass_sim"))
  LandRBiomass_sim <- loadSimList(file.path(simPaths$outputPath, paste0("simInitList_", runName)))

## load the experiment simLists object
if (!exists("factorialSimulations"))
  factorialSimulations <- qs::qread(file.path(simPaths$outputPath, paste0("simList_factorialSimulations_", runName)))

## get necessary rasters from simList (only rasters aren't kept in memory)
if (!inMemory(LandRBiomass_sim$rawBiomassMap)) {
  rasFilename <- raster::filename(LandRBiomass_sim$rawBiomassMap)
  rasFilename <- sub(inputPath(LandRBiomass_sim), simPaths$inputPath, rasFilename)

  rawBiomassMap <- raster(rasFilename)
} else {
  rawBiomassMap <- LandRBiomass_sim$rawBiomassMap
}

if (!inMemory(LandRBiomass_sim$standAgeMap)) {
  rasFilename <- raster::filename(LandRBiomass_sim$standAgeMap)
  rasFilename <- sub(inputPath(LandRBiomass_sim), simPaths$inputPath, rasFilename)

  standAgeMap <- raster(rasFilename)
} else {
  standAgeMap <- LandRBiomass_sim$standAgeMap
}

if (!inMemory(LandRBiomass_sim$speciesLayers)) {
  speciesLayers <- lapply(unstack(LandRBiomass_sim$speciesLayers), function(x) {
    rasFilename <- raster::filename(x)
    rasFilename <- sub(inputPath(LandRBiomass_sim), simPaths$inputPath, rasFilename)
    layerName <- names(x)
    sppLayer <- raster(rasFilename)

    names(sppLayer) <- layerName
    return(sppLayer)
  })
  speciesLayers <- stack(speciesLayers)
} else {
  speciesLayers <- LandRBiomass_sim$speciesLayers
}


if (!inMemory(LandRBiomass_sim$rasterToMatch)) {
  rasFilename <- raster::filename(LandRBiomass_sim$rasterToMatch)
  rasFilename <- sub(inputPath(LandRBiomass_sim), simPaths$inputPath, rasFilename)

  standAgeMap <- raster(rasterToMatch)
} else {
  rasterToMatch <- LandRBiomass_sim$rasterToMatch
}

rasterToMatch[!is.na(rasterToMatch)] <- 1

rstLCChangeAllbin <- Cache(prepInputs,
                           targetFile = "change_2001_2011.dat",
                           archive = "LCchange_AlldisturbancesBinary.zip",
                           url = "https://drive.google.com/file/d/1xd46zkJRdZVVNK93AUcfqUvtL-OFb4h4/view?usp=sharing",
                           alsoExtract = "similar",
                           fun = "raster::raster",
                           destinationPath = options()$reproducible.inputPaths,
                           studyArea = LandRBiomass_sim$studyArea,
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
                           studyArea = LandRBiomass_sim$studyArea,
                           useSAcrs = FALSE,
                           filename2 = file.path(options()$reproducible.inputPaths, "change_2001_2011_postProcess.tiff"),
                           overwrite = TRUE,
                           cacheRepo = simPaths$cachePath,
                           useCache = TRUE)

## convert to mask
rstLCChangeAllbin[getValues(rstLCChangeAllbin) != 1] <- NA

## compile all simulation output tables and replace output paths
simulationOutputs <- lapply(factorialSimulations, FUN = function(x, localSimPaths) {
  oldPath <- dirname(outputPath(x)) ## exclude sim*_rep* folder
  DT <- as.data.table(outputs(x))
  DT[, file := sub(oldPath, localSimPaths$outputPath, file)]
  DT
  }, localSimPaths = simPaths)
simulationOutputs <- rbindlist(simulationOutputs)

## get biomassMap (should be the same across all reps)


