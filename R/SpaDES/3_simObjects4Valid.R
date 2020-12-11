RTM <- postProcess(simOutSpeciesLayers$rasterToMatchLarge,
                   studyArea = if (grepl("study", runName)) get(runName) else studyAreaS,
                   rasterToMatch = simOutSpeciesLayers$rasterToMatchLarge,
                   useSAcrs = FALSE,
                   maskWithRTM = FALSE,
                   useCache = FALSE,
                   method = "bilinear",
                   filename2 = NULL)
RTM[!is.na(RTM)] <- 1

rstLCChangeAllbin <- Cache(prepInputs,
                           targetFile = "change_2001_2011.dat",
                           archive = "LCchange_AlldisturbancesBinary.zip",
                           url = "https://drive.google.com/file/d/1xd46zkJRdZVVNK93AUcfqUvtL-OFb4h4/view?usp=sharing",
                           alsoExtract = "similar",
                           fun = "raster::raster",
                           destinationPath = options()$reproducible.inputPaths,
                           studyArea = if (grepl("study", runName)) get(runName) else studyAreaS,
                           rasterToMatch = RTM,
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
                           studyArea = if (grepl("study", runName)) get(runName) else studyAreaS,
                           useSAcrs = FALSE,
                           filename2 = file.path(options()$reproducible.inputPaths, "change_2001_2011_postProcess.tiff"),
                           overwrite = TRUE,
                           cacheRepo = simPaths$cachePath,
                           useCache = TRUE)

## convert to mask
rstLCChangeAllbin[getValues(rstLCChangeAllbin) != 1] <- NA

rm(RTM); amc::.gc()

