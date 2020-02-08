RTM <- postProcess(simOutSpeciesLayers$rasterToMatchLarge,
                   studyArea = if (grepl("study", runName)) get(runName) else studyAreaS,
                   rasterToMatch = simOutSpeciesLayers$rasterToMatchLarge,
                   useSAcrs = FALSE,
                   maskWithRTM = FALSE,
                   useCache = FALSE,
                   method = "bilinear",
                   filename2 = NULL)
RTM[!is.na(RTM)] <- 1

rstLCChangeAllbin <- unzip(file.path(options()$reproducible.inputPaths, "LCchange_AlldisturbancesBinary.zip"),
                           exdir = file.path(options()$reproducible.inputPaths))
rstLCChangeAllbin <- raster(file.path(options()$reproducible.inputPaths, "change_2001_2011.dat"))

rstLCChangeAllbin <- postProcess(x = rstLCChangeAllbin,
                                 studyArea = if (grepl("study", runName)) get(runName) else studyAreaS,
                                 rasterToMatch = RTM,
                                 useSAcrs = FALSE,
                                 maskWithRTM = TRUE,
                                 method = "ngb",
                                 filename2 = file.path(options()$reproducible.inputPaths, "change_2001_2011_postProcess.tiff"),
                                 overwrite = TRUE,
                                 cacheRepo = simPaths$cachePath,
                                 useCache = TRUE)

## convert to mask
rstLCChangeAllbin[getValues(rstLCChangeAllbin) != 1] <- NA

rm(RTM); amc::.gc()

