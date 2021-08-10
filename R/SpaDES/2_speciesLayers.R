## -----------------------------------
## LOAD/MAKE SPECIES LAYERS
## -----------------------------------

## this script makes a pre-simulation object that makes species layers
## by running Biomass_speciesData This is the longest module to run and,
## unless the study area or the species needed change, it whould only
## be run once (even if other things change, like the simulation rep,
## or other modules). That's why caching is kept separate from the rest
## of the simulation

if (!runName %in% c("baseCase", "studyAreaChange", "studyAreaS", "studyAreaL")) {
  stop("runName must be one of 'baseCase', 'studyAreaChange', 'altParameters', 'studyAreaS' or 'studyAreaL'")
}

speciesPaths <- list(cachePath = file.path(simPaths$cachePath, "speciesLayers"),
                     modulePath = file.path("R/SpaDES/m"),
                     inputPath = file.path("R/SpaDES/inputs"),
                     outputPath = file.path(simPaths$outputPath, "speciesLayers"))

speciesParameters <- list(
  Biomass_speciesData = list(
    "types" = c("KNN")
    , "sppEquivCol" = sppEquivCol
    , ".useCache" = eventCaching
  )
)


speciesObjects <- list(
  "studyAreaLarge" = if (grepl("studyArea(S|L)$", runName)) get(runName) else studyAreaL
  , "sppEquiv" = sppEquivalencies_CA
  , "sppColorVect" = sppColorVect
)

SAname <- studyAreaName(speciesObjects$studyAreaLarge)
simOutSpeciesLayers <- Cache(simInitAndSpades
                             , times = list(start = 0, end = 1)
                             , params = speciesParameters
                             , modules = "Biomass_speciesData"
                             , objects = speciesObjects
                             , paths = speciesPaths
                             , debug = TRUE
                             , .plotInitialTime = NA
                             , .studyAreaName = SAname
                             , cacheRepo = speciesPaths$cachePath
                             , userTags = c("speciesLayersSim", SAname)
                             , omitArgs = c("userTags", ".plotInitialTime", "debug"))
saveSimList(simOutSpeciesLayers, file.path(simPaths$outputPath, paste0("simList_speciesLayers", runName)))   ## only save in first runs
