
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "Biomass_validationKNN",
  description = "Validation module for LandR Biomass predictions of forest succession", #"insert module description here",
  keywords = c("validation", "ecological simulation model",
               "forest dynamics", "forest succession", "data", "prediction"),
  authors = c(person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre"))),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.6", LandR = "0.0.2.9007",
                 Biomass_validationKNN = "0.0.1",
                 Biomass_core = "1.3.2", Biomass_speciesData = "1.0.0",
                 Biomass_borealDataPrep = "1.4.0.9000"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_validationKNN.Rmd"),
  reqdPkgs = list(),
  parameters = rbind(
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    expectsInput("firePerimeters", "shapefile",
                 desc = paste("A map of fire perimeters in the study area that can be used to exclude pixels",
                              "that have been burnt during the validation period. Defaults to the Canadian",
                              "Wildland Fire Information System 1986-2018 National Burned Area Composite."),
                 sourceURL = "http://cwfis.cfs.nrcan.gc.ca/downloads/nbac/nbac_1986_to_2018_20191017.zip"),
    expectsInput("rawBiomassMap", "RasterLayer",
                 desc = "total biomass raster layer in study area, default is Canada national biomass map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"),
    expectsInput("rstLCChange", "RasterLayer",
                 desc = paste("A map of land cover change types in the study area.",
                              "The default layer used, if not supplied, is Canada's forest change national map between 1985-2011 (CFS)",
                              "See https://opendata.nfis.org/mapserver/nfis-change_eng.html for more information."),
                 sourceURL = "https://opendata.nfis.org/downloads/forest_change/C2C_change_type_1985_2011.zip"),
    expectsInput("rstLCChangeYr", "RasterLayer",
                 desc = paste("A map of land cover change years in the study area.",
                              "The default layer used, if not supplied, is Canada's forest change year national map between 1985-2011 (CFS).",
                              "See https://opendata.nfis.org/mapserver/nfis-change_eng.html for more information."),
                 sourceURL = "https://opendata.nfis.org/downloads/forest_change/C2C_change_year_1985_2011.zip"),
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = NA, objectClass = NA, desc = NA)
  )
))

doEvent.Biomass_validationKNN = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "Biomass_validationKNN", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "Biomass_validationKNN", "save")
    },
    plot = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      #plotFun(sim) # uncomment this, replace with object to plot
      # schedule future event(s)

      # e.g.,
      #sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "Biomass_validationKNN", "plot")

      # ! ----- STOP EDITING ----- ! #
    },
    save = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "Biomass_validationKNN", "save")

      # ! ----- STOP EDITING ----- ! #
    },
    event1 = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "Biomass_validationKNN", "templateEvent")

      # ! ----- STOP EDITING ----- ! #
    },
    event2 = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "Biomass_validationKNN", "templateEvent")

      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

### template for save events
Save <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  #Plot(sim$object)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event1
Event1 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event1Test1 <- " this is test for event 1. " # for dummy unit test
  # sim$event1Test2 <- 999 # for dummy unit test


  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event2
Event2 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event2Test1 <- " this is test for event 2. " # for dummy unit test
  # sim$event2Test2 <- 777  # for dummy unit test


  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  browser()
  cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  if (getOption("LandR.verbose", TRUE) > 0)
    message(currentModule(sim), ": using dataPath '", dPath, "'.")

  ## Study area raster and shapefile -------------------------------------------
  if (!suppliedElsewhere("studyArea", sim)) {
    if (getOption("LandR.verbose", TRUE) > 0)
      message("'studyArea' was not provided by user. Using a polygon (6250000 m^2) in southwestern Alberta, Canada")
    sim$studyArea <- randomStudyArea(seed = 1234, size = (250^2)*100)
  }

  needRTM <- FALSE
  if (is.null(sim$rasterToMatch)) {
    if (!suppliedElsewhere("rasterToMatch", sim)) {
      needRTM <- TRUE
      message("There is no rasterToMatch supplied; will attempt to use rawBiomassMap")
    } else {
      stop("rasterToMatch is going to be supplied, but ", currentModule(sim), " requires it ",
           "as part of its .inputObjects. Please make it accessible to ", currentModule(sim),
           " in the .inputObjects by passing it in as an object in simInit(objects = list(rasterToMatch = aRaster)",
           " or in a module that gets loaded prior to ", currentModule(sim))
    }
  }

  if (!suppliedElsewhere("rawBiomassMap", sim) || needRTM) {
    sim$rawBiomassMap <- Cache(prepInputs,
                               targetFile = asPath(basename(rawBiomassMapFilename)),
                               archive = asPath(c("kNN-StructureBiomass.tar",
                                                  "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.zip")),
                               url = extractURL("rawBiomassMap"),
                               destinationPath = dPath,
                               studyArea = sim$studyArea,
                               rasterToMatch = if (!needRTM) sim$rasterToMatch else NULL,
                               # maskWithRTM = TRUE,    ## if RTM not supplied no masking happens (is this intended?)
                               maskWithRTM = if (!needRTM) TRUE else FALSE,
                               ## TODO: if RTM is not needed use SA CRS? -> this is not correct
                               # useSAcrs = if (!needRTM) TRUE else FALSE,
                               useSAcrs = FALSE,     ## never use SA CRS
                               method = "bilinear",
                               datatype = "INT2U",
                               filename2 = TRUE, overwrite = TRUE, userTags = cacheTags,
                               omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))
  }
  if (needRTM) {
    ## if we need rasterToMatch, that means a) we don't have it, but b) we will have rawBiomassMap
    sim$rasterToMatch <- sim$rawBiomassMap
    RTMvals <- getValues(sim$rasterToMatch)
    sim$rasterToMatch[!is.na(RTMvals)] <- 1

    sim$rasterToMatch <- Cache(writeOutputs, sim$rasterToMatch,
                               filename2 = file.path(cachePath(sim), "rasters", "rasterToMatch.tif"),
                               datatype = "INT2U", overwrite = TRUE,
                               userTags = cacheTags,
                               omitArgs = c("userTags"))

    ## this is old, and potentially not needed anymore
    if (FALSE) {
      studyArea <- sim$studyArea # temporary copy because it will be overwritten if it is suppliedElsewhere
      message("  Rasterizing the studyArea polygon map")
      if (!is(studyArea, "SpatialPolygonsDataFrame")) {
        dfData <- if (is.null(rownames(studyArea))) {
          polyID <- sapply(slot(studyArea, "polygons"), function(x) slot(x, "ID"))
          data.frame("field" = as.character(seq_along(length(studyArea))), row.names = polyID)
        } else {
          polyID <- sapply(slot(studyArea, "polygons"), function(x) slot(x, "ID"))
          data.frame("field" = rownames(studyArea), row.names = polyID)
        }
        studyArea <- SpatialPolygonsDataFrame(studyArea, data = dfData)
      }
      if (!identical(crs(studyArea), crs(sim$rasterToMatch))) {
        studyArea <- spTransform(studyArea, crs(sim$rasterToMatch))
        studyArea <- fixErrors(studyArea)
      }


      #TODO: review whether this is necessary (or will break LandWeb if removed) see Git Issue #22
      # layers provided by David Andison sometimes have LTHRC, sometimes LTHFC ... chose whichever
      LTHxC <- grep("(LTH.+C)", names(studyArea), value = TRUE)
      fieldName <- if (length(LTHxC)) {
        LTHxC
      } else {
        if (length(names(studyArea)) > 1) {
          ## study region may be a simple polygon
          names(studyArea)[1]
        } else NULL
      }

      sim$rasterToMatch <- crop(fasterizeFromSp(studyArea, sim$rasterToMatch, fieldName),
                                studyArea)
      sim$rasterToMatch <- Cache(writeRaster, sim$rasterToMatch,
                                 filename = file.path(dPath, "rasterToMatch.tif"),
                                 datatype = "INT2U", overwrite = TRUE,
                                 userTags = cacheTags,
                                 omitArgs = c("userTags"))
    }
  }


  if (!identical(crs(sim$studyArea), crs(sim$rasterToMatch))) {
    warning(paste0("studyArea and rasterToMatch projections differ.\n",
                   "studyArea will be projected to match rasterToMatch"))
    sim$studyArea <- spTransform(sim$studyArea, crs(sim$rasterToMatch))
    sim$studyArea <- fixErrors(sim$studyArea)
  }

  ## Land cover change (type and year) -------------------------------------------
  LCChangeFilename <- "C2C_change_type_1985_2011.tif"
  LCChangeYrFilename <- "C2C_change_year_1985_2011.tif"

  if (!suppliedElsewhere("rstLCChange", sim)) {
    sim$rstLCChange <- Cache(prepInputs,
                        targetFile = LCChangeFilename,
                        archive = asPath("C2C_change_type_1985_2011.zip"),
                        url = extractURL("rstLCChange"),
                        destinationPath = dPath,
                        studyArea = sim$studyArea,
                        rasterToMatch = sim$rasterToMatch,
                        maskWithRTM = TRUE,
                        method = "bilinear",
                        datatype = "INT2U",
                        filename2 = TRUE, overwrite = TRUE,
                        userTags = c("rstLCChange", cacheTags),
                        omitArgs = c("destinationPath", "targetFile", "userTags"))

    if (!identical(projection(sim$rstLCChange),
                   projection(sim$rasterToMatch)))
      projection(sim$rstLCChange) <- projection(sim$rasterToMatch) ## Ceres: this shouldn't be necessary anymore
  }

  if (!suppliedElsewhere("rstLCChangeYr", sim)) {
    sim$rstLCChangeYr <- Cache(prepInputs,
                             targetFile = LCChangeYrFilename,
                             archive = asPath("C2C_change_year_1985_2011.zip"),
                             url = extractURL("rstLCChangeYr"),
                             destinationPath = dPath,
                             studyArea = sim$studyArea,
                             rasterToMatch = sim$rasterToMatch,
                             maskWithRTM = TRUE,
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = TRUE, overwrite = TRUE,
                             userTags = c("rstLCChangeYr", cacheTags),
                             omitArgs = c("destinationPath", "targetFile", "userTags"))

    if (!identical(projection(sim$rstLCChangeYr),
                   projection(sim$rasterToMatch)))
      projection(sim$rstLCChangeYr) <- projection(sim$rasterToMatch) ## Ceres: this shouldn't be necessary anymore
  }

  ## Fire perimeter data ---------------------------------------------------
  firePerimetersFile <- "nbac_1986_to_2018_20191017.shp"

  if (!suppliedElsewhere("firePerimeters", sim)) {
    sim$firePerimeters <- Cache(prepInputs,
                               targetFile = firePerimetersFile,
                               archive = asPath("nbac_1986_to_2018_20191017.zip"),
                               url = extractURL("firePerimeters"),
                               destinationPath = dPath,
                               studyArea = sim$studyArea,
                               datatype = "INT2U",
                               filename2 = TRUE, overwrite = TRUE,
                               userTags = c("firePerimeters", cacheTags),
                               omitArgs = c("destinationPath", "targetFile", "userTags"))

    if (!identical(projection(sim$firePerimeters),
                   projection(sim$rasterToMatch)))
      projection(sim$firePerimeters) <- projection(sim$rasterToMatch) ## Ceres: this shouldn't be necessary anymore
  }

  ## Species equivalencies table -------------------------------------------
  if (!suppliedElsewhere("sppEquiv", sim)) {
    if (!is.null(sim$sppColorVect))
      stop("If you provide sppColorVect, you MUST also provide sppEquiv")

    data("sppEquivalencies_CA", package = "LandR", envir = environment())
    sim$sppEquiv <- as.data.table(sppEquivalencies_CA)
    ## By default, Abies_las is renamed to Abies_sp
    sim$sppEquiv[KNN == "Abie_Las", LandR := "Abie_sp"]

    ## check spp column to use
    if (P(sim)$sppEquivCol == "Boreal") {
      message(paste("There is no 'sppEquiv' table supplied;",
                    "will attempt to use species listed under 'Boreal'",
                    "in the 'LandR::sppEquivalencies_CA' table"))
    } else {
      if (grepl(P(sim)$sppEquivCol, names(sim$sppEquiv))) {
        message(paste("There is no 'sppEquiv' table supplied,",
                      "will attempt to use species listed under", P(sim)$sppEquivCol,
                      "in the 'LandR::sppEquivalencies_CA' table"))
      } else {
        stop("You changed 'sppEquivCol' without providing 'sppEquiv',",
             "and the column name can't be found in the default table ('LandR::sppEquivalencies_CA').",
             "Please provide conforming 'sppEquivCol', 'sppEquiv' and 'sppColorVect'")
      }
    }

    ## remove empty lines/NAs
    sim$sppEquiv <- sim$sppEquiv[!"", on = P(sim)$sppEquivCol]
    sim$sppEquiv <- na.omit(sim$sppEquiv, P(sim)$sppEquivCol)

    ## add default colors for species used in model
    sim$sppColorVect <- sppColors(sim$sppEquiv, P(sim)$sppEquivCol,
                                  newVals = "Mixed", palette = "Accent")
  } else {
    if (is.null(sim$sppColorVect))
      stop("If you provide 'sppEquiv' you MUST also provide 'sppColorVect'")
  }


  ## Species raster layers -------------------------------------------
  if (!suppliedElsewhere("speciesLayers", sim)) {
    #opts <- options(reproducible.useCache = "overwrite")
    sim$speciesLayers <- Cache(loadkNNSpeciesLayers,
                               dPath = dPath,
                               rasterToMatch = sim$rasterToMatch,
                               studyArea = sim$studyArea,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
                               sppEquiv = sim$sppEquiv,
                               knnNamesCol = "KNN",
                               sppEquivCol = P(sim)$sppEquivCol,
                               thresh = 5,
                               url = extractURL("speciesLayers"),
                               userTags = c(cacheTags, "speciesLayers"),
                               omitArgs = c("userTags"))
  }

  return(invisible(sim))
}
### add additional events as needed by copy/pasting from above
