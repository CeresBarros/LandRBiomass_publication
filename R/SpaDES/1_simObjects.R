## -----------------------------------
## LOAD/MAKE NECESSARY OBJECTS
## -----------------------------------

## this script loads/treats/makes the necessary objects for the simulation
## it should be sourced before running any modules

## STUDY AREA(S) ---------------------------------------

## create a larger study area and create a smaller one (half extent)
## note that projection of the orignal CRS is always necessary
originalcrs <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
largeExtent <- extent(-104.757, -104.2197, 55.68663, 56.20319)
smallExtent <- largeExtent

smallExtent@xmax <- largeExtent@xmin + (largeExtent@xmax - largeExtent@xmin)*0.5   ## this is minimum size (10 000 pix at 250m res) -- need to increase large are to double
smallExtent@ymax <- largeExtent@ymin + (largeExtent@ymax - largeExtent@ymin)*0.5

studyAreaL <- as(largeExtent, "SpatialPolygons")
studyAreaL <-  SpatialPolygonsDataFrame(studyAreaL, data.frame(id = 1:length(studyAreaL)))
crs(studyAreaL) <- originalcrs
studyAreaL <- spTransform(studyAreaL, originalcrs)

studyAreaS <- as(smallExtent, "SpatialPolygons")
studyAreaS <-  SpatialPolygonsDataFrame(studyAreaS, data.frame(id = 1:length(studyAreaS)))
crs(studyAreaS) <- originalcrs
studyAreaS <- spTransform(studyAreaS, originalcrs)

## now reproject to Biomass_core standard
Biomass_corecrs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
studyAreaL <- spTransform(studyAreaL, Biomass_corecrs)
studyAreaS <- spTransform(studyAreaS, Biomass_corecrs)

plot(studyAreaL); plot(studyAreaS, add = TRUE, col = "red")

## Set up sppEquiv  ------------------------------
data("sppEquivalencies_CA", package = "LandR")
sppEquivalencies_CA[grep("Pin", LandR), `:=`(EN_generic_short = "Pine",
                                             EN_generic_full = "Pine",
                                             Leading = "Pine leading")]
sppEquivalencies_CA[grep("Betu_pap", LandR), `:=`(EN_generic_short = "Birch",
                                                  EN_generic_full = "Birch",
                                                  Leading = "Birch leading")]
## all Popu will be merged
sppEquivalencies_CA[grep("Popu_", LandR), `:=`(EN_generic_short = "Poplar",
                                               EN_generic_full = "Poplar",
                                               Leading = "Poplar leading")]

sppEquivalencies_CA[grep("Popu_", LandR), Boreal := "Popu_Spp"]

## define spp column to use for model
sppEquivCol <- "Boreal"
sppEquivalencies_CA <- sppEquivalencies_CA[!"", on = sppEquivCol]
sppEquivalencies_CA <- na.omit(sppEquivalencies_CA, sppEquivCol)

## create color palette for species used in model
sppColorVect <- sppColors(sppEquivalencies_CA, sppEquivCol,
                          newVals = "Mixed", palette = "Accent")


## Set up PSP data for LandR_speciesParameters  ------------------------------
PSPmeasure <- Cache(prepInputs,
                    targetFile = "PSPmeasure.Rdat",
                    archive = "randomized_LandR_speciesParameters_Inputs.zip",
                    url = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing",
                    destinationPath = simPaths$inputPath,
                    fun = "readRDS")

PSPplot <- Cache(prepInputs,
                    targetFile = "PSPplot.Rdat",
                    archive = "randomized_LandR_speciesParameters_Inputs.zip",
                    url = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing",
                    destinationPath = simPaths$inputPath,
                    fun = "readRDS")

PSPgis <- Cache(prepInputs,
                 targetFile = "randomizedPSPdata.Rdat",
                 archive = "randomized_LandR_speciesParameters_Inputs.zip",
                 url = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing",
                 destinationPath = simPaths$inputPath,
                 fun = "readRDS")
