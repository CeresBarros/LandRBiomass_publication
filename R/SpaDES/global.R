## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES
## TESTS
##
## Ceres: Nov 2017
## ------------------------------------------------------

## clean workspace
rm(list=ls()); amc::.gc()

## requires as of June 10th 2019
# loading reproducible     0.2.8.9001
# loading quickPlot        0.1.6.9000
# loading SpaDES.core      0.2.5.9004
# loading SpaDES.tools     0.3.2.9000
# loading SpaDES.addins    0.1.2

library(SpaDES)
library(LandR)


## define paths
setPaths(modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))

## download modules
downloadModule()

## -----------------------------------------------
## SIMULATION SETUP
## -----------------------------------------------

## Set up sppEquiv  ---------------------------
data("sppEquivalencies_CA", package = "LandR")
sppEquivalencies_CA[grep("Pin", LandR), `:=`(EN_generic_short = "Pine",
                                             EN_generic_full = "Pine",
                                             Leading = "Pine leading")]

## define spp column to use for model
sppEquivCol <- "Boreal"
sppEquivalencies_CA <- na.omit(sppEquivalencies_CA, cols = sppEquivCol)

## create color palette for species used in model
sppColorVect <- sppColors(sppEquivalencies_CA, sppEquivCol,
                          newVals = "Mixed", palette = "Accent")
