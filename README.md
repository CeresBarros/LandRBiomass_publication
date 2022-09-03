# LandR Biomass Workflow

The LandR Biomass workflow is the implementation of a PERFICT modelling framework that links all steps associated with running a landscape dynamic vegetation model (data downloading, data treatment, parameterisation, calibration, simulation, model validation, visualisation and analysis of results, and model code testing) in a continuous and reproducible way.

We leverage several R packages (e.g. `SpaDES`, `reproducible`) to do this and present how we implemented this workflow with the LandR Biomass model, a re-implementation of LANDIS-II Biomass Succession Extension model (v3.2) using two examples.

For our examples we use a collection of several LandR Biomass modules (see [SpaDES modules](https://spades-core.predictiveecology.org/articles/i-introduction.html) for more information) each of which is developed collaboratively and has its own open git repository (each module folder is a *git submodule* in this repository). Code that is shared among modules was bundled into R packages (e.g. [`LandR` R package](https://github.com/PredictiveEcology/LandR/)), and hosted in open git repositories. 

If you want to learn more about SpaDES go to [https://spades.predictiveecology.org/](https://spades.predictiveecology.org/).

*Modules*
* [PredictiveEcology/Biomass_borealDataPrep](https://github.com/PredictiveEcology/Biomass_borealDataPrep)
* [PredictiveEcology/Biomass_core](https://github.com/PredictiveEcology/Biomass_core)
* [PredictiveEcology/Biomass_speciesData](https://github.com/PredictiveEcology/Biomass_speciesData)
* [PredictiveEcology/Biomass_validationKNN](https://github.com/PredictiveEcology/Biomass_validationKNN)

The easiest way to obtain all the code used in this workflow is to clone the main repository and each of the sub-modules (step 3 of installation notes below)

## Installation notes

**Disclaimer**
This project was run using R v4.0.3, tested on Windows 10 and 11. We cannot guarantee that package installation steps bellow and in the R scripts will work in other versions of R (or that all packages and package dependencies used here will remain available in CRAN for all eternity!).

1. Install development libraries

Some packages need to be built from source, which requires the appropriate development libraries for your operating system. You should install them first:

* *Windows*: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your R version.
* *macOS*: install Xcode commandline tools from the terminal: `xcode-select --install`.
* *Debian/Ubuntu Linux*: ensure `r-base-dev` is installed.

2. Getting the code:

To get this project's code, open a terminal window (e.g. from RStudio) and clone the project repository and its git submodules:

```bash
git clone --recurse-submodules "https://github.com/CeresBarros/LandRBiomass_publication" LandRBiomass_publication/
```

4. Install R packages and run simulations:

* follow R/SpaDES/global.Rmd to learn how to run an (example) simulation
* run R/SpaDES/global.R to reproduce all the scenario simulations (with replication) used in the publication

## Reporting bugs

Contact us via the project GitHub site: https://github.com/CeresBarros/LandRBiomass_publication/issues
