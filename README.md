# LandR Biomass Workflow
The LandR Biomass workflow is the implementation of a framework linking all steps associated with a landscape dynamic vegetation model (data downloading, data treatment, parametrisation, calibration, simulation, model validation, visualisation and analysis of results and model code testing) in a continuous and reproducible way.
We leverage on several R packages (e.g. `SpaDES`, `reproducible`) to do this and present here an example of how we implemented this workflow with the LandR Biomass model, a re-implementation of LANDIS-II Biomass Succession Extension model (v3.2).

For our examples we use a collection of several LandR Biomass modules (see [SpaDES modules](https://spades-core.predictiveecology.org/articles/i-introduction.html) for more information) each of which is developed collaboratively and has it's own repository, and so are treated as *git submodules* in this repository. Code that is shared among modules was bundled into R packages (e.g. (`LandR` R package)[https://github.com/PredictiveEcology/LandR/]), and hosted in open git repositories. 

*Modules*
* (PredictiveEcology/Biomass_borealDataPrep)[https://github.com/PredictiveEcology/Biomass_borealDataPrep]
* (PredictiveEcology/Biomass_core)[https://github.com/PredictiveEcology/Biomass_core]
* (PredictiveEcology/Biomass_speciesData)[https://github.com/PredictiveEcology/Biomass_speciesData]
* (PredictiveEcology/Biomass_validationKNN)[https://github.com/PredictiveEcology/Biomass_validationKNN]

The easiest way to obtain all the code used in this workflow is to clone the main repository and each of the sub-modules (step 3 of installation notes below)

## Installation notes
The installing packages necessary to run the simulations requires installing development tools for R:

1. Install development libraries: building packages from source requires the appropriate development libraries for your operating system. See here for more details.

* *Windows*: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your R version.
* *macOS*: install Xcode commandline tools from the terminal: `xcode-select --install`.
* *Debian/Ubuntu Linux*: ensure `r-base-dev` is installed.

2. Install `devtools` package:
```
install.packages(devtools)
```

3. Getting the code
```
git clone --recurse-submodules -j8 https://github.com/PredictiveEcology/LandWeb
```

## Reporting bugs
Contact us via the package GitHub site: https://github.com/CeresBarros/LandRBiomass_publication
