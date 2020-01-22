## Create simLinks for large folders

## move folders to storage
# rsync -a /home/cbarros/GitHub/LandRBiomass_publication/data /mnt/storage/cbarros/LandRBiomass_publication
# rsync -a /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/cache /mnt/storage/cbarros/LandRBiomass_publication/R/SpaDES/
# rsync -a /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/inputs /mnt/storage/cbarros/LandRBiomass_publication/R/SpaDES/
# rsync -a /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/outputs /mnt/storage/cbarros/LandRBiomass_publication/R/SpaDES/
# rsync -a /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/m/Biomass_core/data /mnt/storage/cbarros/LandRBiomass_publication/R/SpaDES/m/Biomass_core/
# rsync -a /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/m/Biomass_borealDataPrep/data /mnt/storage/cbarros/LandRBiomass_publication/R/SpaDES/m/Biomass_borealDataPrep/
# rsync -a /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/m/Biomass_speciesData/data /mnt/storage/cbarros/LandRBiomass_publication/R/SpaDES/m/Biomass_speciesData/

## after checking all is good remove the folders before creating the links
# rm -r /home/cbarros/GitHub/LandRBiomass_publication/data
# rm -r /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/cache
# rm -r /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/inputs
# rm -r /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/outputs
# rm -r /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/m/Biomass_core/data
# rm -r /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/m/Biomass_borealDataPrep/data
# rm -r /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/m/Biomass_speciesData/data

## create simLinks
ln -s /mnt/storage/cbarros/LandRBiomass_publication/data /home/cbarros/GitHub/LandRBiomass_publication/data
ln -s /mnt/storage/cbarros/LandRBiomass_publication/R/SpaDES/cache /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/cache
ln -s /mnt/storage/cbarros/LandRBiomass_publication/R/SpaDES/inputs /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/inputs
ln -s /mnt/storage/cbarros/LandRBiomass_publication/R/SpaDES/outputs /home/cbarros/GitHub/LandRBiomass_publication/R/SpaDES/outputs

