# Floral color manuscript code

## Workflow:

*Note 1: Gathering of phylogeny and floral color data not included in this workflow*
*Note 2: To recreate results, download `floral_color_manuscript_data.zip` from releases, unzip, and run #8 from this workflow*

1. Download species occurrences from GBIF, clean them up, and format the data using `Distribution Downloader.Rmd`
2. Clean up and standardize taxonomical nomenclature using `Taxonomy Updater.R`
3. Create shapefiles with geometries for each species range based on occurrences using `Range Maker.R` (ignore `occs.R`)
4. Create shapefiles with rasters for each species range based on occurrences using `raster_maker.R`
5. Calculate range overlap for all species pairs using `overlap stuff again.R` or `Range Overlap Calculator.R` (depending on raster or geometry representation)
6. Gather and process phenology data using `Phenology.R`
7. Select focal sites for analysis with `site_picker.R`
8. Run analysis for each major question using `New Question 1.R`, `New Question 2.R`, and/or `New Question 3.R`

- To create a phylogeny of species with tips colored by floral color, use `megatree.R`
- To recreate other figures, use `Figures.R`
