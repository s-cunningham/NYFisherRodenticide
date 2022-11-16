
library(tidyverse)
library(blockCV)
library(sf)
library(raster)

### Read in raster layers ###
aea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

## NLCD
nlcd <- raster("data/rasters/nybuffnlcd.tif")

# Reclassify 0 (unclassified) and 128 (?) to no data
m <- rbind(c(0, NA), c(128, NA))
nlcd <- reclassify(nlcd, m)
unique(nlcd)

# Set up NLCD classes and class values
nlcd_values <- c(11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95)
nlcd_class <- c("Open Water", "Developed, Open Space", "Developed, Low Intensity", 
                "Developed, Medium Intensity", "Developed, High Intensity", "Barren", 
                "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub",
                "Grassland/Herbaceous", "Pasture/Hay", "Cultivated Crops", "Woody Wetlands", 
                "Emergent Herbaceous Wetlands")

# Add class names and numbers to the raster
levels(nlcd) <- list(data.frame(ID = nlcd_values,
                                landcov = nlcd_class))
nlcd <- projectRaster(nlcd, crs=aea)

## Wildland-urban interface
# 100m radius
wui100 <- raster("data/rasters/WUI100mNY.tif")
wui100 <- projectRaster(wui100, nlcd) # Match projection to nlcd

wui250 <- raster("data/rasters/WUI250mNY.tif")
wui250 <- projectRaster(wui250, nlcd) # Match projection to nlcd

wui500 <- raster("data/rasters/WUI500mNY.tif")
wui500 <- projectRaster(wui500, nlcd) # Match projection to nlcd

# Set up WUI classes and values
wui_values <- c(0, 1, 2)
wui_class <- c("not WUI", "intermix WUI", "interface WUI")

# Add class names and numbers to the raster
levels(wui100) <- list(data.frame(ID = wui_values,
                                  landcov = wui_class))

levels(wui250) <- list(data.frame(ID = wui_values,
                                  landcov = wui_class))

levels(wui500) <- list(data.frame(ID = wui_values,
                                  landcov = wui_class))


## Beech basal area
beech <- raster("data/rasters/baa250_masked.tif")
beech <- projectRaster(beech, nlcd)

## Combine rasters
layers <- stack(nlcd, wui100, wui250, wui500, beech)
layers <- brick(layers)

## Read in rodenticide locations
aea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat1 <- dat %>% filter(pt_index==1 & buffsize==30 & radius==100) %>% st_as_sf(coords=c("rand_x", "rand_y"), crs=crs(layers))

## Spatial clustering
sb <- spatialBlock(speciesData=dat1,
                   rasterLayer=layers,
                   rows=8, 
                   cols=8,
                   k=5,
                   selection="random")

## Evaluating models

# Loop over each point set
for (j in 1:10) {
  
  iter <- dat %>% filter(pt_index==j & buffsize==30 & radius==100)
    
  # Save folds (list)
  folds <- sb$folds
  
  # loop over folds
  for (i in seq_len(folds)) {
    
    trainSet <- iter[unlist(folds[[i]][1]),]
    testSet <- iter[unlist(folds[[i]][2])]
    
    
    
  }
  
  
  
}

