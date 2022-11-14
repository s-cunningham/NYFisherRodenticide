

library(blockCV)
library(sf)
library(raster)


### Read in raster layers ###
## NLCD
nlcd <- raster("data/rasters/nybuffnlcd.tif")

# Reclassify 0 (unclassified) and 128 (?) to no data
m <- rbind(c(0, NA), c(128, NA))
nlcd <- classify(nlcd, m)
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

## Wildland-urban interface
# 100m radius
wui100 <- raster("data/rasters/WUI100mNY.tif")
wui100 <- project(wui100, nlcd) # Match projection to nlcd

wui250 <- raster("data/rasters/WUI250mNY.tif")
wui250 <- project(wui250, nlcd) # Match projection to nlcd

wui500 <- raster("data/rasters/WUI500mNY.tif")
wui500 <- project(wui500, nlcd) # Match projection to nlcd

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
beech <- rast("data/rasters/baa250_masked.tif")



