library(terra)
library(sf)
library(tidyverse)

# Read woods and woody wetlands raster
www <- rast("data/rasters/NYwoods_woodywetlands.tif")
plot(www)

# Reclassify 0 (unclassified) and 128 (?) to no data
m <- rbind(c(1, 1), c(128, NA))
www <- classify(www, m)
plot(www)

## Read shapefile
twmu <- st_read("data/spatial/AR_towns_WMUs.shp")
# plot(twmu)

# read streets
streets <- st_read("data/spatial/streets2kmbuff.shp")

# Intersect to clip
twmu <- st_intersection(twmu, streets)

# "Dissolve" because for some reason there are many polygon parts after the intersect
twmu <- twmu %>% group_by(key) %>% summarize()
plot(st_geometry(twmu))

twmu <- vect(twmu)
twmu <- project(twmu, crs(www))

## Extract raster by mask
www_ext <- mask(www, twmu)
plot(www_ext)

# Convert to polygon
www_ext <- as.polygons(www_ext)
plot(www_ext)

# Intersect twmu & forest/woody wetlands
twmu <- st_as_sf(twmu)
www_ext <- st_as_sf(www_ext)
twmu <- st_intersection(twmu, www_ext)
st_write(twmu, "data/spatial/www_AR_townWMUs.shp")
