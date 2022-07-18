## Set up spatial layers for rodenticide analysis
## 2022-06-16

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)

set.seed(123)

## Read data frame with town/wmu location
loc <- read.csv("data/analysis-ready/ar_locations_only.csv")
loc <- loc[,-1]
loc <- loc[loc$RegionalID!="2018-9211",] # Seems wrong

## Read in polygon layer that has union of towns and WMUs
# twmu <- st_read("data/spatial", "WMUtown_union_Harv")
aea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
# twmu <- st_transform(twmu, aea)
# twmu <- unite(twmu, "key", 2:3, sep="-", remove=FALSE)
# Remove town/WMU combos that don't have a fisher liver from them
# twmu <- twmu[twmu$key %in% loc$key,]
# plot(st_geometry(twmu))

# Write to shapefile for erasing in ArcMap
# st_write(twmu, "output/AR_towns.shp", layer_options="SHPT=POLYGON")

## Read edited polygon layer back in
twmu <- st_read("data/spatial", "ARland")
st_crs(twmu)
# Make sure projection is what it's supposed to be, and as a proj4 string
twmu <- st_transform(twmu, aea)

## Select random points for multiple imputation
# How many fishers per polygon
keycount <- loc %>% group_by(key) %>% count() %>% as.data.frame()

# Select random points per polygon
samples_per_polygon <- 10*keycount$n
samples <- st_sample(twmu, samples_per_polygon)
samples <- st_as_sf(samples)

# Add names to points to associate with a liver ID
N.order <- order(loc$key)
loc <- loc[N.order,]
ids <- data.frame(id=rep(loc$RegionalID, each=10), buffno=rep(1:10, length(unique(loc$RegionalID))))
ids <- unite(ids, "id_index", 1:2, sep="_", remove=TRUE)

samples$name <- ids$id_index

# Create buffer for 15km2 area
buff15 <- st_buffer(samples, 2185.1)
buff30 <- st_buffer(samples, 3090.19)
buff60 <- st_buffer(samples, 4370.194)

# Plot and example to see what 
ggplot() + 
  geom_sf(data=twmu) +
  geom_sf(data=samples, shape=20, color="blue", size=3) +
  # geom_sf(data=buff60, fill=NA, color="blue") +
  coord_sf(xlim=c(1583308.486, 1625741.123), ylim=c(861590.893, 888677.666)) +
  theme_bw()

# convert buffers to sf objects
buff15 <- st_as_sf(buff15)
buff30 <- st_as_sf(buff30)
buff60 <- st_as_sf(buff60)

## Load raster layers
nlcd <- rast("data/rasters/nybuffnlcd.tif")

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

## Extract values from NLCD raster based on buffer using exactextractr
landcov_fracs60 <- exact_extract(nlcd, buff60, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs60$buffsize <- 60

landcov_fracs15 <- exact_extract(nlcd, buff15, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs15$buffsize <- 15

landcov_fracs30 <- exact_extract(nlcd, buff30, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs30$buffsize <- 30

# Combine into single data frame
landcov_frac <- bind_rows(landcov_fracs60, landcov_fracs15, landcov_fracs30)

# Remove everything that is not forest or ag
keep_cov <- c(41, 42, 43, 81, 82)
landcov_frac <- landcov_frac[landcov_frac$value %in% keep_cov,]

write.csv(landcov_frac, "data/analysis-ready/nlcd_pct.csv")

#### Read in WUI layers to calculate ####
# 100m radius
wui100 <- rast("data/rasters/WUI100mNY.tif")
wui100 <- project(wui100, nlcd) # Match projection to nlcd

wui250 <- rast("data/rasters/WUI250mNY.tif")
wui250 <- project(wui250, nlcd) # Match projection to nlcd

wui500 <- rast("data/rasters/WUI500mNY.tif")
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

#### Extract WUI

### 100m radius

## 60 km2 buffer
# Fraction
wui100_fracs60 <- exact_extract(wui100, buff60, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui100_fracs60$buffsize <- 60

## 15 km2 buffer
wui100_fracs15 <- exact_extract(wui100, buff15, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui100_fracs15$buffsize <- 15

## 30 km2 buffer
wui100_fracs30 <- exact_extract(wui100, buff30, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui100_fracs30$buffsize <- 30

wui100_fracs <- rbind(wui100_fracs60, wui100_fracs15, wui100_fracs30)
write.csv(wui100_fracs, "data/analysis-ready/wui100_frac.csv")


### 250 m radius
## 60 km2 buffer
wui250_fracs60 <- exact_extract(wui250, buff60, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui250_fracs60$buffsize <- 60

## 15 km2 buffer
wui250_fracs15 <- exact_extract(wui250, buff15, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui250_fracs15$buffsize <- 15

## 30 km2 buffer
wui250_fracs30 <- exact_extract(wui250, buff30, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui250_fracs30$buffsize <- 30

wui250_fracs <- rbind(wui250_fracs60, wui250_fracs15, wui250_fracs30)
write.csv(wui250_fracs, "data/analysis-ready/wui250_frac.csv")

## 500 m radius
wui500_fracs60 <- exact_extract(wui500, buff60, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui500_fracs60$buffsize <- 60

wui500_fracs15 <- exact_extract(wui500, buff15, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui500_fracs15$buffsize <- 15

wui500_fracs30 <- exact_extract(wui500, buff30, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui500_fracs30$buffsize <- 30
wui500_fracs <- rbind(wui500_fracs60, wui500_fracs15, wui500_fracs30)
write.csv(wui500_fracs, "data/analysis-ready/wui500_frac.csv")


#### Read in housing layer to calculate density and extract values ####
build100 <- rast("data/rasters/BuildCount100m.tif")
build100 <- project(build100, nlcd) # Match projection to nlcd

build250 <- rast("data/rasters/BuildCount250m.tif")
build250 <- project(build250, nlcd) # Match projection to nlcd

build500 <- rast("data/rasters/BuildCount500m.tif")
build500 <- project(build500, nlcd) # Match projection to nlcd

buff60 <- st_transform(buff60, crs(nlcd))
buff15 <- st_transform(buff15, crs(nlcd)) 
buff30 <- st_transform(buff30, crs(nlcd))

### Calculate average # buildings 
build100_avg60 <- exact_extract(build100, buff60, 'mean')
build100_avg60 <- data.frame(avg=build100_avg60)
build100_avg60$buffsize <- 60
build100_avg15 <- exact_extract(build100, buff15, 'mean')
build100_avg15 <- data.frame(avg=build100_avg15)
build100_avg15$buffsize <- 15
build100_avg30 <- exact_extract(build100, buff30, 'mean')
build100_avg30 <- data.frame(avg=build100_avg30)
build100_avg30$buffsize <- 30
build100_avg <- rbind(build100_avg60, build100_avg15, build100_avg30)
build100_avg$radius <- 100

build250_avg60 <- exact_extract(build250, buff60, 'mean')
build250_avg60 <- data.frame(avg=build250_avg60)
build250_avg60$buffsize <- 60
build250_avg15 <- exact_extract(build250, buff15, 'mean')
build250_avg15 <- data.frame(avg=build250_avg15)
build250_avg15$buffsize <- 15
build250_avg30 <- exact_extract(build250, buff30, 'mean')
build250_avg30 <- data.frame(avg=build250_avg30)
build250_avg30$buffsize <- 30
build250_avg <- rbind(build250_avg60, build250_avg15, build250_avg30)
build250_avg$radius <- 250

build500_avg60 <- exact_extract(build500, buff60, 'mean')
build500_avg60 <- data.frame(avg=build500_avg60)
build500_avg60$buffsize <- 60
build500_avg15 <- exact_extract(build500, buff15, 'mean')
build500_avg15 <- data.frame(avg=build500_avg15)
build500_avg15$buffsize <- 15
build500_avg30 <- exact_extract(build500, buff30, 'mean')
build500_avg30 <- data.frame(avg=build500_avg30)
build500_avg30$buffsize <- 30
build500_avg <- rbind(build500_avg60, build500_avg15, build500_avg30)
build500_avg$radius <- 500

build_avg <- rbind(build100_avg, build250_avg, build500_avg)
build_avg$name <- rep(buff60$name, 9)
write.csv(build_avg, "data/analysis-ready/build_avg.csv")

#### Read in predicted beech layer ####

## Load raster layers
beech <- rast("data/rasters/baa250_masked.tif")

beech_mean15 <- exact_extract(beech, buff15, 'mean')
beech_mean15 <- data.frame(name=buff15$name, baa=beech_mean15, buffsize=15)

beech_mean30 <- exact_extract(beech, buff30, 'mean')
beech_mean30 <- data.frame(name=buff30$name, baa=beech_mean30, buffsize=30)

beech_mean60 <- exact_extract(beech, buff60, 'mean')
beech_mean60 <- data.frame(name=buff60$name, baa=beech_mean60, buffsize=60)

beech_mean <- bind_rows(beech_mean15, beech_mean30, beech_mean60)
write.csv(beech_mean, "data/analysis-ready/baa_mean.csv")


