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
buff4_5 <- st_buffer(samples, 1196.83)
buff15 <- st_buffer(samples, 2185.1)
buff30 <- st_buffer(samples, 3090.19)

# Plot and example to see what 
ggplot() + 
  geom_sf(data=twmu) +
  geom_sf(data=samples, shape=20, color="blue", size=3) +
  # geom_sf(data=buff4_5, fill=NA, color="blue") +
  coord_sf(xlim=c(1583308.486, 1625741.123), ylim=c(861590.893, 888677.666)) +
  theme_bw()

# convert buffers to sf objects
buff4_5 <- st_as_sf(buff4_5)
buff15 <- st_as_sf(buff15)
buff30 <- st_as_sf(buff30)

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
landcov_fracs4_5 <- exact_extract(nlcd, buff4_5, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs4_5$buffsize <- 4.5

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
landcov_frac <- bind_rows(landcov_fracs4_5, landcov_fracs15, landcov_fracs30)

# Remove everything that is not forest or ag
keep_cov <- c(41, 42, 43, 81, 82)
landcov_frac <- landcov_frac[landcov_frac$value %in% keep_cov,]

write.csv(landcov_frac, "data/analysis-ready/nlcd_pct.csv")

## Read in WUI layers to calculate 
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

# Extract WUI
wui100_fracs4_5 <- exact_extract(wui100, buff4_5, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui100_fracs4_5$buffsize <- 4.5

wui100_sum4_5 <- exact_extract(rast, poly, function(values, coverage_fraction)
  sum(values * coverage_fraction, na.rm=TRUE))

# 250 m radius



# 500 m radius




## Read in housing layer to calculate density and extract values







