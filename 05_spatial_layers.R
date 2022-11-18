## Set up spatial layers for rodenticide analysis
## 2022-06-16

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(landscapemetrics)

set.seed(123)

#### Read point locations ####
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
twmu <- st_read("data/spatial", "ARerase")
st_crs(twmu)
# Make sure projection is what it's supposed to be, and as a proj4 string
twmu <- st_transform(twmu, aea)

## Select random points for multiple imputation
# How many fishers per polygon
keycount <- loc %>% group_by(key) %>% count() %>% as.data.frame()

# Select random points per polygon
samples_per_polygon <- 10*keycount$n
samples <- st_sample(twmu, samples_per_polygon)
samples <- st_as_sf(samples) %>% 
  st_transform(crs=aea)
# st_write(samples, "data/spatial/random_samples.shp", layer_options="SHPT=POINT")

# Add names to points to associate with a liver ID
N.order <- order(loc$key)
loc <- loc[N.order,]
ids <- data.frame(id=rep(loc$RegionalID, each=10), buffno=rep(1:10, length(unique(loc$RegionalID))))
ids <- unite(ids, "id_index", 1:2, sep="_", remove=TRUE)

samples$name <- ids$id_index

pts <- st_coordinates(samples)
pts <- cbind(ids$id_index, pts) |> as.data.frame()
names(pts) <- c("pt_name", "x", "y")
# write_csv(pts, "output/random_point_locs.csv")

# Create buffer for 15km2 area
buff15 <- st_buffer(samples, 2185.1)
buff30 <- st_buffer(samples, 3090.19)
buff60 <- st_buffer(samples, 4370.194)

# Plot and example to see what 
ggplot() + 
  geom_sf(data=twmu) +
  geom_sf(data=samples, shape=20, color="blue", size=3) +
  # geom_sf(data=buff4p5, fill=NA, color="blue") +
  # geom_sf(data=buff60, fill=NA, color="green") +
  coord_sf(xlim=c(1583308.486, 1625741.123), ylim=c(861590.893, 888677.666)) +
  theme_bw()

#### Load NLCD layer ####
nlcd <- rast("data/rasters/nybuffnlcd.tif")

# convert buffers to sf objects
buff15 <- st_as_sf(buff15)%>% 
  st_transform(crs=crs(nlcd))
buff30 <- st_as_sf(buff30)%>% 
  st_transform(crs=crs(nlcd))
buff60 <- st_as_sf(buff60)%>% 
  st_transform(crs=crs(nlcd))

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

landcov_fracs60 <- exact_extract(nlcd, buff60, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs60$buffsize <- 60

# Combine into single data frame
landcov_frac <- bind_rows(landcov_fracs15, landcov_fracs30, landcov_fracs60)

# Remove everything that is not forest or ag
keep_cov <- c(41, 42, 43, 81, 82)
landcov_frac <- landcov_frac[landcov_frac$value %in% keep_cov,]

write_csv(landcov_frac, "data/analysis-ready/nlcd_pct.csv")

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

## 4.5 km2 buffer
wui100_fracs4p5 <- exact_extract(wui100, buff4p5, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui100_fracs4p5$buffsize <- 4.5

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

## 60 km2 buffer
wui100_fracs60 <- exact_extract(wui100, buff60, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui100_fracs60$buffsize <- 60

wui100_fracs <- rbind(wui100_fracs4p5, wui100_fracs15, wui100_fracs30, wui100_fracs60)
write_csv(wui100_fracs, "data/analysis-ready/wui100_frac.csv")


### 250 m radius

## 4.5 km2 buffer
wui250_fracs4p5 <- exact_extract(wui250, buff4p5, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui250_fracs4p5$buffsize <- 4.5

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

## 60 km2 buffer
wui250_fracs60 <- exact_extract(wui250, buff60, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui250_fracs60$buffsize <- 60

# Save as .csv
wui250_fracs <- rbind(wui250_fracs4p5, wui250_fracs15, wui250_fracs30, wui250_fracs60)
write_csv(wui250_fracs, "data/analysis-ready/wui250_frac.csv")

## 500 m radius

## 4.5 km2 buffer
wui500_fracs4p5 <- exact_extract(wui500, buff4p5, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui500_fracs4p5$buffsize <- 4.5

## 15 km2
wui500_fracs15 <- exact_extract(wui500, buff15, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui500_fracs15$buffsize <- 15

## 30 km2
wui500_fracs30 <- exact_extract(wui500, buff30, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui500_fracs30$buffsize <- 30

## 60 km2 radius
wui500_fracs60 <- exact_extract(wui500, buff60, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui500_fracs60$buffsize <- 60

# Save to csv
wui500_fracs <- rbind(wui500_fracs4p5, wui500_fracs15, wui500_fracs30, wui500_fracs60)
write_csv(wui500_fracs, "data/analysis-ready/wui500_frac.csv")

#### Read in predicted beech layer ####

## Load raster layers
beech <- rast("data/rasters/baa250_masked.tif")

# Do math to get total beech mass (sq.ft) instead of basal area per acre
# Cells are 65200 m^2, there are 15.444 acres in 65200 m^2
beech <- beech * 15.444

# annual beech mast index
mast <- read_csv("data/analysis-ready/ALTEMP26_beech-data.csv")
mast_mean <- mean(mast$Total_Beechnuts)
mast_median <- median(mast$Total_Beechnuts)
mast_max <- max(mast$Total_Beechnuts)
mast$devMean <- mast$Total_Beechnuts - mast_mean
mast$devMedian <- mast$Total_Beechnuts - mast_median
mast <- mast[mast$year>2016 & mast$year<=2020,]

# for each year, create a spatially-weighted beech layer
beech17 <- beech * mast$Total_Beechnuts[1]
beech18 <- beech * mast$Total_Beechnuts[2]
beech19 <- beech * mast$Total_Beechnuts[3]
beech20 <- beech * mast$Total_Beechnuts[4]

# put rasters in a list to loop over them
beech_list <- list(beech17, beech18, beech19, beech20)
years <- c(2017, 2018, 2019, 2020)

beech_sum <- data.frame()
for (i in 1:4) {
  beech_sum15 <- exact_extract(beech_list[[i]], buff15, 'sum')
  beech_sum15 <- data.frame(name=buff15$name, year=years[i], baa=beech_sum15, buffsize=15)
  
  beech_sum30 <- exact_extract(beech_list[[i]], buff30, 'sum')
  beech_sum30 <- data.frame(name=buff30$name, year=years[i], baa=beech_sum30, buffsize=30)
  
  beech_sum60 <- exact_extract(beech_list[[i]], buff60, 'sum')
  beech_sum60 <- data.frame(name=buff60$name, year=years[i], baa=beech_sum60, buffsize=60)
  
  beech_sum <- bind_rows(beech_sum, beech_sum15, beech_sum30, beech_sum60)
}

write_csv(beech_sum, "data/analysis-ready/baa_sum.csv")

#### Landscape metrics for forest cover ####

## Load rasters reclassified in ArcMap
tforest <- rast("data/rasters/nybuff_totalforest.tif")

nlcd_values <- c(1,2)

# Add class names and numbers to the raster
nlcd_class <- c("not forest", "total forest")
levels(tforest) <- list(data.frame(ID=nlcd_values, landcov=nlcd_class))

# Reproject samples
samples <- samples %>% st_transform(crs=crs(mixed))

## Landscape metrics...use points and have LSM create buffer
sizes <- c(4370.2, 6180.38, 8740.388)

# total forest
lsm_tforest_output <- sizes %>%
  set_names() %>%
  map_dfr(~sample_lsm(tforest, 
                      y=samples, 
                      plot_id=samples$name, 
                      what=c("lsm_c_ed",
                             "lsm_c_ai",
                             "lsm_c_contig_mn",
                             "lsm_c_cohesion",
                             "lsm_c_cpland",
                             "lsm_c_dcad",
                             "lsm_l_pladj"), 
                      shape="circle", size=.), .id="buffer")

# Save only metrics for total forest (class = s)
lsm_tforest_output <- lsm_tforest_output %>%
                        filter(class==2) %>%
                        select(plot_id, buffer, metric, value)

write_csv(lsm_tforest_output, "data/analysis-ready/forest_lsm.csv")


