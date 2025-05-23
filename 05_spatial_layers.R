## Set up spatial layers for rodenticide analysis
## 2022-06-16

library(tidyverse)
library(sf)
library(terra)
library(exactextractr)
library(landscapemetrics)

set.seed(1)

# Function
add_sub_key <- function(df) {
  df$key_sub <- 1
  key <- unique(df$key)
  for (i in 1:length(key)) {
    x <- df[df$key==key[i],]
    nrowx <- nrow(x)
    if (nrowx>10) {
      nsets <- nrowx/10
      df$key_sub[df$key==key[i]] <- rep(1:nsets, each=10)
    }
  }
  return(df)
}

add_sub_key2 <- function(df) {
  df$key_sub <- 1
  key <- unique(df$key)
  for (i in 1:length(key)) {
    x <- df[df$key==key[i],]
    nrowx <- nrow(x)
    if (nrowx>1) {
      nsets <- nrowx
      df$key_sub[df$key==key[i]] <- rep(1:nsets)
    }
  }
  return(df)
}

#### Read point locations ####
## Read data frame with town-WMU key
loc <- read_csv("output/ncompounds_trace.csv")
loc <- loc[,-1]

## Read in polygon layer that has union of towns and WMUs
aea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

## Read edited polygon layer back in
twmu <- st_read("data/spatial", "www_AR_townWMUs")
st_crs(twmu)

# Make sure projection is what it's supposed to be, and as a proj4 string
twmu <- st_transform(twmu, aea)

# Remove unneeded columns from twmu
twmu <- twmu %>% select(key, geometry) 

# calculate area of town/WMU unioned polygons
sqm_ar <- st_area(twmu)
sqkm_ar <- sqm_ar / 1e6

sort(sqkm_ar)

mean(sqkm_ar)
sd(sqkm_ar)

## Select random points for multiple imputation
# How many fishers per polygon
keycount <- loc %>% group_by(key) %>% 
  count() %>% as.data.frame()

# Reorder according to twmu 
first <- unique(twmu$key)
second <- unique(keycount$key)

reorder_idx <- match(first, second)

keycount <- keycount[reorder_idx,] 
keycount <- keycount[complete.cases(keycount),]

# Select random points per polygon
samples_per_polygon <- 10*keycount$n
samples <- st_sample(twmu, samples_per_polygon)
samples <- st_as_sf(samples) %>% 
  st_transform(crs=aea)
samples$key <- map2(keycount$key, keycount$n*10, rep) %>% unlist()
# st_write(samples, "data/spatial/random_samples.shp", layer_options="SHPT=POINT", append=FALSE)

# Add names to points to associate with a liver ID
sdf <- st_coordinates(samples)
sdf <- bind_cols(sdf, map2(keycount$key, keycount$n*10, rep) %>% unlist()) %>%
  rename(key=`...3`, x=X, y=Y) %>%
  select(key, x, y)
sdf$pt_index <- rep(1:10, 338)
sdf <- add_sub_key(sdf)


loc <- add_sub_key2(loc)
loc <- left_join(loc, sdf, by=c("key", "key_sub"))

ggplot(loc, aes(x=x, y=y, color=factor(Region))) + geom_point()

# convert back to sf
loc <- loc %>% unite("name", c(1,12), sep="_", remove=FALSE) %>%
  select(RegionalID, name, pt_index, key, key_sub, year, Region, x, y)
write_csv(loc, "output/random_point_locs.csv")
samples <- st_as_sf(loc, coords=c("x","y"), crs=aea)
# st_write(samples, "data/spatial/df_random_samples.shp", layer_options="SHPT=POINT", append=FALSE)

# Create buffers
buff15 <- st_buffer(samples, 2185.0969) # 15: 2185.0969
buff30 <- st_buffer(samples, 3090.1936) # 30: 3090.1936
buff45 <- st_buffer(samples, 3784.6988)
# buff45 <- st_buffer(samples, 4370.194)

# Plot and example to see what 
# ggplot() +
#   geom_sf(data=twmu, aes(color=key, fill=key)) +
#   geom_sf(data=samples, shape=20, color="blue", size=3) +
#   geom_sf(data=buff15, fill=NA, color="blue") +
#   geom_sf(data=buff45, fill=NA, color="green") +
#   # coord_sf(xlim=c(1583308.486, 1625741.123), ylim=c(861590.893, 888677.666)) +
#   coord_sf(xlim=c(1668479, 1719120), ylim=c(819906, 894757)) +
#   theme_bw() +
#   theme(legend.position="none")

#### Load NLCD layer ####
nlcd <- rast("data/rasters/nybuffnlcd.tif")

# convert buffers to sf objects
buff15 <- st_as_sf(buff15)%>% 
  st_transform(crs=crs(nlcd))
buff30 <- st_as_sf(buff30)%>% 
  st_transform(crs=crs(nlcd))
buff45 <- st_as_sf(buff45)%>% 
  st_transform(crs=crs(nlcd))

# st_write(buff45, "data/spatial/AR_locations_45km_buffer.shp")


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

landcov_fracs45 <- exact_extract(nlcd, buff45, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs45$buffsize <- 45

# Combine into single data frame
landcov_frac <- bind_rows(landcov_fracs15, landcov_fracs30, landcov_fracs45)

# Remove everything that is not forest or ag
keep_cov <- c(41, 42, 43, 81, 82)
landcov_frac <- landcov_frac[landcov_frac$value %in% keep_cov,]

write_csv(landcov_frac, "data/analysis-ready/nlcd_pct.csv")

#### WUI layer ####

#### Read in WUI layers to calculate ####
wui500 <- rast("data/rasters/WUI500mNY.tif")
wui500 <- project(wui500, nlcd) # Match projection to nlcd

# Set up WUI classes and values
wui_values <- c(0, 1, 2)
wui_class <- c("not WUI", "intermix WUI", "interface WUI")

levels(wui500) <- list(data.frame(ID = wui_values,
                                  landcov = wui_class))

## Extract WUI 
## 15 km2
wui500_fracs15 <- exact_extract(wui500, buff15, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
wui500_fracs15$buffsize <- 15

# Save to csv
write_csv(wui500_fracs15, "data/analysis-ready/wui500_frac15.csv")

#### Read in predicted beech layer ####

## Load raster layers
beech <- rast("data/rasters/baa250_masked.tif")

# Do math to get total beech mass (sq.ft) instead of basal area per acre
# Cells are 65200 m^2, there are 15.444 acres in 65200 m^2
beech <- beech * 15.444

# Extract sum beech mast
beech_sum15 <- exact_extract(beech, buff15, 'sum')
beech_sum15 <- data.frame(name=buff15$name, baa=beech_sum15, buffsize=15)

beech_sum30 <- exact_extract(beech, buff30, 'sum')
beech_sum30 <- data.frame(name=buff30$name, baa=beech_sum30, buffsize=30)

beech_sum45 <- exact_extract(beech, buff45, 'sum')
beech_sum45 <- data.frame(name=buff45$name, baa=beech_sum45, buffsize=45)

beech_sum_single <- bind_rows(beech_sum15, beech_sum30, beech_sum45)
write_csv(beech_sum_single, "data/analysis-ready/baa_sum_single_raster.csv")

# annual beech mast index
mast <- read_csv("data/analysis-ready/ALTEMP26_beech-data.csv")
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
  
  beech_sum45 <- exact_extract(beech_list[[i]], buff45, 'sum')
  beech_sum45 <- data.frame(name=buff45$name, year=years[i], baa=beech_sum45, buffsize=45)
  
  beech_sum <- bind_rows(beech_sum, beech_sum15, beech_sum30, beech_sum45)
}

write_csv(beech_sum, "data/analysis-ready/baa_sum.csv")

#### Building layer raster ####
build_cntroid <- rast("data/rasters/NewYork_centroids.tif")
build_cntroid <- project(build_cntroid, nlcd) # Match projection to nlcd

build_sum15 <- exact_extract(build_cntroid, buff15, 'sum')
build_sum15 <- data.frame(name=buff15$name, nbuildings=build_sum15, buffsize=15)

build_sum30 <- exact_extract(build_cntroid, buff30, 'sum')
build_sum30 <- data.frame(name=buff30$name, nbuildings=build_sum30, buffsize=30)

build_sum45 <- exact_extract(build_cntroid, buff45, 'sum')
build_sum45 <- data.frame(name=buff45$name, nbuildings=build_sum45, buffsize=45)

build_sum <- bind_rows(build_sum15, build_sum30, build_sum45)
write_csv(build_sum, "data/analysis-ready/building-centroid_sum.csv")


#### Stand age ####
stnd <- rast("data/rasters/stand_age_buffer10km.tif")
stnd <- project(stnd, nlcd) # Match projection to nlcd

# mean
stnd_mean15 <- exact_extract(stnd, buff15, 'mean')
stnd_mean15 <- data.frame(name=buff15$name, stand_age_mean=stnd_mean15, buffsize=15)

stnd_mean30 <- exact_extract(stnd, buff30, 'mean')
stnd_mean30 <- data.frame(name=buff30$name, stand_age_mean=stnd_mean30, buffsize=30)

stnd_mean45 <- exact_extract(stnd, buff45, 'mean')
stnd_mean45 <- data.frame(name=buff45$name, stand_age_mean=stnd_mean45, buffsize=45)

stnd_mean <- bind_rows(stnd_mean15, stnd_mean30, stnd_mean45)
write_csv(stnd_mean, "data/analysis-ready/stand-age_mean.csv")

# standard deviation
stnd_stdev15 <- exact_extract(stnd, buff15, 'stdev')
stnd_stdev15 <- data.frame(name=buff15$name, stand_age_sd=stnd_stdev15, buffsize=15)

stnd_stdev30 <- exact_extract(stnd, buff30, 'stdev')
stnd_stdev30 <- data.frame(name=buff30$name, stand_age_sd=stnd_stdev30, buffsize=30)

stnd_stdev45 <- exact_extract(stnd, buff45, 'stdev')
stnd_stdev45 <- data.frame(name=buff45$name, stand_age_sd=stnd_stdev45, buffsize=45)

stnd_stdev <- bind_rows(stnd_stdev15, stnd_stdev30, stnd_stdev45)
write_csv(stnd_stdev, "data/analysis-ready/stand-age_stdev.csv")

#### Existing vegetation cover
evt <- rast("data/rasters/LANDFIRE/ny10kmbuff_lf2020_evt_220.tif")
evt <- as.factor(evt)

evt_count15 <- exact_extract(evt, buff15, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    reframe(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
# evt_count15 <- evt_count15 %>% mutate(freq=freq/min(freq)) %>%
#                 group_by(name) %>% slice_max(freq,n=2) %>%
#                 group_by(value) %>% count() %>%
#                 arrange(desc(n))


evt_count30 <- exact_extract(evt, buff30, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    reframe(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
# evt_count30 <- evt_count30 %>% mutate(freq=freq/min(freq)) %>%
#   group_by(name) %>% slice_max(freq,n=2) %>%
#   group_by(value) %>% count() %>%
#   arrange(desc(n))

evt_count45 <- exact_extract(evt, buff45, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    reframe(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
# evt_count45 <- evt_count45 %>% mutate(freq=freq/min(freq)) %>%
#   group_by(name) %>% slice_max(freq,n=2) %>%
#   group_by(value) %>% count() %>%
#   arrange(desc(n))

# evt_tab <- read_csv("output/EVT_LF220_table.csv")

evt_count15 <- evt_count15 %>% filter(value==7302 | value==7511 | value==7977 | value==9312 | value==7512 |
                                        value==7366 | value==9315 | value==7373) %>% mutate(buffsize=15)
evt_count30 <- evt_count30 %>% filter(value==7302 | value==7511 | value==7977 | value==9312 | value==7512 |
                                        value==7366 | value==9315 | value==7373) %>% mutate(buffsize=30)
evt_count45 <- evt_count45 %>% filter(value==7302 | value==7511 | value==7977 | value==9312 | value==7512 |
                                        value==7366 | value==9315 | value==7373) %>% mutate(buffsize=45)

evt_count <- bind_rows(evt_count15, evt_count30, evt_count45)
write_csv(evt_count, "data/analysis-ready/evt_pct_filtered.csv")

## Canopy cover
cc <- rast("data/rasters/nlcd21_tcc_ny10kbuff.tif")
cc <- as.numeric(cc)

# convert buffers to sf objects
buff15 <- st_as_sf(buff15)%>% 
  st_transform(crs=crs(cc))
buff30 <- st_as_sf(buff30)%>% 
  st_transform(crs=crs(cc))
buff45 <- st_as_sf(buff45)%>% 
  st_transform(crs=crs(cc))

# mean
cc_mean15 <- exact_extract(cc, buff15, 'mean')
cc_mean15 <- data.frame(name=buff15$name,cc_mean=cc_mean15, buffsize=15)

cc_mean30 <- exact_extract(cc, buff30, 'mean')
cc_mean30 <- data.frame(name=buff30$name, cc_mean=cc_mean30, buffsize=30)

cc_mean45 <- exact_extract(cc, buff45, 'mean')
cc_mean45 <- data.frame(name=buff45$name, cc_mean=cc_mean45, buffsize=45)

cc_mean <- bind_rows(cc_mean15, cc_mean30, cc_mean45)
write_csv(cc_mean, "data/analysis-ready/canopy-cover_mean.csv")

# standard deviation
cc_stdev15 <- exact_extract(cc, buff15, 'stdev')
cc_stdev15 <- data.frame(name=buff15$name, cc_sd=cc_stdev15, buffsize=15)

cc_stdev30 <- exact_extract(cc, buff30, 'stdev')
cc_stdev30 <- data.frame(name=buff30$name, cc_sd=cc_stdev30, buffsize=30)

cc_stdev45 <- exact_extract(cc, buff45, 'stdev')
cc_stdev45 <- data.frame(name=buff45$name, cc_sd=cc_stdev45, buffsize=45)

cc_stdev <- bind_rows(cc_stdev15, cc_stdev30, cc_stdev45)
write_csv(cc_stdev, "data/analysis-ready/canopy-cover_stdev.csv")


#### Landscape metrics for forest cover ####

## Load rasters reclassified in ArcMap
tforest <- rast("data/rasters/allNLCD_reclassFOR.tif")
tforest <- project(tforest, nlcd) # Match projection to nlcd

nlcd_values <- c(1,2)

# Add class names and numbers to the raster
nlcd_class <- c("not forest", "forest")
levels(tforest) <- list(data.frame(ID=nlcd_values, landcov=nlcd_class))

# Reproject samples
samples <- samples %>% st_transform(crs=crs(tforest))

## Landscape metrics...use points and have LSM create buffer
sizes <- c(2185.0969, 3090.1936, 3784.6988)  # 15, 30, 45

# total forest
lsm_tforest_output <- sizes %>%
  set_names() %>%
  map_dfr(~sample_lsm(tforest, 
                      y=samples, 
                      plot_id=samples$name, 
                      what=c("lsm_c_ed"),
                      shape="circle", size=.), .id="buffer")

# Save only metrics for total forest (class = s)
lsm_tforest_output <- lsm_tforest_output %>%
  filter(class==2) %>%
  select(plot_id, buffer, metric, value)
  
write_csv(lsm_tforest_output, "data/analysis-ready/forest_edge_density.csv")


#### Load Cropland data layer ####
cldl18 <- rast("F:/PhDwork/map_data/NASSCropScape/CropScape2018_NY.tif")
cldl19 <- rast("F:/PhDwork/map_data/NASSCropScape/CropScape2019_NY.tif")
cldl20 <- rast("F:/PhDwork/map_data/NASSCropScape/CropScape2020_NY.tif")

un.cldl <- unique(cldl18$CropScape2018_NY)
woody <- c(66,67,68,69,70,71,74,75,76,77)
not_woody <- un.cldl$CropScape2018_NY[!(un.cldl$CropScape2018_NY %in% woody)]
not_woody <- not_woody[2:length(not_woody)]
# Reclassify 0 (unclassified) and 128 (?) to no data
m <- matrix(data=c(0,woody,not_woody, NA, rep(1, length(woody)), rep(2, length(not_woody))),
            ncol=2,nrow=69, byrow=FALSE)
cldl18<- classify(cldl18, m)
unique(cldl18)

cldl19 <- classify(cldl19, m)
cldl20 <- classify(cldl20, m)

# Set up NLCD classes and class values
cldl_values <- c(1, 2)
cldl_class <- c("Tree Crops", "Other")

# Add class names and numbers to the raster
levels(cldl18) <- list(data.frame(ID = cldl_values,
                                landcov = cldl_class))
levels(cldl19) <- list(data.frame(ID = cldl_values,
                                  landcov = cldl_class))
levels(cldl20) <- list(data.frame(ID = cldl_values,
                                  landcov = cldl_class))

buff15 <- st_as_sf(buff15)%>% 
  st_transform(crs=crs(cldl18))
buff30 <- st_as_sf(buff30)%>% 
  st_transform(crs=crs(cldl18))
buff45 <- st_as_sf(buff45)%>% 
  st_transform(crs=crs(cldl18))

## Extract values from Cropland data layer raster based on buffer using exactextractr
landcov_fracs15 <- exact_extract(cldl18, buff15, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs15$buffsize <- 15

landcov_fracs30 <- exact_extract(cldl18, buff30, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs30$buffsize <- 30

landcov_fracs45 <- exact_extract(cldl18, buff45, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs45$buffsize <- 45

# Combine into single data frame
landcov_frac <- bind_rows(landcov_fracs15, landcov_fracs30, landcov_fracs45)

# Remove everything that is not forest or ag
keep_cov <- c(1)
landcov_frac <- landcov_frac[landcov_frac$value %in% keep_cov,]

write_csv(landcov_frac, "data/analysis-ready/cldl_woody_pct_2018.csv")

## Extract values from Cropland data layer raster based on buffer using exactextractr
landcov_fracs15 <- exact_extract(cldl19, buff15, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs15$buffsize <- 15

landcov_fracs30 <- exact_extract(cldl19, buff30, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs30$buffsize <- 30

landcov_fracs45 <- exact_extract(cldl19, buff45, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs45$buffsize <- 45

# Combine into single data frame
landcov_frac <- bind_rows(landcov_fracs15, landcov_fracs30, landcov_fracs45)

# Remove everything that is not forest or ag
keep_cov <- c(1)
landcov_frac <- landcov_frac[landcov_frac$value %in% keep_cov,]

write_csv(landcov_frac, "data/analysis-ready/cldl_woody_pct_2019.csv")

## Extract values from Cropland data layer raster based on buffer using exactextractr
landcov_fracs15 <- exact_extract(cldl20, buff15, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs15$buffsize <- 15

landcov_fracs30 <- exact_extract(cldl20, buff30, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs30$buffsize <- 30

landcov_fracs45 <- exact_extract(cldl20, buff45, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)
landcov_fracs45$buffsize <- 45

# Combine into single data frame
landcov_frac <- bind_rows(landcov_fracs15, landcov_fracs30, landcov_fracs45)

# Remove everything that is not forest or ag
keep_cov <- c(1)
landcov_frac <- landcov_frac[landcov_frac$value %in% keep_cov,]

write_csv(landcov_frac, "data/analysis-ready/cldl_woody_pct_2020.csv")


ggplot(landcov_frac) +
  geom_density(aes(x=freq, color=factor(buffsize), fill=factor(buffsize)), alpha=0.6)
