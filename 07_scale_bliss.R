## Simulating spatial points
library(tidyverse)
library(sf)
library(terra)
library(exactextractr)

M <- 400  # Number of samples

# Read in polygon
ny <- st_read("spatial", "forest_2kmroad_WMUtown_union") %>% 
  unite("key", c(2, 3), sep="_") %>% select(key, geometry)

# Generate random points ("true" locations)
set.seed(10)
pts <- st_sample(ny, M)
pts <- vect(pts)
pts$Id <- 1:M
names(pts)[1] <- "name"
pts <- st_as_sf(pts) %>% 
  st_transform(crs=crs(ny))

# Join Town/WMU key to points
pts <- st_join(pts, left=FALSE, ny["key"])

# Create a buffer around random points
truebuff <- st_buffer(pts, 3090.19)

# Select only polygons that had random points
pt_polys <- ny[pts,]

# How manu points per polygon?
pt_count <- pts %>% group_by(key) %>% count() %>% as_tibble() %>% select(key, n)

# reorder according to polygons
first <- unique(ny$key)
second <- unique(pt_count$key)

reorder_idx <- match(first, second)

pt_count <- pt_count[reorder_idx,] 
pt_count <- pt_count[complete.cases(pt_count),]

## 
pts2 <- st_drop_geometry(pts)

# Select random points per polygon ("jittered" points)
set.seed(11)
samples_per_polygon <- 10*pt_count$n
samples <- st_sample(pt_polys, samples_per_polygon)
samples <- st_as_sf(samples) %>% 
  st_transform(crs=crs(ny))
samples <- st_join(samples, left=FALSE, ny["key"])
samples <- samples %>% mutate(region=gsub("([0-9]+).*$", "\\1",key))
samples$pt_index <- rep(1:10, M)
samples <- samples[,1:4]
stest <- left_join(pts2,samples, by="key") %>% as_tibble() %>%
  select(name, pt_index, key, region) %>% 
  distinct() %>%
  unite("pt_name", 1:2, sep="_", remove=FALSE)
samples <- left_join(samples, stest, by=c("key", "pt_index")) %>% distinct(across(-geometry),.keep_all=TRUE) %>%
  select(key:name, geometry) %>% 
  rename(region=region.x)

ggplot(samples, aes(color=factor(region))) + geom_sf()
# st_write(samples, "spatial/simulation_iter_points.shp", layer_options="SHPT=POINT", delete_dsn=TRUE)

# Create a buffer around random points
randbuff <- st_buffer(samples, 3090.19)

#### Read in rasters and extract values ####

## NLCD
nlcd <- rast("rasters/nybuffnlcd.tif")

pts <- st_transform(pts, crs=crs(nlcd))

# Reclassify 0 (unclassified) and 128 (?) to no data
m <- rbind(c(0, NA), c(128, NA))
nlcd <- classify(nlcd, m)

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
landcov_fracs30 <- exact_extract(nlcd, truebuff, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)

# Organize data
forest_cov <- landcov_fracs30 %>% filter(value==41 | value==42 | value==43)   
forest_cov$value <- factor(forest_cov$value, levels=c(41, 42, 43), labels=c("deciduous", "evergreen", "mixed")) |> as.character()
forest_cov <- forest_cov %>% select(name,value,freq) %>%
  pivot_wider(names_from=value, values_from=freq) %>%
  mutate(tforest = deciduous + evergreen + mixed)

past <- landcov_fracs30 %>% filter(value==81) %>%
  select(name, freq) %>%
  rename(pasture=freq)

lc <- left_join(forest_cov, past, by="name")

## WUI
wui250 <- rast("rasters/WUI250mNY.tif")
wui250 <- project(wui250, nlcd) # Match projection to nlcd

# Set up WUI classes and values
m <- rbind(c(NaN, NA))
wui250 <- classify(wui250, m)
wui_values <- c(0, 1, 2)
wui_class <- c("not WUI", "intermix WUI", "interface WUI")
# Add class names and numbers to the raster
levels(wui250) <- list(data.frame(ID = wui_values,
                                  landcov = wui_class))

wui250_fracs <- exact_extract(wui250, truebuff, function(df) {
  df %>%
    mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(name, value) %>%
    summarize(freq = sum(frac_total))
}, summarize_df = TRUE, include_cols = 'name', progress = FALSE)

wui250_fracs <- pivot_wider(wui250_fracs, names_from=value, values_from=freq) %>%
  rename(intermix=`1`) %>%
  select(name, intermix)

lc <- left_join(lc, wui250_fracs, by="name")

sdf <- st_coordinates(pts)
lc <- bind_cols(lc, sdf)

# Fill in NAs with 0s
lc <- lc %>% mutate(pasture = coalesce(pasture, 0),
                    intermix = coalesce(intermix, 0)) 

##### Replace this with actual values from randomly generated points across NYS, then scale####
# wui <- runif(n=M, -1, 1)       # Scaled WUI of a site
wui <- lc$intermix |> scale()
# pasture <- runif(n=M, -1, 1)   # Scaled % pasture at each site
pasture <- lc$pasture |> scale()
# forest <- runif(n=M, -1, 1)    # Scaled % forest at each site
# forest <- lc$tforest |> scale()
# deciduous <- runif(n=M, -1, 1)    # Scaled % forest at each site
deciduous <- lc$deciduous |> scale()

# Randomly select points on a real plot 
# lat <- runif(n=M, -1, 1)       # Scaled location
lat <- lc$Y |> scale()
# lon <- runif(n=M, -1, 1)       # Scaled location
lon <- lc$X |> scale()

# Set relationships with covariates
mean.count <- 1.757396
beta0 <- log(mean.count) # Same on log scale 
beta1 <- 0.42     # Effect (slope) of WUI
beta2 <- -0.18    # Effect (slope) of % pasture
beta3 <- -0.06     # Effect (slope) of forest cover
beta4 <- 0.3      # Effect (slope) of latitude (but more like UTM coords)
beta5 <- -0.4     # Effect (slope) of longitude (but more like UTM coords)

log.count <- beta0 + beta1*wui + beta2*pasture + beta3*deciduous + beta4*lat + beta5*lon
count <- exp(log.count)   # Inverse link transformation
range(count)

par(mfrow = c(3, 2), mar = c(5,4,2,2), cex.main = 1)
plot(wui, count, ylim = c(0, 6), xlab = "WUI", ylab = "Count", pch=16)
curve(exp(beta0 + beta1*x), range(wui)[1], range(wui)[2], col = "red", xlab = "WUI", ylab = "prob", ylim=c(0,1), lwd = 2, add=TRUE)

plot(pasture, count, ylim = c(0, 6), xlab = "Pasture", ylab = "Count", pch=16)
curve(exp(beta0 + beta2*x), range(pasture)[1], range(pasture)[2], col = "red", ylim = c(0, 1), xlab = "Pasture", ylab = "prob", lwd = 2, add=TRUE)

plot(deciduous, count, ylim = c(0, 6), xlab = "deciduous", ylab = "Count", pch=16)
curve(exp(beta0 + beta3*x), range(deciduous)[1], range(deciduous)[2], col = "red", ylim = c(0, 1), xlab = "deciduous", ylab = "prob", lwd = 2, add=TRUE)

plot(lat, count, ylim = c(0, 6), xlab = "Latitude", ylab = "Count", pch=16)
curve(exp(beta0 + beta4*x), range(lat)[1], range(lat)[2], col = "red", ylim = c(0, 1), xlab = "Latitude", ylab = "prob", lwd = 2, add=TRUE)

plot(lon, count, ylim = c(0, 6), xlab = "Longitude", ylab = "Probability", pch=16)
curve(exp(beta0 + beta5*x), range(lon)[1], range(lon)[2], col = "red", ylim = c(0, 1), xlab = "Longitude", ylab = "prob", lwd = 2, add=TRUE)

hist(count)
par(mfrow=c(1,1))

# Convert to tibble and join response variable
true_locs <- st_drop_geometry(pts)
true_locs <- bind_cols(true_locs, st_coordinates(pts))
true_locs$arcount <- floor(count)
true_locs <- true_locs %>% left_join(lc, by="name") %>%
  as_tibble() %>% select(name,key,arcount,X.x,Y.x:intermix) %>%
  rename(x=X.x, y=Y.x)

# Plot to see where exposure is predicted
ggplot(true_locs, aes(x=x, y=y, color=factor(arcount))) + geom_point() + theme_bw() +theme(legend.position="bottom")


