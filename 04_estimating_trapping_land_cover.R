library(tidyverse)
library(terra)
library(sf)

dat <- read_csv("data/JensenHumphries2019/Fisher_RSPF Occurrence & Covariate Data.csv") %>%
          select(OBJECTID, POINT_X, POINT_Y, occur) %>%
          rename(pt_id=OBJECTID, utm_x=POINT_X, utm_y=POINT_Y) %>%
          filter(occur==1) %>%
          mutate(pt_id=1:1078)
          
datsf <- st_as_sf(dat, coords=c("utm_x", "utm_y"), crs=26918)

## Load raster (NLCD)

nlcd <- rast("data/rasters/ny2008nlcd.tif")
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

datsf <- st_transform(datsf, crs=crs(nlcd))


pt_lc <- extract(nlcd, datsf, method="simple")
pt_lc <- pt_lc %>% rename(pt_id=ID)


datsf <- left_join(datsf, pt_lc, by="pt_id")
dat <- bind_cols(st_drop_geometry(datsf), dat[,2:3])

ggplot(dat, aes(x=utm_x, y=utm_y, color=landcov)) + geom_point() +theme_bw()

dat %>% group_by(landcov) %>% count()

st_write(datsf, "data/JensenHumphries2019/harvest_locations_landcov2008.shp")

