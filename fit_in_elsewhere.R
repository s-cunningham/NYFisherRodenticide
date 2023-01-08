library(tidyverse)
# library(sf)

dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
# dat <- dat %>% select(RegionalID:year)

# aea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
# rp <- st_as_sf(dat, coords=c("rand_x", "rand_y"), crs=aea)
# st_write(rp, "data/spatial/random_samples_data20221220.shp", layer_options="SHPT=POINT")


baa <- read_csv("data/analysis-ready/baa_sum_single_raster.csv") %>%
        rename(pt_name=name, BBA=baa)

dat <- left_join(dat, baa, by=c("pt_name", "buffsize")) 
dat <- dat %>% distinct()

mast <- read_csv("data/analysis-ready/ALTEMP26_beech-data.csv")
mast <- mast[mast$year>2014 & mast$year<=2021,1:2]

dat <- left_join(dat, mast, by="year")

names(dat)[45] <- "beechnuts"

mast$year <- mast$year + 1

dat <- left_join(dat, mast, by="year")
names(dat)[46] <- "lag_beechnuts"


write_csv(dat, "data/analysis-ready/combined_AR_covars_new12-2022.csv")

