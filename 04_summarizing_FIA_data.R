## Set up spatial layers for rodenticide analysis
## 2022-06-16

library(rFIA)
library(sf)
library(tidyverse)
library(randomForest)
library(terra)

#### FIA data ####

# Import from csv files
nydb <- readFIA("data/ny_fia/", inMemory=FALSE)

# Estimate biomass
bm <- biomass(db=nydb, byPlot=TRUE, returnSpatial=TRUE, bySpecies=TRUE, totals=TRUE)

bm <- bm[(bm$YEAR==2018 | bm$YEAR==2019 | bm$YEAR==2020) & (bm$SCIENTIFIC_NAME=="Fagus grandifolia"),]


plot(bm)

# Estimate basal area per acre
baa <- tpa(db=nydb, byPlot=TRUE, returnSpatial=TRUE, bySpecies=TRUE)

baa <- baa[(baa$YEAR==2017 | baa$YEAR==2018 | baa$YEAR==2019 | baa$YEAR==2020) & (baa$SCIENTIFIC_NAME=="Fagus grandifolia"),]

plot(st_geometry(baa))

ggplot() +
  geom_sf(data=baa, aes(color=BAA), size=4)

#### Land cover data and covariates ####

## LANDFIRE layers
slp <- rast("data/rasters/LANDFIRE/LF20_SlpD220.tif")
asp <- rast("data/rasters/LANDFIRE/Asp22_NY.tif")
elev <- rast("data/rasters/LANDFIRE/LF20_Elev220_NY.tif")
evc <- rast("data/rasters/LANDFIRE/LC16_EVC200_NY.tif")
evh <- rast("data/rasters/LANDFIRE/LF16_EVH200_NY.tif")
evt <- rast("data/rasters/LANDFIRE/LF20_EVT200_NY.tif")

## PRISM layers
tmean <- rast("data/rasters/PRISM_normals/tmean_NY.tif")
tmin <- rast("data/rasters/PRISM_normals/tmin_NY.tif")
tmax <- rast("data/rasters/PRISM_normals/tmax_NY.tif")
dp <- rast("data/rasters/PRISM_normals/dewpoint_NY.tif")
ppt <- rast("data/rasters/PRISM_normals/ppt_NY.tif")
sol <- rast("data/rasters/PRISM_normals/soltotal_NY.tif")
vpdmax <- rast("data/rasters/PRISM_normals/vpdmax_NY.tif")
vpdmin <- rast("data/rasters/PRISM_normals/vpdmin_NY.tif")

# Calculate relative humidity from dewpoint and tmean
relh <- 100 * ((exp(17.625*(dp/(243.04+dp))))/(exp(17.625*(tmean/(243.04+tmean)))))

# Average vpd
vpd <- (vpdmax + vpdmin)/2

## Resample LANDFIRE to PRISM
slp <- resample(slp, tmin, method='bilinear', filename="data/rasters/LANDFIRE/slp_resamp.tif")
elev <- resample(elev, tmin, method='bilinear', filename="data/rasters/LANDFIRE/elev_resamp.tif")
asp <- resample(asp, tmin, method='bilinear', filename="data/rasters/LANDFIRE/asp_resamp.tif")






#### Extract values at points ####




#### Run random forest model ####


















