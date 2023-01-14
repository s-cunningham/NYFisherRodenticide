## Set up spatial layers for rodenticide analysis
## 2022-06-16

library(rFIA)
library(sf)
library(tidyverse)
library(ranger)
library(terra)

set.seed(123)

#### FIA data ####
# Import from csv files
nydb <- readFIA("data/ny_fia/", inMemory=FALSE)

# Estimate biomass
bm <- biomass(db=nydb, byPlot=TRUE, returnSpatial=TRUE, bySpecies=TRUE, totals=TRUE)

bm <- bm[(bm$YEAR==2018 | bm$YEAR==2019 | bm$YEAR==2020) & (bm$SCIENTIFIC_NAME=="Fagus grandifolia"),]

plot(bm)

# Estimate basal area per acre
baa <- tpa(db=nydb, byPlot=TRUE, returnSpatial=TRUE, bySpecies=TRUE)

baa <- baa[(baa$YEAR==2017 | baa$YEAR==2018 | baa$YEAR==2019 | baa$YEAR==2020) & 
             (baa$SCIENTIFIC_NAME=="Fagus grandifolia") &
             (baa$PROP_FOREST==1),]

plot(st_geometry(baa))

ggplot() +
  geom_sf(data=baa, aes(color=BAA), size=4)

# Convert baa to SpatVect format
baaSV <- vect(baa)

#### Land cover data and covariates ####

## Beech layer (2000-2009) (for comparison)
beech <- rast("data/rasters/NY_Fgrandifolia.tif")

## LANDFIRE layers
slp <- rast("data/rasters/LANDFIRE/LF20_SlpD220.tif")
asp <- rast("data/rasters/LANDFIRE/Asp22_NY.tif")
elev <- rast("data/rasters/LANDFIRE/LF20_Elev220_NY.tif")
evc <- rast("data/rasters/LANDFIRE/EVC_resamp2.tif")
evh <- rast("data/rasters/LANDFIRE/EVH_resamp2.tif")
evt <- rast("data/rasters/LANDFIRE/EVT_resamp2.tif")

# Rename the veg layers so they're not confusing
names(evc) <- "EVC_CLASS"
names(evh) <- "EVH_CLASS"
names(slp) <- "slope"
names(elev) <- "elevation"
names(asp) <- "aspect"

## PRISM layers
tmean <- rast("data/rasters/PRISM_normals/tmean_NY.tif")
tmin <- rast("data/rasters/PRISM_normals/tmin_NY.tif")
tmax <- rast("data/rasters/PRISM_normals/tmax_NY.tif")
dp <- rast("data/rasters/PRISM_normals/dewpoint_NY.tif")
ppt <- rast("data/rasters/PRISM_normals/ppt_NY.tif")
sol <- rast("data/rasters/PRISM_normals/soltotal_NY.tif")
vpdmax <- rast("data/rasters/PRISM_normals/vpdmax_NY.tif")
vpdmin <- rast("data/rasters/PRISM_normals/vpdmin_NY.tif")

# rename climate layers
names(tmin) <- "min_temp"
names(tmax) <- "max_temp"
names(ppt) <- "precip"
names(sol) <- "solar_rad"

# Calculate relative humidity from dewpoint and tmean
relh <- 100 * ((exp(17.625*(dp/(243.04+dp))))/(exp(17.625*(tmean/(243.04+tmean)))))
names(relh) <- "rel_humidity"

# Average vpd
vpd <- (vpdmax + vpdmin)/2
names(vpd) <- "vapor_press_deficit"

## Resample LANDFIRE to beech
slp <- resample(slp, beech, method='bilinear')
elev <- resample(elev, beech, method='bilinear')
asp <- resample(asp, beech, method='bilinear')

## Resample PRISM to beech
tmin <- resample(tmin, beech, method='bilinear')
tmax <- resample(tmax, beech, method='bilinear')
relh <- resample(relh, beech, method='bilinear')
ppt <- resample(ppt, beech, method='bilinear')
sol <- resample(sol, beech, method='bilinear')
vpd <- resample(vpd, beech, method='bilinear')

## Mask landfire veg 
ext(evc) <- ext(slp)
evc <- mask(evc, slp)
ext(evh) <- ext(slp)
evh <- mask(evh, slp)
ext(evt) <- ext(slp)
evt <- mask(evt, slp)

#### Extract values at points ####

# Reproject FIA points
baaSV <- project(baaSV, crs(slp))

# Create raster stacks
stack <- c(slp, elev, asp, tmin, tmax, ppt, relh, sol, vpd, evc, evh, evt) 

# Extract from raster stacks
extr_vals <- extract(stack, baaSV)

#### Run random forest model ####

# Combine covariates with baa
plt_dat <- as.data.frame(baaSV[,c(8)])
plt_dat <- cbind(plt_dat, extr_vals)
plt_dat <- plt_dat[,-2]

# Divide into testing and training
row_idx <- sample(seq_len(nrow(plt_dat)), nrow(plt_dat))
training <- plt_dat[row_idx < nrow(plt_dat) * 0.8, ]
testing <- plt_dat[row_idx >= nrow(plt_dat) * 0.8, ]

# Build random forest model
baa_rf <- ranger(BAA ~ ., data=training, num.trees=1300, mtry=4, min.node.size = 1, 
                 replace=FALSE, sample.fraction=1)

# Predict on trainin data
rf_pred <- predictions(predict(baa_rf, testing))

# Calculate RMSE
sqrt(mean((rf_pred - testing$BAA)^2))

## Tune the random forest
# Mike Mahoney's calc rmse function
calc_rmse <- function(rf_model, data) {
  rf_predictions <- predictions(predict(rf_model, data))
  sqrt(mean((rf_predictions - data$BAA)^2))
}

# Mike Mahoney's k-fold cross validation function
k_fold_cv <- function(data, k, ...) {
  per_fold <- floor(nrow(data) / k)
  fold_order <- sample(seq_len(nrow(data)),
                       size = per_fold * k)
  fold_rows <- split(
    fold_order,
    rep(1:k, each = per_fold)
  )
  vapply(
    fold_rows,
    \(fold_idx) {
      fold_test <- data[fold_idx, ]
      fold_train <- data[-fold_idx, ]
      
      fold_rf <- ranger(BAA ~ ., fold_train, ...)
      calc_rmse(fold_rf, fold_test)
    },
    numeric(1)
  ) |>
    mean()
}

# Create tuning grid
tuning_grid <- expand.grid(
  mtry=floor(ncol(training) * c(0.3,  0.6, 0.9)),
  min.node.size=c(1,3,5),
  replace=c(TRUE, FALSE),
  sample.fraction = c(0.5, 0.63, 0.8),
  rmse=NA
)

# run grid search
for (i in seq_len(nrow(tuning_grid))) {
  tuning_grid$rmse[i] <- k_fold_cv(
    training,
    k = 5,
    mtry = tuning_grid$mtry[i],
    min.node.size = tuning_grid$min.node.size[i],
    replace = tuning_grid$replace[i],
    sample.fraction = tuning_grid$sample.fraction[i]
  )
}

# see how the tuning did
head(tuning_grid[order(tuning_grid$rmse), ])

# Build random forest model with tuned parameters
baa_rf <- ranger(BAA ~ ., data=training, num.trees=1300, mtry=3, min.node.size = 3, 
                 replace=TRUE, sample.fraction=0.5)

# k-fold cross validation
cv_results <- k_fold_cv(plt_dat, k=5, num.trees=1300, mtry=3, min.node.size=3, replace=TRUE, sample.fraction=0.5)

# Predict on training data
rf_pred <- predictions(predict(baa_rf, testing))

# Calculate RMSE of tuned model
sqrt(mean((rf_pred - testing$BAA)^2))

#### Create raster of BAA predictions ####

# Read in raster point centers and convert to SpatVector
pts <- st_read("data/rasters/Fgrandi_pts.shp")
pts <- vect(pts)
pts <- project(pts, crs(slp))
pts_sf <- st_as_sf(pts)

# Extract values from raster stack
vars <- extract(stack, pts)
vars$EVC_CLASS[vars$EVC_CLASS=="NoData"] <- NA
vars <- vars[complete.cases(vars),]

# predict
rf_rast <- predictions(predict(baa_rf, vars))

# add predictions to sf object
pts_sf <- pts_sf[vars$ID,]
pts_sf$BAA_pred <- rf_rast

# Write point dataset
st_write(pts_sf, "output/beech_baa_pred_250.shp")


#### standStructure ####
# ss <- standStruct(db=nydb, grpBy=PLOT, returnSpatial=TRUE, landType='forest', method='LMA', 
#                   totals=TRUE, variance=TRUE, byPlot=TRUE, nCores=2)
# 
# ss <- ss[(ss$YEAR==2017 | ss$YEAR==2018 | ss$YEAR==2019 | ss$YEAR==2020),]
# 
# ggplot() +
#   geom_sf(data=ss, aes(color=STAGE), size=3)
# 




