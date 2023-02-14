
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(maptiles)



longlat="+proj=longlat +datum=WGS84"

## Read WUI raster
wui <- rast("data/rasters/WUI250m_NZmask_epsg4326.tif")
ggplot() + geom_spatraster(data=wui, maxcell=1000000) + 
  scale_fill_continuous(na.value=NA) +
  geom_vline(aes(xintercept=c(-75.5, -75))) +
  geom_hline(aes(yintercept=c(43.3, 43.5)))


ggplot() + geom_spatraster(data=wui, maxcell=1000000) + 
  scale_fill_continuous(na.value=NA) +
  coord_sf(xlim=c(-75.5, -75), ylim=c(43.3, 43.5))

p <- st_point(c(mean(c(-75.5, -75)), mean(c(43.3, 43.5)))) %>%
      st_sfc(crs=4326) %>%
      st_buffer(0.2)
tile1 <- get_tiles(p, provider="Stamen.Terrain", zoom=10, cachedir=".")

ggplot() +
  # geom_spatraster_rgb(data=tile1, maxcell=Inf) +
  geom_spatraster(data=wui, maxcell=Inf, alpha=0.2) +
  geom_sf(data=p, fill=NA) +
  coord_sf(crs=3857,xlim=c(-75.5, -75), ylim=c(43.3, 43.5))


ggplot() + geom_spatraster(data=wui) +
              # geom_vline(aes(xintercept=-76), color="red") +
              coord_sf(xlim=c(-76,-74), ylim=c(43.5,44))

ggplot() + geom_spatraster(data=wui) +
  # geom_vline(aes(xintercept=-76), color="red") +
  coord_sf(xlim=c(-76,-74), ylim=c(43.5,44))


ny_utm <- st_read("data/spatial", "NY_outline_UTM18N")
ny_utm <- vect(ny_utm)
ny_ll <- project(ny_utm, longlat)
