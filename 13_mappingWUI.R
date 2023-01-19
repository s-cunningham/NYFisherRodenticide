
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(maptiles)



longlat="+proj=longlat +datum=WGS84"

## Read WUI raster
wui <- rast("data/rasters/WUI250m_NZmask_epsg4326.tif")

nz <- st_read("data/spatial", "northern_zone_NAD83") %>% st_union()
plot(st_geometry(nz))
nz <- vect(nz)
nz <- project(nz, crs(wui))

ext(wui) <- ext(nz)
wui <- mask(wui, nz)
wui <- project(wui, longlat)

p <- st_point(c(-75.07020555, 43.582844444)) %>%
      st_sfc(crs=4326) %>%
      st_buffer(750)
tile1 <- get_tiles(p, provider="Stamen.Terrain", zoom=10, cachedir=".")

ggplot() +
  geom_spatraster_rgb(data=tile1, maxcell=Inf) +
  geom_spatraster(data=wui, maxcell=Inf, alpha=0.2) +
  geom_sf(data=p, fill=NA) +
  coord_sf(crs=3857)


ggplot() + geom_spatraster(data=wui) +
              # geom_vline(aes(xintercept=-76), color="red") +
              coord_sf(xlim=c(-76,-74), ylim=c(43.5,44))

ggplot() + geom_spatraster(data=wui) +
  # geom_vline(aes(xintercept=-76), color="red") +
  coord_sf(xlim=c(-76,-74), ylim=c(43.5,44))


ny_utm <- st_read("data/spatial", "NY_outline_UTM18N")
ny_utm <- vect(ny_utm)
ny_ll <- project(ny_utm, longlat)
