library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(maptiles)
library(patchwork)

# Load raster
# wui <- rast("data/rasters/NZ_WUI250m.tif")
# 
# # Set up WUIclasses and class values
# wui_values <- c(0, 1, 2)
# wui_class <- c("Not WUI", "Intermix", "Interface")
# 
# # Add class names and numbers to the raster
# levels(wui) <- list(data.frame(ID = wui_values,
#                                wuiclass = wui_class))

interface <- rast("data/rasters/interface_box.tif")
intermix <- rast("data/rasters/intermix_box.tif")
notwui <- rast("data/rasters/notwui_box.tif")

mbox <- st_read("data/spatial", "NZ_map_box")
roads <- st_read("data/spatial", "NZ_box_roads")
build <- st_read("data/spatial", "NZ_box_buildings")



# get tile from a point in sf format
p1 <- st_point(c(-75.43141418299894, 43.77274929367046)) %>%
      st_sfc(crs=4326) %>%
      st_buffer(750)

p2 <- st_point(c(-75.30716109326626, 43.773150761278686)) %>%
  st_sfc(crs=4326) %>%
  st_buffer(750)

tile1 <- get_tiles(p1, provider="Stamen.TerrainBackground", zoom=11, cachedir=".")
tile2 <-  get_tiles(p2, provider="Stamen.TerrainBackground", zoom=11, cachedir=".")

mix <- ggplot() +
  geom_spatraster_rgb(data=tile1, maxcell=Inf) +
  geom_spatraster_rgb(data=tile2, maxcell=Inf) +
  geom_spatraster(data=intermix, maxcell=1000000, alpha=0.8) +
  scale_fill_gradient(low="gray60", high="gray60", na.value=NA) +
  geom_sf(data=roads) +
  geom_sf(data=build, shape=20, size=0.5) +
  coord_sf(xlim=c(-75.49316986833249, -75.24503389115974), ylim=c(43.72828172346491, 43.818319582625065)) +
  theme_void() +
  theme(legend.position="none")

face <- ggplot() +
  geom_spatraster_rgb(data=tile1, maxcell=Inf) +
  geom_spatraster_rgb(data=tile2, maxcell=Inf) +
  geom_spatraster(data=interface, maxcell=1000000, alpha=0.8) +
  scale_fill_gradient(low="gray60", high="gray60", na.value=NA) +
  geom_sf(data=roads) +
  geom_sf(data=build, shape=20, size=0.5) +
  coord_sf(xlim=c(-75.49316986833249, -75.24503389115974), ylim=c(43.72828172346491, 43.818319582625065)) +
  theme_void() +
  theme(legend.position="none")

not <- ggplot() +
  geom_spatraster_rgb(data=tile1, maxcell=Inf) +
  geom_spatraster_rgb(data=tile2, maxcell=Inf) +
  geom_spatraster(data=notwui, maxcell=1000000, alpha=0.8) +
  scale_fill_gradient(low="gray60", high="gray60", na.value=NA) +
  geom_sf(data=roads) +
  geom_sf(data=build, shape=20, size=0.5) +
  coord_sf(xlim=c(-75.49316986833249, -75.24503389115974), ylim=c(43.72828172346491, 43.818319582625065)) +
  theme_void() +
  theme(legend.position="none")

face / mix / not + plot_annotation(tag_levels="a", tag_prefix="(", tag_suffix=")")

  # geom_spatraster(data=wui, alpha=0.5, maxcell=2000000) +
  # scale_fill_discrete(na.value=NA) +
  # geom_sf(data=p1) +
  # geom_sf(data=p2) +
  # coord_sf(crs=4326, xlim=c(-75.55, -75.25), ylim=c(43.6, 43.7))

# xlim=c(-75.1, -74.9), ylim=c(43.55, 43.65), 
# geom_vline(aes(xintercept=c(-75.55, -75.25)))
# geom_hline(aes(yintercept=c(43.6, 43.7))) 