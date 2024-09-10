library(terra)
library(tidyterra)
library(ggspatial)
library(ggplot2)

## Load raster layers
beech <- rast("data/rasters/baa250_masked.tif")

# Do math to get total beech mass (sq.ft) instead of basal area per acre
# Cells are 65200 m^2, there are 15.444 acres in 65200 m^2
beech <- beech * 15.444 / 10.764

writeRaster(beech, "data/rasters/baa250_m2.tif")




beech <- rast("data/rasters/baa250m2_utm18wgs84.tif")

ny <- vect("data/spatial/NYS_outline_albers.shp")
ny <- project(ny, beech)



map <- ggplot() +
  geom_spatraster(data=beech, na.rm=FALSE) +
  geom_spatvector(data=ny, fill=NA, color="black", linewidth=1) +
  scale_fill_gradient(low = '#ccece6',high = '#005824', na.value = "transparent",
                      name = expression("Basal area ("~m^2~")" )) +
  theme_void() +
  annotation_north_arrow(
    which_north = TRUE,
    pad_x = unit(0.75, "npc"),
    pad_y = unit(0.75, "npc"),
    style = north_arrow_fancy_orienteering()) +
  annotation_scale(
    height = unit(0.015, "npc"),
    width_hint = 0.45,
    pad_x = unit(0.07, "npc"),
    pad_y = unit(0.2, "npc"),
    text_cex = .8) +
  theme(legend.position = c(.2,0.8),
        legend.justification=c(0,1))

ggsave("figs/beech_basal_area_885x735.svg", plot=map)
