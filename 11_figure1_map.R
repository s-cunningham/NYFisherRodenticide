library(tidyverse)
library(terra)
library(tidyterra)
library(patchwork)
library(ggspatial)
library(ggnewscale)

## Read layers
# NY outline
ny <- vect("F:/PhDwork/map_data/NY_ne_10m_utm18nad83.shp")

# Harvest area WMUAs
wmuas <- vect("C:/Users/steph/Documents/harvest_area_WMUAs.shp")

# Clipped hydro area dissolve
water <- vect("F:/PhDwork/map_data/hydrogeography/clipped_hyrdo_dissolve.shp")

# Vegetation
wildveg <- rast("F:/PhDwork/map_data/NLCD/nlcd_2019_land_cover_l48_20210604/wildveg.tif")

# Points
pts <- vect("output/random_points_idx1__epsg26918.shp") 

# Hillshade
hs <- rast("F:/PhDwork/map_data/elevation/SRTM/clipped_hillshsrtm3_epsg26918.tif")

## Create map background
# states / provinces
states <- vect("F:/PhDwork/map_data/Northeast_ne_10m_UTM18NAD82.shp") %>%
              aggregate()

# create NY bounding box + buffer
bbox <- c(xmin=105571.4216, ymin=4483089.5, xmax=764142.8229, ymax=4985443.8894)
sq_poly <- sf::st_as_sfc(sf::st_bbox(bbox)) %>% sf::st_set_crs(sf::st_crs(ny))
sq_poly <- vect(sq_poly) 
sq_poly <- buffer(sq_poly, width=10000)

# crop states to bbox
states <- crop(states, sq_poly)

## Reproject vegetation raster
wildveg <- project(wildveg, hs)

## Palette for hillshade
# Hillshading, but we need a palette
pal_greys <- hcl.colors(1000, "Grays")


# Add color table for land cover
coltab(wildveg) <- data.frame(value=c(0,1), 
                              col=c("#dbdbd8", "#518e52"))

# Set up NLCD classes and class values
veg_values <- c(0,1)
veg_class <- c("Non-forest", "Forest")

# Add class names and numbers to the raster
levels(wildveg) <- list(data.frame(ID = veg_values,
                                   NLCD_land = veg_class))

## Plot

ggplot() +
  geom_spatvector(data=states, fill="gray95", color="gray") +
  geom_spatraster(data=wildveg) +
  scale_fill_coltab(data=wildveg, na.value="transparent", alpha=0.4, name="Landcover") +
  geom_spatvector(data=water, color="#2f6e9e", fill="#44a1e3", linewidth=0.01) +
  geom_spatvector(data=wmuas, fill=NA, linewidth=0.75, color="black") +
  geom_spatvector(data=ny, fill=NA, color="black") +
  new_scale_fill() +
  geom_spatvector(data=pts, aes(fill=factor(year)), shape=21, color="black", size=2.5) +
  scale_fill_manual(values=c("#40004b","#9970ab","#e7d4e8"), name="Year") +
  theme_void() +
  theme(legend.position=c(0.92,0.8),
        legend.justification=c(1,1), 
        legend.title=element_text(size=12),
        legend.text=element_text(size=11)) +
  annotation_north_arrow(
    which_north = TRUE,
    pad_x = unit(0.05, "npc"),
    pad_y = unit(0.7, "npc"),
    style = north_arrow_minimal()) +
  annotation_scale(
    height = unit(0.02, "npc"),
    width_hint = 0.4,
    pad_x = unit(0.07, "npc"),
    pad_y = unit(0.07, "npc"),
    text_cex = 1)

ggsave("figs/figure_1_map.svg")
# Saving 9.1 x 7.38 in image