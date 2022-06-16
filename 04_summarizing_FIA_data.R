library(rFIA)
library(sf)

#### Set up spatial area for extracting FIA data ####
# Read in NY shapefile
ny <- st_read("data/spatial", "NYS_outline_albers")
ny <- st_as_sf(ny)

# create grid
spacing = 15000 # 15 km height & width
ny_grid <- st_make_grid(ny, cellsize = c(spacing, spacing), square=FALSE,
                          what = 'polygons') 

# Plot
plot(ny$geometry)
plot(ny_grid, col=NA, add=TRUE)

# Get rid of extra hexagons
ny_grid <- ny_grid[ny]

plot(ny$geometry)
plot(ny_grid, add=TRUE)


#### FIA data ####
nydb <- readFIA("data/ny_fia/")

biomass(db=nydb, byPlot=TRUE)
