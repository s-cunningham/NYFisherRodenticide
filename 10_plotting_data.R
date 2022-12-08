## Plotting covariates


library(tidyverse)
library(sf)
library(stars)
library(viridis)
library(ggspatial)
library(patchwork)

theme_set(theme_classic())

#### Covariate plots (boxplots, xy plots) ####
# read data
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

# Categorical number of compounds (collapsing higher numbers)
dat$catNcompMO <- ifelse(dat$n.compounds.MO>=2, "2+", as.character(dat$n.compounds.MO))
dat$catNcompT <- ifelse(dat$n.compounds.T>=3, "3+", as.character(dat$n.compounds.T))

dat$catAge[dat$Age>=3.5] <- "adult"
dat$catAge[dat$Age==2.5] <- "subadult"
dat$catAge[dat$Ag<2.5] <- "juvenile"


#### Exploratory ####
# Break down by age and sex
dat1 <- dat[dat$pt_index==1 & dat$buffsize==60 & dat$radius==100,]
dat1 %>% group_by(Age, Sex) %>% count()

ggplot(dat1, aes(x=rand_x, y=rand_y, color=Age)) + geom_point(size=3)
ggplot(dat1, aes(x=rand_x, y=rand_y, color=Sex)) + geom_point(size=2)

# Plot by number of compounds
cols4 <- c("#9ebcda", "#8c96c6", "#8856a7", "#810f7c")
ggplot(dat1, aes(x=rand_x, y=rand_y, color=factor(n.compounds.MO))) + 
  geom_point(size=4) + scale_color_manual(values=cols4, name="# compounds") + 
  ggtitle("Measured-only") + theme(legend.position="bottom")

cols6 <- c("#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#6e016b")
ggplot(dat1, aes(x=rand_x, y=rand_y, color=factor(n.compounds.T))) + 
  geom_point(size=4) + scale_color_manual(values=cols6, name="# compounds") + 
  ggtitle("Including trace") + theme(legend.position="bottom")

dat1 %>% group_by(n.compounds.MO, year) %>% count()
dat1 %>% group_by(n.compounds.T, year) %>% count()

dat$catNcompT <- ifelse(dat$n.compounds.T>=3, "3+", as.character(dat$n.compounds.T))
dat$catNcompT  <- ordered(dat$catNcompT , levels=c("0", "1", "2", "3+"))
dat$catNcompMO <- ifelse(dat$n.compounds.MO>=2, "2+", as.character(dat$n.compounds.MO))
dat$catNcompMO <- ordered(dat$catNcompMO, levels=c("0", "1", "2+"))

# Histogram

ggplot(dat1, aes(x=n.compounds.T)) + geom_bar() + 
  facet_grid(.~year) + 
  xlab("Number of Compounds Detected per Individual") +
  theme_bw() 

# Age and sex plots
ggplot(dat1, aes(x=catAge, y=n.compounds.T, fill=Sex)) + 
  geom_boxplot()

ggplot(dat1, aes(x=catAge, y=n.compounds.MO, fill=Sex)) + 
  geom_boxplot()

## Box plots for covariates (according to number of compounds)
# Agriculture
pctAG <- dat[, c(1:16,32,33,18:20)]
pctAG <- distinct(pctAG)
pctAG <- pctAG %>% group_by(RegionalID) %>% 
  pivot_longer(19:21, names_to="AG_category", values_to="value") %>% as.data.frame()

ggplot(pctAG, aes(y=catNcompT, x=value, fill=AG_category)) + 
  geom_boxplot() + facet_grid(buffsize~.) + 
  coord_cartesian(ylim=c(0,1)) + ylab("Percent coverage in buffer area") + xlab("Number of Compounds Detected")
  
# Beech basal area per acre
baa <- dat[, c(1:16,32,33,25)]
baa <- distinct(baa)

ggplot(baa, aes(x=catNcompT, y=baa, fill=factor(buffsize))) + 
  geom_boxplot() + xlab("Number of Compounds Detected") +
  ylab("Average beech basal area per acre")

# Wildland-urban interface
wui <- dat[,c(1:17,32,33,26:28)]
wui <- wui %>% group_by(RegionalID) %>% 
  pivot_longer(20:22, names_to="WUI_category", values_to="value") %>% as.data.frame()

ggplot(wui, aes(x=catNcompT, y=value, fill=factor(WUI_category))) + 
  geom_boxplot() + xlab("Number of Compounds Detected") +
  ylab("WUI %") + facet_grid(factor(buffsize)~factor(radius))


## Plot points colored by covariate value
ggplot(dat, aes(x=rand_x, y=rand_y, color=totalWUI)) + geom_point() +
  facet_grid(~buffsize) + theme(legend.position="bottom")

ggplot(dat, aes(x=rand_x, y=rand_y, color=intermix)) + geom_point() +
  facet_grid(~buffsize) + theme(legend.position="bottom")

ggplot(dat, aes(x=rand_x, y=rand_y, color=interface)) + geom_point() +
  facet_grid(~buffsize) + theme(legend.position="bottom")

ggplot(dat, aes(x=rand_x, y=rand_y, color=pasture)) + geom_point() +
  facet_grid(~buffsize) + theme(legend.position="bottom")

ggplot(dat, aes(x=rand_x, y=rand_y, color=crops)) + geom_point() +
  facet_grid(~buffsize) + theme(legend.position="bottom")

ggplot(dat, aes(x=rand_x, y=rand_y, color=(nbuildings/buffsize))) + geom_point() +
  facet_grid(~buffsize) + theme(legend.position="bottom")

dat60 <- dat[dat$buffsize==60,]
ggplot(dat60, aes(x=rand_x, y=rand_y, color=intermix)) + geom_point(size=2) +
  theme(legend.position="bottom")

ggplot(dat60, aes(x=rand_x, y=rand_y, color=interface)) + geom_point(size=2) +
  theme(legend.position="bottom")

ggplot(dat60, aes(x=rand_x, y=rand_y, color=totalWUI)) + geom_point(size=2) +
  theme(legend.position="bottom")

ggplot(dat60, aes(x=rand_x, y=rand_y, color=pasture)) + geom_point(size=2) +
  theme(legend.position="bottom")

ggplot(dat60, aes(x=rand_x, y=rand_y, color=crops)) + geom_point(size=2) +
  theme(legend.position="bottom")

ggplot(dat60, aes(x=rand_x, y=rand_y, color=totalag)) + geom_point(size=2) +
  theme(legend.position="bottom")

ggplot(dat60, aes(x=rand_x, y=rand_y, color=baa)) + geom_point(size=2) +
  theme(legend.position="bottom")


## Subset to only the most common compounds ##

# Read data
dets <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dets <- as.data.frame(dets)
dets <- dets[,c(1:5, 7:9, 13, 16:20, 25:29 )]

dat <- read_csv("output/summarized_AR_results.csv")
dat <- as.data.frame(dat[,2:23])

# Subset to diphacinone, bromadiolone, brodifacoum
dat <- dat[dat$compound=="Diphacinone" | dat$compound=="Brodifacoum" |
             dat$compound=="Bromadiolone", c(1,18,20) ]

dat <- left_join(dets, dat, by="RegionalID")

dat <- pivot_wider(dat, values_from="exposure", names_from="compound") %>% as.data.frame()


ggplot(dat, aes(x=rand_x, y=rand_y, color=factor(Diphacinone))) + geom_point() +
  theme(legend.position="bottom")

ggplot(dat, aes(x=rand_x, y=rand_y, color=factor(Bromadiolone))) + geom_point() +
  theme(legend.position="bottom")

ggplot(dat, aes(x=rand_x, y=rand_y, color=factor(Brodifacoum))) + geom_point() +
  theme(legend.position="bottom")

# Boxplots

dat$diphac_binary <- ifelse(dat$Diphacinone=="ND", "not detected", "detected")
dat$bromod_binary <- ifelse(dat$Bromadiolone=="ND", "not detected", "detected")
dat$brodif_binary <- ifelse(dat$Brodifacoum=="ND", "not detected", "detected")

ggplot(dat, aes(x=diphac_binary, y=totalag)) + geom_boxplot() + 
  theme(legend.position="bottom") + facet_grid(buffsize~radius)

ggplot(dat, aes(x=bromod_binary, y=totalag)) + geom_boxplot() + 
  theme(legend.position="bottom") + facet_grid(buffsize~radius)

ggplot(dat, aes(x=brodif_binary, y=totalag)) + geom_boxplot() + 
  theme(legend.position="bottom") + facet_grid(buffsize~radius)

# crops
ggplot(dat, aes(x=diphac_binary, y=crops)) + geom_boxplot() + 
  theme(legend.position="bottom") + facet_grid(buffsize~radius)

ggplot(dat, aes(x=bromod_binary, y=crops)) + geom_boxplot() + 
  theme(legend.position="bottom") + facet_grid(buffsize~radius)

ggplot(dat, aes(x=brodif_binary, y=crops)) + geom_boxplot() + 
  theme(legend.position="bottom") + facet_grid(buffsize~radius)



#### Figure 1 for manuscript ####
## Load raster layer
beech <- read_stars("data/rasters/baa250_masked.tif")

# New york state
aea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
utm <- "+proj=utm +zone=18 +datum=WGS84"
nys <- st_read("data/spatial", "NYS_outline_albers")
nys <- st_transform(nys, utm)

# Fisher harvest area
fha <- st_read("data/spatial", "FisherHarvestArea_dissolved")
fha <- st_transform(fha, utm)

# Water
water <- st_read("data/spatial", "water_over_1sq_km")
water <- st_transform(water, utm)

# make data sf object
dat <- dat[dat$buffsize==15,]
datsf <- st_as_sf(dat, coords=c(4,5), crs=aea)
datsf <- st_transform(datsf, utm)
dat2 <- st_coordinates(datsf)
dat <- bind_cols(dat, dat2)

# Read in polygon layer that has union of towns and WMUs
twmu <- st_read("data/spatial", "ARtowns2")
twmu <- st_transform(twmu, utm)

# Plots
ggplot() + geom_sf(data=nys, fill="gray80", color="gray20") + 
  geom_sf(data=fha, fill="white", color="gray20") + 
  geom_point(data=dat, aes(x=X, y=Y, color=factor(year)), shape=16, size=0.8) +
  geom_sf(data=twmu, fill=NA) + 
  coord_sf(xlim=c(100000, 649000), ylim=c(4595000, 4980000)) + 
  scale_color_manual(values=c("#1b9e77", "#d95f02", "#7570b3"), name="Year") +
  guides(colour=guide_legend(override.aes=list(size=3))) +
  annotation_scale(location="bl", text_cex=1, style="bar") +
  theme(legend.position=c(0,1), legend.justification=c(0,1),
        legend.background=element_rect(fill=NA),
        axis.title=element_blank(),
        legend.title=element_text(size=12, face="bold"), 
        legend.text=element_text(size=12),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.border=element_rect(color="black", fill=NA, size=0.5))




