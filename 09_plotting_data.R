## Plotting covariates


library(tidyverse)
library(sf)

theme_set(theme_bw())

# read data
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

# Categorical number of compounds (collapsing higher numbers)
dat$catNcompMO <- ifelse(dat$n.compounds.MO>=2, "2+", as.character(dat$n.compounds.MO))
dat$catNcompT <- ifelse(dat$n.compounds.T>=3, "3+", as.character(dat$n.compounds.T))

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

## Box plots for covariates (according to number of compounds)
# Agriculture
pctAG <- dat[, c(1:16,32,33,18:20)]
pctAG <- distinct(pctAG)
pctAG <- pctAG %>% group_by(RegionalID) %>% 
  pivot_longer(19:21, names_to="AG_category", values_to="value") %>% as.data.frame()

ggplot(pctAG, aes(x=catNcompT, y=value, fill=AG_category)) + 
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

ggplot(dat, aes(x=rand_x, y=rand_y, color=totalag)) + geom_point() +
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

## Read data frame with town/wmu location
loc <- read.csv("data/analysis-ready/ar_locations_only.csv")
loc <- loc[,-1]
loc <- loc[loc$RegionalID!="2018-9211",] # Seems wrong

## Read in polygon layer that has union of towns and WMUs
twmu <- st_read("data/spatial", "WMUtown_union_Harv")
aea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
twmu <- st_transform(twmu, aea)
twmu <- unite(twmu, "key", 2:3, sep="-", remove=FALSE)
# Remove town/WMU combos that don't have a fisher liver from them
twmu <- twmu[twmu$key %in% loc$key,]

# Count how many individual samples from each town-WMU combo
nloc <- loc %>% group_by(key) %>% count()
twmu <- left_join(twmu, nloc, by="key")

ggplot(twmu, aes(fill=n)) + geom_sf() 


