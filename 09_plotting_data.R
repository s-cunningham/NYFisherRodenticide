## Plotting covariates


library(tidyverse)
library(sf)

theme_set(theme_bw())

# read data
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

# Break down by age and sex
dat1 <- dat[dat$pt_index==1 & dat$buffsize==60 & dat$radius==100,]
dat1 %>% group_by(Age, Sex) %>% count()

ggplot(dat1, aes(x=rand_x, y=rand_y, color=Age)) + geom_point(size=3)
ggplot(dat1, aes(x=rand_x, y=rand_y, color=Sex)) + geom_point(size=2)

# Plot by number of compounds
cols4 <- c("#9ebcda", "#8c96c6", "#8856a7", "#810f7c")
ggplot(dat1, aes(x=rand_x, y=rand_y, color=factor(n.compounds.MO))) + 
  geom_point(size=4) + scale_color_manual(values=cols4, name="# compounds") + ggtitle("Measured-only")

cols6 <- c("#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#6e016b")
ggplot(dat1, aes(x=rand_x, y=rand_y, color=factor(n.compounds.T))) + 
  geom_point(size=4) + scale_color_manual(values=cols6, name="# compounds") + ggtitle("Including trace")

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


