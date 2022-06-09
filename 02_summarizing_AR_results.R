## Summarizing rodenticide
## 2022-04-26

library(tidyverse)
library(sf)

## Read AR screening results
b1ar <- read.csv("data/analysis-ready/screening-results_batch1.csv")
b1ar[b1ar=="<MDL"] <- "Traces"
b1ar$ID <- paste0("2018-", b1ar$ID)

b2ar <- read.csv("data/analysis-ready/screening-results_batch2.csv")
b2ar <- b2ar[,-1]

ar <- rbind(b1ar, b2ar)
names(ar)[1] <- "RegionalID"

## Read sample details 
b1det <- read.csv("data/analysis-ready/liver_samples_20201021_batch1.csv")
b1det <- b1det[,c(3,9,7,8,6,5,12:18)]
b1det$RegionalID <- paste0("2018-", b1det$RegionalID)

b2det <- read.csv("data/analysis-ready/liver_samples_20211207_batch2.csv")
b2det <- b2det[,c(2,5,6,7,9,8,11,12,14,15,17:19)]

b3det <- read.csv("data/analysis-ready/liver_samples_20220310_batch3.csv")
b3det <- b3det[,c(1:4,6,5,7,8,17,18)]

dets <- rbind(b1det, b2det)

## Join data to screening results
dat <- left_join(dets, ar, by="RegionalID")
dat <- dat[dat$RegionalID!="2018-9211",]

dat$HarvestDate <- as.Date(dat$HarvestDate, format="%m/%d/%Y")
dat$year <- as.numeric(format(dat$HarvestDate, "%Y"))
dat$year[is.na(dat$year)] <- c(2018, 2020, 2020, 2020, 2020)

yr <- dat[,c(1,5,25)]

# Reformat
dat[14:24] <- lapply(dat[14:24], function(x) replace(x, x=="ND", NA))
dat[14:24] <- lapply(dat[14:24], function(x) replace(x, x=="Traces", 0.000001))

# convert to numeric
dat[,14:24] <- sapply(dat[,14:24], as.numeric)

# Check if there are compounds not detected
for (i in 14:24) {
  print(sum(is.na(dat[,i])))
}

# Convert to long format for ggplotting
datl <- as.data.frame(pivot_longer(dat, cols=14:24, names_to="compound", values_to="ppm"))

# Remove compounds with no detections
datl <- datl[datl$compound!="Pindone" & datl$compound!="Coumafuryl" & 
             datl$compound!="Coumachlor",]

# create column for describing exposure by factor level
datl$exposure <- "ND"
datl$exposure[datl$ppm==0.000001] <- "trace"
datl$exposure[datl$ppm>=0.001] <- "measured"

# How many compounds?
datl$bin.exp <- ifelse(datl$exposure=="ND", 0, 1)
datl$bin.exp.ntr <- ifelse(datl$exposure=="measured", 1, 0)

dat2 <- datl %>% group_by(RegionalID) %>% summarize(n.compounds=sum(bin.exp))
dat2 <- as.data.frame(dat2)
dat2 <- left_join(dat2, yr, by="RegionalID")

dat3 <- datl %>% group_by(RegionalID) %>% summarize(n.compounds=sum(bin.exp.ntr))
dat3 <- as.data.frame(dat3)
dat3 <- left_join(dat3, yr, by="RegionalID")

dat2s <- dat2 %>% group_by(n.compounds, year) %>% count()
dat2s$Trace <- "yes"

dat3s <- dat3 %>% group_by(n.compounds, year) %>% count()
dat3s$Trace <- "no"

datbar <- rbind(dat2s, dat3s)

ggplot(datbar) +
  geom_bar(aes(x=n.compounds, y=n, fill=Trace), stat="identity", position=position_dodge()) + 
  facet_grid(.~year) + theme_bw() + ylab("Count") + xlab("Number of Compounds") +
  theme(strip.text=element_text(size=14),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14, face="bold"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14, face="bold")) 

ggplot(dat2) +
  geom_tile(aes(x=RegionalID, y=n.compounds, fill=n.compounds), color="black") +
  facet_grid(year~., space="free_y", scales="free")

ggplot(dat3) +
  geom_tile(aes(x=RegionalID, y=n.compounds, fill=n.compounds), color="black") +
  facet_grid(year~., space="free_y", scales="free")

# plot by concentration
ggplot(datl) +
  geom_tile(aes(x=RegionalID, y=compound, fill=ppm), color="gray80") +
  scale_fill_gradient(low="#ffffcc", high="#081d58", space="Lab", na.value="gray80", limits=c(0,1)) +
  theme_bw() +
  theme(axis.text.x=element_blank())


# Read in state polygon
nys <- st_read("data/spatial", "NYS_outline_albers")
nys <- st_transform(nys, "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
# nys <- st_transform(nys, "+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs")

# Read in town shapefile
twn <- st_read("data/spatial", "Cities_Towns")
twn <- st_transform(twn, "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
cents <- st_centroid(twn)
centid <- twn$NAME
centlab <- data.frame(id=centid, longitude=cents[,1], latitude=cents[,2])
twn2 <- fortify(twn, region="NAME")

# Read in WMU shapefile
wmu <- st_read("data/spatial", "wmus")
wmu <- st_transform(wmu, CRS="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
wmu.c <- st_centroid(wmu)
wmu.cid <- wmu$Name
wmu.cl <- data.frame(id=wmu.cid, longitude=wmu.c[,1], latitude=wmu.c[,2])
wmu2 <- fortify(wmu, region="Name")

# Read in fisher harvest area
fha <- st_read("data/spatial", "FisherHarvestArea")
fha <- st_transform(fha, CRS="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
fha2 <- fortify(fha, region="Name")

names(wmu.cl)[1] <- "WMU"

# Add town coordinates
dat <- left_join()

ggplot() + 
  geom_sf(data=nys, fill="gray40") +
  geom_point(data=dat, aes(x=longitude, y=latitude)) +facet_grid(.~year)
