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

b3ar <- read.csv("data/analysis-ready/screening-results_batch3.csv")
b3ar <- b3ar[,c("ID", names(b2ar[2:12]))]

ar <- rbind(b1ar, b2ar, b3ar)
names(ar)[1] <- "RegionalID"

## Read sample details 
# First batch, 2018 samples
b1det <- read.csv("data/analysis-ready/liver_samples_20201021_batch1.csv")
b1det <- b1det[,c(3,9,7,8,6,5,12:18)]
b1det$RegionalID <- paste0("2018-", b1det$RegionalID)

# b1det$RegionalID[!(b1det$RegionalID %in% b1ar$ID)]

# 2nd batch, 2020 samples
b2det <- read.csv("data/analysis-ready/liver_samples_20211207_batch2.csv")
b2det <- b2det[,c(2,5,6,7,9,8,11,12,14,15,17:19)]

# 3rd batch, mainly 2019, some 2020 samples
b3det <- read.csv("data/analysis-ready/liver_samples_20220310_batch3.csv")
b3det <- b3det[,c(1:4,6,5,7,8,17,18)]
fg19 <- read.csv("data/2019_forest_groups.csv")
fg19 <- fg19[,c(2,16:18)]
fg20 <- read.csv("data/2020_forest_groups.csv")
fg20 <- fg20[,c(2,17:19)]
b3det <- b3det %>% left_join(fg19, by="RegionalID") %>% left_join(fg20, by="RegionalID") %>%
  mutate(PctForest = coalesce(PctForest.x, PctForest.y)) %>%
  mutate(ForestGroup = coalesce(ForestGroup.x, ForestGroup.y)) %>%
  mutate(AG = coalesce(AG.x, AG.y)) 
b3det <- b3det[,-c(11:16)]

dets <- bind_rows(b1det, b2det, b3det)
dets <- dets[dets$RegionalID!="2018-6108",]

dets$PctForest[dets$RegionalID=="2018-4033"] <- 54.76190 
dets$ForestGroup[dets$RegionalID=="2018-4033"] <- 2 
dets$Town[dets$RegionalID=="2018-6086"] <- "Pierrepont"
dets$Town[dets$RegionalID=="2018-5215"] <- "Moira"
dets$WMU[dets$RegionalID=="2018-5215"] <- "6C"
dets$WMU[dets$RegionalID=="2018-5384"] <- "5G"
dets$WMU[dets$RegionalID=="2019-5142"] <- "5S"
dets$WMU[dets$RegionalID=="2019-8179"] <- "8T"
dets$Town[dets$RegionalID=="2020-4064"] <- "Pittsfield"
dets$Town[dets$RegionalID=="2019-6207"] <- "Sangerfield"
dets$Town[dets$RegionalID=="2019-6209"] <- "Augusta"
dets$County[dets$County=="St. Lawrsnce" | dets$County=="St.Lawrence" | 
              dets$County=="St. Lawrence" | dets$County=="St. Lawrence "] <- "St Lawrence"
dets$County[dets$County=="oneida"] <- "Oneida"
dets$County[dets$RegionalID=="2018-9211"] <- "Cattaraugus"
dets$County[dets$RegionalID=="2018-6016"] <- "Lewis"
dets$County[dets$RegionalID=="2020-6267"] <- "Oswego"
dets$County[dets$RegionalID=="2020-6140"] <- "Oneida"
dets$WMU[dets$RegionalID=="2018-5318"] <- "6J"
dets$Region[dets$RegionalID=="2018-5318"] <- 6

# Update centroid locations based on town AND wmu
twmu <- st_read("data/spatial", "WMUtown_union_Harv")
twmu <- st_transform(twmu, "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
twmu <- st_drop_geometry(twmu)
names(twmu)[c(1,2,3,5)] <- c("FID", "WMU", "Town", "County")
twmu <- twmu[!(twmu$Town=="Corning" & twmu$MUNI_TYPE=="city"),]
twmu <- twmu[!(twmu$Town=="Little Falls" & twmu$MUNI_TYPE=="city"),]
lymehouns <- twmu[twmu$Town=="Lyme" | twmu$Town=="Hounsfield",]
twmu <- twmu[twmu$Town!="Lyme" & twmu$Town!="Hounsfield",]
lymehouns <- lymehouns[(lymehouns$Town=="Lyme" & lymehouns$x_coord==1570473.02421) |
                         (lymehouns$Town=="Hounsfield" & lymehouns$x_coord==1583912.96229),]
twmu <- bind_rows(twmu, lymehouns)
twmu <- twmu[,-c(1,4,6)]
twmu <- unite(twmu, "key", 1:2, sep="-", remove=FALSE)

dets <- left_join(dets, twmu, by=c("Town", "WMU", "County"))
write.csv(dets, "data/analysis-ready/ar_locations_only.csv")

# Save locations as point file
dets.sf <- st_as_sf(dets, coords=c("x_coord", "y_coord"), crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
st_write(dets.sf, "output/liver_pts.shp")

## Join data to screening results
dat <- left_join(dets, ar, by="RegionalID")
dat <- dat[dat$RegionalID!="2018-9211",] # this one seems like it is wrong but don't know which to fix to make it right

### Check what fishers had fur extracted
# fur <- read.csv("data/harvest_hair.csv")
# 
# sum(fur$Fisher.ID %in% dat$RegionalID)
# furliv <- dat[dat$RegionalID %in% fur$Fisher.ID,]
###

# save as shapefile
# dat.sf <- st_as_sf(dat, coords=c("x_coord", "y_coord"), crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# st_write(dat.sf, "output/liver_pts.shp")

# Reformat date 
dat$HarvestDate <- as.Date(dat$HarvestDate, format="%m/%d/%Y")
dat$year <- as.numeric(format(dat$HarvestDate, "%Y"))
dat$year[is.na(dat$year)] <- c(2018, 2020, 2020, 2020, 2020, 2020) 

# look at years
yr <- dat[,c(1,5,28)]  

# Reformat screening results
dat[17:27] <- lapply(dat[17:27], function(x) replace(x, x=="ND", NA))
dat[17:27] <- lapply(dat[17:27], function(x) replace(x, x=="Traces", 0.000001))

# convert to numeric
dat[,17:27] <- sapply(dat[,17:27], as.numeric)

# Check if there are compounds not detected
for (i in 17:27) {
  print(sum(is.na(dat[,i])))
}

# Convert to long format for ggplotting
datl <- as.data.frame(pivot_longer(dat, cols=17:27, names_to="compound", values_to="ppm"))

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

# write.csv(datl, "output/summarized_AR_results.csv")

dat2 <- datl %>% group_by(RegionalID) %>% summarize(n.compounds=sum(bin.exp))
dat2 <- as.data.frame(dat2)
dat2 <- left_join(dat2, yr, by="RegionalID")
# write.csv(dat2, "output/ncompounds_trace.csv")

dat3 <- datl %>% group_by(RegionalID) %>% summarize(n.compounds=sum(bin.exp.ntr))
dat3 <- as.data.frame(dat3)
dat3 <- left_join(dat3, yr, by="RegionalID")
# write.csv(dat3, "output/ncompounds_notrace.csv")

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

datsf <- st_as_sf(dat, coords=c("x_coord", "y_coord"), crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
datsf <- st_transform(datsf, "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

ggplot() + 
  geom_sf(data=nys, fill="gray40") +
  geom_sf(data=datsf, aes(color=factor(year)), size=2) +
  guides(color=guide_legend(title="Year"))


             