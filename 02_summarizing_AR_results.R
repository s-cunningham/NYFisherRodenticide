library(tidyverse)
library(sf)

## Read AR screening results
b1ar <- read_csv("data/analysis-ready/screening-results_batch1.csv")
b1ar[b1ar=="<MDL"] <- "Traces"
b1ar$ID <- paste0("2018-", b1ar$ID)
b1ar$ID[b1ar$ID=="2018-6000"] <- "2018-6001" 

b2ar <- read_csv("data/analysis-ready/screening-results_batch2.csv")
b2ar <- b2ar[,-1]

b3ar <- read_csv("data/analysis-ready/screening-results_batch3.csv")
b3ar <- b3ar[,c("ID", names(b2ar[2:12]))]

ar <- rbind(b1ar, b2ar, b3ar)
names(ar)[1] <- "RegionalID"

## Read in liver data
dets <- read_csv("data/analysis-ready/2016-2020_ages_data.csv")

## Combine location and age data with rodenticide results
dat <- left_join(ar, dets, by="RegionalID") %>% 
  select(RegionalID, HarvestYear:Town, Warfarin:Coumachlor ) %>% distinct()

# Reformat screening results
dat[10:20] <- lapply(dat[10:20], function(x) replace(x, x=="ND", NA))
dat[10:20] <- lapply(dat[10:20], function(x) replace(x, x=="Traces", 0.000001))

dat$HarvestYear <- ifelse(is.na(dat$HarvestYear), str_split(dat$RegionalID, "-")[[1]][1], dat$HarvestYear)

dat$Town[dat$RegionalID=="2018-6086"] <- "Pierrepont"
dat$Town[dat$RegionalID=="2018-5215"] <- "Moira"
dat$WMU[dat$RegionalID=="2018-5215"] <- "6C"
dat$WMU[dat$RegionalID=="2018-5384"] <- "5G"
dat$WMU[dat$RegionalID=="2019-5142"] <- "5S"
dat$WMU[dat$RegionalID=="2019-8179"] <- "8T"
dat$Town[dat$RegionalID=="2020-4064"] <- "Pittsfield"
dat$Town[dat$RegionalID=="2019-6207"] <- "Sangerfield"
dat$Town[dat$RegionalID=="2019-6209"] <- "Augusta"
dat$County[dat$County=="oneida"] <- "Oneida"
dat$County[dat$RegionalID=="2018-9211"] <- "Cattaraugus"
dat$County[dat$RegionalID=="2018-6016"] <- "Lewis"
dat$County[dat$RegionalID=="2020-6267"] <- "Oswego"
dat$County[dat$RegionalID=="2020-6140"] <- "Oneida"
dat$WMU[dat$RegionalID=="2018-5318"] <- "6J"
dat$Region[dat$RegionalID=="2018-5318"] <- 6
dat$WMU[dat$Town=="Argyle" & dat$WMU=="5A"] <- "5S"
dat$WMU[dat$Town=="Freedom" & dat$WMU=="9W"] <- "9N"
dat$WMU[dat$Town=="Western" & dat$WMU=="5H"] <- "6K"
dat$HarvestYear[dat$RegionalID=="2020-4018"] <- "2020"
dat$Sex[dat$RegionalID=="2019-7709" | dat$RegionalID=="2020-70001"] <- "F"

dat <- dat %>% mutate_at(c(10:20), as.numeric) %>% 
  unite("key", c(WMU,Town), sep="-", remove=FALSE) %>%
  arrange(key) %>%
  rename(year=HarvestYear)

unique(dat$key)

# Check if there are compounds not detected
for (i in 10:20) {
  print(sum(is.na(dat[,i])))
}

## write wide data
# write_csv(dat, "output/AR_results_wide.csv")

ordin <- dat %>% mutate(ar_cat=case_when(is.na(Warfarin) & is.na(Coumafuryl) & is.na(Diphacinone) & is.na(Pindone) &
                                         is.na(Brodifacoum) & is.na(Difenacoum) & is.na(Bromadiolone) & is.na(Coumachlor) &
                                         is.na(Chlorophacinone) & is.na(Difethialone) & is.na(Dicoumarol) ~ "none",
                                         !is.na(Warfarin) | !is.na(Coumafuryl) | !is.na(Diphacinone) | !is.na(Coumachlor) |
                                           !is.na(Pindone) | !is.na(Chlorophacinone) | !is.na(Chlorophacinone) |
                                         !is.na(Dicoumarol) ~ "FGAR",
                                         !is.na(Brodifacoum) | !is.na(Difenacoum) | 
                                           !is.na(Bromadiolone) | !is.na(Difethialone) ~ "SGAR",
                                         (!is.na(Pindone) |  !is.na(Chlorophacinone)| !is.na(Warfarin) | 
                                           !is.na(Coumafuryl) | !is.na(Diphacinone) | !is.na(Coumachlor) |
                                            !is.na(Dicoumarol)) &
                                         (!is.na(Brodifacoum) | !is.na(Difenacoum) | 
                                           !is.na(Bromadiolone) | !is.na(Difethialone)) ~ "both")) %>%
        select(RegionalID, ar_cat)


ggplot(ordin) +
  geom_bar(aes(x=ar_cat))

# Separate dicoumarol
dicoum <- dat %>% select(RegionalID:Town, Dicoumarol)

## Convert to long format for ggplotting and later analysis
datl <- dat %>% pivot_longer(cols=c(11:21), names_to="compound", values_to="ppm")

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

# write to file
# write_csv(datl, "output/summarized_AR_results.csv")

brod <- datl %>% filter(compound=="Brodifacoum") %>% select(RegionalID:Town, bin.exp, bin.exp.ntr)
brom <- datl %>% filter(compound=="Bromadiolone") %>% select(RegionalID:Town, bin.exp, bin.exp.ntr)
diph <- datl %>% filter(compound=="Diphacinone") %>% select(RegionalID:Town, bin.exp, bin.exp.ntr)
dico <- datl %>% filter(compound=="Dicoumarol") %>% select(RegionalID:Town, bin.exp, bin.exp.ntr)

# write_csv(brod,"output/binary_brodifacoum.csv")
# write_csv(brom,"output/binary_bromadiolone.csv")
# write_csv(diph,"output/binary_diphacinone.csv")
# write_csv(dico,"output/binary_dicoumarol.csv")

## Summarize by number of compounds
# look at years
yr <- dat[,c(1:7)]

# with trace
dat2 <- datl %>% group_by(RegionalID) %>%
          filter(compound!="Dicoumarol") %>%
          summarize(n.compounds=sum(bin.exp))
dat2 <- left_join(dat2, yr, by="RegionalID") %>%
          select(RegionalID, year:key, n.compounds)
# write.csv(dat2, "output/ncompounds_trace.csv")

# without trace
dat3 <- datl %>% group_by(RegionalID) %>% 
          filter(compound!="Dicoumarol") %>%
          summarize(n.compounds=sum(bin.exp.ntr))
dat3 <- left_join(dat3, yr, by="RegionalID")%>%
          select(RegionalID, year:key, n.compounds)
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

ggplot(dat2) +
  geom_bar(aes(x=n.compounds), fill="darkgreen") +
  ylab("Number of fishers") + xlab("Number of compounds detected") +
  theme_classic(base_size=12) +
  theme(panel.border=element_rect(linewidth=0.5, fill=NA, color="black"))

dat2 %>% group_by(n.compounds) %>% count()

# plot by concentration
ggplot(datl) +
  geom_tile(aes(x=RegionalID, y=compound, fill=ppm), color="gray80") +
  scale_fill_gradient(low="#ffffcc", high="#081d58", space="Lab", na.value="gray80", limits=c(0,1)) +
  theme_bw() +
  theme(axis.text.x=element_blank())

dat2 <- dat2 %>% mutate(binexp=if_else(n.compounds==0, 0, 1))

dat2 %>% summarize(sum(binexp))
dat2 %>% group_by(year) %>% count()

dat2 %>% group_by(Sex) %>% summarize(sum(binexp))
dat2 %>% group_by(Sex) %>% count()

#### Subset polygon ####
twmu <- st_read("data/spatial/WMUtown_union_Harv.shp")
twmu <- unite(twmu, "key", c("Name", "NAME_1"), sep="-")

twmu <- twmu %>% filter(key %in% unique(dat$key)) %>% select(key, x_coord, y_coord, geometry)

# twmu <- st_transform(twmu, 32618)
# xy <- st_centroid(twmu)
# xy <- st_coordinates(xy)
# 
# twmu <- st_drop_geometry(twmu)
# twmu <- bind_cols(twmu, xy)
# 
# plot(st_geometry(twmu))

# Remove islands
twmu <- twmu[-c(104:109,130,145),]

twmu <- st_zm(twmu)
# twmu <- st_transform(twmu, crs=32632)

# st_write(twmu, "data/spatial/AR_towns_WMUs.shp", driver = "ESRI Shapefile", append=FALSE)


towns <- st_read("data/spatial/Cities_Towns.shp")

towns <- st_centroid(towns)
