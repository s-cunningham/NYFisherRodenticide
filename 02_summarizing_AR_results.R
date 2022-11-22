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
dat$WMU[dat$Town=="Freedom" & dat$WMU=="9W"] <- "9T"
dat$WMU[dat$Town=="Western" & dat$WMU=="5H"] <- "6K"

dat <- dat %>% mutate_at(c(10:20), as.numeric) %>% 
  unite("key", c(WMU,Town), sep="-", remove=FALSE) %>%
  arrange(key)

unique(dat$key)

# Check if there are compounds not detected
for (i in 10:20) {
  print(sum(is.na(dat[,i])))
}

# Convert to long format for ggplotting
datl <- as.data.frame(pivot_longer(dat, cols=11:21, names_to="compound", values_to="ppm"))

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

