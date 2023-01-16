## Set up rodenticide covariates
## 2022-06-27, updated 2023-01-12

library(tidyverse)

#### Read in data ####
dat <- read_csv("output/summarized_AR_results.csv")
trace_y <- read_csv("output/ncompounds_trace.csv") %>% 
  select(RegionalID, n.compounds) %>%
  rename(n.compounds.T=n.compounds)
trace_n <- read_csv("output/ncompounds_notrace.csv") %>% 
  select(RegionalID, n.compounds) %>%
  rename(n.compounds.MO=n.compounds)
wui100 <- read_csv("data/analysis-ready/wui100_frac.csv") %>%
  filter(value==1 | value==2) %>%
  filter(complete.cases(.))
wui250 <- read_csv("data/analysis-ready/wui250_frac.csv") %>%
  filter(value==1 | value==2) %>%
  filter(complete.cases(.))
wui500 <- read_csv("data/analysis-ready/wui500_frac.csv") %>%
  filter(value==1 | value==2) %>%
  filter(complete.cases(.))
ag <- read_csv("data/analysis-ready/nlcd_pct.csv")
bmi <- read_csv("data/analysis-ready/baa_sum.csv")
baa <- read_csv("data/analysis-ready/baa_sum_single_raster.csv")
pts <- read_csv("output/random_point_locs.csv")
wmua <- read_csv("data/analysis-ready/wmuas.csv")
# lsm <- read_csv("data/analysis-ready/forest_lsm_2.csv")
build <- read_csv("data/analysis-ready/building-centroid_sum.csv") %>%
            rename(pt_name=name) 
mast <- read_csv("data/analysis-ready/ALTEMP26_beech-data.csv")

#### Combine data ####
## Number of compounds detected
ncomps <- left_join(trace_n, trace_y, by="RegionalID")

## Save details of each sample
dets <- dat %>% select(RegionalID:WMU,Town) %>% distinct()

## Add a column to points with just sample ID
pts <- pts %>% select(RegionalID,name,x,y) %>%
          rename(rand_x=x, rand_y=y, pt_name=name)

## Subset forest and add together
forest <- ag %>% filter(value==41 | value==42 | value==43)   
forest$value <- factor(forest$value, levels=c(41, 42, 43), labels=c("deciduous", "evergreen", "mixed")) |> as.character()
forest <- forest %>% select(name,buffsize,value,freq) %>%
           pivot_wider(names_from=value, values_from=freq) %>%
           rename(pt_name=name)
forest$deciduous[is.na(forest$deciduous)] <- 0
forest$evergreen[is.na(forest$evergreen)] <- 0
forest$mixed[is.na(forest$mixed)] <- 0
forest <- mutate(forest, totalforest=deciduous + evergreen + mixed)

## Subset ag and add together
ag <- ag %>% filter(value==81 | value==82)
ag$value <- factor(ag$value, levels=c(81, 82), labels=c("pasture", "crops")) %>% as.character()
ag <- ag %>% pivot_wider(names_from=value, values_from=freq) %>%
        rename(pt_name=name)
ag$crops[is.na(ag$crops)] <- 0
ag$pasture[is.na(ag$pasture)] <- 0
ag <- mutate(ag, totalag=crops + pasture)

## reorganize beech
# beech mast index
names(bmi)[c(1,3)] <- c("pt_name", "bmi")
bmi$bmi[is.na(bmi$bmi)] <- 0

# beech basal area
names(baa)[1] <- "pt_name"
baa$baa[is.na(baa$baa)] <- 0
baa <- baa %>% select(pt_name, buffsize, baa)

## reorganize WUI
wui100$radius <- 100
wui250$radius <- 250
wui500$radius <- 500
wui <- bind_rows(wui100, wui250, wui500)
wui <- wui[complete.cases(wui),]
wui <- pivot_wider(wui, names_from=value, values_from=freq)
names(wui)[4:5] <- c("intermix", "interface")
wui[is.na(wui)] <- 0
wui <- mutate(wui, totalWUI=intermix + interface)
names(wui)[1] <- "pt_name"

# reorganize landscape metrics
lsm <- lsm %>% 
        filter(metric!="pafrac") %>%
        mutate(buffer=replace(buffer, buffer==4370.200, 15)) %>%
        mutate(buffer=replace(buffer, buffer==4370.200, 15)) %>%
        mutate(buffer=replace(buffer, buffer==6180.38, 30)) %>%
        mutate(buffer=replace(buffer, buffer==8740.388, 60)) %>%
        pivot_wider(names_from=metric, values_from=value) %>%
        rename(buffsize=buffer, pt_name=plot_id) 
        
# Create categorical variable for buildings
build <- build %>% mutate(build_cat=case_when(
                             nbuildings<1 ~ "None",
                             (nbuildings>=1 & nbuildings<=60) & buffsize==15 ~ "1stQuart",
                             (nbuildings>=1 & nbuildings<=149) & buffsize==30 ~ "1stQuart",
                             (nbuildings>=1 & nbuildings<=339) & buffsize==60 ~ "1stQuart",
                             (nbuildings>60 & nbuildings<=132) & buffsize==15 ~ "2ndQuart",
                             (nbuildings>149 & nbuildings<=287) & buffsize==30 ~ "2ndQuart",
                             (nbuildings>339 & nbuildings<=610) & buffsize==60 ~ "2ndQuart",
                             (nbuildings>132 & nbuildings<=241) & buffsize==15 ~ "3rdQuart",
                             (nbuildings>287 & nbuildings<=503) & buffsize==30 ~ "3rdQuart",
                             (nbuildings>610 & nbuildings<=1094) & buffsize==60 ~ "3rdQuart",
                             (nbuildings>241 & nbuildings<=5464) & buffsize==15 ~ "4thQuart",
                             (nbuildings>503 & nbuildings<=8309) & buffsize==30 ~ "4thQuart",
                             (nbuildings>1094 & nbuildings<=10798) & buffsize==60 ~ "4thQuart"))

## Joining data
# Join location, age & sex details to random points
dets <- left_join(pts, dets, by="RegionalID")
dat <- left_join(dets, ncomps, by="RegionalID")
dat <- separate(dat, 2, into=c("id", "pt_index"), sep="_", remove=FALSE) 
dat <- dat[,-c(3)]
dat$pt_index <- as.numeric(dat$pt_index)
dat <- dat %>% distinct()

# Add columns for buffer and radius
dat <- bind_rows(dat, dat, dat)
dat$buffsize <- rep(c(15,30,60), each=3380) # buffer sizes
dat <- bind_rows(dat, dat, dat)
dat$radius <- rep(c(100,250,500), each=10140) # WUI radius sizes
dat <- dat %>% select(RegionalID:n.compounds.MO, n.compounds.T,buffsize,radius)

# join covariate data
dat <- left_join(dat, ag, by=c("pt_name", "buffsize")) %>%
  left_join(forest, by=c("pt_name", "buffsize")) %>%
  left_join(wui, by=c("pt_name", "buffsize", "radius")) %>%
  # left_join(lsm, by=c("pt_name", "buffsize")) %>%
  left_join(build, by=c("pt_name", "buffsize")) %>%
  left_join(baa, by=c("pt_name", "buffsize"))

## Add beech mast index
dat <- left_join(dat, bmi, by=c("pt_name", "buffsize", "year")) %>%
          rename(BMI=bmi)

# add 1 to year to get lagged
bmi$year <- bmi$year + 1
dat <- left_join(dat, bmi, by=c("pt_name", "buffsize", "year")) %>%
  rename(laggedBMI=bmi)

## Fill in 0 for missing values (WUI)
dat <- dat %>% mutate(interface = coalesce(interface, 0),
                      intermix = coalesce(intermix, 0),
                      totalWUI = coalesce(totalWUI, 0),
                      crops = coalesce(crops, 0),
                      pasture = coalesce(pasture, 0),
                      totalag = coalesce(totalag, 0))

# Add WMUA 
dat <- left_join(dat, wmua, by="WMU")

# Reorder columns
dat <- dat %>% select(RegionalID:AgeClass,key,Region,WMUA_code,WMU,Town:laggedBMI)

# Add beechnut counts
dat <- dat %>% mutate(beechnuts=case_when(
                          year==2018 ~ 6,
                          year==2019 ~ 295,
                          year==2020 ~ 14),
                      lag_beechnuts=case_when(
                          year==2018 ~ 145,
                          year==2019 ~ 6,
                          year==2020 ~ 295))


# mast_mean <- mean(mast$Total_Beechnuts)
# mast_median <- median(mast$Total_Beechnuts)
# mast_max <- max(mast$Total_Beechnuts)
# mast$devMean <- mast$Total_Beechnuts - mast_mean
# mast$devMedian <- mast$Total_Beechnuts - mast_median

### Save data to file ####
write_csv(dat, "data/analysis-ready/combined_AR_covars.csv")


ggplot(dat, aes(x=rand_x, y=rand_y, color=factor(Region))) + geom_point() + theme_bw()


#### Semivariogram ####
dat1 <- dat %>% select(RegionalID,rand_x,rand_y, Region, n.compounds.T)
dat1 <- unique(dat1)
sp::coordinates(dat1) <- ~rand_x+rand_y
vario <- gstat::variogram(n.compounds.T~1, data=dat1)
plot(vario)

## Moran's I
midat <- dat[dat$pt_index==1,c(1,4:10,16)]
midat$binary.T <- ifelse(midat$n.compounds.T==0, 0, 1)
midat <- distinct(midat)

# Create distance matrix
ar.dists <- as.matrix(dist(cbind(midat$rand_x, midat$rand_y)))

# create inverse distance matrix
ar.dists.inv <- 1/ar.dists
diag(ar.dists.inv) <- 0
ar.dists.inv[1:5, 1:5]

# Moran's I
ape::Moran.I(midat$n.compounds.T, ar.dists.inv)

# Binary distance matrix, where d=50 km
ar.dists.bin <- (ar.dists >0 & ar.dists <= 50000)
ape::Moran.I(midat$n.compounds.T, ar.dists.bin)

#### Plot covariate values ####
dat100 <- dat[dat$radius==100,]
dat250 <- dat[dat$radius==250,]
dat500 <- dat[dat$radius==500,]

ggplot(dat100, aes(x=factor(buffsize), y=bmi, fill=factor(Region))) + geom_boxplot() + theme_bw() 


ggplot(dat100, aes(x=factor(buffsize), y=pct_ag, fill=factor(Region))) + geom_boxplot() + theme_bw()
ggplot(dat100, aes(x=factor(buffsize), y=intermix, fill=factor(Region))) + geom_boxplot() + theme_bw()
ggplot(dat100, aes(x=factor(buffsize), y=interface, fill=factor(Region))) + geom_boxplot() + theme_bw()

ggplot(dat100, aes(x=factor(buffsize), y=totalWUI)) + geom_boxplot() + theme_bw()
ggplot(dat250, aes(x=factor(buffsize), y=totalWUI)) + geom_boxplot() + theme_bw()
ggplot(dat500, aes(x=factor(buffsize), y=totalWUI, fill=factor(Region))) + 
  geom_boxplot() +  theme_bw()
