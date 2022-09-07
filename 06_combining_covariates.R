## Set up rodenticide covariates
## 2022-06-27

library(tidyverse)

#### Read in data ####
dat <- read_csv("output/summarized_AR_results.csv")
trace_y <- read_csv("output/ncompounds_trace.csv")
trace_y <- trace_y[,-c(1,4,5)]
trace_n <- read_csv("output/ncompounds_notrace.csv")
trace_n <- trace_n[,-c(1,4,5)]
wui100 <- read_csv("data/analysis-ready/wui100_frac.csv")
wui100 <- wui100[wui100$value==1 | wui100$value==2, -1]
wui100 <- wui100[complete.cases(wui100),]
wui250 <- read_csv("data/analysis-ready/wui250_frac.csv")
wui250 <- wui250[wui250$value==1 | wui250$value==2, -1]
wui250 <- wui250[complete.cases(wui250),]
wui500 <- read_csv("data/analysis-ready/wui500_frac.csv")
wui500 <- wui500[wui500$value==1 | wui500$value==2, -1]
wui500 <- wui500[complete.cases(wui500),]
ag <- read_csv("data/analysis-ready/nlcd_pct.csv")
baa <- read_csv("data/analysis-ready/baa_mean.csv")
pts <- read_csv("output/random_point_locs.csv")
wmua <- read_csv("data/analysis-ready/wmuas.csv")

#### Combine data ####
## Number of compounds detected
names(trace_y)[2] <- "n.compounds.T"
names(trace_n)[2] <- "n.compounds.MO"
ncomps <- left_join(trace_n, trace_y, by="RegionalID")

## Save details of each sample
dets <- dat[,c(2,6:9,15:18)]
dets <- unique(dets)

## Add a column to points with just sample ID
pts <- separate(pts, 1, into=c("RegionalID", "throw"), sep="_", remove=FALSE)
pts <- pts[,c(2,1,4,5)]
names(pts)[3:4] <- c("rand_x", "rand_y")

## Subset forest and add together
ag <- ag[,-1]
forest <- ag[ag$value==41 | ag$value==42 | ag$value==43,]
forest$value <- factor(forest$value, levels=c(41, 42, 43), labels=c("deciduous", "evergreen", "mixed")) |> as.character()
forest <- forest[forest$buffsize>10,]
forest <- forest[,c(1,4,2,3)]
forest <- pivot_wider(forest, names_from=value, values_from=freq)
names(forest)[1] <- "pt_name"
forest$deciduous[is.na(forest$deciduous)] <- 0
forest$evergreen[is.na(forest$evergreen)] <- 0
forest$mixed[is.na(forest$mixed)] <- 0
forest <- mutate(forest, totalforest=deciduous + evergreen + mixed)

## Subset ag and add together
ag <- ag[ag$value==81 | ag$value==82, ]
ag$value <- factor(ag$value, levels=c(81, 82), labels=c("pasture", "crops")) |> as.character()
ag <- ag[ag$buffsize>10,]
names(ag)[1] <- "pt_name"
ag <- pivot_wider(ag, names_from=value, values_from=freq)
ag$crops[is.na(ag$crops)] <- 0
ag$pasture[is.na(ag$pasture)] <- 0
ag <- mutate(ag, totalag=crops + pasture)

## reorganize beech
names(baa)[1] <- "pt_name"
baa$baa[is.na(baa$baa)] <- 0

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
wui <- wui[wui$buffsize>10,]

## Joining data
# Join location, age & sex details to random points
dets <- left_join(pts, dets, by="RegionalID")
dat <- left_join(dets, ncomps, by="RegionalID")
# dat1 <- left_join(dets, ncomps, by="RegionalID")
dat <- separate(dat, 2, into=c("id", "pt_index"), sep="_", remove=FALSE)
dat <- dat[,-c(3)]
dat$pt_index <- as.numeric(dat$pt_index)

# Add columns for buffer and radius
dat <- bind_rows(dat, dat, dat)
dat$buffsize <- rep(c(15,30,60), each=3370) # buffer sizes
dat <- bind_rows(dat, dat, dat)
dat$radius <- rep(c(100,250,500), each=10110) # WUI radius sizes

# join covariate data
dat <- left_join(dat, ag, by=c("pt_name", "buffsize"))
dat <- left_join(dat, forest, by=c("pt_name", "buffsize"))
dat <- left_join(dat, wui, by=c("pt_name", "buffsize", "radius"))

# Add beech mast index
dat$laggedYear <- dat$year - 1
baa$laggedYear <- baa$year - 1

dat <- left_join(dat, baa[,c(1,3,4,5)], by=c("pt_name", "buffsize", "laggedYear"))
names(dat)[29] <- "laggedBMI"

dat <- left_join(dat, baa[,c(1,2,3,4)], by=c("pt_name", "buffsize", "year"))
names(dat)[30] <- "BMI"

# Fill in 0 for missing values (WUI)
dat <- dat %>% mutate(interface = coalesce(interface, 0),
                      intermix = coalesce(intermix, 0),
                      totalWUI = coalesce(totalWUI, 0),
                      crops = coalesce(crops, 0),
                      pasture = coalesce(pasture, 0),
                      totalag = coalesce(totalag, 0))

# Add WMUA 
dat <- left_join(dat, wmua, by="WMU")

# Reorder columns
dat <- dat[,c(1:6,32,7:9,13:21,25:27,29,30)]

### Save data to file ####
write_csv(dat, "data/analysis-ready/combined_AR_covars.csv")


#### Semivariogram ####
dat1 <- dat1[,-c(2:4)]
dat1 <- unique(dat1)
sp::coordinates(dat1) <- ~x_coord+y_coord
vario <- gstat::variogram(n.compounds.MO~1, data=dat1)
plot(vario)

## Moran's I
midat <- dat[dat$pt_index==1,c(1,4:10,14:15)]
midat$binary.T <- ifelse(midat$n.compounds.T==0, 0, 1)
midat$binary.MO <- ifelse(midat$n.compounds.MO==0, 0, 1)
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

ggplot(dat100, aes(x=factor(buffsize), y=baa, fill=factor(Region))) + geom_boxplot() + theme_bw() 


ggplot(dat100, aes(x=factor(buffsize), y=pct_ag, fill=factor(Region))) + geom_boxplot() + theme_bw()
ggplot(dat100, aes(x=factor(buffsize), y=intermix, fill=factor(Region))) + geom_boxplot() + theme_bw()
ggplot(dat100, aes(x=factor(buffsize), y=interface, fill=factor(Region))) + geom_boxplot() + theme_bw()

ggplot(dat100, aes(x=factor(buffsize), y=totalWUI)) + geom_boxplot() + theme_bw()
ggplot(dat250, aes(x=factor(buffsize), y=totalWUI)) + geom_boxplot() + theme_bw()
ggplot(dat500, aes(x=factor(buffsize), y=totalWUI, fill=factor(Region))) + 
  geom_boxplot() +  theme_bw()
