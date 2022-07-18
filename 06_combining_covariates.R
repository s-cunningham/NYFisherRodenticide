## Set up rodenticide analysis
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

## Subset ag and add together
ag <- ag[ag$value==81 | ag$value==82, -1]
c_ag <- ag %>% group_by(name, buffsize) %>% summarise(pct_ag=sum(freq))
names(c_ag)[1] <- "pt_name"

## reorganize beech
baa <- baa[,-1]
names(baa)[1] <- "pt_name"

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

# annual beech mast index
mast <- read_csv("data/analysis-ready/ALTEMP26_beech-data.csv")
mast_mean <- mean(mast$Total_Beechnuts)
mast_median <- median(mast$Total_Beechnuts)
mast_max <- max(mast$Total_Beechnuts)
mast$devMean <- mast$Total_Beechnuts - mast_mean

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
dat <- left_join(dat, c_ag, by=c("pt_name", "buffsize"))
dat <- left_join(dat, baa, by=c("pt_name", "buffsize"))
dat <- left_join(dat, wui, by=c("pt_name", "buffsize", "radius"))

# Add beech mast index
dat$laggedBMI <- NA
dat$laggedBMI[dat$year==2018] <- mast$devMean[mast$year==2017]
dat$laggedBMI[dat$year==2019] <- mast$devMean[mast$year==2018]
dat$laggedBMI[dat$year==2020] <- mast$devMean[mast$year==2019]

# Fill in 0 for missing values (ag, WUI)
dat <- dat %>% mutate(pct_ag = coalesce(pct_ag, 0),
                      interface = coalesce(interface, 0),
                      intermix = coalesce(intermix, 0),
                      totalWUI = coalesce(totalWUI, 0))

# Add WMUA 
dat <- left_join(dat, wmua, by="WMU")

### Save data to file ####
write_csv(dat, "data/analysis-ready/combined_AR_covars.csv")


#### Semivariogram ####
dat1 <- dat1[,-c(2:4)]
dat1 <- unique(dat1)
sp::coordinates(dat1) <- ~x_coord+y_coord
vario <- gstat::variogram(n.compounds.MO~1, data=dat1)
plot(vario)

# create correlologram


