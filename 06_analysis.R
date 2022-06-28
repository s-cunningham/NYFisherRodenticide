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
wui250 <- read_csv("data/analysis-ready/wui250_frac.csv")
wui250 <- wui250[wui250$value==1 | wui250$value==2, -1]
wui500 <- read_csv("data/analysis-ready/wui500_frac.csv")
wui500 <- wui500[wui500$value==1 | wui500$value==2, -1]
ag <- read_csv("data/analysis-ready/nlcd_pct.csv")
baa <- read_csv("data/analysis-ready/baa-predicted_avg.csv")
pts <- read_csv("output/random_point_locs.csv")

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
baa <- baa[,c(4,2,3)]
names(baa)[1] <- "pt_name"

## reorganize WUI
wui100$radius <- 100
wui250$radius <- 250
wui500$radius <- 500
wui <- bind_rows(wui100, wui250, wui500)
wui <- wui[complete.cases(wui),]
totwui <- wui %>% group_by(name, buffsize, radius) %>% summarise(total_wui=sum(freq))
wui <- pivot_wider(wui, names_from=value, values_from=freq)

## Joining data
# Join location, age & sex details to random points
dets <- left_join(pts, dets, by="RegionalID")
dat <- left_join(dets, ncomps, by="RegionalID")
dat <- separate(dat, 2, into=c("id", "pt_index"), sep="_", remove=FALSE)
dat <- dat[,-c(3)]
dat$pt_index <- as.numeric(dat$pt_index)

# join covariate data
dat <- left_join(dat, c_ag, by="pt_name")
dat <- left_join(dat, baa, by=c("pt_name", "buffsize"))



#### Semivariogram ####




