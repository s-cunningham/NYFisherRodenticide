library(tidyverse)
library(lme4)
library(broom.mixed)
library(MuMIn)

set.seed(123)

#### Read in data ####
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

#### Analysis ####

## set up binary variable
dat$binary.MO <- ifelse(dat$n.compounds.MO==0, 0, 1)
dat$binary.T <- ifelse(dat$n.compounds.T==0, 0, 1)

## Scale and center variables
dat[,c(9,18:23)] <- scale(dat[,c(9,18:23)])
dat$fBMI <- ordered(dat$year, levels=c(2019, 2018, 2020))

## Subset data by covariates
dat$buffsize <- as.character(dat$buffsize)
dat$buffsize <- gsub("4.5", "4_5", dat$buffsize)

## Use all the data

## Percent AG
pctAG1 <- dat[, c(1:16, 18, 25:27)]
pctAG1 <- pctAG1 %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=pct_ag, values_fn=unique) %>% as.data.frame()
names(pctAG1)[19:22] <- c("4_5km2", "15km2", "30km2", "60km2") 

pctAG_45 <- glmer(binary.T ~ `4_5km2` + (Region|WMU), family=binomial(link="logit"), data=pctAG1)
pctAG_15 <- glmer(binary.T ~ `15km2` + (Region|WMU), family=binomial(link="logit"), data=pctAG1)
pctAG_30 <- glmer(binary.T ~ `30km2` + (Region|WMU), family=binomial(link="logit"), data=pctAG1)
pctAG_60 <- glmer(binary.T ~ `60km2` + (Region|WMU), family=binomial(link="logit"), data=pctAG1)
pctAG_sel <- model.sel(pctAG_45, pctAG_15, pctAG_30, pctAG_60)
pctAG_sel

## Beech basal area
baa1 <- dat[, c(1:16, 19, 25:27)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=baa, values_fn=unique) %>% as.data.frame()
names(baa1)[19:22] <- c("4_5km2", "15km2", "30km2", "60km2") 

baa_45 <- glmer(binary.T ~ `4_5km2` + (Region|WMU), family=binomial(link="logit"), data=baa1)
baa_15 <- glmer(binary.T ~ `15km2` + (Region|WMU), family=binomial(link="logit"), data=baa1)
baa_30 <- glmer(binary.T ~ `30km2` + (Region|WMU), family=binomial(link="logit"), data=baa1)
baa_60 <- glmer(binary.T ~ `60km2` + (Region|WMU), family=binomial(link="logit"), data=baa1)

# Model selection with MuMIn
baa_sel <- model.sel(baa_45, baa_15, baa_30, baa_60)
baa_sel

## Intermix WUI
intermix1 <- dat[, c(1:17, 20, 25:27)]
intermix1 <- unite(intermix1, "buffrad", 16:17, sep="_")
intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[19:30] <- c("4_5km2_100m", "15km2_100m", "30km2_100m", "60km2_100m",
                             "4_5km2_250m", "15km2_250m", "30km2_250m", "60km2_250m",
                             "4_5km2_500m", "15km2_500m", "30km2_500m", "60km2_500m") 

mix_45100 <- glmer(binary.T ~ `4_5km2_100m` + (Region|WMU), family=binomial(link="logit"), data=intermix1)
mix_45250 <- glmer(binary.T ~ `4_5km2_250m` + (Region|WMU), family=binomial(link="logit"), data=intermix1)
mix_45500 <- glmer(binary.T ~ `4_5km2_500m` + (Region|WMU), family=binomial(link="logit"), data=intermix1)

mix_15100 <- glmer(binary.T ~ `15km2_100m` + (Region|WMU), family=binomial(link="logit"), data=intermix1)
mix_15250 <- glmer(binary.T ~ `15km2_250m` + (Region|WMU), family=binomial(link="logit"), data=intermix1)
mix_15500 <- glmer(binary.T ~ `15km2_500m` + (Region|WMU), family=binomial(link="logit"), data=intermix1)

mix_30100 <- glmer(binary.T ~ `30km2_100m` + (Region|WMU), family=binomial(link="logit"), data=intermix1)
mix_30250 <- glmer(binary.T ~ `30km2_250m` + (Region|WMU), family=binomial(link="logit"), data=intermix1)
mix_30500 <- glmer(binary.T ~ `30km2_500m` + (Region|WMU), family=binomial(link="logit"), data=intermix1)

mix_60100 <- glmer(binary.T ~ `60km2_100m` + (Region|WMU), family=binomial(link="logit"), data=intermix1)
mix_60250 <- glmer(binary.T ~ `60km2_250m` + (Region|WMU), family=binomial(link="logit"), data=intermix1)
mix_60500 <- glmer(binary.T ~ `60km2_500m` + (Region|WMU), family=binomial(link="logit"), data=intermix1)

intermix_sel <- model.sel(mix_45100, mix_15100, mix_30100,mix_60100,
                           mix_45250, mix_15250, mix_30250,mix_60250,
                           mix_45500, mix_15500, mix_30500,mix_60500)
intermix_sel

## Interface WUI
interface1 <- dat[, c(1:17, 21, 25:27)]
interface1 <- unite(interface1, "buffrad", 16:17, sep="_")
interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[19:30] <- c("4_5km2_100m", "15km2_100m", "30km2_100m", "60km2_100m",
                              "4_5km2_250m", "15km2_250m", "30km2_250m", "60km2_250m",
                              "4_5km2_500m", "15km2_500m", "30km2_500m", "60km2_500m") 

face_45100 <- glmer(binary.T ~ `4_5km2_100m` + (Region|WMU), family=binomial(link="logit"), data=interface1)
face_45250 <- glmer(binary.T ~ `4_5km2_250m` + (Region|WMU), family=binomial(link="logit"), data=interface1)
face_45500 <- glmer(binary.T ~ `4_5km2_500m` + (Region|WMU), family=binomial(link="logit"), data=interface1)

face_15100 <- glmer(binary.T ~ `15km2_100m` + (Region|WMU), family=binomial(link="logit"), data=interface1)
face_15250 <- glmer(binary.T ~ `15km2_250m` + (Region|WMU), family=binomial(link="logit"), data=interface1)
face_15500 <- glmer(binary.T ~ `15km2_500m` + (Region|WMU), family=binomial(link="logit"), data=interface1)

face_30100 <- glmer(binary.T ~ `30km2_100m` + (Region|WMU), family=binomial(link="logit"), data=interface1)
face_30250 <- glmer(binary.T ~ `30km2_250m` + (Region|WMU), family=binomial(link="logit"), data=interface1)
face_30500 <- glmer(binary.T ~ `30km2_500m` + (Region|WMU), family=binomial(link="logit"), data=interface1)

face_60100 <- glmer(binary.T ~ `60km2_100m` + (Region|WMU), family=binomial(link="logit"), data=interface1)
face_60250 <- glmer(binary.T ~ `60km2_250m` + (Region|WMU), family=binomial(link="logit"), data=interface1)
face_60500 <- glmer(binary.T ~ `60km2_500m` + (Region|WMU), family=binomial(link="logit"), data=interface1)

interface_sel <- model.sel(face_45100, face_15100, face_30100, face_60100,
                            face_45250, face_15250, face_30250, face_60250,
                            face_45500, face_15500, face_30500, face_60500)
interface_sel

## Total WUI
wui1 <- dat[, c(1:17, 22, 25:27)]
wui1 <- unite(wui1, "buffrad", 16:17, sep="_")
wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[19:30] <- c("4_5km2_100m", "15km2_100m", "30km2_100m", "60km2_100m",
                        "4_5km2_250m", "15km2_250m", "30km2_250m", "60km2_250m",
                        "4_5km2_500m", "15km2_500m", "30km2_500m", "60km2_500m") 

# wui_45100 <- glmer(binary.T ~ `4_5km2_100m` + (1|Region), family=binomial(link="logit"), data=wui1)
# wui_45250 <- glmer(binary.T ~ `4_5km2_250m` + (1|Region), family=binomial(link="logit"), data=wui1)
# wui_45500 <- glmer(binary.T ~ `4_5km2_500m` + (1|Region), family=binomial(link="logit"), data=wui1)
# 
# wui_15100 <- glmer(binary.T ~ `15km2_100m` + (1|Region), family=binomial(link="logit"), data=wui1)
# wui_15250 <- glmer(binary.T ~ `15km2_250m` + (1|Region), family=binomial(link="logit"), data=wui1)
# wui_15500 <- glmer(binary.T ~ `15km2_500m` + (1|Region), family=binomial(link="logit"), data=wui1)
# 
# wui_30100 <- glmer(binary.T ~ `30km2_100m` + (1|Region), family=binomial(link="logit"), data=wui1)
# wui_30250 <- glmer(binary.T ~ `30km2_250m` + (1|Region), family=binomial(link="logit"), data=wui1)
# wui_30500 <- glmer(binary.T ~ `30km2_500m` + (1|Region), family=binomial(link="logit"), data=wui1)
# 
# wui_60100 <- glmer(binary.T ~ `60km2_100m` + (1|Region), family=binomial(link="logit"), data=wui1)
# wui_60250 <- glmer(binary.T ~ `60km2_250m` + (1|Region), family=binomial(link="logit"), data=wui1)
# wui_60500 <- glmer(binary.T ~ `60km2_500m` + (1|Region), family=binomial(link="logit"), data=wui1)
# 
# wui_sel <- model.sel(wui_45100, wui_15100, wui_30100, wui_60100, 
#                       wui_45250, wui_15250, wui_30250, wui_60250,  
#                       wui_45500, wui_15500, wui_30500, wui_60500)
# wui_sel

#### Set up data to run for each combination of covariates ####
dat1 <- dat[,c(1:15,23:27)]
dat1 <- distinct(dat1)

# Join percent agriculture
pctAG1 <- pctAG1[,c(1:3, 18:21)]
dat1 <- left_join(dat1, pctAG1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[21:24] <- c("pctAG_4_5", "pctAG_15", "pctAG_30", "pctAG_60")

# Join beech basal area
baa1 <- baa1[,c(1:3, 18:21)]
dat1 <- left_join(dat1, baa1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[25:28] <- c("baa_4_5", "baa_15", "baa_30", "baa_60")

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 18:29)]
dat1 <- left_join(dat1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[29:40] <- c("mix_4_5km2_100m", "mix_15km2_100m", "mix_30km2_100m", "mix_60km2_100m",
                        "mix_4_5km2_250m", "mix_15km2_250m", "mix_30km2_250m", "mix_60km2_250m",
                        "mix_4_5km2_500m", "mix_15km2_500m", "mix_30km2_500m", "mix_60km2_500m") 

# Join interface WUI
interface1 <- interface1[, c(1:3, 18:29)]
dat1 <- left_join(dat1, interface1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[41:52] <- c("face_4_5km2_100m", "face_15km2_100m", "face_30km2_100m", "face_60km2_100m",
                        "face_4_5km2_250m", "face_15km2_250m", "face_30km2_250m", "face_60km2_250m",
                        "face_4_5km2_500m", "face_15km2_500m", "face_30km2_500m", "face_60km2_500m") 

# Join total WUI
wui1 <- wui1[,c(1:3, 18:29)]
dat1 <- left_join(dat1, wui1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[53:64] <- c("wui_4_5km2_100m", "wui_15km2_100m", "wui_30km2_100m", "wui_60km2_100m",
                        "wui_4_5km2_250m", "wui_15km2_250m", "wui_30km2_250m", "wui_60km2_250m",
                        "wui_4_5km2_500m", "wui_15km2_500m", "wui_30km2_500m", "wui_60km2_500m") 

## Dropping 4.5 km2 buffer
dat1 <- dat1[,-c(21,25,29,33,37,41,45,49,53,57,61)]

## square variables
dat1 <- mutate(dat1, pctAG_15sq=pctAG_15^2)
dat1 <- mutate(dat1, pctAG_30sq=pctAG_30^2)
dat1 <- mutate(dat1, pctAG_60sq=pctAG_60^2)

dat1 <- mutate(dat1, mix_15_100sq=mix_15km2_100m^2)
dat1 <- mutate(dat1, mix_30_100sq=mix_30km2_100m^2)
dat1 <- mutate(dat1, mix_60_100sq=mix_60km2_100m^2)

dat1 <- mutate(dat1, mix_15_250sq=mix_15km2_250m^2)
dat1 <- mutate(dat1, mix_30_250sq=mix_30km2_250m^2)
dat1 <- mutate(dat1, mix_60_250sq=mix_60km2_250m^2)

dat1 <- mutate(dat1, mix_15_500sq=mix_15km2_500m^2)
dat1 <- mutate(dat1, mix_30_500sq=mix_30km2_500m^2)
dat1 <- mutate(dat1, mix_60_500sq=mix_60km2_500m^2)

dat1 <- mutate(dat1, face_15_100sq=face_15km2_100m^2)
dat1 <- mutate(dat1, face_30_100sq=face_30km2_100m^2)
dat1 <- mutate(dat1, face_60_100sq=face_60km2_100m^2)

dat1 <- mutate(dat1, face_15_250sq=face_15km2_250m^2)
dat1 <- mutate(dat1, face_30_250sq=face_30km2_250m^2)
dat1 <- mutate(dat1, face_60_250sq=face_60km2_250m^2)

dat1 <- mutate(dat1, face_15_500sq=face_15km2_500m^2)
dat1 <- mutate(dat1, face_30_500sq=face_30km2_500m^2)
dat1 <- mutate(dat1, face_60_500sq=face_60km2_500m^2)

dat1 <- mutate(dat1, wui_15_100sq=wui_15km2_100m^2)
dat1 <- mutate(dat1, wui_30_100sq=wui_30km2_100m^2)
dat1 <- mutate(dat1, wui_60_100sq=wui_60km2_100m^2)

dat1 <- mutate(dat1, wui_15_250sq=wui_15km2_250m^2)
dat1 <- mutate(dat1, wui_30_250sq=wui_30km2_250m^2)
dat1 <- mutate(dat1, wui_60_250sq=wui_60km2_250m^2)

dat1 <- mutate(dat1, wui_15_500sq=wui_15km2_500m^2)
dat1 <- mutate(dat1, wui_30_500sq=wui_30km2_500m^2)
dat1 <- mutate(dat1, wui_60_500sq=wui_60km2_500m^2)

#### Run models ####

# Set up lists to hold data
pt1 <- pt2 <- pt3 <- pt4 <- pt5 <- pt6 <- pt7 <- pt8 <- pt9 <- pt10 <- list()

for (i in 1:10) {
  
  pt_subset <- dat1[dat1$pt_index==i,]
  
  # Set up model with all covariates
  global <- glmer(binary.T ~ Sex*Age +
                 baa_15*laggedBMI + baa_30*laggedBMI + baa_60*laggedBMI + 
                 pctAG_15 + pctAG_15sq + 
                 pctAG_30 + pctAG_30sq + 
                 pctAG_60 + pctAG_60sq + 
                 # mix_15km2_100m + mix_30km2_100m + mix_60km2_100m +
                 # mix_15km2_250m + mix_30km2_250m + mix_60km2_250m +
                 # mix_15km2_500m + mix_30km2_500m + mix_60km2_500m +
                 # face_15km2_100m + face_30km2_100m + face_60km2_100m +
                 # face_15km2_250m + face_30km2_250m + face_60km2_250m +
                 # face_15km2_500m + face_30km2_500m + face_60km2_500m +
                 wui_15km2_100m + wui_30km2_100m + wui_60km2_100m +
                 wui_15km2_250m + wui_30km2_250m + wui_60km2_250m +
                 wui_15km2_500m + wui_30km2_500m + wui_60km2_500m +
                 # mix_15_100sq + mix_30_100sq + mix_60_100sq +
                 # mix_15_250sq + mix_30_250sq + mix_60_250sq +
                 # mix_15_500sq + mix_30_500sq + mix_60_500sq +
                 # face_15_100sq + face_30_100sq + face_60_100sq +
                 # face_15_250sq + face_30_250sq + face_60_250sq +
                 # face_15_500sq + face_30_500sq + face_60_500sq +
                 # wui_15_100sq + wui_30_100sq + wui_60_100sq +
                 # wui_15_250sq + wui_30_250sq + wui_60_250sq +
                 # wui_15_500sq + wui_30_500sq + wui_60_500sq +
                 (1|Region), family=binomial(link="logit"), data=pt_subset)
  
  # Use dredge function to build models with some serious subsetting
  options(na.action = "na.fail")
  models <- dredge(global)
  subset(models, delta < 6)
}


