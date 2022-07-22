library(tidyverse)
# library(brms)
# library(lme4)
# library(broom.mixed)
library(MuMIn)
library(ordinal)


set.seed(123)

#### Read in data ####
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

#### Analysis Set-up ####

## Scale and center variables
dat[,c(9,18:23)] <- scale(dat[,c(9,18:23)])

# Order response
dat$n.compounds.MO <- ordered(dat$n.compounds.MO, levels=c(0,1,2,3))
dat$n.compounds.T <- ordered(dat$n.compounds.T, levels=c(0,1,2,3,4,5))

## Percent AG
pctAG1 <- dat[, c(1:16, 18, 25)]
pctAG1 <- pctAG1 %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=pct_ag, values_fn=unique) %>% as.data.frame()
names(pctAG1)[17:20] <- c("km4_5", "km15", "km30", "km60") 

ag4_5 <- clmm(n.compounds.T ~ km4_5 + (1|Region), data=pctAG1)
summary(ag4_5)
ag15 <- clmm(n.compounds.T ~ km15 + (1|Region), data=pctAG1)
summary(ag15)
ag30 <- clmm(n.compounds.T ~ km30 + (1|Region), data=pctAG1)
summary(ag30)
ag60 <- clmm(n.compounds.T ~ km60 + (1|Region), data=pctAG1)
summary(ag60)

pctAG_sel <- model.sel(ag4_5, ag15, ag30, ag60)
pctAG_sel

## Beech basal area
baa1 <- dat[, c(1:16, 19, 25)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=baa, values_fn=unique) %>% as.data.frame()
names(baa1)[17:20] <- c("km4_5", "km15", "km30", "km60") 

baa4_5 <- clmm(n.compounds.T ~ km4_5 + (1|Region), data=baa1)
summary(baa4_5)
baa15 <- clmm(n.compounds.T ~ km15 + (1|Region), data=baa1)
summary(baa15)
baa30 <- clmm(n.compounds.T ~ km30 + (1|Region), data=baa1)
summary(baa30)
baa60 <- clmm(n.compounds.T ~ km60 + (1|Region), data=baa1)
summary(baa60)

baa_sel <- model.sel(baa4_5, baa15, baa30, baa60)
baa_sel


## Intermix WUI
intermix1 <- dat[, c(1:17, 20, 25)]
intermix1 <- unite(intermix1, "buffrad", 16:17, sep="_")
intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[17:28] <- c("b4_5km2_100m", "b15km2_100m", "b30km2_100m", "b60km2_100m",
                             "b4_5km2_250m", "b15km2_250m", "b30km2_250m", "b60km2_250m",
                             "b4_5km2_500m", "b15km2_500m", "b30km2_500m", "b60km2_500m") 

mix_45100 <- clmm(n.compounds.T ~ b4_5km2_100m + (1|Region), data=intermix1)
mix_45250 <- clmm(n.compounds.T ~ b4_5km2_250m + (1|Region), data=intermix1)
mix_45500 <- clmm(n.compounds.T ~ b4_5km2_500m + (1|Region), data=intermix1)

mix_15100 <- clmm(n.compounds.T ~ b15km2_100m + (1|Region), data=intermix1)
mix_15250 <- clmm(n.compounds.T ~ b15km2_250m + (1|Region), data=intermix1)
mix_15500 <- clmm(n.compounds.T ~ b15km2_500m + (1|Region), data=intermix1)

mix_30100 <- clmm(n.compounds.T ~ b30km2_100m + (1|Region), data=intermix1)
mix_30250 <- clmm(n.compounds.T ~ b30km2_250m + (1|Region), data=intermix1)
mix_30500 <- clmm(n.compounds.T ~ b30km2_500m + (1|Region), data=intermix1)

mix_60100 <- clmm(n.compounds.T ~ b60km2_100m + (1|Region), data=intermix1)
mix_60250 <- clmm(n.compounds.T ~ b60km2_250m + (1|Region), data=intermix1)
mix_60500 <- clmm(n.compounds.T ~ b60km2_500m + (1|Region), data=intermix1)

intermix_sel <- model.sel(mix_45100, mix_15100, mix_30100,mix_60100,
                          mix_45250, mix_15250, mix_30250,mix_60250,
                          mix_45500, mix_15500, mix_30500,mix_60500)
intermix_sel

## Interface WUI
interface1 <- dat[, c(1:17, 21, 25)]
interface1 <- unite(interface1, "buffrad", 16:17, sep="_")
interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[17:28] <- c("b4_5km2_100m", "b15km2_100m", "b30km2_100m", "b60km2_100m",
                              "b4_5km2_250m", "b15km2_250m", "b30km2_250m", "b60km2_250m",
                              "b4_5km2_500m", "b15km2_500m", "b30km2_500m", "b60km2_500m") 

face_45100 <- clmm(n.compounds.T ~ b4_5km2_100m + (1|Region), data=interface1)
face_45250 <- clmm(n.compounds.T ~ b4_5km2_250m + (1|Region), data=interface1)
face_45500 <- clmm(n.compounds.T ~ b4_5km2_500m + (1|Region), data=interface1)

face_15100 <- clmm(n.compounds.T ~ b15km2_100m + (1|Region), data=interface1)
face_15250 <- clmm(n.compounds.T ~ b15km2_250m + (1|Region), data=interface1)
face_15500 <- clmm(n.compounds.T ~ b15km2_500m + (1|Region), data=interface1)

face_30100 <- clmm(n.compounds.T ~ b30km2_100m + (1|Region), data=interface1)
face_30250 <- clmm(n.compounds.T ~ b30km2_250m + (1|Region), data=interface1)
face_30500 <- clmm(n.compounds.T ~ b30km2_500m + (1|Region), data=interface1)

face_60100 <- clmm(n.compounds.T ~ b60km2_100m + (1|Region), data=interface1)
face_60250 <- clmm(n.compounds.T ~ b60km2_250m + (1|Region), data=interface1)
face_60500 <- clmm(n.compounds.T ~ b60km2_500m + (1|Region), data=interface1)

interface_sel <- model.sel(face_45100, face_15100, face_30100, face_60100,
                           face_45250, face_15250, face_30250, face_60250,
                           face_45500, face_15500, face_30500, face_60500)
interface_sel

## Total WUI
wui1 <- dat[, c(1:17, 22, 25)]
wui1 <- unite(wui1, "buffrad", 16:17, sep="_")
wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[17:28] <- c("b4_5km2_100m", "b15km2_100m", "b30km2_100m", "b60km2_100m",
                        "b4_5km2_250m", "b15km2_250m", "b30km2_250m", "b60km2_250m",
                        "b4_5km2_500m", "b15km2_500m", "b30km2_500m", "b60km2_500m") 

wui_45100 <- clmm(n.compounds.T ~ b4_5km2_100m + (1|Region), data=wui1)
wui_45250 <- clmm(n.compounds.T ~ b4_5km2_250m + (1|Region), data=wui1)
wui_45500 <- clmm(n.compounds.T ~ b4_5km2_500m + (1|Region), data=wui1)

wui_15100 <- clmm(n.compounds.T ~ b15km2_100m + (1|Region), data=wui1)
wui_15250 <- clmm(n.compounds.T ~ b15km2_250m + (1|Region), data=wui1)
wui_15500 <- clmm(n.compounds.T ~ b15km2_500m + (1|Region), data=wui1)

wui_30100 <- clmm(n.compounds.T ~ b30km2_100m + (1|Region), data=wui1)
wui_30250 <- clmm(n.compounds.T ~ b30km2_250m + (1|Region), data=wui1)
wui_30500 <- clmm(n.compounds.T ~ b30km2_500m + (1|Region), data=wui1)

wui_60100 <- clmm(n.compounds.T ~ b60km2_100m + (1|Region), data=wui1)
wui_60250 <- clmm(n.compounds.T ~ b60km2_250m + (1|Region), data=wui1)
wui_60500 <- clmm(n.compounds.T ~ b60km2_500m + (1|Region), data=wui1)

wui_sel <- model.sel(wui_45100, wui_15100, wui_30100, wui_60100,
                      wui_45250, wui_15250, wui_30250, wui_60250,
                      wui_45500, wui_15500, wui_30500, wui_60500)
wui_sel

#### Set up data to run for each combination of covariates ####
dat1 <- dat[,c(1:15,23:25)]
dat1 <- distinct(dat1)

# Join percent agriculture
pctAG1 <- pctAG1[,c(1:3, 17:20)]
dat1 <- left_join(dat1, pctAG1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[19:22] <- c("pctAG_4_5", "pctAG_15", "pctAG_30", "pctAG_60")

# Join beech basal area
baa1 <- baa1[,c(1:3, 17, 20)]
dat1 <- left_join(dat1, baa1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[23:24] <- c("baa_4", "baa_60")

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 28)]
dat1 <- left_join(dat1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[25] <- c("intermix_60km2_500m") 

# Join interface WUI
interface1 <- interface1[, c(1:3, 20)]
dat1 <- left_join(dat1, interface1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[26] <- c("interface_60km2_100m") 

# Join total WUI
wui1 <- wui1[,c(1:3, 28)]
dat1 <- left_join(dat1, wui1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[27] <- c("wui_60km2_500m") 


## Set up models

for (i in 1:10) {
  
  m1g <- clmm(n.compounds.T ~ baa_4 + (1|Region), data=dat1, na.action="na.fail")
  summary(m1g)
  m2g <- clmm(n.compounds.T ~ Sex*Age + pctAG_30 + baa_4*laggedBMI + (1|Region), data=dat1, na.action="na.fail")
  summary(m2g)
}



baa_60*year +
Sex*Age + pctAG_15 + pctAG_4_5 + pctAG_30 + pctAG_60 + 
  intermix_60km2_500m + interface_60km2_100m + wui_60km2_500m +
  
  
  
  
  