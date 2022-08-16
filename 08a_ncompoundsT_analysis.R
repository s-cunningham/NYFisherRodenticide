library(tidyverse)
library(MuMIn)
library(ordinal)
library(parallel)

set.seed(123)

#### Set up parallel processing ####

# Determine cluster type (mine is a PSOCK)
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"

# Detect number of cores and create cluster (leave a couple out to not overwhelm pc)
nCores <- detectCores() - 2
cl <- makeCluster(nCores, type=clusterType)

clusterEvalQ(cl,library(ordinal))
clusterEvalQ(cl,library(MuMIn))

#### Read in data ####
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

#### Analysis Set-up ####

# Order response
dat$n.compounds.MO <- ordered(dat$n.compounds.MO, levels=c(0,1,2,3))
dat$n.compounds.T <- ordered(dat$n.compounds.T, levels=c(0,1,2,3,4,5))

# Resort columns
dat <- dat[,c(1:7,30,31,10:12,8,9,13:29)]

## Percent AG
pctAG1 <- dat[, c(1:18, 20:22)]
pctAG1 <- distinct(pctAG1)
pctAG1 <- pctAG1 %>% group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

## Scale and center variables
pctAG1[,c(18:26)] <- scale(pctAG1[,c(18:26)])

# Run models
ag15 <- clmm(n.compounds.T ~ totalag_15 + (1|Region) + (1|WMU), data=pctAG1)
ag15sq <- clmm(n.compounds.T ~ totalag_15 + I(totalag_15^2) + (1|Region) + (1|WMU), data=pctAG1)
ag30 <- clmm(n.compounds.T ~ totalag_30 + (1|Region) + (1|WMU), data=pctAG1)
ag30sq <- clmm(n.compounds.T ~ totalag_30 + I(totalag_30^2) +(1|Region) + (1|WMU), data=pctAG1)
ag60 <- clmm(n.compounds.T ~ totalag_60 + (1|Region) + (1|WMU), data=pctAG1)
ag60sq <- clmm(n.compounds.T ~ totalag_60 + I(totalag_60^2) + (1|Region) + (1|WMU), data=pctAG1)

crop15 <- clmm(n.compounds.T ~ crops_15 + (1|Region) + (1|WMU), data=pctAG1)
crop15sq <- clmm(n.compounds.T ~ crops_15 + I(crops_15^2) + (1|Region) + (1|WMU), data=pctAG1)
crop30 <- clmm(n.compounds.T ~ crops_30 + (1|Region) + (1|WMU), data=pctAG1)
crop30sq <- clmm(n.compounds.T ~ crops_30 + I(crops_30^2) +(1|Region) + (1|WMU), data=pctAG1)
crop60 <- clmm(n.compounds.T ~ crops_60 + (1|Region) + (1|WMU), data=pctAG1)
crop60sq <- clmm(n.compounds.T ~ crops_60 + I(crops_60^2) + (1|Region) + (1|WMU), data=pctAG1)

past15 <- clmm(n.compounds.T ~ pasture_15 + (1|Region) + (1|WMU), data=pctAG1)
past15sq <- clmm(n.compounds.T ~ pasture_15 + I(pasture_15^2) + (1|Region) + (1|WMU), data=pctAG1)
past30 <- clmm(n.compounds.T ~ pasture_30 + (1|Region) + (1|WMU), data=pctAG1)
past30sq <- clmm(n.compounds.T ~ pasture_30 + I(pasture_30^2) +(1|Region) + (1|WMU), data=pctAG1)
past60 <- clmm(n.compounds.T ~ pasture_60 + (1|Region) + (1|WMU), data=pctAG1)
past60sq <- clmm(n.compounds.T ~ pasture_60 + I(pasture_60^2) + (1|Region) + (1|WMU), data=pctAG1)


pctAG_sel <- model.sel(ag15, ag30, ag60, ag15sq, ag30sq, ag60sq,
                       crop15, crop30, crop60, crop15sq, crop30sq, crop60sq,
                       past15, past30, past60, past15sq, past30sq, past60sq)
pctAG_sel

#####
## Percent forest
# pctFOR1 <- dat[, c(1:18, 23:26)]
# pctFOR1 <- distinct(pctFOR1)
# pctFOR1 <- pctFOR1 %>% group_by(RegionalID) %>% 
#   pivot_wider(names_from=buffsize, values_from=c(deciduous, evergreen, mixed, totalforest)) %>% as.data.frame()
# 
# ## Scale and center variables
# pctFOR1[,c(18:29)] <- scale(pctFOR1[,c(18:29)])
# 
# # Run models
# for15 <- clmm(n.compounds.T ~ totalforest_15 + (1|Region) + (1|WMU), data=pctFOR1)
# for15sq <- clmm(n.compounds.T ~ totalforest_15 + I(totalforest_15^2) + (1|Region) + (1|WMU), data=pctFOR1)
# for30 <- clmm(n.compounds.T ~ totalforest_30 + (1|Region) + (1|WMU), data=pctFOR1)
# for30sq <- clmm(n.compounds.T ~ totalforest_30 + I(totalforest_30^2) +(1|Region) + (1|WMU), data=pctFOR1)
# for60 <- clmm(n.compounds.T ~ totalforest_60 + (1|Region) + (1|WMU), data=pctFOR1)
# for60sq <- clmm(n.compounds.T ~ totalforest_60 + I(totalforest_60^2) + (1|Region) + (1|WMU), data=pctFOR1)
# 
# decid15 <- clmm(n.compounds.T ~ deciduous_15 + (1|Region) + (1|WMU), data=pctFOR1)
# decid15sq <- clmm(n.compounds.T ~ deciduous_15 + I(deciduous_15^2) + (1|Region) + (1|WMU), data=pctFOR1)
# decid30 <- clmm(n.compounds.T ~ deciduous_30 + (1|Region) + (1|WMU), data=pctFOR1)
# decid30sq <- clmm(n.compounds.T ~ deciduous_30 + I(deciduous_30^2) +(1|Region) + (1|WMU), data=pctFOR1)
# decid60 <- clmm(n.compounds.T ~ deciduous_60 + (1|Region) + (1|WMU), data=pctFOR1)
# decid60sq <- clmm(n.compounds.T ~ deciduous_60 + I(deciduous_60^2) + (1|Region) + (1|WMU), data=pctFOR1)
# 
# ever15 <- clmm(n.compounds.T ~ evergreen_15 + (1|Region) + (1|WMU), data=pctFOR1)
# ever15sq <- clmm(n.compounds.T ~ evergreen_15 + I(evergreen_15^2) + (1|Region) + (1|WMU), data=pctFOR1)
# ever30 <- clmm(n.compounds.T ~ evergreen_30 + (1|Region) + (1|WMU), data=pctFOR1)
# ever30sq <- clmm(n.compounds.T ~ evergreen_30 + I(evergreen_30^2) +(1|Region) + (1|WMU), data=pctFOR1)
# ever60 <- clmm(n.compounds.T ~ evergreen_60 + (1|Region) + (1|WMU), data=pctFOR1)
# ever60sq <- clmm(n.compounds.T ~ evergreen_60 + I(evergreen_60^2) + (1|Region) + (1|WMU), data=pctFOR1)
# 
# mfor15 <- clmm(n.compounds.T ~ mixed_15 + (1|Region) + (1|WMU), data=pctFOR1)
# mfor15sq <- clmm(n.compounds.T ~ mixed_15 + I(mixed_15^2) + (1|Region) + (1|WMU), data=pctFOR1)
# mfor30 <- clmm(n.compounds.T ~ mixed_30 + (1|Region) + (1|WMU), data=pctFOR1)
# mfor30sq <- clmm(n.compounds.T ~ mixed_30 + I(mixed_30^2) +(1|Region) + (1|WMU), data=pctFOR1)
# mfor60 <- clmm(n.compounds.T ~ mixed_60 + (1|Region) + (1|WMU), data=pctFOR1)
# mfor60sq <- clmm(n.compounds.T ~ mixed_60 + I(mixed_60^2) + (1|Region) + (1|WMU), data=pctFOR1)
# 
# pctFOR_sel <- model.sel(for15, for30, for60, for15sq, for30sq, for60sq,
#                        decid15, decid30, decid60, decid15sq, decid30sq, decid60sq,
#                        ever15, ever30, ever60, ever15sq, ever30sq, ever60sq,
#                        mfor15, mfor30, mfor60, mfor15sq, mfor30sq, mfor60sq)
# pctFOR_sel
#####

## Beech basal area
baa1 <- dat[, c(1:18, 27)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=baa, values_fn=unique) %>% as.data.frame()
names(baa1)[18:20] <- c("baa_15", "baa_30", "baa_60") 

## Scale and center variables
baa1[,c(18:20)] <- scale(baa1[,c(18:20)])

baa15 <- clmm(n.compounds.T ~ baa_15 + (1|Region) + (1|WMU), data=baa1)
baa15sq <- clmm(n.compounds.T ~ baa_15 + I(baa_15^2) + (1|Region) + (1|WMU), data=baa1)
baa30 <- clmm(n.compounds.T ~ baa_30 + (1|Region) + (1|WMU), data=baa1)
baa30sq <- clmm(n.compounds.T ~ baa_30 + I(baa_30^2) + (1|Region) + (1|WMU), data=baa1)
baa60 <- clmm(n.compounds.T ~ baa_60 + (1|Region) + (1|WMU), data=baa1)
baa60sq <- clmm(n.compounds.T ~ baa_60 + I(baa_60^2) + (1|Region) + (1|WMU), data=baa1)

baa_sel <- model.sel(baa15, baa30, baa60, baa15sq, baa30sq, baa60sq)
baa_sel

## Intermix WUI
intermix1 <- dat[, c(1:19, 28)]
intermix1 <- unite(intermix1, "buffrad", 18:19, sep="_")
intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[18:26] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                             "mix_15_250", "mix_30_250", "mix_60_250",
                             "mix_15_500", "mix_30_500", "mix_60_500") 

## Scale and center variables
intermix1[,c(18:26)] <- scale(intermix1[,c(18:26)])

mix_15100 <- clmm(n.compounds.T ~ mix_15_100 + (1|Region) + (1|WMU), data=intermix1)
mix_15250 <- clmm(n.compounds.T ~ mix_15_250 + (1|Region) + (1|WMU), data=intermix1)
mix_15500 <- clmm(n.compounds.T ~ mix_15_500 + (1|Region) + (1|WMU), data=intermix1)
mix_15100sq <- clmm(n.compounds.T ~ mix_15_100 + I(mix_15_100^2) + (1|Region) + (1|WMU), data=intermix1)
mix_15250sq <- clmm(n.compounds.T ~ mix_15_250 + I(mix_15_250^2) + (1|Region) + (1|WMU), data=intermix1)
mix_15500sq <- clmm(n.compounds.T ~ mix_15_500 + I(mix_15_500^2) + (1|Region) + (1|WMU), data=intermix1)

mix_30100 <- clmm(n.compounds.T ~ mix_30_100 + (1|Region) + (1|WMU), data=intermix1)
mix_30250 <- clmm(n.compounds.T ~ mix_30_250 + (1|Region) + (1|WMU), data=intermix1)
mix_30500 <- clmm(n.compounds.T ~ mix_30_500 + (1|Region) + (1|WMU), data=intermix1)
mix_30100sq <- clmm(n.compounds.T ~ mix_30_100 + I(mix_30_100^2) + (1|Region) + (1|WMU), data=intermix1)
mix_30250sq <- clmm(n.compounds.T ~ mix_30_250 + I(mix_30_250^2) + (1|Region) + (1|WMU), data=intermix1)
mix_30500sq <- clmm(n.compounds.T ~ mix_30_500 + I(mix_30_500^2) + (1|Region) + (1|WMU), data=intermix1)

mix_60100 <- clmm(n.compounds.T ~ mix_60_100 + (1|Region) + (1|WMU), data=intermix1)
mix_60250 <- clmm(n.compounds.T ~ mix_60_250 + (1|Region) + (1|WMU), data=intermix1)
mix_60500 <- clmm(n.compounds.T ~ mix_60_500 + (1|Region) + (1|WMU), data=intermix1)
mix_60100sq <- clmm(n.compounds.T ~ mix_60_100 + I(mix_60_100^2) + (1|Region) + (1|WMU), data=intermix1)
mix_60250sq <- clmm(n.compounds.T ~ mix_60_250 + I(mix_60_250^2) + (1|Region) + (1|WMU), data=intermix1)
mix_60500sq <- clmm(n.compounds.T ~ mix_60_500 + I(mix_60_500^2) + (1|Region) + (1|WMU), data=intermix1)

intermix_sel <- model.sel(mix_15100, mix_30100, mix_60100, mix_15100sq, mix_30100sq, mix_60100sq,
                          mix_15250, mix_30250, mix_60250, mix_15250sq, mix_30250sq, mix_60250sq,
                          mix_15500, mix_30500, mix_60500, mix_15500sq, mix_30500sq, mix_60500sq)
intermix_sel

## Interface WUI
interface1 <- dat[, c(1:19, 29)]
interface1 <- unite(interface1, "buffrad", 18:19, sep="_")
interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[18:26] <- c("face_15_100", "face_30_100", "face_60_100",
                              "face_15_250", "face_30_250", "face_60_250",
                              "face_15_500", "face_30_500", "face_60_500") 

## Scale and center variables
interface1[,c(18:26)] <- scale(interface1[,c(18:26)])

face_15100 <- clmm(n.compounds.T ~ face_15_100 + (1|Region) + (1|WMU), data=interface1)
face_15250 <- clmm(n.compounds.T ~ face_15_250 + (1|Region) + (1|WMU), data=interface1)
face_15500 <- clmm(n.compounds.T ~ face_15_500 + (1|Region) + (1|WMU), data=interface1)
face_15100sq <- clmm(n.compounds.T ~ face_15_100 + I(face_15_100^2) + (1|Region) + (1|WMU), data=interface1)
face_15250sq <- clmm(n.compounds.T ~ face_15_250 + I(face_15_250^2) + (1|Region) + (1|WMU), data=interface1)
face_15500sq <- clmm(n.compounds.T ~ face_15_500 + I(face_15_500^2) + (1|Region) + (1|WMU), data=interface1)

face_30100 <- clmm(n.compounds.T ~ face_30_100 + (1|Region) + (1|WMU), data=interface1)
face_30250 <- clmm(n.compounds.T ~ face_30_250 + (1|Region) + (1|WMU), data=interface1)
face_30500 <- clmm(n.compounds.T ~ face_30_500 + (1|Region) + (1|WMU), data=interface1)
face_30100sq <- clmm(n.compounds.T ~ face_30_100 + I(face_30_100^2) + (1|Region) + (1|WMU), data=interface1)
face_30250sq <- clmm(n.compounds.T ~ face_30_250 + I(face_30_250^2) + (1|Region) + (1|WMU), data=interface1)
face_30500sq <- clmm(n.compounds.T ~ face_30_500 + I(face_30_500^2) + (1|Region) + (1|WMU), data=interface1)

face_60100 <- clmm(n.compounds.T ~ face_60_100 + (1|Region) + (1|WMU), data=interface1)
face_60250 <- clmm(n.compounds.T ~ face_60_250 + (1|Region) + (1|WMU), data=interface1)
face_60500 <- clmm(n.compounds.T ~ face_60_500 + (1|Region) + (1|WMU), data=interface1)
face_60100sq <- clmm(n.compounds.T ~ face_60_100 + I(face_60_100^2) + (1|Region) + (1|WMU), data=interface1)
face_60250sq <- clmm(n.compounds.T ~ face_60_250 + I(face_60_250^2) + (1|Region) + (1|WMU), data=interface1)
face_60500sq <- clmm(n.compounds.T ~ face_60_500 + I(face_60_500^2) + (1|Region) + (1|WMU), data=interface1)

interface_sel <- model.sel(face_15100, face_30100, face_60100, face_15100sq, face_30100sq, face_60100sq,  
                           face_15250, face_30250, face_60250, face_15250sq, face_30250sq, face_60250sq,  
                           face_15500, face_30500, face_60500, face_15500sq, face_30500sq, face_60500sq)  
interface_sel

## Total WUI
wui1 <- dat[, c(1:19, 30)]
wui1 <- unite(wui1, "buffrad", 18:19, sep="_")
wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[18:26] <- c("wui_15_100", "wui_30k_100", "wui_60_100",
                        "wui_15_250", "wui_30k_250", "wui_60_250",
                        "wui_15_500", "wui_30k_500", "wui_60_500") 

## Scale and center variables
wui1[,c(18:26)] <- scale(wui1[,c(18:26)])

wui_15100 <- clmm(n.compounds.T ~ wui_15_100 + (1|Region) + (1|WMU), data=wui1)
wui_15250 <- clmm(n.compounds.T ~ wui_15_250 + (1|Region) + (1|WMU), data=wui1)
wui_15500 <- clmm(n.compounds.T ~ wui_15_500 + (1|Region) + (1|WMU), data=wui1)
wui_15100sq <- clmm(n.compounds.T ~ wui_15_100 + I(wui_15_100^2) + (1|Region) + (1|WMU), data=wui1)
wui_15250sq <- clmm(n.compounds.T ~ wui_15_250 + I(wui_15_250^2) + (1|Region) + (1|WMU), data=wui1)
wui_15500sq <- clmm(n.compounds.T ~ wui_15_500 + I(wui_15_500^2) + (1|Region) + (1|WMU), data=wui1)

wui_30100 <- clmm(n.compounds.T ~ wui_30k_100 + (1|Region) + (1|WMU), data=wui1)
wui_30250 <- clmm(n.compounds.T ~ wui_30k_250 + (1|Region) + (1|WMU), data=wui1)
wui_30500 <- clmm(n.compounds.T ~ wui_30k_500 + (1|Region) + (1|WMU), data=wui1)
wui_30100sq <- clmm(n.compounds.T ~ wui_30k_100 + I(wui_30k_100^2) + (1|Region) + (1|WMU), data=wui1)
wui_30250sq <- clmm(n.compounds.T ~ wui_30k_250 + I(wui_30k_250^2) + (1|Region) + (1|WMU), data=wui1)
wui_30500sq <- clmm(n.compounds.T ~ wui_30k_500 + I(wui_30k_500^2) + (1|Region) + (1|WMU), data=wui1)

wui_60100 <- clmm(n.compounds.T ~ wui_60_100 + (1|Region) + (1|WMU), data=wui1)
wui_60250 <- clmm(n.compounds.T ~ wui_60_250 + (1|Region) + (1|WMU), data=wui1)
wui_60500 <- clmm(n.compounds.T ~ wui_60_500 + (1|Region) + (1|WMU), data=wui1)
wui_60100sq <- clmm(n.compounds.T ~ wui_60_100 + I(wui_60_100^2) + (1|Region) + (1|WMU), data=wui1)
wui_60250sq <- clmm(n.compounds.T ~ wui_60_250 + I(wui_60_250^2) + (1|Region) + (1|WMU), data=wui1)
wui_60500sq <- clmm(n.compounds.T ~ wui_60_500 + I(wui_60_500^2) + (1|Region) + (1|WMU), data=wui1)

wui_sel <- model.sel(wui_15100, wui_30100, wui_60100, wui_15100sq, wui_30100sq, wui_60100sq, 
                     wui_15250, wui_30250, wui_60250, wui_15250sq, wui_30250sq, wui_60250sq, 
                     wui_15500, wui_30500, wui_60500, wui_15500sq, wui_30500sq, wui_60500sq)
wui_sel

#### Set up data to run for each combination of covariates ####
dat1 <- dat[,c(1:17,31)]
dat1 <- distinct(dat1)

# Join percent agriculture
pctAG1 <- pctAG1[,c(1:3, 20, 26)]
dat1 <- left_join(dat1, pctAG1, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 20)]
dat1 <- left_join(dat1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 26)]
dat1 <- left_join(dat1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

# Join interface WUI
interface1 <- interface1[, c(1:3, 20)]
dat1 <- left_join(dat1, interface1, by=c("RegionalID", "pt_name", "pt_index"))

# Join total WUI
wui1 <- wui1[,c(1:3, 23, 25, 26)]
dat1 <- left_join(dat1, wui1, by=c("RegionalID", "pt_name", "pt_index"))

# join forest
# pctFOR1 <- pctFOR1[,c(1:3, 24, 29)]
# dat1 <- left_join(dat1, pctFOR1, by=c("RegionalID", "pt_name", "pt_index"))

## Scale age and laggedBMI
dat1$Age <- scale(dat1$Age)
dat1$laggedBMI <- scale(dat1$laggedBMI)

## Check correlation matrix
cor(dat1[,18:26])

## Different way of accounting for year/mast cycle
dat1$fBMI <- ordered(dat1$year, levels=c(2019, 2018, 2020))
dat1$year <- as.factor(dat1$year)

## Set up global models
# Human-driven hypothesis
g1 <- clmm(n.compounds.T ~ Sex*Age + 
                pasture_60 + I(pasture_60^2) + 
                totalag_60 + I(totalag_60^2) +
                mix_60_500 + I(mix_60_500 ^2) +
                face_60_100 + I(face_60_100^2) +
                wui_60_250 + I(wui_60_250^2) +
                wui_30k_500 + I(wui_30k_500^2) +
                wui_60_500 +I(wui_60_500^2) +
                baa_60 + laggedBMI + I(baa_60^2) +
                baa_60:laggedBMI + I(baa_60^2):laggedBMI +
                (1|Region) + (1|WMU), data=dat1, na.action="na.fail")


# Export data and model into the cluster worker nodes
clusterExport(cl, c("dat1","g1"))

g1_dredge <- MuMIn:::.dredge.par(g1, cluster=cl, trace=2, subset=!(pasture_60 && totalag_60) &&
                                                                          !(wui_60_250 && wui_30k_500) &&
                                                                          !(wui_60_250 && wui_60_500) &&
                                                                          !(wui_60_250 && mix_60_500) &&
                                                                          !(wui_60_250 && face_60_100) &&
                                                                          !(mix_60_500 && face_60_100) && 
                                                                          !(wui_30k_500 && face_60_100) &&
                                                                          !(wui_60_500 && face_60_100) &&
                                                                          !(wui_60_500 && mix_60_500) &&
                                                                          !(wui_30k_500 && mix_60_500) &&
                                                                          dc(face_60_100, I(face_60_100^2)) &&
                                                                          dc(mix_60_500, I(mix_60_500^2)) &&
                                                                          dc(wui_60_250, I(wui_60_250^2)) &&
                                                                          dc(wui_60_500, I(wui_60_500^2)) &&
                                                                          dc(wui_30k_500, I(wui_30k_500^2)) &&
                                                                          dc(pasture_60, I(pasture_60^2)) &&
                                                                          dc(totalag_60, I(totalag_60^2)) &&
                                                                          dc(baa_60:laggedBMI, I(baa_60^2):laggedBMI) &&
                                                                          dc(baa_60, I(baa_60^2), I(baa_60^2):laggedBMI))

# Save dredge tables
save(file="output/dredge_tables_humansT.Rdata", list="g1_dredge")

dh <- as.data.frame(g1_dredge)
write_csv(dh, "output/human-models_ncompT.csv")

## Read table back in...hopefully simple data frame
dh <- read_csv("output/human-models_ncompT.csv")
dh <- as.data.frame(dh)

# Remove rows where the quadratic total ag was included without the linear
dh <- dh[!(is.na(dh$totalag_60) & dh$`I(totalag_60^2)`>0),]
dh <- dh[!is.na(dh$AICc),]

# Remove rows where quadratic interaction was included but not the linear interaction
dh <- dh[!(is.na(dh$`baa_60:laggedBMI`) & dh$`I(baa_60^2):laggedBMI`>0),]
dh <- dh[!is.na(dh$AICc),]

# Clear deltaAIC and model weights
dh$delta <- NA
dh$weight <- NA

# Recalculate deltaAIC and model weights
dh$delta <- dh$AICc - dh$AICc[1]
w <- qpcR::akaike.weights(dh$AICc)
dh$weight <- w$weights

#### Running final models ####

# See summary for top models (deltaAICc < 2)
m1 <- clmm(n.compounds.T ~ Sex*Age + baa_60 + laggedBMI + I(baa_60^2) +
             baa_60:laggedBMI + I(baa_60^2):laggedBMI +
             mix_60_500 + I(mix_60_500^2) +
             totalag_60 + (1|Region) + (1|WMU), data=dat1)

m2 <- clmm(n.compounds.T ~ Sex*Age + baa_60 + laggedBMI + I(baa_60^2) +
             baa_60:laggedBMI + I(baa_60^2):laggedBMI +
             mix_60_500 + I(mix_60_500^2) +
             totalag_60 + I(totalag_60^2) +
             (1|Region) + (1|WMU), data=dat1)

## Loop over each set of random points

m_est <- data.frame()
pct2.5 <- data.frame()
pct97.5 <- data.frame()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- dat1[dat1$pt_index==i,]
  
  # Run both models with deltaAICc < 2
  m1_pt <- clmm(n.compounds.T ~ Sex*Age + baa_60 + laggedBMI + I(baa_60^2) +
               baa_60:laggedBMI + I(baa_60^2):laggedBMI + mix_60_500 + I(mix_60_500^2) +
               totalag_60 + (1|Region), data=pt)
  
  m2_pt <- clmm(n.compounds.T ~ Sex*Age + baa_60 + laggedBMI + I(baa_60^2) +
               baa_60:laggedBMI + I(baa_60^2):laggedBMI +
               mix_60_500 + I(mix_60_500^2) + totalag_60 + I(totalag_60^2) +
               (1|Region), data=pt)
  
  # Create model selection object for averaging
  mod_pt <- model.sel(m1_pt, m2_pt)
  
  # Model average estimates
  avgd <- model.avg(mod_pt)
  avg_smry <- summary(avgd)
  
  # save averaged confidence intervals
  pct2.5 <- rbind(pct2.5, t(confint(avgd, full=TRUE))[1,])
  pct97.5 <- rbind(pct97.5, t(confint(avgd, full=TRUE))[2,])

  # Save point set estimates
  m_est <- rbind(m_est, avgd$coefficients[row.names(avgd$coefficients)=="full"])
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- c(colnames(avgd$coefficients))
    names(pct2.5) <- c(colnames(avgd$coefficients))
    names(pct97.5) <- c(colnames(avgd$coefficients))
    }
  
}

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)])
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)])
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)])


# Combine and clean up data frame
coef_summary <- bind_rows(coef_avg, pct2.5_avg, pct97.5_avg)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("coef_mean", "2.5CI", "97.5CI")
coef_summary <- data.frame(coef=coefs, coef_summary)

# Write to file
write_csv(coef_summary, "results/ncompT_coef-summary.csv")


















