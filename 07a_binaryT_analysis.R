library(tidyverse)
library(lme4)
library(parallel)
library(MuMIn)

set.seed(123)
#### Set up parallel processing ####

# Determine cluster type (mine is a PSOCK)
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"

# Detect number of cores and create cluster (leave a couple out to not overwhelm pc)
nCores <- detectCores() - 4
cl <- makeCluster(nCores, type=clusterType)

clusterEvalQ(cl,library(lme4))
clusterEvalQ(cl,library(MuMIn))

#### Read in data ####
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

#### Analysis ####

## set up binary variable
dat$binary.MO <- ifelse(dat$n.compounds.MO==0, 0, 1)
dat$binary.T <- ifelse(dat$n.compounds.T==0, 0, 1)

## Reorder columns
dat <- dat[,c(1:7, 31, 8:13, 32, 33, 16:29)]

## Scale and center variables
dat[,c(10,19:30)] <- scale(dat[,c(10,19:30)])

# Change how beech mast is incorporated
dat$mast <- NA 
dat$mast[dat$year==2018] <- "intermediate"
dat$mast[dat$year==2019] <- "fail"
dat$mast[dat$year==2020] <- "high"
dat$mast <-  ordered(dat$mast, levels=c("fail", "intermediate", "high"))

dat$fyear <- factor(dat$year)

## Use pooled data to determine scale

## Percent AG
pctAG1 <- dat[, c(1:3, 6:8, 16, 17, 19:21)]
pctAG1 <- distinct(pctAG1)
pctAG1 <- pctAG1 %>% group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

ag15 <- glmer(binary.T ~ totalag_15 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
ag30 <- glmer(binary.T ~ totalag_30 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
ag60 <- glmer(binary.T ~ totalag_60 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
ag15sq <- glmer(binary.T ~ totalag_15 + I(totalag_15^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
ag30sq <- glmer(binary.T ~ totalag_30 + I(totalag_30^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
ag60sq <- glmer(binary.T ~ totalag_60 + I(totalag_60^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)

crops15 <- glmer(binary.T ~ crops_15 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
crops30 <- glmer(binary.T ~ crops_30 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
crops60 <- glmer(binary.T ~ crops_60 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
crops15sq <- glmer(binary.T ~ crops_15 + I(crops_15^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
crops30sq <- glmer(binary.T ~ crops_30 + I(crops_30^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
crops60sq <- glmer(binary.T ~ crops_60 + I(crops_60^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)

past15 <- glmer(binary.T ~ pasture_15 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
past30 <- glmer(binary.T ~ pasture_30 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
past60 <- glmer(binary.T ~ pasture_60 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
past15sq <- glmer(binary.T ~ pasture_15 + I(pasture_15^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
past30sq <- glmer(binary.T ~ pasture_30 + I(pasture_30^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)
past60sq <- glmer(binary.T ~ pasture_60 + I(pasture_60^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=pctAG1)

pctAG_sel <- model.sel(ag15, ag15sq, crops15, crops15sq, past15, past15sq,
                       ag30, ag30sq, crops30, crops30sq, past30, past30sq,
                       ag60, ag60sq, crops60, crops60sq, past60, past60sq)
pctAG_sel

## Beech basal area
baa1 <- dat[, c(1:3, 6:8, 16, 17, 26, 30:32)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=baa, values_fn=unique) %>% as.data.frame()
names(baa1)[11:13] <- c("baa_15", "baa_30", "baa_60") 

baa15 <- glmer(binary.T ~ baa_15 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa30 <- glmer(binary.T ~ baa_30 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa60 <- glmer(binary.T ~ baa_60 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa15sq <- glmer(binary.T ~ baa_15 + I(baa_15^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa30sq <- glmer(binary.T ~ baa_30 + I(baa_30^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa60sq <- glmer(binary.T ~ baa_60 + I(baa_60^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)

baa15lBMI <- glmer(binary.T ~ baa_15*laggedBMI + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa15sqlBMI <- glmer(binary.T ~ baa_15*laggedBMI + I(baa_15^2)*laggedBMI + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa30lBMI <- glmer(binary.T ~ baa_30*laggedBMI + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa30sqlBMI <- glmer(binary.T ~ baa_30*laggedBMI + I(baa_30^2)*laggedBMI + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa60lBMI <- glmer(binary.T ~ baa_60*laggedBMI + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa60sqlBMI <- glmer(binary.T ~ baa_60*laggedBMI + I(baa_60^2)*laggedBMI + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)

baa15M <- glmer(binary.T ~ baa_15*mast + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa15sqM <- glmer(binary.T ~ baa_15*mast + I(baa_15^2)*mast + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa30M <- glmer(binary.T ~ baa_30*mast + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa30sqM <- glmer(binary.T ~ baa_30*mast + I(baa_30^2)*mast + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa60M <- glmer(binary.T ~ baa_60*mast + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa60sqM <- glmer(binary.T ~ baa_60*mast + I(baa_60^2)*mast + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)

baa15y <- glmer(binary.T ~ baa_15*fyear + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa15sqy <- glmer(binary.T ~ baa_15*fyear + I(baa_15^2)*fyear + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa30y <- glmer(binary.T ~ baa_30*fyear + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa30sqy <- glmer(binary.T ~ baa_30*fyear + I(baa_30^2)*fyear + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa60y <- glmer(binary.T ~ baa_60*fyear + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)
baa60sqy <- glmer(binary.T ~ baa_60*fyear + I(baa_60^2)*fyear + (1|WMUA_code/WMU), family=binomial(link="logit"), data=baa1)

baa_sel <- model.sel(baa15, baa30, baa60, baa15sq, baa30sq, baa60sq, 
                     baa15lBMI, baa15sqlBMI, baa30lBMI, baa30sqlBMI, baa60lBMI, baa60sqlBMI,
                     baa15M, baa15sqM, baa30M, baa30sqM, baa60M, baa60sqM,
                     baa15y, baa15sqy, baa30y, baa30sqy, baa60y, baa60sqy)
baa_sel

## Wildland-urban interface
intermix1 <- dat[, c(1:3, 6:8, 16:18, 27)]
intermix1 <- unite(intermix1, "buffrad", 8:9, sep="_")
intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[8:16] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                            "mix_15_250", "mix_30_250", "mix_60_250",
                            "mix_15_500", "mix_30_500", "mix_60_500")  


interface1 <- dat[, c(1:3, 6:8, 16:18, 28)]
interface1 <- unite(interface1, "buffrad", 8:9, sep="_")
interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[8:16] <- c("face_15_100", "face_30_100", "face_60_100",
                             "face_15_250", "face_30_250", "face_60_250",
                             "face_15_500", "face_30_500", "face_60_500") 

wui1 <- dat[, c(1:3, 6:8, 16:18, 29)]
wui1 <- unite(wui1, "buffrad", 8:9, sep="_")
wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[8:16] <- c("wui_15_100", "wui_30k_100", "wui_60_100",
                       "wui_15_250", "wui_30k_250", "wui_60_250",
                       "wui_15_500", "wui_30k_500", "wui_60_500")

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 8:16)]
wui1 <- left_join(wui1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

# Join interface WUI
interface1 <- interface1[, c(1:3, 8:16)]
wui1 <- left_join(wui1, interface1, by=c("RegionalID", "pt_name", "pt_index"))

# Intermix
mix_15100 <- glmer(binary.T ~ mix_15_100 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_15250 <- glmer(binary.T ~ mix_15_250 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_15500 <- glmer(binary.T ~ mix_15_500 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_15100sq <- glmer(binary.T ~ mix_15_100 + I(mix_15_100^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_15250sq <- glmer(binary.T ~ mix_15_250 + I(mix_15_250^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_15500sq <- glmer(binary.T ~ mix_15_500 + I(mix_15_500^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)

mix_30100 <- glmer(binary.T  ~ mix_30_100 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_30250 <- glmer(binary.T  ~ mix_30_250 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_30500 <- glmer(binary.T  ~ mix_30_500 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_30100sq <- glmer(binary.T ~ mix_30_100 + I(mix_30_100^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_30250sq <- glmer(binary.T ~ mix_30_250 + I(mix_30_250^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_30500sq <- glmer(binary.T ~ mix_30_500 + I(mix_30_500^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)

mix_60100 <- glmer(binary.T ~ mix_60_100 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_60250 <- glmer(binary.T ~ mix_60_250 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_60500 <- glmer(binary.T ~ mix_60_500 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_60100sq <- glmer(binary.T ~ mix_60_100 + I(mix_60_100^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_60250sq <- glmer(binary.T ~ mix_60_250 + I(mix_60_250^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
mix_60500sq <- glmer(binary.T ~ mix_60_500 + I(mix_60_500^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)

# Interface
face_15100 <- glmer(binary.T ~ face_15_100 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_15250 <- glmer(binary.T ~ face_15_250 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_15500 <- glmer(binary.T ~ face_15_500 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_15100sq <- glmer(binary.T ~ face_15_100 + I(face_15_100^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_15250sq <- glmer(binary.T ~ face_15_250 + I(face_15_250^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_15500sq <- glmer(binary.T ~ face_15_500 + I(face_15_500^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)

face_30100 <- glmer(binary.T  ~ face_30_100 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_30250 <- glmer(binary.T  ~ face_30_250 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_30500 <- glmer(binary.T  ~ face_30_500 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_30100sq <- glmer(binary.T ~ face_30_100 + I(face_30_100^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_30250sq <- glmer(binary.T ~ face_30_250 + I(face_30_250^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_30500sq <- glmer(binary.T ~ face_30_500 + I(face_30_500^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)

face_60100 <- glmer(binary.T ~ face_60_100 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_60250 <- glmer(binary.T ~ face_60_250 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_60500 <- glmer(binary.T ~ face_60_500 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_60100sq <- glmer(binary.T ~ face_60_100 + I(face_60_100^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_60250sq <- glmer(binary.T ~ face_60_250 + I(face_60_250^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
face_60500sq <- glmer(binary.T ~ face_60_500 + I(face_60_500^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)

# Total WUI
wui_15100 <- glmer(binary.T ~ wui_15_100 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_15250 <- glmer(binary.T ~ wui_15_250 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_15500 <- glmer(binary.T ~ wui_15_500 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_15100sq <- glmer(binary.T ~ wui_15_100 + I(wui_15_100^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_15250sq <- glmer(binary.T ~ wui_15_250 + I(wui_15_250^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_15500sq <- glmer(binary.T ~ wui_15_500 + I(wui_15_500^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)

wui_30100 <- glmer(binary.T ~ wui_30k_100 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_30250 <- glmer(binary.T ~ wui_30k_250 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_30500 <- glmer(binary.T ~ wui_30k_500 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_30100sq <- glmer(binary.T ~ wui_30k_100 + I(wui_30k_100^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_30250sq <- glmer(binary.T ~ wui_30k_250 + I(wui_30k_250^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_30500sq <- glmer(binary.T ~ wui_30k_500 + I(wui_30k_500^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)

wui_60100 <- glmer(binary.T ~ wui_60_100 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_60250 <- glmer(binary.T ~ wui_60_250 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_60500 <- glmer(binary.T ~ wui_60_500 + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_60100sq <- glmer(binary.T ~ wui_60_100 + I(wui_60_100^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_60250sq <- glmer(binary.T ~ wui_60_250 + I(wui_60_250^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)
wui_60500sq <- glmer(binary.T ~ wui_60_500 + I(wui_60_500^2) + (1|WMUA_code/WMU), family=binomial(link="logit"), data=wui1)

wui_sel <- model.sel(wui_15100, wui_30100, wui_60100, wui_15100sq, wui_30100sq, wui_60100sq, 
                     wui_15250, wui_30250, wui_60250, wui_15250sq, wui_30250sq, wui_60250sq, 
                     wui_15500, wui_30500, wui_60500, wui_15500sq, wui_30500sq, wui_60500sq,
                     mix_15100, mix_30100, mix_60100, mix_15100sq, mix_30100sq, mix_60100sq,
                     mix_15250, mix_30250, mix_60250, mix_15250sq, mix_30250sq, mix_60250sq,
                     mix_15500, mix_30500, mix_60500, mix_15500sq, mix_30500sq, mix_60500sq,
                     face_15100, face_30100, face_60100, face_15100sq, face_30100sq, face_60100sq,
                     face_15250, face_30250, face_60250, face_15250sq, face_30250sq, face_60250sq,
                     face_15500, face_30500, face_60500, face_15500sq, face_30500sq, face_60500sq)
wui_sel

#### Set up data to run for each combination of covariates ####
dat1 <- dat[,c(1:14,16,31,32)]
dat1 <- distinct(dat1)

# Join percent agriculture
pctAG1 <- pctAG1[,c(1:3, 10, 16)]
dat1 <- left_join(dat1, pctAG1, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 13)]
dat1 <- left_join(dat1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join WUI
wui1 <- wui1[,c(1:3, 25)]
dat1 <- left_join(dat1, wui1, by=c("RegionalID", "pt_name", "pt_index"))

#### Run models ####

## Check correlation matrix
cor(dat1[,18:21])

## Set up global model
g1 <- glmer(binary.T ~ Sex*Age + mast + 
              mix_60_500 + I(mix_60_500^2) +
              totalag_60 + I(totalag_60^2) +
              pasture_60 + I(pasture_60^2) +
              baa_60 + I(baa_60^2) +
              baa_60:mast + I(baa_60^2):mast + 
              (1|WMUA_code/WMU), family=binomial(link="logit"), data=dat1, na.action="na.fail")

g2 <- glmer(binary.T ~ Sex*Age + fyear + 
              mix_60_500 + I(mix_60_500^2) +
              totalag_60 + I(totalag_60^2) +
              pasture_60 + I(pasture_60^2) +
              baa_60 + I(baa_60^2) +
              baa_60:fyear + I(baa_60^2):fyear + 
              (1|WMUA_code/WMU), family=binomial(link="logit"), data=dat1, na.action="na.fail")

# Export data and model into the cluster worker nodes
clusterExport(cl, c("dat1","g1","g2"))

## Build model sets
g1_dredge <- MuMIn:::.dredge.par(g1, cluster=cl, trace=2,  subset=!(pasture_60 && totalag_60) &&
                                                                    dc(mix_60_500, I(mix_60_500^2)) &&
                                                                    dc(pasture_60, I(pasture_60^2)) &&
                                                                    dc(totalag_60, I(totalag_60^2)) &&
                                                                    dc(baa_60:mast, I(baa_60^2):mast) &&
                                                                    dc(baa_60, I(baa_60^2), I(baa_60^2):mast))

g2_dredge <- MuMIn:::.dredge.par(g2, cluster=cl, trace=2,  subset=!(pasture_60 && totalag_60) &&
                                                                    dc(mix_60_500, I(mix_60_500^2)) &&
                                                                    dc(pasture_60, I(pasture_60^2)) &&
                                                                    dc(totalag_60, I(totalag_60^2)) &&
                                                                    dc(baa_60:fyear, I(baa_60^2):fyear) &&
                                                                    dc(baa_60, I(baa_60^2), I(baa_60^2):fyear))


# Save dredge tables
save(file="output/dredge_tables_binaryT.Rdata", list="g1_dredge")
save(file="output/dredge_tables_binaryT2.Rdata", list="g2_dredge")

# Save as csv file
dh <- as.data.frame(g1_dredge)
write_csv(dh, "output/models_binaryT.csv")
dh2 <- as.data.frame(g2_dredge)
write_csv(dh2, "output/models_binaryT2.csv")

#### Read tables back in
dh <- read_csv("output/models_binaryT.csv")
dh <- as.data.frame(dh)

dh2 <- read_csv("output/models_binaryT2.csv")
dh2 <- as.data.frame(dh2)

# Combine and reorder
dh <- bind_rows(dh, dh2)

# Reorder columns
dh <- dh[,c(1:15, 21:23, 16:20)]

# Clear deltaAIC and model weights
dh$delta <- NA
dh$weight <- NA

# Reorder according to AICc
N.order <- order(dh$AICc)
dh <- dh[N.order,]

# Recalculate deltaAIC and model weights
dh$delta <- dh$AICc - dh$AICc[1]
w <- qpcR::akaike.weights(dh$AICc)
dh$weight <- w$weights

#### Running final models ####

# See summary for top models (deltaAICc < 2)
m1 <- glmer(binary.T ~ Sex*Age + baa_60 + laggedBMI + baa_60:laggedBMI +
             I(baa_60^2) + totalag_60 + wui_60_100 + (1|Region), data=dat1)

m2 <- glmer(binary.T ~ Sex*Age + baa_30 + I(baa_30^2) + 
             wui_60_100 + I(wui_60_100^2) +
             crops_60 + I(crops_60^2) + (1|Region), data=dat1)

m3 <- glmer(binary.T ~ Sex*Age + baa_60 + I(baa_60^2) + 
             wui_60_100 + I(wui_60_100^2) +
             crops_60 + I(crops_60^2) + (1|Region), data=dat1)

m4 <- glmer(binary.T ~ Sex*Age + baa_60 + I(baa_60^2) + 
             wui_60_100 + I(wui_60_100^2) +
             crops_60 + I(crops_60^2) + (1|Region), data=dat1)
