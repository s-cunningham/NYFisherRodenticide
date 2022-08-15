library(tidyverse)
library(lme4)
library(parallel)
library(MuMIn)

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

#### Analysis ####

## set up binary variable
dat$binary.MO <- ifelse(dat$n.compounds.MO==0, 0, 1)
dat$binary.T <- ifelse(dat$n.compounds.T==0, 0, 1)

## Reorder columns
dat <- dat[,c(1:7, 31, 8:13, 32, 33, 16:29)]

## Scale and center variables
dat[,c(10,19:30)] <- scale(dat[,c(10,19:30)])

## Use all the data

## Percent AG
pctAG1 <- dat[, c(1:3, 6, 7, 15, 17, 19:21)]
pctAG1 <- distinct(pctAG1)
pctAG1 <- pctAG1 %>% group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

ag15 <- glmer(binary.MO ~ totalag_15 + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
ag30 <- glmer(binary.MO ~ totalag_30 + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
ag60 <- glmer(binary.MO ~ totalag_60 + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
ag15sq <- glmer(binary.MO ~ totalag_15 + I(totalag_15^2) + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
ag30sq <- glmer(binary.MO ~ totalag_30 + I(totalag_30^2) + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
ag60sq <- glmer(binary.MO ~ totalag_60 + I(totalag_60^2) + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)

crops15 <- glmer(binary.MO ~ crops_15 + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
crops30 <- glmer(binary.MO ~ crops_30 + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
crops60 <- glmer(binary.MO ~ crops_60 + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
crops15sq <- glmer(binary.MO ~ crops_15 + I(crops_15^2) + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
crops30sq <- glmer(binary.MO ~ crops_30 + I(crops_30^2) + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
crops60sq <- glmer(binary.MO ~ crops_60 + I(crops_60^2) + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)

past15 <- glmer(binary.MO ~ pasture_15 + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
past30 <- glmer(binary.MO ~ pasture_30 + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
past60 <- glmer(binary.MO ~ pasture_60 + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
past15sq <- glmer(binary.MO ~ pasture_15 + I(pasture_15^2) + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
past30sq <- glmer(binary.MO ~ pasture_30 + I(pasture_30^2) + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)
past60sq <- glmer(binary.MO ~ pasture_60 + I(pasture_60^2) + (1|Region:WMU), family=binomial(link="logit"), data=pctAG1)

pctAG_sel <- model.sel(ag15, ag15sq, crops15, crops15sq, past15, past15sq,
                       ag30, ag30sq, crops30, crops30sq, past30, past30sq,
                       ag60, ag60sq, crops60, crops60sq, past60, past60sq)
pctAG_sel

## Beech basal area
baa1 <- dat[, c(1:3, 6, 7, 15, 17, 26)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=baa, values_fn=unique) %>% as.data.frame()
names(baa1)[7:9] <- c("baa_15", "baa_30", "baa_60") 

baa15 <- glmer(binary.MO ~ baa_15 + (1|Region:WMU), family=binomial(link="logit"), data=baa1)
baa30 <- glmer(binary.MO ~ baa_30 + (1|Region:WMU), family=binomial(link="logit"), data=baa1)
baa60 <- glmer(binary.MO ~ baa_60 + (1|Region:WMU), family=binomial(link="logit"), data=baa1)
baa15sq <- glmer(binary.MO ~ baa_15 + I(baa_15^2) + (1|Region:WMU), family=binomial(link="logit"), data=baa1)
baa30sq <- glmer(binary.MO ~ baa_30 + I(baa_30^2) + (1|Region:WMU), family=binomial(link="logit"), data=baa1)
baa60sq <- glmer(binary.MO ~ baa_60 + I(baa_60^2) + (1|Region:WMU), family=binomial(link="logit"), data=baa1)

# Model selection with MuMIn
baa_sel <- model.sel(baa15, baa30, baa60, baa15sq, baa30sq, baa60sq)
baa_sel

## Intermix WUI
intermix1 <- dat[, c(1:3, 6, 7, 15, 17, 18, 27)]
intermix1 <- unite(intermix1, "buffrad", 7:8, sep="_")
intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[7:15] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                            "mix_15_250", "mix_30_250", "mix_60_250",
                            "mix_15_500", "mix_30_500", "mix_60_500")  

mix_15100 <- glmer(binary.MO ~ mix_15_100 + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_15250 <- glmer(binary.MO ~ mix_15_250 + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_15500 <- glmer(binary.MO ~ mix_15_500 + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_15100sq <- glmer(binary.MO ~ mix_15_100 + I(mix_15_100^2) + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_15250sq <- glmer(binary.MO ~ mix_15_250 + I(mix_15_250^2) + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_15500sq <- glmer(binary.MO ~ mix_15_500 + I(mix_15_500^2) + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)

mix_30100 <- glmer(binary.MO  ~ mix_30_100 + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_30250 <- glmer(binary.MO  ~ mix_30_250 + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_30500 <- glmer(binary.MO  ~ mix_30_500 + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_30100sq <- glmer(binary.MO ~ mix_30_100 + I(mix_30_100^2) + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_30250sq <- glmer(binary.MO ~ mix_30_250 + I(mix_30_250^2) + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_30500sq <- glmer(binary.MO ~ mix_30_500 + I(mix_30_500^2) + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)

mix_60100 <- glmer(binary.MO ~ mix_60_100 + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_60250 <- glmer(binary.MO ~ mix_60_250 + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_60500 <- glmer(binary.MO ~ mix_60_500 + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_60100sq <- glmer(binary.MO ~ mix_60_100 + I(mix_60_100^2) + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_60250sq <- glmer(binary.MO ~ mix_60_250 + I(mix_60_250^2) + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)
mix_60500sq <- glmer(binary.MO ~ mix_60_500 + I(mix_60_500^2) + (1|Region:WMU), family=binomial(link="logit"), data=intermix1)

intermix_sel <- model.sel(mix_15100, mix_30100, mix_60100, mix_15100sq, mix_30100sq, mix_60100sq,
                          mix_15250, mix_30250, mix_60250, mix_15250sq, mix_30250sq, mix_60250sq,
                          mix_15500, mix_30500, mix_60500, mix_15500sq, mix_30500sq, mix_60500sq)
intermix_sel

## Interface WUI
interface1 <- dat[, c(1:3, 6, 7, 15, 17, 18, 28)]
interface1 <- unite(interface1, "buffrad", 7:8, sep="_")
interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[7:15] <- c("face_15_100", "face_30_100", "face_60_100",
                              "face_15_250", "face_30_250", "face_60_250",
                              "face_15_500", "face_30_500", "face_60_500") 

face_15100 <- glmer(binary.MO ~ face_15_100 + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_15250 <- glmer(binary.MO ~ face_15_250 + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_15500 <- glmer(binary.MO ~ face_15_500 + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_15100sq <- glmer(binary.MO ~ face_15_100 + I(face_15_100^2) + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_15250sq <- glmer(binary.MO ~ face_15_250 + I(face_15_250^2) + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_15500sq <- glmer(binary.MO ~ face_15_500 + I(face_15_500^2) + (1|Region:WMU), family=binomial(link="logit"), data=interface1)

face_30100 <- glmer(binary.MO  ~ face_30_100 + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_30250 <- glmer(binary.MO  ~ face_30_250 + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_30500 <- glmer(binary.MO  ~ face_30_500 + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_30100sq <- glmer(binary.MO ~ face_30_100 + I(face_30_100^2) + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_30250sq <- glmer(binary.MO ~ face_30_250 + I(face_30_250^2) + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_30500sq <- glmer(binary.MO ~ face_30_500 + I(face_30_500^2) + (1|Region:WMU), family=binomial(link="logit"), data=interface1)

face_60100 <- glmer(binary.MO ~ face_60_100 + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_60250 <- glmer(binary.MO ~ face_60_250 + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_60500 <- glmer(binary.MO ~ face_60_500 + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_60100sq <- glmer(binary.MO ~ face_60_100 + I(face_60_100^2) + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_60250sq <- glmer(binary.MO ~ face_60_250 + I(face_60_250^2) + (1|Region:WMU), family=binomial(link="logit"), data=interface1)
face_60500sq <- glmer(binary.MO ~ face_60_500 + I(face_60_500^2) + (1|Region:WMU), family=binomial(link="logit"), data=interface1)

interface_sel <- model.sel(face_15100, face_30100, face_60100, face_15100sq, face_30100sq, face_60100sq,
                          face_15250, face_30250, face_60250, face_15250sq, face_30250sq, face_60250sq,
                          face_15500, face_30500, face_60500, face_15500sq, face_30500sq, face_60500sq)
interface_sel

## Total WUI
wui1 <- dat[, c(1:3, 6, 7, 15, 17, 18, 29)]
wui1 <- unite(wui1, "buffrad", 7:8, sep="_")
wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[7:15] <- c("wui_15_100", "wui_30k_100", "wui_60_100",
                        "wui_15_250", "wui_30k_250", "wui_60_250",
                        "wui_15_500", "wui_30k_500", "wui_60_500")

wui_15100 <- glmer(binary.MO ~ wui_15_100 + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_15250 <- glmer(binary.MO ~ wui_15_250 + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_15500 <- glmer(binary.MO ~ wui_15_500 + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_15100sq <- glmer(binary.MO ~ wui_15_100 + I(wui_15_100^2) + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_15250sq <- glmer(binary.MO ~ wui_15_250 + I(wui_15_250^2) + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_15500sq <- glmer(binary.MO ~ wui_15_500 + I(wui_15_500^2) + (1|Region:WMU), family=binomial(link="logit"), data=wui1)

wui_30100 <- glmer(binary.MO ~ wui_30k_100 + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_30250 <- glmer(binary.MO ~ wui_30k_250 + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_30500 <- glmer(binary.MO ~ wui_30k_500 + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_30100sq <- glmer(binary.MO ~ wui_30k_100 + I(wui_30k_100^2) + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_30250sq <- glmer(binary.MO ~ wui_30k_250 + I(wui_30k_250^2) + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_30500sq <- glmer(binary.MO ~ wui_30k_500 + I(wui_30k_500^2) + (1|Region:WMU), family=binomial(link="logit"), data=wui1)

wui_60100 <- glmer(binary.MO ~ wui_60_100 + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_60250 <- glmer(binary.MO ~ wui_60_250 + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_60500 <- glmer(binary.MO ~ wui_60_500 + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_60100sq <- glmer(binary.MO ~ wui_60_100 + I(wui_60_100^2) + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_60250sq <- glmer(binary.MO ~ wui_60_250 + I(wui_60_250^2) + (1|Region:WMU), family=binomial(link="logit"), data=wui1)
wui_60500sq <- glmer(binary.MO ~ wui_60_500 + I(wui_60_500^2) + (1|Region:WMU), family=binomial(link="logit"), data=wui1)

wui_sel <- model.sel(wui_15100, wui_30100, wui_60100, wui_15100sq, wui_30100sq, wui_60100sq, 
                     wui_15250, wui_30250, wui_60250, wui_15250sq, wui_30250sq, wui_60250sq, 
                     wui_15500, wui_30500, wui_60500, wui_15500sq, wui_30500sq, wui_60500sq)
wui_sel


#### Set up data to run for each combination of covariates ####
dat1 <- dat[,c(1:15,30)]
dat1 <- distinct(dat1)

# Join percent agriculture
pctAG1 <- pctAG1[,c(1:3, 9)]
dat1 <- left_join(dat1, pctAG1, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 8:9)]
dat1 <- left_join(dat1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 9)]
dat1 <- left_join(dat1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

# Join interface WUI
interface1 <- interface1[, c(1:3, 9)]
dat1 <- left_join(dat1, interface1, by=c("RegionalID", "pt_name", "pt_index"))

# Join total WUI
wui1 <- wui1[,c(1:3, 9)]
dat1 <- left_join(dat1, wui1, by=c("RegionalID", "pt_name", "pt_index"))

#### Run models ####

## Check correlation matrix
cor(dat1[,16:22])

## Set up global model
g1 <- glmer(binary.MO ~ Sex*Age + laggedBMI + 
                        mix_60_100 + I(mix_60_100^2) +
                        face_60_100 + I(face_60_100^2) +
                        baa_30 + I(baa_30^2) +
                        baa_30:laggedBMI + I(baa_30^2):laggedBMI + 
                        (1|Region:WMU), family=binomial(link="logit"), data=dat1, na.action="na.fail")

g2 <- glmer(binary.MO ~ Sex*Age + laggedBMI + 
                        pasture_60 + I(pasture_60^2) +
                        wui_60_100 + I(wui_60_100^2) + 
                        baa_60 + I(baa_60^2) +
                        baa_60:laggedBMI + I(baa_60^2):laggedBMI +
                        (1|Region:WMU), family=binomial(link="logit"), data=dat1, na.action="na.fail")

## Build model sets
g1_dredge <- dredge(g1, subset=dc(baa_30, laggedBMI, baa_30:laggedBMI) &&
                               dc(baa_30, I(baa_30^2), I(baa_30^2):laggedBMI) &&
                               dc(mix_60_100, I(mix_60_100^2)) &&
                               dc(face_60_100, I(face_60_100^2)))

g2_dredge <- dredge(g1, subset=dc(baa_60, laggedBMI, baa_60:laggedBMI) &&
                               dc(baa_60, I(baa_60^2), I(baa_60^2):laggedBMI) &&
                               dc(wui_60_100, I(wui_60_100^2)) &&
                               dc(pasture_60, I(pasture_60^2)))

# Save dredge tables
save(file="output/dredge_tables_binaryMO.Rdata", list="g1_dredge")
save(file="output/dredge_tables_binaryMO2.Rdata", list="g2_dredge")

# Save as csv file
dh <- as.data.frame(g1_dredge)
write_csv(dh, "output/models_binaryMO.csv")
dh2 <- as.data.frame(g2_dredge)
write_csv(dh2, "output/models_binaryMO2.csv")
