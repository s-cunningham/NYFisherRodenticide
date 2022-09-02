library(tidyverse)
library(lme4)
library(parallel)
library(MuMIn)
library(broom.mixed)

options(scipen=999, digits=3)
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
# Landscape covariates
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

# Individual compounds
dat2 <- read_csv("output/summarized_AR_results.csv")
dat2 <- as.data.frame(dat2)

dat2 <- dat2[dat2$compound=="Diphacinone" | dat2$compound=="Brodifacoum" |
             dat2$compound=="Bromadiolone", c(2,19,21:22) ]

dat <- left_join(dat, dat2, by="RegionalID")

#### Analysis ####

# Change how beech mast is incorporated
dat$mast <- NA 
dat$mast[dat$year==2018] <- "mast"
dat$mast[dat$year==2019] <- "fail"
dat$mast[dat$year==2020] <- "mast"
dat$mast <- as.factor(dat$mast)

# Make age a categorical variable
dat$catAge[dat$Age>=3.5] <- "adult"
dat$catAge[dat$Age==2.5] <- "subadult"
dat$catAge[dat$Ag<2.5] <- "juvenile"

# Make random effects factors
dat$WMU <- as.factor(dat$WMU)
dat$WMUA_code <- as.factor(dat$WMUA_code)
dat$year <- factor(dat$year)

## Reorder columns
dat <- dat[,c(1:7,31,8,9,36,10:13,32:34,16,17,35,29,25,18:20,26:28)]

## Scale and center variables
dat[,c(10,22:29)] <- scale(dat[,c(10,22:29)])

## Use pooled data to determine scale

# Subset by compound
brod <- dat[dat$compound=="Brodifacoum",]
brom <- dat[dat$compound=="Bromadiolone",]
diph <- dat[dat$compound=="Diphacinone",]

#### Brodifacoum ####
## Percent AG
pctAG_brod <- brod[, c(1:3, 6:8,15, 18:19, 24:26)]
pctAG_brod <- distinct(pctAG_brod)
pctAG_brod <- pctAG_brod %>% group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

ag15 <- glmer(bin.exp ~ totalag_15 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
ag30 <- glmer(bin.exp ~ totalag_30 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
ag60 <- glmer(bin.exp ~ totalag_60 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
ag15sq <- glmer(bin.exp ~ totalag_15 + I(totalag_15^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
ag30sq <- glmer(bin.exp ~ totalag_30 + I(totalag_30^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
ag60sq <- glmer(bin.exp ~ totalag_60 + I(totalag_60^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)

crops15 <- glmer(bin.exp ~ crops_15 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
crops30 <- glmer(bin.exp ~ crops_30 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
crops60 <- glmer(bin.exp ~ crops_60 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
crops15sq <- glmer(bin.exp ~ crops_15 + I(crops_15^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
crops30sq <- glmer(bin.exp ~ crops_30 + I(crops_30^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
crops60sq <- glmer(bin.exp ~ crops_60 + I(crops_60^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)

past15 <- glmer(bin.exp ~ pasture_15 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
past30 <- glmer(bin.exp ~ pasture_30 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
past60 <- glmer(bin.exp ~ pasture_60 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
past15sq <- glmer(bin.exp ~ pasture_15 + I(pasture_15^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
past30sq <- glmer(bin.exp ~ pasture_30 + I(pasture_30^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)
past60sq <- glmer(bin.exp ~ pasture_60 + I(pasture_60^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brod)

pctAG_sel <- model.sel(ag15, ag15sq, crops15, crops15sq, past15, past15sq,
                       ag30, ag30sq, crops30, crops30sq, past30, past30sq,
                       ag60, ag60sq, crops60, crops60sq, past60, past60sq)
pctAG_sel

## Beech basal area
baa1 <- brod[, c(1:3, 6:8,15, 18:19, 21:23)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=baa, values_fn=unique) %>% as.data.frame()
names(baa1)[11:13] <- c("baa_15", "baa_30", "baa_60") 

baa15 <- glmer(bin.exp ~ baa_15 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa15sq <- glmer(bin.exp ~ baa_15 + I(baa_15^2) + (1|WMUA_code/WMU) + (1|year),family=binomial(link="logit"),  data=baa1)
baa30 <- glmer(bin.exp ~ baa_30 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa30sq <- glmer(bin.exp ~ baa_30 + I(baa_30^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa60 <- glmer(bin.exp ~ baa_60 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa60sq <- glmer(bin.exp ~ baa_60 + I(baa_60^2) + (1|WMUA_code/WMU) + (1|year),family=binomial(link="logit"),  data=baa1)

baa15lBMI <- glmer(bin.exp ~ baa_15*laggedBMI + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa30lBMI <- glmer(bin.exp ~ baa_30*laggedBMI + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa60lBMI <- glmer(bin.exp ~ baa_60*laggedBMI + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)

baa15M <- glmer(bin.exp ~ baa_15*mast + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa30M <- glmer(bin.exp ~ baa_30*mast + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa60M <- glmer(bin.exp ~ baa_60*mast + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)

baa_sel <- model.sel(baa15, baa30, baa60, baa15sq, baa30sq, baa60sq, 
                     baa15lBMI, baa30lBMI, baa60lBMI, 
                     baa15M, baa30M, baa60M)
baa_sel

## Wildland-urban interface
# Intermix
intermix1 <- brod[, c(1:3, 6:8,15, 18:20, 27)]
intermix1 <- unite(intermix1, "buffrad", 9:10, sep="_")
intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[9:17] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                            "mix_15_250", "mix_30_250", "mix_60_250",
                            "mix_15_500", "mix_30_500", "mix_60_500")  

# Interface
interface1 <- brod[, c(1:3, 6:8,15, 18:20,28)]
interface1 <- unite(interface1, "buffrad", 9:10, sep="_")
interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[9:17] <- c("face_15_100", "face_30_100", "face_60_100",
                             "face_15_250", "face_30_250", "face_60_250",
                             "face_15_500", "face_30_500", "face_60_500") 

# WUI total
wui1 <- brod[, c(1:3, 6:8,15, 18:20, 29)]
wui1 <- unite(wui1, "buffrad", 9:10, sep="_")
wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[9:17] <- c("wui_15_100", "wui_30k_100", "wui_60_100",
                       "wui_15_250", "wui_30k_250", "wui_60_250",
                       "wui_15_500", "wui_30k_500", "wui_60_500")

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 9:17)]
wui1 <- left_join(wui1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

# Join interface WUI
interface1 <- interface1[, c(1:3, 9:17)]
wui1 <- left_join(wui1, interface1, by=c("RegionalID", "pt_name", "pt_index"))

# Intermix
mix_15100 <- glmer(bin.exp ~ mix_15_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15250 <- glmer(bin.exp ~ mix_15_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15500 <- glmer(bin.exp ~ mix_15_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15100sq <- glmer(bin.exp ~ mix_15_100 + I(mix_15_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15250sq <- glmer(bin.exp ~ mix_15_250 + I(mix_15_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15500sq <- glmer(bin.exp ~ mix_15_500 + I(mix_15_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

mix_30100 <- glmer(bin.exp  ~ mix_30_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30250 <- glmer(bin.exp  ~ mix_30_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30500 <- glmer(bin.exp  ~ mix_30_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30100sq <- glmer(bin.exp ~ mix_30_100 + I(mix_30_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30250sq <- glmer(bin.exp ~ mix_30_250 + I(mix_30_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30500sq <- glmer(bin.exp ~ mix_30_500 + I(mix_30_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

mix_60100 <- glmer(bin.exp ~ mix_60_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60250 <- glmer(bin.exp ~ mix_60_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60500 <- glmer(bin.exp ~ mix_60_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60100sq <- glmer(bin.exp ~ mix_60_100 + I(mix_60_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60250sq <- glmer(bin.exp ~ mix_60_250 + I(mix_60_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60500sq <- glmer(bin.exp ~ mix_60_500 + I(mix_60_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

# Interface
face_15100 <- glmer(bin.exp ~ face_15_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15250 <- glmer(bin.exp ~ face_15_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15500 <- glmer(bin.exp ~ face_15_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15100sq <- glmer(bin.exp ~ face_15_100 + I(face_15_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15250sq <- glmer(bin.exp ~ face_15_250 + I(face_15_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15500sq <- glmer(bin.exp ~ face_15_500 + I(face_15_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

face_30100 <- glmer(bin.exp  ~ face_30_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30250 <- glmer(bin.exp  ~ face_30_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30500 <- glmer(bin.exp  ~ face_30_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30100sq <- glmer(bin.exp ~ face_30_100 + I(face_30_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30250sq <- glmer(bin.exp ~ face_30_250 + I(face_30_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30500sq <- glmer(bin.exp ~ face_30_500 + I(face_30_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

face_60100 <- glmer(bin.exp ~ face_60_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60250 <- glmer(bin.exp ~ face_60_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60500 <- glmer(bin.exp ~ face_60_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60100sq <- glmer(bin.exp ~ face_60_100 + I(face_60_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60250sq <- glmer(bin.exp ~ face_60_250 + I(face_60_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60500sq <- glmer(bin.exp ~ face_60_500 + I(face_60_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

# Total WUI
wui_15100 <- glmer(bin.exp ~ wui_15_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15250 <- glmer(bin.exp ~ wui_15_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15500 <- glmer(bin.exp ~ wui_15_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15100sq <- glmer(bin.exp ~ wui_15_100 + I(wui_15_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15250sq <- glmer(bin.exp ~ wui_15_250 + I(wui_15_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15500sq <- glmer(bin.exp ~ wui_15_500 + I(wui_15_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

wui_30100 <- glmer(bin.exp ~ wui_30k_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30250 <- glmer(bin.exp ~ wui_30k_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30500 <- glmer(bin.exp ~ wui_30k_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30100sq <- glmer(bin.exp ~ wui_30k_100 + I(wui_30k_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30250sq <- glmer(bin.exp ~ wui_30k_250 + I(wui_30k_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30500sq <- glmer(bin.exp ~ wui_30k_500 + I(wui_30k_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

wui_60100 <- glmer(bin.exp ~ wui_60_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60250 <- glmer(bin.exp ~ wui_60_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60500 <- glmer(bin.exp ~ wui_60_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60100sq <- glmer(bin.exp ~ wui_60_100 + I(wui_60_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60250sq <- glmer(bin.exp ~ wui_60_250 + I(wui_60_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60500sq <- glmer(bin.exp ~ wui_60_500 + I(wui_60_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

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

## Set up data to run for each combination of covariates ##
brod1 <- brod[,c(1:18)]
brod1 <- distinct(brod1)

# Join percent agriculture
pctAG_brod <- pctAG_brod[,c(1:3, 11)]
brod1 <- left_join(brod1, pctAG_brod, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 10, 13)]
brod1 <- left_join(brod1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join WUI
wui1 <- wui1[,c(1:3, 26, 30)]
brod1 <- left_join(brod1, wui1, by=c("RegionalID", "pt_name", "pt_index"))

## Run models ##

## Check correlation matrix
cor(brod1[,19:23])

## Set up global model
g1 <- glmer(bin.exp ~ Sex*catAge + laggedBMI + 
              mix_60_500 + I(mix_60_500^2) +
              face_15_250 + I(face_15_250^2) +
              pasture_60 + I(pasture_60^2) +
              baa_60 + I(baa_60^2) + baa_60:laggedBMI +  
              (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=brod1, na.action="na.fail")


# Export data and model into the cluster worker nodes
clusterExport(cl, c("brod1","g1"))

## Build model sets
g1_dredge <- MuMIn:::.dredge.par(g1, cluster=cl, trace=2,  subset=dc(mix_60_500, I(mix_60_500^2)) &&
                                                                    dc(pasture_60, I(pasture_60^2)) &&
                                                                    dc(face_15_250, I(face_15_250^2)) &&
                                                                    dc(baa_60:laggedBMI, I(baa_60^2):laggedBMI))


# Save dredge tables
save(file="output/dredge_tables_binaryT.Rdata", list="g1_dredge")

# Save as csv file
dh <- as.data.frame(g1_dredge)
write_csv(dh, "output/models_binaryT.csv")

#### Read tables back in
dh <- read_csv("output/models_binaryT.csv")
dh <- as.data.frame(dh)

## Running final models ##

## Loop over each set of random points

m_est <- data.frame()
m_stderr <- data.frame()
pct2.5 <- data.frame()
pct97.5 <- data.frame()
set.seed(123)

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- brod1[brod1$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- glmer(bin.exp ~ Sex + catAge + pasture_60 + I(pasture_60^2) + 
                  mix_60_500 + I(mix_60_500^2) + face_15_250 + (1|year), 
                  family=binomial(link="logit"), data=pt, nAGQ=0,
                  control=glmerControl(optimizer=c("nloptwrap", "nloptwrap")))
  
  m1s <- summary(m1_pt)
  m1sdf <- m1s$coefficients
  
  # save averaged confidence intervals
  ci <- confint.merMod(m1_pt, method="Wald")
  pct2.5 <- rbind(pct2.5, t(ci)[1,2:10])
  pct97.5 <- rbind(pct97.5, t(ci)[2,2:10])
  
  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- row.names(m1sdf)
    names(m_stderr) <- row.names(m1sdf)
    names(pct2.5) <- row.names(m1sdf)
    names(pct97.5) <- row.names(m1sdf)
  }
  
}

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)

# One-sample t-test to determine "significance"
pvalue <- c()
for (i in 1:ncol(m_est)) {
  
  tresult <- t.test(m_est[,i], mu=0, alternative="two.sided")
  pvalue <- c(pvalue, tresult$p.value)
}
pvalue <- as.data.frame(pvalue)
pvalue <- cbind(names(coef_avg), pvalue)
pvalue <- pivot_wider(pvalue, names_from="names(coef_avg)", values_from="pvalue") %>% as.data.frame()

# Combine and clean up data frame
coef_summary <- bind_rows(coef_avg, stderr_avg, pct2.5_avg, pct97.5_avg, pvalue)
names(coef_summary) <- row.names(m1sdf)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "2.5CI", "97.5CI", "P-value")
coef_summary <- data.frame(coef=coefs, coef_summary)

# Write to file
write_csv(coef_summary, "results/binaryTbrodifacoum_coef-summary.csv")

#### Bromadiolone ####

## Percent AG
pctAG_brom <- brom[, c(1:3, 6:8,15, 18:19, 24:26)]
pctAG_brom <- distinct(pctAG_brom)
pctAG_brom <- pctAG_brom %>% group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

ag15 <- glmer(bin.exp ~ totalag_15 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
ag30 <- glmer(bin.exp ~ totalag_30 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
ag60 <- glmer(bin.exp ~ totalag_60 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
ag15sq <- glmer(bin.exp ~ totalag_15 + I(totalag_15^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
ag30sq <- glmer(bin.exp ~ totalag_30 + I(totalag_30^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
ag60sq <- glmer(bin.exp ~ totalag_60 + I(totalag_60^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)

crops15 <- glmer(bin.exp ~ crops_15 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
crops30 <- glmer(bin.exp ~ crops_30 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
crops60 <- glmer(bin.exp ~ crops_60 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
crops15sq <- glmer(bin.exp ~ crops_15 + I(crops_15^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
crops30sq <- glmer(bin.exp ~ crops_30 + I(crops_30^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
crops60sq <- glmer(bin.exp ~ crops_60 + I(crops_60^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)

past15 <- glmer(bin.exp ~ pasture_15 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
past30 <- glmer(bin.exp ~ pasture_30 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
past60 <- glmer(bin.exp ~ pasture_60 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
past15sq <- glmer(bin.exp ~ pasture_15 + I(pasture_15^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
past30sq <- glmer(bin.exp ~ pasture_30 + I(pasture_30^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)
past60sq <- glmer(bin.exp ~ pasture_60 + I(pasture_60^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_brom)

pctAG_sel <- model.sel(ag15, ag15sq, crops15, crops15sq, past15, past15sq,
                       ag30, ag30sq, crops30, crops30sq, past30, past30sq,
                       ag60, ag60sq, crops60, crops60sq, past60, past60sq)
pctAG_sel

## Beech basal area
baa1 <- brom[, c(1:3, 6:8,15, 18:19, 21:23)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=baa, values_fn=unique) %>% as.data.frame()
names(baa1)[11:13] <- c("baa_15", "baa_30", "baa_60") 

baa15 <- glmer(bin.exp ~ baa_15 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa15sq <- glmer(bin.exp ~ baa_15 + I(baa_15^2) + (1|WMUA_code/WMU) + (1|year),family=binomial(link="logit"),  data=baa1)
baa30 <- glmer(bin.exp ~ baa_30 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa30sq <- glmer(bin.exp ~ baa_30 + I(baa_30^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa60 <- glmer(bin.exp ~ baa_60 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa60sq <- glmer(bin.exp ~ baa_60 + I(baa_60^2) + (1|WMUA_code/WMU) + (1|year),family=binomial(link="logit"),  data=baa1)

baa15lBMI <- glmer(bin.exp ~ baa_15*laggedBMI + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa30lBMI <- glmer(bin.exp ~ baa_30*laggedBMI + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa60lBMI <- glmer(bin.exp ~ baa_60*laggedBMI + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)

baa15M <- glmer(bin.exp ~ baa_15*mast + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa30M <- glmer(bin.exp ~ baa_30*mast + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa60M <- glmer(bin.exp ~ baa_60*mast + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)

baa_sel <- model.sel(baa15, baa30, baa60, baa15sq, baa30sq, baa60sq, 
                     baa15lBMI, baa30lBMI, baa60lBMI, 
                     baa15M, baa30M, baa60M)
baa_sel

## Wildland-urban interface
# Intermix
intermix1 <- brom[, c(1:3, 6:8,15, 18:20, 27)]
intermix1 <- unite(intermix1, "buffrad", 9:10, sep="_")
intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[9:17] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                            "mix_15_250", "mix_30_250", "mix_60_250",
                            "mix_15_500", "mix_30_500", "mix_60_500")  

# Interface
interface1 <- brom[, c(1:3, 6:8,15, 18:20,28)]
interface1 <- unite(interface1, "buffrad", 9:10, sep="_")
interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[9:17] <- c("face_15_100", "face_30_100", "face_60_100",
                             "face_15_250", "face_30_250", "face_60_250",
                             "face_15_500", "face_30_500", "face_60_500") 

# WUI total
wui1 <- brom[, c(1:3, 6:8,15, 18:20, 29)]
wui1 <- unite(wui1, "buffrad", 9:10, sep="_")
wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[9:17] <- c("wui_15_100", "wui_30k_100", "wui_60_100",
                       "wui_15_250", "wui_30k_250", "wui_60_250",
                       "wui_15_500", "wui_30k_500", "wui_60_500")

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 9:17)]
wui1 <- left_join(wui1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

# Join interface WUI
interface1 <- interface1[, c(1:3, 9:17)]
wui1 <- left_join(wui1, interface1, by=c("RegionalID", "pt_name", "pt_index"))

# Intermix
mix_15100 <- glmer(bin.exp ~ mix_15_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15250 <- glmer(bin.exp ~ mix_15_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15500 <- glmer(bin.exp ~ mix_15_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15100sq <- glmer(bin.exp ~ mix_15_100 + I(mix_15_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15250sq <- glmer(bin.exp ~ mix_15_250 + I(mix_15_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15500sq <- glmer(bin.exp ~ mix_15_500 + I(mix_15_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

mix_30100 <- glmer(bin.exp  ~ mix_30_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30250 <- glmer(bin.exp  ~ mix_30_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30500 <- glmer(bin.exp  ~ mix_30_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30100sq <- glmer(bin.exp ~ mix_30_100 + I(mix_30_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30250sq <- glmer(bin.exp ~ mix_30_250 + I(mix_30_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30500sq <- glmer(bin.exp ~ mix_30_500 + I(mix_30_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

mix_60100 <- glmer(bin.exp ~ mix_60_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60250 <- glmer(bin.exp ~ mix_60_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60500 <- glmer(bin.exp ~ mix_60_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60100sq <- glmer(bin.exp ~ mix_60_100 + I(mix_60_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60250sq <- glmer(bin.exp ~ mix_60_250 + I(mix_60_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60500sq <- glmer(bin.exp ~ mix_60_500 + I(mix_60_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

# Interface
face_15100 <- glmer(bin.exp ~ face_15_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15250 <- glmer(bin.exp ~ face_15_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15500 <- glmer(bin.exp ~ face_15_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15100sq <- glmer(bin.exp ~ face_15_100 + I(face_15_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15250sq <- glmer(bin.exp ~ face_15_250 + I(face_15_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15500sq <- glmer(bin.exp ~ face_15_500 + I(face_15_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

face_30100 <- glmer(bin.exp  ~ face_30_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30250 <- glmer(bin.exp  ~ face_30_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30500 <- glmer(bin.exp  ~ face_30_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30100sq <- glmer(bin.exp ~ face_30_100 + I(face_30_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30250sq <- glmer(bin.exp ~ face_30_250 + I(face_30_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30500sq <- glmer(bin.exp ~ face_30_500 + I(face_30_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

face_60100 <- glmer(bin.exp ~ face_60_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60250 <- glmer(bin.exp ~ face_60_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60500 <- glmer(bin.exp ~ face_60_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60100sq <- glmer(bin.exp ~ face_60_100 + I(face_60_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60250sq <- glmer(bin.exp ~ face_60_250 + I(face_60_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60500sq <- glmer(bin.exp ~ face_60_500 + I(face_60_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

# Total WUI
wui_15100 <- glmer(bin.exp ~ wui_15_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15250 <- glmer(bin.exp ~ wui_15_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15500 <- glmer(bin.exp ~ wui_15_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15100sq <- glmer(bin.exp ~ wui_15_100 + I(wui_15_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15250sq <- glmer(bin.exp ~ wui_15_250 + I(wui_15_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15500sq <- glmer(bin.exp ~ wui_15_500 + I(wui_15_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

wui_30100 <- glmer(bin.exp ~ wui_30k_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30250 <- glmer(bin.exp ~ wui_30k_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30500 <- glmer(bin.exp ~ wui_30k_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30100sq <- glmer(bin.exp ~ wui_30k_100 + I(wui_30k_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30250sq <- glmer(bin.exp ~ wui_30k_250 + I(wui_30k_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30500sq <- glmer(bin.exp ~ wui_30k_500 + I(wui_30k_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

wui_60100 <- glmer(bin.exp ~ wui_60_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60250 <- glmer(bin.exp ~ wui_60_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60500 <- glmer(bin.exp ~ wui_60_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60100sq <- glmer(bin.exp ~ wui_60_100 + I(wui_60_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60250sq <- glmer(bin.exp ~ wui_60_250 + I(wui_60_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60500sq <- glmer(bin.exp ~ wui_60_500 + I(wui_60_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

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

## Set up data to run for each combination of covariates ##
brom1 <- brom[,c(1:18)]
brom1 <- distinct(brom1)

# Join percent agriculture
pctAG_brom <- pctAG_brom[,c(1:3, 11)]
brom1 <- left_join(brom1, pctAG_brom, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 13)]
brom1 <- left_join(brom1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join WUI
wui1 <- wui1[,c(1:3, 26)]
brom1 <- left_join(brom1, wui1, by=c("RegionalID", "pt_name", "pt_index"))

## Run models ##

## Check correlation matrix
cor(brom1[,19:21])

## Set up global model
g1 <- glmer(bin.exp ~ Sex*catAge + 
              mix_60_500 + I(mix_60_500^2) +
              pasture_60 + I(pasture_60^2) +
              baa_60 + I(baa_60^2) + 
              (1|year), family=binomial(link="logit"), 
              data=brom1, na.action="na.fail", nAGQ=0,
              control=glmerControl(optimizer=c("nloptwrap", "nloptwrap")))


# Export data and model into the cluster worker nodes
clusterExport(cl, c("brom1","g1"))

## Build model sets
g1_dredge <- MuMIn:::.dredge.par(g1, cluster=cl, trace=2,  subset=dc(mix_60_500, I(mix_60_500^2)) &&
                                   dc(pasture_60, I(pasture_60^2)))


# Save dredge tables
save(file="output/dredge_tables_binaryTbroma.Rdata", list="g1_dredge")

# Save as csv file
dh <- as.data.frame(g1_dredge)
write_csv(dh, "output/models_binaryTbroma.csv")

#### Read tables back in
dh <- read_csv("output/models_binaryTbroma.csv")
dh <- as.data.frame(dh)

## Running final models ##

## Loop over each set of random points

m_est <- data.frame()
m_stderr <- data.frame()
pct2.5 <- data.frame()
pct97.5 <- data.frame()
set.seed(123)

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- brom1[brom1$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- glmer(bin.exp ~ Sex*catAge + pasture_60 + 
                   mix_60_500 + I(mix_60_500^2) + 
                   (1|WMUA_code), family=binomial(link="logit"),data=pt, nAGQ=0,
                 control=glmerControl(optimizer=c("nloptwrap", "nloptwrap")))
  
  m1s <- summary(m1_pt)
  m1sdf <- m1s$coefficients
  
  # save averaged confidence intervals
  ci <- confint.merMod(m1_pt, method="Wald")
  pct2.5 <- rbind(pct2.5, t(ci)[1,2:10])
  pct97.5 <- rbind(pct97.5, t(ci)[2,2:10])
  
  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- row.names(m1sdf)
    names(m_stderr) <- row.names(m1sdf)
    names(pct2.5) <- row.names(m1sdf)
    names(pct97.5) <- row.names(m1sdf)
  }
  
}

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)

# One-sample t-test to determine "significance"
pvalue <- c()
for (i in 1:ncol(m_est)) {
  
  tresult <- t.test(m_est[,i], mu=0, alternative="two.sided")
  pvalue <- c(pvalue, tresult$p.value)
}
pvalue <- as.data.frame(pvalue)
pvalue <- cbind(names(coef_avg), pvalue)
pvalue <- pivot_wider(pvalue, names_from="names(coef_avg)", values_from="pvalue") %>% as.data.frame()

# Combine and clean up data frame
coef_summary <- bind_rows(coef_avg, stderr_avg, pct2.5_avg, pct97.5_avg, pvalue)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "2.5CI", "97.5CI", "P-value")
coef_summary <- data.frame(coef=coefs, coef_summary)

# Write to file
write_csv(coef_summary, "results/binaryTbromadiolone_coef-summary.csv")


#### Diphacinone ####

## Percent AG
pctAG_diph <- diph[, c(1:3, 6:8,15, 18:19, 24:26)]
pctAG_diph <- distinct(pctAG_diph)
pctAG_diph <- pctAG_diph %>% group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

ag15 <- glmer(bin.exp ~ totalag_15 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
ag30 <- glmer(bin.exp ~ totalag_30 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
ag60 <- glmer(bin.exp ~ totalag_60 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
ag15sq <- glmer(bin.exp ~ totalag_15 + I(totalag_15^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
ag30sq <- glmer(bin.exp ~ totalag_30 + I(totalag_30^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
ag60sq <- glmer(bin.exp ~ totalag_60 + I(totalag_60^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)

crops15 <- glmer(bin.exp ~ crops_15 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
crops30 <- glmer(bin.exp ~ crops_30 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
crops60 <- glmer(bin.exp ~ crops_60 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
crops15sq <- glmer(bin.exp ~ crops_15 + I(crops_15^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
crops30sq <- glmer(bin.exp ~ crops_30 + I(crops_30^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
crops60sq <- glmer(bin.exp ~ crops_60 + I(crops_60^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)

past15 <- glmer(bin.exp ~ pasture_15 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
past30 <- glmer(bin.exp ~ pasture_30 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
past60 <- glmer(bin.exp ~ pasture_60 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
past15sq <- glmer(bin.exp ~ pasture_15 + I(pasture_15^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
past30sq <- glmer(bin.exp ~ pasture_30 + I(pasture_30^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)
past60sq <- glmer(bin.exp ~ pasture_60 + I(pasture_60^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=pctAG_diph)

pctAG_sel <- model.sel(ag15, ag15sq, crops15, crops15sq, past15, past15sq,
                       ag30, ag30sq, crops30, crops30sq, past30, past30sq,
                       ag60, ag60sq, crops60, crops60sq, past60, past60sq)
pctAG_sel

## Beech basal area
baa1 <- diph[, c(1:3, 6:8,15, 18:19, 21:23)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=baa, values_fn=unique) %>% as.data.frame()
names(baa1)[11:13] <- c("baa_15", "baa_30", "baa_60") 

baa15 <- glmer(bin.exp ~ baa_15 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa15sq <- glmer(bin.exp ~ baa_15 + I(baa_15^2) + (1|WMUA_code/WMU) + (1|year),family=binomial(link="logit"),  data=baa1)
baa30 <- glmer(bin.exp ~ baa_30 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa30sq <- glmer(bin.exp ~ baa_30 + I(baa_30^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa60 <- glmer(bin.exp ~ baa_60 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa60sq <- glmer(bin.exp ~ baa_60 + I(baa_60^2) + (1|WMUA_code/WMU) + (1|year),family=binomial(link="logit"),  data=baa1)

baa15lBMI <- glmer(bin.exp ~ baa_15*laggedBMI + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa30lBMI <- glmer(bin.exp ~ baa_30*laggedBMI + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa60lBMI <- glmer(bin.exp ~ baa_60*laggedBMI + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)

baa15M <- glmer(bin.exp ~ baa_15*mast + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa30M <- glmer(bin.exp ~ baa_30*mast + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)
baa60M <- glmer(bin.exp ~ baa_60*mast + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=baa1)

baa_sel <- model.sel(baa15, baa30, baa60, baa15sq, baa30sq, baa60sq, 
                     baa15lBMI, baa30lBMI, baa60lBMI, 
                     baa15M, baa30M, baa60M)
baa_sel

## Wildland-urban interface
# Intermix
intermix1 <- diph[, c(1:3, 6:8,15, 18:20, 27)]
intermix1 <- unite(intermix1, "buffrad", 9:10, sep="_")
intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[9:17] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                            "mix_15_250", "mix_30_250", "mix_60_250",
                            "mix_15_500", "mix_30_500", "mix_60_500")  

# Interface
interface1 <- diph[, c(1:3, 6:8,15, 18:20,28)]
interface1 <- unite(interface1, "buffrad", 9:10, sep="_")
interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[9:17] <- c("face_15_100", "face_30_100", "face_60_100",
                             "face_15_250", "face_30_250", "face_60_250",
                             "face_15_500", "face_30_500", "face_60_500") 

# WUI total
wui1 <- diph[, c(1:3, 6:8,15, 18:20, 29)]
wui1 <- unite(wui1, "buffrad", 9:10, sep="_")
wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[9:17] <- c("wui_15_100", "wui_30k_100", "wui_60_100",
                       "wui_15_250", "wui_30k_250", "wui_60_250",
                       "wui_15_500", "wui_30k_500", "wui_60_500")

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 9:17)]
wui1 <- left_join(wui1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

# Join interface WUI
interface1 <- interface1[, c(1:3, 9:17)]
wui1 <- left_join(wui1, interface1, by=c("RegionalID", "pt_name", "pt_index"))

# Intermix
mix_15100 <- glmer(bin.exp ~ mix_15_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15250 <- glmer(bin.exp ~ mix_15_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15500 <- glmer(bin.exp ~ mix_15_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15100sq <- glmer(bin.exp ~ mix_15_100 + I(mix_15_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15250sq <- glmer(bin.exp ~ mix_15_250 + I(mix_15_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_15500sq <- glmer(bin.exp ~ mix_15_500 + I(mix_15_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

mix_30100 <- glmer(bin.exp  ~ mix_30_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30250 <- glmer(bin.exp  ~ mix_30_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30500 <- glmer(bin.exp  ~ mix_30_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30100sq <- glmer(bin.exp ~ mix_30_100 + I(mix_30_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30250sq <- glmer(bin.exp ~ mix_30_250 + I(mix_30_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_30500sq <- glmer(bin.exp ~ mix_30_500 + I(mix_30_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

mix_60100 <- glmer(bin.exp ~ mix_60_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60250 <- glmer(bin.exp ~ mix_60_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60500 <- glmer(bin.exp ~ mix_60_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60100sq <- glmer(bin.exp ~ mix_60_100 + I(mix_60_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60250sq <- glmer(bin.exp ~ mix_60_250 + I(mix_60_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
mix_60500sq <- glmer(bin.exp ~ mix_60_500 + I(mix_60_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

# Interface
face_15100 <- glmer(bin.exp ~ face_15_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15250 <- glmer(bin.exp ~ face_15_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15500 <- glmer(bin.exp ~ face_15_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15100sq <- glmer(bin.exp ~ face_15_100 + I(face_15_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15250sq <- glmer(bin.exp ~ face_15_250 + I(face_15_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_15500sq <- glmer(bin.exp ~ face_15_500 + I(face_15_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

face_30100 <- glmer(bin.exp  ~ face_30_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30250 <- glmer(bin.exp  ~ face_30_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30500 <- glmer(bin.exp  ~ face_30_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30100sq <- glmer(bin.exp ~ face_30_100 + I(face_30_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30250sq <- glmer(bin.exp ~ face_30_250 + I(face_30_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_30500sq <- glmer(bin.exp ~ face_30_500 + I(face_30_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

face_60100 <- glmer(bin.exp ~ face_60_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60250 <- glmer(bin.exp ~ face_60_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60500 <- glmer(bin.exp ~ face_60_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60100sq <- glmer(bin.exp ~ face_60_100 + I(face_60_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60250sq <- glmer(bin.exp ~ face_60_250 + I(face_60_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
face_60500sq <- glmer(bin.exp ~ face_60_500 + I(face_60_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

# Total WUI
wui_15100 <- glmer(bin.exp ~ wui_15_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15250 <- glmer(bin.exp ~ wui_15_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15500 <- glmer(bin.exp ~ wui_15_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15100sq <- glmer(bin.exp ~ wui_15_100 + I(wui_15_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15250sq <- glmer(bin.exp ~ wui_15_250 + I(wui_15_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_15500sq <- glmer(bin.exp ~ wui_15_500 + I(wui_15_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

wui_30100 <- glmer(bin.exp ~ wui_30k_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30250 <- glmer(bin.exp ~ wui_30k_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30500 <- glmer(bin.exp ~ wui_30k_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30100sq <- glmer(bin.exp ~ wui_30k_100 + I(wui_30k_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30250sq <- glmer(bin.exp ~ wui_30k_250 + I(wui_30k_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_30500sq <- glmer(bin.exp ~ wui_30k_500 + I(wui_30k_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

wui_60100 <- glmer(bin.exp ~ wui_60_100 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60250 <- glmer(bin.exp ~ wui_60_250 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60500 <- glmer(bin.exp ~ wui_60_500 + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60100sq <- glmer(bin.exp ~ wui_60_100 + I(wui_60_100^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60250sq <- glmer(bin.exp ~ wui_60_250 + I(wui_60_250^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)
wui_60500sq <- glmer(bin.exp ~ wui_60_500 + I(wui_60_500^2) + (1|WMUA_code/WMU) + (1|year), family=binomial(link="logit"), data=wui1)

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

## Set up data to run for each combination of covariates ##
diph1 <- diph[,c(1:18)]
diph1 <- distinct(diph1)

# Join percent agriculture
pctAG_diph <- pctAG_diph[,c(1:3, 11)]
diph1 <- left_join(diph1, pctAG_diph, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 10, 13)]
diph1 <- left_join(diph1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join WUI
wui1 <- wui1[,c(1:3, 9)]
diph1 <- left_join(diph1, wui1, by=c("RegionalID", "pt_name", "pt_index"))

## Run models ##

## Check correlation matrix
cor(diph1[,19:22])

## Set up global model
g1 <- glmer(bin.exp ~ Sex*catAge + laggedBMI + 
              wui_15_100 + I(wui_15_100^2) +
              pasture_60 + I(pasture_60^2) +
              baa_60 + I(baa_60^2) + baa_60:laggedBMI +  
              (1|year), family=binomial(link="logit"), data=diph1, na.action="na.fail",
              nAGQ=0, control=glmerControl(optimizer=c("nloptwrap", "nloptwrap")))


# Export data and model into the cluster worker nodes
clusterExport(cl, c("diph1","g1"))

## Build model sets
g1_dredge <- MuMIn:::.dredge.par(g1, cluster=cl, trace=2,  subset=dc(wui_15_100, I(wui_15_100^2)) &&
                                   dc(pasture_60, I(pasture_60^2)) &&
                                   dc(baa_60:laggedBMI, I(baa_60^2):laggedBMI))


# Save dredge tables
save(file="output/dredge_tables_binaryTdiphacinone.Rdata", list="g1_dredge")

# Save as csv file
dh <- as.data.frame(g1_dredge)
write_csv(dh, "output/models_binaryTdipha.csv")

#### Read tables back in
dh <- read_csv("output/models_binaryTdipha.csv")
dh <- as.data.frame(dh)

# Clear deltaAIC and model weights
dh$delta <- NA
dh$weight <- NA

# Reorder according to AICc
N.order <- order(dh$AICc)
dh <- dh[N.order,]

# Remove some that didn't end up in subset of dredge
dh <- dh[!(is.na(dh$baa_60) & !is.na(dh$`I(baa_60^2)`)),]

# Recalculate deltaAIC and model weights
dh$delta <- dh$AICc - dh$AICc[1]
w <- qpcR::akaike.weights(dh$AICc)
dh$weight <- w$weights

## Running final models ##

## Loop over each set of random points

m_est <- data.frame()
m_stderr <- data.frame()
pct2.5 <- data.frame()
pct97.5 <- data.frame()
set.seed(123)

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- diph1[diph1$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- glm(bin.exp ~ Sex + catAge + pasture_60 + wui_15_100 + baa_60, family=binomial(link="logit"), data=pt) #,
                    # control=glmerControl(optimizer=c("nloptwrap", "nloptwrap")))
  
  m1s <- summary(m1_pt)
  m1sdf <- m1s$coefficients
  
  # save averaged confidence intervals
  ci <- confint.merMod(m1_pt, method="Wald")
  pct2.5 <- rbind(pct2.5, t(ci)[1,2:10])
  pct97.5 <- rbind(pct97.5, t(ci)[2,2:10])
  
  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- row.names(m1sdf)
    names(m_stderr) <- row.names(m1sdf)
    names(pct2.5) <- row.names(m1sdf)
    names(pct97.5) <- row.names(m1sdf)
  }
  
}

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)

# One-sample t-test to determine "significance"
pvalue <- c()
for (i in 1:ncol(m_est)) {
  
  tresult <- t.test(m_est[,i], mu=0, alternative="two.sided")
  pvalue <- c(pvalue, tresult$p.value)
}
pvalue <- as.data.frame(pvalue)
pvalue <- cbind(names(coef_avg), pvalue)
pvalue <- pivot_wider(pvalue, names_from="names(coef_avg)", values_from="pvalue") %>% as.data.frame()

# Combine and clean up data frame
coef_summary <- bind_rows(coef_avg, stderr_avg, pct2.5_avg, pct97.5_avg, pvalue)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "2.5CI", "97.5CI", "P-value")
coef_summary <- data.frame(coef=coefs, coef_summary)

# Write to file
write_csv(coef_summary, "results/binaryTdiphacinone_coef-summary.csv")


