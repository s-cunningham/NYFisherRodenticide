library(tidyverse)
library(MuMIn)
library(ordinal)


set.seed(123)

#### Read in data ####
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)
dat <- separate(dat, 10, into=c("del", "Town"), sep="-", remove=FALSE)
dat <- dat[,-11]

#### Analysis Set-up ####

# Order response
dat$n.compounds.MO <- ordered(dat$n.compounds.MO, levels=c(0,1,2,3))
dat$n.compounds.T <- ordered(dat$n.compounds.T, levels=c(0,1,2,3,4,5))

dat <- dat[dat$buffsize>10,]

## Percent AG
pctAG1 <- dat[, c(1:17, 19, 26)]
pctAG1 <- pctAG1 %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=pct_ag, values_fn=unique) %>% as.data.frame()
names(pctAG1)[18:20] <- c("km15", "km30", "km60") 

## Scale and center variables
pctAG1[,c(18:20)] <- scale(pctAG1[,c(18:20)])

# Run models
ag15 <- clmm(n.compounds.T ~ km15 + (1|Region) + (1|WMU), data=pctAG1)
ag15sq <- clmm(n.compounds.T ~ km15 + I(km15^2) + (1|Region) + (1|WMU), data=pctAG1)
ag30 <- clmm(n.compounds.T ~ km30 + (1|Region) + (1|WMU), data=pctAG1)
ag30sq <- clmm(n.compounds.T ~ km30 + I(km30^2) +(1|Region) + (1|WMU), data=pctAG1)
ag60 <- clmm(n.compounds.T ~ km60 + (1|Region) + (1|WMU), data=pctAG1)
ag60sq <- clmm(n.compounds.T ~ km60 + I(km60^2) + (1|Region) + (1|WMU), data=pctAG1)

pctAG_sel <- model.sel(ag15, ag30, ag60, ag15sq, ag30sq, ag60sq)
pctAG_sel

## Beech basal area
baa1 <- dat[, c(1:17, 20, 26)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=baa, values_fn=unique) %>% as.data.frame()
names(baa1)[18:20] <- c("km15", "km30", "km60") 

## Scale and center variables
baa1[,c(18:20)] <- scale(baa1[,c(18:20)])

baa15 <- clmm(n.compounds.T ~ km15 + (1|Region) + (1|WMU), data=baa1)
baa15sq <- clmm(n.compounds.T ~ km15 + I(km15^2) + (1|Region) + (1|WMU), data=baa1)
baa30 <- clmm(n.compounds.T ~ km30 + (1|Region) + (1|WMU), data=baa1)
baa30sq <- clmm(n.compounds.T ~ km30 + I(km30^2) + (1|Region) + (1|WMU), data=baa1)
baa60 <- clmm(n.compounds.T ~ km60 + (1|Region) + (1|WMU), data=baa1)
baa60sq <- clmm(n.compounds.T ~ km60 + I(km60^2) + (1|Region) + (1|WMU), data=baa1)

baa_sel <- model.sel(baa15, baa30, baa60, baa15sq, baa30sq, baa60sq)
baa_sel


## Intermix WUI
intermix1 <- dat[, c(1:18, 21, 26)]
intermix1 <- unite(intermix1, "buffrad", 17:18, sep="_")
intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[18:26] <- c("b15km2_100m", "b30km2_100m", "b60km2_100m",
                             "b15km2_250m", "b30km2_250m", "b60km2_250m",
                             "b15km2_500m", "b30km2_500m", "b60km2_500m") 

## Scale and center variables
intermix1[,c(18:26)] <- scale(intermix1[,c(18:26)])

mix_15100 <- clmm(n.compounds.T ~ b15km2_100m + (1|Region) + (1|WMU), data=intermix1)
mix_15250 <- clmm(n.compounds.T ~ b15km2_250m + (1|Region) + (1|WMU), data=intermix1)
mix_15500 <- clmm(n.compounds.T ~ b15km2_500m + (1|Region) + (1|WMU), data=intermix1)
mix_15100sq <- clmm(n.compounds.T ~ b15km2_100m + I(b15km2_100m^2) + (1|Region) + (1|WMU), data=intermix1)
mix_15250sq <- clmm(n.compounds.T ~ b15km2_250m + I(b15km2_250m^2) + (1|Region) + (1|WMU), data=intermix1)
mix_15500sq <- clmm(n.compounds.T ~ b15km2_500m + I(b15km2_500m^2) + (1|Region) + (1|WMU), data=intermix1)

mix_30100 <- clmm(n.compounds.T ~ b30km2_100m + (1|Region) + (1|WMU), data=intermix1)
mix_30250 <- clmm(n.compounds.T ~ b30km2_250m + (1|Region) + (1|WMU), data=intermix1)
mix_30500 <- clmm(n.compounds.T ~ b30km2_500m + (1|Region) + (1|WMU), data=intermix1)
mix_30100sq <- clmm(n.compounds.T ~ b30km2_100m + I(b30km2_100m^2) + (1|Region) + (1|WMU), data=intermix1)
mix_30250sq <- clmm(n.compounds.T ~ b30km2_250m + I(b30km2_250m^2) + (1|Region) + (1|WMU), data=intermix1)
mix_30500sq <- clmm(n.compounds.T ~ b30km2_500m + I(b30km2_500m^2) + (1|Region) + (1|WMU), data=intermix1)

mix_60100 <- clmm(n.compounds.T ~ b60km2_100m + (1|Region) + (1|WMU), data=intermix1)
mix_60250 <- clmm(n.compounds.T ~ b60km2_250m + (1|Region) + (1|WMU), data=intermix1)
mix_60500 <- clmm(n.compounds.T ~ b60km2_500m + (1|Region) + (1|WMU), data=intermix1)
mix_60100sq <- clmm(n.compounds.T ~ b60km2_100m + I(b60km2_100m^2) + (1|Region) + (1|WMU), data=intermix1)
mix_60250sq <- clmm(n.compounds.T ~ b60km2_250m + I(b60km2_250m^2) + (1|Region) + (1|WMU), data=intermix1)
mix_60500sq <- clmm(n.compounds.T ~ b60km2_500m + I(b60km2_500m^2) + (1|Region) + (1|WMU), data=intermix1)

intermix_sel <- model.sel(mix_15100, mix_30100, mix_60100, mix_15100sq, mix_30100sq, mix_60100sq,
                          mix_15250, mix_30250, mix_60250, mix_15250sq, mix_30250sq, mix_60250sq,
                          mix_15500, mix_30500, mix_60500, mix_15500sq, mix_30500sq, mix_60500sq)
intermix_sel

## Interface WUI
interface1 <- dat[, c(1:18, 22, 26)]
interface1 <- unite(interface1, "buffrad", 17:18, sep="_")
interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[18:26] <- c("b15km2_100m", "b30km2_100m", "b60km2_100m",
                              "b15km2_250m", "b30km2_250m", "b60km2_250m",
                              "b15km2_500m", "b30km2_500m", "b60km2_500m") 

## Scale and center variables
interface1[,c(18:26)] <- scale(interface1[,c(18:26)])

face_15100 <- clmm(n.compounds.T ~ b15km2_100m + (1|Region) + (1|WMU), data=interface1)
face_15250 <- clmm(n.compounds.T ~ b15km2_250m + (1|Region) + (1|WMU), data=interface1)
face_15500 <- clmm(n.compounds.T ~ b15km2_500m + (1|Region) + (1|WMU), data=interface1)
face_15100sq <- clmm(n.compounds.T ~ b15km2_100m + I(b15km2_100m^2) + (1|Region) + (1|WMU), data=interface1)
face_15250sq <- clmm(n.compounds.T ~ b15km2_250m + I(b15km2_250m^2) + (1|Region) + (1|WMU), data=interface1)
face_15500sq <- clmm(n.compounds.T ~ b15km2_500m + I(b15km2_500m^2) + (1|Region) + (1|WMU), data=interface1)

face_30100 <- clmm(n.compounds.T ~ b30km2_100m + (1|Region) + (1|WMU), data=interface1)
face_30250 <- clmm(n.compounds.T ~ b30km2_250m + (1|Region) + (1|WMU), data=interface1)
face_30500 <- clmm(n.compounds.T ~ b30km2_500m + (1|Region) + (1|WMU), data=interface1)
face_30100sq <- clmm(n.compounds.T ~ b30km2_100m + I(b30km2_100m^2) + (1|Region) + (1|WMU), data=interface1)
face_30250sq <- clmm(n.compounds.T ~ b30km2_250m + I(b30km2_250m^2) + (1|Region) + (1|WMU), data=interface1)
face_30500sq <- clmm(n.compounds.T ~ b30km2_500m + I(b30km2_500m^2) + (1|Region) + (1|WMU), data=interface1)

face_60100 <- clmm(n.compounds.T ~ b60km2_100m + (1|Region) + (1|WMU), data=interface1)
face_60250 <- clmm(n.compounds.T ~ b60km2_250m + (1|Region) + (1|WMU), data=interface1)
face_60500 <- clmm(n.compounds.T ~ b60km2_500m + (1|Region) + (1|WMU), data=interface1)
face_60100sq <- clmm(n.compounds.T ~ b60km2_100m + I(b60km2_100m^2) + (1|Region) + (1|WMU), data=interface1)
face_60250sq <- clmm(n.compounds.T ~ b60km2_250m + I(b60km2_250m^2) + (1|Region) + (1|WMU), data=interface1)
face_60500sq <- clmm(n.compounds.T ~ b60km2_500m + I(b60km2_500m^2) + (1|Region) + (1|WMU), data=interface1)

interface_sel <- model.sel(face_15100, face_30100, face_60100, face_15100sq, face_30100sq, face_60100sq,  
                           face_15250, face_30250, face_60250, face_15250sq, face_30250sq, face_60250sq,  
                           face_15500, face_30500, face_60500, face_15500sq, face_30500sq, face_60500sq)  
interface_sel

## Total WUI
wui1 <- dat[, c(1:18, 23, 26)]
wui1 <- unite(wui1, "buffrad", 17:18, sep="_")
wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[18:26] <- c("b15km2_100m", "b30km2_100m", "b60km2_100m",
                        "b15km2_250m", "b30km2_250m", "b60km2_250m",
                        "b15km2_500m", "b30km2_500m", "b60km2_500m") 

# Set up WMU as group number
wui1$WMUnum <- as.numeric(factor(wui1$WMU, labels=1:54))

## Scale and center variables
wui1[,c(18:26)] <- scale(wui1[,c(18:26)])

wui_15100 <- clmm2(n.compounds.T ~ b15km2_100m + (1|Region), data=wui1)
wui_15250 <- clmm2(n.compounds.T ~ b15km2_250m + (1|Region) + (1|Region:WMU), data=wui1)
wui_15500 <- clmm2(n.compounds.T ~ b15km2_500m + (1|Region) + (1|Region:WMU), data=wui1)
wui_15100sq <- clmm(n.compounds.T ~ b15km2_100m + I(b15km2_100m^2) + (1|Region) + (1|WMU), data=wui1)
wui_15250sq <- clmm(n.compounds.T ~ b15km2_250m + I(b15km2_250m^2) + (1|Region) + (1|WMU), data=wui1)
wui_15500sq <- clmm(n.compounds.T ~ b15km2_500m + I(b15km2_500m^2) + (1|Region) + (1|WMU), data=wui1)

wui_30100 <- clmm(n.compounds.T ~ b30km2_100m + (1|Region) + (1|WMU), data=wui1)
wui_30250 <- clmm(n.compounds.T ~ b30km2_250m + (1|Region) + (1|WMU), data=wui1)
wui_30500 <- clmm(n.compounds.T ~ b30km2_500m + (1|Region) + (1|WMU), data=wui1)
wui_30100sq <- clmm(n.compounds.T ~ b30km2_100m + I(b30km2_100m^2) + (1|Region) + (1|WMU), data=wui1)
wui_30250sq <- clmm(n.compounds.T ~ b30km2_250m + I(b30km2_250m^2) + (1|Region) + (1|WMU), data=wui1)
wui_30500sq <- clmm(n.compounds.T ~ b30km2_500m + I(b30km2_500m^2) + (1|Region) + (1|WMU), data=wui1)

wui_60100 <- clmm(n.compounds.T ~ b60km2_100m + (1|Region) + (1|WMU), data=wui1)
wui_60250 <- clmm(n.compounds.T ~ b60km2_250m + (1|Region) + (1|WMU), data=wui1)
wui_60500 <- clmm(n.compounds.T ~ b60km2_500m + (1|Region) + (1|WMU), data=wui1)
wui_60100sq <- clmm(n.compounds.T ~ b60km2_100m + I(b60km2_100m^2) + (1|Region) + (1|WMU), data=wui1)
wui_60250sq <- clmm(n.compounds.T ~ b60km2_250m + I(b60km2_250m^2) + (1|Region) + (1|WMU), data=wui1)
wui_60500sq <- clmm(n.compounds.T ~ b60km2_500m + I(b60km2_500m^2) + (1|Region) + (1|WMU), data=wui1)

wui_sel <- model.sel(wui_15100, wui_30100, wui_60100, wui_15100sq, wui_30100sq, wui_60100sq, 
                     wui_15250, wui_30250, wui_60250, wui_15250sq, wui_30250sq, wui_60250sq, 
                     wui_15500, wui_30500, wui_60500, wui_15500sq, wui_30500sq, wui_60500sq)
wui_sel

#### Set up data to run for each combination of covariates ####
dat1 <- dat[,c(1:15,23:25)]
dat1 <- distinct(dat1)

# Join percent agriculture
pctAG1 <- pctAG1[,c(1:3, 17:20)]
dat1 <- left_join(dat1, pctAG1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[19:21] <- c("pctAG_15", "pctAG_30", "pctAG_60")

# Join beech basal area
baa1 <- baa1[,c(1:3, 20)]
dat1 <- left_join(dat1, baa1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[22] <- "baa_60"

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 28)]
dat1 <- left_join(dat1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[23] <- "intermix_60km2_500m" 

# Join interface WUI
interface1 <- interface1[, c(1:3, 20)]
dat1 <- left_join(dat1, interface1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[24] <- "interface_60km2_100m" 

# Join total WUI
wui1 <- wui1[,c(1:3, 28)]
dat1 <- left_join(dat1, wui1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[26] <- "wui_60km2_500m" 

## Check correlation matrix
cor(dat1[,19:26])

# Scale age
dat1[,c()] <- scale(dat1[,c()])

## Different way of accounting for year/mast cycle
dat1$fBMI <- ordered(dat1$year, levels=c(2019, 2018, 2020))
dat1$year <- as.factor(dat1$year)

## Set up models

m1g <- clmm(n.compounds.T ~ Sex*Age +
                baa_60*laggedBMI + baa_60*fBMI * baa_60*year +
                pctAG_15 + pctAG_4_5 + pctAG_30 + pctAG_60 + 
                intermix_60km2_500m + interface_60km2_100m + wui_60km2_500m +
                (1|Region), data=pt_dat, na.action="na.fail")
  summary(m1g)
  
T_dredge <- dredge(m1g, subset=!(pctAG_15 && pctAG_4_5) &&
                                   !(pctAG_30 && pctAG_4_5) &&
                                   !(pctAG_60 && pctAG_4_5) &&
                                   !(pctAG_15 && pctAG_30) &&
                                   !(pctAG_15 && pctAG_60) &&
                                   !(pctAG_30 && pctAG_60) && 
                                   !(laggedBMI && fBMI) &&
                                   !(laggedBMI && year) &&
                                   !(fBMI && year) &&
                                   !(intermix_60km2_500m && interface_60km2_100m) &&
                                   !(intermix_60km2_500m && wui_60km2_500m) &&
                                   !(wui_60km2_500m && interface_60km2_100m))
  



# Save dredge tables
save(file="output/dredge_tables_ncompT.Rdata", list="pt_spec")

#### Do this whole thing again but for measured only ####

## Percent AG
pctAG1 <- dat[, c(1:17, 19, 26)]
pctAG1 <- pctAG1 %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=pct_ag, values_fn=unique) %>% as.data.frame()
names(pctAG1)[18:20] <- c("km15", "km30", "km60") 

## Scale and center variables
pctAG1[,c(18:20)] <- scale(pctAG1[,c(18:20)])

# Run models
ag15 <- clmm(n.compounds.MO ~ km15 + (1|Region) + (1|WMU), data=pctAG1)
ag15sq <- clmm(n.compounds.MO ~ km15 + I(km15^2) + (1|Region) + (1|WMU), data=pctAG1)
ag30 <- clmm(n.compounds.MO ~ km30 + (1|Region) + (1|WMU), data=pctAG1)
ag30sq <- clmm(n.compounds.MO ~ km30 + I(km30^2) +(1|Region) + (1|WMU), data=pctAG1)
ag60 <- clmm(n.compounds.MO ~ km60 + (1|Region) + (1|WMU), data=pctAG1)
ag60sq <- clmm(n.compounds.MO ~ km60 + I(km60^2) + (1|Region) + (1|WMU), data=pctAG1)

pctAG_sel <- model.sel(ag15, ag30, ag60, ag15sq, ag30sq, ag60sq)
pctAG_sel

## Beech basal area
baa1 <- dat[, c(1:17, 20, 26)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=baa, values_fn=unique) %>% as.data.frame()
names(baa1)[18:20] <- c("km15", "km30", "km60") 

## Scale and center variables
baa1[,c(18:20)] <- scale(baa1[,c(18:20)])

baa15 <- clmm(n.compounds.MO ~ km15 + (1|Region) + (1|WMU), data=baa1)
baa15sq <- clmm(n.compounds.MO ~ km15 + I(km15^2) + (1|Region) + (1|WMU), data=baa1)
baa30 <- clmm(n.compounds.MO ~ km30 + (1|Region) + (1|WMU), data=baa1)
baa30sq <- clmm(n.compounds.MO ~ km30 + I(km30^2) + (1|Region) + (1|WMU), data=baa1)
baa60 <- clmm(n.compounds.MO ~ km60 + (1|Region) + (1|WMU), data=baa1)
baa60sq <- clmm(n.compounds.MO ~ km60 + I(km60^2) + (1|Region) + (1|WMU), data=baa1)

baa_sel <- model.sel(baa15, baa30, baa60, baa15sq, baa30sq, baa60sq)
baa_sel


## Intermix WUI
intermix1 <- dat[, c(1:18, 21, 26)]
intermix1 <- unite(intermix1, "buffrad", 17:18, sep="_")
intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[18:26] <- c("b15km2_100m", "b30km2_100m", "b60km2_100m",
                             "b15km2_250m", "b30km2_250m", "b60km2_250m",
                             "b15km2_500m", "b30km2_500m", "b60km2_500m") 

## Scale and center variables
intermix1[,c(18:26)] <- scale(intermix1[,c(18:26)])

mix_15100 <- clmm(n.compounds.MO ~ b15km2_100m + (1|Region) + (1|WMU), data=intermix1)
mix_15250 <- clmm(n.compounds.MO ~ b15km2_250m + (1|Region) + (1|WMU), data=intermix1)
mix_15500 <- clmm(n.compounds.MO ~ b15km2_500m + (1|Region) + (1|WMU), data=intermix1)
mix_15100sq <- clmm(n.compounds.MO ~ b15km2_100m + I(b15km2_100m^2) + (1|Region) + (1|WMU), data=intermix1)
mix_15250sq <- clmm(n.compounds.MO ~ b15km2_250m + I(b15km2_250m^2) + (1|Region) + (1|WMU), data=intermix1)
mix_15500sq <- clmm(n.compounds.MO ~ b15km2_500m + I(b15km2_500m^2) + (1|Region) + (1|WMU), data=intermix1)

mix_30100 <- clmm(n.compounds.MO ~ b30km2_100m + (1|Region) + (1|WMU), data=intermix1)
mix_30250 <- clmm(n.compounds.MO ~ b30km2_250m + (1|Region) + (1|WMU), data=intermix1)
mix_30500 <- clmm(n.compounds.MO ~ b30km2_500m + (1|Region) + (1|WMU), data=intermix1)
mix_30100sq <- clmm(n.compounds.MO ~ b30km2_100m + I(b30km2_100m^2) + (1|Region) + (1|WMU), data=intermix1)
mix_30250sq <- clmm(n.compounds.MO ~ b30km2_250m + I(b30km2_250m^2) + (1|Region) + (1|WMU), data=intermix1)
mix_30500sq <- clmm(n.compounds.MO ~ b30km2_500m + I(b30km2_500m^2) + (1|Region) + (1|WMU), data=intermix1)

mix_60100 <- clmm(n.compounds.MO ~ b60km2_100m + (1|Region) + (1|WMU), data=intermix1)
mix_60250 <- clmm(n.compounds.MO ~ b60km2_250m + (1|Region) + (1|WMU), data=intermix1)
mix_60500 <- clmm(n.compounds.MO ~ b60km2_500m + (1|Region) + (1|WMU), data=intermix1)
mix_60100sq <- clmm(n.compounds.MO ~ b60km2_100m + I(b60km2_100m^2) + (1|Region) + (1|WMU), data=intermix1)
mix_60250sq <- clmm(n.compounds.MO ~ b60km2_250m + I(b60km2_250m^2) + (1|Region) + (1|WMU), data=intermix1)
mix_60500sq <- clmm(n.compounds.MO ~ b60km2_500m + I(b60km2_500m^2) + (1|Region) + (1|WMU), data=intermix1)

intermix_sel <- model.sel(mix_15100, mix_30100, mix_60100, mix_15100sq, mix_30100sq, mix_60100sq,
                          mix_15250, mix_30250, mix_60250, mix_15250sq, mix_30250sq, mix_60250sq,
                          mix_15500, mix_30500, mix_60500, mix_15500sq, mix_30500sq, mix_60500sq)
intermix_sel

## Interface WUI
interface1 <- dat[, c(1:18, 22, 26)]
interface1 <- unite(interface1, "buffrad", 17:18, sep="_")
interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[18:26] <- c("b15km2_100m", "b30km2_100m", "b60km2_100m",
                              "b15km2_250m", "b30km2_250m", "b60km2_250m",
                              "b15km2_500m", "b30km2_500m", "b60km2_500m") 

## Scale and center variables
interface1[,c(18:26)] <- scale(interface1[,c(18:26)])

face_15100 <- clmm(n.compounds.MO ~ b15km2_100m + (1|Region) + (1|WMU), data=interface1)
face_15250 <- clmm(n.compounds.MO ~ b15km2_250m + (1|Region) + (1|WMU), data=interface1)
face_15500 <- clmm(n.compounds.MO ~ b15km2_500m + (1|Region) + (1|WMU), data=interface1)
face_15100sq <- clmm(n.compounds.MO ~ b15km2_100m + I(b15km2_100m^2) + (1|Region) + (1|WMU), data=interface1)
face_15250sq <- clmm(n.compounds.MO ~ b15km2_250m + I(b15km2_250m^2) + (1|Region) + (1|WMU), data=interface1)
face_15500sq <- clmm(n.compounds.MO ~ b15km2_500m + I(b15km2_500m^2) + (1|Region) + (1|WMU), data=interface1)

face_30100 <- clmm(n.compounds.MO ~ b30km2_100m + (1|Region) + (1|WMU), data=interface1)
face_30250 <- clmm(n.compounds.MO ~ b30km2_250m + (1|Region) + (1|WMU), data=interface1)
face_30500 <- clmm(n.compounds.MO ~ b30km2_500m + (1|Region) + (1|WMU), data=interface1)
face_30100sq <- clmm(n.compounds.MO ~ b30km2_100m + I(b30km2_100m^2) + (1|Region) + (1|WMU), data=interface1)
face_30250sq <- clmm(n.compounds.MO ~ b30km2_250m + I(b30km2_250m^2) + (1|Region) + (1|WMU), data=interface1)
face_30500sq <- clmm(n.compounds.MO ~ b30km2_500m + I(b30km2_500m^2) + (1|Region) + (1|WMU), data=interface1)

face_60100 <- clmm(n.compounds.MO ~ b60km2_100m + (1|Region) + (1|WMU), data=interface1)
face_60250 <- clmm(n.compounds.MO ~ b60km2_250m + (1|Region) + (1|WMU), data=interface1)
face_60500 <- clmm(n.compounds.MO ~ b60km2_500m + (1|Region) + (1|WMU), data=interface1)
face_60100sq <- clmm(n.compounds.MO ~ b60km2_100m + I(b60km2_100m^2) + (1|Region) + (1|WMU), data=interface1)
face_60250sq <- clmm(n.compounds.MO ~ b60km2_250m + I(b60km2_250m^2) + (1|Region) + (1|WMU), data=interface1)
face_60500sq <- clmm(n.compounds.MO ~ b60km2_500m + I(b60km2_500m^2) + (1|Region) + (1|WMU), data=interface1)

interface_sel <- model.sel(face_15100, face_30100, face_60100, face_15100sq, face_30100sq, face_60100sq,  
                           face_15250, face_30250, face_60250, face_15250sq, face_30250sq, face_60250sq,  
                           face_15500, face_30500, face_60500, face_15500sq, face_30500sq, face_60500sq)  
interface_sel

## Total WUI
wui1 <- dat[, c(1:18, 23, 26)]
wui1 <- unite(wui1, "buffrad", 17:18, sep="_")
wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[18:26] <- c("b15km2_100m", "b30km2_100m", "b60km2_100m",
                        "b15km2_250m", "b30km2_250m", "b60km2_250m",
                        "b15km2_500m", "b30km2_500m", "b60km2_500m") 

# Set up WMU as group number
wui1$WMUnum <- as.numeric(factor(wui1$WMU, labels=1:54))

## Scale and center variables
wui1[,c(18:26)] <- scale(wui1[,c(18:26)])

wui_15100 <- clmm2(n.compounds.MO ~ b15km2_100m + (1|Region) + (1|WMU), data=wui1)
wui_15250 <- clmm2(n.compounds.MO ~ b15km2_250m + (1|Region) + (1|WMU), data=wui1)
wui_15500 <- clmm2(n.compounds.MO ~ b15km2_500m + (1|Region) + (1|WMU), data=wui1)
wui_15100sq <- clmm(n.compounds.MO ~ b15km2_100m + I(b15km2_100m^2) + (1|Region) + (1|WMU), data=wui1)
wui_15250sq <- clmm(n.compounds.MO ~ b15km2_250m + I(b15km2_250m^2) + (1|Region) + (1|WMU), data=wui1)
wui_15500sq <- clmm(n.compounds.MO ~ b15km2_500m + I(b15km2_500m^2) + (1|Region) + (1|WMU), data=wui1)

wui_30100 <- clmm(n.compounds.MO ~ b30km2_100m + (1|Region) + (1|WMU), data=wui1)
wui_30250 <- clmm(n.compounds.MO ~ b30km2_250m + (1|Region) + (1|WMU), data=wui1)
wui_30500 <- clmm(n.compounds.MO ~ b30km2_500m + (1|Region) + (1|WMU), data=wui1)
wui_30100sq <- clmm(n.compounds.MO ~ b30km2_100m + I(b30km2_100m^2) + (1|Region) + (1|WMU), data=wui1)
wui_30250sq <- clmm(n.compounds.MO ~ b30km2_250m + I(b30km2_250m^2) + (1|Region) + (1|WMU), data=wui1)
wui_30500sq <- clmm(n.compounds.MO ~ b30km2_500m + I(b30km2_500m^2) + (1|Region) + (1|WMU), data=wui1)

wui_60100 <- clmm(n.compounds.MO ~ b60km2_100m + (1|Region) + (1|WMU), data=wui1)
wui_60250 <- clmm(n.compounds.MO ~ b60km2_250m + (1|Region) + (1|WMU), data=wui1)
wui_60500 <- clmm(n.compounds.MO ~ b60km2_500m + (1|Region) + (1|WMU), data=wui1)
wui_60100sq <- clmm(n.compounds.MO ~ b60km2_100m + I(b60km2_100m^2) + (1|Region) + (1|WMU), data=wui1)
wui_60250sq <- clmm(n.compounds.MO ~ b60km2_250m + I(b60km2_250m^2) + (1|Region) + (1|WMU), data=wui1)
wui_60500sq <- clmm(n.compounds.MO ~ b60km2_500m + I(b60km2_500m^2) + (1|Region) + (1|WMU), data=wui1)

wui_sel <- model.sel(wui_15100, wui_30100, wui_60100, wui_15100sq, wui_30100sq, wui_60100sq, 
                     wui_15250, wui_30250, wui_60250, wui_15250sq, wui_30250sq, wui_60250sq, 
                     wui_15500, wui_30500, wui_60500, wui_15500sq, wui_30500sq, wui_60500sq)
wui_sel

#### Set up data to run for each combination of covariates ####
dat1 <- dat[,c(1:15,23:25)]
dat1 <- distinct(dat1)

# Join percent agriculture
pctAG1 <- pctAG1[,c(1:3, 20)]
dat1 <- left_join(dat1, pctAG1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[19] <-  "pctAG_60"
  
# beech basal area
baa1 <- baa1[,c(1:3, 17:20)]
dat1 <- left_join(dat1, baa1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[20:23] <- c("baa_4", "baa_15", "baa_30", "baa_60")

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 24, 28)]
dat1 <- left_join(dat1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[24:25] <- c("intermix_60km2_250m", "intermix_60km2_500m")

# Join interface WUI
interface1 <- interface1[, c(1:3, 20, 25)]
dat1 <- left_join(dat1, interface1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[26:27] <- c("interface_4_5km2_500m", "interface_60km2_100m")

# Join total WUI
wui1 <- wui1[,c(1:3, 28)]
dat1 <- left_join(dat1, wui1, by=c("RegionalID", "pt_name", "pt_index"))
names(dat1)[28] <- "wui_60km2_500m" 

## Check correlation matrix
cor(dat1[,19:28])

## Different way of accounting for year/mast cycle
dat1$fBMI <- ordered(dat1$year, levels=c(2019, 2018, 2020))
dat1$year <- as.factor(dat1$year)

## Set up models

pt_spec <- list()

for (i in 1:10) {
  
  pt_dat <- dat1[dat1$pt_index==i,]
  
  m1g <- clmm(n.compounds.T ~ Sex*Age + pctAG_60 + 
                baa_15*laggedBMI + baa_4*laggedBMI + baa_30*laggedBMI +baa_60*laggedBMI +
                baa_15*year + baa_4*year + baa_30*year +baa_60*year +
                baa_15*fBMI + baa_4*fBMI + baa_30*fBMI +baa_60*fBMI +
                intermix_60km2_500m + interface_60km2_100m + wui_60km2_500m +
                intermix_60km2_250m + interface_4_5km2_500m +
                (1|Region), data=pt_dat, na.action="na.fail")
  summary(m1g)
  
  pt_dredge <- dredge(m1g, subset=!(baa_4 && baa_30) &&
                                  !(baa_4 && baa_15) &&
                                  !(baa_4 && baa_60) &&
                                  !(baa_15 && baa_30) &&
                                  !(baa_15 && baa_60) &&
                                  !(baa_30 && baa_60) && 
                                  !(laggedBMI && fBMI) &&
                                  !(laggedBMI && year) &&
                                  !(fBMI && year) &&
                                  !(interface_60km2_100m && interface_4_5km2_500m) &&
                                  !(intermix_60km2_500m && intermix_60km2_250m) &&
                                  !(intermix_60km2_500m && wui_60km2_500m) &&
                                  !(intermix_60km2_500m && interface_60km2_100m) &&
                                  !(intermix_60km2_500m && interface_4_5km2_500m) &&
                                  !(intermix_60km2_250m && wui_60km2_500m) &&
                                  !(intermix_60km2_250m && interface_60km2_100m) &&
                                  !(intermix_60km2_250m && interface_4_5km2_500m) &&                        
                                  !(interface_60km2_100m && wui_60km2_500m) &&
                                  !(interface_4_5km2_500m && wui_60km2_500m))
  
  pt_spec[[i]] <- pt_dredge
  
}

# Save dredge tables
save(file="output/dredge_tables_ncompMO.Rdata", list="pt_spec")

