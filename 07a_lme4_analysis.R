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

## Create lists to store results
pctAG <- list()
baa <- list()
intermix <- list()
interface <- list()
wui <- list()

## Look over each point set
for (i in 1:10) {
  
  ## Percent AG
  pctAG1 <- dat[dat$pt_index==i, c(1:16, 18, 26, 27)]
  pctAG1 <- pctAG1 %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=pct_ag, values_fn=unique) %>% as.data.frame()
  names(pctAG1)[18:21] <- c("4_5km2", "15km2", "30km2", "60km2") 
  
  pctAG_45 <- glmer(binary.T ~ `4_5km2` + Age*Sex + (1|Region), family=binomial(link="logit"), data=pctAG1)
  pctAG_15 <- glmer(binary.T ~ `15km2` + Age*Sex + (1|Region), family=binomial(link="logit"), data=pctAG1)
  pctAG_30 <- glmer(binary.T ~ `30km2` + Age*Sex + (1|Region), family=binomial(link="logit"), data=pctAG1)
  pctAG_60 <- glmer(binary.T ~ `60km2` + Age*Sex + (1|Region), family=binomial(link="logit"), data=pctAG1)
  
  pctAG[[i]] <- model.sel(pctAG_45, pctAG_15, pctAG_30, pctAG_60)
  
  ## Beech basal area
  baa1 <- dat[dat$pt_index==i, c(1:16, 19, 26, 27)]
  baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=baa, values_fn=unique) %>% as.data.frame()
  names(baa1)[18:21] <- c("4_5km2", "15km2", "30km2", "60km2") 
  
  baa_45 <- glmer(binary.T ~ `4_5km2` + Age*Sex +(1|Region), family=binomial(link="logit"), data=baa1)
  baa_15 <- glmer(binary.T ~ `15km2` + Age*Sex + (1|Region), family=binomial(link="logit"), data=baa1)
  baa_30 <- glmer(binary.T ~ `30km2` + Age*Sex + (1|Region), family=binomial(link="logit"), data=baa1)
  baa_60 <- glmer(binary.T ~ `60km2` + Age*Sex + (1|Region), family=binomial(link="logit"), data=baa1)
  
  # Model selection with MuMIn
  baa[[i]] <- model.sel(baa_45, baa_15, baa_30, baa_60)
  
  ## Intermix WUI
  intermix1 <- dat[dat$pt_index==i, c(1:17, 20, 26, 27)]
  intermix1 <- unite(intermix1, "buffrad", 16:17, sep="_")
  intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
  names(intermix1)[18:29] <- c("4_5km2_100m", "15km2_100m", "30km2_100m", "60km2_100m",
                               "4_5km2_250m", "15km2_250m", "30km2_250m", "60km2_250m",
                               "4_5km2_500m", "15km2_500m", "30km2_500m", "60km2_500m") 
  
  mix_45100 <- glmer(binary.T ~ `4_5km2_100m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=intermix1)
  mix_45250 <- glmer(binary.T ~ `4_5km2_250m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=intermix1)
  mix_45500 <- glmer(binary.T ~ `4_5km2_500m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=intermix1)
  
  mix_15100 <- glmer(binary.T ~ `15km2_100m` +Age*Sex + (1|Region), family=binomial(link="logit"), data=intermix1)
  mix_15250 <- glmer(binary.T ~ `15km2_250m` +Age*Sex + (1|Region), family=binomial(link="logit"), data=intermix1)
  mix_15500 <- glmer(binary.T ~ `15km2_500m` +Age*Sex + (1|Region), family=binomial(link="logit"), data=intermix1)
  
  mix_30100 <- glmer(binary.T ~ `30km2_100m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=intermix1)
  mix_30250 <- glmer(binary.T ~ `30km2_250m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=intermix1)
  mix_30500 <- glmer(binary.T ~ `30km2_500m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=intermix1)
  
  mix_60100 <- glmer(binary.T ~ `60km2_100m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=intermix1)
  mix_60250 <- glmer(binary.T ~ `60km2_250m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=intermix1)
  mix_60500 <- glmer(binary.T ~ `60km2_500m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=intermix1)
  
  intermix[[i]] <- model.sel(mix_45100, mix_15100, mix_30100,mix_60100, 
                             mix_45250, mix_15250, mix_30250,mix_60250,  
                             mix_45500, mix_15500, mix_30500,mix_60500)
  
  ## Interface WUI
  interface1 <- dat[dat$pt_index==i, c(1:17, 21, 26, 27)]
  interface1 <- unite(interface1, "buffrad", 16:17, sep="_")
  interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
  names(interface1)[18:29] <- c("4_5km2_100m", "15km2_100m", "30km2_100m", "60km2_100m",
                                "4_5km2_250m", "15km2_250m", "30km2_250m", "60km2_250m",
                                "4_5km2_500m", "15km2_500m", "30km2_500m", "60km2_500m") 
  
  face_45100 <- glmer(binary.T ~ `4_5km2_100m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=interface1)
  face_45250 <- glmer(binary.T ~ `4_5km2_250m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=interface1)
  face_45500 <- glmer(binary.T ~ `4_5km2_500m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=interface1)
  
  face_15100 <- glmer(binary.T ~ `15km2_100m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=interface1)
  face_15250 <- glmer(binary.T ~ `15km2_250m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=interface1)
  face_15500 <- glmer(binary.T ~ `15km2_500m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=interface1)
  
  face_30100 <- glmer(binary.T ~ `30km2_100m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=interface1)
  face_30250 <- glmer(binary.T ~ `30km2_250m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=interface1)
  face_30500 <- glmer(binary.T ~ `30km2_500m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=interface1)
  
  face_60100 <- glmer(binary.T ~ `60km2_100m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=interface1)
  face_60250 <- glmer(binary.T ~ `60km2_250m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=interface1)
  face_60500 <- glmer(binary.T ~ `60km2_500m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=interface1)
  
  interface[[i]] <- model.sel(face_45100, face_15100, face_30100, face_60100, 
                              face_45250, face_15250, face_30250, face_60250,  
                              face_45500, face_15500, face_30500, face_60500)
  
  ## Total WUI
  wui1 <- dat[dat$pt_index==i, c(1:17, 22, 26, 27)]
  wui1 <- unite(wui1, "buffrad", 16:17, sep="_")
  wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
  names(wui1)[18:29] <- c("4_5km2_100m", "15km2_100m", "30km2_100m", "60km2_100m",
                          "4_5km2_250m", "15km2_250m", "30km2_250m", "60km2_250m",
                          "4_5km2_500m", "15km2_500m", "30km2_500m", "60km2_500m") 
  
  wui_45100 <- glmer(binary.T ~ `4_5km2_100m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=wui1)
  wui_45250 <- glmer(binary.T ~ `4_5km2_250m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=wui1)
  wui_45500 <- glmer(binary.T ~ `4_5km2_500m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=wui1)
  
  wui_15100 <- glmer(binary.T ~ `15km2_100m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=wui1)
  wui_15250 <- glmer(binary.T ~ `15km2_250m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=wui1)
  wui_15500 <- glmer(binary.T ~ `15km2_500m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=wui1)
  
  wui_30100 <- glmer(binary.T ~ `30km2_100m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=wui1)
  wui_30250 <- glmer(binary.T ~ `30km2_250m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=wui1)
  wui_30500 <- glmer(binary.T ~ `30km2_500m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=wui1)
  
  wui_60100 <- glmer(binary.T ~ `60km2_100m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=wui1)
  wui_60250 <- glmer(binary.T ~ `60km2_250m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=wui1)
  wui_60500 <- glmer(binary.T ~ `60km2_500m` + Age*Sex + (1|Region), family=binomial(link="logit"), data=wui1)
  
  wui[[i]] <- model.sel(wui_45100, wui_15100, wui_30100, wui_60100, 
                        wui_45250, wui_15250, wui_30250, wui_60250,  
                        wui_45500, wui_15500, wui_30500, wui_60500)
  
}




#### Set up data to have covariates and scales of interest ####



