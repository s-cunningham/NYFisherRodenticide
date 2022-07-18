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
dat$binary.T <- ifelse(dat$n.compounds.T==0, 0, 1)
dat$binary.MO <- ifelse(dat$n.compounds.MO==0, 0, 1)

## Scale and center variables
dat[,18:23] <- scale(dat[,18:23])
dat$fBMI <- ordered(dat$year, levels=c(2019, 2018, 2020))

cor(dat[,18:23])

## Subset data by buffers and radii
dat_r100_b15 <- dat[dat$buffsize==15 & dat$radius==100 & dat$pt_index==1,]
cor(dat_r100_b15[,18:23])
dat_r100_b30 <- dat[dat$buffsize==30 & dat$radius==100 & dat$pt_index==1,]
cor(dat_r100_b30[,18:23])
dat_r100_b60 <- dat[dat$buffsize==60 & dat$radius==100 & dat$pt_index==1,]
cor(dat_r100_b60[,18:23])

dat_r250_b15 <- dat[dat$buffsize==15 & dat$radius==250 & dat$pt_index==1,]
cor(dat_r250_b15[,18:23])
dat_r250_b30 <- dat[dat$buffsize==30 & dat$radius==250 & dat$pt_index==1,]
cor(dat_r250_b30[,18:23])
dat_r250_b60 <- dat[dat$buffsize==60 & dat$radius==250 & dat$pt_index==1,]
cor(dat_r250_b60[,18:23])

dat_r500_b15 <- dat[dat$buffsize==15 & dat$radius==500 & dat$pt_index==1,]
cor(dat_r500_b15[,18:23])
dat_r500_b30 <- dat[dat$buffsize==30 & dat$radius==500 & dat$pt_index==1,]
cor(dat_r500_b30[,18:23])
dat_r500_b60 <- dat[dat$buffsize==60 & dat$radius==500 & dat$pt_index==1,]
cor(dat_r500_b60[,18:23])


## Run models

m1_r100_b15 <- glmer(binary.T ~ pct_ag + totalWUI + (1|WMUA_code), family=binomial(link="logit"), data=dat_r100_b15)
m1_r100_b30 <- glmer(binary.T ~ pct_ag + totalWUI + (1|WMUA_code), family=binomial(link="logit"), data=dat_r100_b30)
m1_r100_b60 <- glmer(binary.T ~ pct_ag + totalWUI + (1|WMUA_code), family=binomial(link="logit"), data=dat_r100_b60)

m1_r250_b15 <- glmer(binary.T ~ pct_ag + totalWUI + (1|WMUA_code), family=binomial(link="logit"), data=dat_r250_b15)
m1_r250_b30 <- glmer(binary.T ~ pct_ag + totalWUI + (1|WMUA_code), family=binomial(link="logit"), data=dat_r250_b30)
m1_r250_b60 <- glmer(binary.T ~ pct_ag + totalWUI + (1|WMUA_code), family=binomial(link="logit"), data=dat_r250_b60)

summary(m1_r100_b15)
summary(m1_r100_b30)
summary(m1_r100_b60)

summary(m1_r250_b15)
summary(m1_r250_b30)
summary(m1_r250_b60)








