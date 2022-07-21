library(tidyverse)
library(brms)

set.seed(123)

#### Read in data ####
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

#### Analysis Set-up ####

## Scale and center variables
dat[,c(9,18:23)] <- scale(dat[,c(9,18:23)])
dat$fBMI <- factor(dat$year, levels=c(2019, 2018, 2020), labels=c("failure", "intermediate", "high"))

cor(dat[,18:23])

#### Subset data by buffers and radii ####

dat_r100_b4p5 <- dat[dat$buffsize==4.5 & dat$radius==100 & dat$pt_index==1,]
cor(dat_r100_b4p5[,18:23])
dat_r100_b15 <- dat[dat$buffsize==15 & dat$radius==100 & dat$pt_index==1,]
cor(dat_r100_b15[,18:23])
dat_r100_b30 <- dat[dat$buffsize==30 & dat$radius==100 & dat$pt_index==1,]
cor(dat_r100_b30[,18:23])
dat_r100_b60 <- dat[dat$buffsize==60 & dat$radius==100 & dat$pt_index==1,]
cor(dat_r100_b60[,18:23])

dat_r250_b4p5 <- dat[dat$buffsize==4.5 & dat$radius==250 & dat$pt_index==1,]
cor(dat_r250_b4p5[,18:23])
dat_r250_b15 <- dat[dat$buffsize==15 & dat$radius==250 & dat$pt_index==1,]
cor(dat_r250_b15[,18:23])
dat_r250_b30 <- dat[dat$buffsize==30 & dat$radius==250 & dat$pt_index==1,]
cor(dat_r250_b30[,18:23])
dat_r250_b60 <- dat[dat$buffsize==60 & dat$radius==250 & dat$pt_index==1,]
cor(dat_r250_b60[,18:23])

dat_r500_b4p5 <- dat[dat$buffsize==4.5 & dat$radius==500 & dat$pt_index==1,]
cor(dat_r500_b4p5[,18:23])
dat_r500_b15 <- dat[dat$buffsize==15 & dat$radius==500 & dat$pt_index==1,]
cor(dat_r500_b15[,18:23])
dat_r500_b30 <- dat[dat$buffsize==30 & dat$radius==500 & dat$pt_index==1,]
cor(dat_r500_b30[,18:23])
dat_r500_b60 <- dat[dat$buffsize==60 & dat$radius==500 & dat$pt_index==1,]
cor(dat_r500_b60[,18:23])