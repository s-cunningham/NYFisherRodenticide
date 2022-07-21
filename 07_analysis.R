library(tidyverse)
library(brms)

set.seed(123)

#### Read in data ####
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

#### Analysis Set-up ####

## set up binary variable
dat$binary.T <- ifelse(dat$n.compounds.T==0, 0, 1)
dat$binary.MO <- ifelse(dat$n.compounds.MO==0, 0, 1)

## Scale and center variables
dat[,c(9,18:23)] <- scale(dat[,c(9,18:23)])
dat$fBMI <- factor(dat$year, levels=c(2019, 2018, 2020), labels=c("failure", "intermediate", "high"))

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

#### Test covariates individually ####

## Run test model

pctAG15 <- brm(binary.T ~ pct_ag + (1|WMUA_code), family=bernoulli(link = "logit"), 
               data=dat_r100_b15, chains=3, iter=50000, backend="cmdstanr", cores=3)

# Look at priors
prior_summary(pctAG15)

# Try different priors

new_priors <- c(prior_string("normal(0,10)", class="b"),
                prior_string("normal(0,10)", class="Intercept"))


pctAG15 <- brm(binary.T ~ pct_ag + (1|Region), family=bernoulli(link = "logit"), 
               data=dat_r100_b15, chains=3, iter=50000, backend="cmdstanr", cores=3,
               prior=new_priors)





## Run models

m1_r100_b4p5 <- brm(binary.T ~ baa*fBMI + Sex*Age + (1|WMUA_code), family=bernoulli(link = "logit"), data=dat_r100_b15, chains=3, iter=50000, backend="cmdstanr", cores=3)
m1_r100_b30 <- brm(binary.T ~ pct_ag + intermix + (1|Region), family=bernoulli(link = "logit"), data=dat_r100_b30, chains=3, iter=50000, backend="cmdstanr", cores=3)

summary(m1_r100_b4p5)
summary(m1_r100_b30)

m1_r100_b15 <- add_criterion(m1_r100_b15, c("loo", "waic"))
m1_r100_b30 <- add_criterion(m1_r100_b30, c("loo", "waic"))


loo_compare(m1_r100_b15, m1_r100_b30, criterion="waic")






