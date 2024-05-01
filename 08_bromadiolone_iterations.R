library(tidyverse)
library(nimble)
library(MCMCvis)
library(HDInterval)

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-edge_density_15, -edge_density_30, -edge_density_45, -build_cat_15, -build_cat_30,-build_cat_45) %>%
  mutate(mast_year=if_else(year==2019, 2, 1), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat[,c(8,16:42)] <- scale(dat[,c(8,16:42)])

# read data for individual compounds
brom <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Bromadiolone") %>%
  select(RegionalID,bin.exp)
dat <- left_join(dat, brom, by="RegionalID")
dat <- dat %>% select(RegionalID:Town,bin.exp,deciduous_15:mast_year)

# Build model in BUGS language
bromadiolone_code <- nimbleCode({
  
  ## Priors
  # beta coefficient priors
  beta_age ~ dnorm(0, sd=1.4)
  beta_age2 ~ dnorm(0, sd=1.4)
  for (k in 1:2) {
    beta_sex[k] ~ dnorm(0, sd=1.4)
  }
  beta_mast ~ dnorm(0, sd=1.4)
  beta_build ~ dnorm(0, sd=1.4)
  beta_decid ~ dnorm(0, sd=1.4)
  beta_evrgrn ~ dnorm(0, sd=1.4)
  beta_standm ~ dnorm(0, sd=1.4)
  beta_standsd ~ dnorm(0, sd=1.4)
  
  ## random intercepts
  # WMU
  for (k in 1:nWMU) {
    alpha[k] ~ dnorm(mu.alpha, sd=sigma.alpha)
  }
  mu.alpha ~ dnorm(0, 0.001)
  sigma.alpha ~ dunif(0, 100)
  
  ## Likelihood
  for (i in 1:N) {
    
    logit(p[i]) <- alpha[WMU[i]] + beta_age*age[i] + beta_age2*age2[i] + beta_sex[sex[i]] + 
      beta_mast*covars[i,1] + beta_decid*covars[i,2] + beta_evrgrn*covars[i,3] +
      beta_build*covars[i,4] + beta_standm*covars[i,5] + beta_standsd*covars[i,6]
    
    y[i] ~ dbern(p[i])
    
  }
  
})

# parameters to monitor
params <- c("beta_age","beta_age2","beta_sex","beta_mast","beta_decid","beta_evrgrn",
            "beta_build","beta_standm", "beta_standsd", "alpha", "mu.alpha", "sigma.alpha")  

# MCMC options
nt <- 1
ni <- 75000
nb <- 35000
nc <- 3

set.seed(1)
Inits <- list(sigma.alpha=1, mu.alpha=1,
              beta_mast=rnorm(1), beta_decid=rnorm(1), beta_evrgrn=rnorm(1), 
              beta_build=rnorm(1), beta_standm=rnorm(1), beta_standsd=rnorm(1), 
              beta_age=rnorm(1), beta_age2=rnorm(1), beta_sex=rnorm(2)) 

#### Loop over random locations ####
## Iteration 1
dat1 <- dat %>% filter(pt_index==1)

## Set up data
brom1 <- dat1$bin.exp
wmu1 <- as.numeric(factor(dat1$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars1 <- matrix(NA, nrow=nrow(dat1),ncol=6)
covars1[1:nrow(dat1),1] <- dat1$lag_beechnuts
covars1[1:nrow(dat1),2] <- dat1$deciduous_45
covars1[1:nrow(dat1),3] <- dat1$evergreen_45
covars1[1:nrow(dat1),4] <- dat1$nbuildings_45
covars1[1:nrow(dat1),5] <- dat1$stand_age_mean_45
covars1[1:nrow(dat1),6] <- dat1$stand_age_sd_45

## prep fof nimble model
Constants1 <- list(N=nrow(dat1),
                   sex=dat1$Sex,
                   WMU=wmu1, # random intercept
                   nWMU=length(unique(dat1$WMU)))

DataBundle1 <- list(y=brom1, # response
                    covars=covars1, # covariates 
                    age=dat1$Age,
                    age2=dat1$age2) 

set.seed(1)
brom.out1 <- nimbleMCMC(code=bromadiolone_code, constants=Constants1, data=DataBundle1, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum1 <- MCMCsummary(brom.out1)
brom.sum1 <- rownames_to_column(brom.sum1, "parameter")
range(brom.sum1$Rhat)

## Iteration 2
dat2 <- dat %>% filter(pt_index==2)

brom2 <- dat2$bin.exp
wmu2 <- as.numeric(factor(dat2$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars2 <- matrix(NA, nrow=nrow(dat2),ncol=6)
covars2[1:nrow(dat2),1] <- dat2$lag_beechnuts
covars2[1:nrow(dat2),2] <- dat2$deciduous_45
covars2[1:nrow(dat2),3] <- dat2$evergreen_45
covars2[1:nrow(dat2),4] <- dat2$nbuildings_45
covars2[1:nrow(dat2),5] <- dat2$stand_age_mean_45
covars2[1:nrow(dat2),6] <- dat2$stand_age_sd_45

## prep fof nimble model
Constants2 <- list(N=nrow(dat2),
                   sex=dat2$Sex,
                   WMU=wmu2, # random intercept
                   nWMU=length(unique(dat2$WMU)))

DataBundle2 <- list(y=brom2, # response
                    covars=covars2, # covariates 
                    age=dat2$Age,
                    age2=dat2$age2) 

set.seed(1)
brom.out2 <- nimbleMCMC(code=bromadiolone_code, constants=Constants2, data=DataBundle2, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum2 <- MCMCsummary(brom.out2)
brom.sum2 <- rownames_to_column(brom.sum2, "parameter")
range(brom.sum2$Rhat)

## Iteration 3
dat3 <- dat %>% filter(pt_index==3)

brom3 <- dat3$bin.exp
wmu3 <- as.numeric(factor(dat3$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars3 <- matrix(NA, nrow=nrow(dat3),ncol=6)
covars3[1:nrow(dat3),1] <- dat3$lag_beechnuts
covars3[1:nrow(dat3),2] <- dat3$deciduous_45
covars3[1:nrow(dat3),3] <- dat3$evergreen_45
covars3[1:nrow(dat3),4] <- dat3$nbuildings_45
covars3[1:nrow(dat3),5] <- dat3$stand_age_mean_45
covars3[1:nrow(dat3),6] <- dat3$stand_age_sd_45

## prep fof nimble model
Constants3 <- list(N=nrow(dat3),
                   sex=dat3$Sex,
                   WMU=wmu3, # random intercept
                   nWMU=length(unique(dat3$WMU)))

DataBundle3 <- list(y=brom3, # response
                    covars=covars3, # covariates 
                    age=dat3$Age,
                    age2=dat3$age2) 

set.seed(1)
brom.out3 <- nimbleMCMC(code=bromadiolone_code, constants=Constants3, data=DataBundle3, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum3 <- MCMCsummary(brom.out3)
brom.sum3 <- rownames_to_column(brom.sum3, "parameter")
range(brom.sum3$Rhat)

## Iteration 4
dat4 <- dat %>% filter(pt_index==4)

brom4 <- dat4$bin.exp
wmu4 <- as.numeric(factor(dat4$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars4 <- matrix(NA, nrow=nrow(dat4),ncol=6)
covars4[1:nrow(dat4),1] <- dat4$lag_beechnuts
covars4[1:nrow(dat4),2] <- dat4$deciduous_45
covars4[1:nrow(dat4),3] <- dat4$evergreen_45
covars4[1:nrow(dat4),4] <- dat4$nbuildings_45
covars4[1:nrow(dat4),5] <- dat4$stand_age_mean_45
covars4[1:nrow(dat4),6] <- dat4$stand_age_sd_45

## prep fof nimble model
Constants4 <- list(N=nrow(dat4),
                   sex=dat4$Sex,
                   WMU=wmu4, # random intercept
                   nWMU=length(unique(dat4$WMU)))

DataBundle4 <- list(y=brom4, # response
                    covars=covars4, # covariates 
                    age=dat4$Age,
                    age2=dat4$age2) 

set.seed(1)
brom.out4 <- nimbleMCMC(code=bromadiolone_code, constants=Constants4, data=DataBundle4, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum4 <- MCMCsummary(brom.out4)
brom.sum4 <- rownames_to_column(brom.sum4, "parameter")
range(brom.sum4$Rhat)

## Iteration 5
dat5 <- dat %>% filter(pt_index==5)

brom5 <- dat5$bin.exp

# create array for covariate data (column for each covariate)
covars5 <- matrix(NA, nrow=nrow(dat5),ncol=6)
covars5[1:nrow(dat5),1] <- dat5$lag_beechnuts
covars5[1:nrow(dat5),2] <- dat5$deciduous_45
covars5[1:nrow(dat5),3] <- dat5$evergreen_45
covars5[1:nrow(dat5),4] <- dat5$nbuildings_45
covars5[1:nrow(dat5),5] <- dat5$stand_age_mean_45
covars5[1:nrow(dat5),6] <- dat5$stand_age_sd_45

## prep fof nimble model
Constants5 <- list(N=nrow(dat5),
                   sex=dat5$Sex,
                   WMU=as.numeric(factor(dat5$WMU, labels=1:55)), # random intercept
                   nWMU=length(unique(dat5$WMU)))

DataBundle5 <- list(y=brom5, # response
                    covars=covars5, # covariates 
                    age=dat5$Age,
                    age2=dat5$age2) 

set.seed(1)
brom.out5 <- nimbleMCMC(code=bromadiolone_code, constants=Constants5, data=DataBundle5, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum5 <- MCMCsummary(brom.out5)
brom.sum5 <- rownames_to_column(brom.sum5, "parameter")
range(brom.sum5$Rhat)

## Iteration 6
dat6<- dat %>% filter(pt_index==6)

brom6 <- dat6$bin.exp

# create array for covariate data (column for each covariate)
covars6 <- matrix(NA, nrow=nrow(dat6),ncol=6)
covars6[1:nrow(dat6),1] <- dat6$lag_beechnuts
covars6[1:nrow(dat6),2] <- dat6$deciduous_45
covars6[1:nrow(dat6),3] <- dat6$evergreen_45
covars6[1:nrow(dat6),4] <- dat6$nbuildings_45
covars6[1:nrow(dat6),5] <- dat6$stand_age_mean_45
covars6[1:nrow(dat6),6] <- dat6$stand_age_sd_45

## prep fof nimble model
Constants6 <- list(N=nrow(dat6),
                   sex=dat6$Sex,
                   WMU=as.numeric(factor(dat6$WMU, labels=1:55)), # random intercept
                   nWMU=length(unique(dat6$WMU)))

DataBundle6 <- list(y=brom6, # response
                    covars=covars6, # covariates 
                    age=dat6$Age,
                    age2=dat6$age2) 

set.seed(1)
brom.out6 <- nimbleMCMC(code=bromadiolone_code, constants=Constants6, data=DataBundle6, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum6 <- MCMCsummary(brom.out6)
brom.sum6 <- rownames_to_column(brom.sum6, "parameter")
range(brom.sum6$Rhat)

## Iteration 7
dat7 <- dat %>% filter(pt_index==7)

brom7 <- dat7$bin.exp

# create array for covariate data (column for each covariate)
covars7 <- matrix(NA, nrow=nrow(dat7),ncol=6)
covars7[1:nrow(dat7),1] <- dat7$lag_beechnuts
covars7[1:nrow(dat7),2] <- dat7$deciduous_45
covars7[1:nrow(dat7),3] <- dat7$evergreen_45
covars7[1:nrow(dat7),4] <- dat7$nbuildings_45
covars7[1:nrow(dat7),5] <- dat7$stand_age_mean_45
covars7[1:nrow(dat7),6] <- dat7$stand_age_sd_45

## prep fof nimble model
Constants7 <- list(N=nrow(dat7),
                   sex=dat7$Sex,
                   WMU=as.numeric(factor(dat7$WMU, labels=1:55)), # random intercept
                   nWMU=length(unique(dat7$WMU)))

DataBundle7 <- list(y=brom7, # response
                    covars=covars7, # covariates 
                    age=dat7$Age,
                    age2=dat7$age2) 

set.seed(1)
brom.out7 <- nimbleMCMC(code=bromadiolone_code, constants=Constants5, data=DataBundle5, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum7 <- MCMCsummary(brom.out7)
brom.sum7 <- rownames_to_column(brom.sum7, "parameter")
range(brom.sum7$Rhat)


## Iteration 8
dat8 <- dat %>% filter(pt_index==8)

brom8 <- dat8$bin.exp

# create array for covariate data (column for each covariate)
covars8 <- matrix(NA, nrow=nrow(dat8),ncol=6)
covars8[1:nrow(dat8),1] <- dat8$lag_beechnuts
covars8[1:nrow(dat8),2] <- dat8$deciduous_45
covars8[1:nrow(dat8),3] <- dat8$evergreen_45
covars8[1:nrow(dat8),4] <- dat8$nbuildings_45
covars8[1:nrow(dat8),5] <- dat8$stand_age_mean_45
covars8[1:nrow(dat8),6] <- dat8$stand_age_sd_45

## prep fof nimble model
Constants8 <- list(N=nrow(dat8),
                   sex=dat8$Sex,
                   WMU=as.numeric(factor(dat8$WMU, labels=1:55)), # random intercept
                   nWMU=length(unique(dat8$WMU)))

DataBundle8 <- list(y=brom8, # response
                    covars=covars8, # covariates 
                    age=dat8$Age,
                    age2=dat8$age2) 

set.seed(1)
brom.out8 <- nimbleMCMC(code=bromadiolone_code, constants=Constants8, data=DataBundle8, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum8 <- MCMCsummary(brom.out8)
brom.sum8 <- rownames_to_column(brom.sum8, "parameter")
range(brom.sum8$Rhat)

## Iteration 9
dat9 <- dat %>% filter(pt_index==9)

brom9 <- dat9$bin.exp

# create array for covariate data (column for each covariate)
covars9 <- matrix(NA, nrow=nrow(dat9),ncol=6)
covars9[1:nrow(dat9),1] <- dat9$lag_beechnuts
covars9[1:nrow(dat9),2] <- dat9$deciduous_45
covars9[1:nrow(dat9),3] <- dat9$evergreen_45
covars9[1:nrow(dat9),4] <- dat9$nbuildings_45
covars9[1:nrow(dat9),5] <- dat9$stand_age_mean_45
covars9[1:nrow(dat9),6] <- dat9$stand_age_sd_45

## prep fof nimble model
Constants9 <- list(N=nrow(dat9),
                   sex=dat9$Sex,
                   WMU=as.numeric(factor(dat9$WMU, labels=1:55)), # random intercept
                   nWMU=length(unique(dat9$WMU)))

DataBundle9 <- list(y=brom9, # response
                    covars=covars9, # covariates 
                    age=dat9$Age,
                    age2=dat9$age2) 

set.seed(1)
brom.out9 <- nimbleMCMC(code=bromadiolone_code, constants=Constants9, data=DataBundle9, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum9 <- MCMCsummary(brom.out9)
brom.sum9 <- rownames_to_column(brom.sum9, "parameter")
range(brom.sum9$Rhat)

## Iteration 10
dat10 <- dat %>% filter(pt_index==10)

brom10 <- dat10$bin.exp

# create array for covariate data (column for each covariate)
covars10 <- matrix(NA, nrow=nrow(dat10),ncol=6)
covars10[1:nrow(dat10),1] <- dat10$lag_beechnuts
covars10[1:nrow(dat10),2] <- dat10$deciduous_45
covars10[1:nrow(dat10),3] <- dat10$evergreen_45
covars10[1:nrow(dat10),4] <- dat10$nbuildings_45
covars10[1:nrow(dat10),5] <- dat10$stand_age_mean_45
covars10[1:nrow(dat10),6] <- dat10$stand_age_sd_45

## prep fof nimble model
Constants10 <- list(N=nrow(dat10),
                    sex=dat10$Sex,
                    WMU=as.numeric(factor(dat10$WMU, labels=1:55)), # random intercept
                    nWMU=length(unique(dat10$WMU)))

DataBundle10 <- list(y=brom10, # response
                     covars=covars10, # covariates 
                     age=dat10$Age,
                     age2=dat10$age2) 

set.seed(1)
brom.out10 <- nimbleMCMC(code=bromadiolone_code, constants=Constants10, data=DataBundle10, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum10 <- MCMCsummary(brom.out10)
brom.sum10 <- rownames_to_column(brom.sum10, "parameter")
range(brom.sum10$Rhat)


# Combine samples

age <- c(brom.out1$chain1[,56],brom.out2$chain1[,56],brom.out3$chain1[,56],brom.out4$chain1[,56],brom.out5$chain1[,56],
         brom.out6$chain1[,56],brom.out7$chain1[,56],brom.out8$chain1[,56],brom.out9$chain1[,56],brom.out10$chain1[,56],
         brom.out1$chain2[,56],brom.out2$chain2[,56],brom.out3$chain2[,56],brom.out4$chain2[,56],brom.out5$chain2[,56],
         brom.out6$chain2[,56],brom.out7$chain2[,56],brom.out8$chain2[,56],brom.out9$chain2[,56],brom.out10$chain2[,56],
         brom.out1$chain3[,56],brom.out2$chain3[,56],brom.out3$chain3[,56],brom.out4$chain3[,56],brom.out5$chain3[,56],
         brom.out6$chain3[,56],brom.out7$chain3[,56],brom.out8$chain3[,56],brom.out9$chain3[,56],brom.out10$chain3[,56])

# Calculate HDI and quantiles
hdi(age)
quantile(age, probs=c(0.5,0.025,0.975))


age2 <- c(brom.out1$chain1[,57],brom.out2$chain1[,57],brom.out3$chain1[,57],brom.out4$chain1[,57],brom.out5$chain1[,57],
          brom.out6$chain1[,57],brom.out7$chain1[,57],brom.out8$chain1[,57],brom.out9$chain1[,57],brom.out10$chain1[,57],
          brom.out1$chain2[,57],brom.out2$chain2[,57],brom.out3$chain2[,57],brom.out4$chain2[,57],brom.out5$chain2[,57],
          brom.out6$chain2[,57],brom.out7$chain2[,57],brom.out8$chain2[,57],brom.out9$chain2[,57],brom.out10$chain2[,57],
          brom.out1$chain3[,57],brom.out2$chain3[,57],brom.out3$chain3[,57],brom.out4$chain3[,57],brom.out5$chain3[,57],
          brom.out6$chain3[,57],brom.out7$chain3[,57],brom.out8$chain3[,57],brom.out9$chain3[,57],brom.out10$chain3[,57])

# Calculate HDI and quantiles
hdi(age2)
quantile(age2, probs=c(0.5,0.025,0.975))

build <- c(brom.out1$chain1[,58],brom.out2$chain1[,58],brom.out3$chain1[,58],brom.out4$chain1[,58],brom.out5$chain1[,58],
           brom.out6$chain1[,58],brom.out7$chain1[,58],brom.out8$chain1[,58],brom.out9$chain1[,58],brom.out10$chain1[,58],
           brom.out1$chain2[,58],brom.out2$chain2[,58],brom.out3$chain2[,58],brom.out4$chain2[,58],brom.out5$chain2[,58],
           brom.out6$chain2[,58],brom.out7$chain2[,58],brom.out8$chain2[,58],brom.out9$chain2[,58],brom.out10$chain2[,58],
           brom.out1$chain3[,58],brom.out2$chain3[,58],brom.out3$chain3[,58],brom.out4$chain3[,58],brom.out5$chain3[,58],
           brom.out6$chain3[,58],brom.out7$chain3[,58],brom.out8$chain3[,58],brom.out9$chain3[,58],brom.out10$chain3[,58])

# Calculate HDI and quantiles
hdi(build)
quantile(build, probs=c(0.5,0.025,0.975))

decid <- c(brom.out1$chain1[,59],brom.out2$chain1[,59],brom.out3$chain1[,59],brom.out4$chain1[,59],brom.out5$chain1[,59],
           brom.out6$chain1[,59],brom.out7$chain1[,59],brom.out8$chain1[,59],brom.out9$chain1[,59],brom.out10$chain1[,59],
           brom.out1$chain2[,59],brom.out2$chain2[,59],brom.out3$chain2[,59],brom.out4$chain2[,59],brom.out5$chain2[,59],
           brom.out6$chain2[,59],brom.out7$chain2[,59],brom.out8$chain2[,59],brom.out9$chain2[,59],brom.out10$chain2[,59],
           brom.out1$chain3[,59],brom.out2$chain3[,59],brom.out3$chain3[,59],brom.out4$chain3[,59],brom.out5$chain3[,59],
           brom.out6$chain3[,59],brom.out7$chain3[,59],brom.out8$chain3[,59],brom.out9$chain3[,59],brom.out10$chain3[,59])

# Calculate HDI and quantiles
hdi(decid)
quantile(decid, probs=c(0.5,0.025,0.975))

evrgrn <- c(brom.out1$chain1[,60],brom.out2$chain1[,60],brom.out3$chain1[,60],brom.out4$chain1[,60],brom.out5$chain1[,60],
            brom.out6$chain1[,60],brom.out7$chain1[,60],brom.out8$chain1[,60],brom.out9$chain1[,60],brom.out10$chain1[,60],
            brom.out1$chain2[,60],brom.out2$chain2[,60],brom.out3$chain2[,60],brom.out4$chain2[,60],brom.out5$chain2[,60],
            brom.out6$chain2[,60],brom.out7$chain2[,60],brom.out8$chain2[,60],brom.out9$chain2[,60],brom.out10$chain2[,60],
            brom.out1$chain3[,60],brom.out2$chain3[,60],brom.out3$chain3[,60],brom.out4$chain3[,60],brom.out5$chain3[,60],
            brom.out6$chain3[,60],brom.out7$chain3[,60],brom.out8$chain3[,60],brom.out9$chain3[,60],brom.out10$chain3[,60])

# Calculate HDI and quantiles
hdi(evrgrn)
quantile(evrgrn, probs=c(0.5,0.025,0.975))

mast <- c(brom.out1$chain1[,61],brom.out2$chain1[,61],brom.out3$chain1[,61],brom.out4$chain1[,61],brom.out5$chain1[,61],
          brom.out6$chain1[,61],brom.out7$chain1[,61],brom.out8$chain1[,61],brom.out9$chain1[,61],brom.out10$chain1[,61],
          brom.out1$chain2[,61],brom.out2$chain2[,61],brom.out3$chain2[,61],brom.out4$chain2[,61],brom.out5$chain2[,61],
          brom.out6$chain2[,61],brom.out7$chain2[,61],brom.out8$chain2[,61],brom.out9$chain2[,61],brom.out10$chain2[,61],
          brom.out1$chain3[,61],brom.out2$chain3[,61],brom.out3$chain3[,61],brom.out4$chain3[,61],brom.out5$chain3[,61],
          brom.out6$chain3[,61],brom.out7$chain3[,61],brom.out8$chain3[,61],brom.out9$chain3[,61],brom.out10$chain3[,61])

# Calculate HDI and quantiles
hdi(mast)
quantile(mast, probs=c(0.5,0.025,0.975))

sex1 <- c(brom.out1$chain1[,62],brom.out2$chain1[,62],brom.out3$chain1[,62],brom.out4$chain1[,62],brom.out5$chain1[,62],
          brom.out6$chain1[,62],brom.out7$chain1[,62],brom.out8$chain1[,62],brom.out9$chain1[,62],brom.out10$chain1[,62],
          brom.out1$chain2[,62],brom.out2$chain2[,62],brom.out3$chain2[,62],brom.out4$chain2[,62],brom.out5$chain2[,62],
          brom.out6$chain2[,62],brom.out7$chain2[,62],brom.out8$chain2[,62],brom.out9$chain2[,62],brom.out10$chain2[,62],
          brom.out1$chain3[,62],brom.out2$chain3[,62],brom.out3$chain3[,62],brom.out4$chain3[,62],brom.out5$chain3[,62],
          brom.out6$chain3[,62],brom.out7$chain3[,62],brom.out8$chain3[,62],brom.out9$chain3[,62],brom.out10$chain3[,62])

# Calculate HDI and quantiles
hdi(sex1)
quantile(sex1, probs=c(0.5,0.025,0.975))

sex2 <- c(brom.out1$chain1[,63],brom.out2$chain1[,63],brom.out3$chain1[,63],brom.out4$chain1[,63],brom.out5$chain1[,63],
          brom.out6$chain1[,63],brom.out7$chain1[,63],brom.out8$chain1[,63],brom.out9$chain1[,63],brom.out10$chain1[,63],
          brom.out1$chain2[,63],brom.out2$chain2[,63],brom.out3$chain2[,63],brom.out4$chain2[,63],brom.out5$chain2[,63],
          brom.out6$chain2[,63],brom.out7$chain2[,63],brom.out8$chain2[,63],brom.out9$chain2[,63],brom.out10$chain2[,63],
          brom.out1$chain3[,63],brom.out2$chain3[,63],brom.out3$chain3[,63],brom.out4$chain3[,63],brom.out5$chain3[,63],
          brom.out6$chain3[,63],brom.out7$chain3[,63],brom.out8$chain3[,63],brom.out9$chain3[,63],brom.out10$chain3[,63])

# Calculate HDI and quantiles
hdi(sex2)
quantile(sex2, probs=c(0.5,0.025,0.975))

standmn <- c(brom.out1$chain1[,64],brom.out2$chain1[,64],brom.out3$chain1[,64],brom.out4$chain1[,64],brom.out5$chain1[,64],
             brom.out6$chain1[,64],brom.out7$chain1[,64],brom.out8$chain1[,64],brom.out9$chain1[,64],brom.out10$chain1[,64],
             brom.out1$chain2[,64],brom.out2$chain2[,64],brom.out3$chain2[,64],brom.out4$chain2[,64],brom.out5$chain2[,64],
             brom.out6$chain2[,64],brom.out7$chain2[,64],brom.out8$chain2[,64],brom.out9$chain2[,64],brom.out10$chain2[,64],
             brom.out1$chain3[,64],brom.out2$chain3[,64],brom.out3$chain3[,64],brom.out4$chain3[,64],brom.out5$chain3[,64],
             brom.out6$chain3[,64],brom.out7$chain3[,64],brom.out8$chain3[,64],brom.out9$chain3[,64],brom.out10$chain3[,64])

# Calculate HDI and quantiles
hdi(standmn)
quantile(standmn, probs=c(0.5,0.025,0.975))

standsd <- c(brom.out1$chain1[,65],brom.out2$chain1[,65],brom.out3$chain1[,65],brom.out4$chain1[,65],brom.out5$chain1[,65],
             brom.out6$chain1[,65],brom.out7$chain1[,65],brom.out8$chain1[,65],brom.out9$chain1[,65],brom.out10$chain1[,65],
             brom.out1$chain2[,65],brom.out2$chain2[,65],brom.out3$chain2[,65],brom.out4$chain2[,65],brom.out5$chain2[,65],
             brom.out6$chain2[,65],brom.out7$chain2[,65],brom.out8$chain2[,65],brom.out9$chain2[,65],brom.out10$chain2[,65],
             brom.out1$chain3[,65],brom.out2$chain3[,65],brom.out3$chain3[,65],brom.out4$chain3[,65],brom.out5$chain3[,65],
             brom.out6$chain3[,65],brom.out7$chain3[,65],brom.out8$chain3[,65],brom.out9$chain3[,65],brom.out10$chain3[,65])

# Calculate HDI and quantiles
hdi(standsd)
quantile(standsd, probs=c(0.5,0.025,0.975))

