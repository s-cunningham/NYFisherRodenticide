library(tidyverse)
library(nimble)
library(MCMCvis)

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-build_cat) %>%
  mutate(mast_year=if_else(year==2019, 2, 1), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat[,c(8,17:41)] <- scale(dat[,c(8,17:41)])

# read data for individual compounds
brom <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Bromadiolone") %>%
  select(RegionalID,bin.exp)
dat <- left_join(dat, brom, by="RegionalID")
dat <- dat %>% select(RegionalID:Town,bin.exp,deciduous:mast_year)

# Build model in BUGS language
Bromadiolone_code <- nimbleCode({
  
  ## Priors
  # beta coefficient priors
  beta_age ~ dnorm(0, sd=1.4)
  beta_age2 ~ dnorm(0, sd=1.4)
  beta_sex[1] <- 0
  beta_sex[2] ~ dnorm(0, 0.001)
  beta_build ~ dnorm(0, sd=1.4)
  beta_mast ~ dnorm(0, sd=1.4)
  beta_standm ~ dnorm(0, sd=1.4)
  
  ## random intercepts
  # WMU
  for (k in 1:nWMUA) {
    alpha[k] ~ dnorm(mu.alpha, sd=sigma.alpha)
  }
  mu.alpha ~ dnorm(0, 0.001)
  sigma.alpha ~ dunif(0, 1000)
  
  ## Likelihood
  for (i in 1:N) {
    
    logit(p[i]) <- alpha[WMUA[i]] + beta_age*age[i] + beta_age2*age2[i] + beta_sex[sex[i]] + 
      beta_mast*covars[i,1] + beta_build*covars[i,2] + beta_standm*covars[i,3]
    
    y[i] ~ dbern(p[i])
    
  }
  
})

# parameters to monitor
params <- c("beta_age","beta_age2","beta_sex","beta_mast",
            "beta_build","beta_standm",  "alpha", "mu.alpha", "sigma.alpha")  

# MCMC options
nt <- 1
ni <- 25000
nb <- 12000
nc <- 3

set.seed(1)
Inits <- list(sigma.alpha=1, mu.alpha=1, alpha=rnorm(18), #beta0=rnorm(1),
              beta_mast=rnorm(1), beta_build=rnorm(1), beta_standm=rnorm(1), 
              beta_age=rnorm(1), beta_age2=rnorm(1), beta_sex=rnorm(2)) 

#### Loop over random locations ####
## Iteration 1
dat1 <- dat %>% filter(pt_index==1)

## Set up data
brom1 <- dat1$bin.exp
wmua1 <- as.numeric(factor(dat1$WMUA_code, labels=1:18))

# create array for covariate data (column for each covariate)
covars1 <- matrix(NA, nrow=nrow(dat1),ncol=3)
covars1[1:nrow(dat1),1] <- dat1$lag_beechnuts
covars1[1:nrow(dat1),2] <- dat1$nbuildings
covars1[1:nrow(dat1),3] <- dat1$stand_age_mean

## prep fof nimble model
Constants1 <- list(N=nrow(dat1),
                   sex=dat1$Sex,
                   WMUA=wmua1, # random intercept
                   nWMUA=length(unique(dat1$WMUA_code)))

DataBundle1 <- list(y=brom1, # response
                    covars=covars1, # covariates 
                    age=dat1$Age,
                    age2=dat1$age2) 

set.seed(1)
brom.out1 <- nimbleMCMC(code=Bromadiolone_code, constants=Constants1, data=DataBundle1, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum1 <- MCMCsummary(brom.out1)
brom.sum1 <- rownames_to_column(brom.sum1, "parameter")
range(brom.sum1$Rhat, na.rm=TRUE)

plot(brom.out1)

## Iteration 2
dat2 <- dat %>% filter(pt_index==2)

brom2 <- dat2$bin.exp
wmua2 <- as.numeric(factor(dat2$WMUA_code, labels=1:18))

# create array for covariate data (column for each covariate)
covars2 <- matrix(NA, nrow=nrow(dat2),ncol=3)
covars2[1:nrow(dat2),1] <- dat2$lag_beechnuts
covars2[1:nrow(dat2),2] <- dat2$nbuildings
covars2[1:nrow(dat2),3] <- dat2$stand_age_mean

## prep fof nimble model
Constants2 <- list(N=nrow(dat2),
                   sex=dat2$Sex,
                   WMUA=wmua2, # random intercept
                   nWMUA=length(unique(dat2$WMUA_code)))

DataBundle2 <- list(y=brom2, # response
                    covars=covars2, # covariates 
                    age=dat2$Age,
                    age2=dat2$age2) 

set.seed(1)
brom.out2 <- nimbleMCMC(code=Bromadiolone_code, constants=Constants2, data=DataBundle2, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum2 <- MCMCsummary(brom.out2)
brom.sum2 <- rownames_to_column(brom.sum2, "parameter")
range(brom.sum2$Rhat, na.rm=TRUE)

## Iteration 3
dat3 <- dat %>% filter(pt_index==3)

brom3 <- dat3$bin.exp
wmua3 <- as.numeric(factor(dat3$WMUA_code, labels=1:18))

# create array for covariate data (column for each covariate)
covars3 <- matrix(NA, nrow=nrow(dat3),ncol=3)
covars3[1:nrow(dat3),1] <- dat3$lag_beechnuts
covars3[1:nrow(dat3),2] <- dat3$nbuildings
covars3[1:nrow(dat3),3] <- dat3$stand_age_mean

## prep fof nimble model
Constants3 <- list(N=nrow(dat3),
                   sex=dat3$Sex,
                   WMUA=wmua3, # random intercept
                   nWMUA=length(unique(dat3$WMUA_code)))

DataBundle3 <- list(y=brom3, # response
                    covars=covars3, # covariates 
                    age=dat3$Age,
                    age2=dat3$age2) 

set.seed(1)
brom.out3 <- nimbleMCMC(code=Bromadiolone_code, constants=Constants3, data=DataBundle3, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum3 <- MCMCsummary(brom.out3)
brom.sum3 <- rownames_to_column(brom.sum3, "parameter")
range(brom.sum3$Rhat, na.rm=TRUE)

## Iteration 4
dat4 <- dat %>% filter(pt_index==4)

brom4 <- dat4$bin.exp
wmua4 <- as.numeric(factor(dat1$WMUA_code, labels=1:18))

# create array for covariate data (column for each covariate)
covars4 <- matrix(NA, nrow=nrow(dat4),ncol=3)
covars4[1:nrow(dat4),1] <- dat4$lag_beechnuts
covars4[1:nrow(dat4),2] <- dat4$nbuildings
covars4[1:nrow(dat4),3] <- dat4$stand_age_mean

## prep fof nimble model
Constants4 <- list(N=nrow(dat4),
                   sex=dat4$Sex,
                   WMUA=wmua4, # random intercept
                   nWMUA=length(unique(dat4$WMUA_code)))

DataBundle4 <- list(y=brom4, # response
                    covars=covars4, # covariates 
                    age=dat4$Age,
                    age2=dat4$age2) 

set.seed(1)
brom.out4 <- nimbleMCMC(code=Bromadiolone_code, constants=Constants4, data=DataBundle4, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum4 <- MCMCsummary(brom.out4)
brom.sum4 <- rownames_to_column(brom.sum4, "parameter")
range(brom.sum4$Rhat, na.rm=TRUE)

## Iteration 5
dat5 <- dat %>% filter(pt_index==5)

brom5 <- dat5$bin.exp
wmua5 <- as.numeric(factor(dat1$WMUA_code, labels=1:18))

# create array for covariate data (column for each covariate)
covars5 <- matrix(NA, nrow=nrow(dat5),ncol=3)
covars5[1:nrow(dat5),1] <- dat5$lag_beechnuts
covars5[1:nrow(dat5),2] <- dat5$nbuildings
covars5[1:nrow(dat5),3] <- dat5$stand_age_mean

## prep fof nimble model
Constants5 <- list(N=nrow(dat5),
                   sex=dat5$Sex,
                   WMUA=wmua5, # random intercept
                   nWMUA=length(unique(dat5$WMUA_code)))

DataBundle5 <- list(y=brom5, # response
                    covars=covars5, # covariates 
                    age=dat5$Age,
                    age2=dat5$age2) 

set.seed(1)
brom.out5 <- nimbleMCMC(code=Bromadiolone_code, constants=Constants5, data=DataBundle5, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum5 <- MCMCsummary(brom.out5)
brom.sum5 <- rownames_to_column(brom.sum5, "parameter")
range(brom.sum5$Rhat, na.rm=TRUE)

## Iteration 6
dat6<- dat %>% filter(pt_index==6)

brom6 <- dat6$bin.exp
wmua6 <- as.numeric(factor(dat6$WMUA_code, labels=1:18))

# create array for covariate data (column for each covariate)
covars6 <- matrix(NA, nrow=nrow(dat6),ncol=3)
covars6[1:nrow(dat6),1] <- dat6$lag_beechnuts
covars6[1:nrow(dat6),2] <- dat6$nbuildings
covars6[1:nrow(dat6),3] <- dat6$stand_age_mean

## prep fof nimble model
Constants6 <- list(N=nrow(dat6),
                   sex=dat6$Sex,
                   WMUA=wmu6, # random intercept
                   nWMUA=length(unique(dat6$WMUA_code)))

DataBundle6 <- list(y=brom6, # response
                    covars=covars6, # covariates 
                    age=dat6$Age,
                    age2=dat6$age2) 

set.seed(1)
brom.out6 <- nimbleMCMC(code=Bromadiolone_code, constants=Constants6, data=DataBundle6, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum6 <- MCMCsummary(brom.out6)
brom.sum6 <- rownames_to_column(brom.sum6, "parameter")
range(brom.sum6$Rhat, na.rm=TRUE)

## Iteration 7
dat7 <- dat %>% filter(pt_index==7)

brom7 <- dat7$bin.exp
wmua7 <- as.numeric(factor(dat7$WMUA_code, labels=1:18))

# create array for covariate data (column for each covariate)
covars7 <- matrix(NA, nrow=nrow(dat7),ncol=3)
covars7[1:nrow(dat7),1] <- dat7$lag_beechnuts
covars7[1:nrow(dat7),2] <- dat7$nbuildings
covars7[1:nrow(dat7),3] <- dat7$stand_age_mean

## prep fof nimble model
Constants7 <- list(N=nrow(dat7),
                   sex=dat7$Sex,
                   WMUA=wmua7, # random intercept
                   nWMUA=length(unique(dat7$WMUA_code)))

DataBundle7 <- list(y=brom7, # response
                    covars=covars7, # covariates 
                    age=dat7$Age,
                    age2=dat7$age2) 

set.seed(1)
brom.out7 <- nimbleMCMC(code=Bromadiolone_code, constants=Constants5, data=DataBundle5, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum7 <- MCMCsummary(brom.out7)
brom.sum7 <- rownames_to_column(brom.sum7, "parameter")
range(brom.sum7$Rhat, na.rm=TRUE)


## Iteration 8
dat8 <- dat %>% filter(pt_index==8)

brom8 <- dat8$bin.exp
wmua8 <- as.numeric(factor(dat8$WMUA_code, labels=1:18))

# create array for covariate data (column for each covariate)
covars8 <- matrix(NA, nrow=nrow(dat8),ncol=3)
covars8[1:nrow(dat8),1] <- dat8$lag_beechnuts
covars8[1:nrow(dat8),2] <- dat8$nbuildings
covars8[1:nrow(dat8),3] <- dat8$stand_age_mean

## prep fof nimble model
Constants8 <- list(N=nrow(dat8),
                   sex=dat8$Sex,
                   WMUA=wmua8, # random intercept
                   nWMUA=length(unique(dat8$WMUA_code)))

DataBundle8 <- list(y=brom8, # response
                    covars=covars8, # covariates 
                    age=dat8$Age,
                    age2=dat8$age2) 

set.seed(1)
brom.out8 <- nimbleMCMC(code=Bromadiolone_code, constants=Constants8, data=DataBundle8, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum8 <- MCMCsummary(brom.out8)
brom.sum8 <- rownames_to_column(brom.sum8, "parameter")
range(brom.sum8$Rhat, na.rm=TRUE)

## Iteration 9
dat9 <- dat %>% filter(pt_index==9)

brom9 <- dat9$bin.exp
wmua9 <- as.numeric(factor(dat9$WMUA_code, labels=1:18))

# create array for covariate data (column for each covariate)
covars9 <- matrix(NA, nrow=nrow(dat9),ncol=3)
covars9[1:nrow(dat9),1] <- dat9$lag_beechnuts
covars9[1:nrow(dat9),2] <- dat9$nbuildings
covars9[1:nrow(dat9),3] <- dat9$stand_age_mean

## prep fof nimble model
Constants9 <- list(N=nrow(dat9),
                   sex=dat9$Sex,
                   WMUA=wmua9, # random intercept
                   nWMUA=length(unique(dat9$WMUA_code)))

DataBundle9 <- list(y=brom9, # response
                    covars=covars9, # covariates 
                    age=dat9$Age,
                    age2=dat9$age2) 

set.seed(1)
brom.out9 <- nimbleMCMC(code=Bromadiolone_code, constants=Constants9, data=DataBundle9, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum9 <- MCMCsummary(brom.out9)
brom.sum9 <- rownames_to_column(brom.sum9, "parameter")
range(brom.sum9$Rhat, na.rm=TRUE)

## Iteration 10
dat10 <- dat %>% filter(pt_index==10)

brom10 <- dat10$bin.exp
wmua10 <- as.numeric(factor(dat1$WMUA_code, labels=1:18))

# create array for covariate data (column for each covariate)
covars10 <- matrix(NA, nrow=nrow(dat10),ncol=3)
covars10[1:nrow(dat10),1] <- dat10$lag_beechnuts
covars10[1:nrow(dat10),2] <- dat10$nbuildings
covars10[1:nrow(dat10),3] <- dat10$stand_age_mean

## prep fof nimble model
Constants10 <- list(N=nrow(dat10),
                    sex=dat10$Sex,
                    WMUA=wmua10, # random intercept
                    nWMUA=length(unique(dat10$WMUA_code)))

DataBundle10 <- list(y=brom10, # response
                     covars=covars10, # covariates 
                     age=dat10$Age,
                     age2=dat10$age2) 

set.seed(1)
brom.out10 <- nimbleMCMC(code=Bromadiolone_code, constants=Constants10, data=DataBundle10, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

brom.sum10 <- MCMCsummary(brom.out10)
brom.sum10 <- rownames_to_column(brom.sum10, "parameter")
range(brom.sum10$Rhat, na.rm=TRUE)



## Save model results
saveRDS(brom.sum1, "output/model_output/brom.sum1.rds")
saveRDS(brom.out1, "output/model_output/brom.out1.rds")

saveRDS(brom.sum2, "output/model_output/brom.sum2.rds")
saveRDS(brom.out2, "output/model_output/brom.out2.rds")

saveRDS(brom.sum3, "output/model_output/brom.sum3.rds")
saveRDS(brom.out3, "output/model_output/brom.out3.rds")

saveRDS(brom.sum4, "output/model_output/brom.sum4.rds")
saveRDS(brom.out4, "output/model_output/brom.out4.rds")

saveRDS(brom.sum5, "output/model_output/brom.sum5.rds")
saveRDS(brom.out5, "output/model_output/brom.out5.rds")

saveRDS(brom.sum6, "output/model_output/brom.sum6.rds")
saveRDS(brom.out6, "output/model_output/brom.out6.rds")

saveRDS(brom.sum7, "output/model_output/brom.sum7.rds")
saveRDS(brom.out7, "output/model_output/brom.out7.rds")

saveRDS(brom.sum8, "output/model_output/brom.sum8.rds")
saveRDS(brom.out8, "output/model_output/brom.out8.rds")

saveRDS(brom.sum9, "output/model_output/brom.sum9.rds")
saveRDS(brom.out9, "output/model_output/brom.out9.rds")

saveRDS(brom.sum10, "output/model_output/brom.sum10.rds")
saveRDS(brom.out10, "output/model_output/brom.out10.rds")



