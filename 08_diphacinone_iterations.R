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
diph <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Diphacinone") %>%
  select(RegionalID,bin.exp)
dat <- left_join(dat, diph, by="RegionalID")
dat <- dat %>% select(RegionalID:Town,bin.exp,deciduous:mast_year)

# Build model in BUGS language
diphifacoum_code <- nimbleCode({
  
  ## Priors
  # beta coefficient priors
  beta_age ~ dnorm(0, sd=1.4)
  beta_age2 ~ dnorm(0, sd=1.4)
  beta_sex[1] <- 0
  beta_sex[2] ~ dnorm(0, sd=1.4)
  beta_wui ~ dnorm(0, sd=1.4)
  beta_mast ~ dnorm(0, sd=1.4)
  beta_intx ~ dnorm(0, sd=1.4)
  
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
      beta_mast*beechnuts[i] + beta_wui*wui[i] + beta_intx*wui[i]*beechnuts[i]
    
    y[i] ~ dbern(p[i])
    
  }
  
})

# parameters to monitor
params <- c("beta_age","beta_age2","beta_sex","beta_mast","beta_wui","beta_intx", 
            "alpha", "mu.alpha", "sigma.alpha")  

# MCMC options
nt <- 1
ni <- 50000
nb <- 25000
nc <- 3

set.seed(1)
Inits <- list(sigma.alpha=1, mu.alpha=1, alpha=rnorm(18), #beta0=rnorm(1),
              beta_mast=0, beta_wui=0, beta_intx=rnorm(1), 
              beta_age=0, beta_age2=0, beta_sex=rep(0,2)) 

#### Loop over random locations ####
## Iteration 1
dat1 <- dat %>% filter(pt_index==1)

## Set up data
diph1 <- dat1$bin.exp
wmua1 <- as.numeric(factor(dat1$WMUA_code, labels=1:18))

## prep fof nimble model
Constants1 <- list(N=nrow(dat1),
                   sex=dat1$Sex,
                   WMUA=wmua1, # random intercept
                   nWMUA=length(unique(dat1$WMUA_code)))

DataBundle1 <- list(y=diph1, # response
                    beechnuts=dat1$lag_beechnuts,
                    wui=dat1$intermix,
                    age=dat1$Age,
                    age2=dat1$age2) 

set.seed(1)
diph.out1 <- nimbleMCMC(code=diphifacoum_code, constants=Constants1, data=DataBundle1, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

diph.sum1 <- MCMCsummary(diph.out1)
diph.sum1 <- rownames_to_column(diph.sum1, "parameter")
range(diph.sum1$Rhat, na.rm=TRUE)

plot(diph.out1)

## Iteration 2
dat2 <- dat %>% filter(pt_index==2)

diph2 <- dat2$bin.exp
wmua2 <- as.numeric(factor(dat2$WMUA_code, labels=1:18))

## prep fof nimble model
Constants2 <- list(N=nrow(dat2),
                   sex=dat2$Sex,
                   WMUA=wmua2, # random intercept
                   nWMUA=length(unique(dat2$WMUA_code)))

DataBundle2 <- list(y=diph2, # response
                    beechnuts=dat2$lag_beechnuts,
                    wui=dat2$intermix,
                    age=dat2$Age,
                    age2=dat2$age2) 

set.seed(1)
diph.out2 <- nimbleMCMC(code=diphifacoum_code, constants=Constants2, data=DataBundle2, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

diph.sum2 <- MCMCsummary(diph.out2)
diph.sum2 <- rownames_to_column(diph.sum2, "parameter")
range(diph.sum2$Rhat, na.rm=TRUE)

## Iteration 3
dat3 <- dat %>% filter(pt_index==3)

diph3 <- dat3$bin.exp
wmua3 <- as.numeric(factor(dat3$WMUA_code, labels=1:18))

## prep fof nimble model
Constants3 <- list(N=nrow(dat3),
                   sex=dat3$Sex,
                   WMUA=wmua3, # random intercept
                   nWMUA=length(unique(dat3$WMUA_code)))

DataBundle3 <- list(y=diph3, # response
                    beechnuts=dat3$lag_beechnuts,
                    wui=dat3$intermix,
                    age=dat3$Age,
                    age2=dat3$age2) 

set.seed(1)
diph.out3 <- nimbleMCMC(code=diphifacoum_code, constants=Constants3, data=DataBundle3, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

diph.sum3 <- MCMCsummary(diph.out3)
diph.sum3 <- rownames_to_column(diph.sum3, "parameter")
range(diph.sum3$Rhat, na.rm=TRUE)

## Iteration 4
dat4 <- dat %>% filter(pt_index==4)

diph4 <- dat4$bin.exp
wmua4 <- as.numeric(factor(dat1$WMUA_code, labels=1:18))

## prep fof nimble model
Constants4 <- list(N=nrow(dat4),
                   sex=dat4$Sex,
                   WMUA=wmua4, # random intercept
                   nWMUA=length(unique(dat4$WMUA_code)))

DataBundle4 <- list(y=diph4, # response
                    beechnuts=dat4$lag_beechnuts,
                    wui=dat4$intermix,
                    age=dat4$Age,
                    age2=dat4$age2) 

set.seed(1)
diph.out4 <- nimbleMCMC(code=diphifacoum_code, constants=Constants4, data=DataBundle4, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

diph.sum4 <- MCMCsummary(diph.out4)
diph.sum4 <- rownames_to_column(diph.sum4, "parameter")
range(diph.sum4$Rhat, na.rm=TRUE)

## Iteration 5
dat5 <- dat %>% filter(pt_index==5)

diph5 <- dat5$bin.exp
wmua5 <- as.numeric(factor(dat1$WMUA_code, labels=1:18))

## prep fof nimble model
Constants5 <- list(N=nrow(dat5),
                   sex=dat5$Sex,
                   WMUA=wmua5, # random intercept
                   nWMUA=length(unique(dat5$WMUA_code)))

DataBundle5 <- list(y=diph5, # response
                    beechnuts=dat5$lag_beechnuts,
                    wui=dat5$intermix,
                    age=dat5$Age,
                    age2=dat5$age2) 

set.seed(1)
diph.out5 <- nimbleMCMC(code=diphifacoum_code, constants=Constants5, data=DataBundle5, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

diph.sum5 <- MCMCsummary(diph.out5)
diph.sum5 <- rownames_to_column(diph.sum5, "parameter")
range(diph.sum5$Rhat, na.rm=TRUE)

## Iteration 6
dat6<- dat %>% filter(pt_index==6)

diph6 <- dat6$bin.exp
wmua6 <- as.numeric(factor(dat6$WMUA_code, labels=1:18))

## prep fof nimble model
Constants6 <- list(N=nrow(dat6),
                   sex=dat6$Sex,
                   WMUA=wmua6, # random intercept
                   nWMUA=length(unique(dat6$WMUA_code)))

DataBundle6 <- list(y=diph6, # response
                    beechnuts=dat6$lag_beechnuts,
                    wui=dat6$intermix,
                    age=dat6$Age,
                    age2=dat6$age2) 

set.seed(1)
diph.out6 <- nimbleMCMC(code=diphifacoum_code, constants=Constants6, data=DataBundle6, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

diph.sum6 <- MCMCsummary(diph.out6)
diph.sum6 <- rownames_to_column(diph.sum6, "parameter")
range(diph.sum6$Rhat, na.rm=TRUE)

## Iteration 7
dat7 <- dat %>% filter(pt_index==7)

diph7 <- dat7$bin.exp
wmua7 <- as.numeric(factor(dat7$WMUA_code, labels=1:18))

## prep fof nimble model
Constants7 <- list(N=nrow(dat7),
                   sex=dat7$Sex,
                   WMUA=wmua7, # random intercept
                   nWMUA=length(unique(dat7$WMUA_code)))

DataBundle7 <- list(y=diph7, # response
                    beechnuts=dat7$lag_beechnuts,
                    wui=dat7$intermix,
                    age=dat7$Age,
                    age2=dat7$age2) 

set.seed(1)
diph.out7 <- nimbleMCMC(code=diphifacoum_code, constants=Constants5, data=DataBundle5, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

diph.sum7 <- MCMCsummary(diph.out7)
diph.sum7 <- rownames_to_column(diph.sum7, "parameter")
range(diph.sum7$Rhat, na.rm=TRUE)


## Iteration 8
dat8 <- dat %>% filter(pt_index==8)

diph8 <- dat8$bin.exp
wmua8 <- as.numeric(factor(dat8$WMUA_code, labels=1:18))

## prep fof nimble model
Constants8 <- list(N=nrow(dat8),
                   sex=dat8$Sex,
                   WMUA=wmua8, # random intercept
                   nWMUA=length(unique(dat8$WMUA_code)))

DataBundle8 <- list(y=diph8, # response
                    beechnuts=dat8$lag_beechnuts,
                    wui=dat8$intermix,
                    age=dat8$Age,
                    age2=dat8$age2) 

set.seed(1)
diph.out8 <- nimbleMCMC(code=diphifacoum_code, constants=Constants8, data=DataBundle8, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

diph.sum8 <- MCMCsummary(diph.out8)
diph.sum8 <- rownames_to_column(diph.sum8, "parameter")
range(diph.sum8$Rhat, na.rm=TRUE)

## Iteration 9
dat9 <- dat %>% filter(pt_index==9)

diph9 <- dat9$bin.exp
wmua9 <- as.numeric(factor(dat9$WMUA_code, labels=1:18))

## prep fof nimble model
Constants9 <- list(N=nrow(dat9),
                   sex=dat9$Sex,
                   WMUA=wmua9, # random intercept
                   nWMUA=length(unique(dat9$WMUA_code)))

DataBundle9 <- list(y=diph9, # response
                    beechnuts=dat9$lag_beechnuts,
                    wui=dat9$intermix,
                    age=dat9$Age,
                    age2=dat9$age2) 

set.seed(1)
diph.out9 <- nimbleMCMC(code=diphifacoum_code, constants=Constants9, data=DataBundle9, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

diph.sum9 <- MCMCsummary(diph.out9)
diph.sum9 <- rownames_to_column(diph.sum9, "parameter")
range(diph.sum9$Rhat, na.rm=TRUE)

## Iteration 10
dat10 <- dat %>% filter(pt_index==10)

diph10 <- dat10$bin.exp
wmua10 <- as.numeric(factor(dat1$WMUA_code, labels=1:18))

## prep fof nimble model
Constants10 <- list(N=nrow(dat10),
                    sex=dat10$Sex,
                    WMUA=wmua10, # random intercept
                    nWMUA=length(unique(dat10$WMUA_code)))

DataBundle10 <- list(y=diph10, # response
                     beechnuts=dat10$lag_beechnuts,
                     wui=dat10$intermix,
                     age=dat10$Age,
                     age2=dat10$age2) 

set.seed(1)
diph.out10 <- nimbleMCMC(code=diphifacoum_code, constants=Constants10, data=DataBundle10, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

diph.sum10 <- MCMCsummary(diph.out10)
diph.sum10 <- rownames_to_column(diph.sum10, "parameter")
range(diph.sum10$Rhat, na.rm=TRUE)



## Save model results
saveRDS(diph.sum1, "output/model_output/diph.sum1.rds")
saveRDS(diph.out1, "output/model_output/diph.out1.rds")

saveRDS(diph.sum2, "output/model_output/diph.sum2.rds")
saveRDS(diph.out2, "output/model_output/diph.out2.rds")

saveRDS(diph.sum3, "output/model_output/diph.sum3.rds")
saveRDS(diph.out3, "output/model_output/diph.out3.rds")

saveRDS(diph.sum4, "output/model_output/diph.sum4.rds")
saveRDS(diph.out4, "output/model_output/diph.out4.rds")

saveRDS(diph.sum5, "output/model_output/diph.sum5.rds")
saveRDS(diph.out5, "output/model_output/diph.out5.rds")

saveRDS(diph.sum6, "output/model_output/diph.sum6.rds")
saveRDS(diph.out6, "output/model_output/diph.out6.rds")

saveRDS(diph.sum7, "output/model_output/diph.sum7.rds")
saveRDS(diph.out7, "output/model_output/diph.out7.rds")

saveRDS(diph.sum8, "output/model_output/diph.sum8.rds")
saveRDS(diph.out8, "output/model_output/diph.out8.rds")

saveRDS(diph.sum9, "output/model_output/diph.sum9.rds")
saveRDS(diph.out9, "output/model_output/diph.out9.rds")

saveRDS(diph.sum10, "output/model_output/diph.sum10.rds")
saveRDS(diph.out10, "output/model_output/diph.out10.rds")





