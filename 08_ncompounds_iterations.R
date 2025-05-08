## Nimble model after variable and scale selection
library(tidyverse)
library(nimble)
library(COMPoissonReg)
library(MCMCvis)
library(HDInterval)


## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(RegionalID:ncomp,nbuildings,notWUI:totalWUI,beechnuts:age2) %>%
  mutate(mast_year=if_else(year==2019, 1, 2), # mast years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat[,c(8,16:23)] <- scale(dat[,c(8,16:23)])

## Nimble-ize Conway-Maxwell Poisson functions
# Random values function
rCOMP <- nimbleRcall(
  function(n=double(0), lambda=double(0), nu=double(0)){},
  Rfun = 'rcmp',
  returnType=double()) 

# Density function
dCOMP <- nimbleRcall(
  function(x=double(0), lambda=double(0), nu=double(0), log=double(0, default=0)){},
  Rfun = 'dcmp',
  returnType=double()) 

# define in global environment
assign('dCOMP', dCOMP, envir=.GlobalEnv)
assign('rCOMP', rCOMP, envir=.GlobalEnv)

# Build model in BUGS language
ncompounds_code <- nimbleCode({
  
  ## Priors
  nu ~ dunif(1,2.5) # prior for CMP dispersion parameter
  
  # beta coefficient priors
  beta_age ~ dnorm(0, 0.01)
  beta_age2 ~ dnorm(0, 0.01)
  beta_sex[1] <- 0
  beta_sex[2] ~ dnorm(0, 0.01)
  beta_wui ~ dnorm(0, 0.01)
  beta_mast ~ dnorm(0, 0.01)
  beta_intx ~ dnorm(0, 0.01)
  
  ## random intercepts - WMUA
  for (k in 1:nWMUA) {
    alpha[k] ~ dnorm(mu.alpha, tau.alpha)
  }
  mu.alpha ~ dnorm(0, 0.001)
  sigma.alpha ~ dunif(0, 10)
  tau.alpha <- 1/(sigma.alpha*sigma.alpha)
  
  ## Likelihood
  for (i in 1:N) {
    
    lambda[i] <- exp(alpha[WMUA[i]] + beta_age*age[i] + beta_age2*age2[i] + beta_sex[sex[i]] + 
                       beta_wui*wui[i] + beta_mast*beechnuts[i] + beta_intx*wui[i]*beechnuts[i]) 

    ncomp[i] ~ dCOMP(lambda[i], nu)
    
  }
  
})

# parameters to monitor
params <- c("beta_age","beta_age2","beta_sex", "beta_wui", "beta_mast", "beta_intx",
            "nu", "alpha", "mu.alpha", "sigma.alpha")  

# MCMC options
nt <- 1
ni <- 150000
nb <- 75000
nc <- 3

set.seed(1)
Inits <- list(nu=1.5, sigma.alpha=1, alpha=rnorm(18), mu.alpha=rnorm(1), beta_mast=0, beta_intx=0, 
              beta_sex=rep(0,2), beta_age=0, beta_age2=0, beta_wui=0) #

#### Loop over random locations ####
## Iteration 1
dat1 <- dat %>% filter(pt_index==1)

## Set up data
ncomp1 <- dat1$ncomp
wmua1 <- as.numeric(factor(dat1$WMUA_code, labels=1:18))

## prep fof nimble model
Constants1 <- list(N=nrow(dat1),
                  sex=dat1$Sex,
                   WMUA=wmua1, # random intercept
                   nWMUA=length(unique(dat1$WMUA_code)))

DataBundle1 <- list(ncomp=ncomp1, # response
                    wui=dat1$intermix, # covariates (buildings)
                    beechnuts=dat1$lag_beechnuts,   
                    age=dat1$Age,
                    age2=dat1$age2) 

set.seed(1)
  ncomp.out1 <- nimbleMCMC(code=ncompounds_code, constants=Constants1, data=DataBundle1, 
                           inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                           nchains=nc, check=FALSE, samples=TRUE,
                           samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum1 <- MCMCsummary(ncomp.out1)
ncomp.sum1 <- rownames_to_column(ncomp.sum1, "parameter")
range(ncomp.sum1$Rhat, na.rm=TRUE)

MCMCtrace(ncomp.out1, pdf = FALSE)


saveRDS(ncomp.sum1, "output/model_output/ncomp.sum1.rds")
saveRDS(ncomp.out1, "output/model_output/ncomp.out1.rds")


## Iteration 2
dat2 <- dat %>% filter(pt_index==2)

## Set up data
ncomp2 <- dat2$ncomp
wmua2 <- as.numeric(factor(dat2$WMUA_code, labels=1:18))

## prep fof nimble model
Constants2 <- list(N=nrow(dat2),
                   sex=dat2$Sex,
                   WMUA=wmua2, # random intercept
                   nWMUA=length(unique(dat2$WMUA_code)))

DataBundle2 <- list(ncomp=ncomp2, # response
                    wui=dat2$intermix, # covariates (buildings)
                    beechnuts=dat2$lag_beechnuts,   
                    age=dat2$Age,
                    age2=dat2$age2) 

set.seed(1)
ncomp.out2 <- nimbleMCMC(code=ncompounds_code, constants=Constants2, data=DataBundle2, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum2 <- MCMCsummary(ncomp.out2)
ncomp.sum2 <- rownames_to_column(ncomp.sum2, "parameter")
range(ncomp.sum2$Rhat, na.rm=TRUE)


saveRDS(ncomp.sum2, "output/model_output/ncomp.sum2.rds")
saveRDS(ncomp.out2, "output/model_output/ncomp.out2.rds")


## Iteration 3
dat3 <- dat %>% filter(pt_index==3)

## Set up data
ncomp3 <- dat3$ncomp
wmua3 <- as.numeric(factor(dat3$WMUA_code, labels=1:18))

## prep fof nimble model
Constants3 <- list(N=nrow(dat3),
                   sex=dat3$Sex,
                   WMUA=wmua3, # random intercept
                   nWMUA=length(unique(dat3$WMUA_code)))

DataBundle3 <- list(ncomp=ncomp3, # response
                    wui=dat3$intermix, 
                    beechnuts=dat3$lag_beechnuts,   
                    age=dat3$Age,
                    age2=dat3$age2) 

set.seed(1)
ncomp.out3 <- nimbleMCMC(code=ncompounds_code, constants=Constants3, data=DataBundle3, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum3 <- MCMCsummary(ncomp.out3)
ncomp.sum3 <- rownames_to_column(ncomp.sum3, "parameter")
range(ncomp.sum3$Rhat, na.rm=TRUE)


saveRDS(ncomp.sum3, "output/model_output/ncomp.sum3.rds")
saveRDS(ncomp.out3, "output/model_output/ncomp.out3.rds")


## Iteration 4
dat4 <- dat %>% filter(pt_index==4)

## Set up data
ncomp4 <- dat4$ncomp
wmua4 <- as.numeric(factor(dat4$WMUA_code, labels=1:18))

## prep fof nimble model
Constants4 <- list(N=nrow(dat4),
                   sex=dat4$Sex,
                   WMUA=wmua4, # random intercept
                   nWMUA=length(unique(dat4$WMUA_code)))

DataBundle4 <- list(ncomp=ncomp4, # response
                    wui=dat4$intermix, 
                    beechnuts=dat4$lag_beechnuts,   
                    age=dat4$Age,
                    age2=dat4$age2) 

set.seed(1)
ncomp.out4 <- nimbleMCMC(code=ncompounds_code, constants=Constants4, data=DataBundle4, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum4 <- MCMCsummary(ncomp.out4)
ncomp.sum4 <- rownames_to_column(ncomp.sum4, "parameter")
range(ncomp.sum4$Rhat, na.rm=TRUE)

saveRDS(ncomp.sum4, "output/model_output/ncomp.sum4.rds")
saveRDS(ncomp.out4, "output/model_output/ncomp.out4.rds")


## Iteration 5
dat5 <- dat %>% filter(pt_index==5)

## Set up data
ncomp5 <- dat5$ncomp
wmua5 <- as.numeric(factor(dat5$WMUA_code, labels=1:18))

## prep fof nimble model
Constants5 <- list(N=nrow(dat5),
                   sex=dat5$Sex,
                   WMUA=wmua5, # random intercept
                   nWMUA=length(unique(dat5$WMUA_code)))

DataBundle5 <- list(ncomp=ncomp5, # response
                    wui=dat5$intermix, 
                    beechnuts=dat5$lag_beechnuts,   
                    age=dat5$Age,
                    age2=dat5$age2) 

set.seed(1)
ncomp.out5 <- nimbleMCMC(code=ncompounds_code, constants=Constants5, data=DataBundle5, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum5 <- MCMCsummary(ncomp.out5)
ncomp.sum5 <- rownames_to_column(ncomp.sum5, "parameter")
range(ncomp.sum5$Rhat, na.rm=TRUE)

saveRDS(ncomp.sum5, "output/model_output/ncomp.sum5.rds")
saveRDS(ncomp.out5, "output/model_output/ncomp.out5.rds")


## Iteration 6
dat6 <- dat %>% filter(pt_index==6)

## Set up data
ncomp6 <- dat6$ncomp
wmua6 <- as.numeric(factor(dat6$WMUA_code, labels=1:18))

## prep fof nimble model
Constants6 <- list(N=nrow(dat6),
                   sex=dat6$Sex,
                   WMUA=wmua6, # random intercept
                   nWMUA=length(unique(dat6$WMUA_code)))

DataBundle6 <- list(ncomp=ncomp6, # response
                    wui=dat6$intermix, 
                    beechnuts=dat6$lag_beechnuts,   
                    age=dat6$Age,
                    age2=dat6$age2) 

set.seed(1)
ncomp.out6 <- nimbleMCMC(code=ncompounds_code, constants=Constants6, data=DataBundle6, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum6 <- MCMCsummary(ncomp.out6)
ncomp.sum6 <- rownames_to_column(ncomp.sum6, "parameter")
range(ncomp.sum6$Rhat, na.rm=TRUE)

## Iteration 7
dat7 <- dat %>% filter(pt_index==7)

## Set up data
ncomp7 <- dat7$ncomp
wmua7 <- as.numeric(factor(dat7$WMUA_code, labels=1:18))

## prep fof nimble model
Constants7 <- list(N=nrow(dat7),
                   sex=dat7$Sex,
                   WMUA=wmua7, # random intercept
                   nWMUA=length(unique(dat7$WMUA_code)))

DataBundle7 <- list(ncomp=ncomp7, # response
                    wui=dat7$intermix, 
                    beechnuts=dat7$lag_beechnuts,   
                    age=dat7$Age,
                    age2=dat7$age2) 

set.seed(1)
ncomp.out7 <- nimbleMCMC(code=ncompounds_code, constants=Constants7, data=DataBundle7, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum7 <- MCMCsummary(ncomp.out7)
ncomp.sum7 <- rownames_to_column(ncomp.sum7, "parameter")
range(ncomp.sum7$Rhat, na.rm=TRUE)

## Iteration 8
dat8 <- dat %>% filter(pt_index==8)

## Set up data
ncomp8 <- dat8$ncomp
wmua8 <- as.numeric(factor(dat8$WMUA_code, labels=1:18))

## prep fof nimble model
Constants8 <- list(N=nrow(dat8),
                   sex=dat8$Sex,
                   WMUA=wmua8, # random intercept
                   nWMUA=length(unique(dat8$WMUA_code)))

DataBundle8 <- list(ncomp=ncomp8, # response
                    wui=dat8$intermix, 
                    beechnuts=dat8$lag_beechnuts,   
                    age=dat8$Age,
                    age2=dat8$age2) 

set.seed(1)
ncomp.out8 <- nimbleMCMC(code=ncompounds_code, constants=Constants8, data=DataBundle8, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum8 <- MCMCsummary(ncomp.out8)
ncomp.sum8 <- rownames_to_column(ncomp.sum8, "parameter")
range(ncomp.sum8$Rhat, na.rm=TRUE)

saveRDS(ncomp.sum8, "output/model_output/ncomp.sum8.rds")
saveRDS(ncomp.out8, "output/model_output/ncomp.out8.rds")

## Iteration 9
dat9 <- dat %>% filter(pt_index==9)

## Set up data
ncomp9 <- dat9$ncomp
wmua9 <- as.numeric(factor(dat9$WMUA_code, labels=1:18))

## prep fof nimble model
Constants9 <- list(N=nrow(dat9),
                   sex=dat9$Sex,
                   WMUA=wmua9, # random intercept
                   nWMUA=length(unique(dat9$WMUA_code)))

DataBundle9 <- list(ncomp=ncomp9, # response
                    wui=dat9$intermix, 
                    beechnuts=dat9$lag_beechnuts,   
                    age=dat9$Age,
                    age2=dat9$age2) 

set.seed(1)
ncomp.out9 <- nimbleMCMC(code=ncompounds_code, constants=Constants9, data=DataBundle9, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum9 <- MCMCsummary(ncomp.out9)
ncomp.sum9 <- rownames_to_column(ncomp.sum9, "parameter")
range(ncomp.sum9$Rhat, na.rm=TRUE)

saveRDS(ncomp.sum9, "output/model_output/ncomp.sum9.rds")
saveRDS(ncomp.out9, "output/model_output/ncomp.out9.rds")

## Iteration 10
dat10 <- dat %>% filter(pt_index==10)

## Set up data
ncomp10 <- dat10$ncomp
wmua10 <- as.numeric(factor(dat10$WMUA_code, labels=1:18))

## prep fof nimble model
Constants10 <- list(N=nrow(dat10),
                   sex=dat10$Sex,
                   WMUA=wmua10, # random intercept
                   nWMUA=length(unique(dat10$WMUA_code)))

DataBundle10 <- list(ncomp=ncomp10, # response
                    wui=dat10$intermix, 
                    beechnuts=dat10$lag_beechnuts,   
                    age=dat10$Age,
                    age2=dat10$age2) 

set.seed(1)
ncomp.out10 <- nimbleMCMC(code=ncompounds_code, constants=Constants10, data=DataBundle10, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum10 <- MCMCsummary(ncomp.out10)
ncomp.sum10 <- rownames_to_column(ncomp.sum10, "parameter")
range(ncomp.sum10$Rhat, na.rm=TRUE)

### save the output data

saveRDS(ncomp.sum1, "output/model_output/ncomp.sum1.rds")
saveRDS(ncomp.out1, "output/model_output/ncomp.out1.rds")

saveRDS(ncomp.sum2, "output/model_output/ncomp.sum2.rds")
saveRDS(ncomp.out3, "output/model_output/ncomp.out2.rds")

saveRDS(ncomp.sum3, "output/model_output/ncomp.sum3.rds")
saveRDS(ncomp.out3, "output/model_output/ncomp.out3.rds")

saveRDS(ncomp.sum4, "output/model_output/ncomp.sum4.rds")
saveRDS(ncomp.out4, "output/model_output/ncomp.out4.rds")

saveRDS(ncomp.sum5, "output/model_output/ncomp.sum5.rds")
saveRDS(ncomp.out5, "output/model_output/ncomp.out5.rds")

saveRDS(ncomp.sum6, "output/model_output/ncomp.sum6.rds")
saveRDS(ncomp.out6, "output/model_output/ncomp.out6.rds")

saveRDS(ncomp.sum7, "output/model_output/ncomp.sum7.rds")
saveRDS(ncomp.out7, "output/model_output/ncomp.out7.rds")

saveRDS(ncomp.sum8, "output/model_output/ncomp.sum8.rds")
saveRDS(ncomp.out8, "output/model_output/ncomp.out8.rds")

saveRDS(ncomp.sum9, "output/model_output/ncomp.sum9.rds")
saveRDS(ncomp.out9, "output/model_output/ncomp.out9.rds")

saveRDS(ncomp.sum10, "output/model_output/ncomp.sum10.rds")
saveRDS(ncomp.out10, "output/model_output/ncomp.out10.rds")

# Combine samples (still need to check column numbers)

age <- c(ncomp.out1$chain1[,19],ncomp.out2$chain1[,19],ncomp.out3$chain1[,19],ncomp.out4$chain1[,19],ncomp.out5$chain1[,19],
         ncomp.out6$chain1[,19],ncomp.out7$chain1[,19],ncomp.out8$chain1[,19],ncomp.out9$chain1[,19],ncomp.out10$chain1[,19],
         ncomp.out1$chain2[,19],ncomp.out2$chain2[,19],ncomp.out3$chain2[,19],ncomp.out4$chain2[,19],ncomp.out5$chain2[,19],
         ncomp.out6$chain2[,19],ncomp.out7$chain2[,19],ncomp.out8$chain2[,19],ncomp.out9$chain2[,19],ncomp.out10$chain2[,19],
         ncomp.out1$chain3[,19],ncomp.out2$chain3[,19],ncomp.out3$chain3[,19],ncomp.out4$chain3[,19],ncomp.out5$chain3[,19],
         ncomp.out6$chain3[,19],ncomp.out7$chain3[,19],ncomp.out8$chain3[,19],ncomp.out9$chain3[,19],ncomp.out10$chain3[,19])

# Calculate HDI and quantiles
hdi(age)
quantile(age, probs=c(0.5,0.025,0.975))


age2 <- c(ncomp.out1$chain1[,20],ncomp.out2$chain1[,20],ncomp.out3$chain1[,20],ncomp.out4$chain1[,20],ncomp.out5$chain1[,20],
          ncomp.out6$chain1[,20],ncomp.out7$chain1[,20],ncomp.out8$chain1[,20],ncomp.out9$chain1[,20],ncomp.out10$chain1[,20],
          ncomp.out1$chain2[,20],ncomp.out2$chain2[,20],ncomp.out3$chain2[,20],ncomp.out4$chain2[,20],ncomp.out5$chain2[,20],
          ncomp.out6$chain2[,20],ncomp.out7$chain2[,20],ncomp.out8$chain2[,20],ncomp.out9$chain2[,20],ncomp.out10$chain2[,20],
          ncomp.out1$chain3[,20],ncomp.out2$chain3[,20],ncomp.out3$chain3[,20],ncomp.out4$chain3[,20],ncomp.out5$chain3[,20],
          ncomp.out6$chain3[,20],ncomp.out7$chain3[,20],ncomp.out8$chain3[,20],ncomp.out9$chain3[,20],ncomp.out10$chain3[,20])

# Calculate HDI and quantiles
hdi(age2)
quantile(age2, probs=c(0.5,0.025,0.975))

build <- c(ncomp.out1$chain1[,21],ncomp.out2$chain1[,21],ncomp.out3$chain1[,21],ncomp.out4$chain1[,21],ncomp.out5$chain1[,21],
           ncomp.out6$chain1[,21],ncomp.out7$chain1[,21],ncomp.out8$chain1[,21],ncomp.out9$chain1[,21],ncomp.out10$chain1[,21],
           ncomp.out1$chain2[,21],ncomp.out2$chain2[,21],ncomp.out3$chain2[,21],ncomp.out4$chain2[,21],ncomp.out5$chain2[,21],
           ncomp.out6$chain2[,21],ncomp.out7$chain2[,21],ncomp.out8$chain2[,21],ncomp.out9$chain2[,21],ncomp.out10$chain2[,21],
           ncomp.out1$chain3[,21],ncomp.out2$chain3[,21],ncomp.out3$chain3[,21],ncomp.out4$chain3[,21],ncomp.out5$chain3[,21],
           ncomp.out6$chain3[,21],ncomp.out7$chain3[,21],ncomp.out8$chain3[,21],ncomp.out9$chain3[,21],ncomp.out10$chain3[,21])

# Calculate HDI and quantiles
hdi(build)
quantile(build, probs=c(0.5,0.025,0.975))

mast <- c(ncomp.out1$chain1[,22],ncomp.out2$chain1[,22],ncomp.out3$chain1[,22],ncomp.out4$chain1[,22],ncomp.out5$chain1[,22],
          ncomp.out6$chain1[,22],ncomp.out7$chain1[,22],ncomp.out8$chain1[,22],ncomp.out9$chain1[,22],ncomp.out10$chain1[,22],
          ncomp.out1$chain2[,22],ncomp.out2$chain2[,22],ncomp.out3$chain2[,22],ncomp.out4$chain2[,22],ncomp.out5$chain2[,22],
          ncomp.out6$chain2[,22],ncomp.out7$chain2[,22],ncomp.out8$chain2[,22],ncomp.out9$chain2[,22],ncomp.out10$chain2[,22],
          ncomp.out1$chain3[,22],ncomp.out2$chain3[,22],ncomp.out3$chain3[,22],ncomp.out4$chain3[,22],ncomp.out5$chain3[,22],
          ncomp.out6$chain3[,22],ncomp.out7$chain3[,22],ncomp.out8$chain3[,22],ncomp.out9$chain3[,22],ncomp.out10$chain3[,22])

# Calculate HDI and quantiles
hdi(mast)
quantile(mast, probs=c(0.5,0.025,0.975))

sex1 <- c(ncomp.out1$chain1[,23],ncomp.out2$chain1[,23],ncomp.out3$chain1[,23],ncomp.out4$chain1[,23],ncomp.out5$chain1[,23],
          ncomp.out6$chain1[,23],ncomp.out7$chain1[,23],ncomp.out8$chain1[,23],ncomp.out9$chain1[,23],ncomp.out10$chain1[,23],
          ncomp.out1$chain2[,23],ncomp.out2$chain2[,23],ncomp.out3$chain2[,23],ncomp.out4$chain2[,23],ncomp.out5$chain2[,23],
          ncomp.out6$chain2[,23],ncomp.out7$chain2[,23],ncomp.out8$chain2[,23],ncomp.out9$chain2[,23],ncomp.out10$chain2[,23],
          ncomp.out1$chain3[,23],ncomp.out2$chain3[,23],ncomp.out3$chain3[,23],ncomp.out4$chain3[,23],ncomp.out5$chain3[,23],
          ncomp.out6$chain3[,23],ncomp.out7$chain3[,23],ncomp.out8$chain3[,23],ncomp.out9$chain3[,23],ncomp.out10$chain3[,23])

# Calculate HDI and quantiles
hdi(sex1)
quantile(sex1, probs=c(0.5,0.025,0.975))

sex2 <- c(ncomp.out1$chain1[,24],ncomp.out2$chain1[,24],ncomp.out3$chain1[,24],ncomp.out4$chain1[,24],ncomp.out5$chain1[,24],
          ncomp.out6$chain1[,24],ncomp.out7$chain1[,24],ncomp.out8$chain1[,24],ncomp.out9$chain1[,24],ncomp.out10$chain1[,24],
          ncomp.out1$chain2[,24],ncomp.out2$chain2[,24],ncomp.out3$chain2[,24],ncomp.out4$chain2[,24],ncomp.out5$chain2[,24],
          ncomp.out6$chain2[,24],ncomp.out7$chain2[,24],ncomp.out8$chain2[,24],ncomp.out9$chain2[,24],ncomp.out10$chain2[,24],
          ncomp.out1$chain3[,24],ncomp.out2$chain3[,24],ncomp.out3$chain3[,24],ncomp.out4$chain3[,24],ncomp.out5$chain3[,24],
          ncomp.out6$chain3[,24],ncomp.out7$chain3[,24],ncomp.out8$chain3[,24],ncomp.out9$chain3[,24],ncomp.out10$chain3[,24])

# Calculate HDI and quantiles
hdi(sex2)
quantile(sex2, probs=c(0.5,0.025,0.975))



#### 



