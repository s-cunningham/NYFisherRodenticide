library(tidyverse)
library(nimble)
library(COMPoissonReg)
library(MCMCvis)

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-edge_density_15, -edge_density_30, -edge_density_45, -build_cat_15, -build_cat_30,-build_cat_45) %>%
  mutate(mast_year=if_else(year==2019, 1, 0), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat[,c(8,16:42)] <- scale(dat[,c(8,16:42)])

## Set up data
numScaleVars <- 5
nScales <- 3
nonScaleVars <- 3
ncomp <- dat$ncomp
wmu <- as.numeric(factor(dat$WMU, labels=1:55))
ids <- as.numeric(factor(dat$RegionalID, labels=1:length(unique(dat$RegionalID))))

# create array for covariate data (slice for each covariate)
scale_covars <- array(NA, dim=c(nrow(dat), nScales, numScaleVars))
scale_covars[1:nrow(dat),1:3,1] <- as.matrix(dat[1:nrow(dat),16:18]) # % deciduous
scale_covars[1:nrow(dat),1:3,2] <- as.matrix(dat[1:nrow(dat),19:21]) # % evergreen
scale_covars[1:nrow(dat),1:3,3] <- as.matrix(dat[1:nrow(dat),28:30]) # number of buildings
scale_covars[1:nrow(dat),1:3,4] <- as.matrix(dat[1:nrow(dat),34:36]) # stand age mean
scale_covars[1:nrow(dat),1:3,5] <- as.matrix(dat[1:nrow(dat),37:39]) # stand age standard deviation

covars <- matrix(NA, nrow=nrow(dat),ncol=nonScaleVars)
covars[1:nrow(dat),1] <- dat$beechnuts
covars[1:nrow(dat),2] <- dat$lag_beechnuts
covars[1:nrow(dat),3] <- as.numeric(dat$mast_year)

# prep fof nimble model
vsConstants <- list(N=nrow(dat),
                    sex=dat$Sex,
                    nsamples=length(unique(dat$RegionalID)), 
                    numVars=numScaleVars + nonScaleVars, 
                    numScaleVars=numScaleVars,
                    WMU=wmu, # random intercept
                    sampleID=ids, # random intercept
                    nWMU=length(unique(dat$WMU)))

vsDataBundle <- list(ncomp=ncomp, # response
                     covars=covars, # covariates (no scale)
                     scale_covars=scale_covars,
                     age=dat$Age,
                     age2=dat$age2) # covariates (scale)

vsInits <- list(sigma.alpha=1, mu.alpha=1, sigma.eta=1, mu.eta=1, nu=1.5, 
                beta=rnorm(vsConstants$numVars), 
                x_scale=rep(1, numScaleVars),
                beta_age=rnorm(1), beta_age2=rnorm(1), beta_sex=rnorm(2), #beta_mast=rnorm(2),
                z=sample(0:1,(vsConstants$numVars), 0.5), cat_prob=rep(1/3,3)) 

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
var_scale_code <- nimbleCode({
  
  ## Priors
  nu ~ dunif(1,2.5) # prior for CMP dispersion parameter
  # V ~ dgamma(3.29, 7.8) # total beta variance
  # 
  # beta_var <- V/numVars
  
  # beta coefficient priors
  beta_age ~ dnorm(0, sd=10)
  beta_age2 ~ dnorm(0, sd=10)
  for (k in 1:2) {
    beta_sex[k] ~ dnorm(0, sd=10)
  }
  
  for (k in 1:numVars) {
    beta[k] ~ dnorm(0, sd=10)
    z[k] ~ dbern(0.5) # indicator for each coefficient
  }
  
  cat_prob[1:3] <- c(1/3, 1/3, 1/3)
  for (k in 1:numScaleVars) {
    x_scale[k] ~ dcat(cat_prob[1:3])
  }
  
  ## random intercepts
  # WMU
  for (k in 1:nWMU) {
    alpha[k] ~ dnorm(mu.alpha, sd=sigma.alpha)
  }
  mu.alpha ~ dnorm(0, 0.001)
  sigma.alpha ~ dunif(0, 100)
  
  # sample
  for (k in 1:nsamples) {
    eta[k] ~ dnorm(mu.eta, sd=sigma.eta)
  }
  mu.eta ~ dnorm(0, 0.001)
  sigma.eta ~ dunif(0, 100)
  
  ## Likelihood
  for (i in 1:N) {
    
    temp1[i] <- beta_age*age[i] + beta_age2*age2[i] + beta_sex[sex[i]] + alpha[WMU[i]] + eta[sampleID[i]] +
      z[1]*beta[1]*covars[i,1] + z[2]*beta[2]*covars[i,2] + z[3]*beta[3]*covars[i,3] 
    
    temp2[i] <- beta[4]*scale_covars[i, x_scale[1], 1] + beta[5]*scale_covars[i, x_scale[2], 2] +
                beta[6]*scale_covars[i, x_scale[3], 3] + beta[7]*scale_covars[i, x_scale[4], 4] +
                beta[8]*scale_covars[i, x_scale[5], 5] 
    
    lambda[i] <- exp(temp1[i] + temp2[i])
    
    ncomp[i] ~ dCOMP(lambda[i], nu)
    
  }
  
})

## Set up the model.
vsModel <- nimbleModel(code=var_scale_code, constants=vsConstants,
                       inits=vsInits, data=vsDataBundle)

vsIndicatorConf <- configureMCMC(vsModel)
vsIndicatorConf$addMonitors('z')

configureRJ(vsIndicatorConf,
            targetNodes = 'beta',
            indicatorNodes = 'z',  
            control = list(mean = 0, scale = .2))

## Check the assigned samplers
vsIndicatorConf$printSamplers(c("z[1]", "beta[1]"))

## Build and run MCMC
mcmcIndicatorRJ <- buildMCMC(vsIndicatorConf)

cIndicatorModel <- compileNimble(vsModel)

CMCMCIndicatorRJ <- compileNimble(mcmcIndicatorRJ, project = vsModel)

set.seed(1)
system.time(samplesIndicator <- runMCMC(CMCMCIndicatorRJ, niter=1000, nburnin=500))