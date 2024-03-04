## see code:
# https://figshare.com/articles/dataset/A_Bayesian_approach_for_multiscale_modeling_of_the_influence_of_seasonal_and_annual_habitat_variation_on_relative_abundance_of_ring-necked_pheasant_roosters/21901740

library(tidyverse)
library(nimble)
library(COMPoissonReg)
library(MCMCvis)

## Read data
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")

nWMU <- length(unique(dat$WMU))

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
  psi ~ dunif(0,1) # prior for inclusion probability
  nu ~ dunif(1,2) # prior for CMP dispersion parameter
  
  # beta coefficient priors
  for (k in 1:numVars) {
    z[k] ~ dbern(psi) # indicator for each coefficient
    beta[k] ~ dnorm(0, sd=10)
    zbeta[k] <- z[k] * beta[k]
  }
  
  # random intercept
  for (k in 1:nWMU) {
    alpha[k] ~ dnorm(mu.alpha, sd=sigma.alpha)
  }
  mu.alpha ~ dnorm(0, 0.001)
  sigma.alpha ~ dunif(0, 100)
  
  ## Likelihood
  for (i in 1:N) {
    
    log(lambda[i]) <- alpha[WMU[i]] + inprod()
    ncomp ~ dCOMP(lambda[i], nu)
    
  }
  
})

## Set up the model.
vsIndicatorConstants <- list(N=100, numVars=15, nWMU=55)
vsIndicatorInits <- list(sigma.alpha=1, psi=0.5, nu=1.5, 
                         beta=rnorm(vsIndicatorConstants$numVars),
                         z=sample(0:1, vsIndicatorConstants$numVars, 0.5))

vsIndicatorData  <- list(y = y, X = X)
vsIndicatorModel <- nimbleModel(code=var_scale_code, constants=vsIndicatorConstants,
                                inits=vsIndicatorInits, data=vsIndicatorData)

vsIndicatorConf <- configureMCMC(vsIndicatorModel)
vsIndicatorConf$addMonitors('z') # need to add other variables?

configureRJ(vsIndicatorConf,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2))

## Check the assigned samplers
vsIndicatorConf$printSamplers(c("z[1]", "beta[1]"))

## Build and run MCMC
mcmcIndicatorRJ <- buildMCMC(vsIndicatorConf)

cIndicatorModel <- compileNimble(vsIndicatorModel)

CMCMCIndicatorRJ <- compileNimble(mcmcIndicatorRJ, project = vsIndicatorModel)
