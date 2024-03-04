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
bliss_nimble <- nimbleCode({
  
  ## Priors
  psi ~ dunif(0,1) # prior for inclusion probability
  nu ~ dunif(1,2) # prior for CMP dispersion parameter
  
  # beta coefficient priors
  for (k in 1:numVars) {
    z[k] ~ dbern(psi) # indicator for each coefficient
    beta[k] ~ dnorm(0, sd=100)
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
    
    lambda <- alpha[WMU[i]] + inprod()
    ncomp ~ dCOMP(lambda, nu)
    
  }
  
})

## Set up the model.
lmIndicatorConstants <- list(N=100, numVars=15, nWMU=55)
lmIndicatorInits <- list(sigma.alpha=1, psi=0.5, nu=1.21, 
                         beta=rnorm(lmIndicatorConstants$numVars),
                         z=sample(0:1, lmIndicatorConstants$numVars, 0.5))

lmIndicatorData  <- list(y = y, X = X)
lmIndicatorModel <- nimbleModel(code=bliss_nimble, constants=lmIndicatorConstants,
                                inits=lmIndicatorInits, data=lmIndicatorData)

lmIndicatorConf <- configureMCMC(lmIndicatorModel)
lmIndicatorConf$addMonitors('z') # need to add other variables?

configureRJ(lmIndicatorConf,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2))

## Check the assigned samplers
lmIndicatorConf$printSamplers(c("z[1]", "beta[1]"))

## Build and run MCMC
mcmcIndicatorRJ <- buildMCMC(lmIndicatorConf)

cIndicatorModel <- compileNimble(lmIndicatorModel)

CMCMCIndicatorRJ <- compileNimble(mcmcIndicatorRJ, project = lmIndicatorModel)

# *probably in wrong spot*
# Individual inclusion proportion
# We can calculate the proportion of times each coefficient is included in the model.
betaCols <- grep("beta\\[", colnames(samplesNoIndicator))
posterior_inclusion_proportions <- colMeans(apply(samplesNoIndicator[, betaCols],
                                                  2, function(x) x != 0))
posterior_inclusion_proportions
