library(tidyverse)
library(nimble)
library(COMPoissonReg)
library(MCMCvis)

## Read data




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
  
  ## Priords
  psi ~ dunif(0,1) # prior for inclusion probability
  
  # beta coefficient priors
  for (k in 1:numVars) {
    z[k] ~ dbern(psi) # indicator for each coefficient
    beta[k] ~ dnorm(0, sd=100)
    zbeta[k] <- z[k] * beta[k]
  }
  
  # random intercept
  for (k in 1:nWMU) {
    
    
  }
  
  ## Likelihood
  for (i in 1:N) {
    
    
  }
  
  
})


