## Nimble model after variable and scale selection
library(tidyverse)
library(nimble)
library(COMPoissonReg)


## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-edge_density_15, -edge_density_30, -edge_density_45, -build_cat_15, -build_cat_30,-build_cat_45) %>%
  mutate(mast_year=if_else(year==2019, 2, 1), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat[,c(8,16:42)] <- scale(dat[,c(8,16:42)])

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
  beta_age ~ dnorm(0, sd=10)
  beta_age2 ~ dnorm(0, sd=10)
  for (k in 1:2) {
    beta_sex[k] ~ dnorm(0, sd=10)
  }
  for (k in 1:2) {
    beta_mast[k] ~ dnorm(0, sd=10)
  }
  
  beta_build ~ dnorm(0, sd=10)
  beta_decid ~ dnorm(0, sd=10)
  beta_evrgrn ~ dnorm(0, sd=10)
  beta_standm ~ dnorm(0, sd=10)
  beta_standsd ~ dnorm(0, sd=10)
  
  ## random intercepts
  # WMU
  for (k in 1:nWMU) {
    alpha[k] ~ dnorm(mu.alpha, sd=sigma.alpha)
  }
  mu.alpha ~ dnorm(0, 0.001)
  sigma.alpha ~ dunif(0, 100)
  
  ## Likelihood
  for (i in 1:N) {
    
    lambda[i] <- alpha[WMU[i]] + beta_age*age[i] + beta_age2*age2[i] + beta_sex[sex[i]] + 
                  beta_mast[mast[i]] + beta_decid*covars[i,1] + beta_evrgrn*covars[i,2] +
                  beta_build*covars[i,3] + beta_standm*covars[i,4] + beta_standsd*covars[i,5]
    
    ncomp[i] ~ dCOMP(lambda[i], nu)
    
  }
  
})

# parameters to monitor
params <- c("beta_age","beta_age2","beta_sex","beta_mast","beta_decid","beta_evrgrn",
            "beta_build","beta_standm", "beta_standsd", "nu", "alpha", "mu.alpha", "sigma.alpha")  

# MCMC options
nt <- 1
ni <- 75000
nb <- 35000
nc <- 3

set.seed(1)
Inits <- list(sigma.alpha=1, mu.alpha=1, nu=1.5,
              beta_mast=rnorm(2), beta_decid=rnorm(1), beta_evrgrn=rnorm(1), 
              beta_build=rnorm(1), beta_standm=rnorm(1), beta_standsd=rnorm(1), 
              beta_age=rnorm(1), beta_age2=rnorm(1), beta_sex=rnorm(2)) 

#### Loop over random locations ####
## Iteration 1
dat1 <- dat %>% filter(pt_index==1)

## Set up data
ncomp1 <- dat1$ncomp
wmu1 <- as.numeric(factor(dat1$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars1 <- matrix(NA, nrow=nrow(dat1),ncol=5)
covars1[1:nrow(dat1),1] <- dat1$deciduous_45
covars1[1:nrow(dat1),2] <- dat1$evergreen_45
covars1[1:nrow(dat1),3] <- dat1$nbuildings_45
covars1[1:nrow(dat1),4] <- dat1$stand_age_mean_45
covars1[1:nrow(dat1),5] <- dat1$stand_age_sd_45

## prep fof nimble model
Constants1 <- list(N=nrow(dat1),
                   sex=dat1$Sex,
                   mast=dat1$mast_year,
                   WMU=wmu1, # random intercept
                   nWMU=length(unique(dat1$WMU)))

DataBundle1 <- list(ncomp=ncomp1, # response
                    covars=covars1, # covariates 
                    age=dat1$Age,
                    age2=dat1$age2) 

set.seed(1)
ncomp.out1 <- nimbleMCMC(code=ncompounds_code, constants=Constants1, data=DataBundle1, 
                        inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                        nchains=nc, check=FALSE, samples=TRUE,
                        samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum1 <- MCMCsummary(ncomp.out1)
ncomp.sum1 <- rownames_to_column(ncomp.sum1, "parameter")
range(ncomp.sum1$Rhat)
