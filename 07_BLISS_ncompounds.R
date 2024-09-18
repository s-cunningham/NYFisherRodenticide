library(tidyverse)
library(nimble)
library(COMPoissonReg)

# https://r-nimble.org/html_manual/cha-mcmc.html#sec:rjmcmc

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  select(-build_cat) %>%
  mutate(age2=Age^2) %>%
  mutate(mast_year=if_else(year==2019, 1, 2), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat[,c(8,17:(ncol(dat)-1))] <- scale(dat[,c(8,17:(ncol(dat)-1))])

## Set up data
nVars <- 3
ncomp <- dat$ncomp
wmua <- as.numeric(factor(dat$WMUA_code, labels=1:18))
ids <- as.numeric(factor(dat$RegionalID, labels=1:length(unique(dat$RegionalID))))

# create array for covariate data (slice for each covariate)
scale_covars <- array(NA, dim=c(nrow(dat), 4, 2))
scale_covars[1:nrow(dat),1:4,1] <- as.matrix(dat[1:nrow(dat),37:40]) # mast
scale_covars[1:nrow(dat),1:4,2] <- as.matrix(dat[1:nrow(dat),c(21,34:36)]) # WUI

scale_covars2 <- array(NA, dim=c(nrow(dat), 5, 1))
scale_covars2[1:nrow(dat),1:5,1] <- as.matrix(dat[1:nrow(dat),c(18,23:26)]) # forest structure

# prep fof nimble model
vsConstants <- list(N=nrow(dat),
                    sex=dat$Sex,
                    nsamples=length(unique(dat$RegionalID)), 
                    nVars=nVars,
                    nWMUA=length(unique(wmua)),
                    WMUA=wmua,
                    sampleID=ids) # Random intercepts

vsDataBundle <- list(ncomp=ncomp, # response
                     scale_covars=scale_covars,
                     scale_covars2=scale_covars2,
                     age=dat$Age,
                     age2=dat$age2) 

vsInits <- list(sigma.eta=1, mu.eta=1, eta=rnorm(length(unique(ids))), 
                sigma.eps=1, mu.eps=1, eps=rnorm(length(unique(wmua))), 
                nu=1.5, mast_scale=1, wui_scale=1, fstruct_scale=1,
                beta=rnorm(vsConstants$nVars), beta0=rnorm(1),
                beta_age=rnorm(1), beta_age2=rnorm(1),
                beta_sex=rnorm(2), cat_prob=rep(1/4,4), cat_prob2=rep(1/5,5)) 

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

  # beta coefficient priors
  beta_age ~ dnorm(0,  0.001)
  beta_age2 ~ dnorm(0,  0.001)
  beta_sex[1] <- 0
  beta_sex[2] ~ dnorm(0,  0.001)
  beta0 ~ dnorm(0,  0.001)
  
  for (j in 1:nVars) {
    beta[j] ~ dnorm(0,  0.001)
  }
  
  # Scale variables
  cat_prob[1:4] <- c(1/4, 1/4, 1/4, 1/4)
  mast_scale ~ dcat(cat_prob[1:4])
  wui_scale ~ dcat(cat_prob[1:4])
  cat_prob2[1:5] <- c(1/5, 1/5, 1/5, 1/5, 1/5)
  fstruct_scale ~ dcat(cat_prob2[1:5])
  
  # sample
  for (k in 1:nsamples) {
    eta[k] ~ dnorm(mu.eta, sd=sigma.eta)
  }
  mu.eta ~ dnorm(0, 0.001)
  sigma.eta ~ dunif(0, 100)
  
  # random intercept for WMUA
  for (j in 1:nWMUA) {
    eps[j] ~ dnorm(mu.eps, sd=sigma.eps)
  }
  mu.eps ~ dnorm(0, 0.001)
  sigma.eps ~ dunif(0, 100)
  
  ## Likelihood
  for (i in 1:N) {
    
    lambda[i] <- exp(eta[sampleID[i]] + eps[WMUA[i]] + beta0 +
                       beta_age*age[i] + beta_age2*age2[i] + beta_sex[sex[i]] + 
                       beta[1]*scale_covars[i, mast_scale, 1] +
                       beta[2]*scale_covars[i, wui_scale, 2] +
                       beta[3]*scale_covars2[i, fstruct_scale, 1])
    
    ncomp[i] ~ dCOMP(lambda[i], nu)
    
  }
  
})

params <- c("beta0", "beta_age", "beta_age2", "beta_sex", "beta", "eta", "eps", "mast_scale", "wui_scale", "fstruct_scale")

samples <- nimbleMCMC(
  code = var_scale_code,
  constants = vsConstants, 
  data =vsDataBundle, 
  inits = vsInits,
  monitors = params,
  niter = 10000,
  nburnin = 5000,
  thin = 1)


saveRDS(samples, file = "results/ncomp_indicators.rds")
# samplesIndicator <- readRDS("results/ncomp_indicators.rds")

## Looking at results
plot(samples[,'beta[1]'], type="l", cex = 0.4, main = "beta traceplot")

## Plot scale probabilities
sCols <- grep("_scale", colnames(samples))
posterior_scales <- samples[, sCols]

posterior_scales <- as.data.frame(posterior_scales)

names(posterior_scales) <- c("fstruct_scale", "mast_scale", "wui_scale")

posterior_scales <- posterior_scales %>% pivot_longer(1:3, names_to="covar", values_to="scale")


ggplot(posterior_scales) +
  geom_bar(aes(x=scale)) +
  facet_wrap(vars(covar))
  
posterior_scales %>% group_by(covar, scale) %>% count()
