library(tidyverse)
library(nimble)
library(COMPoissonReg)


## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-edge_density_15, -edge_density_30, -edge_density_45, -build_cat_15, -build_cat_30,-build_cat_45) %>%
  mutate(mast_year=if_else(year==2019, 1, 0), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat[,c(8,16:42)] <- scale(dat[,c(8,16:42)])

## Set up data
numScaleVars <- 1
nScales <- 3
nonScaleVars <- 1
ncomp <- dat$ncomp
wmu <- as.numeric(factor(dat$WMU, labels=1:55))
ids <- as.numeric(factor(dat$RegionalID, labels=1:length(unique(dat$RegionalID))))

# create array for covariate data (slice for each covariate)
scale_covars <- array(NA, dim=c(nrow(dat), nScales, numScaleVars))
scale_covars[1:nrow(dat),1:3,1] <- as.matrix(dat[1:nrow(dat),28:30]) # number of buildings

covars <- matrix(NA, nrow=nrow(dat),ncol=nonScaleVars)
covars[1:nrow(dat),1] <- dat$lag_beechnuts

# prep fof nimble model
vsConstants <- list(N=nrow(dat),
                    sex=dat$Sex,
                    nsamples=length(unique(dat$RegionalID)), 
                    numVars=numScaleVars + nonScaleVars, 
                    numScaleVars=numScaleVars,
                    # nonScaleVars=nonScaleVars,
                    sampleID=ids) # Random intercepts

vsDataBundle <- list(ncomp=ncomp, # response
                     covars=covars, # covariates (no scale)
                     scale_covars=scale_covars,
                     age=dat$Age,
                     age2=dat$age2) # covariates (scale)

vsInits <- list( sigma.eta=1, mu.eta=1, nu=1.5, 
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
  
  # sample
  for (k in 1:nsamples) {
    eta[k] ~ dnorm(mu.eta, sd=sigma.eta)
  }
  mu.eta ~ dnorm(0, 0.001)
  sigma.eta ~ dunif(0, 100)
  
  ## Likelihood
  for (i in 1:N) {
    
    lambda[i] <- exp(beta_age*age[i] + beta_age2*age2[i] + beta_sex[sex[i]] + 
                       eta[sampleID[i]] + 
                       z[1]*beta[1]*covars[i,1] + 
                       z[2]*beta[2]*scale_covars[i, x_scale[1], 1])
    
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
system.time(samplesIndicator <- runMCMC(CMCMCIndicatorRJ, thin=1, niter=100000, nburnin=50000))

saveRDS(samplesIndicator, file = "results/ncomp_indicators.rds")
# samplesIndicator <- readRDS("results/ncomp_indicators.rds")

## Looking at results
par(mfrow = c(2, 1))
plot(samplesIndicator[,'beta[3]'], pch = 16, cex = 0.4, main = "beta[3] traceplot")
plot(samplesIndicator[,'z[3]'], pch = 16, cex = 0.4, main = "z[3] traceplot")

# Individual inclusion probabilities
par(mfrow = c(1, 1))
zCols <- grep("z\\[", colnames(samplesIndicator))
posterior_inclusion_prob <- colMeans(samplesIndicator[, zCols])
plot(1:3, posterior_inclusion_prob, ylim=c(0,1),
     xlab = "beta", ylab = "inclusion probability",
     main = "Inclusion probabilities for each beta")

## Plot scale probabilities
sCols <- grep("x_scale\\[", colnames(samplesIndicator))
posterior_scales <- samplesIndicator[, sCols]

posterior_scales <- as.data.frame(posterior_scales)

names(posterior_scales) <- c("pct_decid", "pct_evrgrn", "nbuildings", "stand_mean", "stand_sd")

posterior_scales <- posterior_scales %>% pivot_longer(1:5, names_to="covar", values_to="scale")


ggplot(posterior_scales) +
  geom_bar(aes(x=scale)) +
  facet_wrap(vars(covar))
  

