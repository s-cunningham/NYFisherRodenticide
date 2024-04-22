library(tidyverse)
library(nimble)
library(MCMCvis)

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-edge_density_15, -edge_density_30, -edge_density_45, -build_cat_15, -build_cat_30,-build_cat_45) %>%
  mutate(mast_year=if_else(year==2019, 2, 1), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat[,c(8,16:42)] <- scale(dat[,c(8,16:42)])

# read data for individual compounds
diph <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Diphacinone") %>%
            select(RegionalID,bin.exp)
dat <- left_join(dat, diph, by="RegionalID")
dat <- dat %>% select(RegionalID:Town,bin.exp,deciduous_15:mast_year)

## Set up data
numScaleVars <- 5
nScales <- 3
nonScaleVars <- 3
diph <- dat$bin.exp
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

## prep fof nimble model
vsConstants <- list(N=nrow(dat),
                    sex=dat$Sex,
                    nsamples=length(unique(dat$RegionalID)), 
                    numVars=nonScaleVars, 
                    numScaleVars=numScaleVars,
                    WMU=wmu, # random intercept
                    sampleID=ids, # random intercept
                    nWMU=length(unique(dat$WMU)))

vsDataBundle <- list(y=diph, # response
                     covars=covars, # covariates (no scale)
                     scale_covars=scale_covars,
                     age=dat$Age,
                     age2=dat$age2) # covariates (scale)

vsInits <- list(sigma.alpha=1, mu.alpha=1, sigma.eta=1, mu.eta=1,
                beta=rnorm(vsConstants$numVars), sc_beta=rnorm(vsConstants$numScaleVars), 
                x_scale=rep(1, numScaleVars),
                beta_age=rnorm(1), beta_age2=rnorm(1), beta_sex=rnorm(2), #beta_mast=rnorm(2),
                z=sample(0:1,(vsConstants$numVars), 0.5), cat_prob=rep(1/3,3)) 



# Build model in BUGS language
var_scale_code <- nimbleCode({
  
  ## Priors
  
  # beta coefficient priors
  beta_age ~ dnorm(0, sd=10)
  beta_age2 ~ dnorm(0, sd=10)
  for (k in 1:2) {
    beta_sex[k] ~ dnorm(0, sd=10)
  }
  
  # Indicator betas
  for (k in 1:numVars) {
    beta[k] ~ dnorm(0, sd=10)
    z[k] ~ dbern(0.5) # indicator for each coefficient
    zbeta[k] <- z[k]*beta[k]
  }
  
  # Scale betas
  for (k in 1:numScaleVars) {
    sc_beta[k] ~ dnorm(0, sd=10)
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
    
    logit(p[i]) <- beta_age*age[i] + beta_age2*age2[i] + beta_sex[sex[i]] + alpha[WMU[i]] + eta[sampleID[i]] +
                    zbeta[1]*covars[i,1] + zbeta[2]*covars[i,2] + zbeta[3]*covars[i,3] +
                    sc_beta[1]*scale_covars[i, x_scale[1], 1] + sc_beta[2]*scale_covars[i, x_scale[2], 2] +
                    sc_beta[3]*scale_covars[i, x_scale[3], 3] + sc_beta[4]*scale_covars[i, x_scale[4], 4] +
                    sc_beta[5]*scale_covars[i, x_scale[5], 5] 
    
    y[i] ~ dbern(p[i])
    
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
system.time(samplesIndicator <- runMCMC(CMCMCIndicatorRJ, niter=10000, nburnin=4000))

saveRDS(samplesIndicator, file = "results/diphacinone_indicators.rds")

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
