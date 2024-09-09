library(tidyverse)
library(nimble)

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-build_cat_15, -build_cat_30,-build_cat_45) %>%
  mutate(mast_year=if_else(year==2019, 2, 1), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat[,c(8,16:(ncol(dat)-1))] <- scale(dat[,c(8,16:(ncol(dat)-1))])

# read data for individual compounds
diph <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="diphacinone") %>%
  select(RegionalID,bin.exp)
dat <- left_join(dat, diph, by="RegionalID")
dat <- dat %>% select(RegionalID:Town,bin.exp,deciduous_15:mast_year)

## Set up data
numScaleVars <- 10
nScales <- 3
nonScaleVars <- 4
diph <- dat$bin.exp
wmua <- as.numeric(factor(dat$WMUA_code, labels=1:18))
ids <- as.numeric(factor(dat$RegionalID, labels=1:length(unique(dat$RegionalID))))

# create array for covariate data (slice for each covariate)
scale_covars <- array(NA, dim=c(nrow(dat), nScales, numScaleVars))
scale_covars[1:nrow(dat),1:3,1] <- as.matrix(dat[1:nrow(dat),34:36]) # edge density
scale_covars[1:nrow(dat),1:3,2] <- as.matrix(dat[1:nrow(dat),28:30]) # number of buildings
scale_covars[1:nrow(dat),1:3,3] <- as.matrix(dat[1:nrow(dat),37:39]) # stand age mean
scale_covars[1:nrow(dat),1:3,4] <- as.matrix(dat[1:nrow(dat),40:42]) # stand age standard deviation
scale_covars[1:nrow(dat),1:3,5] <- as.matrix(dat[1:nrow(dat),43:45]) # LaurentianAcadianNorthernHardwoods
scale_covars[1:nrow(dat),1:3,6] <- as.matrix(dat[1:nrow(dat),46:48]) # LaurentianAcadianPinesHemlock
scale_covars[1:nrow(dat),1:3,7] <- as.matrix(dat[1:nrow(dat),49:51]) # NortheasternNATemperatePlantation
scale_covars[1:nrow(dat),1:3,8] <- as.matrix(dat[1:nrow(dat),52:54]) # NorthCentralRuderalForest
scale_covars[1:nrow(dat),1:3,9] <- as.matrix(dat[1:nrow(dat),55:57]) # AcadianLowElevationSpruceFir
scale_covars[1:nrow(dat),1:3,10] <- as.matrix(dat[1:nrow(dat),58:60]) # AppalachianHardwoodsHemlocks

## prep fof nimble model
vsConstants <- list(N=nrow(dat),
                    sex=dat$Sex,
                    nsamples=length(unique(dat$RegionalID)), 
                    numScaleVars=numScaleVars,
                    sampleID=ids) # random intercept

vsDataBundle <- list(y=diph, # response
                     beechnuts=dat$lag_beechnuts, 
                     scale_covars=scale_covars,
                     age=dat$Age,
                     age2=dat$age2) # covariates (scale)

vsInits <- list(sigma.eta=1, mu.eta=1, eta=rnorm(length(unique(ids))), 
                beta=rnorm(vsConstants$numScaleVars), 
                x_scale=rep(1, numScaleVars),
                beta_age=rnorm(1), beta_age2=rnorm(1), beta_sex=rnorm(2), beta_mast=rnorm(1),
                z=sample(0:1,(vsConstants$numScaleVars), 0.5), cat_prob=rep(1/3,3)) 


# Build model in BUGS language
var_scale_code <- nimbleCode({
  
  ## Priors
  
  # beta coefficient priors
  beta_age ~ dnorm(0, sd=1.4)
  beta_age2 ~ dnorm(0, sd=1.4)
  for (k in 1:2) {
    beta_sex[k] ~ dnorm(0, sd=1.4)
  }
  
  # Indicator betas
  for (k in 1:numScaleVars) {
    beta[k] ~ dnorm(0, sd=1.4)
    z[k] ~ dbern(0.5) # indicator for each coefficient
  }
  
  # Scale variables
  cat_prob[1:3] <- c(1/3, 1/3, 1/3)
  for (k in 1:numScaleVars) {
    x_scale[k] ~ dcat(cat_prob[1:3])
  }
  
  # random intercept for sample
  for (k in 1:nsamples) {
    eta[k] ~ dnorm(mu.eta, sd=sigma.eta)
  }
  mu.eta ~ dnorm(0, 0.001)
  sigma.eta ~ dunif(0, 100)
  
  ## Likelihood
  for (i in 1:N) {
    
    logit(p[i]) <- eta[sampleID[i]] + beta_age*age[i] + beta_age2*age2[i] + beta_sex[sex[i]] + beta_mast*beechnuts[i] +
      z[1]*beta[1]*scale_covars[i, x_scale[1], 1] + z[2]*beta[2]*scale_covars[i, x_scale[2], 2] +
      z[3]*beta[3]*scale_covars[i, x_scale[3], 3] + z[4]*beta[4]*scale_covars[i, x_scale[4], 4] +
      z[5]*beta[5]*scale_covars[i, x_scale[5], 5] + z[6]*beta[6]*scale_covars[i, x_scale[6], 6] +
      z[7]*beta[7]*scale_covars[i, x_scale[7], 7] + z[8]*beta[8]*scale_covars[i, x_scale[8], 8] +
      z[9]*beta[9]*scale_covars[i, x_scale[9], 9] + z[10]*beta[10]*scale_covars[i, x_scale[10], 10]
    
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
system.time(samplesIndicator <- runMCMC(CMCMCIndicatorRJ, niter=100000, nburnin=50000))

saveRDS(samplesIndicator, file = "results/diphacinone_indicators.rds")
# samplesIndicator <- readRDS("results/diphacinone_indicators.rds")

## Looking at results
par(mfrow = c(2, 1))
plot(samplesIndicator[,'beta[1]'], pch = 16, cex = 0.4, main = "beta[1] traceplot")
plot(samplesIndicator[,'z[1]'], pch = 16, cex = 0.4, main = "z[1] traceplot")

# Individual inclusion probabilities
par(mfrow = c(1, 1))
zCols <- grep("z\\[", colnames(samplesIndicator))
posterior_inclusion_prob <- colMeans(samplesIndicator[, zCols])
plot(1:9, posterior_inclusion_prob, ylim=c(0,1),
     xlab = "beta", ylab = "inclusion probability",
     main = "Inclusion probabilities for each beta")
posterior_inclusion_prob

## Plot scale probabilities
sCols <- grep("x_scale\\[", colnames(samplesIndicator))
posterior_scales <- samplesIndicator[, sCols]

posterior_scales <- as.data.frame(posterior_scales)

names(posterior_scales) <- c("pct_evergrn", "pct_decid", "nbuildings", "stand_mean", "stand_sd")

posterior_scales <- posterior_scales %>% pivot_longer(1:5, names_to="covar", values_to="scale")

ggplot(posterior_scales) +
  geom_bar(aes(x=scale)) +
  facet_wrap(vars(covar))
