library(tidyverse)
library(nimble)

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-build_cat) %>%
  mutate(mast_year=if_else(year==2019, 2, 1), # failure years are reference
         Sex=if_else(Sex=='F',0,1)) # Females are reference

# Scale variables
dat[,c(8,17:41)] <- scale(dat[,c(8,17:41)])

# read data for individual compounds
brod <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Brodifacoum") %>%
  select(RegionalID,bin.exp)
dat <- left_join(dat, brod, by="RegionalID")
dat <- dat %>% select(RegionalID:Town,bin.exp,deciduous:mast_year)

## Set up data
numScaleVars <- 4
nVars <- 6
brod <- dat$bin.exp
wmua <- as.numeric(factor(dat$WMUA_code, labels=1:18))
ids <- as.numeric(factor(dat$RegionalID, labels=1:length(unique(dat$RegionalID))))

# create array for covariate data (slice for each covariate)
scale_covars <- array(NA, dim=c(nrow(dat), 4, 3))
scale_covars[1:nrow(dat),1:4,1] <- as.matrix(dat[1:nrow(dat),36:39]) # mast
scale_covars[1:nrow(dat),1:4,2] <- as.matrix(dat[1:nrow(dat),c(20,33,34,35)]) # WUI
scale_covars[1:nrow(dat),1:4,3] <- as.matrix(dat[1:nrow(dat),22:25]) # forest structure

## prep fof nimble model
vsConstants <- list(N=nrow(dat),
                    sex=dat$Sex,
                    nsamples=length(unique(dat$RegionalID)), 
                    nWMUA=length(unique(wmua)),
                    nVars=nVars,
                    sampleID=ids,
                    WMUA=wmua) # random intercept

vsDataBundle <- list(y=brod, # response
                     scale_covars=scale_covars,
                     age=dat$Age,
                     age2=dat$age2) # covariates (scale)

set.seed(1)
vsInits <- list(sigma.eta=1, mu.eta=1, eta=rnorm(length(unique(ids))), 
                sigma.eps=1, mu.eps=1, eps=rnorm(length(unique(wmua))), 
                beta=rnorm(nVars), cat_prob=rep(1/4,4), 
                psi=1, mast_scale=1, wui_scale=1, 
                fstruct_scale=1) 

# Build model in BUGS language
var_scale_code <- nimbleCode({
  
  ## Priors
  psi ~ dgamma()   ## prior on inclusion probability

  # Indicator betas
  for (k in 1:nVars) {
    beta[k] ~ dnorm(0, sd=1.4)
    z[k] ~ dbern(psi) # indicator for each coefficient
  }
  
  # Scale variables
  cat_prob[1:4] <- c(1/4, 1/4, 1/4, 1/4)
  mast_scale ~ dcat(cat_prob[1:4])
  wui_scale ~ dcat(cat_prob[1:4])
  fstruct_scale ~ dcat(cat_prob[1:4])
  
  # random intercept for sample
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
    
    logit(p[i]) <- eta[sampleID[i]] + eps[WMUA[i]] +
                    z[1]*beta[1]*age[i] +
                    z[2]*beta[2]*age2[i] +
                    z[3]*beta[3]*sex[i] +
                    z[4]*beta[4]*scale_covars[i, mast_scale, 1] +
                    z[5]*beta[5]*scale_covars[i, wui_scale, 2] +
                    z[6]*beta[6]*scale_covars[i, fstruct_scale, 3] 
      
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
system.time(samplesIndicator <- runMCMC(CMCMCIndicatorRJ, niter=25000, nburnin=12000))

saveRDS(samplesIndicator, file = "results/brodifacoum_indicators.rds")
# samplesIndicator <- readRDS("results/brodifacoum_indicators.rds")

## Looking at results
par(mfrow = c(2, 1))
plot(samplesIndicator[,'beta[3]'], pch = 16, cex = 0.4, main = "beta[3] traceplot")
plot(samplesIndicator[,'z[3]'], pch = 16, cex = 0.4, main = "z[3] traceplot")

# Individual inclusion probabilities
par(mfrow = c(1, 1))
zCols <- grep("z\\[", colnames(samplesIndicator))
posterior_inclusion_prob <- colMeans(samplesIndicator[, zCols])
plot(1:6, posterior_inclusion_prob, ylim=c(0,1),
     xlab = "beta", ylab = "inclusion probability",
     main = "Inclusion probabilities for each beta")
abline(h=0.5)

## Plot scale probabilities
sCols <- grep("_scale", colnames(samplesIndicator))
posterior_scales <- samplesIndicator[, sCols]

posterior_scales <- as.data.frame(posterior_scales)

names(posterior_scales) <- c( "fstruct_vars", "mast_vars", "wui_vars")

posterior_scales <- posterior_scales %>% pivot_longer(1:4, names_to="covar", values_to="scale")

ggplot(posterior_scales) +
  geom_bar(aes(x=scale)) +
  facet_wrap(vars(covar))
