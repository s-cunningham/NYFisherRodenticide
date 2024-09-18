library(tidyverse)
library(nimble)

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-build_cat) %>%
  mutate(mast_year=if_else(year==2019, 2, 1), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat[,c(8,17:41)] <- scale(dat[,c(8,17:41)])

# read data for individual compounds
diph <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Diphacinone") %>%
  select(RegionalID,bin.exp)
dat <- left_join(dat, diph, by="RegionalID")
dat <- dat %>% select(RegionalID:Town,bin.exp,deciduous:mast_year)

## Set up data
numScaleVars <- 4
nVars <- 6
diph <- dat$bin.exp
wmua <- as.numeric(factor(dat$WMUA_code, labels=1:18))
ids <- as.numeric(factor(dat$RegionalID, labels=1:length(unique(dat$RegionalID))))

# create array for covariate data (slice for each covariate)
scale_covars <- array(NA, dim=c(nrow(dat), 4, 2))
scale_covars[1:nrow(dat),1:4,1] <- as.matrix(dat[1:nrow(dat),36:39]) # mast
scale_covars[1:nrow(dat),1:4,2] <- as.matrix(dat[1:nrow(dat),c(20,33,34,35)]) # WUI

scale_covars2 <- array(NA, dim=c(nrow(dat), 5, 1))
scale_covars2[1:nrow(dat),1:5,1] <- as.matrix(dat[1:nrow(dat),c(17,22:25)]) # forest structure

## prep fof nimble model
vsConstants <- list(N=nrow(dat),
                    sex=dat$Sex,
                    nsamples=length(unique(dat$RegionalID)), 
                    nWMUA=length(unique(wmua)),
                    nVars=nVars,
                    sampleID=ids,
                    WMUA=wmua) # random intercept

vsDataBundle <- list(y=diph, # response
                     scale_covars=scale_covars,
                     scale_covars2=scale_covars2,
                     age=dat$Age,
                     age2=dat$age2) # covariates (scale)

set.seed(1)
vsInits <- list(sigma.eta=1, mu.eta=1, eta=rnorm(length(unique(ids))), 
                sigma.eps=1, mu.eps=1, eps=rnorm(length(unique(wmua))), 
                beta=rnorm(nVars), cat_prob=rep(1/4,4), beta_sex=rnorm(2),
                cat_prob2=rep(1/5,5), 
                mast_scale=1, wui_scale=1, 
                fstruct_scale=1) 

# Build model in BUGS language
var_scale_code <- nimbleCode({
  
  ## Priors
  # betas
  for (k in 1:(nVars-1)) {
    beta[k] ~ dnorm(0, sd=1.4)
  }
  beta_sex[1] <- 0
  beta_sex[2] ~ dnorm(0, sd=1.4)
  
  # Scale variables
  cat_prob[1:4] <- c(1/4, 1/4, 1/4, 1/4)
  mast_scale ~ dcat(cat_prob[1:4])
  wui_scale ~ dcat(cat_prob[1:4])
  
  cat_prob2[1:5] <- c(1/5, 1/5, 1/5, 1/5, 1/5)
  fstruct_scale ~ dcat(cat_prob2[1:5])
  
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
      beta_sex[sex[i]] +
      beta[1]*age[i] +
      beta[2]*age2[i] +
      beta[3]*scale_covars[i, mast_scale, 1] +
      beta[4]*scale_covars[i, wui_scale, 2] +
      beta[5]*scale_covars2[i, fstruct_scale, 1] 
    
    y[i] ~ dbern(p[i])
    
  }
  
})

params <- c("beta", "beta_sex", "eta", "eps", "mast_scale", "wui_scale", "fstruct_scale")

samples <- nimbleMCMC(
  code = var_scale_code,
  constants = vsConstants, 
  data =vsDataBundle, 
  inits = vsInits,
  monitors = params,
  niter = 12000,
  nburnin = 6000,
  thin = 1)

# Save samples
saveRDS(samples, file = "results/diphacinone_indicators.rds")

## Looking at results
plot(samples[,'beta[1]'], type="l", cex = 0.4, main = "beta[3] traceplot")

## Plot scale probabilities
sCols <- grep("_scale", colnames(samples))
posterior_scales <- samples[, sCols]

posterior_scales <- as.data.frame(posterior_scales)

names(posterior_scales) <- c( "fstruct_vars", "mast_vars", "wui_vars")

posterior_scales <- posterior_scales %>% pivot_longer(1:3, names_to="covar", values_to="scale")

ggplot(posterior_scales) +
  geom_bar(aes(x=scale)) +
  facet_wrap(vars(covar), scales="free_x")
