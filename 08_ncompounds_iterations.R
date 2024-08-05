## Nimble model after variable and scale selection
library(tidyverse)
library(nimble)
library(COMPoissonReg)
library(MCMCvis)
library(HDInterval)


## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-edge_density_15, -edge_density_30, -edge_density_45, -build_cat_15, -build_cat_30,-build_cat_45) %>%
  mutate(mast_year=if_else(year==2019, 1, 2), # mast years are reference
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
  beta0 ~ dnorm(0, sd=10)
  beta_age ~ dnorm(0, sd=10)
  beta_age2 ~ dnorm(0, sd=10)
  beta_build ~ dnorm(0, sd=10)
  beta_sex[1] <- 0
  beta_sex[2] ~ dnorm(0, sd=10)
  beta_mast[1] <- 0
  beta_mast[2] ~ dnorm(0, sd=10)

  ## random intercepts
  # WMU
  for (k in 1:nWMU) {
    alpha[k] ~ dnorm(mu.alpha, sd=sigma.alpha)
  }
  mu.alpha ~ dnorm(0, 0.01)
  sigma.alpha ~ dunif(0, 100)
  
  ## Likelihood
  for (i in 1:N) {
    
    log(lambda[i]) <- beta0 + beta_age*age[i] + beta_age2*age2[i] + beta_sex[sex[i]] + beta_build*covars[i,1] + beta_mast[mast[i]] + alpha[WMU[i]] 
    
    ncomp[i] ~ dCOMP(lambda[i], nu)
    
  }
  
})

# parameters to monitor
params <- c("beta0","beta_age","beta_age2","beta_mast","beta_sex","beta_build",
            "nu", "alpha", "mu.alpha", "sigma.alpha")  

# MCMC options
nt <- 1
ni <- 50000
nb <- 20000
nc <- 3

set.seed(1)
Inits <- list(sigma.alpha=1, mu.alpha=1, nu=1.5, beta0=rnorm(1), beta_age=1,beta_build=1,beta_mast=rep(1,2), beta_age2=1)
              # 
              #  beta_sex=rep(1,2),

#### Loop over random locations ####
## Iteration 1
dat1 <- dat %>% filter(pt_index==1)

## Set up data
ncomp1 <- dat1$ncomp
wmu1 <- as.numeric(factor(dat1$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars1 <- matrix(NA, nrow=nrow(dat1),ncol=1)
covars1[1:nrow(dat1),1] <- dat1$nbuildings_15

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
system.time(
  ncomp.out1 <- nimbleMCMC(code=ncompounds_code, constants=Constants1, data=DataBundle1, 
                           inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                           nchains=nc, check=FALSE, samples=TRUE,
                           samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)
)


ncomp.sum1 <- MCMCsummary(ncomp.out1)
ncomp.sum1 <- rownames_to_column(ncomp.sum1, "parameter")
range(ncomp.sum1$Rhat, na.rm=TRUE)

MCMCtrace(ncomp.out1, 
          params = c('beta_mast[1]', 'beta_mast[2]'), 
          ISB = FALSE, 
          exact = TRUE,
          pdf = FALSE)

## Iteration 2
dat2 <- dat %>% filter(pt_index==2)

## Set up data
ncomp2 <- dat2$ncomp
wmu2 <- as.numeric(factor(dat2$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars2 <- matrix(NA, nrow=nrow(dat2),ncol=1)
covars2[1:nrow(dat2),1] <- dat2$nbuildings_15

## prep fof nimble model
Constants2 <- list(N=nrow(dat2),
                   sex=dat2$Sex,
                   mast=dat2$mast_year,
                   WMU=wmu2, # random intercept
                   nWMU=length(unique(dat2$WMU)))

DataBundle2 <- list(ncomp=ncomp2, # response
                    covars=covars2, # covariates 
                    age=dat2$Age,
                    age2=dat2$age2) 

set.seed(1)
ncomp.out2 <- nimbleMCMC(code=ncompounds_code, constants=Constants2, data=DataBundle2, 
                         inits=Inits, monitors=params, thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum2 <- MCMCsummary(ncomp.out2)
ncomp.sum2 <- rownames_to_column(ncomp.sum2, "parameter")
range(ncomp.sum2$Rhat)

## Iteration 3
dat3 <- dat %>% filter(pt_index==3)

## Set up data
ncomp3 <- dat3$ncomp
wmu3 <- as.numeric(factor(dat3$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars3 <- matrix(NA, nrow=nrow(dat3),ncol=1)
covars3[1:nrow(dat3),1] <- dat3$nbuildings_15

## prep fof nimble model
Constants3 <- list(N=nrow(dat3),
                   sex=dat3$Sex,
                   mast=dat3$mast_year,
                   WMU=wmu3, # random intercept
                   nWMU=length(unique(dat3$WMU)))

DataBundle3 <- list(ncomp=ncomp3, # response
                    covars=covars3, # covariates 
                    age=dat3$Age,
                    age2=dat3$age2) 

set.seed(1)
ncomp.out3 <- nimbleMCMC(code=ncompounds_code, constants=Constants3, data=DataBundle3, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum3 <- MCMCsummary(ncomp.out3)
ncomp.sum3 <- rownames_to_column(ncomp.sum3, "parameter")
range(ncomp.sum3$Rhat)

## Iteration 4
dat4 <- dat %>% filter(pt_index==4)

## Set up data
ncomp4 <- dat4$ncomp
wmu4 <- as.numeric(factor(dat4$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars4 <- matrix(NA, nrow=nrow(dat4),ncol=1)
covars4[1:nrow(dat4),1] <- dat4$nbuildings_15

## prep fof nimble model
Constants4 <- list(N=nrow(dat4),
                   sex=dat4$Sex,
                   mast=dat4$mast_year,
                   WMU=wmu4, # random intercept
                   nWMU=length(unique(dat4$WMU)))

DataBundle4 <- list(ncomp=ncomp4, # response
                    covars=covars4, # covariates 
                    age=dat4$Age,
                    age2=dat4$age2) 

set.seed(1)
ncomp.out4 <- nimbleMCMC(code=ncompounds_code, constants=Constants4, data=DataBundle4, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum4 <- MCMCsummary(ncomp.out4)
ncomp.sum4 <- rownames_to_column(ncomp.sum4, "parameter")
range(ncomp.sum4$Rhat)

## Iteration 5
dat5 <- dat %>% filter(pt_index==5)

## Set up data
ncomp5 <- dat5$ncomp
wmu5 <- as.numeric(factor(dat5$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars5 <- matrix(NA, nrow=nrow(dat5),ncol=1)
covars5[1:nrow(dat5),1] <- dat5$nbuildings_15

## prep fof nimble model
Constants5 <- list(N=nrow(dat5),
                   sex=dat5$Sex,
                   mast=dat5$mast_year,
                   WMU=wmu5, # random intercept
                   nWMU=length(unique(dat5$WMU)))

DataBundle5 <- list(ncomp=ncomp5, # response
                    covars=covars5, # covariates 
                    age=dat5$Age,
                    age2=dat5$age2) 

set.seed(1)
ncomp.out5 <- nimbleMCMC(code=ncompounds_code, constants=Constants5, data=DataBundle5, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum5 <- MCMCsummary(ncomp.out5)
ncomp.sum5 <- rownames_to_column(ncomp.sum5, "parameter")
range(ncomp.sum5$Rhat)

## Iteration 6
dat6 <- dat %>% filter(pt_index==6)

## Set up data
ncomp6 <- dat6$ncomp
wmu6 <- as.numeric(factor(dat6$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars6 <- matrix(NA, nrow=nrow(dat6),ncol=1)
covars6[1:nrow(dat6),1] <- dat6$nbuildings_15

## prep fof nimble model
Constants6 <- list(N=nrow(dat6),
                   sex=dat6$Sex,
                   mast=dat6$mast_year,
                   WMU=wmu6, # random intercept
                   nWMU=length(unique(dat6$WMU)))

DataBundle6 <- list(ncomp=ncomp6, # response
                    covars=covars6, # covariates 
                    age=dat6$Age,
                    age2=dat6$age2) 

set.seed(1)
ncomp.out6 <- nimbleMCMC(code=ncompounds_code, constants=Constants6, data=DataBundle6, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum6 <- MCMCsummary(ncomp.out6)
ncomp.sum6 <- rownames_to_column(ncomp.sum6, "parameter")
range(ncomp.sum6$Rhat)

## Iteration 7
dat7 <- dat %>% filter(pt_index==7)

## Set up data
ncomp7 <- dat7$ncomp
wmu7 <- as.numeric(factor(dat7$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars7 <- matrix(NA, nrow=nrow(dat7),ncol=1)
covars7[1:nrow(dat7),1] <- dat7$nbuildings_15

## prep fof nimble model
Constants7 <- list(N=nrow(dat7),
                   sex=dat7$Sex,
                   mast=dat7$mast_year,
                   WMU=wmu7, # random intercept
                   nWMU=length(unique(dat7$WMU)))

DataBundle7 <- list(ncomp=ncomp7, # response
                    covars=covars7, # covariates 
                    age=dat7$Age,
                    age2=dat7$age2) 

set.seed(1)
ncomp.out7 <- nimbleMCMC(code=ncompounds_code, constants=Constants7, data=DataBundle7, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum7 <- MCMCsummary(ncomp.out7)
ncomp.sum7 <- rownames_to_column(ncomp.sum7, "parameter")
range(ncomp.sum7$Rhat)

## Iteration 8
dat8 <- dat %>% filter(pt_index==8)

## Set up data
ncomp8 <- dat8$ncomp
wmu8 <- as.numeric(factor(dat8$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars8 <- matrix(NA, nrow=nrow(dat8),ncol=1)
covars8[1:nrow(dat8),1] <- dat8$nbuildings_15

## prep fof nimble model
Constants8 <- list(N=nrow(dat8),
                   sex=dat8$Sex,
                   mast=dat8$mast_year,
                   WMU=wmu8, # random intercept
                   nWMU=length(unique(dat8$WMU)))

DataBundle8 <- list(ncomp=ncomp8, # response
                    covars=covars8, # covariates 
                    age=dat8$Age,
                    age2=dat8$age2) 

set.seed(1)
ncomp.out8 <- nimbleMCMC(code=ncompounds_code, constants=Constants8, data=DataBundle8, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum8 <- MCMCsummary(ncomp.out8)
ncomp.sum8 <- rownames_to_column(ncomp.sum8, "parameter")
range(ncomp.sum8$Rhat)

## Iteration 9
dat9 <- dat %>% filter(pt_index==9)

## Set up data
ncomp9 <- dat9$ncomp
wmu9 <- as.numeric(factor(dat9$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars9 <- matrix(NA, nrow=nrow(dat9),ncol=1)
covars9[1:nrow(dat9),1] <- dat9$nbuildings_15

## prep fof nimble model
Constants9 <- list(N=nrow(dat9),
                   sex=dat9$Sex,
                   mast=dat9$mast_year,
                   WMU=wmu9, # random intercept
                   nWMU=length(unique(dat9$WMU)))

DataBundle9 <- list(ncomp=ncomp9, # response
                    covars=covars9, # covariates 
                    age=dat9$Age,
                    age2=dat9$age2) 

set.seed(1)
ncomp.out9 <- nimbleMCMC(code=ncompounds_code, constants=Constants9, data=DataBundle9, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum9 <- MCMCsummary(ncomp.out9)
ncomp.sum9 <- rownames_to_column(ncomp.sum9, "parameter")
range(ncomp.sum9$Rhat)

## Iteration 10
dat10 <- dat %>% filter(pt_index==10)

## Set up data
ncomp10 <- dat1$ncomp
wmu10 <- as.numeric(factor(dat10$WMU, labels=1:55))

# create array for covariate data (column for each covariate)
covars10 <- matrix(NA, nrow=nrow(dat10),ncol=1)
covars10[1:nrow(dat10),1] <- dat10$nbuildings_15

## prep fof nimble model
Constants10 <- list(N=nrow(dat10),
                   sex=dat10$Sex,
                   mast=dat10$mast_year,
                   WMU=wmu10, # random intercept
                   nWMU=length(unique(dat10$WMU)))

DataBundle10 <- list(ncomp=ncomp10, # response
                    covars=covars10, # covariates 
                    age=dat10$Age,
                    age2=dat10$age2) 

set.seed(1)
ncomp.out10 <- nimbleMCMC(code=ncompounds_code, constants=Constants10, data=DataBundle10, 
                         inits=Inits, monitors=params,thin=nt, niter=ni, nburnin=nb, 
                         nchains=nc, check=FALSE, samples=TRUE,
                         samplesAsCodaMCMC=TRUE, summary=FALSE, WAIC=FALSE)

ncomp.sum10 <- MCMCsummary(ncomp.out10)
ncomp.sum10 <- rownames_to_column(ncomp.sum10, "parameter")
range(ncomp.sum10$Rhat)


# Combine samples (still need to check column numbers)

age <- c(ncomp.out1$chain1[,56],ncomp.out2$chain1[,56],ncomp.out3$chain1[,56],ncomp.out4$chain1[,56],ncomp.out5$chain1[,56],
         ncomp.out6$chain1[,56],ncomp.out7$chain1[,56],ncomp.out8$chain1[,56],ncomp.out9$chain1[,56],ncomp.out10$chain1[,56],
         ncomp.out1$chain2[,56],ncomp.out2$chain2[,56],ncomp.out3$chain2[,56],ncomp.out4$chain2[,56],ncomp.out5$chain2[,56],
         ncomp.out6$chain2[,56],ncomp.out7$chain2[,56],ncomp.out8$chain2[,56],ncomp.out9$chain2[,56],ncomp.out10$chain2[,56],
         ncomp.out1$chain3[,56],ncomp.out2$chain3[,56],ncomp.out3$chain3[,56],ncomp.out4$chain3[,56],ncomp.out5$chain3[,56],
         ncomp.out6$chain3[,56],ncomp.out7$chain3[,56],ncomp.out8$chain3[,56],ncomp.out9$chain3[,56],ncomp.out10$chain3[,56])

# Calculate HDI and quantiles
hdi(age)
quantile(age, probs=c(0.5,0.025,0.975))


age2 <- c(ncomp.out1$chain1[,57],ncomp.out2$chain1[,57],ncomp.out3$chain1[,57],ncomp.out4$chain1[,57],ncomp.out5$chain1[,57],
          ncomp.out6$chain1[,57],ncomp.out7$chain1[,57],ncomp.out8$chain1[,57],ncomp.out9$chain1[,57],ncomp.out10$chain1[,57],
          ncomp.out1$chain2[,57],ncomp.out2$chain2[,57],ncomp.out3$chain2[,57],ncomp.out4$chain2[,57],ncomp.out5$chain2[,57],
          ncomp.out6$chain2[,57],ncomp.out7$chain2[,57],ncomp.out8$chain2[,57],ncomp.out9$chain2[,57],ncomp.out10$chain2[,57],
          ncomp.out1$chain3[,57],ncomp.out2$chain3[,57],ncomp.out3$chain3[,57],ncomp.out4$chain3[,57],ncomp.out5$chain3[,57],
          ncomp.out6$chain3[,57],ncomp.out7$chain3[,57],ncomp.out8$chain3[,57],ncomp.out9$chain3[,57],ncomp.out10$chain3[,57])

# Calculate HDI and quantiles
hdi(age2)
quantile(age2, probs=c(0.5,0.025,0.975))

build <- c(ncomp.out1$chain1[,58],ncomp.out2$chain1[,58],ncomp.out3$chain1[,58],ncomp.out4$chain1[,58],ncomp.out5$chain1[,58],
           ncomp.out6$chain1[,58],ncomp.out7$chain1[,58],ncomp.out8$chain1[,58],ncomp.out9$chain1[,58],ncomp.out10$chain1[,58],
           ncomp.out1$chain2[,58],ncomp.out2$chain2[,58],ncomp.out3$chain2[,58],ncomp.out4$chain2[,58],ncomp.out5$chain2[,58],
           ncomp.out6$chain2[,58],ncomp.out7$chain2[,58],ncomp.out8$chain2[,58],ncomp.out9$chain2[,58],ncomp.out10$chain2[,58],
           ncomp.out1$chain3[,58],ncomp.out2$chain3[,58],ncomp.out3$chain3[,58],ncomp.out4$chain3[,58],ncomp.out5$chain3[,58],
           ncomp.out6$chain3[,58],ncomp.out7$chain3[,58],ncomp.out8$chain3[,58],ncomp.out9$chain3[,58],ncomp.out10$chain3[,58])

# Calculate HDI and quantiles
hdi(build)
quantile(build, probs=c(0.5,0.025,0.975))

decid <- c(ncomp.out1$chain1[,59],ncomp.out2$chain1[,59],ncomp.out3$chain1[,59],ncomp.out4$chain1[,59],ncomp.out5$chain1[,59],
           ncomp.out6$chain1[,59],ncomp.out7$chain1[,59],ncomp.out8$chain1[,59],ncomp.out9$chain1[,59],ncomp.out10$chain1[,59],
           ncomp.out1$chain2[,59],ncomp.out2$chain2[,59],ncomp.out3$chain2[,59],ncomp.out4$chain2[,59],ncomp.out5$chain2[,59],
           ncomp.out6$chain2[,59],ncomp.out7$chain2[,59],ncomp.out8$chain2[,59],ncomp.out9$chain2[,59],ncomp.out10$chain2[,59],
           ncomp.out1$chain3[,59],ncomp.out2$chain3[,59],ncomp.out3$chain3[,59],ncomp.out4$chain3[,59],ncomp.out5$chain3[,59],
           ncomp.out6$chain3[,59],ncomp.out7$chain3[,59],ncomp.out8$chain3[,59],ncomp.out9$chain3[,59],ncomp.out10$chain3[,59])

# Calculate HDI and quantiles
hdi(decid)
quantile(decid, probs=c(0.5,0.025,0.975))

evrgrn <- c(ncomp.out1$chain1[,60],ncomp.out2$chain1[,60],ncomp.out3$chain1[,60],ncomp.out4$chain1[,60],ncomp.out5$chain1[,60],
            ncomp.out6$chain1[,60],ncomp.out7$chain1[,60],ncomp.out8$chain1[,60],ncomp.out9$chain1[,60],ncomp.out10$chain1[,60],
            ncomp.out1$chain2[,60],ncomp.out2$chain2[,60],ncomp.out3$chain2[,60],ncomp.out4$chain2[,60],ncomp.out5$chain2[,60],
            ncomp.out6$chain2[,60],ncomp.out7$chain2[,60],ncomp.out8$chain2[,60],ncomp.out9$chain2[,60],ncomp.out10$chain2[,60],
            ncomp.out1$chain3[,60],ncomp.out2$chain3[,60],ncomp.out3$chain3[,60],ncomp.out4$chain3[,60],ncomp.out5$chain3[,60],
            ncomp.out6$chain3[,60],ncomp.out7$chain3[,60],ncomp.out8$chain3[,60],ncomp.out9$chain3[,60],ncomp.out10$chain3[,60])

# Calculate HDI and quantiles
hdi(evrgrn)
quantile(evrgrn, probs=c(0.5,0.025,0.975))

mast <- c(ncomp.out1$chain1[,61],ncomp.out2$chain1[,61],ncomp.out3$chain1[,61],ncomp.out4$chain1[,61],ncomp.out5$chain1[,61],
          ncomp.out6$chain1[,61],ncomp.out7$chain1[,61],ncomp.out8$chain1[,61],ncomp.out9$chain1[,61],ncomp.out10$chain1[,61],
          ncomp.out1$chain2[,61],ncomp.out2$chain2[,61],ncomp.out3$chain2[,61],ncomp.out4$chain2[,61],ncomp.out5$chain2[,61],
          ncomp.out6$chain2[,61],ncomp.out7$chain2[,61],ncomp.out8$chain2[,61],ncomp.out9$chain2[,61],ncomp.out10$chain2[,61],
          ncomp.out1$chain3[,61],ncomp.out2$chain3[,61],ncomp.out3$chain3[,61],ncomp.out4$chain3[,61],ncomp.out5$chain3[,61],
          ncomp.out6$chain3[,61],ncomp.out7$chain3[,61],ncomp.out8$chain3[,61],ncomp.out9$chain3[,61],ncomp.out10$chain3[,61])

# Calculate HDI and quantiles
hdi(mast)
quantile(mast, probs=c(0.5,0.025,0.975))

sex1 <- c(ncomp.out1$chain1[,62],ncomp.out2$chain1[,62],ncomp.out3$chain1[,62],ncomp.out4$chain1[,62],ncomp.out5$chain1[,62],
          ncomp.out6$chain1[,62],ncomp.out7$chain1[,62],ncomp.out8$chain1[,62],ncomp.out9$chain1[,62],ncomp.out10$chain1[,62],
          ncomp.out1$chain2[,62],ncomp.out2$chain2[,62],ncomp.out3$chain2[,62],ncomp.out4$chain2[,62],ncomp.out5$chain2[,62],
          ncomp.out6$chain2[,62],ncomp.out7$chain2[,62],ncomp.out8$chain2[,62],ncomp.out9$chain2[,62],ncomp.out10$chain2[,62],
          ncomp.out1$chain3[,62],ncomp.out2$chain3[,62],ncomp.out3$chain3[,62],ncomp.out4$chain3[,62],ncomp.out5$chain3[,62],
          ncomp.out6$chain3[,62],ncomp.out7$chain3[,62],ncomp.out8$chain3[,62],ncomp.out9$chain3[,62],ncomp.out10$chain3[,62])

# Calculate HDI and quantiles
hdi(sex1)
quantile(sex1, probs=c(0.5,0.025,0.975))

sex2 <- c(ncomp.out1$chain1[,63],ncomp.out2$chain1[,63],ncomp.out3$chain1[,63],ncomp.out4$chain1[,63],ncomp.out5$chain1[,63],
          ncomp.out6$chain1[,63],ncomp.out7$chain1[,63],ncomp.out8$chain1[,63],ncomp.out9$chain1[,63],ncomp.out10$chain1[,63],
          ncomp.out1$chain2[,63],ncomp.out2$chain2[,63],ncomp.out3$chain2[,63],ncomp.out4$chain2[,63],ncomp.out5$chain2[,63],
          ncomp.out6$chain2[,63],ncomp.out7$chain2[,63],ncomp.out8$chain2[,63],ncomp.out9$chain2[,63],ncomp.out10$chain2[,63],
          ncomp.out1$chain3[,63],ncomp.out2$chain3[,63],ncomp.out3$chain3[,63],ncomp.out4$chain3[,63],ncomp.out5$chain3[,63],
          ncomp.out6$chain3[,63],ncomp.out7$chain3[,63],ncomp.out8$chain3[,63],ncomp.out9$chain3[,63],ncomp.out10$chain3[,63])

# Calculate HDI and quantiles
hdi(sex2)
quantile(sex2, probs=c(0.5,0.025,0.975))

standmn <- c(ncomp.out1$chain1[,64],ncomp.out2$chain1[,64],ncomp.out3$chain1[,64],ncomp.out4$chain1[,64],ncomp.out5$chain1[,64],
             ncomp.out6$chain1[,64],ncomp.out7$chain1[,64],ncomp.out8$chain1[,64],ncomp.out9$chain1[,64],ncomp.out10$chain1[,64],
             ncomp.out1$chain2[,64],ncomp.out2$chain2[,64],ncomp.out3$chain2[,64],ncomp.out4$chain2[,64],ncomp.out5$chain2[,64],
             ncomp.out6$chain2[,64],ncomp.out7$chain2[,64],ncomp.out8$chain2[,64],ncomp.out9$chain2[,64],ncomp.out10$chain2[,64],
             ncomp.out1$chain3[,64],ncomp.out2$chain3[,64],ncomp.out3$chain3[,64],ncomp.out4$chain3[,64],ncomp.out5$chain3[,64],
             ncomp.out6$chain3[,64],ncomp.out7$chain3[,64],ncomp.out8$chain3[,64],ncomp.out9$chain3[,64],ncomp.out10$chain3[,64])

# Calculate HDI and quantiles
hdi(standmn)
quantile(standmn, probs=c(0.5,0.025,0.975))

standsd <- c(ncomp.out1$chain1[,65],ncomp.out2$chain1[,65],ncomp.out3$chain1[,65],ncomp.out4$chain1[,65],ncomp.out5$chain1[,65],
             ncomp.out6$chain1[,65],ncomp.out7$chain1[,65],ncomp.out8$chain1[,65],ncomp.out9$chain1[,65],ncomp.out10$chain1[,65],
             ncomp.out1$chain2[,65],ncomp.out2$chain2[,65],ncomp.out3$chain2[,65],ncomp.out4$chain2[,65],ncomp.out5$chain2[,65],
             ncomp.out6$chain2[,65],ncomp.out7$chain2[,65],ncomp.out8$chain2[,65],ncomp.out9$chain2[,65],ncomp.out10$chain2[,65],
             ncomp.out1$chain3[,65],ncomp.out2$chain3[,65],ncomp.out3$chain3[,65],ncomp.out4$chain3[,65],ncomp.out5$chain3[,65],
             ncomp.out6$chain3[,65],ncomp.out7$chain3[,65],ncomp.out8$chain3[,65],ncomp.out9$chain3[,65],ncomp.out10$chain3[,65])

# Calculate HDI and quantiles
hdi(standsd)
quantile(standsd, probs=c(0.5,0.025,0.975))
