## Running models
## 2022-09-08

library(tidyverse)
library(MuMIn)
library(glmmTMB)
library(DHARMa)
library(broom.mixed)

options(scipen=999, digits=3)
set.seed(1)

#### Parallel processing ####
nt <- min(parallel::detectCores(),6)

#### Read in data ####
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

#### Analysis Set-up ####

# Categorical number of compounds (collapsing higher numbers)
dat$catNcompT <- ifelse(dat$n.compounds.T>=3, "3+", as.character(dat$n.compounds.T))
dat$catNcompT  <- ordered(dat$catNcompT , levels=c("0", "1", "2", "3+"))

# Make random effects factors
dat$WMU <- as.factor(dat$WMU)
dat$WMUA_code <- as.factor(dat$WMUA_code)
dat$year <- factor(dat$year)
dat$RegionalID <- factor(dat$RegionalID)

# Make age a categorical variable
dat$catAge[dat$Age>=2.5] <- "adult"
dat$catAge[dat$Age==1.5] <- "subadult"
dat$catAge[dat$Ag==0.5] <- "juvenile"

# Resort columns
dat <- dat[,c(1:10,25,11,13:23)] 

## Scale and center variables
dat[,c(10,16:23)] <- scale(dat[,c(10,16:23)])

## Percent AG
pctAG1 <- dat[, c(1:14, 16:18)]
pctAG1 <- distinct(pctAG1)
pctAG1 <- pctAG1 %>% group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

# Run models with glmmTMB
ag15 <- glmmTMB(n.compounds.T ~ totalag_15 + (1|WMU) + (1|year), data=pctAG1, 
                family=compois, control=glmmTMBControl(parallel=nt))
ag30 <- glmmTMB(n.compounds.T ~ totalag_30 + (1|WMU) + (1|year), data=pctAG1, 
                family=compois, control=glmmTMBControl(parallel=nt))
ag60 <- glmmTMB(n.compounds.T ~ totalag_60 + (1|WMU) + (1|year), data=pctAG1, 
                family=compois, control=glmmTMBControl(parallel=nt))

crop15 <- glmmTMB(n.compounds.T ~ crops_15 + (1|WMU) + (1|year), data=pctAG1, 
                  family=compois, control=glmmTMBControl(parallel=nt))
crop30 <- glmmTMB(n.compounds.T ~ crops_30 + (1|WMU) + (1|year), data=pctAG1, 
                  family=compois, control=glmmTMBControl(parallel=nt))
crop60 <- glmmTMB(n.compounds.T ~ crops_60 + (1|WMU) + (1|year), data=pctAG1, 
                  family=compois, control=glmmTMBControl(parallel=nt))

past15 <- glmmTMB(n.compounds.T ~ pasture_15 + (1|WMU) + (1|year), data=pctAG1, 
                  family=compois, control=glmmTMBControl(parallel=nt))
past30 <- glmmTMB(n.compounds.T ~ pasture_30 + (1|WMU) + (1|year), data=pctAG1, 
                  family=compois, control=glmmTMBControl(parallel=nt))
past60 <- glmmTMB(n.compounds.T ~ pasture_60 + (1|WMU) + (1|year), data=pctAG1, 
                  family=compois, control=glmmTMBControl(parallel=nt))

pctAG_sel <- model.sel(ag15, ag30, ag60, 
                       crop15, crop30, crop60,
                       past15, past30, past60)

pctAG_sel

## Beech basal area
baa1 <- dat[, c(1:14,22,23)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=c(BMI, laggedBMI), values_fn=unique) %>% as.data.frame()

bmi15 <- glmmTMB(n.compounds.T ~ BMI_15 + (1|WMU) + (1|year), data=baa1, 
                 family=compois, control=glmmTMBControl(parallel=nt))
bmi30 <- glmmTMB(n.compounds.T ~ BMI_30 + (1|WMU) + (1|year), data=baa1, 
                 family=compois, control=glmmTMBControl(parallel=nt))
bmi60 <- glmmTMB(n.compounds.T ~ BMI_60 + (1|WMU) + (1|year), data=baa1, 
                 family=compois, control=glmmTMBControl(parallel=nt))

lbmi15 <- glmmTMB(n.compounds.T ~ laggedBMI_15 + (1|WMU) + (1|year), data=baa1, 
                  family=compois, control=glmmTMBControl(parallel=nt))
lbmi30 <- glmmTMB(n.compounds.T ~ laggedBMI_30 + (1|WMU) + (1|year), data=baa1, 
                  family=compois, control=glmmTMBControl(parallel=nt))
lbmi60 <- glmmTMB(n.compounds.T ~ laggedBMI_60 + (1|WMU) + (1|year), data=baa1, 
                  family=compois, control=glmmTMBControl(parallel=nt))

baa_sel <- model.sel(bmi15,lbmi15, bmi30,lbmi30,bmi60,lbmi60)
baa_sel

## Wildland-urban interface
intermix1 <- dat[, c(1:15, 19)]
intermix1 <- unite(intermix1, "buffrad", 14:15, sep="_")
intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[14:22] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                             "mix_15_250", "mix_30_250", "mix_60_250",
                             "mix_15_500", "mix_30_500", "mix_60_500") 

interface1 <- dat[, c(1:15, 20)]
interface1 <- unite(interface1, "buffrad", 14:15, sep="_")
interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[14:22] <- c("face_15_100", "face_30_100", "face_60_100",
                              "face_15_250", "face_30_250", "face_60_250",
                              "face_15_500", "face_30_500", "face_60_500") 

wui1 <- dat[, c(1:15, 21)]
wui1 <- unite(wui1, "buffrad", 14:15, sep="_")
wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[14:22] <- c("wui_15_100", "wui_30_100", "wui_60_100",
                        "wui_15_250", "wui_30_250", "wui_60_250",
                        "wui_15_500", "wui_30_500", "wui_60_500") 

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 14:22)]
wui1 <- left_join(wui1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

# Join interface WUI
interface1 <- interface1[, c(1:3, 14:22)]
wui1 <- left_join(wui1, interface1, by=c("RegionalID", "pt_name", "pt_index"))

# Intermix WUI
mix_15100 <- glmmTMB(n.compounds.T ~ mix_15_100 + (1|WMU) + (1|year), data=wui1, 
                     family=compois, control=glmmTMBControl(parallel=nt))
mix_15250 <- glmmTMB(n.compounds.T ~ mix_15_250 + (1|WMU) + (1|year), data=wui1, 
                     family=compois, control=glmmTMBControl(parallel=nt))
mix_15500 <- glmmTMB(n.compounds.T ~ mix_15_500 + (1|WMU) + (1|year), data=wui1, 
                     family=compois, control=glmmTMBControl(parallel=nt))

mix_30100 <- glmmTMB(n.compounds.T ~ mix_30_100 + (1|WMU) + (1|year), data=wui1,
                     family=compois, control=glmmTMBControl(parallel=nt))
mix_30250 <- glmmTMB(n.compounds.T ~ mix_30_250 + (1|WMU) + (1|year), data=wui1,
                     family=compois, control=glmmTMBControl(parallel=nt))
mix_30500 <- glmmTMB(n.compounds.T ~ mix_30_500 + (1|WMU) + (1|year), data=wui1,
                     family=compois, control=glmmTMBControl(parallel=nt))

mix_60100 <- glmmTMB(n.compounds.T ~ mix_60_100 + (1|WMU) + (1|year), data=wui1,
                     family=compois, control=glmmTMBControl(parallel=nt))
mix_60250 <- glmmTMB(n.compounds.T ~ mix_60_250 + (1|WMU) + (1|year), data=wui1,
                     family=compois, control=glmmTMBControl(parallel=nt))
mix_60500 <- glmmTMB(n.compounds.T ~ mix_60_500 + (1|WMU) + (1|year), data=wui1,
                     family=compois, control=glmmTMBControl(parallel=nt))

# Interface WUI
face_15100 <- glmmTMB(n.compounds.T ~ face_15_100 + (1|WMU) + (1|year), data=wui1, 
                      family=compois, control=glmmTMBControl(parallel=nt))
face_15250 <- glmmTMB(n.compounds.T ~ face_15_250 + (1|WMU) + (1|year), data=wui1, 
                      family=compois, control=glmmTMBControl(parallel=nt))
face_15500 <- glmmTMB(n.compounds.T ~ face_15_500 + (1|WMU) + (1|year), data=wui1, 
                      family=compois, control=glmmTMBControl(parallel=nt))

face_30100 <- glmmTMB(n.compounds.T ~ face_30_100 + (1|WMU) + (1|year), data=wui1, 
                      family=compois, control=glmmTMBControl(parallel=nt))
face_30250 <- glmmTMB(n.compounds.T ~ face_30_250 + (1|WMU) + (1|year), data=wui1, 
                      family=compois, control=glmmTMBControl(parallel=nt))
face_30500 <- glmmTMB(n.compounds.T ~ face_30_500 + (1|WMU) + (1|year), data=wui1, 
                      family=compois, control=glmmTMBControl(parallel=nt))

face_60100 <- glmmTMB(n.compounds.T ~ face_60_100 + (1|WMU) + (1|year), data=wui1, 
                      family=compois, control=glmmTMBControl(parallel=nt))
face_60250 <- glmmTMB(n.compounds.T ~ face_60_250 + (1|WMU) + (1|year), data=wui1, 
                      family=compois, control=glmmTMBControl(parallel=nt))
face_60500 <- glmmTMB(n.compounds.T ~ face_60_500 + (1|WMU) + (1|year), data=wui1, 
                      family=compois, control=glmmTMBControl(parallel=nt))

# Total WUI
wui_15100 <- glmmTMB(n.compounds.T ~ wui_15_100 + (1|WMU) + (1|year), data=wui1, 
                     family=compois, control=glmmTMBControl(parallel=nt))
wui_15250 <- glmmTMB(n.compounds.T ~ wui_15_250 + (1|WMU) + (1|year), data=wui1, 
                     family=compois, control=glmmTMBControl(parallel=nt))
wui_15500 <- glmmTMB(n.compounds.T ~ wui_15_500 + (1|WMU) + (1|year), data=wui1, 
                     family=compois, control=glmmTMBControl(parallel=nt))

wui_30100 <- glmmTMB(n.compounds.T ~ wui_30_100 + (1|WMU) + (1|year), data=wui1, 
                     family=compois, control=glmmTMBControl(parallel=nt))
wui_30250 <- glmmTMB(n.compounds.T ~ wui_30_250 + (1|WMU) + (1|year), data=wui1, 
                     family=compois, control=glmmTMBControl(parallel=nt))
wui_30500 <- glmmTMB(n.compounds.T ~ wui_30_500 + (1|WMU) + (1|year), data=wui1, 
                     family=compois, control=glmmTMBControl(parallel=nt))

wui_60100 <- glmmTMB(n.compounds.T ~ wui_60_100 + (1|WMU) + (1|year), data=wui1, 
                     family=compois, control=glmmTMBControl(parallel=nt))
wui_60250 <- glmmTMB(n.compounds.T ~ wui_60_250 + (1|WMU) + (1|year), data=wui1, 
                     family=compois, control=glmmTMBControl(parallel=nt))
wui_60500 <- glmmTMB(n.compounds.T ~ wui_60_500 + (1|WMU) + (1|year), data=wui1, 
                     family=compois, control=glmmTMBControl(parallel=nt))

wui_sel <- model.sel(wui_15100, wui_30100, wui_60100, face_15100, face_30100, face_60100,mix_15100, mix_30100, mix_60100, 
                     wui_15250, wui_30250, wui_60250, face_15250, face_30250, face_60250,mix_15250, mix_30250, mix_60250, 
                     wui_15500, wui_30500, wui_60500, face_15500, face_30500, face_60500,mix_15500, mix_30500, mix_60500)
wui_sel

#### Set up data to run for each combination of covariates ####
dat1 <- dat[,c(1:13)]
dat1 <- distinct(dat1)

# Join percent agriculture
pctAG1 <- pctAG1[,c(1:3, 15)]
dat1 <- left_join(dat1, pctAG1, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 18)]
dat1 <- left_join(dat1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join total WUI
wui1 <- wui1[,c(1:3, 24)]
dat1 <- left_join(dat1, wui1, by=c("RegionalID", "pt_name", "pt_index"))
# write_csv(dat1, "output/model_data.csv")

## Check correlation matrix
cor(dat1[,14:16])

#### Running final models ####


library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = m1_pt, re.form = NULL)
testDispersion(simulationOutput, alternative = "less", plot = FALSE) # only underdispersion


## Loop over each set of random points
m_est <- m_stderr <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- dat1[dat1$pt_index==i,]
  
  # Run model with deltaAICc < 2

  m1_pt <- glmmTMB(n.compounds.T ~ Sex*Age + pasture_30 + mix_30_100 + 
                      laggedBMI_30 + (1|WMU) + (1|year), data=pt, 
                      family=compois(link = "log"), control=glmmTMBControl(parallel=nt))

  m1s <- summary(m1_pt)

  # save averaged confidence intervals
  pct2.5 <- rbind(pct2.5, t(confint(m1_pt))[1,1:7])
  pct97.5 <- rbind(pct97.5, t(confint(m1_pt))[2,1:7])

  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)$cond[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)$cond[,2])
  m_ranef <- rbind(m_ranef, unlist(VarCorr(m1_pt)))

  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- c(row.names(m1s$coefficients$cond))
    names(m_stderr) <- c(row.names(m1s$coefficients$cond))
    names(pct2.5) <- c(row.names(m1s$coefficients$cond))
    names(pct97.5) <- c(row.names(m1s$coefficients$cond))
    names(m_ranef) <- c("RE_WMU", "RE_year")
  }
  
}

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)
ranef_avg <- as.data.frame(colMeans(m_ranef[sapply(m_ranef, is.numeric)])) %>% rownames_to_column("RE")
names(ranef_avg)[2] <- "variance"

mean(m_est[,2])
sd(m_est[,2])

# One-sample t-test to determine "significance"
pvalue <- c()
for (i in 1:ncol(m_est)) {
  
  tresult <- t.test(m_est[,i], mu=0, alternative="two.sided")
  pvalue <- c(pvalue, tresult$p.value)
}
pvalue <- as.data.frame(pvalue)
pvalue <- cbind(names(coef_avg), pvalue)
pvalue <- pivot_wider(pvalue, names_from="names(coef_avg)", values_from="pvalue") %>% as.data.frame()

# Combine and clean up data frame
coef_summary <- bind_rows(coef_avg, stderr_avg, pct2.5_avg, pct97.5_avg, pvalue)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "2.5CI", "97.5CI", "P-value")
coef_summary <- data.frame(coef=coefs, coef_summary)
names(coef_summary)[2] <- "Intercept"

# Write to file
write_csv(coef_summary, "results/ncompT_coef-summary.csv")
write_csv(ranef_avg, "results/ncompT_coef-random-effects.csv")






