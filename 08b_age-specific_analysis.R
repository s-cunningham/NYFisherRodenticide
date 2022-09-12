library(tidyverse)
library(MuMIn)
library(ordinal)

options(scipen=999, digits=3)
set.seed(123)

#### Read in data ####
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

#### Analysis Set-up ####

# Categorical number of compounds (collapsing higher numbers)
dat$catNcompMO <- ifelse(dat$n.compounds.MO>=2, "2+", as.character(dat$n.compounds.MO))
dat$catNcompMO <- ordered(dat$catNcompMO, levels=c("0", "1", "2+"))

# Order response
dat$n.compounds.MO <- ordered(dat$n.compounds.MO, levels=c(0,1,2,3))

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
dat <- dat[,c(1:10,26,11,24,25,14:23)] 

## Scale and center variables
dat[,c(10,17:24)] <- scale(dat[,c(10,17:24)])

#### Split data, set up covariates ####
subads <- dat %>% filter(catAge=="subadult")

## Percent AG
pctAG1 <- subads[, c(1:15, 17:19)]
pctAG1 <- distinct(pctAG1)
pctAG1 <- pctAG1 %>% group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

## Beech basal area
baa1 <- subads[, c(1:15,23,24)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=c(BMI, laggedBMI), values_fn=unique) %>% as.data.frame()

## Wildland-urban interface
intermix1 <- subads[, c(1:16, 20)] 
intermix1 <- intermix1 %>% unite("buffrad", 15:16, sep="_") %>% 
  group_by(RegionalID) %>% 
  pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[15:23] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                             "mix_15_250", "mix_30_250", "mix_60_250",
                             "mix_15_500", "mix_30_500", "mix_60_500") 

interface1 <- dat[, c(1:16, 21)]
interface1 <- interface1  %>% unite("buffrad", 15:16, sep="_") %>%
  group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[15:23] <- c("face_15_100", "face_30_100", "face_60_100",
                              "face_15_250", "face_30_250", "face_60_250",
                              "face_15_500", "face_30_500", "face_60_500") 

wui1 <- dat[, c(1:16, 22)]
wui1 <- wui1 %>% unite("buffrad", 15:16, sep="_") %>%
  group_by(RegionalID) %>% 
  pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[15:23] <- c("wui_15_100", "wui_30k_100", "wui_60_100",
                        "wui_15_250", "wui_30k_250", "wui_60_250",
                        "wui_15_500", "wui_30k_500", "wui_60_500") 

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 14:22)]
wui1 <- left_join(wui1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

# Join interface WUI
interface1 <- interface1[, c(1:3, 14:22)]
wui1 <- left_join(wui1, interface1, by=c("RegionalID", "pt_name", "pt_index"))

## Add covariates
subad1 <- subads %>% select(1:14) %>% distinct()

# Join percent agriculture
pctAG1 <- pctAG1[,c(1:3, 16)]
subad1 <- left_join(subad1, pctAG1, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 19)]
subad1 <- left_join(subad1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join total WUI
wui1 <- wui1[,c(1:3, 26)]
subad1 <- left_join(subad1, wui1, by=c("RegionalID", "pt_name", "pt_index"))

#### With trace ####
## Loop over each set of random points ##
m_est <- m_stderr <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- subad1 %>% filter(pt_index==i)
  
  # Run model with deltaAICc < 2
  m1_pt <- clmm(catNcompT ~ Sex + pasture_30 + mix_30_100 + laggedBMI_30 + (1|WMU) + (1|year), data=pt, na.action="na.fail")
  
  m1s <- summary.merMod(m1_pt)
  
  # save averaged confidence intervals
  ci <- confint.merMod(m1_pt, method="boot")
  pct2.5 <- rbind(pct2.5, t(ci)[1,2:10])
  pct97.5 <- rbind(pct97.5, t(ci)[2,2:10])
  
  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  m_ranef <- rbind(m_ranef, unlist(VarCorr(m1_pt)))
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- c(names(coef(m1_pt)))
    names(m_stderr) <- c(names(coef(m1_pt)))
    names(pct2.5) <- c(names(coef(m1_pt)))
    names(pct97.5) <- c(names(coef(m1_pt)))
    names(m_ranef) <- c("RE_WMU", "RE_year")
  }
  
}

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)
ranef_avg <- as.data.frame(colMeans(m_ranef[sapply(m_ranef, is.numeric)]))

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

# Write to file
write_csv(coef_summary, "results/subadsT.csv")
write_csv(ranef_avg, "results/subadsT-random-effects.csv")

#### Without trace ####
## Loop over each set of random points ##
m_est <- m_stderr <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- subad1 %>% filter(pt_index==i)
  
  # Run model with deltaAICc < 2
  m1_pt <- clmm(catNcompMO ~ Sex + pasture_30 + mix_30_100 + laggedBMI_30 + (1|WMU) + (1|year), data=pt, na.action="na.fail")
  
  m1s <- summary(m1_pt)
  
  # save averaged confidence intervals
  pct2.5 <- rbind(pct2.5, t(confint(m1_pt))[1,])
  pct97.5 <- rbind(pct97.5, t(confint(m1_pt))[2,])
  
  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  m_ranef <- rbind(m_ranef, unlist(VarCorr(m1_pt)))
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- c(names(coef(m1_pt)))
    names(m_stderr) <- c(names(coef(m1_pt)))
    names(pct2.5) <- c(names(coef(m1_pt)))
    names(pct97.5) <- c(names(coef(m1_pt)))
    names(m_ranef) <- c("RE_WMU", "RE_year")
  }
  
}

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)
ranef_avg <- as.data.frame(colMeans(m_ranef[sapply(m_ranef, is.numeric)]))

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

# Write to file
write_csv(coef_summary, "results/subadsMO.csv")
write_csv(ranef_avg, "results/subadsMO-random-effects.csv")




