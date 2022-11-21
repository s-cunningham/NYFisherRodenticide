library(tidyverse)
library(lme4)
library(MuMIn)

options(scipen=999, digits=3)
set.seed(123)

#### Read in data ####
# Landscape covariates
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

# Individual compounds
dat2 <- read_csv("output/summarized_AR_results.csv")
dat2 <- as.data.frame(dat2)

dat2 <- dat2[dat2$compound=="Diphacinone" | dat2$compound=="Brodifacoum" |
             dat2$compound=="Bromadiolone", c(2,19,21:22) ]

dat <- left_join(dat, dat2, by="RegionalID")

# Make age a categorical variable
dat$catAge[dat$Age>=3.5] <- "adult"
dat$catAge[dat$Age==2.5] <- "subadult"
dat$catAge[dat$Ag<2.5] <- "juvenile"

# Make random effects factors
dat$WMU <- as.factor(dat$WMU)
dat$WMUA_code <- as.factor(dat$WMUA_code)
dat$year <- factor(dat$year)

## Reorder columns
dat <- dat[,c(1:10,27,11,24:26,14:23)]

## Scale and center variables
dat[,c(10,18:25)] <- scale(dat[,c(10,18:25)])

#### Analysis ####
## Use pooled data to determine scale

## Pe AG
pctAG <- dat[, c(1:3,6:8,15,16,18:20)]
pctAG <- distinct(pctAG)
pctAG <- pctAG %>% group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

## Beech basal area
baa1 <- dat[, c(1:3,6:8,15,16,24,25)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=c(BMI, laggedBMI), values_fn=unique) %>% as.data.frame()

## Wildland-urban interface
# Intermix
intermix1 <- dat[, c(1:3,6:8,15:17,21)]
intermix1 <- unite(intermix1, "buffrad", 8:9, sep="_") %>% group_by(RegionalID) %>% 
  pivot_wider(names_from=buffrad, values_from=intermix, values_fn=unique) %>% as.data.frame()
names(intermix1)[8:16] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                            "mix_15_250", "mix_30_250", "mix_60_250",
                            "mix_15_500", "mix_30_500", "mix_60_500")  

## Set up data to run for each combination of covariates ##
dat1 <- dat[,c(1:15)]
dat1 <- distinct(dat1)

# Join percent agriculture
pctAG <- pctAG[,c(1:3, 16)]
dat1 <- left_join(dat1, pctAG, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 12)]
dat1 <- left_join(dat1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join WUI
intermix1 <- intermix1[,c(1:3, 8)]
dat1 <- left_join(dat1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

# write_csv(dat1, "output/binary_model_data.csv")

# Subset by compound
brod <- dat1[dat1$compound=="Brodifacoum",]
brom <- dat1[dat1$compound=="Bromadiolone",]
diph <- dat1[dat1$compound=="Diphacinone",]

#### Brodifacoum ####

## Run models ##

## Loop over each set of random points
m_est <- m_stderr <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- brod[brod$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- glmer(bin.exp ~ Sex*Age + totalag_60 + mix_15_100 + 
                   laggedBMI_30 + (1|WMU) + (1|year), 
                  family=binomial(link="logit"), data=pt)
  
  
  
  
  
  
  m1s <- summary(m1_pt)
  m1sdf <- m1s$coefficients
  
  # save averaged confidence intervals
  ci <- confint.merMod(m1_pt, method="Wald")
  ci <- ci[complete.cases(ci),]
  pct2.5 <- rbind(pct2.5, t(ci)[1,1:nrow(ci)])
  pct97.5 <- rbind(pct97.5, t(ci)[2,1:nrow(ci)])
  
  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  m_ranef <- rbind(m_ranef, unlist(VarCorr(m1_pt)))
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- row.names(m1sdf)
    names(m_stderr) <- row.names(m1sdf)
    names(pct2.5) <- row.names(m1sdf)
    names(pct97.5) <- row.names(m1sdf)
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

# Combine and clean up data frame
coef_summary <- bind_rows(coef_avg, stderr_avg, pct2.5_avg, pct97.5_avg)
names(coef_summary) <- row.names(m1sdf)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "2.5CI", "97.5CI")
coef_summary <- data.frame(coef=coefs, coef_summary)

# Write to file
write_csv(coef_summary, "results/binaryTbrodifacoum_coef-summary.csv")
write_csv(ranef_avg, "results/binaryTbrodifacoum_coef-random-effects.csv")

#### Bromadiolone ####
## Running final models ##

## Loop over each set of random points
m_est <- m_stderr <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- brom[brom$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- glmer(bin.exp ~ Sex*Age + totalag_60 + mix_15_100 + laggedBMI_30 + 
                   (1|WMU) + (1|year), family=binomial(link="logit"),data=pt)
  
  m1s <- summary(m1_pt)
  m1sdf <- m1s$coefficients
  
  # save averaged confidence intervals
  ci <- confint.merMod(m1_pt, method="Wald")
  ci <- ci[complete.cases(ci),]
  pct2.5 <- rbind(pct2.5, t(ci)[1,1:nrow(ci)])
  pct97.5 <- rbind(pct97.5, t(ci)[2,1:nrow(ci)])
  
  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  m_ranef <- rbind(m_ranef, unlist(VarCorr(m1_pt)))
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- row.names(m1sdf)
    names(m_stderr) <- row.names(m1sdf)
    names(pct2.5) <- row.names(m1sdf)
    names(pct97.5) <- row.names(m1sdf)
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

# Combine and clean up data frame
coef_summary <- bind_rows(coef_avg, stderr_avg, pct2.5_avg, pct97.5_avg)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "2.5CI", "97.5CI")
coef_summary <- data.frame(coef=coefs, coef_summary)

# Write to file
write_csv(coef_summary, "results/binaryTbromadiolone_coef-summary.csv")
write_csv(ranef_avg, "results/binaryTbromadiolone_coef-random-effects.csv")

#### Diphacinone ####

## Loop over each set of random points
m_est <- m_stderr <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- diph[diph$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- glmer(bin.exp ~ Sex*Age + totalag_60 + mix_15_100 + laggedBMI_30 + (1|WMU), 
               family=binomial(link="logit"), data=pt)
  if (!isSingular(m1_pt)) {
    
    m1s <- summary(m1_pt)
    m1sdf <- m1s$coefficients
    
    # save averaged confidence intervals
    ci <- confint.merMod(m1_pt, method="Wald")
    ci <- ci[complete.cases(ci),]
    pct2.5 <- rbind(pct2.5, t(ci)[1,1:nrow(ci)])
    pct97.5 <- rbind(pct97.5, t(ci)[2,1:nrow(ci)])
    
    # Save point set estimates
    m_est <- rbind(m_est, coef(m1s)[,1])
    m_stderr <- rbind(m_stderr, coef(m1s)[,2])
    m_ranef <- rbind(m_ranef, unlist(VarCorr(m1_pt)))
    
    # Rename (only need to do once)
    if (i==1) {
      names(m_est) <- row.names(m1sdf)
      names(m_stderr) <- row.names(m1sdf)
      names(pct2.5) <- row.names(m1sdf)
      names(pct97.5) <- row.names(m1sdf)
      names(m_ranef) <- c("RE_year")
    }
    
  }
}

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)
ranef_avg <- as.data.frame(colMeans(m_ranef[sapply(m_ranef, is.numeric)])) %>% rownames_to_column("RE")
names(ranef_avg)[2] <- "variance"

# Combine and clean up data frame
coef_summary <- bind_rows(coef_avg, stderr_avg, pct2.5_avg, pct97.5_avg)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "2.5CI", "97.5CI")
coef_summary <- data.frame(coef=coefs, coef_summary)

# Write to file
write_csv(coef_summary, "results/binaryTdiphacinone_coef-summary.csv")
write_csv(ranef_avg, "results/binaryTdiphacinone_coef-random-effects.csv")

