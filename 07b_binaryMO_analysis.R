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

# trace = 0 binary
dat$MObinary <- ifelse(dat$exposure=="measured", 1, 0)

## Reorder columns
dat <- dat[,c(1:10,27,11,24,26,28,14:23)]

## Scale and center variables
dat[,c(10,18:25)] <- scale(dat[,c(10,18:25)])

#### Analysis ####

# Subset by compound
brod <- dat[dat$compound=="Brodifacoum",]
brom <- dat[dat$compound=="Bromadiolone",]
diph <- dat[dat$compound=="Diphacinone",]

#### Brodifacoum ####
## Percent AG
pctAG_brod <- brod %>% select(1:3,6:8,15,16,18:20) %>% 
  distinct() %>% 
  group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

## Beech basal area
baa1 <- brod[, c(1:3,6:8,15,16,24,25)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=c(BMI, laggedBMI), values_fn=unique) %>% as.data.frame()

## Wildland-urban interface
# Intermix
intermix1 <- brod %>% select(1:3,6:8,15:17,21) %>%
  unite("buffrad", 8:9, sep="_") %>%
  group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[8:16] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                            "mix_15_250", "mix_30_250", "mix_60_250",
                            "mix_15_500", "mix_30_500", "mix_60_500")  

## Set up data to run for each combination of covariates ##
brod1 <- brod %>% select(1:15) %>% distinct()

# Join percent agriculture
pctAG_brod <- pctAG_brod[,c(1:3, 9)]
brod1 <- left_join(brod1, pctAG_brod, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 12)]
brod1 <- left_join(brod1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join WUI
intermix1 <- intermix1[,c(1:3, 9)]
brod1 <- left_join(brod1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

## Run models ##

## Check correlation matrix
cor(brod1[,16:18])

## Loop over each set of random points
m_est <- m_stderr <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- brod1[brod1$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- glm(MObinary ~ Sex*catAge + laggedBMI_30 + mix_30_100 + pasture_30, 
                 family=binomial(link="logit"), data=pt)
  m1s <- summary(m1_pt)
  m1sdf <- m1s$coefficients
  
  # save averaged confidence intervals
  ci <- confint(m1_pt) 
  pct2.5 <- rbind(pct2.5, t(ci)[1,1:nrow(ci)])
  pct97.5 <- rbind(pct97.5, t(ci)[2,1:nrow(ci)])
  rsquared <- r.squaredLR(m1_pt)  # MuMIn

  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- row.names(m1sdf)
    names(m_stderr) <- row.names(m1sdf)
    names(pct2.5) <- row.names(m1sdf)
    names(pct97.5) <- row.names(m1sdf)
  }
  
}

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)

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
names(coef_summary) <- row.names(m1sdf)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "2.5CI", "97.5CI", "P-value")
coef_summary <- data.frame(coef=coefs, coef_summary)

# Write to file
write_csv(coef_summary, "results/binaryMObrodifacoum_coef-summary.csv")

#### Bromadiolone ####
## Percent AG
pctAG_brom <- brom %>% select(1:3,6:8,15,16,18:20) %>% 
  distinct() %>% 
  group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

## Beech basal area
baa1 <- brom[, c(1:3,6:8,15,16,24,25)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=c(BMI, laggedBMI), values_fn=unique) %>% as.data.frame()

## Wildland-urban interface
# Intermix
intermix1 <- brom %>% select(1:3,6:8,15:17,21) %>%
  unite("buffrad", 8:9, sep="_") %>%
  group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[8:16] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                            "mix_15_250", "mix_30_250", "mix_60_250",
                            "mix_15_500", "mix_30_500", "mix_60_500")  

## Set up data to run for each combination of covariates ##
brom1 <- brom %>% select(1:15) %>% distinct()

# Join percent agriculture
pctAG_brom <- pctAG_brom[,c(1:3, 9)]
brom1 <- left_join(brom1, pctAG_brom, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 12)]
brom1 <- left_join(brom1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join WUI
intermix1 <- intermix1[,c(1:3, 9)]
brom1 <- left_join(brom1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

## Run models ##

## Check correlation matrix
cor(brom1[,16:18])

## Loop over each set of random points
m_est <- m_stderr <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- brom1[brom1$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- glmer(bin.exp ~ Sex*catAge + pasture_30 + mix_30_100 + laggedBMI_30 + 
                   (1|WMU) + (1|year), family=binomial(link="logit"),data=pt)
  
  m1s <- summary(m1_pt)
  m1sdf <- m1s$coefficients
  
  # save averaged confidence intervals
  ci <- confint.merMod(m1_pt, method="Wald") 
  ci <- ci[complete.cases(ci),]
  pct2.5 <- rbind(pct2.5, t(ci)[1,1:nrow(ci)])
  pct97.5 <- rbind(pct97.5, t(ci)[2,1:nrow(ci)])
  m_ranef <- rbind(m_ranef, unlist(VarCorr(m1_pt)))
  
  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  
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
write_csv(coef_summary, "results/binaryMObromadiolone_coef-summary.csv")


#### Diphacinone ####
## Percent AG
pctAG_diph <- diph %>% select(1:3,6:8,15,16,18:20) %>% 
  distinct() %>% 
  group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

## Beech basal area
baa1 <- diph[, c(1:3,6:8,15,16,24,25)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=c(BMI, laggedBMI), values_fn=unique) %>% as.data.frame()

## Wildland-urban interface
# Intermix
intermix1 <- diph %>% select(1:3,6:8,15:17,21) %>%
  unite("buffrad", 8:9, sep="_") %>%
  group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[8:16] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                            "mix_15_250", "mix_30_250", "mix_60_250",
                            "mix_15_500", "mix_30_500", "mix_60_500")  

## Set up data to run for each combination of covariates ##
diph1 <- diph %>% select(1:15) %>% distinct()

# Join percent agriculture
pctAG_diph <- pctAG_diph[,c(1:3, 9)]
diph1 <- left_join(diph1, pctAG_diph, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 12)]
diph1 <- left_join(diph1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join WUI
intermix1 <- intermix1[,c(1:3, 9)]
diph1 <- left_join(diph1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

## Check correlation matrix
cor(diph1[,16:18])

## Loop over each set of random points
m_est <- m_stderr <- pct2.5 <- pct97.5 <- m_rsq <- data.frame()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- diph1[diph1$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- glmer(MObinary ~ Sex*catAge + pasture_30 + mix_30_100 + laggedBMI_30 + (1|WMUA_code),
                 family=binomial(link="logit"),data=pt)
  # null_mod <- glm(MObinary ~ 1 + (1|WMUA_code), family=binomial(link="logit"),data=pt)
  # rsq <- 1-logLik(m1_pt)/logLik(null_mod)
  
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
  # m_rsq <- rbind(m_rsq, unlist(rsq))
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- row.names(m1sdf)
    names(m_stderr) <- row.names(m1sdf)
    names(pct2.5) <- row.names(m1sdf)
    names(pct97.5) <- row.names(m1sdf)
    # names(m_rsq) <- "pseudoRsq"
  }
  
}

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)
# rsq_avg <- mean(m_rsq$pseudoRsq)

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
write_csv(coef_summary, "results/binaryMOdiphacinone_coef-summary.csv")






