library(tidyverse)
library(lme4)
library(MuMIn)

options(scipen=999, digits=3)
set.seed(123)

#### Read in data ####
# Landscape covariates
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")

# Individual compounds
dat2 <- read_csv("output/summarized_AR_results.csv") %>%
          filter(compound=="Diphacinone" | compound=="Brodifacoum" | compound=="Bromadiolone") %>%
          select(RegionalID,compound,bin.exp,bin.exp.ntr)
dat <- left_join(dat, dat2, by="RegionalID") %>%
        rename(catAge=AgeClass)

# Make random effects factors
dat$WMU <- as.factor(dat$WMU)
dat$WMUA_code <- as.factor(dat$WMUA_code)
dat$year <- factor(dat$year)
dat$RegionalID <- factor(dat$RegionalID)

## Reorder columns
dat <- dat %>% select(RegionalID:Town,compound,bin.exp,bin.exp.ntr,buffsize:laggedBMI)

## Use pooled data to determine scale

## Percent AG
pctAG <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, pasture, crops, totalag) %>% 
  distinct() %>% 
  group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag))

## Beech basal area
baa1 <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, BMI, laggedBMI) %>% 
  distinct() %>% 
  group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(BMI, laggedBMI), values_fn=unique) 

## Percent forest
pctFOR <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, deciduous, evergreen, mixed, totalforest) %>% 
  distinct() %>% 
  group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(deciduous, evergreen, mixed, totalforest)) %>% 
  as.data.frame()

## Number of buildings
build1 <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, nbuildings, build_cat) %>%
  distinct() %>% 
  group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(nbuildings, build_cat), values_fn=unique) %>% 
  as.data.frame()

## Wildland-urban interface
intermix1 <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, radius, intermix) %>%
  unite("buffrad", 4:5, sep="_") %>% 
  group_by(RegionalID) %>% distinct() %>%
  pivot_wider(names_from=buffrad, values_from=intermix) 
names(intermix1)[4:12] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                            "mix_15_250", "mix_30_250", "mix_60_250",
                            "mix_15_500", "mix_30_500", "mix_60_500") 

interface1  <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, radius, interface) %>%
  unite("buffrad", 4:5, sep="_") %>% 
  group_by(RegionalID) %>% distinct() %>%
  pivot_wider(names_from=buffrad, values_from=interface) 
names(interface1)[4:12] <- c("face_15_100", "face_30_100", "face_60_100",
                             "face_15_250", "face_30_250", "face_60_250",
                             "face_15_500", "face_30_500", "face_60_500") 

wui1  <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, radius, totalWUI) %>%
  unite("buffrad", 4:5, sep="_") %>% 
  group_by(RegionalID) %>% distinct() %>%
  pivot_wider(names_from=buffrad, values_from=totalWUI)
names(wui1)[4:12] <- c("wui_15_100", "wui_30_100", "wui_60_100",
                       "wui_15_250", "wui_30_250", "wui_60_250",
                       "wui_15_500", "wui_30_500", "wui_60_500") 

## Landscape metrics
lsm1 <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, ai:shape_mn) %>% 
  distinct() %>% 
  group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(ai:shape_mn), values_fn=unique) 

#### Set up data to run for each combination of covariates ####
dat1 <- dat %>% select(RegionalID:bin.exp) %>% distinct() %>%
  left_join(pctAG, by=c("RegionalID", "pt_name", "pt_index")) %>%
  left_join(pctFOR, by=c("RegionalID", "pt_name", "pt_index")) %>%
  left_join(baa1, by=c("RegionalID", "pt_name", "pt_index")) %>%
  left_join(intermix1, by=c("RegionalID", "pt_name", "pt_index")) %>%
  left_join(interface1, by=c("RegionalID", "pt_name", "pt_index")) %>%
  left_join(wui1, by=c("RegionalID", "pt_name", "pt_index")) %>%
  left_join(lsm1, by=c("RegionalID", "pt_name", "pt_index")) %>%
  left_join(build1, by=c("RegionalID", "pt_name", "pt_index"))

## Scale and center variables
dat1[,c(8,17:106)] <- scale(dat1[,c(8,17:106)])
write_csv(dat1, "output/binary_model_data.csv")

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
  m1_pt <- glmer(bin.exp ~ Sex + Age +  *totalforest_30 + pasture_60 + laggedBMI_30 + (1|WMU) + (1|year), 
                 family=binomial(link="logit"), data=pt)
  
    # 5-fold cross validation
  row_idx <- sample(1:5, nrow(pt), replace=TRUE)
  for (j in 1:5) {
    
    # Split data into train & test
    trainSet <- pt[row_idx[row_idx==j],] %>% as_tibble()
    testSet <- pt[row_idx[row_idx!=j],] %>% as_tibble()
    
    # Fit model on training set
    m1_cv <- glmer(bin.exp ~ Sex*Age + totalag_60 + mix_15_100 + laggedBMI_30 + (1|WMU) + (1|year), 
                   family=binomial(link="logit"), data=pt)
    pred <- predict(m1_cv, newdata=testSet, type="response", se.fit=TRUE)
    
    # Create data frame to store prediction for each fold
    testTable <- testSet %>% select(RegionalID:WMUA_code,n.compounds.T)
    
    testTable$pred <- pred$fit
    testTable$se.fit <- pred$se.fit
    
    testTable <- testTable %>% 
      mutate(lower=pred-(se.fit*1.96), upper=pred+(se.fit*1.96)) %>%
      mutate(correct=if_else(n.compounds.T>=lower & n.compounds.T<=upper, 1, 0))
    
    perf[i,j+6] <- sum(testTable$correct)/nrow(testTable)
    perf[i,j+1] <- sqrt(mean((testTable$pred - testTable$n.compounds.T)^2))
  }
  
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

# Calculate average performance metric
perf <- as.data.frame(perf)
names(perf) <- c("Iter", "RMSE1", "RMSE2", "RMSE3", "RMSE4", "RMSE5", "Accur1", "Accur2", "Accur3", "Accur4", "Accur5")
perf <- perf %>% as_tibble() %>% 
  mutate(RMSE=(RMSE1+RMSE2+RMSE3+RMSE4+RMSE5)/5,
         Accuracy=(Accur1+Accur2+Accur3+Accur4+Accur5)/5) %>%
  select(Iter, RMSE, Accuracy)

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
  
  # 5-fold cross validation
  row_idx <- sample(1:5, nrow(pt), replace=TRUE)
  for (j in 1:5) {
    
    # Split data into train & test
    trainSet <- pt[row_idx[row_idx==j],] %>% as_tibble()
    testSet <- pt[row_idx[row_idx!=j],] %>% as_tibble()
    
    # Fit model on training set
    m1_cv <- glmer(bin.exp ~ Sex*Age + totalag_60 + mix_15_100 + 
                     laggedBMI_30 + (1|WMU) + (1|year), 
                   family=binomial(link="logit"), data=pt)
    pred <- predict(m1_cv, newdata=testSet, type="response", se.fit=TRUE)
    
    # Create data frame to store prediction for each fold
    testTable <- testSet %>% select(RegionalID:WMUA_code,n.compounds.T)
    
    testTable$pred <- pred$fit
    testTable$se.fit <- pred$se.fit
    
    testTable <- testTable %>% 
      mutate(lower=pred-(se.fit*1.96), upper=pred+(se.fit*1.96)) %>%
      mutate(correct=if_else(n.compounds.T>=lower & n.compounds.T<=upper, 1, 0))
    
    perf[i,j+6] <- sum(testTable$correct)/nrow(testTable)
    perf[i,j+1] <- sqrt(mean((testTable$pred - testTable$n.compounds.T)^2))
  }
  
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

# Calculate average performance metric
perf <- as.data.frame(perf)
names(perf) <- c("Iter", "RMSE1", "RMSE2", "RMSE3", "RMSE4", "RMSE5", "Accur1", "Accur2", "Accur3", "Accur4", "Accur5")
perf <- perf %>% as_tibble() %>% 
  mutate(RMSE=(RMSE1+RMSE2+RMSE3+RMSE4+RMSE5)/5,
         Accuracy=(Accur1+Accur2+Accur3+Accur4+Accur5)/5) %>%
  select(Iter, RMSE, Accuracy)

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
    
    # 5-fold cross validation
    row_idx <- sample(1:5, nrow(pt), replace=TRUE)
    for (j in 1:5) {
      
      # Split data into train & test
      trainSet <- pt[row_idx[row_idx==j],] %>% as_tibble()
      testSet <- pt[row_idx[row_idx!=j],] %>% as_tibble()
      
      # Fit model on training set
      m1_cv <- glmer(bin.exp ~ Sex*Age + totalag_60 + mix_15_100 + 
                       laggedBMI_30 + (1|WMU) + (1|year), 
                     family=binomial(link="logit"), data=pt)
      pred <- predict(m1_cv, newdata=testSet, type="response", se.fit=TRUE)
      
      # Create data frame to store prediction for each fold
      testTable <- testSet %>% select(RegionalID:WMUA_code,n.compounds.T)
      
      testTable$pred <- pred$fit
      testTable$se.fit <- pred$se.fit
      
      testTable <- testTable %>% 
        mutate(lower=pred-(se.fit*1.96), upper=pred+(se.fit*1.96)) %>%
        mutate(correct=if_else(n.compounds.T>=lower & n.compounds.T<=upper, 1, 0))
      
      perf[i,j+6] <- sum(testTable$correct)/nrow(testTable)
      perf[i,j+1] <- sqrt(mean((testTable$pred - testTable$n.compounds.T)^2))
    }
    
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

# Calculate average performance metric
perf <- as.data.frame(perf)
names(perf) <- c("Iter", "RMSE1", "RMSE2", "RMSE3", "RMSE4", "RMSE5", "Accur1", "Accur2", "Accur3", "Accur4", "Accur5")
perf <- perf %>% as_tibble() %>% 
  mutate(RMSE=(RMSE1+RMSE2+RMSE3+RMSE4+RMSE5)/5,
         Accuracy=(Accur1+Accur2+Accur3+Accur4+Accur5)/5) %>%
  select(Iter, RMSE, Accuracy)

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

