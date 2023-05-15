library(tidyverse)
library(lme4)
library(MuMIn)
library(caret)
library(broom)
library(performance)
library(partR2)

source("00_AR_functions.R")

wmua <- read_csv("data/wmu_wmua_conversion.csv")

# Load covariate information
covar <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>% 
  rename(BBA=baa) %>%
  select(RegionalID:rand_y,buffsize:lag_beechnuts)

# Load binary for each compound, combine covariates
brod <- read_csv("output/binary_brodifacoum.csv") %>% 
  left_join(covar, by="RegionalID") %>% 
  select(RegionalID, pt_name:rand_y, year:bin.exp.ntr,buffsize:lag_beechnuts) %>%
  covar_org()%>%
  select(-c(BMI_15, BMI_30, BMI_60)) %>%
  left_join(wmua, by="WMU")
brom <- read_csv("output/binary_bromadiolone.csv") %>% 
  left_join(covar, by="RegionalID") %>% 
  select(RegionalID, pt_name:rand_y, year:bin.exp.ntr,buffsize:lag_beechnuts) %>%
  covar_org()%>%
  select(-c(BMI_15, BMI_30, BMI_60)) %>%
  left_join(wmua, by="WMU")
diph <- read_csv("output/binary_diphacinone.csv") %>% 
  left_join(covar, by="RegionalID") %>% 
  select(RegionalID, pt_name:rand_y, year:bin.exp.ntr,buffsize:lag_beechnuts) %>%
  covar_org()%>%
  select(-c(BMI_15, BMI_30, BMI_60)) %>%
  left_join(wmua, by="WMU")
dico <- read_csv("output/binary_dicoumarol.csv") %>% 
  left_join(covar, by="RegionalID") %>% 
  select(RegionalID, pt_name:rand_y, year:bin.exp.ntr,buffsize:lag_beechnuts) %>%
  covar_org() %>%
  select(-c(BMI_15, BMI_30, BMI_60)) %>%
  left_join(wmua, by="WMU")

# Scale
brod[,c(8,17:78)] <- scale(brod[,c(8,17:78)])
brom[,c(8,17:78)] <- scale(brom[,c(8,17:78)])
diph[,c(8,17:78)] <- scale(diph[,c(8,17:78)])
dico[,c(8,17:78)] <- scale(dico[,c(8,17:78)])

#### Brodifacoum ####

## Run models ##

## Loop over each set of random points
m_est <- m_stderr <- m_zscore <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()

ranef_coef <- data.frame()

kappa <- matrix(NA, ncol=6, nrow=10)
kappa[,1] <- 1:10

cmpm <- data.frame()
vif <- list()
r2cm <- list()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- brod[brod$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- lme4::glmer(bin.exp ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * lag_beechnuts + (1|WMU), 
                 family=binomial(link="logit"), data=pt)
  
  m1s <- summary(m1_pt)
  m1sdf <- m1s$coefficients
  vif[[i]] <- check_collinearity(m1_pt)
  
  r2cm[[i]] <- r2(m1_pt, metrics="R2", ci=TRUE)
  
  # save averaged confidence intervals
  ci <- confint.merMod(m1_pt, method="Wald")
  ci <- ci[complete.cases(ci),]
  pct2.5 <- rbind(pct2.5, t(ci)[1,1:nrow(ci)])
  pct97.5 <- rbind(pct97.5, t(ci)[2,1:nrow(ci)])
  
  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  m_zscore <- rbind(m_zscore, coef(m1s)[,3])
  m_ranef <- rbind(m_ranef, unlist(VarCorr(m1_pt)))
  
  # Save the coefficients for each level of the random effect
  ranefc <- as_tibble(ranef(m1_pt)) 
  ranefc$iter <- i
  ranef_coef <- bind_rows(ranef_coef, ranefc)
  
  # 5-fold cross validation
  row_idx <- sample(1:5, nrow(pt), replace=TRUE)
  for (j in 1:5) {
    
    # Split data into train & test
    trainSet <- pt[row_idx==j,] %>% as_tibble()
    testSet <- pt[row_idx!=j,] %>% as_tibble()
    
    # Fit model on training set
    m1_cv <- glmer(bin.exp ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * lag_beechnuts + (1|WMU), 
                   family=binomial(link="logit"), data=pt)
    pred <- predict(m1_cv, newdata=testSet, type="response", re.form=NA)
    
    # Code from lme4 FAQ (Bolker)
    mm <- model.matrix(terms(m1_cv),testSet)
    pvar1 <- diag(mm %*% tcrossprod(vcov(m1_cv),mm))
    tvar1 <- pvar1+VarCorr(m1_cv)$WMU[1]
    cmult <- 1.96

    # Create data frame to store prediction for each fold
    testTable <- testSet %>% select(RegionalID:WMU,bin.exp)
    
    testTable <- data.frame(
      testTable, pred=round(pred)
      , plo = testSet$bin.exp-cmult*sqrt(pvar1)
      , phi = testSet$bin.exp+cmult*sqrt(pvar1)
      , tlo = testSet$bin.exp-cmult*sqrt(tvar1)
      , thi = testSet$bin.exp+cmult*sqrt(tvar1)
    )

    bin_confusion <- confusionMatrix(
      # predictions then true values
      data = factor(testTable$pred, levels=0:1),
      reference = factor(testTable$bin.exp, levels=0:1),
      # what is positive exposure value
      positive = "1"
    )
    
    cm <- as.data.frame(t(bin_confusion$overall))
    cm$iteration <- i
    cm$fold <- j
    cmpm <- bind_rows(cmpm, cm)
  }
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- row.names(m1sdf)
    names(m_stderr) <- row.names(m1sdf)
    names(m_zscore) <- row.names(m1sdf)
    names(pct2.5) <- row.names(m1sdf)
    names(pct97.5) <- row.names(m1sdf)
    names(m_ranef) <- c("RE_WMU")
  }
  
}

# Calculate average performance metric
mean(cmpm$Kappa)
mean(cmpm$Accuracy)

# save overall performance stats
write_csv(cmpm, "results/binary_brodifacoum_performance.csv")

# save VIF values
saveRDS(vif, "results/brod_VIF.rds")

r2cm_out <- plyr::ldply(r2cm)
write_csv(r2cm_out, "results/r2_glmm_brod.csv")

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
zscore_avg <- colMeans(m_zscore[sapply(m_zscore, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)
ranef_avg <- as.data.frame(colMeans(m_ranef[sapply(m_ranef, is.numeric)])) %>% rownames_to_column("RE")
names(ranef_avg)[2] <- "variance"

# Save averaged values for each level of random effect
re <- ranef_coef %>% 
            group_by(grp) %>% 
            summarize(REval=mean(condval), REsd=mean(condsd), 
              RErangeL=range(condval)[1], RErangeH=range(condval)[2]) 
write_csv(re, "results/brodifacoum_random_effects_coefficients.csv")


# Calculate imputation variance and 95% CI
m_var <- m_stderr**2
u_bar <- colMeans(m_var[sapply(m_var, is.numeric)], na.rm=TRUE) # within impuation variance
btwnimpvar <- sapply(m_est, B)

total_var <- data.frame()
for (i in 1:length(u_bar)) {
  
  v <- u_bar[i] + (1 + (1/10))*btwnimpvar[i]
  total_var <- bind_rows(total_var, v)
  
}

# consense, saving only the diagonals (removing NAs)
coef_var <- total_var %>% summarise(across(everything(), ~ na.omit(.x)))

# Calculate confidence intervals
ucl <- coef_avg + 1.96*(sqrt(coef_var)/sqrt(338))
lcl <- coef_avg - 1.96*(sqrt(coef_var)/sqrt(338))

# Combine and clean up data frame for fixed effects
coef_summary <- bind_rows(coef_avg, stderr_avg, zscore_avg, lcl, ucl, pct2.5_avg, pct97.5_avg)
names(coef_summary) <- row.names(m1sdf)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "z-score", "2.5CI", "97.5CI", "2.5CIavg", "97.5CIavg")
coef_summary <- data.frame(coef=coefs, coef_summary)
names(coef_summary)[2] <- "Intercept"

# Write to file
write_csv(coef_summary, "results/binaryTbrodifacoum_coef-summary.csv")
write_csv(ranef_avg, "results/binaryTbrodifacoum_random-effect-var.csv")

#### Bromadiolone ####
## Running final models ##

## Loop over each set of random points
m_est <- m_stderr <- m_zscore <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()

kappa <- matrix(NA, ncol=6, nrow=10)
kappa[,1] <- 1:10

ranef_coef <- data.frame()

cmpm <- data.frame()

vif <- list()
r2cm <- list()
# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- brom[brom$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- glmer(bin.exp ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * lag_beechnuts + (1|WMU),
                 family=binomial(link="logit"),data=pt)
  
  vif[[i]] <- check_collinearity(m1_pt)
  
  r2cm[[i]] <- r2(m1_pt, metrics="R2", ci=TRUE)
  
  m1s <- summary(m1_pt)
  m1sdf <- m1s$coefficients
  
  # Save the coefficients for each level of the random effect
  ranefc <- as_tibble(ranef(m1_pt)) 
  ranefc$iter <- i
  ranef_coef <- bind_rows(ranef_coef, ranefc)
  
  # save averaged confidence intervals
  ci <- confint.merMod(m1_pt, method="Wald")
  ci <- ci[complete.cases(ci),]
  pct2.5 <- rbind(pct2.5, t(ci)[1,1:nrow(ci)])
  pct97.5 <- rbind(pct97.5, t(ci)[2,1:nrow(ci)])
  
  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  m_zscore <- rbind(m_zscore, coef(m1s)[,3])
  m_ranef <- rbind(m_ranef, unlist(VarCorr(m1_pt)))
  
  # 5-fold cross validation
  row_idx <- sample(1:5, nrow(pt), replace=TRUE)
  for (j in 1:5) {
    
    # Split data into train & test
    trainSet <- pt[row_idx==j,] %>% as_tibble()
    testSet <- pt[row_idx!=j,] %>% as_tibble()
    
    # Fit model on training set
    m1_cv <- glmer(bin.exp ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * lag_beechnuts + (1|WMU), 
                   family=binomial(link="logit"), data=pt)
    pred <- predict(m1_cv, newdata=testSet, type="response", re.form=NA)
    
    # Code from lme4 FAQ (Bolker)
    mm <- model.matrix(terms(m1_cv),testSet)
    pvar1 <- diag(mm %*% tcrossprod(vcov(m1_cv),mm))
    tvar1 <- pvar1+VarCorr(m1_cv)$WMU[1]
    cmult <- 1.96
    
    # Create data frame to store prediction for each fold
    testTable <- testSet %>% select(RegionalID:WMU,bin.exp)
    
    testTable <- data.frame(
      testTable, pred=round(pred)
      , plo = testSet$bin.exp-cmult*sqrt(pvar1)
      , phi = testSet$bin.exp+cmult*sqrt(pvar1)
      , tlo = testSet$bin.exp-cmult*sqrt(tvar1)
      , thi = testSet$bin.exp+cmult*sqrt(tvar1)
    )
    
    bin_confusion <- confusionMatrix(
      # predictions then true values
      data = factor(testTable$pred, levels=0:1),
      reference = factor(testTable$bin.exp, levels=0:1),
      # what is positive exposure value
      positive = "1"
    )
    
    cm <- as.data.frame(t(bin_confusion$overall))
    cm$iteration <- i
    cm$fold <- j
    cmpm <- bind_rows(cmpm, cm)
  }

  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- row.names(m1sdf)
    names(m_stderr) <- row.names(m1sdf)
    names(m_zscore) <- row.names(m1sdf)
    names(pct2.5) <- row.names(m1sdf)
    names(pct97.5) <- row.names(m1sdf)
    names(m_ranef) <- c("RE_WMU")
  }
}

# Calculate average performance metric
mean(cmpm$Kappa)
mean(cmpm$Accuracy)

# save overall performance stats
write_csv(cmpm, "results/binary_bromadiolone_performance.csv")

# save VIF values
saveRDS(vif, "results/brom_VIF.rds")

r2cm_out <- plyr::ldply(r2cm)
write_csv(r2cm_out, "results/r2_glmm_brom.csv")

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
zscore_avg <- colMeans(m_zscore[sapply(m_zscore, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)
ranef_avg <- as.data.frame(colMeans(m_ranef[sapply(m_ranef, is.numeric)])) %>% rownames_to_column("RE")
names(ranef_avg)[2] <- "variance"

# Save averaged values for each level of random effect
re <- ranef_coef %>% 
  group_by(grp) %>% 
  summarize(REval=mean(condval), REsd=mean(condsd), 
            RErangeL=range(condval)[1], RErangeH=range(condval)[2]) 
write_csv(re, "results/bromadiolone_random_effects_coefficients.csv")

# Calculate imputation variance and 95% CI
m_var <- m_stderr**2
u_bar <- colMeans(m_var[sapply(m_var, is.numeric)], na.rm=TRUE) # within impuation variance
btwnimpvar <- sapply(m_est, B)

total_var <- data.frame()
for (i in 1:length(u_bar)) {
  
  v <- u_bar[i] + (1 + (1/10))*btwnimpvar[i]
  total_var <- bind_rows(total_var, v)
  
}

# consense, saving only the diagonals (removing NAs)
coef_var <- total_var %>% summarise(across(everything(), ~ na.omit(.x)))

# Calculate confidence intervals
ucl <- coef_avg + 1.96*(sqrt(coef_var)/sqrt(338))
lcl <- coef_avg - 1.96*(sqrt(coef_var)/sqrt(338))

# Combine and clean up data frame for fixed effects
coef_summary <- bind_rows(coef_avg, stderr_avg, zscore_avg, lcl, ucl, pct2.5_avg, pct97.5_avg)
names(coef_summary) <- row.names(m1sdf)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "z-score", "2.5CI", "97.5CI", "2.5CIavg", "97.5CIavg")
coef_summary <- data.frame(coef=coefs, coef_summary)
names(coef_summary)[2] <- "Intercept"

# Write to file
write_csv(coef_summary, "results/binaryTbromadiolone_coef-summary.csv")
write_csv(ranef_avg, "results/binaryTbromadiolone_random-effect-var.csv")

#### Diphacinone ####

## Loop over each set of random points
m_est <- m_stderr <- m_zscore <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()

kappa <- matrix(NA, ncol=6, nrow=10)
kappa[,1] <- 1:10

ranef_coef <- data.frame()
cmpm <- data.frame()
vif <- list()
r2cm <- list()
# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- diph[diph$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- glmer(bin.exp ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * lag_beechnuts + (1|WMU), 
               family=binomial(link="logit"), data=pt)
  
  vif[[i]] <- check_collinearity(m1_pt)
  
  r2cm[[i]] <- r2(m1_pt, metrics="R2", ci=TRUE)
  
    
    m1s <- summary(m1_pt)
    m1sdf <- m1s$coefficients
    
    # Save the coefficients for each level of the random effect
    ranefc <- as_tibble(ranef(m1_pt)) 
    ranefc$iter <- i
    ranef_coef <- bind_rows(ranef_coef, ranefc)
    
    # save averaged confidence intervals
    ci <- confint.merMod(m1_pt, method="Wald")
    ci <- ci[complete.cases(ci),]
    pct2.5 <- rbind(pct2.5, t(ci)[1,1:nrow(ci)])
    pct97.5 <- rbind(pct97.5, t(ci)[2,1:nrow(ci)])
    
    # Save point set estimates
    m_est <- rbind(m_est, coef(m1s)[,1])
    m_stderr <- rbind(m_stderr, coef(m1s)[,2])
    m_zscore <- rbind(m_zscore, coef(m1s)[,3])
    m_ranef <- rbind(m_ranef, unlist(VarCorr(m1_pt)))
    
    # 5-fold cross validation
    row_idx <- sample(1:5, nrow(pt), replace=TRUE)
    for (j in 1:5) {
      
      # Split data into train & test
      trainSet <- pt[row_idx==j,] %>% as_tibble()
      testSet <- pt[row_idx!=j,] %>% as_tibble()
      
      # Fit model on training set
      m1_cv <- glmer(bin.exp ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * lag_beechnuts + (1|WMU), 
                     family=binomial(link="logit"), data=pt)
      pred <- predict(m1_cv, newdata=testSet, type="response", re.form=NA)
      
      # Code from lme4 FAQ (Bolker)
      mm <- model.matrix(terms(m1_cv),testSet)
      pvar1 <- diag(mm %*% tcrossprod(vcov(m1_cv),mm))
      tvar1 <- pvar1+VarCorr(m1_cv)$WMU[1]
      cmult <- 1.96
      
      # Create data frame to store prediction for each fold
      testTable <- testSet %>% select(RegionalID:WMU,bin.exp)
      
      testTable <- data.frame(
        testTable, pred=round(pred)
        , plo = testSet$bin.exp-cmult*sqrt(pvar1)
        , phi = testSet$bin.exp+cmult*sqrt(pvar1)
        , tlo = testSet$bin.exp-cmult*sqrt(tvar1)
        , thi = testSet$bin.exp+cmult*sqrt(tvar1)
      )
      
      bin_confusion <- confusionMatrix(
        # predictions then true values
        data = factor(testTable$pred, levels=0:1),
        reference = factor(testTable$bin.exp, levels=0:1),
        # what is positive exposure value
        positive = "1"
      )
      
      cm <- as.data.frame(t(bin_confusion$overall))
      cm$iteration <- i
      cm$fold <- j
      cmpm <- bind_rows(cmpm, cm)
    }

    # Rename (only need to do once)
    if (i==1) {
      names(m_est) <- row.names(m1sdf)
      names(m_stderr) <- row.names(m1sdf)
      names(m_zscore) <- row.names(m1sdf)
      names(pct2.5) <- row.names(m1sdf)
      names(pct97.5) <- row.names(m1sdf)
      names(m_ranef) <- c("RE_WMU")
    }
    
  # }
}

# Calculate average performance metric
mean(cmpm$Kappa)
mean(cmpm$Accuracy)

# save overall performance stats
write_csv(cmpm, "results/binary_diphacinone_performance.csv")

# save VIF values
saveRDS(vif, "results/diph_VIF.rds")

r2cm_out <- plyr::ldply(r2cm)
write_csv(r2cm_out, "results/r2_glmm_diph.csv")

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
zscore_avg <- colMeans(m_zscore[sapply(m_zscore, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)
ranef_avg <- as.data.frame(colMeans(m_ranef[sapply(m_ranef, is.numeric)])) %>% rownames_to_column("RE")
names(ranef_avg)[2] <- "variance"

# Save averaged values for each level of random effect
re <- ranef_coef %>% 
  group_by(grp) %>% 
  summarize(REval=mean(condval), REsd=mean(condsd), 
            RErangeL=range(condval)[1], RErangeH=range(condval)[2]) 
write_csv(re, "results/diphacinone_random_effects_coefficients.csv")

# Calculate imputation variance and 95% CI
m_var <- m_stderr**2
u_bar <- colMeans(m_var[sapply(m_var, is.numeric)], na.rm=TRUE) # within impuation variance
btwnimpvar <- sapply(m_est, B)

total_var <- data.frame()
for (i in 1:length(u_bar)) {
  
  v <- u_bar[i] + (1 + (1/10))*btwnimpvar[i]
  total_var <- bind_rows(total_var, v)
  
}

# consense, saving only the diagonals (removing NAs)
coef_var <- total_var %>% summarise(across(everything(), ~ na.omit(.x)))

# Calculate confidence intervals
ucl <- coef_avg + 1.96*(sqrt(coef_var)/sqrt(338))
lcl <- coef_avg - 1.96*(sqrt(coef_var)/sqrt(338))

# Combine and clean up data frame for fixed effects
coef_summary <- bind_rows(coef_avg, stderr_avg, zscore_avg, lcl, ucl, pct2.5_avg, pct97.5_avg)
names(coef_summary) <- row.names(m1sdf)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "z-score", "2.5CI", "97.5CI", "2.5CIavg", "97.5CIavg")
coef_summary <- data.frame(coef=coefs, coef_summary)
names(coef_summary)[2] <- "Intercept"

# Write to file
write_csv(coef_summary, "results/binaryTdiphacinone_coef-summary.csv")
write_csv(ranef_avg, "results/binaryTdiphacinone_random-effect-var.csv")


#### Dicoumarol ####

## Loop over each set of random points
m_est <- m_stderr <- m_zscore <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()

kappa <- matrix(NA, ncol=6, nrow=10)
kappa[,1] <- 1:10

ranef_coef <- data.frame()
cmpm <- data.frame()
vif <- list()
r2cm <- list()
# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- dico[dico$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- glmer(bin.exp ~ Sex + Age + pasture_15 + (1|Region), 
                 family=binomial(link="logit"), data=pt)
  
  m1s <- summary(m1_pt)
  m1sdf <- m1s$coefficients
  vif[[i]] <- check_collinearity(m1_pt)
  
  r2cm[[i]] <- r2(m1_pt, metrics="R2", ci=TRUE)
  
  # Save the coefficients for each level of the random effect
  ranefc <- as_tibble(ranef(m1_pt)) 
  ranefc$iter <- i
  ranef_coef <- bind_rows(ranef_coef, ranefc)
  
  # save averaged confidence intervals
  ci <- confint.merMod(m1_pt, method="Wald")
  ci <- ci[complete.cases(ci),]
  pct2.5 <- rbind(pct2.5, t(ci)[1,1:nrow(ci)])
  pct97.5 <- rbind(pct97.5, t(ci)[2,1:nrow(ci)])
  
  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  m_zscore <- rbind(m_zscore, coef(m1s)[,3])
  m_ranef <- rbind(m_ranef, unlist(VarCorr(m1_pt)))
  
  # 5-fold cross validation
  row_idx <- sample(1:5, nrow(pt), replace=TRUE)
  for (j in 1:5) {
    # Split data into train & test
    trainSet <- pt[row_idx==j,] %>% as_tibble()
    testSet <- pt[row_idx!=j,] %>% as_tibble()
    # Fit model on training set
    m1_cv <- glmer(bin.exp ~ Sex + Age + pasture_15 + (1|Region),
                   family=binomial(link="logit"), data=pt)
    pred <- predict(m1_cv, newdata=testSet, type="response", re.form=NA)
    # Code from lme4 FAQ (Bolker)
    mm <- model.matrix(terms(m1_cv),testSet)
    pvar1 <- diag(mm %*% tcrossprod(vcov(m1_cv),mm))
    tvar1 <- pvar1+VarCorr(m1_cv)$Region[1]
    cmult <- 1.96
    # Create data frame to store prediction for each fold
    testTable <- testSet %>% select(RegionalID:WMUA,bin.exp)
    testTable <- data.frame(
      testTable, pred=round(pred)
      , plo = testSet$bin.exp-cmult*sqrt(pvar1)
      , phi = testSet$bin.exp+cmult*sqrt(pvar1)
      , tlo = testSet$bin.exp-cmult*sqrt(tvar1)
      , thi = testSet$bin.exp+cmult*sqrt(tvar1)
    )
    bin_confusion <- confusionMatrix(
      # predictions then true values
      data = factor(testTable$pred, levels=0:1),
      reference = factor(testTable$bin.exp, levels=0:1),
      # what is positive exposure value
      positive = "1"
    )
    cm <- as.data.frame(t(bin_confusion$overall))
    cm$iteration <- i
    cm$fold <- j
    cmpm <- bind_rows(cmpm, cm)
  }
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- row.names(m1sdf)
    names(m_stderr) <- row.names(m1sdf)
    names(m_zscore) <- row.names(m1sdf)
    names(pct2.5) <- row.names(m1sdf)
    names(pct97.5) <- row.names(m1sdf)
    names(m_ranef) <- c("RE_WMU")
  }

  # }
}

# Calculate average performance metric
mean(cmpm$Kappa)
mean(cmpm$Accuracy)

# save overall performance stats
write_csv(cmpm, "results/binary_dico_performance.csv")

# save VIF values
saveRDS(vif, "results/dico_VIF.rds")

r2cm_out <- plyr::ldply(r2cm)
write_csv(r2cm_out, "results/r2_glmm_dico.csv")

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
zscore_avg <- colMeans(m_zscore[sapply(m_zscore, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)
ranef_avg <- as.data.frame(colMeans(m_ranef[sapply(m_ranef, is.numeric)])) %>% rownames_to_column("RE")
names(ranef_avg)[2] <- "variance"

# Save averaged values for each level of random effect
re <- ranef_coef %>% 
  group_by(grp) %>% 
  summarize(REval=mean(condval), REsd=mean(condsd), 
            RErangeL=range(condval)[1], RErangeH=range(condval)[2]) 
write_csv(re, "results/dico_random_effects_coefficients.csv")

# Calculate imputation variance and 95% CI
m_var <- m_stderr**2
u_bar <- colMeans(m_var[sapply(m_var, is.numeric)], na.rm=TRUE) # within impuation variance
btwnimpvar <- sapply(m_est, B)

total_var <- data.frame()
for (i in 1:length(u_bar)) {
  
  v <- u_bar[i] + (1 + (1/10))*btwnimpvar[i]
  total_var <- bind_rows(total_var, v)
  
}

# consense, saving only the diagonals (removing NAs)
coef_var <- total_var %>% summarise(across(everything(), ~ na.omit(.x)))

# Calculate confidence intervals
ucl <- coef_avg + 1.96*(sqrt(coef_var)/sqrt(338))
lcl <- coef_avg - 1.96*(sqrt(coef_var)/sqrt(338))

# Combine and clean up data frame for fixed effects
coef_summary <- bind_rows(coef_avg, stderr_avg, zscore_avg, lcl, ucl, pct2.5_avg, pct97.5_avg)
names(coef_summary) <- row.names(m1sdf)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "z-score", "2.5CI", "97.5CI", "2.5CIavg", "97.5CIavg")
coef_summary <- data.frame(coef=coefs, coef_summary)
names(coef_summary)[2] <- "Intercept"

# Write to file
write_csv(coef_summary, "results/binaryTdicoumarol_coef-summary.csv")
write_csv(ranef_avg, "results/binaryTdicoumarol_random-effect-var.csv")