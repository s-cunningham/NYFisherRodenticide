## Functions for rodenticide

library(tidyverse)

# Predict exposure by age
logistic_pred_age <- function(fixed, random, compound, sex, meanWUI, meanPasture, meanBBA, meanMast, ageStart, ageEnd, lo) {
  
  # create a sequence of values to estimate for age
  age_iter <- seq(ageStart, ageEnd, length.out=lo)
  
  # create empty data frame to store
  full_vals <- data.frame()
  
  # Loop over each level of random effect
  for (j in 1:nrow(random)) {
    
    # Empty vector for each level of random effect
    exp_val <- c()
    
    # Loop over values of age
    for (i in 1:length(age_iter)) {
      
      exp_val[i] <- inv.logit(fixed$intercept[1] + 
                                fixed$SexM[1]*sex + 
                                fixed$Age[1]*age_iter[i] + 
                                fixed$Age2[1]*(age_iter[i]^2) +
                                fixed$WUI[1]*meanWUI + 
                                fixed$pasture[1]*meanPasture + 
                                fixed$basalarea[1]*meanBBA +
                                fixed$mast[1]*meanMast + 
                                fixed$intx_beech[1]*meanBBA*meanMast +
                                random$REval[j])
      
      # Does RE need to be weighted somehow?
      # level_probs <- random %>% group_by(grp) %>% summarize(wmu_prop=n()/nrow(random))
      
    }
    
    # Save the expected values to data frame
    full_vals <- bind_rows(full_vals, as.data.frame(t(exp_val)))
    
  }
  
  names(full_vals) <- 1:lo
  full_vals$level <- random$grp
  
  full_vals <- pivot_longer(full_vals, 1:100, names_to="index", values_to="exp_val") %>%
    mutate_at(c('index', 'exp_val'), as.numeric)
  
  full_vals$Age <- rep(age_iter, 55)
  
  # add column for sex
  full_vals$compound <- compound
  full_vals$sex <- ifelse(sex==1, "Male", "Female")
  
  full_vals <- full_vals %>% select(compound, sex, level, Age, index, exp_val)
  
  return(full_vals)
  
}

poisson_pred_age <- function(fixed, random, sex, meanWUI, meanPasture, meanBBA, meanMast, ageStart, ageEnd, lo) {
  
  # crete a sequence of values to estimate for age
  age_iter <- seq(ageStart, ageEnd, length.out=lo)
  
  # create empty data frame to store
  full_vals <- data.frame()
  
  # Loop over each level of random effect
  for (j in 1:nrow(random)) {
    
    # Empty vector for each level of random effect
    exp_val <- c()
    
    # Loop over values of age
    for (i in 1:length(age_iter)) {
      
      exp_val[i] <- exp(fixed$intercept[1] + 
                          fixed$SexM[1]*sex + 
                          fixed$Age[1]*age_iter[i] + 
                          fixed$Age2[1]*(age_iter[i]^2) +
                          fixed$WUI[1]*meanWUI + 
                          fixed$pasture[1]*meanPasture + 
                          fixed$basalarea[1]*meanBBA +
                          fixed$mast[1]*meanMast + 
                          fixed$intx_beech[1]*meanBBA*meanMast +
                          random$REval[j])
    }
    
    # Save the expected values to data frame
    full_vals <- bind_rows(full_vals, as.data.frame(t(exp_val)))
    
  }
  
  # Does RE need to be weighted somehow? maybe in calculating the mean?
  # level_probs <- random %>% 
  #                   group_by(grp) %>% 
  #                   summarize(wmu_prop=n()/nrow(random)) %>%
  #                   rename(level=grp)
  
  names(full_vals) <- 1:lo
  full_vals$level <- random$grp
  
  full_vals <- pivot_longer(full_vals, 1:100, names_to="index", values_to="exp_val") %>%
    mutate_at(c('index', 'exp_val'), as.numeric)
  
  full_vals$Age <- rep(age_iter, 55)
  
  # add column for sex
  full_vals$sex <- ifelse(sex==1, "Male", "Female")
  
  full_vals <- full_vals %>% select(sex, level, Age, index, exp_val)
  # full_vals <- left_join(full_vals, level_probs, by="level")
  
  return(full_vals)
  
}

# Predict beech mast count
logistic_pred_mast <- function(fixed, random, compound, sex, meanWUI, meanPasture, meanBBA, meanAge, mastStart, mastEnd, lo) {
  
  # crete a sequence of values to estimate for age
  mast_iter <- seq(mastStart, mastEnd, length.out=lo)
  
  # create empty data frame to store
  full_vals <- data.frame()
  
  # Loop over each level of random effect
  for (j in 1:nrow(random)) {
    
    # Empty vector for each level of random effect
    exp_val <- c()
    
    # Loop over values of age
    for (i in 1:length(mast_iter)) {
      
      exp_val[i] <- inv.logit(fixed$intercept[1] + 
                                fixed$SexM[1]*sex + 
                                fixed$Age[1]*meanAge + 
                                fixed$Age2[1]*(meanAge^2) +
                                fixed$WUI[1]*meanWUI + 
                                fixed$pasture[1]*meanPasture + 
                                fixed$basalarea[1]*meanBBA +
                                fixed$mast[1]*mast_iter[i] + 
                                fixed$intx_beech[1]*meanBBA*mast_iter[i] +
                                random$REval[j])
      
      # Does RE need to be weighted somehow?
      # level_probs <- random %>% group_by(grp) %>% summarize(wmu_prop=n()/nrow(random))
      
    }
    
    # Save the expected values to data frame
    full_vals <- bind_rows(full_vals, as.data.frame(t(exp_val)))
    
  }
  
  names(full_vals) <- 1:lo
  full_vals$level <- random$grp
  
  full_vals <- pivot_longer(full_vals, 1:100, names_to="index", values_to="exp_val") %>%
    mutate_at(c('index', 'exp_val'), as.numeric)
  
  full_vals$Mast <- rep(mast_iter, 55)
  
  # add column for sex
  full_vals$compound <- compound
  full_vals$sex <- ifelse(sex==1, "Male", "Female")
  
  full_vals <- full_vals %>% select(compound, sex, level, Mast, index, exp_val)
  
  return(full_vals)
  
}

poisson_pred_mast <- function(fixed, random, sex, meanWUI, meanPasture, meanBBA, meanAge, mastStart, mastEnd, lo) {
  
  # create a sequence of values to estimate for age
  mast_iter <- seq(mastStart, mastEnd, length.out=lo)
  
  # create empty data frame to store
  full_vals <- data.frame()
  
  # Loop over each level of random effect
  for (j in 1:nrow(random)) {
    
    # Empty vector for each level of random effect
    exp_val <- c()
    
    # Loop over values of age
    for (i in 1:length(mast_iter)) {
      
      exp_val[i] <- exp(fixed$intercept[1] + 
                          fixed$SexM[1]*sex + 
                          fixed$Age[1]*meanAge + 
                          fixed$Age2[1]*(meanAge^2) +
                          fixed$WUI[1]*meanWUI +
                          fixed$pasture[1]*meanPasture + 
                          fixed$basalarea[1]*meanBBA +
                          fixed$mast[1]*mast_iter[i] + 
                          fixed$intx_beech[1]*meanBBA*mast_iter[i] +
                          random$REval[j])
      
    }
    
    # Save the expected values to data frame
    full_vals <- bind_rows(full_vals, as.data.frame(t(exp_val)))
    
  }
  
  names(full_vals) <- 1:lo
  full_vals$level <- random$grp
  
  full_vals <- pivot_longer(full_vals, 1:100, names_to="index", values_to="exp_val") %>%
    mutate_at(c('index', 'exp_val'), as.numeric)
  
  full_vals$Mast <- rep(mast_iter, 55)
  
  # add column for sex
  full_vals$sex <- ifelse(sex==1, "Male", "Female")
  
  full_vals <- full_vals %>% select(sex, level, Mast, index, exp_val)
  # full_vals <- left_join(full_vals, level_probs, by="level")
  
  return(full_vals)
  
}

# Predict proportion interface
logistic_pred_wui <- function(fixed, random, compound, sex, meanAge, meanPasture, meanBBA, meanMast, wuiStart, wuiEnd, lo) {
  
  # crete a sequence of values to estimate for age
  wui_iter <- seq(wuiStart, wuiEnd, length.out=lo)
  
  # create empty data frame to store
  full_vals <- data.frame()
  
  # Loop over each level of random effect
  for (j in 1:nrow(random)) {
    
    # Empty vector for each level of random effect
    exp_val <- c()
    
    # Loop over values of age
    for (i in 1:length(wui_iter)) {
      
      exp_val[i] <- inv.logit(fixed$intercept[1] + 
                                fixed$SexM[1]*sex + 
                                fixed$Age[1]*meanAge + 
                                fixed$Age2[1]*(meanAge^2) +
                                fixed$WUI[1]*wui_iter[i] + 
                                fixed$pasture[1]*meanPasture + 
                                fixed$basalarea[1]*meanBBA +
                                fixed$mast[1]*meanMast + 
                                fixed$intx_beech[1]*meanBBA*meanMast +
                                random$REval[j])
      
      # Does RE need to be weighted somehow?
      # level_probs <- random %>% group_by(grp) %>% summarize(wmu_prop=n()/nrow(random))
      
    }
    
    # Save the expected values to data frame
    full_vals <- bind_rows(full_vals, as.data.frame(t(exp_val)))
    
  }
  
  names(full_vals) <- 1:lo
  full_vals$level <- random$grp
  
  full_vals <- pivot_longer(full_vals, 1:100, names_to="index", values_to="exp_val") %>%
    mutate_at(c('index', 'exp_val'), as.numeric)
  
  full_vals$WUI <- rep(wui_iter, 55)
  
  # add column for sex
  full_vals$compound <- compound
  full_vals$sex <- ifelse(sex==1, "Male", "Female")
  
  full_vals <- full_vals %>% select(compound, sex, level, WUI, index, exp_val)
  
  return(full_vals)
  
}

poisson_pred_wui <- function(fixed, random, sex, meanAge, meanPasture, meanBBA, meanMast, wuiStart, wuiEnd, lo) {
  
  # create a sequence of values to estimate for age
  wui_iter <- seq(wuiStart, wuiEnd, length.out=lo)
  
  # create empty data frame to store
  full_vals <- data.frame()
  
  # Loop over each level of random effect
  for (j in 1:nrow(random)) {
    
    # Empty vector for each level of random effect
    exp_val <- c()
    
    # Loop over values of age
    for (i in 1:length(wui_iter)) {
      
      exp_val[i] <- exp(fixed$intercept[1] + 
                          fixed$SexM[1]*sex + 
                          fixed$Age[1]*meanAge + 
                          fixed$Age2[1]*(meanAge^2) +
                          fixed$WUI[1]*wui_iter[i] +
                          fixed$pasture[1]*meanPasture + 
                          fixed$basalarea[1]*meanBBA +
                          fixed$mast[1]*meanMast + 
                          fixed$intx_beech[1]*meanBBA*meanMast +
                          random$REval[j])
      
    }
    
    # Save the expected values to data frame
    full_vals <- bind_rows(full_vals, as.data.frame(t(exp_val)))
    
  }
  
  names(full_vals) <- 1:lo
  full_vals$level <- random$grp
  
  full_vals <- pivot_longer(full_vals, 1:100, names_to="index", values_to="exp_val") %>%
    mutate_at(c('index', 'exp_val'), as.numeric)
  
  full_vals$WUI <- rep(wui_iter, 55)
  
  # add column for sex
  full_vals$sex <- ifelse(sex==1, "Male", "Female")
  
  full_vals <- full_vals %>% select(sex, level, WUI, index, exp_val)
  # full_vals <- left_join(full_vals, level_probs, by="level")
  
  return(full_vals)
  
}

# Predict proportion pasture
logistic_pred_wui <- function(fixed, random, compound, sex, meanAge, meanWUI, meanBBA, meanMast, agStart, agEnd, lo) {
  
  # crete a sequence of values to estimate for age
  ag_iter <- seq(agStart, agEnd, length.out=lo)
  
  # create empty data frame to store
  full_vals <- data.frame()
  
  # Loop over each level of random effect
  for (j in 1:nrow(random)) {
    
    # Empty vector for each level of random effect
    exp_val <- c()
    
    # Loop over values of age
    for (i in 1:length(ag_iter)) {
      
      exp_val[i] <- inv.logit(fixed$intercept[1] + 
                                fixed$SexM[1]*sex + 
                                fixed$Age[1]*meanAge + 
                                fixed$Age2[1]*(meanAge^2) +
                                fixed$WUI[1]*meanWUI + 
                                fixed$pasture[1]*ag_iter[i] + 
                                fixed$basalarea[1]*meanBBA +
                                fixed$mast[1]*meanMast + 
                                fixed$intx_beech[1]*meanBBA*meanMast +
                                random$REval[j])
    
    }
    
    # Save the expected values to data frame
    full_vals <- bind_rows(full_vals, as.data.frame(t(exp_val)))
    
  }
  
  names(full_vals) <- 1:lo
  full_vals$level <- random$grp
  
  full_vals <- pivot_longer(full_vals, 1:100, names_to="index", values_to="exp_val") %>%
    mutate_at(c('index', 'exp_val'), as.numeric)
  
  full_vals$pasture <- rep(ag_iter, 55)
  
  # add column for sex
  full_vals$compound <- compound
  full_vals$sex <- ifelse(sex==1, "Male", "Female")
  
  full_vals <- full_vals %>% select(compound, sex, level, pasture, index, exp_val)
  
  return(full_vals)
  
}


## reorganize data for anlaysis
covar_org <- function(dat) {
  # Make random effects factors
  dat$WMU <- as.factor(dat$WMU)
  dat$year <- factor(dat$year)
  dat$RegionalID <- factor(dat$RegionalID)
  
  ## Percent agriculture
  pctAG <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, pasture, crops, totalag) %>% 
    distinct() %>% 
    group_by(RegionalID) %>% 
    pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% 
    as.data.frame()
  
  ## Beech basal area
  baa1 <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, BMI, laggedBMI, BBA, beechnuts, lag_beechnuts) %>% 
    distinct() %>% 
    group_by(RegionalID) %>% 
    pivot_wider(names_from=buffsize, values_from=c(BMI, laggedBMI, BBA, beechnuts, lag_beechnuts), values_fn=unique) %>% 
    as.data.frame()
  
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
    group_by(RegionalID) %>% 
    pivot_wider(names_from=buffrad, values_from=intermix) %>% 
    as.data.frame()
  names(intermix1)[4:12] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                              "mix_15_250", "mix_30_250", "mix_60_250",
                              "mix_15_500", "mix_30_500", "mix_60_500") 
  
  interface1  <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, radius, interface) %>%
    unite("buffrad", 4:5, sep="_") %>% 
    group_by(RegionalID) %>% 
    pivot_wider(names_from=buffrad, values_from=interface) %>% 
    as.data.frame()
  names(interface1)[4:12] <- c("face_15_100", "face_30_100", "face_60_100",
                               "face_15_250", "face_30_250", "face_60_250",
                               "face_15_500", "face_30_500", "face_60_500") 
  
  wui1  <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, radius, totalWUI) %>%
    unite("buffrad", 4:5, sep="_") %>% 
    group_by(RegionalID) %>% 
    pivot_wider(names_from=buffrad, values_from=totalWUI) %>% 
    as.data.frame()
  names(wui1)[4:12] <- c("wui_15_100", "wui_30_100", "wui_60_100",
                         "wui_15_250", "wui_30_250", "wui_60_250",
                         "wui_15_500", "wui_30_500", "wui_60_500") 
  
  #### Set up data to run for each combination of covariates ####
  dat1 <- dat %>% select(RegionalID:Town,bin.exp,bin.exp.ntr,beechnuts,lag_beechnuts) %>% distinct() %>%
    left_join(pctAG, by=c("RegionalID", "pt_name", "pt_index")) %>%
    left_join(pctFOR, by=c("RegionalID", "pt_name", "pt_index")) %>%
    left_join(baa1, by=c("RegionalID", "pt_name", "pt_index")) %>%
    left_join(intermix1, by=c("RegionalID", "pt_name", "pt_index")) %>%
    left_join(interface1, by=c("RegionalID", "pt_name", "pt_index")) %>%
    left_join(wui1, by=c("RegionalID", "pt_name", "pt_index")) %>%
    left_join(build1, by=c("RegionalID", "pt_name", "pt_index"))
  dat1$Sex[dat1$RegionalID=="2019-7709" | dat1$RegionalID=="2020-70001"] <- "F"
  return(dat1)
}


## Function for calculating variances for parameters (based on math in Blanchong et al 2006 WSB 34(3))

# m is the number of imputations (10), U should be a vector of estimated variances from each imputation
u_bar <- function(U, m) {
  ubar <- (1/m) * sum(U)
  return(ubar)
}

# Between-imputation variance
# q is vector of sample variances 
B <- function(m, q) {
  Q <- mean(q)
  b <- (1/(m-1))*sum((q-Q)^2)
  return(b)
}

# Total variance
total_var <- function(u_bar, m, B) {
  
  Tvar <- u_bar + (1+(1/m))*B
  return(Tvar)
  
}