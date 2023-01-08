## Functions for rodenticide

library(tidyverse)

# Predict exposure by age
logistic_pred_age <- function(fixed, random, compound, sex, meanBuild, meanDCAD, meanPasture, meanBBA, meanMast, ageStart, ageEnd, lo) {
  
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
                                fixed$nbuildings[1]*meanBuild + 
                                fixed$DCAD[1]*meanDCAD +
                                fixed$pasture[1]*meanPasture + 
                                fixed$intx_encroach[1]*meanDCAD*meanBuild +
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

# Predict beech mast count
logistic_pred_mast <- function(fixed, random, compound, sex, meanBuild, meanDCAD, meanPasture, meanBBA, meanAge, mastStart, mastEnd, lo) {
  
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
                                fixed$nbuildings[1]*meanBuild + 
                                fixed$DCAD[1]*meanDCAD +
                                fixed$intx_encroach[1]*meanDCAD*meanBuild +
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


logistic_pred_intx <- function(fixed, random, compound, sex, meanWUI, meanPasture, meanBBA, meanAge, mastStart, mastEnd, lo) {
  
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
                                fixed$bba[1]*meanBBA +
                                fixed$mast[1]*mast_iter[i] + 
                                fixed$bba_mast[1]*meanBBA*mast_iter[i] +
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




poisson_pred_age <- function(fixed, random, sex, meanWUI, meanPasture, meanMast, meanBBA, ageStart, ageEnd, lo) {
  
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
      
      exp_val[i] <- exp(fixed$intercept[1] + fixed$SexM[1]*sex + fixed$Age[1]*age_iter[i] + fixed$Age2[1]*(age_iter[i]^2) +
                          fixed$WUI[1]*meanWUI + fixed$pasture[1]*meanPasture + fixed$basalarea[1]*meanBBA + 
                          fixed$mast[1]*meanMast + fixed$interaction[1]*meanMast*meanBBA + random$REval[j])
      
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

poisson_pred_mast <- function(fixed, random, sex, meanWUI, meanPasture, meanAge, meanBBA, mastStart, mastEnd, lo) {
  
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
                          fixed$interaction[1]*meanBBA*mast_iter[i] +
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
  
  full_vals$Mast <- rep(mast_iter, 55)
  
  # add column for sex
  full_vals$sex <- ifelse(sex==1, "Male", "Female")
  
  full_vals <- full_vals %>% select(sex, level, Mast, index, exp_val)
  # full_vals <- left_join(full_vals, level_probs, by="level")
  
  return(full_vals)
  
}











# create empty data frame to store
# full_vals <- matrix(NA, ncol=100+1, nrow=nrow(diph_re))
# full_vals <- data.frame()
# full_vals <- bind_rows(full_vals, as.data.frame(t(1:100)))
# 
# full_vals$level <- diph_re$grp[1:43]
# 
# # crete a sequence of values to estimate for age
# age_iter <- seq(0.5, 8.5, length.out=100)
