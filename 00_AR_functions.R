## Functions for rodenticide

library(tidyverse)


exp_val_calc <- function(fixed, random, compound, meanWUI, meanPasture, meanMast, ageStart, ageEnd, lo) {
  
  # crete a sequence of values to estimate for age
  age_iter <- seq(ageStart, ageEnd, length.out=lo)
  
  # create empty data frame to store
  full_vals <- matrix(NA, ncol=lo+1, nrow=nrow(random))
  full_vals[,1] <- random$grp
  
  # Loop over each level of random effect
  for (j in 1:nrow(random)) {
    
    # Empty vector for each level of random effect
    exp_val <- c()
    
    # Loop over values of age
    for (i in 1:length(age_iter)) {
      
      exp_val[i] <- inv.logit(fixed$intercept[1] + fixed$SexM[1]*sex + fixed$Age[1]*age_iter[i] + fixed$Age2[1]*(age_iter[i]^2) +
                                fixed$WUI[1]*meanWUI + fixed$pasture[1]*meanPasture + fixed$mast[1]*meanMast + REval[j])
      
      # Does RE need to be weighted somehow?
      # level_probs <- random %>% group_by(grp) %>% summarize(wmu_prop=n()/nrow(random))
      
    }
    
    # Save the expected values to data frame
    full_vals[,i] <- exp_val
    
  }
  
  # Reorganize data into data frame
  full_vals <- full_vals |> as.data.frame()
  names(full_vals)[2:ncol(full_vals)] <- 1:(ncol(full_vals)-1)
  names(full_vals)[1] <- "level"
  full_vals <- pivot_longer(full_vals, 2:ncol(full_vals), names_to="index", values_to="exp_val") %>%
    mutate_at(c('index', 'exp_val'), as.numeric)
  
  # add column for sex
  full_vals$compound <- compound
  full_vals$sex <- ifelse(sex==1, "Male", "Female")
  
  full_vals <- full_vals %>% select(compound, sex, level:exp_val)
  
  return(full_vals)
  
}

