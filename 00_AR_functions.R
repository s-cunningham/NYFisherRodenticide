## Functions for rodenticide

library(tidyverse)


exp_val_calc <- function(fixed, random, compound, sex, meanWUI, meanPasture, meanMast, ageStart, ageEnd, lo) {
  
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
      
      exp_val[i] <- inv.logit(fixed$intercept[1] + fixed$SexM[1]*sex + fixed$Age[1]*age_iter[i] + fixed$Age2[1]*(age_iter[i]^2) +
                                fixed$WUI[1]*meanWUI + fixed$pasture[1]*meanPasture + fixed$mast[1]*meanMast + random$REval[j])
      
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



# create empty data frame to store
# full_vals <- matrix(NA, ncol=100+1, nrow=nrow(diph_re))
# full_vals <- data.frame()
# full_vals <- bind_rows(full_vals, as.data.frame(t(1:100)))
# 
# full_vals$level <- diph_re$grp[1:43]
# 
# # crete a sequence of values to estimate for age
# age_iter <- seq(0.5, 8.5, length.out=100)
