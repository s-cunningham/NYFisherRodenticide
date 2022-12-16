## Plotting prediction curves
## 2022-09-14

library(boot)
library(tidyverse)
library(patchwork)

#### Number of compounds ####
dat <- read_csv("output/model_data.csv")

#### Binary analyses #####
dat <- read_csv("output/binary_model_data_unscaled.csv") %>%
            select(RegionalID:bin.exp, pasture_60, laggedBMI_30, wui_60_100)

dat[,c(8,17:19)] <- scale(dat[,c(8,17:19)])

# Subset by compound
brod <- dat[dat$compound=="Brodifacoum",]
brom <- dat[dat$compound=="Bromadiolone",]
diph <- dat[dat$compound=="Diphacinone",]

# Read random effects coefficients
brod_re <- read_csv("results/brodifacoum_random_effects_coefficients.csv") 
brom_re <- read_csv("results/bromadiolone_random_effects_coefficients.csv")
diph_re <- read_csv("results/diphacinone_random_effects_coefficients.csv")

# Read fixed effects coefficients
brod_fe <- read_csv("results/binaryTbrodifacoum_coef-summary.csv") %>% rename(intercept=X.Intercept., 
                                                                              Age2=I.Age.2.,
                                                                              WUI=wui_60_100,
                                                                              pasture=pasture_60,
                                                                              mast=laggedBMI_30)
brom_fe <- read_csv("results/binaryTbromadiolone_coef-summary.csv") %>% rename(intercept=X.Intercept., 
                                                                               Age2=I.Age.2.,
                                                                               WUI=wui_60_100,
                                                                               pasture=pasture_60,
                                                                               mast=laggedBMI_30)
diph_fe <- read_csv("results/binaryTdiphacinone_coef-summary.csv") %>% rename(intercept=X.Intercept., 
                                                                              Age2=I.Age.2.,
                                                                              WUI=wui_60_100,
                                                                              pasture=pasture_60,
                                                                              mast=laggedBMI_30)

# Create function to build data for prediction curve

exp_val_calc <- function(fixed, random, sex, meanWUI, meanPasture, meanMast, ageStart, ageEnd, inc) {}

age_iter <- seq(ageStart, ageEnd, inc)

level_probs <- random %>% group_by(grp) %>% summarize(wmu_prop=n()/nrow(random))

for (i in 1:length(age_iter)) {
  
  exp_val <- inv.logit(fixed$intercept[1] + fixed$SexM[1]*sex + fixed$Age[1]*age_iter[i] + fixed$Age2[1]*(age_iter[i]^2) +
                         fixed$WUI[1]*meanWUI + fixed$pasture[1]*meanPasture + fixed$mast[1]*meanMast + REval[j]*)
  
}



