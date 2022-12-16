## Plotting prediction curves
## 2022-09-14

library(boot)
library(tidyverse)
library(patchwork)

## Source function files
source("00_AR_functions.R")

## Set ggplot theme
theme_set(theme_classic())

#### Number of compounds ####
# dat <- read_csv("output/model_data.csv")

#### Binary analyses #####
dat <- read_csv("output/binary_model_data_unscaled.csv") %>%
            select(RegionalID:bin.exp, pasture_60, laggedBMI_30, wui_60_100)

# Scale data
dat <- within(dat, Age.s <- scale(Age))
dat[,c(8,17:19)] <- scale(dat[,c(8,17:19)])

# Caluculate means
meanWUI <- mean(dat$wui_60_100)
meanPasture <- mean(dat$pasture_60)
meanMast <- mean(dat$laggedBMI_30)

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

## Calculate expected values
brod_m <- exp_val_calc(brod_fe, brod_re, compound="brodifacoum", sex=1, meanWUI, meanPasture, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)
brod_f <- exp_val_calc(brod_fe, brod_re, compound="brodifacoum", sex=0, meanWUI, meanPasture, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)

brom_m <- exp_val_calc(brom_fe, brom_re, compound="bromadiolone", sex=1, meanWUI, meanPasture, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)
brom_f <- exp_val_calc(brom_fe, brom_re, compound="bromadiolone", sex=0, meanWUI, meanPasture, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)

diph_m <- exp_val_calc(diph_fe, diph_re, compound="diphacinone", sex=1, meanWUI, meanPasture, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)
diph_f <- exp_val_calc(diph_fe, diph_re, compound="diphacinone", sex=0, meanWUI, meanPasture, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)

# Combine into single data frame
pred_vals <- list(brod_m, brod_f, brom_m, brom_f, diph_m, diph_f) %>%
                reduce(full_join) %>%
                select(compound:Age, exp_val) %>%
                mutate(unscAge=Age * attr(dat$Age.s, 'scaled:scale') + attr(dat$Age.s, 'scaled:center'))
pred_vals 

pred_vals %>% group_by(compound, sex, level, Age) %>% count()

# mean line
pred_mean <- pred_vals %>% group_by(compound, sex, unscAge) %>% summarize(mval = mean(exp_val))


## Plot

ggplot() + 
  geom_line(data=pred_vals, 
              aes(x=unscAge, y=exp_val, group=interaction(level, sex), color=factor(sex)), alpha=0.15) + 
  geom_line(data=pred_mean, aes(x=unscAge, y=mval, color=sex), size=1) +   
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  facet_grid(compound~.) +
  ylab("Probability of exposure") + xlab("Age (years)") +
  theme(legend.position=c(0,0),
        legend.justification=c(0,0),
        legend.background=element_rect(fill=NA),
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        axis.text=element_text(size=11),
        strip.text=element_text(size=11),
        axis.title=element_text(size=11),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10))


ggplot() +geom_line(data=pred_mean, aes(x=unscAge, y=mval, color=sex)) + facet_grid(.~compound)

