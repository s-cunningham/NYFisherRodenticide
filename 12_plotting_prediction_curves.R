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
dat <- read_csv("output/model_data_notscaled.csv") %>%
          select(RegionalID:Town, n.compounds.T, pasture_60, laggedBMI_30, wui_60_100) %>%
          rename(WUI=wui_60_100,
                 pasture=pasture_60, 
                 mast=laggedBMI_30)
# scale data
dat <- within(dat, mast.s <- scale(mast))
dat <- within(dat, Age.s <- scale(Age))
dat[,c(8,16:18)] <- scale(dat[,c(8,16:18)])

# Calculate means
meanWUI <- mean(dat$WUI)
meanPasture <- mean(dat$pasture)
meanMast <- mean(dat$mast)
meanAge <- mean(dat$Age)

# read random effects
ncomp_re <- read_csv("results/ncompounds_random_effects_coefficients.csv")

# read fixed effects
ncomp_fe <- read_csv("results/ncompT_coef-summary.csv") %>%
                rename(WUI=wui_60_100,
                       pasture=pasture_60, 
                       mast=laggedBMI_30,
                       intercept=Intercept, 
                       Age2=I.Age.2.)
## effects of age
age_m <- poisson_pred_age(fixed=ncomp_fe, 
                            random=ncomp_re, 
                            sex=1, 
                            meanWUI, 
                            meanPasture, 
                            meanMast, 
                            ageStart=min(dat$Age), 
                            ageEnd=max(dat$Age), 
                            lo=100)
age_f <- poisson_pred_age(fixed=ncomp_fe, 
                          random=ncomp_re, 
                          sex=0, 
                          meanWUI, 
                          meanPasture, 
                          meanMast, 
                          ageStart=min(dat$Age), 
                          ageEnd=max(dat$Age), 
                          lo=100)

# Combine into single data frame
pred_vals <- list(age_m, age_f) %>%
  reduce(full_join) %>%
  select(level, sex, Age, exp_val) %>%
  mutate(unscAge=Age * attr(dat$Age.s, 'scaled:scale') + attr(dat$Age.s, 'scaled:center'))
pred_vals 

# mean line
pred_mean <- pred_vals %>% group_by(sex, unscAge) %>% summarize(mval = mean(exp_val))

# Plot age
age_plot <- ggplot() + 
                geom_line(data=pred_vals, 
                          aes(x=unscAge, y=exp_val, group=interaction(level, sex), color=factor(sex)), alpha=0.1) +
                geom_line(data=pred_mean, aes(x=unscAge, y=mval, color=sex), size=1) +
                scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
                ylab("Expected number of compounds") + xlab("Age (years)") +
                theme(legend.position=c(0,1),
                      legend.justification=c(0,1),
                      legend.background=element_rect(fill=NA),
                      panel.border=element_rect(color="black", fill=NA, size=0.5),
                      axis.text=element_text(size=11),
                      axis.title=element_text(size=11),
                      legend.title=element_text(size=11),
                      legend.text=element_text(size=10))

## effect of beech mast

mast_m <- poisson_pred_mast(ncomp_fe, 
                            ncomp_re, 
                            sex=1, 
                            meanWUI, 
                            meanPasture, 
                            meanAge, 
                            mastStart=min(dat$mast), 
                            mastEnd=max(dat$mast), 
                            lo=100)

mast_f <- poisson_pred_mast(ncomp_fe, 
                            ncomp_re, 
                            sex=0, 
                            meanWUI, 
                            meanPasture, 
                            meanAge, 
                            mastStart=min(dat$mast), 
                            mastEnd=max(dat$mast), 
                            lo=100)

pred_mast <- list(mast_m, mast_f) %>%
  reduce(full_join) %>%
  select(level, sex, Mast, exp_val) %>%
  mutate(unscMast=Mast * attr(dat$mast.s, 'scaled:scale') + attr(dat$mast.s, 'scaled:center'))
pred_mast 

# mean line
mast_mean <- pred_mast %>% group_by(sex, unscMast) %>% summarize(mval = mean(exp_val))

mast_plot <- ggplot() + 
                geom_line(data=pred_mast, 
                          aes(x=unscMast, y=exp_val, group=interaction(level, sex), color=factor(sex)), alpha=0.1) +
                geom_line(data=mast_mean, aes(x=unscMast, y=mval, color=sex), size=1) +
                scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
                ylab("Expected number of compounds") + xlab("Spatially weighted beech mast index") +
                theme(legend.position="none",
                      panel.border=element_rect(color="black", fill=NA, size=0.5),
                      axis.text.x=element_text(size=11),
                      axis.text.y=element_blank(),
                      axis.title.x=element_text(size=11),
                      axis.title.y=element_blank())

age_plot + mast_plot

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
brod_m <- logistic_pred(brod_fe, 
                        brod_re, 
                        compound="brodifacoum", 
                        sex=1, 
                        meanWUI, 
                        meanPasture, 
                        meanMast, 
                        ageStart=min(dat$Age), 
                        ageEnd=max(dat$Age), 
                        lo=100)
brod_f <- logistic_pred(brod_fe, brod_re, compound="brodifacoum", sex=0, meanWUI, meanPasture, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)

brom_m <- logistic_pred(brom_fe, brom_re, compound="bromadiolone", sex=1, meanWUI, meanPasture, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)
brom_f <- logistic_pred(brom_fe, brom_re, compound="bromadiolone", sex=0, meanWUI, meanPasture, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)

diph_m <- logistic_pred(diph_fe, diph_re, compound="diphacinone", sex=1, meanWUI, meanPasture, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)
diph_f <- logistic_pred(diph_fe, diph_re, compound="diphacinone", sex=0, meanWUI, meanPasture, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)

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

