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
          select(RegionalID:Town, n.compounds.T, pasture_60, dcad_15, BBA_15, nbuildings_60, lag_beechnuts, mix_60_250) %>%
          rename(WUI=mix_60_250,
                 DCAD=dcad_15,
                 nbuildings=nbuildings_60,
                 pasture=pasture_60,
                 basalarea=BBA_15,
                 mast=lag_beechnuts)
# scale data
dat <- within(dat, mast.s <- scale(mast))
dat <- within(dat, Age.s <- scale(Age))
dat <- within(dat, DCAD.s <- scale(DCAD))
dat[,c(8,16:21)] <- scale(dat[,c(8,16:21)])

# Calculate means
meanWUI <- mean(dat$WUI)
meanPasture <- mean(dat$pasture)
meanMast <- mean(dat$mast)
meanAge <- mean(dat$Age)
meanBBA <- mean(dat$basalarea)
meanBuild <- mean(dat$nbuildings)
meanDCAD <- mean(dat$DCAD)

# read random effects
ncomp_re <- read_csv("results/ncompounds_random_effects_coefficients.csv")

# read fixed effects
ncomp_fe <- read_csv("results/ncompT_coef-summary.csv") %>%
                rename(intercept=Intercept, 
                       Age2=I.Age.2.,
                       nbuildings=nbuildings_60,
                       DCAD=dcad_15,
                       pasture=pasture_60, 
                       basalarea=BBA_15,
                       mast=lag_beechnuts,
                       intx_encroach=nbuildings_60.dcad_15,
                       intx_beech=BBA_15.lag_beechnuts)


## effects of age
age_m <- poisson_pred_age(fixed=ncomp_fe, 
                          random=ncomp_re, 
                          sex=1, 
                          meanBuild=meanBuild, 
                          meanDCAD=meanDCAD,
                          meanPasture=meanPasture, 
                          meanMast=meanMast, 
                          meanBBA=meanBBA,
                          ageStart=min(dat$Age), 
                          ageEnd=max(dat$Age), 
                          lo=100)
age_f <- poisson_pred_age(fixed=ncomp_fe, 
                          random=ncomp_re, 
                          sex=0, 
                          meanBuild=meanBuild, 
                          meanDCAD=meanDCAD,
                          meanPasture=meanPasture, 
                          meanMast=meanMast, 
                          meanBBA=meanBBA,
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
level_freq <- dat %>% 
  select(RegionalID:Town) %>% 
  filter(pt_index==1) %>%
  distinct() %>%
  group_by(WMU) %>%
  summarize(n=n()) %>%
  mutate(freq = n / sum(n))

# mean line
pred_mean <- pred_vals %>% group_by(sex, unscAge) %>% summarize(mval = weighted.mean(exp_val, level_freq$freq))
# pred_mean <- pred_vals %>% group_by(sex, unscAge) %>% summarize(mval = mean(exp_val))

# Plot age
age_plot <- ggplot() + 
                geom_line(data=pred_vals, 
                          aes(x=unscAge, y=exp_val, group=interaction(level, sex), color=factor(sex)), alpha=0.2) +
                geom_line(data=pred_mean, aes(x=unscAge, y=mval, color=sex), size=1) +
                scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
                ylim(0,3.5) +
                ylab("Expected number of compounds") + xlab("Age (years)") +
                theme(legend.position=c(0,1),
                      legend.justification=c(0,1),
                      legend.background=element_rect(fill=NA),
                      panel.border=element_rect(color="black", fill=NA, size=0.5),
                      axis.text=element_text(size=12),
                      axis.title=element_text(size=12),
                      legend.title=element_text(size=12),
                      legend.text=element_text(size=11))

## effect of beech mast
mast_m <- poisson_pred_mast(ncomp_fe, 
                            ncomp_re, 
                            sex=1, 
                            meanBuild=meanBuild, 
                            meanDCAD=meanDCAD,
                            meanPasture=meanPasture, 
                            meanAge=meanAge, 
                            meanBBA=meanBBA,
                            mastStart=min(dat$mast), 
                            mastEnd=max(dat$mast), 
                            lo=100)

mast_f <- poisson_pred_mast(ncomp_fe, 
                            ncomp_re, 
                            sex=0, 
                            meanBuild=meanBuild, 
                            meanDCAD=meanDCAD,
                            meanPasture=meanPasture, 
                            meanAge=meanAge, 
                            meanBBA=meanBBA,
                            mastStart=min(dat$mast), 
                            mastEnd=max(dat$mast), 
                            lo=100)

pred_mast <- list(mast_m, mast_f) %>%
  reduce(full_join) %>%
  select(level, sex, Mast, exp_val) %>%
  mutate(unscMast=Mast * attr(dat$mast.s, 'scaled:scale') + attr(dat$mast.s, 'scaled:center'))
pred_mast 

# mean line
mast_mean <- pred_mast %>% group_by(sex, unscMast) %>% summarize(mval = weighted.mean(exp_val, level_freq$freq))

mast_plot <- ggplot() + 
                geom_line(data=pred_mast, 
                          aes(x=unscMast, y=exp_val, group=interaction(level, sex), color=factor(sex)), alpha=0.15) +
                geom_line(data=mast_mean, aes(x=unscMast, y=mval, color=sex), size=1) +
                scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
                ylim(0, 3.5) +
                ylab("Expected number of compounds") + xlab("Beechnut count") +
                theme(legend.position="none",
                      panel.border=element_rect(color="black", fill=NA, size=0.5),
                      axis.text.x=element_text(size=12),
                      axis.text.y=element_blank(),
                      axis.title.x=element_text(size=12),
                      axis.title.y=element_blank())

age_plot + mast_plot + plot_annotation(tag_levels="a", tag_prefix="(", tag_suffix=")")

#### Binary analyses #####
dat <- read_csv("output/binary_model_data_unscaled.csv") %>%
            select(RegionalID:bin.exp, pasture_60, dcad_15, BBA_15, nbuildings_60, lag_beechnuts, mix_60_250) %>%
                  rename(WUI=mix_60_250,
                         DCAD=dcad_15,
                         nbuildings=nbuildings_60,
                         pasture=pasture_60,
                         basalarea=BBA_15,
                         mast=lag_beechnuts)

# Scale data
dat <- within(dat, Age.s <- scale(Age))
dat <- within(dat, mast.s <- scale(mast))
dat <- within(dat, dcad.s <- scale(DCAD))

dat[,c(8,17:22)] <- scale(dat[,c(8,17:22)])

# Caluculate means
meanWUI <- mean(dat$WUI)
meanPasture <- mean(dat$pasture)
meanMast <- mean(dat$mast)
meanAge <- mean(dat$Age)
meanBBA <- mean(dat$basalarea)
meanBuild <- mean(dat$nbuildings)
meanDCAD <- mean(dat$DCAD)

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
                                                                              nbuildings=nbuildings_60,
                                                                              DCAD=dcad_15,
                                                                              pasture=pasture_60, 
                                                                              basalarea=BBA_15,
                                                                              mast=lag_beechnuts,
                                                                              intx_encroach=nbuildings_60.dcad_15,
                                                                              intx_beech=BBA_15.lag_beechnuts)
brom_fe <- read_csv("results/binaryTbromadiolone_coef-summary.csv") %>% rename(intercept=X.Intercept., 
                                                                                Age2=I.Age.2.,
                                                                                nbuildings=nbuildings_60,
                                                                                DCAD=dcad_15,
                                                                                pasture=pasture_60, 
                                                                                basalarea=BBA_15,
                                                                                mast=lag_beechnuts,
                                                                                intx_encroach=nbuildings_60.dcad_15,
                                                                                intx_beech=BBA_15.lag_beechnuts)
diph_fe <- read_csv("results/binaryTdiphacinone_coef-summary.csv") %>% rename(intercept=X.Intercept., 
                                                                              Age2=I.Age.2.,
                                                                              nbuildings=nbuildings_60,
                                                                              DCAD=dcad_15,
                                                                              pasture=pasture_60, 
                                                                              basalarea=BBA_15,
                                                                              mast=lag_beechnuts,
                                                                              intx_encroach=nbuildings_60.dcad_15,
                                                                              intx_beech=BBA_15.lag_beechnuts)
 
## Calculate expected values
brod_m <- logistic_pred_age(fixed=brod_fe, 
                        random=brod_re, 
                        compound="brodifacoum", 
                        sex=1, 
                        meanBuild=meanBuild,
                        meanDCAD=meanDCAD, 
                        meanPasture=meanPasture, 
                        meanBBA=meanBBA,
                        meanMast=meanMast, 
                        ageStart=min(dat$Age), 
                        ageEnd=max(dat$Age), 
                        lo=100)
brod_f <- logistic_pred_age(brod_fe, brod_re, compound="brodifacoum", sex=0, meanBuild=meanBuild,meanDCAD=meanDCAD, meanPasture, meanBBA, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)

brom_m <- logistic_pred_age(brom_fe, brom_re, compound="bromadiolone", sex=1, meanBuild=meanBuild,meanDCAD=meanDCAD,  meanPasture, meanBBA,meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)
brom_f <- logistic_pred_age(brom_fe, brom_re, compound="bromadiolone", sex=0, meanBuild=meanBuild,meanDCAD=meanDCAD,  meanPasture, meanBBA, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)

diph_m <- logistic_pred_age(diph_fe, diph_re, compound="diphacinone", sex=1, meanBuild=meanBuild,meanDCAD=meanDCAD,  meanPasture, meanBBA, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)
diph_f <- logistic_pred_age(diph_fe, diph_re, compound="diphacinone", sex=0, meanBuild=meanBuild,meanDCAD=meanDCAD,  meanPasture, meanBBA, meanMast, ageStart=min(dat$Age), ageEnd=max(dat$Age), lo=100)

# Combine into single data frame
pred_vals <- list(brod_m, brod_f, brom_m, brom_f, diph_m, diph_f) %>%
                reduce(full_join) %>%
                select(compound:Age, exp_val) %>%
                mutate(unscAge=Age * attr(dat$Age.s, 'scaled:scale') + attr(dat$Age.s, 'scaled:center'))
pred_vals 

level_freq <- dat %>% select(RegionalID:Town) %>% 
                 filter(pt_index==1) %>%
                 distinct() %>%
                 group_by(WMU) %>%
                 summarize(n=n()) %>%
                 mutate(freq = n / sum(n))

# mean line
# pred_mean <- pred_vals %>% group_by(compound, sex, unscAge) %>% summarize(mval = mean(exp_val))
pred_mean <- pred_vals %>% group_by(compound, sex, unscAge) %>% summarize(mval = weighted.mean(exp_val, level_freq$freq))

pred_mean %>% group_by(compound, sex) %>% summarize(maxprob=max(mval))
pred_mean %>% group_by(compound, sex) %>% summarize(probrange=range(mval))
pred_mean %>% group_by(compound, sex) %>% summarize(probsd=sd(mval))
pred_mean %>% group_by(compound, sex) %>% summarize(meanprob=mean(mval))

pred_mean %>% filter(compound=="brodifacoum" & sex=="Female" & mval>=0.840)

## Plot

ggplot() + 
  geom_line(data=pred_vals, 
              aes(x=unscAge, y=exp_val, group=interaction(level, sex), color=factor(sex)), alpha=0.1) + 
  geom_line(data=pred_mean, aes(x=unscAge, y=mval, color=sex), size=1) +   
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  facet_grid(compound~.) +
  ylim(0,1) +
  ylab("Probability of exposure") + xlab("Age (years)") +
  theme(legend.position="bottom",
        # legend.position=c(0,0),
        # legend.justification=c(0,0),
        legend.background=element_rect(fill=NA),
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        axis.text=element_text(size=12),
        strip.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))

# ggplot() +geom_line(data=pred_mean, aes(x=unscAge, y=mval, group=interaction(compound,sex), color=sex)) #+ facet_grid(.~compound)
brod_m <- logistic_pred_dcad(fixed=brod_fe, 
                             random=brod_re, 
                             compound="brodifacoum",
                             sex=1, 
                             meanWUI, 
                             meanPasture,
                             meanBBA,
                             meanAge, 
                             mastStart=min(dat$lag_beechnuts), 
                             mastEnd=max(dat$lag_beechnuts), 
                             lo=100)

brod_f <- logistic_pred_dcad(fixed=brod_fe, 
                             random=brod_re, 
                             compound="brodifacoum",
                             sex=0, 
                             meanBuild=meanBuild,
                             meanDCAD=meanDCAD, 
                             meanPasture=meanPasture, 
                             meanBBA=meanBBA,
                             meanAge=meanAge, 
                             mastStart=min(dat$mast), 
                             mastEnd=max(dat$mast), 
                             lo=100)


brom_m <- logistic_pred_mast(fixed=brom_fe, 
                            random=brom_re, 
                            compound="bromadiolone",
                            sex=1, 
                            meanBuild=meanBuild,
                            meanDCAD=meanDCAD, 
                            meanPasture=meanPasture, 
                            meanBBA=meanBBA,
                            meanAge=meanAge, 
                            mastStart=min(dat$mast), 
                            mastEnd=max(dat$mast), 
                            lo=100)

brom_f <- logistic_pred_mast(fixed=brom_fe, 
                             random=brom_re, 
                             compound="bromadiolone",
                             sex=0, 
                             meanBuild=meanBuild,
                             meanDCAD=meanDCAD, 
                             meanPasture=meanPasture, 
                             meanBBA=meanBBA,
                             meanAge=meanAge, 
                             mastStart=min(dat$mast), 
                             mastEnd=max(dat$mast), 
                             lo=100)

diph_m <- logistic_pred_mast(fixed=diph_fe, 
                             random=diph_re, 
                             compound="diphacinone",
                             sex=1, 
                             meanWUI, 
                             meanPasture,
                             meanBBA,
                             meanAge, 
                             mastStart=min(dat$lag_beechnuts), 
                             mastEnd=max(dat$lag_beechnuts), 
                             lo=100)

diph_f <- logistic_pred_mast(fixed=diph_fe, 
                             random=diph_re, 
                             compound="diphacinone",
                             sex=0, 
                             meanWUI, 
                             meanPasture,
                             meanBBA,
                             meanAge, 
                             mastStart=min(dat$lag_beechnuts), 
                             mastEnd=max(dat$lag_beechnuts), 
                             lo=100)


pred_mast <- list(brod_m, brod_f, brom_m, brom_f, diph_m, diph_f) %>%
  reduce(full_join) %>%
  select(compound, level, sex, Mast, exp_val) %>%
  mutate(unscMast=Mast * attr(dat$mast.s, 'scaled:scale') + attr(dat$mast.s, 'scaled:center'))
pred_mast 

# mean line
mast_mean <- pred_mast %>% group_by(sex, compound, unscMast) %>% summarize(mval=mean(exp_val))

ggplot() + 
  geom_line(data=pred_mast, 
            aes(x=unscMast, y=exp_val, group=interaction(level, sex), color=factor(sex)), alpha=0.2) +
  geom_line(data=mast_mean, aes(x=unscMast, y=mval, color=sex), size=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  ylab("Probability of exposure") + xlab("Total beech nuts (previous year)") +
  facet_grid(compound~.) +
  theme()


