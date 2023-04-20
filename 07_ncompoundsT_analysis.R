## Running models
## 2022-09-08, updated 2023-01-11

library(tidyverse)
library(MuMIn)
library(glmmTMB)
library(DHARMa)
library(caret)
library(performance)

options(scipen=999, digits=3)
set.seed(123)

#### Parallel processing ####
nt <- min(parallel::detectCores(),4)

#### Read in data ####
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>% rename(BBA=baa)

#### Analysis Set-up ####

# Make random effects factors
dat$WMU <- as.factor(dat$WMU)
dat$WMUA_code <- as.factor(dat$WMUA_code)
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

## Landscape metrics
lsm1 <- dat %>% select(RegionalID, pt_name, pt_index, buffsize, ai:lsi) %>% 
  distinct() %>% 
  group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(ai:lsi), values_fn=unique) 

#### Set up data to run for each combination of covariates ####
dat1 <- dat %>% select(RegionalID:Town,n.compounds.T,beechnuts,lag_beechnuts) %>% distinct() %>%
          left_join(pctAG, by=c("RegionalID", "pt_name", "pt_index")) %>%
          left_join(pctFOR, by=c("RegionalID", "pt_name", "pt_index")) %>%
          left_join(baa1, by=c("RegionalID", "pt_name", "pt_index")) %>%
          left_join(intermix1, by=c("RegionalID", "pt_name", "pt_index")) %>%
          left_join(interface1, by=c("RegionalID", "pt_name", "pt_index")) %>%
          left_join(wui1, by=c("RegionalID", "pt_name", "pt_index")) %>%
          left_join(lsm1, by=c("RegionalID", "pt_name", "pt_index")) %>%
          left_join(build1, by=c("RegionalID", "pt_name", "pt_index"))
dat1$Sex[dat1$RegionalID=="2019-7709" | dat1$RegionalID=="2020-70001"] <- "F"
# write_csv(dat1, "output/model_data_notscaled.csv")

## Scale and center variables
dat1[,c(8,16:107)] <- scale(dat1[,c(8,16:107)])
dat1$mast <- as.factor(ifelse(dat1$year==2018 | dat1$year==2020, "fail", "mast"))

# correlation coefficient
cormat <- cor(dat1[,c(16:104)]) |> as.data.frame()
# write_csv(cormat, "output/correlation_matrix.csv")

#### Read in formulas ####
source("00_model_lists.R")

#### Run models with glmmTMB ####
## Model selection for scale
# Age and sex
agesex.models <- lapply(agesex_formulae, FUN=glmmTMB, data=dat1, 
                        family=compois(link = "log"), control=glmmTMBControl(parallel=nt, 
                                                                             profile=TRUE, 
                                                                             optCtrl=list(iter.max=1e20,eval.max=1e20),
                                                                             optimizer=optim, 
                                                                             optArgs=list(method="BFGS"))) 
model.list <- model.sel(agesex.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
write_csv(model_tab, "output/model_selection/agesex_model_selection_table.csv")


all_ag <- c(ag_formulae, crops_formulae, pasture_form)

# Agriculture
ag.models <- lapply(all_ag, FUN=glmmTMB, data=dat1, 
                    family=compois(link = "log"), 
                    control=glmmTMBControl(parallel=nt, 
                                           profile=TRUE, 
                                           optCtrl=list(iter.max=1e11,eval.max=1e11),
                                           optimizer=optim, 
                                           optArgs=list(method="BFGS")))
model.list <- model.sel(ag.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/totalag_selection_table.csv")

# Beech basal area
beech.models <- lapply(beech_formulae, FUN=glmmTMB, data=dat1, 
                       family=compois(link = "log"), control=glmmTMBControl(parallel=nt,
                                                                            profile=TRUE, 
                                                                            optCtrl=list(iter.max=1e11,eval.max=1e11),
                                                                            optimizer=optim, 
                                                                            optArgs=list(method="BFGS"))) 
model.list <- model.sel(beech.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/beech_model_selection_table.csv")

# Forest metrics
forest.models <- lapply(forest_formulae, FUN=glmmTMB, data=dat1, 
                        family=compois(link = "log"), control=glmmTMBControl(parallel=nt, 
                                                                             profile=TRUE, 
                                                                             optCtrl=list(iter.max=1e20,eval.max=1e20),
                                                                             optimizer=optim, 
                                                                             optArgs=list(method="BFGS"))) 
model.list <- model.sel(forest.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/forest_model_selection_table.csv")

# Buildings
build.models <- lapply(build_formulae, FUN=glmmTMB, data=dat1, 
                       family=compois(link = "log"), control=glmmTMBControl(parallel=nt, 
                                                                            profile=TRUE, 
                                                                            optCtrl=list(iter.max=1e19,eval.max=1e19),
                                                                            optimizer=optim, 
                                                                            optArgs=list(method="BFGS")))
model.list <- model.sel(build.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/building_model_selection_table.csv")


all_wui <- c(interface_formulae, intermix_formulae, wui_formulae)
## WUI
wui.models <- lapply(all_wui, FUN=glmmTMB, data=dat1, 
                     family=compois(link = "log"), control=glmmTMBControl(parallel=nt,
                                                                          profile=TRUE, 
                                                                          optCtrl=list(iter.max=1e20,eval.max=1e20), 
                                                                          optimizer=optim, 
                                                                          optArgs=list(method="BFGS"))) 
model.list <- model.sel(wui.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
write_csv(model_tab, "output/model_selection/wui_model_selection_table.csv")


# Check small set of "global models"
global.models <- lapply(global_formulae, FUN=glmmTMB, data=dat1, 
                       family=compois(link = "log"), 
                       control=glmmTMBControl(parallel=nt, 
                                              profile=TRUE,
                                              optCtrl=list(iter.max=1e11,eval.max=1e11),
                                              optimizer=optim, 
                                              optArgs=list(method="BFGS")))
model.list <- model.sel(global.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/global_model_selection_table.csv")


system.time(m1 <- glmmTMB(n.compounds.T ~ Sex + Age + I(Age^2) + face_15_100 + (1|RegionalID), data=dat1,
        family=compois(link = "log"), 
        control=glmmTMBControl(parallel=nt, 
                               profile=TRUE, 
                               optCtrl=list(iter.max=1e19,eval.max=1e19))))

#### Running iteration models ####

## Loop over each set of random points
m_est <- m_stderr <- m_zscore <- pct2.5 <- pct97.5 <- m_ranef <- data.frame()
perf <- matrix(NA, ncol=11, nrow=10)
perf[,1] <- 1:10

kappa <- matrix(NA, ncol=6, nrow=10)
kappa[,1] <- 1:10

classstats <- data.frame()
overallstats <- data.frame()
ranef_coef <- data.frame()
vif <- list()

system.time(
# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- dat1[dat1$pt_index==i,]
  
  # Run model with deltaAICc < 2

  m1_pt <- glmmTMB(n.compounds.T ~ Sex + Age + I(Age^2) + mix_30_250 + pasture_30 + BBA_15 * lag_beechnuts + (1|WMU), data=pt, 
                      family=compois(link = "log"), control=glmmTMBControl(parallel=nt))
  
  vif[[i]] <- check_collinearity(m1_pt, "count")
  m1s <- summary(m1_pt)
  
  # save averaged confidence intervals
  pct2.5 <- rbind(pct2.5, t(confint(m1_pt))[1,1:(nrow(confint(m1_pt))-1)])
  pct97.5 <- rbind(pct97.5, t(confint(m1_pt))[2,1:(nrow(confint(m1_pt))-1)])
  
  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)$cond[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)$cond[,2])
  m_zscore <- rbind(m_zscore, coef(m1s)$cond[,3])
  
  # broom.mixed::tidy(m1_pt)
  
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
    m1_cv <- glmmTMB(n.compounds.T ~ Sex + Age + I(Age^2) + mix_60_250 + pasture_30 + BBA_15 * lag_beechnuts + (1|WMU), data=pt,
                     family=compois(link = "log"), control=glmmTMBControl(parallel=nt))
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
    
    all_confusion <- confusionMatrix(
      # predictions then true values
      data = factor(round(testTable$pred), levels=0:5),
      reference = factor(testTable$n.compounds.T, levels=0:5),
    )
    
    classcm <- as_tibble(all_confusion$byClass)
    classcm$iteration <- i
    classcm$fold <- j
    classstats <- bind_rows(classstats, classcm)
    
    overall <- as_tibble(t(all_confusion$overall))
    overall$iteration <- i
    overall$fold <- j
    overallstats <- bind_rows(overallstats, overall)
    
    kappa[i,j+1] <- all_confusion$overall[[2]]

  }

  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- c(row.names(m1s$coefficients$cond))
    names(m_stderr) <- c(row.names(m1s$coefficients$cond))
    names(m_zscore) <- c(row.names(m1s$coefficients$cond))
    names(pct2.5) <- c(row.names(m1s$coefficients$cond))
    names(pct97.5) <- c(row.names(m1s$coefficients$cond))
    names(m_ranef) <- c("RE_WMU")
  }
  
})

# Calculate average performance metric
perf <- as.data.frame(perf)
names(perf) <- c("Iter", "RMSE1", "RMSE2", "RMSE3", "RMSE4", "RMSE5", "Accur1", "Accur2", "Accur3", "Accur4", "Accur5")
perf2 <- perf %>% as_tibble() %>% 
          mutate(RMSE=(RMSE1+RMSE2+RMSE3+RMSE4+RMSE5)/5,
                 Accuracy=(Accur1+Accur2+Accur3+Accur4+Accur5)/5) %>%
          select(Iter, RMSE, Accuracy)
mean(perf2$Accuracy)
range(perf2$Accuracy)

# Save averaged values for each level of random effect
re <- ranef_coef %>% 
  group_by(grp) %>% 
  summarize(REval=mean(condval), REsd=mean(condsd), 
            RErangeL=range(condval)[1], RErangeH=range(condval)[2]) 
write_csv(re, "results/ncompounds_random_effects_coefficients.csv")

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
zscore_avg <- colMeans(m_zscore[sapply(m_zscore, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)
ranef_avg <- as.data.frame(colMeans(m_ranef[sapply(m_ranef, is.numeric)])) %>% rownames_to_column("RE")
names(ranef_avg)[2] <- "variance"

mean(m_est[,2])
sd(m_est[,2])

# Combine and clean up data frame
coef_summary <- bind_rows(coef_avg, stderr_avg, zscore_avg, pct2.5_avg, pct97.5_avg)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "z-score", "2.5CI", "97.5CI")
coef_summary <- data.frame(coef=coefs, coef_summary)
names(coef_summary)[2] <- "Intercept"

# Write to file
write_csv(coef_summary, "results/ncompT_coef-summary.csv")
write_csv(ranef_avg, "results/ncompT_coef-random-effects.csv")

# Calculate class accuracy
classstats <- classstats %>% as_tibble %>%
                mutate(Accuracy = Sensitivity*Prevalence + (Specificity*(1-Prevalence)))
classstats$group <- rep(0:5, 10*5)

class_accu <- classstats %>% group_by(iteration, group) %>% 
  summarize(class_acc=mean(`Balanced Accuracy`)) %>% 
  pivot_wider(names_from=group, values_from=class_acc)

write_csv(classstats, "results/ncompounds_class-stats.csv")
write_csv(overallstats,  "results/ncompounds_overall_stats.csv")




