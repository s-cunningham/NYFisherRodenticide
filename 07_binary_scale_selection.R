## Logistic models
library(tidyverse)
library(lme4)
library(MuMIn)
library(caret)
library(broom)

# Load function file
source("00_AR_functions.R")

# Load covariate information
covar <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>% 
            rename(BBA=baa) %>%
            select(RegionalID:rand_y,buffsize:lag_beechnuts)

# Load binary for each compound, combine covariates
brod <- read_csv("output/binary_brodifacoum.csv") %>% 
          left_join(covar, by="RegionalID") %>% 
          select(RegionalID, pt_name:rand_y, year:bin.exp.ntr,buffsize:lag_beechnuts) %>%
          covar_org()%>%
          select(-c(BMI_15, BMI_30, BMI_60))
brom <- read_csv("output/binary_bromadiolone.csv") %>% 
          left_join(covar, by="RegionalID") %>% 
          select(RegionalID, pt_name:rand_y, year:bin.exp.ntr,buffsize:lag_beechnuts) %>%
          covar_org()%>%
          select(-c(BMI_15, BMI_30, BMI_60))
diph <- read_csv("output/binary_diphacinone.csv") %>% 
          left_join(covar, by="RegionalID") %>% 
          select(RegionalID, pt_name:rand_y, year:bin.exp.ntr,buffsize:lag_beechnuts) %>%
          covar_org()%>%
          select(-c(BMI_15, BMI_30, BMI_60))
dico <- read_csv("output/binary_dicoumarol.csv") %>% 
          left_join(covar, by="RegionalID") %>% 
          select(RegionalID, pt_name:rand_y, year:bin.exp.ntr,buffsize:lag_beechnuts) %>%
          covar_org() %>%
          select(-c(BMI_15, BMI_30, BMI_60))

# Scale
brod[,c(8,17:78)] <- scale(brod[,c(8,17:78)])
brom[,c(8,17:78)] <- scale(brom[,c(8,17:78)])
diph[,c(8,17:78)] <- scale(diph[,c(8,17:78)])
dico[,c(8,17:78)] <- scale(dico[,c(8,17:78)])

# Load models
source("00_logistic_models.R")

#### Run models  ####

## Brodifacoum
# Age sex
as.models <- lapply(agesex_forms, FUN=glmer, data=brod, family="binomial",
                    control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))
model.list <- model.sel(as.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brod_age-sex_selection.csv")

# Full models
brd.models <- lapply(brod_models, FUN=glmer, data=brod, family='binomial',
                    control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))
model.list <- model.sel(brd.models )
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brod_full_selection.csv")


## Bromadiolone
# Age sex
as.models <- lapply(agesex_forms, FUN=glmer, data=brom, family="binomial",
                    control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))
model.list <- model.sel(as.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brom_age-sex_selection.csv")

# Full models
brd.models <- lapply(brod_models, FUN=glmer, data=brom, family='binomial',
                     control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))
model.list <- model.sel(brd.models )
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brom_full_selection.csv")



## Diphacinone
# Age sex
as.models <- lapply(agesex_forms, FUN=glmer, data=diph, family="binomial",
                    control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))
model.list <- model.sel(as.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/diph_age-sex_selection.csv")

# Full models
brd.models <- lapply(brod_models, FUN=glmer, data=diph, family='binomial',
                     control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))
model.list <- model.sel(brd.models )
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/diph_full_selection.csv")

## Dicoumarol
# Age sex
as.models <- lapply(agesex_forms, FUN=glmer, data=dico, family="binomial",
                    control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))
model.list <- model.sel(as.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/dico_age-sex_selection.csv")


