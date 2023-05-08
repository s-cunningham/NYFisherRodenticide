## Logistic models
library(tidyverse)
library(lme4)
library(MuMIn)
library(caret)
library(broom)

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
brom <- read_csv("output/binary_bromadiolone.csv")%>% 
          left_join(covar, by="RegionalID") %>% 
          select(RegionalID, pt_name:rand_y, year:bin.exp.ntr,buffsize:lag_beechnuts) %>%
          covar_org()%>%
          select(-c(BMI_15, BMI_30, BMI_60))
diph <- read_csv("output/binary_diphacinone.csv")%>% 
          left_join(covar, by="RegionalID") %>% 
          select(RegionalID, pt_name:rand_y, year:bin.exp.ntr,buffsize:lag_beechnuts) %>%
          covar_org()%>%
          select(-c(BMI_15, BMI_30, BMI_60))
dico <- read_csv("output/binary_dicoumarol.csv")%>% 
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

# agriculture
ag.models <- lapply(ag_formulae, FUN=glmer, data=brod, family='binomial',
                    control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))
model.list <- model.sel(ag.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/diph_ag_selection.csv")

# beech
beech.models <- lapply(beech_formulae,FUN=glmer, data=brod, family='binomial',
                       control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5))) 
model.list <- model.sel(beech.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brod_beech_selection.csv")

# forest
forest.models <- lapply(forest_formulae, FUN=glmer, data=brod,family='binomial',
                        control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5))) 
model.list <- model.sel(forest.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brod_forest_selection.csv")

# building models
build.models <- lapply(build_formulae, FUN=glmer, data=brod, family='binomial',
                       control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5))) 
model.list <- model.sel(build.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brod_building_selection.csv")

# WUI models
wui.models <- lapply(wui_formulae, data=brod, FUN=glmer, family='binomial',
                     control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5))) 
model.list <- model.sel(wui.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brod_wui_selection.csv")




## Bromadiolone






## Diphacinone
# Age sex
as.models <- lapply(agesex_forms, FUN=glmer, data=diph, family="binomial")
model.list <- model.sel(ag.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/diph_age-sex_selection.csv")

# agriculture
ag.models <- lapply(ag_formulae, FUN=glmer, data=diph, family='binomial')
model.list <- model.sel(ag.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/diph_ag_selection.csv")

# beech
beech.models <- lapply(beech_formulae,FUN=glmer, data=diph,family='binomial') 
model.list <- model.sel(beech.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/diph_beech_selection.csv")

# forest
forest.models <- lapply(forest_formulae, FUN=glmer, data=diph,family='binomial') 
model.list <- model.sel(forest.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/diph_forest_logistic_model_selection_table.csv")

# building models
build.models <- lapply(build_formulae, FUN=glmer, data=diph, family='binomial') 
model.list <- model.sel(build.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/diph_building_logistic_model_selection_table.csv")

# WUI models
wui.models <- lapply(wui_formulae, data=diph, FUN=glmer, family='binomial') 
model.list <- model.sel(wui.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/diph_wui_logistic_model_selection_table.csv")

# Global models
full.models <- lapply(full_models, data=diph, FUN=glmer, family='binomial') 
model.list <- model.sel(full.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab


## Dicoumarol


