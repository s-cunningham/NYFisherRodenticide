## Logistic models
library(lme4)

dat1 <- dat1 %>% mutate(bin_exp = if_else(n.compounds.T > 0, 1, 0)) %>%
          select(RegionalID:n.compounds.T,bin_exp,pasture_15:build_cat_60)


source("00_logistic_models.R")

#### Run models with glmmTMB ####
ag.models <- lapply(ag_formulae, FUN=glmer, data=diph,family='binomial')
model.list <- model.sel(ag.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/ag_logistic_model_selection_table.csv")

beech.models <- lapply(beech_formulae,FUN=glmer, data=diph,family='binomial') 
model.list <- model.sel(beech.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/beech_logistic_model_selection_table.csv")

forest.models <- lapply(forest_formulae, FUN=glmer, data=diph,family='binomial') 
model.list <- model.sel(forest.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/forest_logistic_model_selection_table.csv")

build.models <- lapply(build_formulae, FUN=glmer, data=diph, family='binomial') 
model.list <- model.sel(build.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/building_logistic_model_selection_table.csv")

lsm.models <- lapply(lsm_formulae, FUN=glmmTMB, data=diph,FUN=glmer, data=dat1,family='binomial') 
model.list <- model.sel(lsm.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/lsm_logistic_model_selection_table.csv")

wui.models <- lapply(wui_formulae2, FUN=glmmTMB, data=diph, FUN=glmer, data=dat1,family='binomial') 
model.list <- model.sel(wui.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/wui_logistic_model_selection_table.csv")

