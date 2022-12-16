## Logistic models
library(lme4)



source("00_logistic_models.R")

#### Run models with glmmTMB ####
ag.models <- lapply(ag_formulae, FUN=glmer, data=brom, family='binomial')
model.list <- model.sel(ag.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brom_ag_logistic_model_selection_table.csv")

beech.models <- lapply(beech_formulae,FUN=glmer, data=brom,family='binomial') 
model.list <- model.sel(beech.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brom_beech_logistic_model_selection_table.csv")

forest.models <- lapply(forest_formulae, FUN=glmer, data=brom,family='binomial') 
model.list <- model.sel(forest.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brom_forest_logistic_model_selection_table.csv")

build.models <- lapply(build_formulae, FUN=glmer, data=brom, family='binomial') 
model.list <- model.sel(build.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brom_building_logistic_model_selection_table.csv")

lsm.models <- lapply(lsm_formulae, FUN=glmer, data=brom, family='binomial') 
model.list <- model.sel(lsm.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brom_lsm_logistic_model_selection_table.csv")

wui.models <- lapply(wui_formulae, data=brom, FUN=glmer, family='binomial') 
model.list <- model.sel(wui.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab
write_csv(model_tab, "output/model_selection/brom_wui_logistic_model_selection_table.csv")


# Full models
full.models <- lapply(full_models, data=brom, FUN=glmer, family='binomial') 
model.list <- model.sel(full.models)
model_tab <- as.data.frame(model.list)
model_tab <- model_tab %>% select(df:weight) %>% rownames_to_column(var="model") %>% as_tibble()
model_tab



ggplot(brod, aes(x=rand_x, y=rand_y, color=factor(bin.exp))) + geom_point() + theme_bw() + theme(legend.position='bottom') 
