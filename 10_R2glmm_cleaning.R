library(tidyverse)


brod_r2cm <- readRDS("results/brod_R2glmm.rds")
brom_r2cm <- readRDS("results/brom_R2glmm.rds")
diph_r2cm <- readRDS("results/diph_R2glmm.rds")
dico_r2cm <- readRDS("results/dico_R2glmm.rds")

## Conditional R2
brod_cond <- bind_rows(brod_r2cm[[1]]$R2_conditional,
                       brod_r2cm[[2]]$R2_conditional,
                       brod_r2cm[[3]]$R2_conditional,
                       brod_r2cm[[4]]$R2_conditional,
                       brod_r2cm[[5]]$R2_conditional,
                       brod_r2cm[[6]]$R2_conditional,
                       brod_r2cm[[7]]$R2_conditional,
                       brod_r2cm[[8]]$R2_conditional,
                       brod_r2cm[[9]]$R2_conditional,
                       brod_r2cm[[10]]$R2_conditional) %>% 
              rename(ConditionalR2=`Conditional R2`) %>%
              summarize(mean_r2=mean(ConditionalR2), CIL=mean(CI_low), CIH=mean(CI_high))

brod_cond$compound <- "brodifacoum"
brod_cond$R2glmm <- "conditional"

brom_cond <- bind_rows(brom_r2cm[[1]]$R2_conditional,
                       brom_r2cm[[2]]$R2_conditional,
                       brom_r2cm[[3]]$R2_conditional,
                       brom_r2cm[[4]]$R2_conditional,
                       brom_r2cm[[5]]$R2_conditional,
                       brom_r2cm[[6]]$R2_conditional,
                       brom_r2cm[[7]]$R2_conditional,
                       brom_r2cm[[8]]$R2_conditional,
                       brom_r2cm[[9]]$R2_conditional,
                       brom_r2cm[[10]]$R2_conditional) %>% 
  rename(ConditionalR2=`Conditional R2`) %>%
  summarize(mean_r2=mean(ConditionalR2), CIL=mean(CI_low), CIH=mean(CI_high))

brom_cond$compound <- "bromadiolone"
brom_cond$R2glmm <- "conditional"

diph_cond <- bind_rows(diph_r2cm[[1]]$R2_conditional,
                       diph_r2cm[[2]]$R2_conditional,
                       diph_r2cm[[3]]$R2_conditional,
                       diph_r2cm[[4]]$R2_conditional,
                       diph_r2cm[[5]]$R2_conditional,
                       diph_r2cm[[6]]$R2_conditional,
                       diph_r2cm[[7]]$R2_conditional,
                       diph_r2cm[[8]]$R2_conditional,
                       diph_r2cm[[9]]$R2_conditional,
                       diph_r2cm[[10]]$R2_conditional) %>% 
  rename(ConditionalR2=`Conditional R2`) %>%
  summarize(mean_r2=mean(ConditionalR2), CIL=mean(CI_low), CIH=mean(CI_high))            

diph_cond$compound <- "diphacinone"
diph_cond$R2glmm <- "conditional"

dico_cond <- bind_rows(dico_r2cm[[1]]$R2_conditional,
                       dico_r2cm[[2]]$R2_conditional,
                       dico_r2cm[[3]]$R2_conditional,
                       dico_r2cm[[4]]$R2_conditional,
                       dico_r2cm[[5]]$R2_conditional,
                       dico_r2cm[[6]]$R2_conditional,
                       dico_r2cm[[7]]$R2_conditional,
                       dico_r2cm[[8]]$R2_conditional,
                       dico_r2cm[[9]]$R2_conditional,
                       dico_r2cm[[10]]$R2_conditional) %>% 
  rename(ConditionalR2=`Conditional R2`) %>%
  summarize(mean_r2=mean(ConditionalR2), CIL=mean(CI_low), CIH=mean(CI_high)) 

dico_cond$compound <- "dicoumarol"
dico_cond$R2glmm <- "conditional"

condR2 <- bind_rows(brom_cond, brod_cond, diph_cond, dico_cond)


## Marginal R2
brod_mar <- bind_rows(brod_r2cm[[1]]$R2_marginal,
                      brod_r2cm[[2]]$R2_marginal,
                      brod_r2cm[[3]]$R2_marginal,
                      brod_r2cm[[4]]$R2_marginal,
                      brod_r2cm[[5]]$R2_marginal,
                      brod_r2cm[[6]]$R2_marginal,
                      brod_r2cm[[7]]$R2_marginal,
                      brod_r2cm[[8]]$R2_marginal,
                      brod_r2cm[[9]]$R2_marginal,
                      brod_r2cm[[10]]$R2_marginal) %>% 
  rename(MarginalR2=`Marginal R2`) %>%
  summarize(mean_r2=mean(MarginalR2), CIL=mean(CI_low), CIH=mean(CI_high))

brod_mar$compound <- "brodifacoum"
brod_mar$R2glmm <- "marginal"

brom_mar <- bind_rows(brom_r2cm[[1]]$R2_marginal,
                      brom_r2cm[[2]]$R2_marginal,
                      brom_r2cm[[3]]$R2_marginal,
                      brom_r2cm[[4]]$R2_marginal,
                      brom_r2cm[[5]]$R2_marginal,
                      brom_r2cm[[6]]$R2_marginal,
                      brom_r2cm[[7]]$R2_marginal,
                      brom_r2cm[[8]]$R2_marginal,
                      brom_r2cm[[9]]$R2_marginal,
                      brom_r2cm[[10]]$R2_marginal) %>% 
  rename(MarginalR2=`Marginal R2`) %>%
  summarize(mean_r2=mean(MarginalR2), CIL=mean(CI_low), CIH=mean(CI_high))

brom_mar$compound <- "bromadiolone"
brom_mar$R2glmm <- "marginal"

diph_mar <- bind_rows(diph_r2cm[[1]]$R2_marginal,
                      diph_r2cm[[2]]$R2_marginal,
                      diph_r2cm[[3]]$R2_marginal,
                      diph_r2cm[[4]]$R2_marginal,
                      diph_r2cm[[5]]$R2_marginal,
                      diph_r2cm[[6]]$R2_marginal,
                      diph_r2cm[[7]]$R2_marginal,
                      diph_r2cm[[8]]$R2_marginal,
                      diph_r2cm[[9]]$R2_marginal,
                      diph_r2cm[[10]]$R2_marginal) %>% 
  rename(MarginalR2=`Marginal R2`) %>%
  summarize(mean_r2=mean(MarginalR2), CIL=mean(CI_low), CIH=mean(CI_high))            

diph_mar$compound <- "diphacinone"
diph_mar$R2glmm <- "marginal"

dico_mar <- bind_rows(dico_r2cm[[1]]$R2_marginal,
                      dico_r2cm[[2]]$R2_marginal,
                      dico_r2cm[[3]]$R2_marginal,
                      dico_r2cm[[4]]$R2_marginal,
                      dico_r2cm[[5]]$R2_marginal,
                      dico_r2cm[[6]]$R2_marginal,
                      dico_r2cm[[7]]$R2_marginal,
                      dico_r2cm[[8]]$R2_marginal,
                      dico_r2cm[[9]]$R2_marginal,
                      dico_r2cm[[10]]$R2_marginal) %>% 
  rename(MarginalR2=`Marginal R2`) %>%
  summarize(mean_r2=mean(MarginalR2), CIL=mean(CI_low), CIH=mean(CI_high))     

dico_mar$compound <- "dicoumarol"
dico_mar$R2glmm <- "marginal"

marR2 <- bind_rows(brom_mar, brod_mar, diph_mar, dico_mar)

### combine and save
glmmR2 <- bind_rows(marR2, condR2)

write_csv(glmmR2, "results/binary_ncomp_R2glmm.csv")
