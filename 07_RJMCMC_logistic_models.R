library(tidyverse)
library(nimble)
library(MCMCvis)

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-edge_density_15, -edge_density_30, -edge_density_45, -build_cat_15, -build_cat_30,-build_cat_45) %>%
  mutate(mast_year=if_else(year==2019, 2, 1), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat[,c(8,16:42)] <- scale(dat[,c(8,16:42)])

# read data for individual compounds
diph <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Diphacinone") %>%
            select(RegionalID,bin.exp)
diph <- left_join(dat, diph, by="RegionalID")
diph <- diph %>% select(RegionalID:Town,bin.exp,deciduous_15:mast_year)

brom <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Bromadiolone") %>%
  select(RegionalID,bin.exp)
brom <- left_join(dat, brom, by="RegionalID")
brom <- brom %>% select(RegionalID:Town,bin.exp,deciduous_15:mast_year)

brod <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Brodifacoum") %>%
  select(RegionalID,bin.exp)
brod <- left_join(dat, brod, by="RegionalID")
brod <- brod %>% select(RegionalID:Town,bin.exp,deciduous_15:mast_year)




