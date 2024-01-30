library(tidyverse)

# Individual compounds
# Landscape covariates
dat <- read_csv("output/AR_results_wide.csv")


dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")

ggplot(dat) +
  geom_point(aes(x=rand_x, y=rand_y, color=n.compounds.T))
