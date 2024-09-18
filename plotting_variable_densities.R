library(tidyverse)
library(tagger)

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  select(pt_index, evergreen, edge_density, stand_age_mean, stand_age_sd, cc_mean, nbuildings, intermix, interface, totalWUI,
         BMI, laggedBMI) %>%
  pivot_longer(2:12, names_to="variable", values_to="value")
  
dat$variable <- factor(dat$variable, levels=c("nbuildings", "intermix", "interface", "totalWUI",
                                              "BMI", "laggedBMI", "evergreen", "edge_density", 
                                              "stand_age_mean", "stand_age_sd", "cc_mean"), 
                       labels=c("Building count", "% Intermix", "% Interface", "% WUI", "Basal Area x Mast",
                                "Basal Area x Lagged Mast", "% Evergreen", "Forest Edge Density",
                                "Mean Stand Age", "St. Dev. Stand Age", "Mean Canopy Cover"))


ggplot(dat) +
  geom_density(aes(x=value, group=pt_index, color=pt_index, fill=pt_index), alpha=0.1) +
  scale_color_gradient(high="white", low="#5ec962") +
  scale_fill_gradient(high="white", low="#5ec962") +
  facet_wrap(vars(variable), ncol=3, scales="free") +
  ylab("Density") +
  theme_classic() +
  tag_facets(tag_prefix="(", position = "tr",) +
  theme(panel.border=element_rect(fill=NA, color="black"),
        panel.background = element_rect(fill="gray20"),
        legend.position = "none",
        strip.background = element_rect(color=NA),
        axis.title.x=element_blank(),
        strip.text=element_text(size=11),
        tagger.panel.tag.text = element_text(color = "white"),
        tagger.panel.tag.background = element_rect(fill = NA)) 

ggsave("figs/pt_index_densities.svg")
# Saving 7.86 x 9.06 in image