library(tidyverse)


dat <- read_csv("output/AR_results_wide.csv") %>%
        # Fill NAs with 0
        replace(is.na(.), 0) %>%
        # Pivot longer
        pivot_longer(11:21, names_to="compound", values_to="concentration") %>%
        # Remove compounds that were not detected
        filter(compound!="Pindone" & compound!="Coumafuryl" & compound!="Coumachlor") %>%
        # subset to only those above 0
        filter(concentration > 0) %>%
        # add FGAR & SGAR
        mutate(gen=if_else(compound=="Brodifacoum" |
                             compound=="Bromadiolone" |
                             compound=="Difethialone" |
                             compound=="Difenacoum", "SGAR", "FGAR"))

comps <- c("Diphacinone", "Brodifacoum", "Bromadiolone", "Dicoumarol", 
           "Difethialone", "Chlorophacinone", "Difenacoum", "Warfarin")
comps <- rev(comps)

dat$compound <- factor(dat$compound, levels=comps)


ggplot(dat) +
  geom_jitter(aes(x=concentration, y=compound, color=gen), alpha=0.3, size=2) +
  scale_color_manual(values=c("#40004b", "#00441b"), name="Compound type") +
  theme_classic() +
  ylab("Compound") + xlab("Concentration (ppm)") +
  theme(legend.position=c(1,0),
        legend.justification=c(1,0),
        legend.background = element_rect(fill=NA),
        panel.border=element_rect(fill=NA, color="black", linewidth=1),
        axis.title=element_text(size=12),
        axis.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))
