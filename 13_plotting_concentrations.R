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

# Convert ppm to ng/g
dat <- dat %>% mutate(concentration=concentration*1000)

avg <- dat %>% group_by(compound) %>%
  reframe(mean=mean(concentration)) %>%
  filter()
  

acute <- data.frame(FisherID=c("20F509","20F509","20F509","20M608","20M609"), 
                    id=c("Adult Female","Adult Female","Adult Female","Adult Male 1","Adult Male 2"),
                    compound=c("Bromadiolone", "Brodifacoum", "Diphacinone","Bromadiolone","Diphacinone"),
                    concentration=c(2070,1850,696,1780,1060))


conc_plot <- ggplot() +
  geom_jitter(data=dat, aes(x=concentration, y=compound, color=gen), alpha=0.3, size=2) +
  scale_color_manual(values=c("#40004b", "#00441b"), name="Compound type") +
  geom_point(data=acute, aes(x=concentration, y=compound, shape=id), color="black", size=3) +
  guides(shape=guide_legend(title="Fisher")) +
  theme_classic() +
  ylab("Compound") + xlab("Concentration (ng/g)") +
  theme(legend.position=c(1,0),
        legend.justification=c(1,0),
        legend.background = element_rect(fill=NA),
        panel.border=element_rect(fill=NA, color="black", linewidth=1),
        axis.title=element_text(size=12),
        axis.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))


avgs <- dat %>% filter(concentration>0.001) %>% 
  group_by(compound) %>%
  reframe(mean=mean(concentration))
  
  