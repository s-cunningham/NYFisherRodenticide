library(tidyverse)
library(patchwork)

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(RegionalID:ncomp) %>%
  mutate(mast_year=if_else(year==2019, 1, 2), # mast years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

annu_ncomp <- dat %>% filter(pt_index==1) %>% group_by(year, ncomp) %>% count() %>%
                  mutate(pct=if_else(year==2020, n/138, n/100))

p1 <- ggplot(annu_ncomp) +
  coord_cartesian(ylim=c(0,0.6)) +
  geom_col(aes(x=year, y=pct, group=ncomp, color=factor(ncomp), fill=factor(ncomp)), position="dodge") +
  scale_fill_manual(values=alpha(c("#762a83","#af8dc3","#e7d4e8","#d9f0d3","#7fbf7b","#1b7837"), 0.6), name="Number of compounds detected") +
  scale_color_manual(values=c("#40004b","#762a83","#9970ab","#5aae61","#1b7837","#00441b"), name="Number of compounds detected") +
  ylab("Proportion of fishers\nwith x compounds") + xlab("Year") +
  guides(fill=guide_legend(position = "inside", nrow=1), color=guide_legend(position = "inside", nrow=1)) +
  theme_classic() +
  theme(legend.position.inside=c(0.27,0.895),  #c(0.23,0.895)
        legend.background = element_rect(fill=NA),
        panel.border=element_rect(fill=NA, color="black", linewidth=1),
        axis.title.y=element_text(size=12),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))

# read data for individual compounds
diph <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Diphacinone") %>%
  select(RegionalID,bin.exp)
diph <- left_join(dat, diph, by="RegionalID")
diph <- diph %>% filter(pt_index==1) %>% select(RegionalID,year,bin.exp) %>% rename(Diphacinone=bin.exp)

brod <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Brodifacoum") %>%
  select(RegionalID,bin.exp)
brod <- left_join(dat, brod, by="RegionalID")
brod <- brod %>% filter(pt_index==1) %>% select(RegionalID,year,bin.exp) %>% rename(Brodifacoum=bin.exp)

brom <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Bromadiolone") %>%
  select(RegionalID,bin.exp)
brom <- left_join(dat, brom, by="RegionalID")
brom <- brom %>% filter(pt_index==1) %>% select(RegionalID,year,bin.exp) %>% rename(Bromadiolone=bin.exp)


dbb <- left_join(diph, brod, by=c("RegionalID", "year")) 
dbb <- left_join(dbb, brom, by=c("RegionalID", "year"))

dbb <- dbb %>% pivot_longer(3:5, names_to="compound", values_to="detect")

dbb <- dbb %>% group_by(year, compound) %>% reframe(n=sum(detect))

dbb$compound <- factor(dbb$compound, levels=c("Diphacinone", "Brodifacoum", "Bromadiolone"))
dbb <- dbb %>% mutate(pct=if_else(year==2020, n/138, n/100))

p2 <- ggplot(dbb) +
  coord_cartesian(ylim=c(0,1)) +
  geom_col(aes(x=year, y=pct, color=factor(compound), fill=factor(compound)), position="dodge") +
  scale_color_manual(values=c("#40004b","gray50","#00441b"), name="Compound") +
  scale_fill_manual(values=alpha(c("#762a83","white","#1b7837"), 0.5), name="Compound") +
  ylab("Proportion of\nfishers positive") + xlab("Year") +
  guides(fill=guide_legend(position = "inside", nrow=1), color=guide_legend(position = "inside", nrow=1)) +
  theme_classic() +
  theme(legend.position.inside=c(0.42,0.895),  # c(0.32,0.895),
        legend.background = element_rect(fill=NA),
        panel.border=element_rect(fill=NA, color="black", linewidth=1),
        axis.title=element_text(size=12),
        axis.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))


p1 / p2 + plot_annotation(tag_levels=c('a'), tag_prefix='(', tag_suffix=')') 


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


p3 <- ggplot() +
  geom_jitter(data=dat, aes(y=concentration, x=compound, color=gen), alpha=0.3, size=2) +
  scale_color_manual(values=c("#40004b", "#00441b"), name="Compound type") +
  geom_point(data=acute, aes(y=concentration, x=compound, shape=id), color="black", size=3) +
  guides(shape=guide_legend(title="Fisher mortality")) +
  theme_classic() +
  xlab("Compound") + ylab("Concentration (ng/g)") +
  theme(legend.position=c(0,1),
        legend.justification=c(0,1),
        legend.background = element_rect(fill=NA),
        panel.border=element_rect(fill=NA, color="black", linewidth=1),
        axis.title=element_text(size=12),
        axis.text.x=element_text(size=12, angle=45, hjust=0.95),
        axis.text.y=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))  #,legend.box="horizontal"


# free(p1) / free(p2) / p3 + plot_annotation(tag_levels=c('a'), tag_prefix='(', tag_suffix=')') 

free(p1 / p2) | p3 + plot_layout(widths=c(3,2)) + plot_annotation(tag_levels=c('a'), tag_prefix='(', tag_suffix=')') 

