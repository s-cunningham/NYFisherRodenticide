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
  coord_cartesian(ylim=c(0,0.5)) +
  geom_col(aes(x=year, y=pct, group=ncomp, color=factor(ncomp), fill=factor(ncomp)), position="dodge") +
  scale_color_manual(values=c("#762a83","#af8dc3","#e7d4e8","#d9f0d3","#7fbf7b","#1b7837"), name="Number of compounds detected") +
  scale_fill_manual(values=c("#762a83","#af8dc3","#e7d4e8","#d9f0d3","#7fbf7b","#1b7837"), name="Number of compounds detected") +
  ylab("Proportion of fishers with x compounds") + xlab("Year") +
  guides(fill=guide_legend(position = "inside", nrow=1), color=guide_legend(position = "inside", nrow=1)) +
  theme_classic() +
  theme(legend.position.inside=c(0.23,0.895),
        panel.border=element_rect(fill=NA, color="black", linewidth=1),
        axis.title=element_text(size=14),
        axis.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))

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
  scale_color_manual(values=c("#762a83","gray80","#1b7837"), name="Compound") +
  scale_fill_manual(values=c("#762a83","gray80","#1b7837"), name="Compound") +
  ylab("Proportion of fishers positive") + xlab("Year") +
  guides(fill=guide_legend(position = "inside", nrow=1), color=guide_legend(position = "inside", nrow=1)) +
  theme_classic() +
  theme(legend.position.inside=c(0.32,0.895),
        panel.border=element_rect(fill=NA, color="black", linewidth=1),
        axis.title=element_text(size=14),
        axis.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))


p1 / p2 + plot_annotation(tag_levels=c('a'), tag_prefix='(', tag_suffix=')') + plot_layout(axes='collect')


