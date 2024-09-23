library(tidyverse)
library(tagger)

theme_set(theme_classic())

load(file="data/diphacinone_preds.Rdata")
load(file="data/brodifacoum_preds.Rdata")
load(file="data/bromadiolone_preds.Rdata")
# 
# 
# intx.all <- bind_rows(diph.data$intx.qt.diph, brom.data$intx.qt.brom, brod.data$intx.qt.brod)
# 
# ggplot(intx.all) +
#   coord_cartesian(ylim=c(0, 1)) +
#   geom_ribbon(aes(x=median, ymin=lci, ymax=uci, color=Sex, fill=Sex, linetype=Mast), alpha=.1) +
#   geom_line(aes(x=median, y=median, color=Sex,  linetype=Mast), linewidth=1) +
#   scale_color_manual(values=c("#1b7837","#762a83"), name="Sex") +
#   scale_fill_manual(values=c("#1b7837","#762a83"), name="Sex") +
#   facet_grid(compound~Sex) + guides(colour="none", fill="none") +
#   ylab("Probability of exposure") + 
#   theme(panel.border=element_rect(fill=NA, color="black"),
#         legend.position=c(0,0),
#         legend.justification = c(0,0), 
#         legend.background = element_rect(fill=NA))

qt.diph <- bind_rows(diph.data$age.qt.diph, diph.data$mast.qt.diph, diph.data$wui.qt.diph)
qt.brom <- bind_rows(brom.data$age.qt.brom, brom.data$mast.qt.brom, brom.data$wui.qt.brom)
qt.brod <- bind_rows(brod.data$age.qt.brod, brod.data$mast.qt.brod, brod.data$wui.qt.brod)

all.qt <- bind_rows(qt.diph, qt.brom, qt.brod)
all.qt$compound <- factor(all.qt$compound, levels=c("Diphacinone", "Brodifacoum", "Bromadiolone"))

all.qt <- all.qt %>% mutate(x=case_when(x=="Age" ~ "Age (years)",
                                        x=="Beechnuts" ~ "Lagged beechnut Count",
                                        x=="Intermix" ~ "% Wildland-urban intermix"))
all.qt$x <- factor(all.qt$x, levels=c("Age (years)", "Lagged beechnut Count", "% Wildland-urban intermix"))

ggplot(all.qt) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=x_val, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=x_val, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  facet_grid(compound~x, scales="free_x", switch="both", axes = "all", axis.labels = "margins") +
  ylab("Probability of exposure to:") +
  theme(legend.position="bottom",
        panel.border=element_rect(fill=NA, color="black"),
        strip.placement = "outside",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text=element_text(size=11),
        strip.text.x=element_text(size=12),
        strip.text.y=element_text(size=13),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11),
        strip.background = element_rect(color=NA, fill=NA),
        axis.ticks.length=unit(-0.1, "cm")) + 
  tag_facets(tag_prefix="    (")


# probably need to remove the left and top borders

ggsave("figs/prob_exp_marginal.svg")
# Saving 10.5 x 7.98 in image


