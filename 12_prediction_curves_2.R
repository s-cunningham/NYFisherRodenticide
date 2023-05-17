library(tidyverse)
library(glmmTMB)
library(sjPlot)
library(sjmisc)
library(lme4)
library(patchwork)

options(scipen=999, digits=3)
set.seed(123)
theme_set(theme_classic())

#### Parallel processing ####
nt <- min(parallel::detectCores(),4)

## Read data
dat1 <- read_csv("output/model_data_notscaled.csv")

# scale data
dat1 <- within(dat1, mast.s <- scale(lag_beechnuts))
dat1 <- within(dat1, Age.s <- scale(Age))
dat1 <- within(dat1, intermix.s <- scale(mix_15_100))
dat1 <- within(dat1, pasture.s <- scale(pasture_15))

# Scale and center variables
dat1[,c(8,16:83)] <- scale(dat1[,c(8,16:83)])

dfp_mix <- list()
dfp_lbn <- list()
dfp_past <- list()
dfp_as <- list()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- dat1[dat1$pt_index==i,]
  
  # Run model 
  m1_pt <- glmmTMB(n.compounds.T ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * lag_beechnuts + (1|WMU), data=pt, 
                   family=compois(link = "log"), control=glmmTMBControl(parallel=nt))
  

  p_mix <- plot_model(m1_pt, type="pred", terms=c("mix_15_100 [-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3]", 
                                                  "Sex [all]"))
  dfp_mix[[i]] <- as.data.frame(p_mix$data)
  
  p_past <- plot_model(m1_pt, type="pred", terms=c("pasture_15 [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3]",
                                                   "Sex [all]" ))
  dfp_past[[i]] <- as.data.frame(p_past$data)
  
  
  p_lbn <- plot_model(m1_pt, type="pred", terms=c("lag_beechnuts [-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]", 
                                                  "Sex [all]"), pred.type="fe")
  dfp_lbn[[i]] <- as.data.frame(p_lbn$data)
  

  p_as <- plot_model(m1_pt, type="pred", terms=c("Age [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]", 
                                                 "Sex [all]"), pred.type="fe")
  dfp_as[[i]] <- as.data.frame(p_as$data)
  
}

## Average predictions
df_mix <- dfp_mix %>% plyr::ldply(data.frame) %>% as_tibble() %>%
                select(-group_col) %>% 
                group_by(group, x) %>%
                summarize(pred=mean(predicted), se=mean(std.error), ci_low=mean(conf.low), ci_high=mean(conf.high)) %>%
                ungroup() %>%
                mutate(unscWUI=x * attr(dat1$intermix.s, 'scaled:scale') + attr(dat1$intermix.s, 'scaled:center')) %>%
                rename(Sex=group)

df_lbn <- dfp_lbn %>% plyr::ldply(data.frame) %>% as_tibble() %>%
                select(-group_col) %>% 
                group_by(group, x) %>%
                summarize(pred=mean(predicted), se=mean(std.error), ci_low=mean(conf.low), ci_high=mean(conf.high)) %>%
                ungroup() %>%
                mutate(unscBeechnuts=x * attr(dat1$mast.s, 'scaled:scale') + attr(dat1$mast.s, 'scaled:center')) %>%
                rename(Sex=group)

df_past <- dfp_past %>% plyr::ldply(data.frame) %>% as_tibble() %>%
                select(-group_col) %>% 
                group_by(group, x) %>%
                summarize(pred=mean(predicted), se=mean(std.error), ci_low=mean(conf.low), ci_high=mean(conf.high)) %>%
                ungroup() %>%
                mutate(unscPasture=x * attr(dat1$pasture.s, 'scaled:scale') + attr(dat1$pasture.s, 'scaled:center')) %>%
                rename(Sex=group)

df_as <- dfp_as %>% plyr::ldply(data.frame) %>% as_tibble() %>%
                select(-group_col) %>% 
                group_by(group, x) %>%
                summarize(pred=mean(predicted), se=mean(std.error), ci_low=mean(conf.low), ci_high=mean(conf.high)) %>%
                ungroup() %>%
                mutate(unscAge=x * attr(dat1$Age.s, 'scaled:scale') + attr(dat1$Age.s, 'scaled:center')) %>%
                rename(Sex=group)

ggplot(df_mix) + 
  geom_ribbon(aes(x=unscWUI, ymin=ci_low, ymax=ci_high, group=Sex, fill=Sex), alpha=0.2) +
  geom_line(aes(x=unscWUI, y=pred, group=Sex, color=Sex), linewidth=1) 

ggplot(df_lbn) + 
  geom_ribbon(aes(x=unscBeechnuts, ymin=ci_low, ymax=ci_high, group=Sex, fill=Sex), alpha=0.2) +
  geom_line(aes(x=unscBeechnuts, y=pred, group=Sex, color=Sex), linewidth=1) 



age_plot <- ggplot(data=df_as) + 
  geom_ribbon(aes(x=unscAge, ymin=ci_low, ymax=ci_high, color=Sex, fill=Sex), alpha=0.3) +
  geom_line(aes(x=unscAge, y=pred, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), limits=c(0, 5.5)) +
  ylab("Expected number of compounds") + xlab("Age (years)") +
  theme(legend.position=c(0,1),
        legend.justification=c(0,1),
        legend.background=element_rect(fill=NA),
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11))

mix_plot <- ggplot(data=df_mix) + 
  geom_ribbon(aes(x=unscWUI, ymin=ci_low, ymax=ci_high, color=Sex, fill=Sex), alpha=0.3) +
  geom_line(aes(x=unscWUI, y=pred, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), limits=c(0, 5.5)) +
  xlab("Proportion intermix") +
  theme(legend.position="none",
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        axis.title.y=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11))

beech_plot <- ggplot(data=df_lbn) + 
  geom_ribbon(aes(x=unscBeechnuts, ymin=ci_low, ymax=ci_high, color=Sex, fill=Sex), alpha=0.3) +
  geom_line(aes(x=unscBeechnuts, y=pred, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), limits=c(0, 5.5)) +
  ylab("Expected number of compounds") + xlab("Beechnut count (1-yr lag)") +
  theme(legend.position="none",
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11))


past_plot <- ggplot(data=df_past) + 
  geom_ribbon(aes(x=unscPasture, ymin=ci_low, ymax=ci_high, color=Sex, fill=Sex), alpha=0.3) +
  geom_line(aes(x=unscPasture, y=pred, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), limits=c(0, 5.5)) +
  xlab("Proportion pasture") +
  theme(legend.position="none",
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        axis.title.y=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11))


(age_plot + mix_plot) / (beech_plot + past_plot) + plot_annotation(tag_levels="a", tag_prefix="(", tag_suffix=")")


###### Binary models #####

# Load binary for each compound, combine covariates
## subset to only the covariates needed and save file
brod <- read_csv("output/brod_full_data_unsc.csv") %>% 
            select(RegionalID:Town, WMUA, bin.exp, mix_15_100, pasture_15, lag_beechnuts, BBA_15)
brom <- read_csv("output/brom_full_data_unsc.csv") %>% 
  select(RegionalID:Town, WMUA, bin.exp, mix_15_100, pasture_15, lag_beechnuts, BBA_15)
diph <- read_csv("output/diph_full_data_unsc.csv") %>% 
  select(RegionalID:Town, WMUA, bin.exp, mix_15_100, pasture_15, lag_beechnuts, BBA_15)
dico <- read_csv("output/dico_full_data_unsc.csv") %>% 
  select(RegionalID:Town, WMUA, bin.exp, mix_15_100, pasture_15, lag_beechnuts, BBA_15)

brod <- within(brod, mast.s <- scale(lag_beechnuts))
brod <- within(brod, Age.s <- scale(Age))
brod <- within(brod, intermix.s <- scale(mix_15_100))
brod <- within(brod, pasture.s <- scale(pasture_15))

brom <- within(brom, mast.s <- scale(lag_beechnuts))
brom <- within(brom, Age.s <- scale(Age))
brom <- within(brom, intermix.s <- scale(mix_15_100))
brom <- within(brom, pasture.s <- scale(pasture_15))

diph <- within(diph, mast.s <- scale(lag_beechnuts))
diph <- within(diph, Age.s <- scale(Age))
diph <- within(diph, intermix.s <- scale(mix_15_100))
diph <- within(diph, pasture.s <- scale(pasture_15))

dico <- within(dico, Age.s <- scale(Age))
dico <- within(dico, pasture.s <- scale(pasture_15))

# Scale
brod[,c(8,17:20)] <- scale(brod[,c(8,17:20)])
brom[,c(8,17:20)] <- scale(brom[,c(8,17:20)])
diph[,c(8,17:20)] <- scale(diph[,c(8,17:20)])
dico[,c(8,17:20)] <- scale(dico[,c(8,17:20)])

## Loop over point sets

brod_mix <- brod_as <- brod_past <- brod_lbn <- list()
brom_mix <- brom_as <- brom_past <- brom_lbn <- list()
diph_mix <- diph_as <- diph_past <- diph_lbn <- list()
dico_past <- dico_as <- list()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  brod_pt <- brod[brod$pt_index==i,]
  brom_pt <- brom[brom$pt_index==i,]
  diph_pt <- diph[diph$pt_index==i,]
  dico_pt <- dico[dico$pt_index==i,]
  
  # Brodifacoum
  m1_brod <- lme4::glmer(bin.exp ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * lag_beechnuts + (1|WMU), 
                       family=binomial(link="logit"), data=brod_pt)
  
  brod_mix_p <- plot_model(m1_brod, type="pred", terms=c("mix_15_100 [-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3]", 
                                                         "Sex [all]"))
  brod_mix[[i]] <- as.data.frame(brod_mix_p$data)
  
  brod_past_p <- plot_model(m1_brod, type="pred", terms=c("pasture_15 [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3]",
                                                          "Sex [all]" ))
  brod_past[[i]] <- as.data.frame(brod_past_p$data)
  
  brod_lbn_p <- plot_model(m1_brod, type="pred", terms=c("lag_beechnuts [-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]", 
                                                         "Sex [all]"), pred.type="fe")
  brod_lbn[[i]] <- as.data.frame(brod_lbn_p$data)
  
  brod_as_p <- plot_model(m1_brod, type="pred", terms=c("Age [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]", 
                                                        "Sex [all]"), pred.type="fe")
  brod_as[[i]] <- as.data.frame(brod_as_p$data)

  
  ## Bromadiolone
  m1_brom <- lme4::glmer(bin.exp ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * lag_beechnuts + (1|WMU), 
                         family=binomial(link="logit"), data=brom_pt)
  
  brom_mix_p <- plot_model(m1_brom, type="pred", terms=c("mix_15_100 [-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3]", 
                                                         "Sex [all]"))
  brom_mix[[i]] <- as.data.frame(brom_mix_p$data)
  
  brom_past_p <- plot_model(m1_brom, type="pred", terms=c("pasture_15 [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3]",
                                                          "Sex [all]" ))
  brom_past[[i]] <- as.data.frame(brom_past_p$data)
  
  brom_lbn_p <- plot_model(m1_brom, type="pred", terms=c("lag_beechnuts [-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]", 
                                                         "Sex [all]"), pred.type="fe")
  brom_lbn[[i]] <- as.data.frame(brom_lbn_p$data)
  
  brom_as_p <- plot_model(m1_brom, type="pred", terms=c("Age [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]", 
                                                        "Sex [all]"), pred.type="fe")
  brom_as[[i]] <- as.data.frame(brom_as_p$data)
  
  
  ## Diphacinone
  m1_diph <- lme4::glmer(bin.exp ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * lag_beechnuts + (1|WMU), 
                         family=binomial(link="logit"), data=diph_pt)
  
  diph_mix_p <- plot_model(m1_diph, type="pred", terms=c("mix_15_100 [-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3]", 
                                                         "Sex [all]"))
  diph_mix[[i]] <- as.data.frame(diph_mix_p$data)
  
  diph_past_p <- plot_model(m1_diph, type="pred", terms=c("pasture_15 [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3]",
                                                          "Sex [all]" ))
  diph_past[[i]] <- as.data.frame(diph_past_p$data)
  
  diph_lbn_p <- plot_model(m1_diph, type="pred", terms=c("lag_beechnuts [-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]", 
                                                         "Sex [all]"), pred.type="fe")
  diph_lbn[[i]] <- as.data.frame(diph_lbn_p$data)
  
  diph_as_p <- plot_model(m1_diph, type="pred", terms=c("Age [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]", 
                                                        "Sex [all]"), pred.type="fe")
  diph_as[[i]] <- as.data.frame(diph_as_p$data)
  
  
  ## Dicoumarol
  m1_dico <- lme4::glmer(bin.exp ~ Sex + Age + pasture_15 + (1|Region), 
                         family=binomial(link="logit"), data=dico_pt)

  dico_past_p <- plot_model(m1_dico, type="pred", terms=c("pasture_15 [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3]",
                                                   "Sex [all]" ))
  dico_past[[i]] <- as.data.frame(dico_past_p$data)
  
  dico_as_p <- plot_model(m1_dico, type="pred", terms=c("Age [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]", 
                                                 "Sex [all]"), pred.type="fe")
  dico_as[[i]] <- as.data.frame(dico_as_p$data)
  
}


### Exposure to individual compounds by age
brod_as <- brod_as %>% plyr::ldply(data.frame) %>% as_tibble() %>%
  select(-group_col) %>% 
  group_by(group, x) %>%
  summarize(pred=mean(predicted), se=mean(std.error), ci_low=mean(conf.low), ci_high=mean(conf.high)) %>%
  ungroup() %>%
  mutate(unscX=x * attr(brod$Age.s, 'scaled:scale') + attr(brod$Age.s, 'scaled:center')) %>%
  rename(Sex=group)
brod_as$compound <- "Brodifacoum"

brom_as <- brom_as %>% plyr::ldply(data.frame) %>% as_tibble() %>%
  select(-group_col) %>% 
  group_by(group, x) %>%
  summarize(pred=mean(predicted), se=mean(std.error), ci_low=mean(conf.low), ci_high=mean(conf.high)) %>%
  ungroup() %>%
  mutate(unscX=x * attr(brom$Age.s, 'scaled:scale') + attr(brom$Age.s, 'scaled:center')) %>%
  rename(Sex=group)
brom_as$compound <- "Bromadiolone"

diph_as <- diph_as %>% plyr::ldply(data.frame) %>% as_tibble() %>%
  select(-group_col) %>% 
  group_by(group, x) %>%
  summarize(pred=mean(predicted), se=mean(std.error), ci_low=mean(conf.low), ci_high=mean(conf.high)) %>%
  ungroup() %>%
  mutate(unscX=x * attr(diph$Age.s, 'scaled:scale') + attr(diph$Age.s, 'scaled:center')) %>%
  rename(Sex=group)
diph_as$compound <- "Diphacinone"

dico_as <- dico_as %>% plyr::ldply(data.frame) %>% as_tibble() %>%
  select(-group_col) %>% 
  group_by(group, x) %>%
  summarize(pred=mean(predicted), se=mean(std.error), ci_low=mean(conf.low), ci_high=mean(conf.high)) %>%
  ungroup() %>%
  mutate(unscX=x * attr(dico$Age.s, 'scaled:scale') + attr(dico$Age.s, 'scaled:center')) %>%
  rename(Sex=group)
dico_as$compound <- "Dicoumarol"

# combine
bin_age <- bind_rows(brod_as, brom_as, diph_as, dico_as)

# plot
ggplot(bin_age) + 
  geom_ribbon(aes(x=unscX, ymin=ci_low, ymax=ci_high, color=Sex, fill=Sex), alpha=0.3) +
  geom_line(aes(x=unscX, y=pred, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83")) +
  scale_fill_manual(values=c("#1b7837", "#762a83")) +
  facet_grid(.~compound) +
  ylim(0,1) +
  ylab("Probability of exposure") + xlab("Age (years)") +
  theme(legend.position=c(0,0),
        legend.justification=c(0,0),
        # legend.position="bottom",
        legend.background=element_rect(fill=NA),
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        axis.text=element_text(size=12),
        strip.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11))




### Exposure to individual compounds by WUI intermix
brod_mix <- brod_mix %>% plyr::ldply(data.frame) %>% as_tibble() %>%
  select(-group_col) %>% 
  group_by(group, x) %>%
  summarize(pred=mean(predicted), se=mean(std.error), ci_low=mean(conf.low), ci_high=mean(conf.high)) %>%
  ungroup() %>%
  mutate(unscX=x * attr(brod$intermix.s, 'scaled:scale') + attr(brod$intermix.s, 'scaled:center')) %>%
  rename(Sex=group)
brod_mix$compound <- "Brodifacoum"

brom_mix <- brom_mix %>% plyr::ldply(data.frame) %>% as_tibble() %>%
  select(-group_col) %>% 
  group_by(group, x) %>%
  summarize(pred=mean(predicted), se=mean(std.error), ci_low=mean(conf.low), ci_high=mean(conf.high)) %>%
  ungroup() %>%
  mutate(unscX=x * attr(brom$intermix.s, 'scaled:scale') + attr(brom$intermix.s, 'scaled:center')) %>%
  rename(Sex=group)
brom_mix$compound <- "Bromadiolone"

diph_mix <- diph_mix %>% plyr::ldply(data.frame) %>% as_tibble() %>%
  select(-group_col) %>% 
  group_by(group, x) %>%
  summarize(pred=mean(predicted), se=mean(std.error), ci_low=mean(conf.low), ci_high=mean(conf.high)) %>%
  ungroup() %>%
  mutate(unscX=x * attr(diph$intermix.s, 'scaled:scale') + attr(diph$intermix.s, 'scaled:center')) %>%
  rename(Sex=group)
diph_mix$compound <- "Diphacinone"


# combine
bin_mix <- bind_rows(brod_mix, brom_mix, diph_mix)

# plot
ggplot(bin_mix) + 
  geom_ribbon(aes(x=unscX, ymin=ci_low, ymax=ci_high, color=Sex, fill=Sex), alpha=0.3) +
  geom_line(aes(x=unscX, y=pred, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83")) +
  scale_fill_manual(values=c("#1b7837", "#762a83")) +
  facet_grid(.~compound) +
  ylim(0,1) +
  ylab("Probability of exposure") + xlab("Proportion intermix") +
  theme(legend.position=c(0,0),
        legend.justification=c(0,0),
        # legend.position="bottom",
        legend.background=element_rect(fill=NA),
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        axis.text=element_text(size=12),
        strip.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11))
