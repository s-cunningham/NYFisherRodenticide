library(tidyverse)
library(boot)
library(tagger)

theme_set(theme_classic())

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-build_cat) %>%
  mutate(mast_year=if_else(year==2019, 2, 1), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# read data for individual compounds
diph <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Diphacinone") %>%
  select(RegionalID,bin.exp)
dat <- left_join(dat, diph, by="RegionalID")
dat <- dat %>% select(RegionalID:Town,bin.exp,deciduous:mast_year)

beech <- dat %>% select(beechnuts, bin.exp)
structures <- dat %>% select(nbuildings, bin.exp)
standage <- dat %>% select(stand_age_mean, bin.exp)

# Scale variables
dat$Age <- scale(dat$Age)
dat$nbuildings <- scale(dat$nbuildings)
dat$lag_beechnuts <- scale(dat$lag_beechnuts)
dat$stand_age_mean <- scale(dat$stand_age_mean)

# Load posterior samples
diph.out1 <- readRDS("output/model_output/diph.out1.rds")
diph.out1 <- do.call("rbind",diph.out1)
diph.out2 <- readRDS("output/model_output/diph.out2.rds")
diph.out2 <- do.call("rbind",diph.out2)
diph.out3 <- readRDS("output/model_output/diph.out3.rds")
diph.out3 <- do.call("rbind",diph.out3)
diph.out4 <- readRDS("output/model_output/diph.out4.rds")
diph.out4 <- do.call("rbind",diph.out4)
diph.out5 <- readRDS("output/model_output/diph.out5.rds")
diph.out5 <- do.call("rbind",diph.out5)
diph.out6 <- readRDS("output/model_output/diph.out6.rds")
diph.out6 <- do.call("rbind",diph.out6)
diph.out7 <- readRDS("output/model_output/diph.out7.rds")
diph.out7 <- do.call("rbind",diph.out7)
diph.out8 <- readRDS("output/model_output/diph.out8.rds")
diph.out8 <- do.call("rbind",diph.out8)
diph.out9 <- readRDS("output/model_output/diph.out9.rds")
diph.out9 <- do.call("rbind",diph.out9)
diph.out10 <- readRDS("output/model_output/diph.out10.rds")
diph.out10 <- do.call("rbind",diph.out10)

## Combine all iterations
beta_age <- c(diph.out1[,19],diph.out2[,19],diph.out3[,19],diph.out4[,19],diph.out5[,19],
              diph.out6[,19],diph.out7[,19],diph.out8[,19],diph.out9[,19],diph.out10[,19])

beta_age2 <- c(diph.out1[,20],diph.out2[,20],diph.out3[,20],diph.out4[,20],diph.out5[,20],
               diph.out6[,20],diph.out7[,20],diph.out8[,20],diph.out9[,20],diph.out10[,20])

beta_build <- c(diph.out1[,21],diph.out2[,21],diph.out3[,21],diph.out4[,21],diph.out5[,21],
                diph.out6[,21],diph.out7[,21],diph.out8[,21],diph.out9[,21],diph.out10[,21])

beta_mast <- c(diph.out1[,22],diph.out2[,22],diph.out3[,22],diph.out4[,22],diph.out5[,22],
               diph.out6[,22],diph.out7[,22],diph.out8[,22],diph.out9[,22],diph.out10[,22])

beta_sex2 <- c(diph.out1[,24],diph.out2[,24],diph.out3[,24],diph.out4[,24],diph.out5[,24],
               diph.out6[,24],diph.out7[,24],diph.out8[,24],diph.out9[,24],diph.out10[,24])

beta_stand <- c(diph.out1[,25],diph.out2[,25],diph.out3[,25],diph.out4[,25],diph.out5[,25],
                diph.out6[,25],diph.out7[,25],diph.out8[,25],diph.out9[,25],diph.out10[,25])

##############

## Look at intercept
alpha <- c(rowMeans(diph.out1[,1:18]), rowMeans(diph.out2[,1:18]), rowMeans(diph.out3[,1:18]), rowMeans(diph.out4[,1:18]), rowMeans(diph.out5[,1:18]),
           rowMeans(diph.out6[,1:18]), rowMeans(diph.out7[,1:18]), rowMeans(diph.out8[,1:18]), rowMeans(diph.out9[,1:18]), rowMeans(diph.out10[,1:18]))
plot(density(alpha))

## Predicting diphifacoum exposure by age (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
age_pred <- seq(min(dat$Age),max(dat$Age),length.out=pred_length)
age_pred2 <- age_pred^2

# Average beechnut counts and number of buildlings
mean_build <- mean(dat$nbuildings)
mean_mast <- mean(dat$lag_beechnuts)
mean_stand <- mean(dat$stand_age_mean)

# Predict
age.diphM <- matrix(, nmcmc, pred_length)
age.diphF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  age.diphF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_build*mean_build + beta_mast*mean_mast + beta_stand*mean_stand) # males
  age.diphM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_build*mean_build + beta_mast*mean_mast + beta_stand*mean_stand) # females
}

# Calculate quantiles
age.diphF.qt <- apply(age.diphF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
age.diphM.qt <- apply(age.diphM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

age.qt <- bind_rows(age.diphF.qt, age.diphM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
age_pred <- age_pred * attr(dat$Age, 'scaled:scale') + attr(dat$Age, 'scaled:center')

age.qt.diph <- age.qt %>% mutate(Age=rep(age_pred, 2))
age.qt.diph$compound <- "Diphacinone"

ggplot(age.qt.diph) +
  coord_cartesian(ylim=c(0, 1), xlim=c(0,8.5)) +
  geom_ribbon(aes(x=Age, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Age, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_x_continuous(breaks=seq(0,9)) +
  ylab("Probability of exposure") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position = "none", 
        legend.background = element_rect(fill=NA))



## Predicting diphifacoum exposure by mast cycles (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
mast_pred <- seq(min(dat$lag_beechnuts),max(dat$lag_beechnuts),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_build <- mean(dat$nbuildings)
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_stand <- mean(dat$stand_age_mean)

# Predict
mast.diphM <- matrix(, nmcmc, pred_length)
mast.diphF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  mast.diphF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mast_pred[j] + beta_stand*mean_stand) # males
  mast.diphM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mast_pred[j] + beta_stand*mean_stand) # females
}

# Calculate quantiles
mast.diphF.qt <- apply(mast.diphF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
mast.diphM.qt <- apply(mast.diphM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

mast.qt <- bind_rows(mast.diphF.qt, mast.diphM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
mast_pred <- mast_pred * attr(dat$lag_beechnuts, 'scaled:scale') + attr(dat$lag_beechnuts, 'scaled:center')

mast.qt.diph  <- mast.qt %>% mutate(Beechnuts=rep(mast_pred, 2))
mast.qt.diph$compound <- "Diphacinone"

ggplot(mast.qt.diph) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=Beechnuts, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Beechnuts, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  # scale_x_continuous(breaks=seq(0,9)) +
  ylab("Probability of exposure") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position="none")



## Predicting diphifacoum exposure by stand age(and sex)
nmcmc <- length(beta_stand)
pred_length <- 100
stand_pred <- seq(min(dat$stand_age_mean),max(dat$stand_age_mean),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_build <- mean(dat$nbuildings)
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_mast <- mean(dat$lag_beechnuts)

# Predict
stand.diphM <- matrix(, nmcmc, pred_length)
stand.diphF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  stand.diphF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mean_mast + beta_stand*stand_pred[j]) # males
  stand.diphM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mean_mast + beta_stand*stand_pred[j]) # females
}

# Calculate quantiles
stand.diphF.qt <- apply(stand.diphF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
stand.diphM.qt <- apply(stand.diphM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

stand.qt <- bind_rows(stand.diphF.qt, stand.diphM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
stand_pred <- stand_pred * attr(dat$stand_age_mean, 'scaled:scale') + attr(dat$stand_age_mean, 'scaled:center')

stand.qt.diph  <- stand.qt %>% mutate(StandAge=rep(stand_pred, 2))
stand.qt.diph$compound <- "Diphacinone"

ggplot(stand.qt.diph) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=StandAge, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=StandAge, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_x_continuous(breaks=seq(0,90,10)) +
  ylab("Probability of exposure") + xlab("Stand age (years)") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position="none")


## Predicting diphifacoum exposure by stand age(and sex)
nmcmc <- length(beta_stand)
pred_length <- 100
build_pred <- seq(min(dat$nbuildings),max(dat$nbuildings),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_mast <- mean(dat$lag_beechnuts)
mean_stand <- mean(dat$stand_age_mean)

# Predict
build.diphM <- matrix(, nmcmc, pred_length)
build.diphF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  build.diphF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*build_pred[j] + beta_mast*mean_mast + beta_stand*mean_stand) # males
  build.diphM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*build_pred[j] + beta_mast*mean_mast + beta_stand*mean_stand) # females
}

# Calculate quantiles
build.diphF.qt <- apply(build.diphF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
build.diphM.qt <- apply(build.diphM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

build.qt <- bind_rows(build.diphF.qt, build.diphM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
build_pred <- build_pred * attr(dat$nbuildings, 'scaled:scale') + attr(dat$nbuildings, 'scaled:center')

build.qt.diph <- build.qt %>% mutate(Buildings=rep(build_pred, 2))
build.qt.diph$compound <- "Diphacinone"

ggplot(build.qt.diph) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=Buildings, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Buildings, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  # scale_x_continuous(breaks=seq(0,9)) +
  ylab("Probability of exposure") + 
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position="none")

## Organize for full plot
age.qt.diph <- age.qt.diph %>% mutate(x="Age") %>% rename(x_val=Age) %>% select(compound, Sex, x, x_val, median:uci)
mast.qt.diph <- mast.qt.diph %>% mutate(x="Beechnuts") %>% rename(x_val=Beechnuts) %>% select(compound, Sex, x, x_val, median:uci)
stand.qt.diph <- stand.qt.diph %>% mutate(x="StandAge") %>% rename(x_val=StandAge) %>% select(compound, Sex, x, x_val, median:uci)
build.qt.diph <- build.qt.diph %>% mutate(x="Buildings") %>% rename(x_val=Buildings) %>% select(compound, Sex, x, x_val, median:uci)

qt.diph <- bind_rows(age.qt.diph, mast.qt.diph, stand.qt.diph, build.qt.diph)

ggplot(qt.diph) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=x_val, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=x_val, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  facet_wrap(vars(x), scales="free_x")


############# Plotting all three compounds

all.qt <- bind_rows(qt.diph, qt.brom, qt.brod)
all.qt$compound <- factor(all.qt$compound, levels=c("Diphacinone", "Brodifacoum", "Bromadiolone"))

all.qt <- all.qt %>% mutate(x=case_when(x=="Age" ~ "Age (years)",
                                        x=="Beechnuts" ~ "Beechnut count",
                                        x=="Buildings" ~ "Building count",
                                        x=="StandAge" ~ "Stand age (years)"))


ggplot(all.qt) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=x_val, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=x_val, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  facet_grid(compound~x, scales="free_x", switch="both", axes = "all", axis.labels = "margins") +
  ylab("Probability of exposure") +
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
# Save 11.4 x 7.38


