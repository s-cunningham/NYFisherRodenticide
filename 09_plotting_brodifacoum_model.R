library(tidyverse)
library(boot)

theme_set(theme_classic())

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-build_cat) %>%
  mutate(mast_year=if_else(year==2019, 2, 1), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# read data for individual compounds
brod <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Brodifacoum") %>%
  select(RegionalID,bin.exp)
dat <- left_join(dat, brod, by="RegionalID")
dat <- dat %>% select(RegionalID:Town,bin.exp,deciduous:mast_year)

beech <- dat %>% select(beechnuts, bin.exp)
structures <- dat %>% select(nbuildings, bin.exp)
standage <- dat %>% select(stand_age_mean, bin.exp)

# Scale variables
dat$Age <- scale(dat$Age)
dat$nbuildings <- scale(dat$nbuildings)
dat$beechnuts <- scale(dat$beechnuts)
dat$stand_age_mean <- scale(dat$stand_age_mean)

# Load posterior samples
brod.out1 <- readRDS("output/model_output/brod.out1.rds")
brod.out1 <- do.call("rbind",brod.out1)
brod.out2 <- readRDS("output/model_output/brod.out2.rds")
brod.out2 <- do.call("rbind",brod.out2)
brod.out3 <- readRDS("output/model_output/brod.out3.rds")
brod.out3 <- do.call("rbind",brod.out3)
brod.out4 <- readRDS("output/model_output/brod.out4.rds")
brod.out4 <- do.call("rbind",brod.out4)
brod.out5 <- readRDS("output/model_output/brod.out5.rds")
brod.out5 <- do.call("rbind",brod.out5)
brod.out6 <- readRDS("output/model_output/brod.out6.rds")
brod.out6 <- do.call("rbind",brod.out6)
brod.out7 <- readRDS("output/model_output/brod.out7.rds")
brod.out7 <- do.call("rbind",brod.out7)
brod.out8 <- readRDS("output/model_output/brod.out8.rds")
brod.out8 <- do.call("rbind",brod.out8)
brod.out9 <- readRDS("output/model_output/brod.out9.rds")
brod.out9 <- do.call("rbind",brod.out9)
brod.out10 <- readRDS("output/model_output/brod.out10.rds")
brod.out10 <- do.call("rbind",brod.out10)

## Combine all iterations
beta_age <- c(brod.out1[,19],brod.out2[,19],brod.out3[,19],brod.out4[,19],brod.out5[,19],
              brod.out6[,19],brod.out7[,19],brod.out8[,19],brod.out9[,19],brod.out10[,19])

# Calculate HDI and quantiles
quantile(beta_age, probs=c(0.025,0.5,0.975))


beta_age2 <- c(brod.out1[,20],brod.out2[,20],brod.out3[,20],brod.out4[,20],brod.out5[,20],
               brod.out6[,20],brod.out7[,20],brod.out8[,20],brod.out9[,20],brod.out10[,20])

# Calculate HDI and quantiles
quantile(beta_age2, probs=c(0.025,0.5,0.975))

beta_build <- c(brod.out1[,21],brod.out2[,21],brod.out3[,21],brod.out4[,21],brod.out5[,21],
                brod.out6[,21],brod.out7[,21],brod.out8[,21],brod.out9[,21],brod.out10[,21])

# Calculate quantiles
quantile(beta_build, probs=c(0.025,0.5,0.975))

beta_mast <- c(brod.out1[,22],brod.out2[,22],brod.out3[,22],brod.out4[,22],brod.out5[,22],
               brod.out6[,22],brod.out7[,22],brod.out8[,22],brod.out9[,22],brod.out10[,22])

# Calculate quantiles
quantile(beta_mast, probs=c(0.025,0.5,0.975))

beta_sex2 <- c(brod.out1[,24],brod.out2[,24],brod.out3[,24],brod.out4[,24],brod.out5[,24],
               brod.out6[,24],brod.out7[,24],brod.out8[,24],brod.out9[,24],brod.out10[,24])

# Calculate quantiles
quantile(beta_sex2, probs=c(0.025,0.5,0.975))


beta_stand <- c(brod.out1[,25],brod.out2[,25],brod.out3[,25],brod.out4[,25],brod.out5[,25],
               brod.out6[,25],brod.out7[,25],brod.out8[,25],brod.out9[,25],brod.out10[,25])

# Calculate quantiles
quantile(beta_stand, probs=c(0.025,0.5,0.975))


##############

## Look at intercept
alpha <- c(rowMeans(brod.out1[,1:18]), rowMeans(brod.out2[,1:18]), rowMeans(brod.out3[,1:18]), rowMeans(brod.out4[,1:18]), rowMeans(brod.out5[,1:18]),
           rowMeans(brod.out6[,1:18]), rowMeans(brod.out7[,1:18]), rowMeans(brod.out8[,1:18]), rowMeans(brod.out9[,1:18]), rowMeans(brod.out10[,1:18]))
plot(density(alpha))

## Predicting brodifacoum exposure by age (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
age_pred <- seq(min(dat$Age),max(dat$Age),length.out=pred_length)
age_pred2 <- age_pred^2

# Average beechnut counts and number of buildlings
mean_build <- mean(dat$nbuildings)
mean_mast <- mean(dat$beechnuts)
mean_stand <- mean(dat$stand_age_mean)

# Predict
age.brodM <- matrix(, nmcmc, pred_length)
age.brodF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  age.brodF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_build*mean_build + beta_mast*mean_mast + beta_stand*mean_stand) # males
  age.brodM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_build*mean_build + beta_mast*mean_mast + beta_stand*mean_stand) # females
}

# Calculate quantiles
age.brodF.qt <- apply(age.brodF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
age.brodM.qt <- apply(age.brodM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

age.qt <- bind_rows(age.brodF.qt, age.brodM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
age_pred <- age_pred * attr(dat$Age, 'scaled:scale') + attr(dat$Age, 'scaled:center')

age.qt.brod <- age.qt %>% mutate(Age=rep(age_pred, 2))
age.qt.brod$compound <- "Brodifacoum"


ggplot(age.qt.brod) +
  coord_cartesian(ylim=c(0, 1), xlim=c(0,8.5)) +
  geom_ribbon(aes(x=Age, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Age, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_x_continuous(breaks=seq(0,9)) +
  ylab("Probability of exposure") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position = c(1,1),
        legend.justification=c(1,1), 
        legend.background = element_rect(fill=NA))



## Predicting brodifacoum exposure by mast cycles (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
mast_pred <- seq(min(dat$beechnuts),max(dat$beechnuts),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_build <- mean(dat$nbuildings)
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_stand <- mean(dat$stand_age_mean)

# Predict
mast.brodM <- matrix(, nmcmc, pred_length)
mast.brodF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  mast.brodF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mast_pred[j] + beta_stand*mean_stand) # males
  mast.brodM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mast_pred[j] + beta_stand*mean_stand) # females
}

# Calculate quantiles
mast.brodF.qt <- apply(mast.brodF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
mast.brodM.qt <- apply(mast.brodM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

mast.qt <- bind_rows(mast.brodF.qt, mast.brodM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
mast_pred <- mast_pred * attr(dat$beechnuts, 'scaled:scale') + attr(dat$beechnuts, 'scaled:center')

mast.qt.brod <- mast.qt %>% mutate(Beechnuts=rep(mast_pred, 2))
mast.qt.brod$compound <- "Brodifacoum"

ggplot(mast.qt.brod) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=Beechnuts, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Beechnuts, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  # scale_x_continuous(breaks=seq(0,9)) +
  ylab("Probability of exposure") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position="none")



## Predicting brodifacoum exposure by stand age(and sex)
nmcmc <- length(beta_stand)
pred_length <- 100
stand_pred <- seq(min(dat$stand_age_mean),max(dat$stand_age_mean),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_build <- mean(dat$nbuildings)
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_mast <- mean(dat$beechnuts)

# Predict
stand.brodM <- matrix(, nmcmc, pred_length)
stand.brodF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  stand.brodF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mean_mast + beta_stand*stand_pred[j]) # males
  stand.brodM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mean_mast + beta_stand*stand_pred[j]) # females
}

# Calculate quantiles
stand.brodF.qt <- apply(stand.brodF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
stand.brodM.qt <- apply(stand.brodM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

stand.qt <- bind_rows(stand.brodF.qt, stand.brodM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
stand_pred <- stand_pred * attr(dat$stand_age_mean, 'scaled:scale') + attr(dat$stand_age_mean, 'scaled:center')

stand.qt.brod <- stand.qt %>% mutate(StandAge=rep(stand_pred, 2))
stand.qt.brod$compound <- "Brodifacoum"

ggplot(stand.qt.brod) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=StandAge, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=StandAge, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_x_continuous(breaks=seq(0,90,10)) +
  ylab("Probability of exposure") + xlab("Stand age (years)") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position="none")


## Predicting brodifacoum exposure by stand age(and sex)
nmcmc <- length(beta_stand)
pred_length <- 100
build_pred <- seq(min(dat$nbuildings),max(dat$nbuildings),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_mast <- mean(dat$beechnuts)
mean_stand <- mean(dat$stand_age_mean)

# Predict
build.brodM <- matrix(, nmcmc, pred_length)
build.brodF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  build.brodF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*build_pred[j] + beta_mast*mean_mast + beta_stand*mean_stand) # males
  build.brodM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*build_pred[j] + beta_mast*mean_mast + beta_stand*mean_stand) # females
}

# Calculate quantiles
build.brodF.qt <- apply(build.brodF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
build.brodM.qt <- apply(build.brodM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

build.qt <- bind_rows(build.brodF.qt, build.brodM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
build_pred <- build_pred * attr(dat$nbuildings, 'scaled:scale') + attr(dat$nbuildings, 'scaled:center')

build.qt.brod <- build.qt %>% mutate(Buildings=rep(build_pred, 2))
build.qt.brod$compound <- "Brodifacoum"


ggplot(build.qt.brod) +
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
age.qt.brod <- age.qt.brod %>% mutate(x="Age") %>% rename(x_val=Age) %>% select(compound, Sex, x, x_val, median:uci)
mast.qt.brod <- mast.qt.brod %>% mutate(x="Beechnuts") %>% rename(x_val=Beechnuts) %>% select(compound, Sex, x, x_val, median:uci)
stand.qt.brod <- stand.qt.brod %>% mutate(x="StandAge") %>% rename(x_val=StandAge) %>% select(compound, Sex, x, x_val, median:uci)
build.qt.brod <- build.qt.brod %>% mutate(x="Buildings") %>% rename(x_val=Buildings) %>% select(compound, Sex, x, x_val, median:uci)

qt.brod <- bind_rows(age.qt.brod, mast.qt.brod, stand.qt.brod, build.qt.brod)

ggplot(qt.brod) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=x_val, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=x_val, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  facet_wrap(vars(x), scales="free_x")

