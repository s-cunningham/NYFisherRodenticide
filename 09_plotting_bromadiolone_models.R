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
brom <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Bromadiolone") %>%
  select(RegionalID,bin.exp)
dat <- left_join(dat, brom, by="RegionalID")
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
brom.out1 <- readRDS("output/model_output/brom.out1.rds")
brom.out1 <- do.call("rbind",brom.out1)
brom.out2 <- readRDS("output/model_output/brom.out2.rds")
brom.out2 <- do.call("rbind",brom.out2)
brom.out3 <- readRDS("output/model_output/brom.out3.rds")
brom.out3 <- do.call("rbind",brom.out3)
brom.out4 <- readRDS("output/model_output/brom.out4.rds")
brom.out4 <- do.call("rbind",brom.out4)
brom.out5 <- readRDS("output/model_output/brom.out5.rds")
brom.out5 <- do.call("rbind",brom.out5)
brom.out6 <- readRDS("output/model_output/brom.out6.rds")
brom.out6 <- do.call("rbind",brom.out6)
brom.out7 <- readRDS("output/model_output/brom.out7.rds")
brom.out7 <- do.call("rbind",brom.out7)
brom.out8 <- readRDS("output/model_output/brom.out8.rds")
brom.out8 <- do.call("rbind",brom.out8)
brom.out9 <- readRDS("output/model_output/brom.out9.rds")
brom.out9 <- do.call("rbind",brom.out9)
brom.out10 <- readRDS("output/model_output/brom.out10.rds")
brom.out10 <- do.call("rbind",brom.out10)

## Combine all iterations
beta_age <- c(brom.out1[,19],brom.out2[,19],brom.out3[,19],brom.out4[,19],brom.out5[,19],
              brom.out6[,19],brom.out7[,19],brom.out8[,19],brom.out9[,19],brom.out10[,19])

# Calculate HDI and quantiles
quantile(beta_age, probs=c(0.025,0.5,0.975))


beta_age2 <- c(brom.out1[,20],brom.out2[,20],brom.out3[,20],brom.out4[,20],brom.out5[,20],
               brom.out6[,20],brom.out7[,20],brom.out8[,20],brom.out9[,20],brom.out10[,20])

# Calculate HDI and quantiles
quantile(beta_age2, probs=c(0.025,0.5,0.975))

beta_build <- c(brom.out1[,21],brom.out2[,21],brom.out3[,21],brom.out4[,21],brom.out5[,21],
                brom.out6[,21],brom.out7[,21],brom.out8[,21],brom.out9[,21],brom.out10[,21])

# Calculate quantiles
quantile(beta_build, probs=c(0.025,0.5,0.975))

beta_mast <- c(brom.out1[,22],brom.out2[,22],brom.out3[,22],brom.out4[,22],brom.out5[,22],
               brom.out6[,22],brom.out7[,22],brom.out8[,22],brom.out9[,22],brom.out10[,22])

# Calculate quantiles
quantile(beta_mast, probs=c(0.025,0.5,0.975))

beta_sex2 <- c(brom.out1[,24],brom.out2[,24],brom.out3[,24],brom.out4[,24],brom.out5[,24],
               brom.out6[,24],brom.out7[,24],brom.out8[,24],brom.out9[,24],brom.out10[,24])

# Calculate quantiles
quantile(beta_sex2, probs=c(0.025,0.5,0.975))


beta_stand <- c(brom.out1[,25],brom.out2[,25],brom.out3[,25],brom.out4[,25],brom.out5[,25],
                brom.out6[,25],brom.out7[,25],brom.out8[,25],brom.out9[,25],brom.out10[,25])

# Calculate quantiles
quantile(beta_stand, probs=c(0.025,0.5,0.975))


##############

## Look at intercept
alpha <- c(rowMeans(brom.out1[,1:18]), rowMeans(brom.out2[,1:18]), rowMeans(brom.out3[,1:18]), rowMeans(brom.out4[,1:18]), rowMeans(brom.out5[,1:18]),
           rowMeans(brom.out6[,1:18]), rowMeans(brom.out7[,1:18]), rowMeans(brom.out8[,1:18]), rowMeans(brom.out9[,1:18]), rowMeans(brom.out10[,1:18]))
plot(density(alpha))

## Predicting bromifacoum exposure by age (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
age_pred <- seq(min(dat$Age),max(dat$Age),length.out=pred_length)
age_pred2 <- age_pred^2

# Average beechnut counts and number of buildlings
mean_build <- mean(dat$nbuildings)
mean_mast <- mean(dat$lag_beechnuts)
mean_stand <- mean(dat$stand_age_mean)

# Predict
age.bromM <- matrix(, nmcmc, pred_length)
age.bromF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  age.bromF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_build*mean_build + beta_mast*mean_mast + beta_stand*mean_stand) # males
  age.bromM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_build*mean_build + beta_mast*mean_mast + beta_stand*mean_stand) # females
}

# Calculate quantiles
age.bromF.qt <- apply(age.bromF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
age.bromM.qt <- apply(age.bromM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

age.qt <- bind_rows(age.bromF.qt, age.bromM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
age_pred <- age_pred * attr(dat$Age, 'scaled:scale') + attr(dat$Age, 'scaled:center')

age.qt.brom <- age.qt %>% mutate(Age=rep(age_pred, 2))
age.qt.brom$compound <- "Bromadiolone"

ggplot(age.qt.brom) +
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



## Predicting bromifacoum exposure by mast cycles (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
mast_pred <- seq(min(dat$lag_beechnuts),max(dat$lag_beechnuts),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_build <- mean(dat$nbuildings)
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_stand <- mean(dat$stand_age_mean)

# Predict
mast.bromM <- matrix(, nmcmc, pred_length)
mast.bromF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  mast.bromF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mast_pred[j] + beta_stand*mean_stand) # males
  mast.bromM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mast_pred[j] + beta_stand*mean_stand) # females
}

# Calculate quantiles
mast.bromF.qt <- apply(mast.bromF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
mast.bromM.qt <- apply(mast.bromM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

mast.qt <- bind_rows(mast.bromF.qt, mast.bromM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
mast_pred <- mast_pred * attr(dat$lag_beechnuts, 'scaled:scale') + attr(dat$lag_beechnuts, 'scaled:center')

mast.qt.brom <- mast.qt %>% mutate(Beechnuts=rep(mast_pred, 2))
mast.qt.brom$compound <- "Bromadiolone"

ggplot(mast.qt.brom) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=Beechnuts, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Beechnuts, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  # scale_x_continuous(breaks=seq(0,9)) +
  ylab("Probability of exposure") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position="none")



## Predicting bromifacoum exposure by stand age(and sex)
nmcmc <- length(beta_stand)
pred_length <- 100
stand_pred <- seq(min(dat$stand_age_mean),max(dat$stand_age_mean),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_build <- mean(dat$nbuildings)
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_mast <- mean(dat$lag_beechnuts)

# Predict
stand.bromM <- matrix(, nmcmc, pred_length)
stand.bromF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  stand.bromF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mean_mast + beta_stand*stand_pred[j]) # males
  stand.bromM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mean_mast + beta_stand*stand_pred[j]) # females
}

# Calculate quantiles
stand.bromF.qt <- apply(stand.bromF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
stand.bromM.qt <- apply(stand.bromM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

stand.qt <- bind_rows(stand.bromF.qt, stand.bromM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
stand_pred <- stand_pred * attr(dat$stand_age_mean, 'scaled:scale') + attr(dat$stand_age_mean, 'scaled:center')

stand.qt.brom <- stand.qt %>% mutate(StandAge=rep(stand_pred, 2))
stand.qt.brom$compound <- "Bromadiolone"

ggplot(stand.qt.brom) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=StandAge, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=StandAge, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_x_continuous(breaks=seq(0,90,10)) +
  ylab("Probability of exposure") + xlab("Stand age (years)") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position="none")


## Predicting bromifacoum exposure by stand age(and sex)
nmcmc <- length(beta_stand)
pred_length <- 100
build_pred <- seq(min(dat$nbuildings),max(dat$nbuildings),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_mast <- mean(dat$lag_beechnuts)
mean_stand <- mean(dat$stand_age_mean)

# Predict
build.bromM <- matrix(, nmcmc, pred_length)
build.bromF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  build.bromF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*build_pred[j] + beta_mast*mean_mast + beta_stand*mean_stand) # males
  build.bromM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*build_pred[j] + beta_mast*mean_mast + beta_stand*mean_stand) # females
}

# Calculate quantiles
build.bromF.qt <- apply(build.bromF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
build.bromM.qt <- apply(build.bromM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

build.qt <- bind_rows(build.bromF.qt, build.bromM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
build_pred <- build_pred * attr(dat$nbuildings, 'scaled:scale') + attr(dat$nbuildings, 'scaled:center')

build.qt.brom <- build.qt %>% mutate(Buildings=rep(build_pred, 2))
build.qt.brom$compound <- "Bromadiolone"

ggplot(build.qt.brom) +
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
age.qt.brom <- age.qt.brom %>% mutate(x="Age") %>% rename(x_val=Age) %>% select(compound, Sex, x, x_val, median:uci)
mast.qt.brom <- mast.qt.brom %>% mutate(x="Beechnuts") %>% rename(x_val=Beechnuts) %>% select(compound, Sex, x, x_val, median:uci)
stand.qt.brom <- stand.qt.brom %>% mutate(x="StandAge") %>% rename(x_val=StandAge) %>% select(compound, Sex, x, x_val, median:uci)
build.qt.brom <- build.qt.brom %>% mutate(x="Buildings") %>% rename(x_val=Buildings) %>% select(compound, Sex, x, x_val, median:uci)

qt.brom <- bind_rows(age.qt.brom, mast.qt.brom, stand.qt.brom, build.qt.brom)

ggplot(qt.brom) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=x_val, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=x_val, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  facet_wrap(vars(x), scales="free_x")


