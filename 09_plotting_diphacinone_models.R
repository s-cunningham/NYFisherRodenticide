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
diph <- read_csv("output/summarized_AR_results.csv") %>% filter(compound=="Diphacinone") %>%
  select(RegionalID,bin.exp)
dat <- left_join(dat, diph, by="RegionalID")
dat <- dat %>% select(RegionalID:Town,bin.exp,deciduous:mast_year)

beech <- dat %>% select(lag_beechnuts, bin.exp)
intermix <- dat %>% select(totalWUI, bin.exp)

# Scale variables
dat$Age <- scale(dat$Age)
dat$intermix <- scale(dat$intermix)
dat$lag_beechnuts <- scale(dat$lag_beechnuts)

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

beta_intx <- c(diph.out1[,21],diph.out2[,21],diph.out3[,21],diph.out4[,21],diph.out5[,21],
               diph.out6[,21],diph.out7[,21],diph.out8[,21],diph.out9[,21],diph.out10[,21])

beta_mast <- c(diph.out1[,22],diph.out2[,22],diph.out3[,22],diph.out4[,22],diph.out5[,22],
               diph.out6[,22],diph.out7[,22],diph.out8[,22],diph.out9[,22],diph.out10[,22])

beta_sex2 <- c(diph.out1[,24],diph.out2[,24],diph.out3[,24],diph.out4[,24],diph.out5[,24],
               diph.out6[,24],diph.out7[,24],diph.out8[,24],diph.out9[,24],diph.out10[,24])

beta_wui <- c(diph.out1[,25],diph.out2[,25],diph.out3[,25],diph.out4[,25],diph.out5[,25],
              diph.out6[,25],diph.out7[,25],diph.out8[,25],diph.out9[,25],diph.out10[,25])

## Look at intercept
alpha <- c(rowMeans(diph.out1[,1:18]), rowMeans(diph.out2[,1:18]), rowMeans(diph.out3[,1:18]), rowMeans(diph.out4[,1:18]), rowMeans(diph.out5[,1:18]),
           rowMeans(diph.out6[,1:18]), rowMeans(diph.out7[,1:18]), rowMeans(diph.out8[,1:18]), rowMeans(diph.out9[,1:18]), rowMeans(diph.out10[,1:18]))
# plot(density(alpha))

## Predicting Diphacinone exposure by age (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
age_pred <- seq(min(dat$Age),max(dat$Age),length.out=pred_length)
age_pred2 <- age_pred^2

# Average beechnut counts and number of buildlings
mean_wui <- mean(dat$intermix)
mean_mast <- mean(dat$lag_beechnuts)

# Predict
age.diphM <- matrix(, nmcmc, pred_length)
age.diphF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  age.diphF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_wui*mean_wui + beta_mast*mean_mast + beta_intx*mean_wui*mean_mast) # males
  age.diphM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_wui*mean_wui + beta_mast*mean_mast + beta_intx*mean_wui*mean_mast) # females
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
        legend.position = c(1,1),
        legend.justification=c(1,1), 
        legend.background = element_rect(fill=NA))

## Predicting Diphacinone exposure by mast cycles (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
mast_pred <- seq(min(dat$lag_beechnuts),max(dat$lag_beechnuts),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_wui <- mean(dat$intermix)
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2

# Predict
mast.diphM <- matrix(, nmcmc, pred_length)
mast.diphF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  mast.diphF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*mean_wui + beta_mast*mast_pred[j] + beta_intx*mean_wui*mast_pred[j]) # males
  mast.diphM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*mean_wui + beta_mast*mast_pred[j] + beta_intx*mean_wui*mast_pred[j]) # females
}

# Calculate quantiles
mast.diphF.qt <- apply(mast.diphF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
mast.diphM.qt <- apply(mast.diphM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

mast.qt <- bind_rows(mast.diphF.qt, mast.diphM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
mast_pred <- mast_pred * attr(dat$lag_beechnuts, 'scaled:scale') + attr(dat$lag_beechnuts, 'scaled:center')

mast.qt.diph <- mast.qt %>% mutate(Beechnuts=rep(mast_pred, 2))
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



## Predicting Diphacinone exposure by stand age(and sex)
nmcmc <- length(beta_wui)
pred_length <- 100
wui_pred <- seq(min(dat$intermix),max(dat$intermix),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_mast <- mean(dat$lag_beechnuts)

# Predict
wui.diphM <- matrix(, nmcmc, pred_length)
wui.diphF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  wui.diphF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_mast*mean_mast + beta_wui*wui_pred[j] + beta_intx*wui_pred[j]*mean_mast) # males
  wui.diphM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_mast*mean_mast + beta_wui*wui_pred[j] + beta_intx*wui_pred[j]*mean_mast) # females
}

# Calculate quantiles
wui.diphF.qt <- apply(wui.diphF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
wui.diphM.qt <- apply(wui.diphM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

wui.qt <- bind_rows(wui.diphF.qt, wui.diphM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
wui_pred <- wui_pred * attr(dat$intermix, 'scaled:scale') + attr(dat$intermix, 'scaled:center')

wui.qt.diph <- wui.qt %>% mutate(Intermix=rep(wui_pred, 2))
wui.qt.diph$compound <- "Diphacinone"

ggplot(wui.qt.diph) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=Intermix, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Intermix, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  # scale_x_continuous(breaks=seq(0,90,10)) +
  ylab("Probability of exposure") + xlab("% Wildland-Urban Intermix") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position="none")


### Interaction plot
nmcmc <- length(beta_wui)
pred_length <- 100
intx_pred <- seq(min(dat$intermix),max(dat$intermix),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2

# Predict
intx.diphMl <- matrix(, nmcmc, pred_length)
intx.diphMm <- matrix(, nmcmc, pred_length)
intx.diphMh <- matrix(, nmcmc, pred_length)
intx.diphFl <- matrix(, nmcmc, pred_length)
intx.diphFm <- matrix(, nmcmc, pred_length)
intx.diphFh <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  intx.diphFl[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-1.321305) + beta_intx*wui_pred[j]*(-1.321305)) # males
  intx.diphFm[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-0.167061) + beta_intx*wui_pred[j]*(-0.167061)) # males
  intx.diphFh[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(1.078526) + beta_intx*wui_pred[j]*(1.078526)) # males
  intx.diphMl[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-1.321305) + beta_intx*wui_pred[j]*(-1.321305)) # females
  intx.diphMm[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-0.167061) + beta_intx*wui_pred[j]*(-0.167061)) # females
  intx.diphMh[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(1.078526) + beta_intx*wui_pred[j]*(1.078526)) # females
}

# Calculate quantiles
intx.diphFl.qt <- apply(intx.diphFl, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female", Mast="Fail")
intx.diphFm.qt <- apply(intx.diphFm, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female", Mast="Normal")
intx.diphFh.qt <- apply(intx.diphFh, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female", Mast="High")
intx.diphMl.qt <- apply(intx.diphMl, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male", Mast="Fail")
intx.diphMm.qt <- apply(intx.diphMm, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male", Mast="Normal")
intx.diphMh.qt <- apply(intx.diphMh, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male", Mast="High")

intx.qt <- bind_rows(intx.diphFl.qt, intx.diphMl.qt,
                     intx.diphFm.qt, intx.diphMm.qt,
                     intx.diphFh.qt, intx.diphMh.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
intx_pred <- intx_pred * attr(dat$intermix, 'scaled:scale') + attr(dat$intermix, 'scaled:center')

intx.qt.diph <- intx.qt %>% mutate(Intermix=rep(intx_pred, 6))
intx.qt.diph$compound <- "Diphacinone"
intx.qt.diph$Mast <- factor(intx.qt.diph$Mast, levels=c("Fail", "Normal", "High"))

ggplot(intx.qt.diph) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=Intermix, ymin=lci, ymax=uci, color=Sex, fill=Sex, linetype=Mast), alpha=.1) +
  geom_line(aes(x=Intermix, y=median, color=Sex,  linetype=Mast), linewidth=1) +
  scale_color_manual(values=c("#1b7837","#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837","#762a83"), name="Sex") +
  facet_grid(Sex~.) + guides(colour="none", fill="none") +
  ylab("Probability of exposure") + 
  theme(panel.border=element_rect(fill=NA, color="black"),
        legend.position=c(0,0),
        legend.justification = c(0,0), 
        legend.background = element_rect(fill=NA))


# Organize for full plot
age.qt.diph <- age.qt.diph %>% mutate(x="Age") %>% rename(x_val=Age) %>% select(compound, Sex, x, x_val, median:uci)
mast.qt.diph <- mast.qt.diph %>% mutate(x="Beechnuts") %>% rename(x_val=Beechnuts) %>% select(compound, Sex, x, x_val, median:uci)
wui.qt.diph <- wui.qt.diph %>% mutate(x="Intermix") %>% rename(x_val=Intermix) %>% select(compound, Sex, x, x_val, median:uci)
intx.qt.diph <- intx.qt.diph %>% mutate(x="%WildlandUrbanIntermix") %>% rename(x_val=Intermix) %>% select(compound, Sex, Mast, x, x_val, median:uci)

diph.data <- list()
diph.data$age.qt.diph <- age.qt.diph
diph.data$mast.qt.diph <- mast.qt.diph
diph.data$wui.qt.diph <- wui.qt.diph
diph.data$intx.qt.diph <- intx.qt.diph

## Save as R data
save(diph.data, file="data/diphacinone_preds.Rdata")



# qt.diph <- bind_rows(age.qt.diph, mast.qt.diph, wui.qt.diph, intx.qt.diph)

ggplot(qt.diph) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=x_val, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=x_val, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  facet_wrap(vars(x), scales="free_x")


############# Plotting all three compounds
