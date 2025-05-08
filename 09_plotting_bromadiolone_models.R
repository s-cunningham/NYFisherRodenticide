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

beech <- dat %>% select(lag_beechnuts, bin.exp)
intermix <- dat %>% select(totalWUI, bin.exp)

# Scale variables
dat$Age <- scale(dat$Age)
dat$intermix <- scale(dat$intermix)
dat$lag_beechnuts <- scale(dat$lag_beechnuts)

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

beta_age2 <- c(brom.out1[,20],brom.out2[,20],brom.out3[,20],brom.out4[,20],brom.out5[,20],
               brom.out6[,20],brom.out7[,20],brom.out8[,20],brom.out9[,20],brom.out10[,20])

beta_intx <- c(brom.out1[,21],brom.out2[,21],brom.out3[,21],brom.out4[,21],brom.out5[,21],
               brom.out6[,21],brom.out7[,21],brom.out8[,21],brom.out9[,21],brom.out10[,21])

beta_mast <- c(brom.out1[,22],brom.out2[,22],brom.out3[,22],brom.out4[,22],brom.out5[,22],
               brom.out6[,22],brom.out7[,22],brom.out8[,22],brom.out9[,22],brom.out10[,22])

beta_sex2 <- c(brom.out1[,24],brom.out2[,24],brom.out3[,24],brom.out4[,24],brom.out5[,24],
               brom.out6[,24],brom.out7[,24],brom.out8[,24],brom.out9[,24],brom.out10[,24])

beta_wui <- c(brom.out1[,25],brom.out2[,25],brom.out3[,25],brom.out4[,25],brom.out5[,25],
              brom.out6[,25],brom.out7[,25],brom.out8[,25],brom.out9[,25],brom.out10[,25])

## Look at intercept
alpha <- c(rowMeans(brom.out1[,1:18]), rowMeans(brom.out2[,1:18]), rowMeans(brom.out3[,1:18]), rowMeans(brom.out4[,1:18]), rowMeans(brom.out5[,1:18]),
           rowMeans(brom.out6[,1:18]), rowMeans(brom.out7[,1:18]), rowMeans(brom.out8[,1:18]), rowMeans(brom.out9[,1:18]), rowMeans(brom.out10[,1:18]))
# plot(density(alpha))

## Predicting Bromadiolone exposure by age (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
age_pred <- seq(min(dat$Age),max(dat$Age),length.out=pred_length)
age_pred2 <- age_pred^2

# Average beechnut counts and number of buildlings
mean_wui <- mean(dat$intermix)
mean_mast <- mean(dat$lag_beechnuts)

# Predict
age.bromM <- matrix(, nmcmc, pred_length)
age.bromF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  age.bromF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_wui*mean_wui + beta_mast*mean_mast + beta_intx*mean_wui*mean_mast) # males
  age.bromM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_wui*mean_wui + beta_mast*mean_mast + beta_intx*mean_wui*mean_mast) # females
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
        legend.position = c(1,1),
        legend.justification=c(1,1), 
        legend.background = element_rect(fill=NA))

## Predicting Bromadiolone exposure by mast cycles (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
mast_pred <- seq(min(dat$lag_beechnuts),max(dat$lag_beechnuts),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_wui <- mean(dat$intermix)
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2

# Predict
mast.bromM <- matrix(, nmcmc, pred_length)
mast.bromF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  mast.bromF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*mean_wui + beta_mast*mast_pred[j] + beta_intx*mean_wui*mast_pred[j]) # males
  mast.bromM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*mean_wui + beta_mast*mast_pred[j] + beta_intx*mean_wui*mast_pred[j]) # females
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



## Predicting Bromadiolone exposure by stand age(and sex)
nmcmc <- length(beta_wui)
pred_length <- 100
wui_pred <- seq(min(dat$intermix),max(dat$intermix),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_mast <- mean(dat$lag_beechnuts)

# Predict
wui.bromM <- matrix(, nmcmc, pred_length)
wui.bromF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  wui.bromF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_mast*mean_mast + beta_wui*wui_pred[j] + beta_intx*wui_pred[j]*mean_mast) # males
  wui.bromM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_mast*mean_mast + beta_wui*wui_pred[j] + beta_intx*wui_pred[j]*mean_mast) # females
}

# Calculate quantiles
wui.bromF.qt <- apply(wui.bromF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
wui.bromM.qt <- apply(wui.bromM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

wui.qt <- bind_rows(wui.bromF.qt, wui.bromM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
wui_pred <- wui_pred * attr(dat$intermix, 'scaled:scale') + attr(dat$intermix, 'scaled:center')

wui.qt.brom <- wui.qt %>% mutate(Intermix=rep(wui_pred, 2))
wui.qt.brom$compound <- "Bromadiolone"

ggplot(wui.qt.brom) +
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
intx.bromMl <- matrix(, nmcmc, pred_length)
intx.bromMm <- matrix(, nmcmc, pred_length)
intx.bromMh <- matrix(, nmcmc, pred_length)
intx.bromFl <- matrix(, nmcmc, pred_length)
intx.bromFm <- matrix(, nmcmc, pred_length)
intx.bromFh <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  intx.bromFl[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-1.321305) + beta_intx*wui_pred[j]*(-1.321305)) # males
  intx.bromFm[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-0.167061) + beta_intx*wui_pred[j]*(-0.167061)) # males
  intx.bromFh[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(1.078526) + beta_intx*wui_pred[j]*(1.078526)) # males
  intx.bromMl[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-1.321305) + beta_intx*wui_pred[j]*(-1.321305)) # females
  intx.bromMm[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-0.167061) + beta_intx*wui_pred[j]*(-0.167061)) # females
  intx.bromMh[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(1.078526) + beta_intx*wui_pred[j]*(1.078526)) # females
}

# Calculate quantiles
intx.bromFl.qt <- apply(intx.bromFl, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female", Mast="Fail")
intx.bromFm.qt <- apply(intx.bromFm, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female", Mast="Normal")
intx.bromFh.qt <- apply(intx.bromFh, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female", Mast="High")
intx.bromMl.qt <- apply(intx.bromMl, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male", Mast="Fail")
intx.bromMm.qt <- apply(intx.bromMm, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male", Mast="Normal")
intx.bromMh.qt <- apply(intx.bromMh, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male", Mast="High")

intx.qt <- bind_rows(intx.bromFl.qt, intx.bromMl.qt,
                     intx.bromFm.qt, intx.bromMm.qt,
                     intx.bromFh.qt, intx.bromMh.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
intx_pred <- intx_pred * attr(dat$intermix, 'scaled:scale') + attr(dat$intermix, 'scaled:center')

intx.qt.brom <- intx.qt %>% mutate(Intermix=rep(intx_pred, 6))
intx.qt.brom$compound <- "Bromadiolone"
intx.qt.brom$Mast <- factor(intx.qt.brom$Mast, levels=c("Fail", "Normal", "High"))

ggplot(intx.qt.brom) +
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
age.qt.brom <- age.qt.brom %>% mutate(x="Age") %>% rename(x_val=Age) %>% select(compound, Sex, x, x_val, median:uci)
mast.qt.brom <- mast.qt.brom %>% mutate(x="Beechnuts") %>% rename(x_val=Beechnuts) %>% select(compound, Sex, x, x_val, median:uci)
wui.qt.brom <- wui.qt.brom %>% mutate(x="Intermix") %>% rename(x_val=Intermix) %>% select(compound, Sex, x, x_val, median:uci)
intx.qt.brom <- intx.qt.brom %>% mutate(x="%WildlandUrbanIntermix") %>% rename(x_val=Intermix) %>% select(compound, Sex, Mast, x, x_val, median:uci)

brom.data <- list()
brom.data$age.qt.brom <- age.qt.brom
brom.data$mast.qt.brom <- mast.qt.brom
brom.data$wui.qt.brom <- wui.qt.brom
brom.data$intx.qt.brom <- intx.qt.brom

## Save as R data
save(brom.data, file="data/bromadiolone_preds.Rdata")


qt.brom <- bind_rows(age.qt.brom, mast.qt.brom, wui.qt.brom, intx.qt.brom)

ggplot(qt.brom) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=x_val, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=x_val, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  facet_wrap(vars(x), scales="free_x")




beech <- qt.brom %>% filter(x=="Beechnuts")


mean(c(0.1810927,0.2255170))
mean(c(0.3335808,0.3972044))
