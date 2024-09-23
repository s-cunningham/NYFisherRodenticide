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

beech <- dat %>% select(lag_beechnuts, bin.exp)
intermix <- dat %>% select(totalWUI, bin.exp)

# Scale variables
dat$Age <- scale(dat$Age)
dat$intermix <- scale(dat$intermix)
dat$lag_beechnuts <- scale(dat$lag_beechnuts)

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

beta_age2 <- c(brod.out1[,20],brod.out2[,20],brod.out3[,20],brod.out4[,20],brod.out5[,20],
               brod.out6[,20],brod.out7[,20],brod.out8[,20],brod.out9[,20],brod.out10[,20])

beta_intx <- c(brod.out1[,21],brod.out2[,21],brod.out3[,21],brod.out4[,21],brod.out5[,21],
                brod.out6[,21],brod.out7[,21],brod.out8[,21],brod.out9[,21],brod.out10[,21])

beta_mast <- c(brod.out1[,22],brod.out2[,22],brod.out3[,22],brod.out4[,22],brod.out5[,22],
               brod.out6[,22],brod.out7[,22],brod.out8[,22],brod.out9[,22],brod.out10[,22])

beta_sex2 <- c(brod.out1[,24],brod.out2[,24],brod.out3[,24],brod.out4[,24],brod.out5[,24],
               brod.out6[,24],brod.out7[,24],brod.out8[,24],brod.out9[,24],brod.out10[,24])

beta_wui <- c(brod.out1[,25],brod.out2[,25],brod.out3[,25],brod.out4[,25],brod.out5[,25],
               brod.out6[,25],brod.out7[,25],brod.out8[,25],brod.out9[,25],brod.out10[,25])

## Look at intercept
alpha <- c(rowMeans(brod.out1[,1:18]), rowMeans(brod.out2[,1:18]), rowMeans(brod.out3[,1:18]), rowMeans(brod.out4[,1:18]), rowMeans(brod.out5[,1:18]),
           rowMeans(brod.out6[,1:18]), rowMeans(brod.out7[,1:18]), rowMeans(brod.out8[,1:18]), rowMeans(brod.out9[,1:18]), rowMeans(brod.out10[,1:18]))
# plot(density(alpha))

## Predicting brodifacoum exposure by age (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
age_pred <- seq(min(dat$Age),max(dat$Age),length.out=pred_length)
age_pred2 <- age_pred^2

# Average beechnut counts and number of buildlings
mean_wui <- mean(dat$intermix)
mean_mast <- mean(dat$lag_beechnuts)

# Predict
age.brodM <- matrix(, nmcmc, pred_length)
age.brodF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  age.brodF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_wui*mean_wui + beta_mast*mean_mast + beta_intx*mean_wui*mean_mast) # males
  age.brodM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_wui*mean_wui + beta_mast*mean_mast + beta_intx*mean_wui*mean_mast) # females
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
mast_pred <- seq(min(dat$lag_beechnuts),max(dat$lag_beechnuts),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_wui <- mean(dat$intermix)
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2

# Predict
mast.brodM <- matrix(, nmcmc, pred_length)
mast.brodF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  mast.brodF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*mean_wui + beta_mast*mast_pred[j] + beta_intx*mean_wui*mast_pred[j]) # males
  mast.brodM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*mean_wui + beta_mast*mast_pred[j] + beta_intx*mean_wui*mast_pred[j]) # females
}

# Calculate quantiles
mast.brodF.qt <- apply(mast.brodF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
mast.brodM.qt <- apply(mast.brodM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

mast.qt <- bind_rows(mast.brodF.qt, mast.brodM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
mast_pred <- mast_pred * attr(dat$lag_beechnuts, 'scaled:scale') + attr(dat$lag_beechnuts, 'scaled:center')

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
nmcmc <- length(beta_wui)
pred_length <- 100
wui_pred <- seq(min(dat$intermix),max(dat$intermix),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_mast <- mean(dat$lag_beechnuts)

# Predict
wui.brodM <- matrix(, nmcmc, pred_length)
wui.brodF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  wui.brodF[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_mast*mean_mast + beta_wui*wui_pred[j] + beta_intx*wui_pred[j]*mean_mast) # males
  wui.brodM[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_mast*mean_mast + beta_wui*wui_pred[j] + beta_intx*wui_pred[j]*mean_mast) # females
}

# Calculate quantiles
wui.brodF.qt <- apply(wui.brodF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
wui.brodM.qt <- apply(wui.brodM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

wui.qt <- bind_rows(wui.brodF.qt, wui.brodM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
wui_pred <- wui_pred * attr(dat$intermix, 'scaled:scale') + attr(dat$intermix, 'scaled:center')

wui.qt.brod <- wui.qt %>% mutate(Intermix=rep(wui_pred, 2))
wui.qt.brod$compound <- "Brodifacoum"

ggplot(wui.qt.brod) +
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
nmcmc <- length(beta_intx)
pred_length <- 100
intx_pred <- seq(min(dat$intermix),max(dat$intermix),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2

# Predict
intx.brodMl <- matrix(, nmcmc, pred_length)
intx.brodMm <- matrix(, nmcmc, pred_length)
intx.brodMh <- matrix(, nmcmc, pred_length)
intx.brodFl <- matrix(, nmcmc, pred_length)
intx.brodFm <- matrix(, nmcmc, pred_length)
intx.brodFh <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  intx.brodFl[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-1.321305) + beta_intx*wui_pred[j]*(-1.321305)) # males
  intx.brodFm[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-0.167061) + beta_intx*wui_pred[j]*(-0.167061)) # males
  intx.brodFh[,j] <- inv.logit(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(1.078526) + beta_intx*wui_pred[j]*(1.078526)) # males
  intx.brodMl[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-1.321305) + beta_intx*wui_pred[j]*(-1.321305)) # females
  intx.brodMm[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-0.167061) + beta_intx*wui_pred[j]*(-0.167061)) # females
  intx.brodMh[,j] <- inv.logit(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(1.078526) + beta_intx*wui_pred[j]*(1.078526)) # females
}

# Calculate quantiles
intx.brodFl.qt <- apply(intx.brodFl, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female", Mast="Fail")
intx.brodFm.qt <- apply(intx.brodFm, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female", Mast="Normal")
intx.brodFh.qt <- apply(intx.brodFh, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female", Mast="High")
intx.brodMl.qt <- apply(intx.brodMl, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male", Mast="Fail")
intx.brodMm.qt <- apply(intx.brodMm, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male", Mast="Normal")
intx.brodMh.qt <- apply(intx.brodMh, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male", Mast="High")

intx.qt <- bind_rows(intx.brodFl.qt, intx.brodMl.qt,
                     intx.brodFm.qt, intx.brodMm.qt,
                     intx.brodFh.qt, intx.brodMh.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
intx_pred <- intx_pred * attr(dat$intermix, 'scaled:scale') + attr(dat$intermix, 'scaled:center')

intx.qt.brod <- intx.qt %>% mutate(Intermix=rep(intx_pred, 6))
intx.qt.brod$compound <- "Brodifacoum"
intx.qt.brod$Mast <- factor(intx.qt.brod$Mast, levels=c("Fail", "Normal", "High"))

ggplot(intx.qt.brod) +
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
age.qt.brod <- age.qt.brod %>% mutate(x="Age") %>% rename(x_val=Age) %>% select(compound, Sex, x, x_val, median:uci)
mast.qt.brod <- mast.qt.brod %>% mutate(x="Beechnuts") %>% rename(x_val=Beechnuts) %>% select(compound, Sex, x, x_val, median:uci)
wui.qt.brod <- wui.qt.brod %>% mutate(x="Intermix") %>% rename(x_val=Intermix) %>% select(compound, Sex, x, x_val, median:uci)
intx.qt.brod <- intx.qt.brod %>% mutate(x="%WildlandUrbanIntermix") %>% rename(x_val=Intermix) %>% select(compound, Sex, Mast, x, x_val, median:uci)


brod.data <- list()
brod.data$age.qt.brod <- age.qt.brod
brod.data$mast.qt.brod <- mast.qt.brod
brod.data$wui.qt.brod <- wui.qt.brod
brod.data$intx.qt.brod <- intx.qt.brod

## Save as R data
save(brod.data, file="data/brodifacoum_preds.Rdata")


# qt.brod <- bind_rows(age.qt.brod, mast.qt.brod, wui.qt.brod, intx.qt.brod)

ggplot(qt.brod) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_ribbon(aes(x=x_val, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=x_val, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  facet_wrap(vars(x),nrow=1, scales="free_x")

