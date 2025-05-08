library(tidyverse)

theme_set(theme_classic())

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-build_cat) %>%
  mutate(mast_year=if_else(year==2019, 2, 1), # failure years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat$Age <- scale(dat$Age)
dat$intermix <- scale(dat$intermix)
dat$lag_beechnuts <- scale(dat$lag_beechnuts)

# Load posterior samples
ncomp.out1 <- readRDS("output/model_output/ncomp.out1.rds")
ncomp.out1 <- do.call("rbind",ncomp.out1)
ncomp.out2 <- readRDS("output/model_output/ncomp.out2.rds")
ncomp.out2 <- do.call("rbind",ncomp.out2)
ncomp.out3 <- readRDS("output/model_output/ncomp.out3.rds")
ncomp.out3 <- do.call("rbind",ncomp.out3)
ncomp.out4 <- readRDS("output/model_output/ncomp.out4.rds")
ncomp.out4 <- do.call("rbind",ncomp.out4)
ncomp.out5 <- readRDS("output/model_output/ncomp.out5.rds")
ncomp.out5 <- do.call("rbind",ncomp.out5)
ncomp.out6 <- readRDS("output/model_output/ncomp.out6.rds")
ncomp.out6 <- do.call("rbind",ncomp.out6)
ncomp.out7 <- readRDS("output/model_output/ncomp.out7.rds")
ncomp.out7 <- do.call("rbind",ncomp.out7)
ncomp.out8 <- readRDS("output/model_output/ncomp.out8.rds")
ncomp.out8 <- do.call("rbind",ncomp.out8)
ncomp.out9 <- readRDS("output/model_output/ncomp.out9.rds")
ncomp.out9 <- do.call("rbind",ncomp.out9)
ncomp.out10 <- readRDS("output/model_output/ncomp.out10.rds")
ncomp.out10 <- do.call("rbind",ncomp.out10)

## Combine all iterations
beta_age <- c(ncomp.out1[,19],ncomp.out2[,19],ncomp.out3[,19],ncomp.out5[,19],ncomp.out4[,19],
              ncomp.out6[,19],ncomp.out7[,19],ncomp.out8[,19],ncomp.out9[,19],ncomp.out10[,19])

beta_age2 <- c(ncomp.out1[,20],ncomp.out2[,20],ncomp.out3[,20],ncomp.out5[,20],ncomp.out4[,20],
               ncomp.out6[,20],ncomp.out7[,20],ncomp.out8[,20],ncomp.out9[,20],ncomp.out10[,20])

beta_intx <- c(ncomp.out1[,21],ncomp.out2[,21],ncomp.out3[,21],ncomp.out5[,21],ncomp.out4[,21],
               ncomp.out6[,21],ncomp.out7[,21],ncomp.out8[,21],ncomp.out9[,21],ncomp.out10[,21])

beta_mast <- c(ncomp.out1[,22],ncomp.out2[,22],ncomp.out3[,22],ncomp.out5[,22],ncomp.out4[,22],
               ncomp.out6[,22],ncomp.out7[,22],ncomp.out8[,22],ncomp.out9[,22],ncomp.out10[,22])

beta_sex2 <- c(ncomp.out1[,24],ncomp.out2[,24],ncomp.out3[,24],ncomp.out5[,24],ncomp.out4[,24],
               ncomp.out6[,24],ncomp.out7[,24],ncomp.out8[,24],ncomp.out9[,24],ncomp.out10[,24])

beta_wui <- c(ncomp.out1[,25],ncomp.out2[,25],ncomp.out3[,25],ncomp.out5[,25],ncomp.out4[,25],
              ncomp.out6[,25],ncomp.out7[,25],ncomp.out8[,25],ncomp.out9[,25],ncomp.out10[,25])

## Look at intercept
alpha <- c(rowMeans(ncomp.out1[,1:18]), rowMeans(ncomp.out2[,1:18]), rowMeans(ncomp.out3[,1:18]), rowMeans(ncomp.out5[,1:18]), rowMeans(ncomp.out4[,1:18]),
           rowMeans(ncomp.out6[,1:18]), rowMeans(ncomp.out7[,1:18]), rowMeans(ncomp.out8[,1:18]), rowMeans(ncomp.out9[,1:18]), rowMeans(ncomp.out10[,1:18]))
# plot(density(alpha))

## Look at nu
nu <- c(ncomp.out1[,27],ncomp.out2[,27],ncomp.out3[,27],ncomp.out5[,27],ncomp.out4[,25],
        ncomp.out6[,27],ncomp.out7[,27],ncomp.out8[,27],ncomp.out9[,27],ncomp.out10[,27])
plot(density(nu))
quantile(nu, probs=c(0.025, 0.5, 0.975))
mean_nu <- mean(nu)


## Predicting ncompifacoum exposure by age (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
age_pred <- seq(min(dat$Age),max(dat$Age),length.out=pred_length)
age_pred2 <- age_pred^2

# Average beechnut counts and number of buildlings
mean_wui <- mean(dat$intermix)
mean_mast <- mean(dat$lag_beechnuts)

# Predict
age.ncompM <- matrix(, nmcmc, pred_length)
age.ncompF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  age.ncompF[,j] <- exp(alpha + 0*1 + beta_sex2*0 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_wui*mean_wui + beta_mast*mean_mast + beta_intx*mean_wui*mean_mast)  # males
  age.ncompM[,j] <- exp(alpha + 0*0 + beta_sex2*1 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_wui*mean_wui + beta_mast*mean_mast + beta_intx*mean_wui*mean_mast)  # females
}

# Calculate quantiles
age.ncompF.qt <- apply(age.ncompF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
age.ncompM.qt <- apply(age.ncompM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

age.qt <- bind_rows(age.ncompF.qt, age.ncompM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
age_pred <- age_pred * attr(dat$Age, 'scaled:scale') + attr(dat$Age, 'scaled:center')

age.qt.ncomp <- age.qt %>% mutate(Age=rep(age_pred, 2))

mean(age.qt.ncomp$median)
var(age.qt.ncomp$median)

ggplot(age.qt.ncomp) +
  coord_cartesian(xlim=c(0,8.5)) +
  geom_hline(yintercept=0) +
  geom_ribbon(aes(x=Age, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Age, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_x_continuous(breaks=seq(0,9)) +
  ylab("Predicted number of compounds") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position = c(1,1),
        legend.justification=c(1,1), 
        legend.background = element_rect(fill=NA))

## Predicting ncompifacoum exposure by mast cycles (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
mast_pred <- seq(min(dat$lag_beechnuts),max(dat$lag_beechnuts),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_wui <- mean(dat$intermix)
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2

# Predict
mast.ncompM <- matrix(, nmcmc, pred_length)
mast.ncompF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  mast.ncompF[,j] <- exp(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*mean_wui + beta_mast*mast_pred[j] + beta_intx*mean_wui*mast_pred[j]) # males
  mast.ncompM[,j] <- exp(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*mean_wui + beta_mast*mast_pred[j] + beta_intx*mean_wui*mast_pred[j]) # females
}

# Calculate quantiles
mast.ncompF.qt <- apply(mast.ncompF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
mast.ncompM.qt <- apply(mast.ncompM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

mast.qt <- bind_rows(mast.ncompF.qt, mast.ncompM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
mast_pred <- mast_pred * attr(dat$lag_beechnuts, 'scaled:scale') + attr(dat$lag_beechnuts, 'scaled:center')

mast.qt.ncomp <- mast.qt %>% mutate(Beechnuts=rep(mast_pred, 2))

ggplot(mast.qt.ncomp) +
  coord_cartesian(ylim=c(0, 8)) +
  geom_ribbon(aes(x=Beechnuts, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Beechnuts, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  # scale_x_continuous(breaks=seq(0,9)) +
  ylab("Predicted number of compounds") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position="none")

## Predicting ncompifacoum exposure by WUI intermix
nmcmc <- length(beta_wui)
pred_length <- 100
wui_pred <- seq(min(dat$intermix),max(dat$intermix),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2
mean_mast <- mean(dat$lag_beechnuts)

# Predict
wui.ncompM <- matrix(, nmcmc, pred_length)
wui.ncompF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  wui.ncompF[,j] <- exp(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_mast*mean_mast + beta_wui*wui_pred[j] + beta_intx*wui_pred[j]*mean_mast) # males
  wui.ncompM[,j] <- exp(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_mast*mean_mast + beta_wui*wui_pred[j] + beta_intx*wui_pred[j]*mean_mast) # females
}

# Calculate quantiles
wui.ncompF.qt <- apply(wui.ncompF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
wui.ncompM.qt <- apply(wui.ncompM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

wui.qt <- bind_rows(wui.ncompF.qt, wui.ncompM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
wui_pred <- wui_pred * attr(dat$intermix, 'scaled:scale') + attr(dat$intermix, 'scaled:center')

wui.qt.ncomp <- wui.qt %>% mutate(Intermix=rep(wui_pred, 2))

ggplot(wui.qt.ncomp) +
  coord_cartesian(ylim=c(0, 10)) +
  geom_ribbon(aes(x=Intermix, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Intermix, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  # scale_x_continuous(breaks=seq(0,90,10)) +
  ylab("Predicted number of compounds") + xlab("% Wildland-Urban Intermix") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position="none")


# Organize for full plot
age.qt.ncomp <- age.qt.ncomp %>% mutate(x="Age") %>% rename(x_val=Age) %>% select(Sex, x, x_val, median:uci)
mast.qt.ncomp <- mast.qt.ncomp %>% mutate(x="Beechnuts") %>% rename(x_val=Beechnuts) %>% select(Sex, x, x_val, median:uci)
wui.qt.ncomp <- wui.qt.ncomp %>% mutate(x="Intermix") %>% rename(x_val=Intermix) %>% select(Sex, x, x_val, median:uci)

qt <- bind_rows(age.qt.ncomp, mast.qt.ncomp, wui.qt.ncomp)

qt <- qt %>% mutate(x=case_when(x=="Age" ~ "Age (years)",
                                        x=="Beechnuts" ~ "Difference in\nbeech seed count (1-yr lag)",
                                        x=="Intermix" ~ "% Wildland-urban intermix"))
qt$x <- factor(qt$x, levels=c("Age (years)", "Difference in\nbeech seed count (1-yr lag)", "% Wildland-urban intermix"))

vline.data <- data.frame(z=c(NA, 6,145,295,NA), x=c("Age (years)", rep("Difference in\nbeech seed count (1-yr lag)",3),"% Wildland-urban intermix"))
vline.data$x <- factor(vline.data$x, levels=c("Age (years)", "Difference in\nbeech seed count (1-yr lag)", "% Wildland-urban intermix"))

ggplot(qt) +
  coord_cartesian(ylim=c(0,8)) +
  # geom_vline(aes(xintercept = z), data=vline.data, colour = "gray85")+
  geom_ribbon(aes(x=x_val, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=x_val, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  facet_grid(.~x, scales="free_x", switch="both", axes = "all", axis.labels = "margins") +
  ylab("Predicted number of compounds") +
  theme(legend.position=c(0.325,1), legend.justification=c(1,1),
        panel.border=element_rect(fill=NA, color="black"),
        strip.placement = "outside",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text=element_text(size=11),
        strip.text.x=element_text(size=12),
        strip.text.y=element_text(size=13),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11),
        legend.background = element_rect(fill=NA),
        strip.background = element_rect(color=NA, fill=NA),
        axis.ticks.length=unit(-0.1, "cm"))

ggsave("figs/ncomp_preds.svg")
# Saving 10.1 x 4.04 in image



mast.qt.ncomp



mean(c(3.309645, 2.538521))

mean(c(3.233165,4.214554))










### Interaction plot
# nmcmc <- length(beta_intx)
# pred_length <- 100
# intx_pred <- seq(min(dat$intermix),max(dat$intermix),length.out=pred_length)
# 
# # Average beechnut counts and number of buildlings
# mean_age <- mean(dat$Age)
# mean_age2 <- mean_age^2
# 
# # Predict
# intx.ncompMl <- matrix(, nmcmc, pred_length)
# intx.ncompMm <- matrix(, nmcmc, pred_length)
# intx.ncompMh <- matrix(, nmcmc, pred_length)
# intx.ncompFl <- matrix(, nmcmc, pred_length)
# intx.ncompFm <- matrix(, nmcmc, pred_length)
# intx.ncompFh <- matrix(, nmcmc, pred_length)
# for (j in 1:pred_length) {
#   intx.ncompFl[,j] <- exp(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-1.321305) + beta_intx*wui_pred[j]*(-1.321305)) # males
#   intx.ncompFm[,j] <- exp(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-0.167061) + beta_intx*wui_pred[j]*(-0.167061)) # males
#   intx.ncompFh[,j] <- exp(alpha + 0*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(1.078526) + beta_intx*wui_pred[j]*(1.078526)) # males
#   intx.ncompMl[,j] <- exp(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-1.321305) + beta_intx*wui_pred[j]*(-1.321305)) # females
#   intx.ncompMm[,j] <- exp(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(-0.167061) + beta_intx*wui_pred[j]*(-0.167061)) # females
#   intx.ncompMh[,j] <- exp(alpha + 0*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_wui*wui_pred[j] + beta_mast*(1.078526) + beta_intx*wui_pred[j]*(1.078526)) # females
# }
# 
# # Calculate quantiles
# intx.ncompFl.qt <- apply(intx.ncompFl, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female", Mast="Fail")
# intx.ncompFm.qt <- apply(intx.ncompFm, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female", Mast="Normal")
# intx.ncompFh.qt <- apply(intx.ncompFh, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female", Mast="High")
# intx.ncompMl.qt <- apply(intx.ncompMl, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male", Mast="Fail")
# intx.ncompMm.qt <- apply(intx.ncompMm, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male", Mast="Normal")
# intx.ncompMh.qt <- apply(intx.ncompMh, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male", Mast="High")
# 
# intx.qt <- bind_rows(intx.ncompFl.qt, intx.ncompMl.qt,
#                      intx.ncompFm.qt, intx.ncompMm.qt,
#                      intx.ncompFh.qt, intx.ncompMh.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 
# 
# # Back-transform age prediction values
# intx_pred <- intx_pred * attr(dat$intermix, 'scaled:scale') + attr(dat$intermix, 'scaled:center')
# 
# intx.qt.ncomp <- intx.qt %>% mutate(Intermix=rep(intx_pred, 6))
# intx.qt.ncomp$Mast <- factor(intx.qt.ncomp$Mast, levels=c("High", "Normal", "Fail"))
# 
# ggplot(intx.qt.ncomp) +
#   coord_cartesian(ylim=c(0, 8)) +
#   geom_ribbon(aes(x=Intermix, ymin=lci, ymax=uci, color=Sex, fill=Sex, linetype=Mast), alpha=.1) +
#   geom_line(aes(x=Intermix, y=median, color=Sex,  linetype=Mast), linewidth=1) +
#   scale_color_manual(values=c("#1b7837","#762a83"), name="Sex") +
#   scale_fill_manual(values=c("#1b7837","#762a83"), name="Sex") +
#   facet_grid(.~Sex) + guides(colour="none", fill="none") +
#   ylab("Predicted number of compounds") + 
#   theme(panel.border=element_rect(fill=NA, color="black"),
#         legend.position=c(0,0),
#         legend.justification = c(0,0), 
#         legend.background = element_rect(fill=NA))


# intx.qt.ncomp <- intx.qt.ncomp %>% mutate(x="%WildlandUrbanIntermix") %>% rename(x_val=Intermix) %>% select(compound, Sex, Mast, x, x_val, median:uci)

