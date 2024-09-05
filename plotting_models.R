library(tidyverse)

theme_set(theme_classic())

## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-edge_density_15, -edge_density_30, -edge_density_45, -build_cat_15, -build_cat_30,-build_cat_45) %>%
  mutate(mast_year=if_else(year==2019, 1, 2), # mast years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

beech <- dat %>% select(lag_beechnuts, ncomp)
structures <- dat %>% select(nbuildings_15, ncomp)

# Scale variables
dat$Age <- scale(dat$Age)
dat$nbuildings_15 <- scale(dat$nbuildings_15)
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
beta_age <- c(ncomp.out1[,19],ncomp.out2[,19],ncomp.out3[,19],ncomp.out4[,19],ncomp.out5[,19],
              ncomp.out6[,19],ncomp.out7[,19],ncomp.out8[,19],ncomp.out9[,19],ncomp.out10[,19])

# Calculate HDI and quantiles
quantile(beta_age, probs=c(0.025,0.5,0.975))


beta_age2 <- c(ncomp.out1[,20],ncomp.out2[,20],ncomp.out3[,20],ncomp.out4[,20],ncomp.out5[,20],
               ncomp.out6[,20],ncomp.out7[,20],ncomp.out8[,20],ncomp.out9[,20],ncomp.out10[,20])

# Calculate HDI and quantiles
quantile(beta_age2, probs=c(0.025,0.5,0.975))

beta_build <- c(ncomp.out1[,21],ncomp.out2[,21],ncomp.out3[,21],ncomp.out4[,21],ncomp.out5[,21],
                ncomp.out6[,21],ncomp.out7[,21],ncomp.out8[,21],ncomp.out9[,21],ncomp.out10[,21])

# Calculate quantiles
quantile(beta_build, probs=c(0.025,0.5,0.975))

beta_mast <- c(ncomp.out1[,22],ncomp.out2[,22],ncomp.out3[,22],ncomp.out4[,22],ncomp.out5[,22],
               ncomp.out6[,22],ncomp.out7[,22],ncomp.out8[,22],ncomp.out9[,22],ncomp.out10[,22])

# Calculate HDI and quantiles
quantile(beta_mast, probs=c(0.025,0.5,0.975))

beta_sex2 <- c(ncomp.out1[,24],ncomp.out2[,24],ncomp.out3[,24],ncomp.out4[,24],ncomp.out5[,24],
               ncomp.out6[,24],ncomp.out7[,24],ncomp.out8[,24],ncomp.out9[,24],ncomp.out10[,24])

# Calculate quantiles
quantile(beta_sex2, probs=c(0.025,0.5,0.975))


## Plot densities of beta coefficients
sex2 <- as_tibble(data.frame(sexM=beta_sex2))
ggplot(sex2) +
  geom_vline(xintercept=0, color="red") +
  geom_density(aes(x=sexM)) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5))

build <- as_tibble(data.frame(buildings=beta_build))
ggplot(build) +
  geom_vline(xintercept=0, color="red") +
  geom_density(aes(x=buildings)) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5))

mast <- as_tibble(data.frame(beechnuts=beta_mast))
ggplot(mast) +
  geom_vline(xintercept=0, color="red") +
  geom_density(aes(x=beechnuts)) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5))


## nu
nu <- c(ncomp.out1[,26],ncomp.out2[,26],ncomp.out3[,26],ncomp.out4[,26],ncomp.out5[,26],
        ncomp.out6[,26],ncomp.out7[,26],ncomp.out8[,26],ncomp.out9[,26],ncomp.out10[,26])
plot(density(nu))

## Look at intercept
alpha <- c(rowMeans(ncomp.out1[,1:18]), rowMeans(ncomp.out2[,1:18]), rowMeans(ncomp.out3[,1:18]), rowMeans(ncomp.out4[,1:18]), rowMeans(ncomp.out5[,1:18]),
           rowMeans(ncomp.out6[,1:18]), rowMeans(ncomp.out7[,1:18]), rowMeans(ncomp.out8[,1:18]), rowMeans(ncomp.out9[,1:18]), rowMeans(ncomp.out10[,1:18]))
plot(density(alpha))

## Predicting number of compounds by age (and sex)
nmcmc <- length(beta_age)
pred_length <- 100
age_pred <- seq(min(dat$Age),max(dat$Age),length.out=pred_length)
age_pred2 <- age_pred^2

# Average beechnut counts and number of buildlings
mean_build <- mean(dat$nbuildings_15)
mean_mast <- mean(dat$lag_beechnuts)

# Predict
age.ncompM <- matrix(, nmcmc, pred_length)
age.ncompF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  age.ncompF[,j] <- (exp(alpha + beta_sex1*1 + beta_sex2*0 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_build*mean_build + beta_mast*mean_mast)) #^(1/nu)) + (1/(2*nu)) #- (0.5)# males
  age.ncompM[,j] <- (exp(alpha + beta_sex1*0 + beta_sex2*1 + beta_age*age_pred[j] + beta_age2*(age_pred2[j]) + beta_build*mean_build + beta_mast*mean_mast)) #^(1/nu)) + (1/(2*nu)) #- (0.5) # females
}

# Calculate quantiles
age.ncompF.qt <- apply(age.ncompF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
age.ncompM.qt <- apply(age.ncompM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

age.qt <- bind_rows(age.ncompF.qt, age.ncompM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
age_pred <- age_pred * attr(dat$Age, 'scaled:scale') + attr(dat$Age, 'scaled:center')

age.qt <- age.qt %>% mutate(Age=rep(age_pred, 2))

ggplot(age.qt) +
  coord_cartesian(ylim=c(0, 6)) +
  geom_hline(yintercept=0) +
  geom_ribbon(aes(x=Age, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Age, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position = c(1,1),
        legend.justification=c(1,1), 
        legend.background = element_rect(fill=NA))


## Predicting number of compounds by beechnut count (and sex)
nmcmc <- length(beta_mast)
pred_length <- 100
mast_pred <- seq(min(dat$lag_beechnuts),max(dat$lag_beechnuts),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_build <- mean(dat$nbuildings_15)
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2

# Predict
mast.ncompM <- matrix(, nmcmc, pred_length)
mast.ncompF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  mast.ncompF[,j] <- (exp(alpha + beta_sex1*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mast_pred[j])) #^(1/nu)) + (1/(2*nu)) - (0.5) # males
  mast.ncompM[,j] <- (exp(alpha + beta_sex1*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*mean_build + beta_mast*mast_pred[j])) #^(1/nu)) + (1/(2*nu)) - (0.5) # females
}

# Calculate quantiles
mast.ncompF.qt <- apply(mast.ncompF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
mast.ncompM.qt <- apply(mast.ncompM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

mast.qt <- bind_rows(mast.ncompF.qt, mast.ncompM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
mast_pred <- mast_pred * attr(dat$lag_beechnuts, 'scaled:scale') + attr(dat$lag_beechnuts, 'scaled:center')

mast.qt <- mast.qt %>% mutate(Beechnuts=rep(mast_pred, 2))

ggplot(mast.qt) +
  coord_cartesian(ylim=c(0,6)) +
  geom_ribbon(aes(x=Beechnuts, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Beechnuts, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  theme(panel.border=element_rect(fill=NA, color="black"), 
        legend.position = "none") 



## Predicting number of compounds by number of buildings in 15 km2 buffer (and sex)
nmcmc <- length(beta_build)
pred_length <- 100
build_pred <- seq(min(dat$nbuildings_15),max(dat$nbuildings_15),length.out=pred_length)

# Average beechnut counts and number of buildlings
mean_mast <- mean(dat$lag_beechnuts)
mean_age <- mean(dat$Age)
mean_age2 <- mean_age^2

# Predict
build.ncompM <- matrix(, nmcmc, pred_length)
build.ncompF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  build.ncompF[,j] <- (exp(alpha + beta_sex1*1 + beta_sex2*0 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*build_pred[j] + beta_mast*mean_mast)^(1/nu)) + (1/(2*nu)) - (0.5) # males
  build.ncompM[,j] <- (exp(alpha + beta_sex1*0 + beta_sex2*1 + beta_age*mean_age + beta_age2*mean_age2 + beta_build*build_pred[j] + beta_mast*mean_mast)^(1/nu)) + (1/(2*nu)) - (0.5) # females
}

# Calculate quantiles
build.ncompF.qt <- apply(build.ncompF, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Female")
build.ncompM.qt <- apply(build.ncompM, 2, quantile, probs=c(.5, .025, .975)) |> t() %>% as_tibble() %>% mutate(Sex="Male")

build.qt <- bind_rows(build.ncompF.qt, build.ncompM.qt) %>% rename(median=`50%`, lci=`2.5%`, uci=`97.5%`) 

# Back-transform age prediction values
build_pred <- build_pred * attr(dat$nbuildings_15, 'scaled:scale') + attr(dat$nbuildings_15, 'scaled:center')

build.qt <- build.qt %>% mutate(Buildings=rep(build_pred, 2))

ggplot(build.qt) +
  geom_ribbon(aes(x=Buildings, ymin=lci, ymax=uci, color=Sex, fill=Sex), alpha=.4) +
  geom_line(aes(x=Buildings, y=median, color=Sex), linewidth=1) +
  scale_color_manual(values=c("#1b7837", "#762a83"), name="Sex") +
  scale_fill_manual(values=c("#1b7837", "#762a83"), name="Sex") 






