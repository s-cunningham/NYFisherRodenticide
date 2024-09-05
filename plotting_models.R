library(tidyverse)


## Read data, remove some columns we don't want to test
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv") %>%
  mutate(age2=Age^2) %>%
  dplyr::select(-edge_density_15, -edge_density_30, -edge_density_45, -build_cat_15, -build_cat_30,-build_cat_45) %>%
  mutate(mast_year=if_else(year==2019, 1, 2), # mast years are reference
         Sex=if_else(Sex=='F',1,2)) # Females are reference

# Scale variables
dat[,c(8,16:42)] <- scale(dat[,c(8,16:42)])

# Load posterior samples
ncomp.out1 <- readRDS("output/model_output/ncomp.out1.rds")
ncomp.out2 <- readRDS("output/model_output/ncomp.out2.rds")
ncomp.out3 <- readRDS("output/model_output/ncomp.out3.rds")
ncomp.out4 <- readRDS("output/model_output/ncomp.out4.rds")
ncomp.out5 <- readRDS("output/model_output/ncomp.out5.rds")
ncomp.out6 <- readRDS("output/model_output/ncomp.out6.rds")
ncomp.out7 <- readRDS("output/model_output/ncomp.out7.rds")
ncomp.out8 <- readRDS("output/model_output/ncomp.out8.rds")
ncomp.out9 <- readRDS("output/model_output/ncomp.out9.rds")
ncomp.out10 <- readRDS("output/model_output/ncomp.out10.rds")


beta_age <- c(ncomp.out1$chain1[,19],ncomp.out2$chain1[,19],ncomp.out3$chain1[,19],ncomp.out4$chain1[,19],ncomp.out5$chain1[,19],
         ncomp.out6$chain1[,19],ncomp.out7$chain1[,19],ncomp.out8$chain1[,19],ncomp.out9$chain1[,19],ncomp.out10$chain1[,19],
         ncomp.out1$chain2[,19],ncomp.out2$chain2[,19],ncomp.out3$chain2[,19],ncomp.out4$chain2[,19],ncomp.out5$chain2[,19],
         ncomp.out6$chain2[,19],ncomp.out7$chain2[,19],ncomp.out8$chain2[,19],ncomp.out9$chain2[,19],ncomp.out10$chain2[,19],
         ncomp.out1$chain3[,19],ncomp.out2$chain3[,19],ncomp.out3$chain3[,19],ncomp.out4$chain3[,19],ncomp.out5$chain3[,19],
         ncomp.out6$chain3[,19],ncomp.out7$chain3[,19],ncomp.out8$chain3[,19],ncomp.out9$chain3[,19],ncomp.out10$chain3[,19])

# Calculate HDI and quantiles
quantile(beta_age, probs=c(0.025,0.5,0.975))


beta_age2 <- c(ncomp.out1$chain1[,20],ncomp.out2$chain1[,20],ncomp.out3$chain1[,20],ncomp.out4$chain1[,20],ncomp.out5$chain1[,20],
          ncomp.out6$chain1[,20],ncomp.out7$chain1[,20],ncomp.out8$chain1[,20],ncomp.out9$chain1[,20],ncomp.out10$chain1[,20],
          ncomp.out1$chain2[,20],ncomp.out2$chain2[,20],ncomp.out3$chain2[,20],ncomp.out4$chain2[,20],ncomp.out5$chain2[,20],
          ncomp.out6$chain2[,20],ncomp.out7$chain2[,20],ncomp.out8$chain2[,20],ncomp.out9$chain2[,20],ncomp.out10$chain2[,20],
          ncomp.out1$chain3[,20],ncomp.out2$chain3[,20],ncomp.out3$chain3[,20],ncomp.out4$chain3[,20],ncomp.out5$chain3[,20],
          ncomp.out6$chain3[,20],ncomp.out7$chain3[,20],ncomp.out8$chain3[,20],ncomp.out9$chain3[,20],ncomp.out10$chain3[,20])

# Calculate HDI and quantiles
quantile(beta_age2, probs=c(0.025,0.5,0.975))

beta_build <- c(ncomp.out1$chain1[,21],ncomp.out2$chain1[,21],ncomp.out3$chain1[,21],ncomp.out4$chain1[,21],ncomp.out5$chain1[,21],
           ncomp.out6$chain1[,21],ncomp.out7$chain1[,21],ncomp.out8$chain1[,21],ncomp.out9$chain1[,21],ncomp.out10$chain1[,21],
           ncomp.out1$chain2[,21],ncomp.out2$chain2[,21],ncomp.out3$chain2[,21],ncomp.out4$chain2[,21],ncomp.out5$chain2[,21],
           ncomp.out6$chain2[,21],ncomp.out7$chain2[,21],ncomp.out8$chain2[,21],ncomp.out9$chain2[,21],ncomp.out10$chain2[,21],
           ncomp.out1$chain3[,21],ncomp.out2$chain3[,21],ncomp.out3$chain3[,21],ncomp.out4$chain3[,21],ncomp.out5$chain3[,21],
           ncomp.out6$chain3[,21],ncomp.out7$chain3[,21],ncomp.out8$chain3[,21],ncomp.out9$chain3[,21],ncomp.out10$chain3[,21])

# Calculate quantiles
quantile(beta_build, probs=c(0.025,0.5,0.975))

beta_mast <- c(ncomp.out1$chain1[,22],ncomp.out2$chain1[,22],ncomp.out3$chain1[,22],ncomp.out4$chain1[,22],ncomp.out5$chain1[,22],
          ncomp.out6$chain1[,22],ncomp.out7$chain1[,22],ncomp.out8$chain1[,22],ncomp.out9$chain1[,22],ncomp.out10$chain1[,22],
          ncomp.out1$chain2[,22],ncomp.out2$chain2[,22],ncomp.out3$chain2[,22],ncomp.out4$chain2[,22],ncomp.out5$chain2[,22],
          ncomp.out6$chain2[,22],ncomp.out7$chain2[,22],ncomp.out8$chain2[,22],ncomp.out9$chain2[,22],ncomp.out10$chain2[,22],
          ncomp.out1$chain3[,22],ncomp.out2$chain3[,22],ncomp.out3$chain3[,22],ncomp.out4$chain3[,22],ncomp.out5$chain3[,22],
          ncomp.out6$chain3[,22],ncomp.out7$chain3[,22],ncomp.out8$chain3[,22],ncomp.out9$chain3[,22],ncomp.out10$chain3[,22])

# Calculate HDI and quantiles
quantile(beta_mast, probs=c(0.025,0.5,0.975))

beta_sex1 <- c(ncomp.out1$chain1[,23],ncomp.out2$chain1[,23],ncomp.out3$chain1[,23],ncomp.out4$chain1[,23],ncomp.out5$chain1[,23],
          ncomp.out6$chain1[,23],ncomp.out7$chain1[,23],ncomp.out8$chain1[,23],ncomp.out9$chain1[,23],ncomp.out10$chain1[,23],
          ncomp.out1$chain2[,23],ncomp.out2$chain2[,23],ncomp.out3$chain2[,23],ncomp.out4$chain2[,23],ncomp.out5$chain2[,23],
          ncomp.out6$chain2[,23],ncomp.out7$chain2[,23],ncomp.out8$chain2[,23],ncomp.out9$chain2[,23],ncomp.out10$chain2[,23],
          ncomp.out1$chain3[,23],ncomp.out2$chain3[,23],ncomp.out3$chain3[,23],ncomp.out4$chain3[,23],ncomp.out5$chain3[,23],
          ncomp.out6$chain3[,23],ncomp.out7$chain3[,23],ncomp.out8$chain3[,23],ncomp.out9$chain3[,23],ncomp.out10$chain3[,23])

# Calculate quantiles
quantile(beta_sex1, probs=c(0.025,0.5,0.975))

beta_sex2 <- c(ncomp.out1$chain1[,24],ncomp.out2$chain1[,24],ncomp.out3$chain1[,24],ncomp.out4$chain1[,24],ncomp.out5$chain1[,24],
          ncomp.out6$chain1[,24],ncomp.out7$chain1[,24],ncomp.out8$chain1[,24],ncomp.out9$chain1[,24],ncomp.out10$chain1[,24],
          ncomp.out1$chain2[,24],ncomp.out2$chain2[,24],ncomp.out3$chain2[,24],ncomp.out4$chain2[,24],ncomp.out5$chain2[,24],
          ncomp.out6$chain2[,24],ncomp.out7$chain2[,24],ncomp.out8$chain2[,24],ncomp.out9$chain2[,24],ncomp.out10$chain2[,24],
          ncomp.out1$chain3[,24],ncomp.out2$chain3[,24],ncomp.out3$chain3[,24],ncomp.out4$chain3[,24],ncomp.out5$chain3[,24],
          ncomp.out6$chain3[,24],ncomp.out7$chain3[,24],ncomp.out8$chain3[,24],ncomp.out9$chain3[,24],ncomp.out10$chain3[,24])

# Calculate quantiles
quantile(beta_sex2, probs=c(0.025,0.5,0.975))



## Plot densities of beta coefficients
sex2 <- as_tibble(data.frame(sexM=sex2))
ggplot(sex2) +
  geom_vline(xintercept=0, color="red") +
  geom_density(aes(x=sexM)) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5))

build <- as_tibble(data.frame(buildings=build))
ggplot(build) +
  geom_vline(xintercept=0, color="red") +
  geom_density(aes(x=buildings)) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5))


mast <- as_tibble(data.frame(beechnuts=mast))
ggplot(mast) +
  geom_vline(xintercept=0, color="red") +
  geom_density(aes(x=beechnuts)) +
  theme_classic() +
  theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5))

## Predicting number of compounds by age (and sex)
nmcmc <- length(mast)
pred_length <- 100
age_pred <- seq(min(dat$Age),max(dat$Age),length.out=pred_length)

age.ncompM <- matrix(, nmcmc, pred_length)
age.ncompF <- matrix(, nmcmc, pred_length)
for (j in 1:pred_length) {
  age.ncompM[,j] <- exp(alpha + beta_sex1*0 +   beta_age*age_pred[j] + beta_age2*age_pred) # males
  age.ncompF[,j] <- exp(alpha + beta_sex2*1 + ) # females
}