library(tidyverse)
library(MCMCvis)

se <- function(x) {
  out <- sd(x)/sqrt(length((x)))
  return(out)
}

#### Diphacinone ####
# Load posterior samples
diph.out1 <- readRDS("output/model_output/diph.out1.rds")
diph1.sum <- MCMCsummary(diph.out1) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=1)

diph.out2 <- readRDS("output/model_output/diph.out2.rds")
diph2.sum <- MCMCsummary(diph.out2) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=2)
diph.out3 <- readRDS("output/model_output/diph.out3.rds")
diph3.sum <- MCMCsummary(diph.out3) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=3)
diph.out4 <- readRDS("output/model_output/diph.out4.rds")
diph4.sum <- MCMCsummary(diph.out4) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=4)
diph.out5 <- readRDS("output/model_output/diph.out5.rds")
diph5.sum <- MCMCsummary(diph.out5) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=5)
diph.out6 <- readRDS("output/model_output/diph.out6.rds")
diph6.sum <- MCMCsummary(diph.out6) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=6)
diph.out7 <- readRDS("output/model_output/diph.out7.rds")
diph7.sum <- MCMCsummary(diph.out7) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=7)
diph.out8 <- readRDS("output/model_output/diph.out8.rds")
diph8.sum <- MCMCsummary(diph.out8) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=8)
diph.out9 <- readRDS("output/model_output/diph.out9.rds")
diph9.sum <- MCMCsummary(diph.out9) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=9)
diph.out10 <- readRDS("output/model_output/diph.out10.rds")
diph10.sum <- MCMCsummary(diph.out10) %>% rownames_to_column("param") %>% as_tibble() %>%
                  mutate(iter=10)

diph.sum <- bind_rows(diph1.sum, diph2.sum, diph3.sum, diph4.sum, diph5.sum, 
                      diph6.sum, diph7.sum, diph8.sum, diph9.sum, diph10.sum)

diph.sum <- diph.sum %>% mutate(compound="diphacinone") %>%
                  select(compound, iter, param:n.eff) %>%
                  # exponentiate for odds ratios 
                  mutate(mean=exp(mean),
                         `2.5%`=exp(`2.5%`),
                         `50%`=exp(`50%`),
                         `97.5%`=exp(`97.5%`))

diph.odds <- diph.sum %>% group_by(param) %>% reframe(mean=mean(mean), `2.5%`=mean(`2.5%`), 
                                         `50%`=mean(`50%`), `97.5%`=mean(`97.5%`)) %>% 
  filter(str_detect(param, 'beta'))

## Combine all iterations
diph.out1 <- do.call("rbind",diph.out1)
diph.out2 <- do.call("rbind",diph.out2)
diph.out3 <- do.call("rbind",diph.out3)
diph.out4 <- do.call("rbind",diph.out4)
diph.out5 <- do.call("rbind",diph.out5)
diph.out6 <- do.call("rbind",diph.out6)
diph.out7 <- do.call("rbind",diph.out7)
diph.out8 <- do.call("rbind",diph.out8)
diph.out9 <- do.call("rbind",diph.out9)
diph.out10 <- do.call("rbind",diph.out10)

## Age
beta_age <- c(diph.out1[,19],diph.out2[,19],diph.out3[,19],diph.out4[,19],diph.out5[,19],
              diph.out6[,19],diph.out7[,19],diph.out8[,19],diph.out9[,19],diph.out10[,19])

# Calculate quantiles
quantile(beta_age, probs=c(0.025,0.5,0.975))
exp(mean(beta_age))
se(beta_age)

exp_age <- exp(beta_age)
sum(exp_age>1)/length(beta_age)

## Age^2
beta_age2 <- c(diph.out1[,20],diph.out2[,20],diph.out3[,20],diph.out4[,20],diph.out5[,20],
               diph.out6[,20],diph.out7[,20],diph.out8[,20],diph.out9[,20],diph.out10[,20])

# Calculate quantiles
quantile(beta_age2, probs=c(0.025,0.5,0.975))
mean(beta_age2)
se(beta_age2)
plot(density(beta_age2))

exp_age2 <- exp(beta_age2)
sum(exp_age2<1)/length(beta_age2)

plot(density(exp_age2))

## Building count
beta_build <- c(diph.out1[,21],diph.out2[,21],diph.out3[,21],diph.out4[,21],diph.out5[,21],
                diph.out6[,21],diph.out7[,21],diph.out8[,21],diph.out9[,21],diph.out10[,21])

# Calculate quantiles
quantile(beta_build, probs=c(0.025,0.5,0.975))
mean(exp(beta_build))
se(beta_build)

beta_build <- exp(beta_build)
sum(beta_build>1)/length(beta_build)

## Beech mast
beta_mast <- c(diph.out1[,22],diph.out2[,22],diph.out3[,22],diph.out4[,22],diph.out5[,22],
               diph.out6[,22],diph.out7[,22],diph.out8[,22],diph.out9[,22],diph.out10[,22])

# Calculate quantiles
quantile(beta_mast, probs=c(0.025,0.5,0.975))
mean(beta_mast)
se(beta_mast)


mean(exp(beta_mast))
beta_mast <- exp(beta_mast)
sum(beta_mast>1)/length(beta_mast)

## Sex (M)
beta_sex2 <- c(diph.out1[,24],diph.out2[,24],diph.out3[,24],diph.out4[,24],diph.out5[,24],
               diph.out6[,24],diph.out7[,24],diph.out8[,24],diph.out9[,24],diph.out10[,24])

# Calculate quantiles
quantile(beta_sex2, probs=c(0.025,0.5,0.975))
mean(beta_sex2)
se(beta_sex2)

mean(exp(beta_sex2))
beta_sex2 <- exp(beta_sex2)
sum(beta_sex2>1)/length(beta_sex2)

## Stand age
beta_stand <- c(diph.out1[,25],diph.out2[,25],diph.out3[,25],diph.out4[,25],diph.out5[,25],
                diph.out6[,25],diph.out7[,25],diph.out8[,25],diph.out9[,25],diph.out10[,25])

# Calculate quantiles
quantile(beta_stand, probs=c(0.025,0.5,0.975))
mean(beta_stand)
se(beta_stand)

mean(exp(beta_stand))
beta_stand <- exp(beta_stand)
sum(beta_stand>1)/length(beta_stand)

#### Brodifacoum ####
# Load posterior samples
brod.out1 <- readRDS("output/model_output/brod.out1.rds")
brod1.sum <- MCMCsummary(brod.out1) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=1)

brod.out2 <- readRDS("output/model_output/brod.out2.rds")
brod2.sum <- MCMCsummary(brod.out2) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=2)
brod.out3 <- readRDS("output/model_output/brod.out3.rds")
brod3.sum <- MCMCsummary(brod.out3) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=3)
brod.out4 <- readRDS("output/model_output/brod.out4.rds")
brod4.sum <- MCMCsummary(brod.out4) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=4)
brod.out5 <- readRDS("output/model_output/brod.out5.rds")
brod5.sum <- MCMCsummary(brod.out5) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=5)
brod.out6 <- readRDS("output/model_output/brod.out6.rds")
brod6.sum <- MCMCsummary(brod.out6) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=6)
brod.out7 <- readRDS("output/model_output/brod.out7.rds")
brod7.sum <- MCMCsummary(brod.out7) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=7)
brod.out8 <- readRDS("output/model_output/brod.out8.rds")
brod8.sum <- MCMCsummary(brod.out8) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=8)
brod.out9 <- readRDS("output/model_output/brod.out9.rds")
brod9.sum <- MCMCsummary(brod.out9) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=9)
brod.out10 <- readRDS("output/model_output/brod.out10.rds")
brod10.sum <- MCMCsummary(brod.out10) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=10)

brod.sum <- bind_rows(brod1.sum, brod2.sum, brod3.sum, brod4.sum, brod5.sum, 
                      brod6.sum, brod7.sum, brod8.sum, brod9.sum, brod10.sum)

brod.sum <- brod.sum %>% mutate(compound="brodifacoum") %>%
  select(compound, iter, param:n.eff) %>%
  # exponentiate for odds ratios 
  mutate(mean=exp(mean),
         `2.5%`=exp(`2.5%`),
         `50%`=exp(`50%`),
         `97.5%`=exp(`97.5%`))

brod.odds <- brod.sum %>% group_by(param) %>% reframe(mean=mean(mean), `2.5%`=mean(`2.5%`), 
                                                 `50%`=mean(`50%`), `97.5%`=mean(`97.5%`)) %>% 
  filter(str_detect(param, 'beta'))

## Combine all iterations
brod.out1 <- do.call("rbind",brod.out1)
brod.out2 <- do.call("rbind",brod.out2)
brod.out3 <- do.call("rbind",brod.out3)
brod.out4 <- do.call("rbind",brod.out4)
brod.out5 <- do.call("rbind",brod.out5)
brod.out6 <- do.call("rbind",brod.out6)
brod.out7 <- do.call("rbind",brod.out7)
brod.out8 <- do.call("rbind",brod.out8)
brod.out9 <- do.call("rbind",brod.out9)
brod.out10 <- do.call("rbind",brod.out10)

### Age
beta_age <- c(brod.out1[,19],brod.out2[,19],brod.out3[,19],brod.out4[,19],brod.out5[,19],
              brod.out6[,19],brod.out7[,19],brod.out8[,19],brod.out9[,19],brod.out10[,19])

mean(exp(beta_age))
beta_age <- exp(beta_age)
sum(beta_age>1)/length(beta_age)

# Calculate quantiles
quantile(beta_age, probs=c(0.025,0.5,0.975))
exp(mean(beta_age))
se(beta_age)

### Age^2
beta_age2 <- c(brod.out1[,20],brod.out2[,20],brod.out3[,20],brod.out4[,20],brod.out5[,20],
               brod.out6[,20],brod.out7[,20],brod.out8[,20],brod.out9[,20],brod.out10[,20])

# Calculate quantiles
quantile(beta_age2, probs=c(0.025,0.5,0.975))
mean(beta_age2)
se(beta_age2)

beta_age2 <- exp(beta_age2)
sum(beta_age2<1)/length(beta_age2)

### Buildings
beta_build <- c(brod.out1[,21],brod.out2[,21],brod.out3[,21],brod.out4[,21],brod.out5[,21],
                brod.out6[,21],brod.out7[,21],brod.out8[,21],brod.out9[,21],brod.out10[,21])

# Calculate quantiles
quantile(beta_build, probs=c(0.025,0.5,0.975))
mean(beta_build)
se(beta_build)

beta_build <- exp(beta_build)
sum(beta_build>1)/length(beta_build)

## Beech mast
beta_mast <- c(brod.out1[,22],brod.out2[,22],brod.out3[,22],brod.out4[,22],brod.out5[,22],
               brod.out6[,22],brod.out7[,22],brod.out8[,22],brod.out9[,22],brod.out10[,22])

# Calculate quantiles
quantile(beta_mast, probs=c(0.025,0.5,0.975))
mean(beta_mast)
se(beta_mast)

beta_mast <- exp(beta_mast)
sum(beta_mast>1)/length(beta_mast)

## Sex (male)
beta_sex2 <- c(brod.out1[,24],brod.out2[,24],brod.out3[,24],brod.out4[,24],brod.out5[,24],
               brod.out6[,24],brod.out7[,24],brod.out8[,24],brod.out9[,24],brod.out10[,24])

# Calculate quantiles
quantile(beta_sex2, probs=c(0.025,0.5,0.975))
mean(beta_sex2)
se(beta_sex2)

beta_sex2 <- exp(beta_sex2)
sum(beta_sex2>1)/length(beta_sex2)

## Stand age
beta_stand <- c(brod.out1[,25],brod.out2[,25],brod.out3[,25],brod.out4[,25],brod.out5[,25],
                brod.out6[,25],brod.out7[,25],brod.out8[,25],brod.out9[,25],brod.out10[,25])

# Calculate quantiles
quantile(beta_stand, probs=c(0.025,0.5,0.975))
mean(beta_stand)
se(beta_stand)

beta_stand <- exp(beta_stand)
sum(beta_stand>1)/length(beta_stand)


#### Bromadiolone ####
# Load posterior samples
brom.out1 <- readRDS("output/model_output/brom.out1.rds")
brom1.sum <- MCMCsummary(brom.out1) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=1)

brom.out2 <- readRDS("output/model_output/brom.out2.rds")
brom2.sum <- MCMCsummary(brom.out2) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=2)
brom.out3 <- readRDS("output/model_output/brom.out3.rds")
brom3.sum <- MCMCsummary(brom.out3) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=3)
brom.out4 <- readRDS("output/model_output/brom.out4.rds")
brom4.sum <- MCMCsummary(brom.out4) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=4)
brom.out5 <- readRDS("output/model_output/brom.out5.rds")
brom5.sum <- MCMCsummary(brom.out5) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=5)
brom.out6 <- readRDS("output/model_output/brom.out6.rds")
brom6.sum <- MCMCsummary(brom.out6) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=6)
brom.out7 <- readRDS("output/model_output/brom.out7.rds")
brom7.sum <- MCMCsummary(brom.out7) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=7)
brom.out8 <- readRDS("output/model_output/brom.out8.rds")
brom8.sum <- MCMCsummary(brom.out8) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=8)
brom.out9 <- readRDS("output/model_output/brom.out9.rds")
brom9.sum <- MCMCsummary(brom.out9) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=9)
brom.out10 <- readRDS("output/model_output/brom.out10.rds")
brom10.sum <- MCMCsummary(brom.out10) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=10)

brom.sum <- bind_rows(brom1.sum, brom2.sum, brom3.sum, brom4.sum, brom5.sum, 
                      brom6.sum, brom7.sum, brom8.sum, brom9.sum, brom10.sum)

brom.sum <- brom.sum %>% mutate(compound="bromacinone") %>%
  select(compound, iter, param:n.eff) %>%
  # exponentiate for odds ratios 
  mutate(mean=exp(mean),
         `2.5%`=exp(`2.5%`),
         `50%`=exp(`50%`),
         `97.5%`=exp(`97.5%`)) 

brom.odds <- brom.sum %>% group_by(param) %>% 
                 reframe(mean=mean(mean), `2.5%`=mean(`2.5%`), 
                         `50%`=mean(`50%`), `97.5%`=mean(`97.5%`)) %>% 
                 filter(str_detect(param, 'beta'))

## Combine all iterations
brom.out1 <- do.call("rbind",brom.out1)
brom.out2 <- do.call("rbind",brom.out2)
brom.out3 <- do.call("rbind",brom.out3)
brom.out4 <- do.call("rbind",brom.out4)
brom.out5 <- do.call("rbind",brom.out5)
brom.out6 <- do.call("rbind",brom.out6)
brom.out7 <- do.call("rbind",brom.out7)
brom.out8 <- do.call("rbind",brom.out8)
brom.out9 <- do.call("rbind",brom.out9)
brom.out10 <- do.call("rbind",brom.out10)

## Age
beta_age <- c(brom.out1[,19],brom.out2[,19],brom.out3[,19],brom.out4[,19],brom.out5[,19],
              brom.out6[,19],brom.out7[,19],brom.out8[,19],brom.out9[,19],brom.out10[,19])

# Calculate quantiles
quantile(beta_age, probs=c(0.025,0.5,0.975))
exp(mean(beta_age))
se(beta_age)

beta_age <- exp(beta_age)
sum(beta_age>1)/length(beta_age)

### Age^2
beta_age2 <- c(brom.out1[,20],brom.out2[,20],brom.out3[,20],brom.out4[,20],brom.out5[,20],
               brom.out6[,20],brom.out7[,20],brom.out8[,20],brom.out9[,20],brom.out10[,20])

# Calculate quantiles
quantile(beta_age2, probs=c(0.025,0.5,0.975))
mean(beta_age2)
se(beta_age2)

beta_age2 <- exp(beta_age2)
sum(beta_age2<1)/length(beta_age2)

### Building count
beta_build <- c(brom.out1[,21],brom.out2[,21],brom.out3[,21],brom.out4[,21],brom.out5[,21],
                brom.out6[,21],brom.out7[,21],brom.out8[,21],brom.out9[,21],brom.out10[,21])

# Calculate quantiles
quantile(beta_build, probs=c(0.025,0.5,0.975))
mean(beta_build)
se(beta_build)

beta_build <- exp(beta_build)
sum(beta_build>1)/length(beta_build)

### Beech mast
beta_mast <- c(brom.out1[,22],brom.out2[,22],brom.out3[,22],brom.out4[,22],brom.out5[,22],
               brom.out6[,22],brom.out7[,22],brom.out8[,22],brom.out9[,22],brom.out10[,22])

# Calculate quantiles
quantile(beta_mast, probs=c(0.025,0.5,0.975))
mean(beta_mast)
se(beta_mast)

beta_mast <- exp(beta_mast)
sum(beta_mast>1)/length(beta_mast)

### Sex (male)
beta_sex2 <- c(brom.out1[,24],brom.out2[,24],brom.out3[,24],brom.out4[,24],brom.out5[,24],
               brom.out6[,24],brom.out7[,24],brom.out8[,24],brom.out9[,24],brom.out10[,24])

# Calculate quantiles
quantile(beta_sex2, probs=c(0.025,0.5,0.975))
mean(beta_sex2)
se(beta_sex2)

beta_sex2 <- exp(beta_sex2)
sum(beta_sex2>1)/length(beta_sex2)

## Stand age
beta_stand <- c(brom.out1[,25],brom.out2[,25],brom.out3[,25],brom.out4[,25],brom.out5[,25],
                brom.out6[,25],brom.out7[,25],brom.out8[,25],brom.out9[,25],brom.out10[,25])

# Calculate quantiles
quantile(beta_stand, probs=c(0.025,0.5,0.975))
mean(beta_stand)
se(beta_stand)

beta_stand <- exp(beta_stand)
sum(beta_stand>1)/length(beta_stand)
