library(tidyverse)
library(MCMCvis)

se <- function(x) {
  out <- sd(x)/sqrt(length((x)))
  return(out)
}

# Load posterior samples
ncomp.out1 <- readRDS("output/model_output/ncomp.out1.rds")
ncomp1.sum <- MCMCsummary(ncomp.out1) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=1)
ncomp.out2 <- readRDS("output/model_output/ncomp.out2.rds")
ncomp2.sum <- MCMCsummary(ncomp.out2) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=2)
ncomp.out3 <- readRDS("output/model_output/ncomp.out3.rds")
ncomp3.sum <- MCMCsummary(ncomp.out3) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=3)
ncomp.out4 <- readRDS("output/model_output/ncomp.out4.rds")
ncomp4.sum <- MCMCsummary(ncomp.out4) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=4)
ncomp.out5 <- readRDS("output/model_output/ncomp.out5.rds")
ncomp5.sum <- MCMCsummary(ncomp.out5) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=5)
ncomp.out6 <- readRDS("output/model_output/ncomp.out6.rds")
ncomp6.sum <- MCMCsummary(ncomp.out6) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=6)
ncomp.out7 <- readRDS("output/model_output/ncomp.out7.rds")
ncomp7.sum <- MCMCsummary(ncomp.out7) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=7)
ncomp.out8 <- readRDS("output/model_output/ncomp.out8.rds")
ncomp8.sum <- MCMCsummary(ncomp.out8) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=8)
ncomp.out9 <- readRDS("output/model_output/ncomp.out9.rds")
ncomp9.sum <- MCMCsummary(ncomp.out9) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=9)
ncomp.out10 <- readRDS("output/model_output/ncomp.out10.rds")
ncomp10.sum <- MCMCsummary(ncomp.out10) %>% rownames_to_column("param") %>% as_tibble() %>%
  mutate(iter=10)

ncomp.sum <- bind_rows(ncomp1.sum, ncomp2.sum, ncomp3.sum,  ncomp5.sum, #ncomp4.sum,
                      ncomp6.sum, ncomp7.sum, ncomp8.sum, ncomp9.sum, ncomp10.sum)

ncomp.sum <- ncomp.sum %>% select(iter, param:n.eff) %>%
  group_by(param) %>% reframe(mean=mean(mean), 
                              `2.5%`=mean(`2.5%`), 
                              # `50%`=mean(`50%`), 
                              `97.5%`=mean(`97.5%`)) %>% 
  filter(str_detect(param, 'beta'))

## Combine all iterations
ncomp.out1 <- do.call("rbind",ncomp.out1)
ncomp.out2 <- do.call("rbind",ncomp.out2)
ncomp.out3 <- do.call("rbind",ncomp.out3)
ncomp.out4 <- do.call("rbind",ncomp.out4)
ncomp.out5 <- do.call("rbind",ncomp.out5)
ncomp.out6 <- do.call("rbind",ncomp.out6)
ncomp.out7 <- do.call("rbind",ncomp.out7)
ncomp.out8 <- do.call("rbind",ncomp.out8)
ncomp.out9 <- do.call("rbind",ncomp.out9)
ncomp.out10 <- do.call("rbind",ncomp.out10)

## Age
beta_age <- c(ncomp.out1[,19],ncomp.out2[,19],ncomp.out3[,19],ncomp.out5[,19],ncomp.out4[,19],
              ncomp.out6[,19],ncomp.out7[,19],ncomp.out8[,19],ncomp.out9[,19],ncomp.out10[,19])

quantile(beta_age, probs=c(0.025,0.5,0.975))
exp(mean(beta_age))
se(beta_age)
sum(beta_age>0)/length(beta_age)

## Age^2
beta_age2 <- c(ncomp.out1[,20],ncomp.out2[,20],ncomp.out3[,20],ncomp.out5[,20],ncomp.out4[,20],
               ncomp.out6[,20],ncomp.out7[,20],ncomp.out8[,20],ncomp.out9[,20],ncomp.out10[,20])

quantile(beta_age2, probs=c(0.025,0.5,0.975))
mean(beta_age2)
se(beta_age2)
sum(beta_age2<0)/length(beta_age2)

# Intermix * mast interaction
beta_intx <- c(ncomp.out1[,21],ncomp.out2[,21],ncomp.out3[,21],ncomp.out5[,21],ncomp.out4[,21],
               ncomp.out6[,21],ncomp.out7[,21],ncomp.out8[,21],ncomp.out9[,21],ncomp.out10[,21])

quantile(beta_intx, probs=c(0.025,0.5,0.975))
mean(exp(beta_intx))
se(beta_intx)
sum(beta_intx>0)/length(beta_intx)

## Beech mast
beta_mast <- c(ncomp.out1[,22],ncomp.out2[,22],ncomp.out3[,22],ncomp.out5[,22],ncomp.out4[,22],
               ncomp.out6[,22],ncomp.out7[,22],ncomp.out8[,22],ncomp.out9[,22],ncomp.out10[,22])

quantile(beta_mast, probs=c(0.025,0.5,0.975))
mean(beta_mast)
se(beta_mast)
sum(beta_mast>0)/length(beta_mast)

## Sex (M)
beta_sex2 <- c(ncomp.out1[,24],ncomp.out2[,24],ncomp.out3[,24],ncomp.out5[,24],ncomp.out4[,24],
               ncomp.out6[,24],ncomp.out7[,24],ncomp.out8[,24],ncomp.out9[,24],ncomp.out10[,24])

quantile(beta_sex2, probs=c(0.025,0.5,0.975))
mean(beta_sex2)
se(beta_sex2)
sum(beta_sex2>0)/length(beta_sex2)

## WUI intermix
beta_wui <- c(ncomp.out1[,25],ncomp.out2[,25],ncomp.out3[,25],ncomp.out5[,25],ncomp.out4[,25],
              ncomp.out6[,25],ncomp.out7[,25],ncomp.out8[,25],ncomp.out9[,25],ncomp.out10[,25])

quantile(beta_wui, probs=c(0.025,0.5,0.975))
mean(beta_wui)
se(beta_wui)
sum(beta_wui>0)/length(beta_wui)

