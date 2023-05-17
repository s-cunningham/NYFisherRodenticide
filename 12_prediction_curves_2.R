library(tidyverse)
library(glmmTMB)
library(sjPlot)
library(sjmisc)
library(ggeffects)

options(scipen=999, digits=3)
set.seed(123)
theme_set(theme_classic())

#### Parallel processing ####
nt <- min(parallel::detectCores(),4)

## Read data
dat1 <- read_csv("output/model_data_notscaled.csv")

# Scale and center variables
dat1[,c(8,16:83)] <- scale(dat1[,c(8,16:83)])

dfp_mix <- list()
dfp_lbn <- list()
dfp_past <- list()
dfp_as <- list()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- dat1[dat1$pt_index==i,]
  
  # Run model 
  m1_pt <- glmmTMB(n.compounds.T ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * lag_beechnuts + (1|WMU), data=pt, 
                   family=compois(link = "log"), control=glmmTMBControl(parallel=nt))
  

  p_mix <- plot_model(m1_pt, type="pred", terms=c("mix_15_100 [-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4]", 
                                                  "Sex [all]"))
  dfp_mix[[i]] <- as.data.frame(p_mix$data)
  
  p_past <- plot_model(m1_pt, type="pred", terms=c("pasture_15 [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3]",
                                                   "Sex [all]" ))
  dfp_past[[i]] <- as.data.frame(p_past$data)
  
  
  p_lbn <- plot_model(m1_pt, type="pred", terms=c("lag_beechnuts [-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]", 
                                                  "Sex [all]"), pred.type="fe")
  dfp_lbn[[i]] <- as.data.frame(p_lbn$data)
  

  p_as <- plot_model(m1_pt, type="pred", terms=c("Age [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]", 
                                                 "Sex [all]"), pred.type="fe")
  dfp_as[[i]] <- as.data.frame(p_as$data)
  
}



