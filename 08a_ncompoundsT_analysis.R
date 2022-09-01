library(tidyverse)
library(MuMIn)
library(ordinal)
library(parallel)

set.seed(123)

#### Set up parallel processing ####

# Determine cluster type (mine is a PSOCK)
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"

# Detect number of cores and create cluster (leave a couple out to not overwhelm pc)
nCores <- detectCores() - 3
cl <- makeCluster(nCores, type=clusterType)

clusterEvalQ(cl,library(ordinal))
clusterEvalQ(cl,library(MuMIn))

#### Read in data ####
dat <- read_csv("data/analysis-ready/combined_AR_covars.csv")
dat <- as.data.frame(dat)

#### Analysis Set-up ####

# Categorical number of compounds (collapsing higher numbers)
dat$catNcompT <- ifelse(dat$n.compounds.T>=3, "3+", as.character(dat$n.compounds.T))
dat$catNcompT  <- ordered(dat$catNcompT , levels=c("0", "1", "2", "3+"))

# Order response
dat$n.compounds.T <- ordered(dat$n.compounds.T, levels=c(0,1,2,3,4,5))

# Make random effects factors
dat$WMU <- as.factor(dat$WMU)
dat$WMUA_code <- as.factor(dat$WMUA_code)

# Change how beech mast is incorporated (lagged values)
dat$mast <- NA 
dat$mast[dat$year==2018] <- "mast"
dat$mast[dat$year==2019] <- "fail"
dat$mast[dat$year==2020] <- "mast"
dat$mast <- as.factor(dat$mast)

# Make age a categorical variable
dat$catAge[dat$Age>=3.5] <- "adult"
dat$catAge[dat$Age==2.5] <- "subadult"
dat$catAge[dat$Ag<2.5] <- "juvenile"

dat$year <- factor(dat$year)

# Resort columns
dat <- dat[,c(1:7,31,10:13,8,9,34,15,32,16:17,29,33,18:28)] 

## Scale and center variables
dat[,c(14,20,22:32)] <- scale(dat[,c(14,20,22:32)])

## Percent AG
pctAG1 <- dat[, c(1:18, 22:24)]
pctAG1 <- distinct(pctAG1)
pctAG1 <- pctAG1 %>% group_by(RegionalID) %>% 
  pivot_wider(names_from=buffsize, values_from=c(pasture, crops, totalag)) %>% as.data.frame()

# Run models
ag15 <- clmm(catNcompT ~ totalag_15 + (1|WMUA_code/WMU) + (1|year), data=pctAG1)
ag15sq <- clmm(catNcompT ~ totalag_15 + I(totalag_15^2) + (1|WMUA_code/WMU) + (1|year), data=pctAG1)
ag30 <- clmm(catNcompT ~ totalag_30 + (1|WMUA_code/WMU) + (1|year), data=pctAG1)
ag30sq <- clmm(catNcompT ~ totalag_30 + I(totalag_30^2) +(1|WMUA_code/WMU) + (1|year), data=pctAG1)
ag60 <- clmm(catNcompT ~ totalag_60 + (1|WMUA_code/WMU) + (1|year), data=pctAG1)
ag60sq <- clmm(catNcompT ~ totalag_60 + I(totalag_60^2) + (1|WMUA_code/WMU) + (1|year), data=pctAG1)

crop15 <- clmm(catNcompT ~ crops_15 + (1|WMUA_code/WMU) + (1|year), data=pctAG1)
crop15sq <- clmm(catNcompT ~ crops_15 + I(crops_15^2) + (1|WMUA_code/WMU) + (1|year), data=pctAG1)
crop30 <- clmm(catNcompT ~ crops_30 + (1|WMUA_code/WMU) + (1|year), data=pctAG1)
crop30sq <- clmm(catNcompT ~ crops_30 + I(crops_30^2) +(1|WMUA_code/WMU) + (1|year), data=pctAG1)
crop60 <- clmm(catNcompT ~ crops_60 + (1|WMUA_code/WMU) + (1|year), data=pctAG1)
crop60sq <- clmm(catNcompT ~ crops_60 + I(crops_60^2) + (1|WMUA_code/WMU) + (1|year), data=pctAG1)

past15 <- clmm(catNcompT ~ pasture_15 + (1|WMUA_code/WMU) + (1|year), data=pctAG1)
past15sq <- clmm(catNcompT ~ pasture_15 + I(pasture_15^2) + (1|WMUA_code/WMU) + (1|year), data=pctAG1)
past30 <- clmm(catNcompT ~ pasture_30 + (1|WMUA_code/WMU) + (1|year), data=pctAG1)
past30sq <- clmm(catNcompT ~ pasture_30 + I(pasture_30^2) +(1|WMUA_code/WMU) + (1|year), data=pctAG1)
past60 <- clmm(catNcompT ~ pasture_60 + (1|WMUA_code/WMU) + (1|year), data=pctAG1)
past60sq <- clmm(catNcompT ~ pasture_60 + I(pasture_60^2) + (1|WMUA_code/WMU) + (1|year), data=pctAG1)

pctAG_sel <- model.sel(ag15, ag30, ag60, ag15sq, ag30sq, ag60sq,
                       crop15, crop30, crop60, crop15sq, crop30sq, crop60sq,
                       past15, past30, past60, past15sq, past30sq, past60sq)
pctAG_sel

## Beech basal area
baa1 <- dat[, c(1:18,20,21,29)]
baa1  <- baa1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffsize, values_from=baa, values_fn=unique) %>% as.data.frame()
names(baa1)[20:22] <- c("baa_15", "baa_30", "baa_60") 

baa15 <- clmm(catNcompT ~ baa_15 + (1|WMUA_code/WMU) + (1|year), data=baa1)
baa15sq <- clmm(catNcompT ~ baa_15 + I(baa_15^2) + (1|WMUA_code/WMU) + (1|year), data=baa1)
baa30 <- clmm(catNcompT ~ baa_30 + (1|WMUA_code/WMU) + (1|year), data=baa1)
baa30sq <- clmm(catNcompT ~ baa_30 + I(baa_30^2) + (1|WMUA_code/WMU) + (1|year), data=baa1)
baa60 <- clmm(catNcompT ~ baa_60 + (1|WMUA_code/WMU) + (1|year), data=baa1)
baa60sq <- clmm(catNcompT ~ baa_60 + I(baa_60^2) + (1|WMUA_code/WMU) + (1|year), data=baa1)

baa15lBMI <- clmm(catNcompT ~ baa_15*laggedBMI + (1|WMUA_code/WMU) + (1|year), data=baa1)
baa30lBMI <- clmm(catNcompT ~ baa_30*laggedBMI + (1|WMUA_code/WMU) + (1|year), data=baa1)
baa60lBMI <- clmm(catNcompT ~ baa_60*laggedBMI + (1|WMUA_code/WMU) + (1|year), data=baa1)

baa15M <- clmm(catNcompT ~ baa_15*mast + (1|WMUA_code/WMU) + (1|year), data=baa1)
baa30M <- clmm(catNcompT ~ baa_30*mast + (1|WMUA_code/WMU) + (1|year), data=baa1)
baa60M <- clmm(catNcompT ~ baa_60*mast + (1|WMUA_code/WMU) + (1|year), data=baa1)

baa_sel <- model.sel(baa15, baa30, baa60, baa15sq, baa30sq, baa60sq, 
                     baa15lBMI, baa30lBMI, baa60lBMI, 
                     baa15M, baa30M, baa60M)
baa_sel

## Wildland-urban interface
intermix1 <- dat[, c(1:19, 30)]
intermix1 <- unite(intermix1, "buffrad", 18:19, sep="_")
intermix1  <- intermix1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=intermix) %>% as.data.frame()
names(intermix1)[18:26] <- c("mix_15_100", "mix_30_100", "mix_60_100",
                             "mix_15_250", "mix_30_250", "mix_60_250",
                             "mix_15_500", "mix_30_500", "mix_60_500") 

interface1 <- dat[, c(1:19, 31)]
interface1 <- unite(interface1, "buffrad", 18:19, sep="_")
interface1  <- interface1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=interface) %>% as.data.frame()
names(interface1)[18:26] <- c("face_15_100", "face_30_100", "face_60_100",
                              "face_15_250", "face_30_250", "face_60_250",
                              "face_15_500", "face_30_500", "face_60_500") 

wui1 <- dat[, c(1:19, 32)]
wui1 <- unite(wui1, "buffrad", 18:19, sep="_")
wui1  <- wui1  %>% group_by(RegionalID) %>% pivot_wider(names_from=buffrad, values_from=totalWUI) %>% as.data.frame()
names(wui1)[18:26] <- c("wui_15_100", "wui_30k_100", "wui_60_100",
                        "wui_15_250", "wui_30k_250", "wui_60_250",
                        "wui_15_500", "wui_30k_500", "wui_60_500") 

# Join intermix WUI
intermix1 <- intermix1[,c(1:3, 18:26)]
wui1 <- left_join(wui1, intermix1, by=c("RegionalID", "pt_name", "pt_index"))

# Join interface WUI
interface1 <- interface1[, c(1:3, 18:26)]
wui1 <- left_join(wui1, interface1, by=c("RegionalID", "pt_name", "pt_index"))

# Intermix WUI
mix_15100 <- clmm(catNcompT ~ mix_15_100 + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_15250 <- clmm(catNcompT ~ mix_15_250 + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_15500 <- clmm(catNcompT ~ mix_15_500 + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_15100sq <- clmm(catNcompT ~ mix_15_100 + I(mix_15_100^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_15250sq <- clmm(catNcompT ~ mix_15_250 + I(mix_15_250^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_15500sq <- clmm(catNcompT ~ mix_15_500 + I(mix_15_500^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)

mix_30100 <- clmm(catNcompT ~ mix_30_100 + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_30250 <- clmm(catNcompT ~ mix_30_250 + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_30500 <- clmm(catNcompT ~ mix_30_500 + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_30100sq <- clmm(catNcompT ~ mix_30_100 + I(mix_30_100^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_30250sq <- clmm(catNcompT ~ mix_30_250 + I(mix_30_250^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_30500sq <- clmm(catNcompT ~ mix_30_500 + I(mix_30_500^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)

mix_60100 <- clmm(catNcompT ~ mix_60_100 + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_60250 <- clmm(catNcompT ~ mix_60_250 + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_60500 <- clmm(catNcompT ~ mix_60_500 + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_60100sq <- clmm(catNcompT ~ mix_60_100 + I(mix_60_100^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_60250sq <- clmm(catNcompT ~ mix_60_250 + I(mix_60_250^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
mix_60500sq <- clmm(catNcompT ~ mix_60_500 + I(mix_60_500^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)

# Interface WUI
face_15100 <- clmm(catNcompT ~ face_15_100 + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_15250 <- clmm(catNcompT ~ face_15_250 + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_15500 <- clmm(catNcompT ~ face_15_500 + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_15100sq <- clmm(catNcompT ~ face_15_100 + I(face_15_100^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_15250sq <- clmm(catNcompT ~ face_15_250 + I(face_15_250^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_15500sq <- clmm(catNcompT ~ face_15_500 + I(face_15_500^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)

face_30100 <- clmm(catNcompT ~ face_30_100 + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_30250 <- clmm(catNcompT ~ face_30_250 + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_30500 <- clmm(catNcompT ~ face_30_500 + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_30100sq <- clmm(catNcompT ~ face_30_100 + I(face_30_100^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_30250sq <- clmm(catNcompT ~ face_30_250 + I(face_30_250^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_30500sq <- clmm(catNcompT ~ face_30_500 + I(face_30_500^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)

face_60100 <- clmm(catNcompT ~ face_60_100 + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_60250 <- clmm(catNcompT ~ face_60_250 + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_60500 <- clmm(catNcompT ~ face_60_500 + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_60100sq <- clmm(catNcompT ~ face_60_100 + I(face_60_100^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_60250sq <- clmm(catNcompT ~ face_60_250 + I(face_60_250^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
face_60500sq <- clmm(catNcompT ~ face_60_500 + I(face_60_500^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)

# Total WUI
wui_15100 <- clmm(catNcompT ~ wui_15_100 + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_15250 <- clmm(catNcompT ~ wui_15_250 + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_15500 <- clmm(catNcompT ~ wui_15_500 + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_15100sq <- clmm(catNcompT ~ wui_15_100 + I(wui_15_100^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_15250sq <- clmm(catNcompT ~ wui_15_250 + I(wui_15_250^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_15500sq <- clmm(catNcompT ~ wui_15_500 + I(wui_15_500^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)

wui_30100 <- clmm(catNcompT ~ wui_30k_100 + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_30250 <- clmm(catNcompT ~ wui_30k_250 + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_30500 <- clmm(catNcompT ~ wui_30k_500 + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_30100sq <- clmm(catNcompT ~ wui_30k_100 + I(wui_30k_100^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_30250sq <- clmm(catNcompT ~ wui_30k_250 + I(wui_30k_250^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_30500sq <- clmm(catNcompT ~ wui_30k_500 + I(wui_30k_500^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)

wui_60100 <- clmm(catNcompT ~ wui_60_100 + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_60250 <- clmm(catNcompT ~ wui_60_250 + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_60500 <- clmm(catNcompT ~ wui_60_500 + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_60100sq <- clmm(catNcompT ~ wui_60_100 + I(wui_60_100^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_60250sq <- clmm(catNcompT ~ wui_60_250 + I(wui_60_250^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)
wui_60500sq <- clmm(catNcompT ~ wui_60_500 + I(wui_60_500^2) + (1|WMUA_code/WMU) + (1|year), data=wui1)

wui_sel <- model.sel(wui_15100, wui_30100, wui_60100, wui_15100sq, wui_30100sq, wui_60100sq, 
                     wui_15250, wui_30250, wui_60250, wui_15250sq, wui_30250sq, wui_60250sq, 
                     wui_15500, wui_30500, wui_60500, wui_15500sq, wui_30500sq, wui_60500sq,
                     face_15100, face_30100, face_60100, face_15100sq, face_30100sq, face_60100sq,  
                     face_15250, face_30250, face_60250, face_15250sq, face_30250sq, face_60250sq,  
                     face_15500, face_30500, face_60500, face_15500sq, face_30500sq, face_60500sq,
                     mix_15100, mix_30100, mix_60100, mix_15100sq, mix_30100sq, mix_60100sq,
                     mix_15250, mix_30250, mix_60250, mix_15250sq, mix_30250sq, mix_60250sq,
                     mix_15500, mix_30500, mix_60500, mix_15500sq, mix_30500sq, mix_60500sq)
wui_sel

#### Set up data to run for each combination of covariates ####
dat1 <- dat[,c(1:17)]
dat1 <- distinct(dat1)

# Join percent agriculture
pctAG1 <- pctAG1[,c(1:3, 20)]
dat1 <- left_join(dat1, pctAG1, by=c("RegionalID", "pt_name", "pt_index"))

# Join beech basal area
baa1 <- baa1[,c(1:3, 22, 18)]
dat1 <- left_join(dat1, baa1, by=c("RegionalID", "pt_name", "pt_index"))

# Join total WUI
wui1 <- wui1[,c(1:3, 35)]
dat1 <- left_join(dat1, wui1, by=c("RegionalID", "pt_name", "pt_index"))

## Check correlation matrix
cor(dat1[,18:21])

## Set up global models
g1 <- clmm(catNcompT ~ Sex*catAge + 
                pasture_60 + I(pasture_60^2) + 
                mix_60_500 + I(mix_60_500 ^2) +
                baa_60 + laggedBMI + baa_60:laggedBMI + 
                (1|WMUA_code/WMU) + (1|year), data=dat1, na.action="na.fail")

# Export data and model into the cluster worker nodes
clusterExport(cl, c("dat1","g1"))

g1_dredge <- MuMIn:::.dredge.par(g1, cluster=cl, trace=2, rank="AIC", 
                                 subset=dc(mix_60_500, I(mix_60_500^2)) &&
                                        dc(pasture_60, I(pasture_60^2)) &&
                                        dc(baa_60, laggedBMI, baa_60:laggedBMI))
  
# Save dredge tables
save(file="output/dredge_tables_ncompT.Rdata", list="g1_dredge")

# Save as csv file
dh <- as.data.frame(g1_dredge)
write_csv(dh, "output/ncompT_dredgetable.csv")

## Read tables back in
dh <- read_csv("output/ncompT_dredgetable.csv")
dh <- as.data.frame(dh)

#### Running final models ####

## Loop over each set of random points

m_est <- data.frame()
m_stderr <- data.frame()
pct2.5 <- data.frame()
pct97.5 <- data.frame()

# Loop over each point set
for (i in 1:10) {
  
  # Subset to one point set at a time
  pt <- dat1[dat1$pt_index==i,]
  
  # Run model with deltaAICc < 2
  m1_pt <- clmm(catNcompT ~ Sex*catAge + pasture_60 + 
                  mix_60_500 + I(mix_60_500^2) + 
                  (1|WMUA_code/WMU) + (1|year), data=pt, na.action="na.fail")
  
  m1s <- summary(m1_pt)
  
  # save averaged confidence intervals
  pct2.5 <- rbind(pct2.5, t(confint(m1_pt))[1,])
  pct97.5 <- rbind(pct97.5, t(confint(m1_pt))[2,])

  # Save point set estimates
  m_est <- rbind(m_est, coef(m1s)[,1])
  m_stderr <- rbind(m_stderr, coef(m1s)[,2])
  
  # Rename (only need to do once)
  if (i==1) {
    names(m_est) <- c(names(coef(m1_pt)))
    names(m_stderr) <- c(names(coef(m1_pt)))
    names(pct2.5) <- c(names(coef(m1_pt)))
    names(pct97.5) <- c(names(coef(m1_pt)))
  }
  
}

# Calculate averages for each coefficient
coef_avg <- colMeans(m_est[sapply(m_est, is.numeric)], na.rm=TRUE)
stderr_avg <- colMeans(m_stderr[sapply(m_stderr, is.numeric)], na.rm=TRUE)
pct2.5_avg <- colMeans(pct2.5[sapply(pct2.5, is.numeric)], na.rm=TRUE)
pct97.5_avg <- colMeans(pct97.5[sapply(pct97.5, is.numeric)], na.rm=TRUE)

# One-sample t-test to determine "significance"
pvalue <- c()
for (i in 1:ncol(m_est)) {
  
  tresult <- t.test(m_est[,i], mu=0, alternative="two.sided")
  pvalue <- c(pvalue, tresult$p.value)
}
pvalue <- as.data.frame(pvalue)
pvalue <- cbind(names(coef_avg), pvalue)
pvalue <- pivot_wider(pvalue, names_from="names(coef_avg)", values_from="pvalue") %>% as.data.frame()

# Combine and clean up data frame
coef_summary <- bind_rows(coef_avg, stderr_avg, pct2.5_avg, pct97.5_avg, pvalue)
coef_summary <- as.data.frame(coef_summary)
coefs <- c("param_est", "std_error", "2.5CI", "97.5CI", "P-value")
coef_summary <- data.frame(coef=coefs, coef_summary)

# Write to file
write_csv(coef_summary, "results/ncompT_coef-summary.csv")











