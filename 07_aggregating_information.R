library(tidyverse)

## Calculating number of individuals with 
# Individual compounds
# Landscape covariates
dat <- read_csv("output/AR_results_wide.csv")

max(dat$Diphacinone, na.rm=TRUE)
max(dat$Brodifacoum, na.rm=TRUE)
max(dat$Bromadiolone, na.rm=TRUE)

# # Reformat screening results
# dat[11:21] <- lapply(dat[11:21], function(x) replace(x, is.na(x), 0))
# dat[11:21] <- lapply(dat[11:21], function(x) replace(x, x>0, 1))
# 
# # sum columns to count how many had at least 1 exposure
# dat <- dat %>% mutate(any=Warfarin+Coumachlor+Diphacinone+Pindone+Brodifacoum+Difenacoum+
#                         Bromadiolone+Difethialone+Dicoumarol+Coumafuryl+Chlorophacinone)
# max(dat$any)
# dat %>% group_by(any) %>% count()
# 
# ggplot(dat, aes(x=any)) + geom_bar() + 
#   ylab("Number of individuals") + xlab("Number of compounds") +
#   theme_bw()
# 
# dat <- dat %>% mutate(any=if_else(any==0, 0, 1))
# 
# sum(dat$any)
# 
# dat %>% group_by(AgeClass,Sex) %>% count()
# 
# # make long format
# dat <- dat %>% pivot_longer(11:21, names_to="compound", values_to="exposed")
# 
# ## Sum values
# dat %>% group_by(compound) %>% summarize(n=sum(exposed))

# Reformat screening results
dat[11:21] <- lapply(dat[11:21], function(x) replace(x, x==1e-6, NA))

# pivot longer
# dat <- dat %>% pivot_longer(11:21, names_to="compound", values_to="ppm")
# dat %>% group_by(compound) %>% summarize(medianPPM=median(ppm, na.rm=TRUE))
# dat %>% group_by(compound) %>% summarize(meanPPM=mean(ppm, na.rm=TRUE))
# dat %>% group_by(compound) %>% summarize(rangePPM=range(ppm, na.rm=TRUE)) %>% print(n=22)


dat[11:21] <- lapply(dat[11:21], function(x) replace(x, x>1e-6, "measured"))
dat[11:21] <- lapply(dat[11:21], function(x) replace(x, x=="1e-06", "trace"))
dat[11:21] <- lapply(dat[11:21], function(x) replace(x, is.na(x), "not detect"))
dat <- dat %>% pivot_longer(11:21, names_to="compound", values_to="level")


dat %>% group_by(compound, level) %>% count() %>% print(n=24)




