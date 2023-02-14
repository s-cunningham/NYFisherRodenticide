library(tidyverse)
library(glmmTMB)
library(MuMIn)

nt <- min(parallel::detectCores(),4)

# read age data
age <- read_csv("data/analysis-ready/2016-2020_ages_data.csv") %>% 
          filter(HarvestYear>=2018) %>%
          select(RegionalID, HarvestYear:Age, Region:Town) %>%
          rename(Year=HarvestYear)

# read CL data
cl <- read_csv("data/analysis-ready/CL_data_20230111.csv") %>% select(RegionalID,`Active CL`) %>%
        rename(ActiveCL=`Active CL`) 
        
cl$ActiveCL[cl$ActiveCL=="ND" | 
            cl$ActiveCL=="NP" | 
            cl$ActiveCL=="autolyed" | 
            cl$ActiveCL=="AUTOLYZED" | 
            cl$ActiveCL=="Autolyzed" | 
            cl$ActiveCL=="I havethis as 2019"] <- NA

dat <- left_join(age, cl, by=c("RegionalID")) %>% filter(!is.na(ActiveCL) & Age>0.5) 
dat$ActiveCL <- as.numeric(dat$ActiveCL)


dat <- dat %>% mutate(Section=case_when(
                    Region==5 ~ "Adirondack",
                    Region==6 ~ "Tug Hill",
                    Region<=4 ~ "Catskills",
                    Region>=7 ~ "southern Tier"))
dat %>% group_by(Section) %>% count()

ggplot(dat) + geom_boxplot(aes(x=Section, y=ActiveCL)) 

# Join CL to rodenticide
dets <- read_csv("output/summarized_AR_results.csv") %>% 
            select(RegionalID:compound,bin.exp) %>%
            pivot_wider(names_from=compound, values_from=bin.exp) %>%
            mutate(FGAR=Warfarin+Diphacinone+Chlorophacinone+Dicoumarol,
                   SGAR=Brodifacoum+Difenacoum+Bromadiolone+Difethialone) %>%
            mutate(ARgroup=case_when(
                   FGAR==0 & SGAR==0 ~ "none",
                   FGAR>=1 & SGAR>=1 ~ "both",
                   FGAR>=1 & SGAR==0 ~ "FGAR",
                   FGAR==0 & SGAR>=1 ~ "SGAR"
                   )) %>%
            select(RegionalID:WMU,Town,ARgroup)

ar <- read_csv("output/ncompounds_trace.csv") %>% 
  select(RegionalID, n.compounds) %>%
  rename(n.compounds=n.compounds)

ar <- left_join(ar, cl, by="RegionalID") 
ar <- left_join(dets, ar, by="RegionalID") %>% filter(!is.na(ActiveCL) & Age>0.5)
ar$Sex[ar$Sex=="M"] <- "F" 
ar$ActiveCL <- as.numeric(ar$ActiveCL)

ar <- ar %>% mutate(Section=case_when(
                    Region==5 ~ "Adirondack",
                    Region==6 ~ "Tug Hill",
                    Region<=4 ~ "Catskills",
                    Region>=7 ~ "Southern Tier"))
ar %>% group_by(Section) %>% summarize(median(ActiveCL))


m1 <- glmmTMB(ActiveCL ~ scale(n.compounds) + (1|WMU), data=ar, family=compois(link = "log"), control=glmmTMBControl(parallel=nt))
m2 <- glmmTMB(ActiveCL ~ scale(n.compounds) + scale(Age) + (1|WMU), data=ar, family=compois(link="log"), control=glmmTMBControl(parallel=nt))
m3 <- glmmTMB(ActiveCL ~ scale(n.compounds) + factor(year) + (1|WMU), data=ar, family=compois(link="log"), control=glmmTMBControl(parallel=nt))
m4 <- glmmTMB(ActiveCL ~ scale(n.compounds) + factor(Section) + (1|WMU), data=ar, family=compois(link="log"), control=glmmTMBControl(parallel=nt))
m5 <- glmmTMB(ActiveCL ~ scale(n.compounds) + scale(Age) + factor(Section) + (1|WMU), data=ar, family=compois(link="log"), control=glmmTMBControl(parallel=nt))
m6 <- glmmTMB(ActiveCL ~ scale(n.compounds) + scale(Age)*factor(Section) + (1|WMU), data=ar, family=compois(link="log"), control=glmmTMBControl(parallel=nt))
m7 <- glmmTMB(ActiveCL ~ scale(n.compounds) * scale(Age) + (1|WMU), data=ar, family=compois(link="log"), control=glmmTMBControl(parallel=nt))
# m7 <- glmmTMB(ActiveCL ~ factor(ARgroup) + (1|WMU), data=ar, family=compois(link="log"), control=glmmTMBControl(parallel=nt))
# m8 <- glmmTMB(ActiveCL ~ factor(ARgroup) + scale(Age) + (1|WMU), data=ar, family=compois(link="log"), control=glmmTMBControl(parallel=nt))
# m9 <- glmmTMB(ActiveCL ~ factor(ARgroup) + factor(Section) + (1|WMU), data=ar, family=compois(link="log"), control=glmmTMBControl(parallel=nt))
# m10 <- glmmTMB(ActiveCL ~ factor(ARgroup) + scale(Age) + factor(Section) + (1|WMU), data=ar, family=compois(link="log"), control=glmmTMBControl(parallel=nt))
# m11 <- glmmTMB(ActiveCL ~ factor(ARgroup) + scale(Age)*factor(Section) + (1|WMU), data=ar, family=compois(link="log"), control=glmmTMBControl(parallel=nt))
# m12 <- glmmTMB(ActiveCL ~ factor(ARgroup) + factor(Section) + (1|WMU), data=ar, family=compois(link="log"), control=glmmTMBControl(parallel=nt))

model.sel(m1,m2,m3,m4,m5,m6,m7)

summary(m4)

new.dat <- data.frame(diameter= d$diameter,
                      plant_density = d$plant_density,
                      plot= d$plot) 

plot(ggpredict(m4))

plot(ggpredict(m1))


new.dat$prediction <- predict(glmm.model, new.data = new.dat, 
                              type = "response", re.form = NA)


# ,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17
# factor(ARgroup) + factor(year) + factor(Section)
# factor(ARgroup) + factor(Section)
# factor(ARgroup) 
# factor(ARgroup) + factor(year)
# factor(ARgroup) + factor(year)
# factor(ARgroup) + scale(Age) + factor(Section)
