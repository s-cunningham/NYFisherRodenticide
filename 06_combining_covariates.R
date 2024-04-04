## Set up rodenticide covariates
## 2022-06-27, updated 2023-01-12

library(tidyverse)

#### Read in data ####
dat <- read_csv("output/summarized_AR_results.csv")
trace_y <- read_csv("output/ncompounds_trace.csv") %>% 
  select(RegionalID, n.compounds) %>%
  rename(n.compounds.T=n.compounds)
trace_n <- read_csv("output/ncompounds_notrace.csv") %>% 
  select(RegionalID, n.compounds) %>%
  rename(n.compounds.MO=n.compounds)
nlcd <- read_csv("data/analysis-ready/nlcd_pct.csv")
bmi <- read_csv("data/analysis-ready/baa_sum.csv")
baa <- read_csv("data/analysis-ready/baa_sum_single_raster.csv")
pts <- read_csv("output/random_point_locs.csv")
wmua <- read_csv("data/analysis-ready/wmuas.csv")
build <- read_csv("data/analysis-ready/building-centroid_sum.csv") %>%
            rename(pt_name=name) 
mast <- read_csv("data/analysis-ready/ALTEMP26_beech-data.csv")

stand_mn <- read_csv("data/analysis-ready/stand-age_mean.csv") %>%
                  rename(pt_name=name) 

stand_sd <- read_csv("data/analysis-ready/stand-age_stdev.csv") %>%
  rename(pt_name=name) 

lsm <- read_csv("data/analysis-ready/forest_edge_density.csv") %>%
            mutate(buffsize=case_when(buffer==2185.0969 ~ 15,
                                      buffer==3090.1936 ~ 30,
                                      buffer==3784.699 ~ 45)) %>%
            rename(edge_density=value, pt_name=plot_id) %>%
            select(pt_name, buffsize, edge_density) %>%
            mutate(buffsize=coalesce(buffsize, 45))

#### Combine data ####
## Number of compounds detected
# ncomps <- left_join(trace_n, trace_y, by="RegionalID")

## Save details of each sample
dets <- dat %>% select(RegionalID:WMU,Town) %>% distinct()

## Add a column to points with just sample ID
pts <- pts %>% select(RegionalID,name,x,y) %>%
          rename(rand_x=x, rand_y=y, pt_name=name)

## Subset forest and add together
forest <- nlcd %>% filter(value==41 | value==42 | value==43)   
forest$value <- factor(forest$value, levels=c(41, 42, 43), labels=c("deciduous", "evergreen", "mixed")) |> as.character()
forest <- forest %>% select(name,buffsize,value,freq) %>%
           pivot_wider(names_from=value, values_from=freq) %>%
           rename(pt_name=name)
forest$deciduous[is.na(forest$deciduous)] <- 0
forest$evergreen[is.na(forest$evergreen)] <- 0
forest$mixed[is.na(forest$mixed)] <- 0
forest <- forest %>% mutate(totalforest=deciduous + evergreen + mixed) 
  
## reorganize beech
# beech mast index
names(bmi)[c(1,3)] <- c("pt_name", "bmi")
bmi$bmi[is.na(bmi$bmi)] <- 0

# beech basal area
names(baa)[1] <- "pt_name"
baa$baa[is.na(baa$baa)] <- 0
baa <- baa %>% select(pt_name, buffsize, baa) %>% rename(bba=baa)

# Create categorical variable for buildings
qntl <- build %>% group_by(buffsize) %>% summarize(first=quantile(nbuildings, probs=0.25), median=quantile(nbuildings, probs=0.5),
                                                   third=quantile(nbuildings, probs=0.75), fourth=quantile(nbuildings, probs=1))

build <- build %>% mutate(build_cat=case_when(
                             nbuildings<1 ~ "None",
                             (nbuildings>=1 & nbuildings<=as.vector(qntl[1,2])) & buffsize==15 ~ "1stQuart",
                             (nbuildings>=1 & nbuildings<=as.vector(qntl[2,2])) & buffsize==30 ~ "1stQuart",
                             (nbuildings>=1 & nbuildings<=as.vector(qntl[3,2])) & buffsize==45 ~ "1stQuart",
                             (nbuildings>as.vector(qntl[1,2]) & nbuildings<=as.vector(qntl[1,3])) & buffsize==15 ~ "2ndQuart",
                             (nbuildings>as.vector(qntl[2,2]) & nbuildings<=as.vector(qntl[2,3])) & buffsize==30 ~ "2ndQuart",
                             (nbuildings>as.vector(qntl[3,2]) & nbuildings<=as.vector(qntl[3,3])) & buffsize==45 ~ "2ndQuart",
                             (nbuildings>as.vector(qntl[1,3]) & nbuildings<=as.vector(qntl[1,4])) & buffsize==15 ~ "3rdQuart",
                             (nbuildings>as.vector(qntl[2,3]) & nbuildings<=as.vector(qntl[2,4])) & buffsize==30 ~ "3rdQuart",
                             (nbuildings>as.vector(qntl[3,3]) & nbuildings<=as.vector(qntl[3,4])) & buffsize==45 ~ "3rdQuart",
                             (nbuildings>as.vector(qntl[1,4]) & nbuildings<=as.vector(qntl[1,5])) & buffsize==15 ~ "4thQuart",
                             (nbuildings>as.vector(qntl[2,4]) & nbuildings<=as.vector(qntl[2,5])) & buffsize==30 ~ "4thQuart",
                             (nbuildings>as.vector(qntl[3,4]) & nbuildings<=as.vector(qntl[3,5])) & buffsize==45 ~ "4thQuart"))

## Joining data
# Join location, age & sex details to random points
dets <- left_join(pts, dets, by="RegionalID")
dat <- left_join(dets, trace_y, by="RegionalID")
dat <- separate(dat, 2, into=c("id", "pt_index"), sep="_", remove=FALSE) 
dat <- dat[,-c(3)]
dat$pt_index <- as.numeric(dat$pt_index)
dat <- dat %>% distinct()

# Add columns for buffer and radius
dat <- bind_rows(dat, dat, dat)
dat$buffsize <- rep(c(15,30,45), each=3380) # buffer sizes
dat <- dat %>% select(RegionalID:n.compounds.T,buffsize)

# join covariate data
dat <- dat %>% left_join(forest, by=c("pt_name", "buffsize")) 
dat <- dat %>% left_join(build, by=c("pt_name", "buffsize")) 
dat <- dat %>% left_join(baa, by=c("pt_name", "buffsize")) 
dat <- dat %>% left_join(lsm, by=c("pt_name", "buffsize"))
dat <- dat %>% left_join(stand_mn, by=c("pt_name", "buffsize"))
dat <- dat %>% left_join(stand_sd, by=c("pt_name", "buffsize"))

## Add beech mast index
dat <- left_join(dat, bmi, by=c("pt_name", "buffsize", "year")) %>%
          rename(BMI=bmi)

# add 1 to year to get lagged
bmi$year <- bmi$year + 1
dat <- left_join(dat, bmi, by=c("pt_name", "buffsize", "year")) %>%
  rename(laggedBMI=bmi)

## Fill in 0 for missing values (WUI)
dat <- dat %>% mutate(evergreen = coalesce(evergreen, 0),
                      mixed = coalesce(mixed, 0),
                      totalforest = coalesce(totalforest, 0),
                      deciduous = coalesce(deciduous, 0))

# Add WMUA 
dat <- left_join(dat, wmua, by="WMU")

# Reorder columns
dat <- dat %>% select(RegionalID:AgeClass,key,Region,WMUA_code,WMU,Town:laggedBMI) %>% 
  rename(ncomp=n.compounds.T) %>%
  select(-BMI, -laggedBMI) %>%
  pivot_wider(id_cols=c(RegionalID:ncomp), names_from=buffsize, values_from=deciduous:stand_age_sd)

# Add beechnut counts
dat <- dat %>% mutate(beechnuts=case_when(
                          year==2018 ~ 6,
                          year==2019 ~ 295,
                          year==2020 ~ 14),
                      lag_beechnuts=case_when(
                          year==2018 ~ 145,
                          year==2019 ~ 6,
                          year==2020 ~ 295))

### Save data to file ####
write_csv(dat, "data/analysis-ready/combined_AR_covars.csv")


cor(dat[,c(16:30,34:45)])

