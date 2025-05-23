## Set up rodenticide covariates
## 2022-06-27, updated 2023-01-12, 2024-09-09

library(tidyverse)

#### Read in data ####
dat <- read_csv("output/summarized_AR_results.csv")
trace_y <- read_csv("output/ncompounds_trace.csv") %>% 
  select(RegionalID, n.compounds) %>%
  rename(n.compounds.T=n.compounds)
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
evt <- read_csv("data/analysis-ready/evt_pct_filtered.csv")
cc_mean <- read_csv("data/analysis-ready/canopy-cover_mean.csv") %>%
  rename(pt_name=name) 
wui <- read_csv("data/analysis-ready/wui500_frac15.csv") %>%
          mutate(value=case_when(value==0 ~ "notWUI",
                                 value==1 ~ "intermix",
                                 value==2 ~ "interface")) %>%
          pivot_wider(names_from="value", values_from = "freq") %>%
          mutate(intermix=coalesce(intermix, 0),
                 interface=coalesce(interface, 0),
                 totalWUI=intermix+interface) %>%
          select(-`NA`) %>%
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

## Reorganize EVT (landfire)
evt$value <- factor(evt$value, levels=c(7302, 7511, 7977, 9312, 7512, 7366, 9315, 7373), 
                    labels=c("LaurentianAcadianNorthernHardwoods", "AppalachianNorthernHardwoods", "EasterCoolTemperatePasture",
                             "NortheasternNATemperatePlantation", "AppalachainHemlockHardwoods", "LaurentianAcadianPinesHemlocks",
                             "NorthCentralRuderalForest", "AcadianLowElevationSpruceFir"))
evt <- evt %>% select(name,buffsize,value,freq) %>%
  pivot_wider(names_from=value, values_from=freq) %>%
  rename(pt_name=name) %>%
  mutate(LaurentianAcadianNorthernHardwoods=coalesce(LaurentianAcadianNorthernHardwoods, 0),
         AppalachianNorthernHardwoods=coalesce(AppalachianNorthernHardwoods, 0),
         EasterCoolTemperatePasture=coalesce(EasterCoolTemperatePasture, 0),
         NortheasternNATemperatePlantation=coalesce(NortheasternNATemperatePlantation, 0),
         AppalachainHemlockHardwoods=coalesce(AppalachainHemlockHardwoods, 0),
         LaurentianAcadianPinesHemlocks=coalesce(LaurentianAcadianPinesHemlocks, 0),
         NorthCentralRuderalForest=coalesce(NorthCentralRuderalForest, 0),
         AcadianLowElevationSpruceFir=coalesce(AcadianLowElevationSpruceFir, 0)) %>%
  mutate(AppalachianHardwoodsHemlocks=AppalachianNorthernHardwoods+AppalachainHemlockHardwoods) %>%
  select(-AppalachianNorthernHardwoods,-AppalachainHemlockHardwoods, -EasterCoolTemperatePasture)

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
dat <- dat %>% left_join(cc_mean, by=c("pt_name", "buffsize"))
dat <- dat %>% left_join(evt, by=c("pt_name", "buffsize"))
dat <- dat %>% left_join(wui, by=c("pt_name", "buffsize"))

## Add beech mast index
dat <- left_join(dat, bmi, by=c("pt_name", "buffsize", "year")) %>%
          rename(BMI=bmi)

# add 1 to year to get lagged
bmi$year <- bmi$year + 1
dat <- left_join(dat, bmi, by=c("pt_name", "buffsize", "year")) %>%
  rename(laggedBMI=bmi)

## Fill in 0 for missing values 
dat <- dat %>% mutate(evergreen = coalesce(evergreen, 0),
                      mixed = coalesce(mixed, 0),
                      totalforest = coalesce(totalforest, 0),
                      deciduous = coalesce(deciduous, 0))

# Add WMUA 
dat <- left_join(dat, wmua, by="WMU")

# Reorder columns
dat <- dat %>% select(RegionalID:AgeClass,key,Region,WMUA_code,WMU,Town:laggedBMI) %>% 
  rename(ncomp=n.compounds.T) %>%
  # select(-BMI, -laggedBMI) %>%
  filter(buffsize==15)
  # pivot_wider(id_cols=c(RegionalID:ncomp), names_from=buffsize, values_from=deciduous:AppalachianHardwoodsHemlocks)

# Add beechnut counts
dat <- dat %>% mutate(beechnuts=case_when(
                          year==2018 ~ -41.3,
                          year==2019 ~ 247.7,
                          year==2020 ~ -33.3),
                      lag_beechnuts=case_when(
                          year==2018 ~ 97.7,
                          year==2019 ~ -41.3,
                          year==2020 ~ 247.7))

### Save data to file ####
write_csv(dat, "data/analysis-ready/combined_AR_covars.csv")


# cor_mat <- cor(dat[,c(16:30,34:63)])
# cor_mat <- as.data.frame(cor_mat) %>% rownames_to_column(var="variable")
# write_csv(cor_mat, "output/correlation_matrix_20240909.csv")

# check <- dat %>% filter(pt_index==1)
# mean(check$ncomp)
# var(check$ncomp)