library(ape)
library(tidyverse)

pts <- read_csv("output/random_point_locs.csv")
dat <- read_csv("output/summarized_AR_results.csv")
trace_y <- read_csv("output/ncompounds_trace.csv") %>% 
  select(RegionalID, n.compounds) %>%
  rename(n.compounds.T=n.compounds)
trace_n <- read_csv("output/ncompounds_notrace.csv") %>% 
  select(RegionalID, n.compounds) %>%
  rename(n.compounds.MO=n.compounds)

## Number of compounds detected
ncomps <- left_join(trace_n, trace_y, by="RegionalID")

## Save details of each sample
dets <- dat %>% select(RegionalID:WMU,Town) %>% distinct()

## Add a column to points with just sample ID
pts <- pts %>% select(RegionalID,name,x,y) %>%
  rename(rand_x=x, rand_y=y, pt_name=name)

dets <- left_join(pts, dets, by="RegionalID")
dat <- left_join(dets, ncomps, by="RegionalID")
dat <- separate(dat, 2, into=c("id", "pt_index"), sep="_", remove=FALSE) 
dat <- dat[,-c(3)]
dat$pt_index <- as.numeric(dat$pt_index)
dat <- dat %>% distinct()

brod <- read_csv("output/binary_brodifacoum.csv") %>%
  left_join(pts, by="RegionalID") %>%
  separate(pt_name, into=c("name", "pt_index"), sep="_", remove=FALSE) %>% select(-name)
brom <- read_csv("output/binary_bromadiolone.csv") %>%
  left_join(pts, by="RegionalID")%>%
  separate(pt_name, into=c("name", "pt_index"), sep="_", remove=FALSE) %>% select(-name)
diph <- read_csv("output/binary_diphacinone.csv") %>%
  left_join(pts, by="RegionalID")%>%
  separate(pt_name, into=c("name", "pt_index"), sep="_", remove=FALSE) %>% select(-name)
dico <- read_csv("output/binary_dicoumarol.csv") %>%
  left_join(pts, by="RegionalID")%>%
  separate(pt_name, into=c("name", "pt_index"), sep="_", remove=FALSE) %>% select(-name)

#### Semivariogram ####
dat1 <- dat %>% select(RegionalID,rand_x,rand_y, Region, n.compounds.T)
dat1 <- unique(dat1)
sp::coordinates(dat1) <- ~rand_x+rand_y
vario <- gstat::variogram(n.compounds.T~1, data=dat1)
plot(vario)

## Moran's I

# Number of compounds
ncomps <- data.frame()
for(i in 1:10) {
  midat <- dat %>% filter(pt_index==i)
  # midat$compo.T <- ifelse(midat$n.compounds.T==0, 0, 1)
  midat <- distinct(midat)

  # Create distance matrix
  ar.dists <- as.matrix(dist(cbind(midat$rand_x, midat$rand_y)))

  # create inverse distance matrix
  ar.dists.inv <- 1/ar.dists
  diag(ar.dists.inv) <- 0
  # ar.dists.inv[1:5, 1:5]

  # Moran's I
  ncomps <- bind_rows(ncomps, as.data.frame(ape::Moran.I(midat$n.compounds.T, ar.dists.inv)))
}

colMeans(ncomps)

# Brodifacoum
brodmi <- data.frame()
for (i in 1:10) {
  
  midat <- brod %>% filter(pt_index==i)
  
  # Create distance matrix
  ar.dists <- as.matrix(dist(cbind(midat$rand_x, midat$rand_y)))
  
  # create inverse distance matrix
  ar.dists.inv <- 1/ar.dists
  diag(ar.dists.inv) <- 0
  
  # Moran's I
  brodmi <- bind_rows(brodmi, as.data.frame(ape::Moran.I(midat$bin.exp, ar.dists.inv)))
  
}
colMeans(brodmi)

# Bromadiolone
brommi <- data.frame()
for (i in 1:10) {
  
  midat <- brom %>% filter(pt_index==i)
  
  # Create distance matrix
  ar.dists <- as.matrix(dist(cbind(midat$rand_x, midat$rand_y)))
  
  # create inverse distance matrix
  ar.dists.inv <- 1/ar.dists
  diag(ar.dists.inv) <- 0
  
  # Moran's I
  brommi <- bind_rows(brommi, as.data.frame(ape::Moran.I(midat$bin.exp, ar.dists.inv)))
  
}
colMeans(brommi)

# Diphacinone
diphmi <- data.frame()
for (i in 1:10) {
  
  midat <- diph %>% filter(pt_index==i)
  
  # Create distance matrix
  ar.dists <- as.matrix(dist(cbind(midat$rand_x, midat$rand_y)))
  
  # create inverse distance matrix
  ar.dists.inv <- 1/ar.dists
  diag(ar.dists.inv) <- 0
  
  # Moran's I
  diphmi <- bind_rows(diphmi, as.data.frame(ape::Moran.I(midat$bin.exp, ar.dists.inv)))
  
}
colMeans(diphmi)

# Dicoumarol
dicomi <- data.frame()
for (i in 1:10) {
  
  midat <- dico %>% filter(pt_index==i)
  
  # Create distance matrix
  ar.dists <- as.matrix(dist(cbind(midat$rand_x, midat$rand_y)))
  
  # create inverse distance matrix
  ar.dists.inv <- 1/ar.dists
  diag(ar.dists.inv) <- 0
  
  # Moran's I
  dicomi <- bind_rows(dicomi, as.data.frame(ape::Moran.I(midat$bin.exp, ar.dists.inv)))
  
}
colMeans(dicomi)
