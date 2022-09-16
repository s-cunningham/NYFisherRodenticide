## Plotting prediction curves
## 2022-09-14

library(tidyverse)
library(patchwork)

#### Number of compounds ####
dat <- read_csv("output/model_data.csv")




#### Age-specific Number of Compounds ####




#### Binary analyses #####
dat2 <- read_csv("output/binary_model_data.csv")

# Subset by compound
brod <- dat2[dat2$compound=="Brodifacoum",]
brom <- dat2[dat2$compound=="Bromadiolone",]
diph <- dat2[dat2$compound=="Diphacinone",]

