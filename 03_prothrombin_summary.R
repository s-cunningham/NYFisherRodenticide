## 

library(tidyverse)
library(scales)
library(colormap)
library(viridis)
library(gridExtra)

## Read PT 
coag <- read.csv("data/analysis-ready/blood_clotting_20220425.csv")
names(coag)[c(1,3:5)] <- c("FisherID", "PT", "PTdilute", "fibrinogen")

coag <- coag[is.na(coag$comments),]
coag$PT <- as.numeric(coag$PT)
coag$PTdilute <- as.numeric(coag$PTdilute)
coag$fibrinogen <- as.numeric(coag$fibrinogen)

## Deaths
dead <- read.csv("data/analysis-ready/deaths_current20210803.csv")

## Join PT to deaths
dead <- left_join(dead, coag, by="FisherID")

arcause <- dead[dead$Cause.of.Death=="Rodenticide" | dead$Cause.of.Death=="Possible rodenticide",]

deadpt <- dead[!is.na(dead$PT), ]
deadpt <- deadpt[is.na(deadpt$comments),]

deadpt$PT <- as.numeric(deadpt$PT)
deadpt$PTdilute <- as.numeric(deadpt$PTdilute)
deadpt$fibrinogen <- as.numeric(deadpt$fibrinogen)

### Make heatmap
coag$StudyArea <- "Tug Hill"
coag$StudyArea[str_sub(coag$FisherID, 4, 4)=="5"] <- "Central Adirondack"

# convert to long
coag <- pivot_longer(coag, 3:5, names_to="testname", values_to="testvalue")



coag1 <- coag[coag$testname=="PT", ]
# coag1$testvalue[coag1$testname=="PTdilute" | coag1$testname=="fibrinogen"] <- NA
coag2 <- coag[coag$testname=="PTdilute", ]
# coag2$testvalue[coag2$testname=="PT" | coag2$testname=="fibrinogen"] <- NA
coag3 <- coag[coag$testname=="fibrinogen", ]
# coag3$testvalue[coag3$testname=="PT" | coag3$testname=="PTdilute"] <- NA

ggplot(coag2) + 
  geom_boxplot(aes(y=testvalue)) + 
  theme_bw() + ylab("Prothrombin Time (dilute)") +
  theme(
    axis.text.y=element_text(size=12), axis.title=element_text(size=14),
    axis.text.x=element_blank(),
    strip.text.x=element_text(size=12, face="bold"))



lowcol <- rgb(68,1,84, max=255)
midcol <- rgb(33,144,141, max=255)
highcol <- rgb(253,231,37, max=255)

show_col(colormap(colormaps$viridis), labels=TRUE)
colormap(colormaps$viridis, nshades=71, format="rgb")

ggplot() +
  geom_tile(data=coag1, aes(x=as.factor(FisherID), y=as.factor(testname), fill=testvalue), colour = "black") +
  scale_fill_gradient2(midpoint=10.51, low=lowcol, mid=midcol, high=highcol, limits=c(8.3, 52.9)) +
  theme_classic() + xlab("Fisher ID") + ylab("Prothrombin Time") +
  theme(#legend.position = 'none',
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text.x=element_text(size=12, face="bold"),
        plot.title=element_text(size=16, face="bold", hjust=0.5)) +
  facet_grid(.~StudyArea, space="free_x", scales="free") 
  
  
 +
  geom_tile(data=coag2, aes(x=as.factor(FisherID), y=as.factor(testname), fill=testvalue), colour = "black") +
  scale_fill_gradient2(midpoint=14, low="#2166ac", mid="white", high="#b2182b",  na.value=mycol) +
  geom_tile(data=coag3, aes(x=as.factor(FisherID), y=as.factor(testname), fill=testvalue), colour = "black") +
  scale_fill_gradient2(midpoint=245, low="#2166ac", mid="white", high="#b2182b", na.value=mycol) +
  theme_classic() + xlab("Test") 
  
  
   # + 
  theme(legend.position = 'none',
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text.y=element_text(size=12, face="bold"),
        strip.text.x=element_blank(),
        plot.title=element_text(size=16, face="bold", hjust=0.5)) +
  facet_grid(.~StudyArea, space="free_x", scales="free",) 


