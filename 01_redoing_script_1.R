library(tidyverse)
library(sf)
library(terra)
library(stars)

#### Step 1: QA/QC on fisher harvest data ####

# Read in age data
age2018 <- read.csv("data/20211111_2018_age_data.csv")
age2019 <- read.csv("data/20211111_2019_age_data.csv")
age2020 <- read.csv("data/20211111_2020_age_data.csv")

# Bind together and rename. Fill missing values with NA. Change date type.
ages <- rbind(age2018,age2019,age2020)
names(ages)[c(1,7,9,11,12)] <- c("HarvestDate", "AgeClass", "AgeRange", "TrapperID", "RegionalSampleID")
ages$AgeRange[ages$AgeRange==""] <- NA
ages$TrapperID[ages$TrapperID==""] <- NA
ages$HarvestDate <- as.Date(ages$HarvestDate, format="%m/%d/%Y")
ages <- ages %>% as_tibble()

# Separate WMU by region
ages <- ages %>% mutate(WMU=str_to_upper(WMU))
ages <- separate(ages, 4, into=c("Region", "x"), sep="[A-Z]", remove=FALSE)
ages <- ages[,-6]

# Add column for submission year
ages$SubmissionYear <- c(rep(2016,1202), rep(2017, 1344), rep(2018,1761), rep(2019, 1448), rep(2020, 1795))

# Add column for harvest year
ages$HarvestYear <- format(ages$HarvestDate, "%Y")

# Remove samples that don't have town/county data
no_loc <- ages[ages$Town=="",]
ages <- ages[ages$Town!="",]

# Remove trailing whitespace
ages$Town <- str_trim(ages$Town, side="right")
ages$County <- str_trim(ages$County, side="right")

# Correct spelling errors
ages$Town[ages$Town=="Dannamora"] <- "Dannemora"
ages$Town[ages$Town=="Loranne" | ages$Town=="Loraine"] <- "Lorraine"
ages$Town[ages$Town=="Bleeker"] <- "Bleecker"
ages$Town[ages$Town=="Housfield"] <- "Hounsfield"
ages$Town[ages$Town=="Montegue"] <- "Montague"
ages$Town[ages$Town=="Lincklean" | ages$Town=="Linklaen"] <- "Lincklaen"
ages$Town[ages$Town=="Pickney" | ages$Town=="Pinckey"] <- "Pinckney"
ages$Town[ages$Town=="Pierpont" | ages$Town=="Pirrpont"] <- "Pierrepont"
ages$Town[ages$Town=="Strattford"] <- "Stratford"
ages$Town[ages$Town=="Ansville"] <- "Annsville"
ages$Town[ages$Town=="Agusta"] <- "Augusta"
ages$Town[ages$Town=="German Flats" | ages$Town=="Germanflat" ] <- "German Flatts"
ages$Town[ages$Town=="Camren"] <- "Cameron"
ages$Town[ages$Town=="tyrone"] <- "Tyrone"
ages$Town[ages$Town=="Middlevileel"] <- "Middleville"
ages$Town[ages$Town=="Schylar" | ages$Town=="Schuylar"] <- "Schuyler"
ages$Town[ages$Town=="Red field" | ages$Town=="Red Field" | ages$Town=="Redfiled"] <- "Redfield"
ages$Town[ages$Town=="Hardenburg"] <- "Hardenburgh"
ages$Town[ages$Town=="Petersburg"] <- "Petersburgh"
ages$Town[ages$Town=="Leray"] <- "Le Ray"
ages$Town[ages$Town=="Poladn"] <- "Poland"
ages$Town[ages$Town=="Petersburg"] <- "Petersburgh"
ages$Town[ages$Town=="West Dale"] <- "Westdale"
ages$Town[ages$Town=="Duanseburg"] <- "Duanesburg"
ages$Town[ages$Town=="Flornce" | ages$Town=="Flourence"] <- "Florence"
ages$Town[ages$Town=="Leigh"] <- "Lee"
ages$Town[ages$Town=="Contantia" | ages$Town=="Consantia"] <- "Constantia"
ages$Town[ages$Town=="Boyelston" | ages$Town=="Boylsotn"] <- "Boylston"
ages$Town[ages$Town=="Edmaston"] <- "Edmeston"
ages$Town[ages$Town=="Ellenbury"] <- "Ellenburg"
ages$Town[ages$Town=="St. Armand"] <- "St Armand"
ages$Town[ages$Town=="Lychfield" | ages$Town=="Lynchfield"] <- "Litchfield"
ages$Town[ages$Town=="Schaticoke"] <- "Schaghticoke"
ages$Town[ages$Town=="Stueben" | ages$Town=="Steben" | ages$Town=="Stuben"] <- "Steuben"
ages$Town[ages$Town=="Tomkins"] <- "Tompkins"
ages$Town[ages$Town=="Lyden"] <- "Leyden"
ages$Town[ages$Town=="Cornith"] <- "Corinth"
ages$Town[ages$Town=="Deruyter"] <- "DeRuyter"
ages$Town[ages$Town=="Sangerfiled"] <- "Sangerfield"
ages$Town[ages$Town=="Charlston"] <- "Charleston"
ages$Town[ages$Town=="Pitcarin"] <- "Pitcairn"
ages$Town[ages$Town=="Otsellic" | ages$Town=="Ostelic"] <- "Otselic"
ages$Town[ages$Town=="Forrestport"] <- "Forestport"
ages$Town[ages$Town=="Au Sable"] <- "AuSable"
ages$Town[ages$Town=="Windfield"] <- "Winfield"
ages$Town[ages$Town=="Remensen"] <- "Remsen"
ages$Town[ages$Town=="Alvion"] <- "Albion"
ages$Town[ages$Town=="Russel"] <- "Russell"
ages$Town[ages$Town=="Lweis"] <- "Lewis"
ages$Town[ages$Town=="Salsbury"] <- "Salisbury"
ages$Town[ages$Town=="OLIVE"] <- "Olive"
ages$Town[ages$Town=="Wiliamstown"] <- "Williamstown"
ages$Town[ages$Town=="Fairfeild"] <- "Fairfield"
ages$Town[ages$Town=="St. Johnsville"] <- "St Johnsville"
ages$Town[ages$Town=="Crohgan" | ages$Town=="Crogan"] <- "Croghan"
ages$Town[ages$Town=="Harriburg" | ages$Town=="harrisburg"] <- "Harrisburg"
ages$Town[ages$Town=="Wadington"] <- "Waddington"
ages$Town[ages$Town=="Forestburg"] <- "Forestburgh"
ages$Town[ages$Town=="Matinsburg"] <- "Martinsburg"
ages$Town[ages$Town=="Grieg" | ages$Town=="Graig"] <- "Greig"
ages$Town[ages$Town=="Marshal"] <- "Marshall"
ages$Town[ages$Town=="Dekalb"] <- "De Kalb"
ages$Town[ages$Town=="SUMMIT"] <- "Summit"
ages$Town[ages$Town=="Palentine"] <- "Palatine"
ages$Town[ages$Town=="Oreleans"] <- "Orleans"
ages$Town[ages$Town=="Dianna"] <- "Diana"
ages$Town[ages$Town=="amboy"] <- "Amboy"
ages$Town[ages$Town=="ATHENS"] <- "Athens"
ages$Town[ages$Town=="Guildford"] <- "Guilford"
ages$Town[ages$Town=="Middleburg"] <- "Middleburgh"
ages$Town[ages$Town=="fine"] <- "Fine"
ages$Town[ages$Town=="Poestenskill"] <- "Poestenkill"
ages$Town[ages$Town=="Coeyman"] <- "Coeymans"
ages$Town[ages$Town=="Wilma"] <- "Wilna"
ages$Town[ages$Town=="Vanetten"] <- "Van Etten"
ages$Town[ages$Town=="Tabers" | ages$Town=="Taburg"] <- "Taberg"
ages$Town[ages$Town=="Moria"] <- "Moriah"
ages$Town[ages$Town=="lumberland"] <- "Lumberland"
ages$Town[ages$Town=="Marblet"] <- "Marbletown"
ages$Town[ages$Town=="Fallsburgh"] <- "Fallsburg"
ages$Town[ages$Town=="Blenhiem"] <- "Blenheim"
ages$Town[ages$Town=="Lake Luzurne"] <- "Lake Luzerne"
ages$Town[ages$Town=="Chenengo"] <- "Chenango"
ages$Town[ages$Town=="Manaheim"] <- "Manheim"
ages$Town[ages$Town=="Turn"] <- "Turin"
ages$Town[ages$Town=="Plattek"] <- "Plattekill"
ages$Town[ages$Town=="Exter"] <- "Exeter"
ages$Town[ages$Town=="Rennslearville"] <- "Rensselaerville"
ages$Town[ages$Town=="St. Vincent"] <- "Cape Vincent"
ages$Town[ages$Town=="Brandalow"] <- "Broadalbin"
ages$County[ages$County=="St. Lawrsnce" | ages$County=="St.Lawrence" | 
              ages$County=="St. Lawrence"] <- "St Lawrence"

# Add column to preserve village/hamlet name
ages$Village <- ages$Town

# Change village/hamlet name to town name
ages$Town[ages$Town=="Mariaville" | ages$Town=="Delanson"] <- "Duanesburg"
ages$Town[ages$Town=="Caroga Lake"] <- "Caroga"
ages$Town[ages$Town=="Schenevus"] <- "Maryland"
ages$Town[ages$Town=="Pierrpont"] <- "Pierrepont"
ages$Town[ages$Town=="Taborton" | ages$Town=="Averill Park" | ages$Town=="Averil Park"] <- "Sand Lake"
ages$Town[ages$Town=="Melrose"] <- "Schaghticoke"
ages$Town[ages$Town=="Taberg"] <- "Annsville"
ages$Town[ages$Town=="Blossvale" | ages$Town=="Blassvale" | 
            ages$Town=="Blossvalle" | ages$Town=="Taeberg" | ages$Town=="Bossvalle" | ages$Town=="Bossvale"] <- "Annsville"
ages$Town[ages$Town=="Brasher Falls"] <- "Brasher"
ages$Town[ages$Town=="West Valley"] <- "Ashford"
ages$Town[ages$Town=="Little Genesee"] <- "Genesee"
ages$Town[ages$Town=="Pattersonville"] <- "Rotterdam"
ages$Town[ages$Town=="Panama" | ages$Town=="Niobe"] <- "Harmony"
ages$Town[ages$Town=="Tannersville"] <- "Hunter"
ages$Town[ages$Town=="Freehold"] <- "Greenville"
ages$Town[ages$Town=="Livingston Manor"] <- "Rockland"
ages$Town[ages$Town=="Richfield Springs"] <- "Richfield"
ages$Town[ages$Town=="Cold Spring"] <- "Philipstown"
ages$Town[ages$Town=="Downsville"] <- "Colchester"
ages$Town[ages$Town=="Cold Brook"] <- "Russia"
ages$Town[ages$Town=="South Corning"] <- "Corning"
ages$Town[ages$Town=="West Clarksville"] <- "Clarksville"
ages$Town[ages$Town=="Cleveland"] <- "Constantia"
ages$Town[ages$Town=="Pointrock"] <- "Lee"
ages$Town[ages$Town=="Modena"] <- "Plattekill"
ages$Town[ages$Town=="Walden"] <- "Montgomery"
ages$Town[ages$Town=="Dalton" | ages$Town=="Dalten"] <- "Nunda"
ages$Town[ages$Town=="Star Lake"] <- "Clifton"
ages$Town[ages$Town=="Livingstonville"] <- "Broome"
ages$Town[ages$Town=="Pulaski"] <- "Richland"
ages$Town[ages$Town=="Belmont"] <- "Amity"
ages$Town[ages$Town=="Whitney Point"] <- "Triangle"
ages$Town[ages$Town=="South Otselic"] <- "Otselic"
ages$Town[ages$Town=="Natural Bridge"] <- "Wilna"
ages$Town[ages$Town=="Fort Drum" | ages$Town=="Evans Mills"] <- "Le Ray"
ages$Town[ages$Town=="Risco" | ages$Town=="Roscoe"] <- "Rockland"
ages$Town[ages$Town=="Garratvilles"] <- "New Lisbon"
ages$Town[ages$Town=="Brushton"] <- "Moira"
ages$Town[ages$Town=="West Winfield"] <- "Winfield"
ages$Town[ages$Town=="Painted Post"] <- "Erwin"
ages$Town[ages$Town=="Cattaraugus"] <- "New Albion"
ages$Town[ages$Town=="West Winfield"] <- "Winfield"
ages$Town[ages$Town=="Vernon Center"] <- "Vernon"
ages$Town[ages$Town=="Parksville"] <- "Liberty"
ages$Town[ages$Town=="Middleville"] <- "Fairfield"
ages$Town[ages$Town=="New Kingston" | ages$Town=="Margaretville"] <- "Middletown"
ages$Town[ages$Town=="Acra"] <- "Cairo"
ages$Town[ages$Town=="Slingerlands"] <- "Bethlehem"
ages$Town[ages$Town=="Harris"] <- "Thompson"
ages$Town[ages$Town=="Starkville"] <- "Stark"
ages$Town[ages$Town=="Smithville Flats"] <- "Smithville"
ages$Town[ages$Town=="Deer Park"] <- "Babylon"
ages$Town[ages$Town=="Harrisville"] <- "Diana"
ages$Town[ages$Town=="Wanakena"] <- "Fine"
ages$Town[ages$Town=="Fort Plain"] <- "Minden"
ages$Town[ages$Town=="Rexville"] <- "West Union"
ages$Town[ages$Town=="Cassadega" | ages$Town=="Cassada"] <- "Stockton"
ages$Town[ages$Town=="Leeds"] <- "Catskill"
ages$Town[ages$Town=="Loon Lake"] <- "Franklin"
ages$Town[ages$Town=="Sheds"] <- "DeRuyter"
ages$Town[ages$Town=="Cherry Plain"] <- "Berlin"
ages$Town[ages$Town=="West Herkimer"] <- "Herkimer"
ages$Town[ages$Town=="Lagrangeville"] <- "La Grange"
ages$Town[ages$Town=="West Fulton"] <- "Fulton"
ages$Town[ages$Town=="Erieville"] <- "Nelson"
ages$Town[ages$Town=="Johnsonville"] <- "Pittstown"
ages$Town[ages$Town=="Medusa"] <- "Rensselaerville"
ages$Town[ages$Town=="North Lawrence"] <- "Lawrence"
ages$Town[ages$Town=="Apalachin"] <- "Owego"
ages$Town[ages$Town=="Lacona"] <- "Sandy Creek"
ages$Town[ages$Town=="Grand Gorge"] <- "Roxbury"
ages$Town[ages$Town=="Black Creek"] <- "New Hudson"
ages$Town[ages$Town=="Mayville"] <- "Chautauqua"
ages$Town[ages$Town=="Maplecrest"] <- "Windham"
ages$Town[ages$Town=="Cooperstown"] <- "Otsego"
ages$Town[ages$Town=="Smithville Flatts"] <- "Smithville"
ages$Town[ages$Town=="Carthage"] <- "Wilna"
ages$Town[ages$Town=="East Pharsalia" | ages$Town=="North Pharsalia"] <- "Pharsalia"
ages$Town[ages$Town=="East Herkimer"] <- "Herkimer"
ages$Town[ages$Town=="Waterville"] <- "Sangerfield"
ages$Town[ages$Town=="Selkirk"] <- "Bethlehem"
ages$Town[ages$Town=="Alcove"] <- "Coeymans"
ages$Town[ages$Town=="Rockhill"] <- "Thompson"
ages$Town[ages$Town=="Hunt"] <- "Portage"
ages$Town[ages$Town=="Deerpar"] <- "Babylon"
ages$Town[ages$Town=="North Creek"] <- "Johnsburg"
ages$Town[ages$Town=="Mowhawk"] <- "German Flatts"
ages$Town[ages$Town=="Chase Mills"] <- "Waddington"
ages$Town[ages$Town=="Westdale"] <- "Camden"
ages$Town[ages$Town=="Selkirk"] <- "Bethlehem"
ages$Town[ages$Town=="South New Berlin"] <- "New Berlin"
ages$Town[ages$Town=="Alexandria Bay"] <- "Alexandria"
ages$Town[ages$Town=="New Woodstock"] <- "Cazenovia"
ages$Town[ages$Town=="Old Forge"] <- "Webb"
ages$Town[ages$Town=="Lockwood"] <- "Barton"
ages$Town[ages$Town=="Lusselville"] <- "Ephratah"
ages$Town[ages$Town=="Wylden"] <- "Ava"
# ages$Town[ages$Town=="West Lydon" & ages$County=="Oneida"] <- "Western"
ages$Town[ages$Town=="Pittfield"] <- "Pittsfield"

# Other corrections
ages$TrapperID[ages$Town=="Monettca"] <- "2666-8000-0482"
ages$Town[ages$Town=="Monettca"] <- "Arietta"
ages$Town[which(ages$RegionalSampleID=="2018-6389" | ages$RegionalSampleID=="2018-6386")] <- "Fairfield"
ages$Town[ages$Town=="Clinton" & ages$County=="St. Lawrence" & ages$WMU=="6F"] <- "Clifton"
ages$County[ages$Town=="Freedom" & ages$County=="Allegany"] <- "Cattaraugus"

#### Step 2: Add spatial context to harvest data ####

# Read in town shapefile
twn <- st_read("data/spatial/Cities_Towns.shp")
twn <- st_transform(twn, crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
cents <- twn %>% st_centroid() %>% st_coordinates()
centid <- twn$NAME
muntype <- twn$MUNI_TYPE
centlab <- data.frame(id=centid, longitude=cents[,1], latitude=cents[,2], munic=muntype)
names(centlab)[1] <- "Town"

cities <- subset(centlab, munic=="city")
centlab <- subset(centlab, munic=="town")

dup <- subset(centlab, duplicated(Town))[,1]
et <- centlab[centlab$Town %in% dup,]

ggplot(et, aes(x=longitude, y=latitude, color=Town)) + geom_point()


dapart <- data.frame()
for (i in 1:length(unique(et$Town))) {
  
  temp <- et[et$Town==unique(et$Town)[i],]
  dist <- st_distance(st_as_sf(temp[1,2:3], coords=c("longitude", "latitude"), crs=26918), 
                      st_as_sf(temp[2,2:3], coords=c("longitude", "latitude"), crs=26918), lonlat=FALSE)/1000
  
  new <- data.frame(Town=temp$Town[1], Dist_km=dist)
  dapart <- rbind(dapart, new)
}

age2 <- ages[ages$Town %in% et$Town,]

N.order <- order(age2$Town, decreasing=FALSE)
age2 <- age2[N.order,]


# Join town coordinates
ages <- left_join(ages, centlab, by="Town", relationship="many-to-many")

# Check for missing towns
miss <- ages[is.na(ages$longitude),]
unique(miss$Town)

ages$longitude[ages$Town=="Rochester"] <- 1776022
ages$latitude[ages$Town=="Rochester"] <- 2293860
ages$longitude[ages$Town=="Albion"] <- 1595844
ages$latitude[ages$Town=="Albion"] <- 2446223
ages$longitude[ages$Town=="Ashland"] <- 1758895
ages$latitude[ages$Town=="Ashland"] <- 2347964
ages$longitude[ages$Town=="Chester"] <- 1761007
ages$latitude[ages$Town=="Chester"] <- 2505715
ages$longitude[ages$Town=="Dickinson"] <- 1681869
ages$latitude[ages$Town=="Dickinson"] <- 2606303
ages$longitude[ages$Town=="Franklin"] <- 1696534
ages$latitude[ages$Town=="Franklin"] <- 2335408
ages$longitude[ages$Town=="Fremont"] <- 1714738
ages$latitude[ages$Town=="Fremont"] <- 2285173
ages$longitude[ages$Town=="Greenville"& ages$County=="Greene"] <- 1782873
ages$latitude[ages$Town=="Greenville" & ages$County=="Greene"] <- 2362152
ages$longitude[ages$Town=="Greenville"& ages$County=="Orange"] <- 1761780
ages$latitude[ages$Town=="Greenville" & ages$County=="Orange"] <- 2240094
ages$longitude[ages$Town=="Lewis"& ages$County=="Lewis"] <- 1633203
ages$latitude[ages$Town=="Lewis" & ages$County=="Lewis"] <- 2452525
ages$longitude[ages$Town=="Lewis"& ages$County=="Essex"] <- 1767925
ages$latitude[ages$Town=="Lewis" & ages$County=="Essex"] <- 2577509
ages$longitude[ages$Town=="Middletown"& ages$County=="Delaware"] <- 1739612
ages$latitude[ages$Town=="Middletown" & ages$County=="Delaware"] <- 2325361
ages$longitude[ages$Town=="Hornell"] <- 1492045
ages$latitude[ages$Town=="Hornell"] <- 2292524
ages$longitude[ages$Town=="Rensselaer"] <- 1798357
ages$latitude[ages$Town=="Rensselaer"] <- 2394953
ages$longitude[ages$Town=="Rome"] <- 1644803
ages$latitude[ages$Town=="Rome"] <- 2426122

# Add column for string length of Regional Sample ID
ages$strl <- str_length(ages$RegionalSampleID)
ages$RegionalSampleID <- ifelse(ages$strl<=5, paste(ages$SubmissionYear, ages$RegionalSampleID, sep="-"), ages$RegionalSampleID)

ages <- ages[,1:18]
names(ages)[13] <- "RegionalID" 

## Save cleaned-up file
ages2 <- ages %>% dplyr::select(RegionalID, TrapperID, HarvestDate, HarvestYear, Sex, Age, AgeClass, Region,
                                WMU, County, Town, Village) 

# Save cleaned file
write_csv(ages2, "data/analysis-ready/2016-2020_ages_data.csv")

#### Step 3: Forest cover grid ####

# load NYS outline
latlon <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
nys <- st_read("data/spatial/NYS_outline_albers.shp")
spacing <- 10000
nys_grid <- st_make_grid(nys, cellsize = c(spacing, spacing),
                         what = 'polygons') 
nys_grid <- nys_grid[nys] %>% st_sf('geometry' = ., data.frame('ID' = 1:length(.)))

# Read percent forest data
plc <- read.csv("data/grid_pct_forest.csv")
plc <- plc[,-1]
plc <- cbind(seq(1,1415,1), plc)
names(plc)[1] <- "ID"

# combine
grid2 <- left_join(nys_grid, plc, by="ID")

ggplot() + 
  geom_sf(data=grid2, aes(fill=Forest)) +
  scale_fill_gradient(low="gray90", high="#006837") +
  theme_classic()

## Convert sf object to stars object
grid.r <- st_rasterize(grid2 %>% dplyr::select(Forest, geometry))

# convert stars object to SpatRast
grid.r <- rast(grid.r)

# plot to check it works
plot(grid.r)

# write to file
writeRaster(grid.r, "data/rasters/forest_grid_sample_select.tif")

## Read in polygon with harvest area
aea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
fha <- st_read("data/spatial/FisherHarvestArea.shp")
fha <- st_transform(fha, crs=st_crs(grid.r))
plot(st_geometry(fha))
# fha2 <- fortify(fha, region="Name")


for.ha <- mask(x=grid.r, mask=fha)
for.ha2 <- as.polygons(for.ha, values=TRUE, na.rm=TRUE, na.all=TRUE)



rfor <- for.ha2@data$Forest 
rfor <- as.data.frame(cbind(seq(1,963,1),rfor))
names(rfor) <- c("id", "Forest")
rfor$id <- as.character(rfor$id)
for.ha2 <- fortify(for.ha2)
for.ha2 <- left_join(for.ha2, rfor, by="id")


#### 2018 Livers ####


#### 2019 Livers ####


# subset to 2019 and non-missing towns
ages <- subset(ages, !is.na(longitude) & HarvestYear==2019) 
ages <- unique(ages)

# Read livers from 2020
liv19 <- read.csv("data/2019livers_456789.csv")
liv19$SubmissionYear <- 2019

liv19$strl <- str_length(liv19$RegionalID)
liv19$RegionalID <- ifelse(liv19$strl<=5, paste(liv19$SubmissionYear, liv19$RegionalID, sep="-"), liv19$RegionalID)

# Join liver info to metadata
dat <- ages[ages$RegionalID %in% liv19$RegionalID,]

lv19 <- dat %>% group_by(Town) %>% count()
lv19 <- left_join(lv19, centlab, by="Town")

ggplot() + 
  geom_polygon(data=nys, aes(x=long, y=lat, group=group), fill="gray40") +
  geom_polygon(data=for.ha2, aes(x=long, y=lat, group=group, fill=Forest), color="gray40") +
  geom_point(data=lv19, aes(x=longitude, y=latitude, size=n), color="gray80", pch=21, fill="red") +
  scale_fill_gradient(low="#ffffe5",high="#004529") +
  theme_void() +
  theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))

plot(dat$longitude, dat$latitude)

dat %>% group_by(Region) %>% count()

ggplot(dat, aes(x=longitude, y=latitude, color=Region)) + geom_point() + theme_bw()

dat2 <- dat %>% group_by(Town) %>% count()
dat2 <- left_join(dat2, centlab, by="Town")

# make new SPDF to extract points
dat <- dat[,c(13:15,1:8,17:18)]

xy <- dat[c(12,13)]
dat.sp <- SpatialPointsDataFrame(xy, dat, proj4string=CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

for.grid <- raster::extract(grid.r, dat.sp, method="simple") 

dat$PctForest <- for.grid

hist(dat$PctForest)

dat$ForestGroup <- 1
dat$ForestGroup[dat$PctForest>48 & dat$PctForest<=69] <- 2
dat$ForestGroup[dat$PctForest>69 & dat$PctForest<=91.8] <- 3
dat$ForestGroup[dat$PctForest>91.8 & dat$PctForest<=100] <- 4

dat$AG <- ifelse(dat$Age < 3.5, "J", "A")
dat <- unite(dat, "AgeSex", c(16,9), sep="", remove=FALSE)

dat %>% group_by(Age) %>% count()

N.order <- order(dat$HarvestDate, decreasing=FALSE)
dat <- dat[N.order,]

# write.csv(dat, "data/2019_forest_groups.csv")

# Select livers
nal <- data.frame()
for (i in 1:4) {
  temp <- dat[dat$ForestGroup==i,]
  
  ag <- temp %>% group_by(AG) %>% count()
  
  if (ag$n[ag$AG=="A"] >= 15) {
    
    juv <- temp[temp$Age < 3.5,]
    jdat <- juv[sample(nrow(juv), 15, replace=FALSE),]
    
    ad <- temp[temp$Age >= 3.5,]
    addat <- ad[sample(nrow(ad), 15, replace=FALSE),]
    
  } else {
    ma <- min(ag[,2])
    
    ad <- temp[temp$Age >= 3.5,]
    addat <- ad[sample(nrow(ad), ma, replace=FALSE),]
    
    juv <- temp[temp$Age < 3.5,]
    jdat <- juv[sample(nrow(juv), (30-ma), replace=FALSE),]
  }
  
  nal <- rbind(nal, jdat, addat)
}
lv <- nal %>% group_by(Town) %>% count()
lv <- left_join(lv, centlab, by="Town")

ggplot() + 
  geom_polygon(data=nys, aes(x=long, y=lat, group=group), fill="gray40") +
  geom_polygon(data=for.ha2, aes(x=long, y=lat, group=group, fill=Forest), color="gray40") +
  geom_point(data=lv, aes(x=longitude, y=latitude, size=n), color="gray80", pch=21, fill="red") +
  scale_fill_gradient(low="#ffffe5",high="#004529") +
  theme_void() +
  theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))
# 
# write.csv(nal, "data/2019liver_samples20220308.csv")

nal %>% group_by(Age) %>% count()
nal %>% group_by(Region) %>% count()

#### 2020 Livers ####

# subset to 2020 and non-missing towns
ages <- subset(ages, !is.na(longitude) & HarvestYear==2020)
ages <- unique(ages)

# Read livers from 2020
liv20 <- read.csv("data/2020-livers-all.csv")

# Join liver info to metadata
ages <- left_join(ages, liv20, by="RegionalID")

# Subset to animals with liver samples
dat <- subset(ages, liver==1 | RegionalID=="2020-6134" | RegionalID=="2020-9481")

# IDK why 4018 is lost, but just adding in manually here
l4018 <- data.frame(HarvestDate=NA, County="Otsego", Town="Hartwick", WMU="4F", Region="4",
                    Sex="F", Age="0.5", AgeClass="Juvenile", CC="A", AgeRange=NA, 
                    Zone="Harvest Expansion Area", TrapperID="2029-8000-0448", RegionalID="2020-4018",
                    SubmissionYear=as.numeric(2020), HarvestYear="2020", Village=NA, 
                    longitude=as.numeric(1696505), latitude=as.numeric(2371246), liver=NA)
dat <- rbind(dat, l4018)

plot(dat$longitude, dat$latitude)

dat %>% group_by(Region) %>% count()

ggplot(dat, aes(x=longitude, y=latitude, color=Region)) + geom_point() + theme_bw()

dat2 <- dat %>% group_by(Town) %>% count()
dat2 <- left_join(dat2, centlab, by="Town")

# make new SPDF to extract points
dat <- dat[,c(13:15,1:8,17:19)]

xy <- dat[c(12,13)]
dat.sp <- SpatialPointsDataFrame(xy, dat, proj4string=CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))

for.grid <- raster::extract(grid.r, dat.sp, method="simple") 

dat$PctForest <- for.grid
# ggplot(dat, aes(x=longitude, y=latitude, color=PctForest)) + geom_point()

hist(dat$PctForest)
# quantile(dat$PctForest, probs=c(0.25, 0.50, 0.75, 1), na.rm=TRUE)

dat$ForestGroup <- 1
dat$ForestGroup[dat$PctForest>48 & dat$PctForest<=69] <- 2
dat$ForestGroup[dat$PctForest>69 & dat$PctForest<=91.8] <- 3
dat$ForestGroup[dat$PctForest>91.8 & dat$PctForest<=100] <- 4

dat$AG <- ifelse(dat$Age < 3.5, "J", "A")
dat <- unite(dat, "AgeSex", c(17,9), sep="", remove=FALSE)

dat %>% group_by(Age) %>% count()

N.order <- order(dat$HarvestDate, decreasing=FALSE)
dat <- dat[N.order,]

# write.csv(dat, "data/2020_forest_groups.csv")

# Select livers
livers <- data.frame()
for (i in 1:4) {
  temp <- dat[dat$ForestGroup==i,]
  
  # Save adults
  adults <- temp[temp$AG=="A",]
  
  # Pull out regions 3 and 4 because there are so few
  reg3 <- temp[(temp$Region==3 | temp$Region==4),]
  
  # Combine adults & regions 3/4, determine how many and any overlap.
  temp2 <- rbind(adults, reg3)
  
  # How many of each age class already selected
  af <- nrow(temp2[temp2$AgeSex=="AF",])
  am <- nrow(temp2[temp2$AgeSex=="AM",])
  jf <- nrow(temp2[temp2$AgeSex=="JF",])
  jm <- nrow(temp2[temp2$AgeSex=="JM",])
  
  n <- af+am+jf+jm
  f <- af+jf
  m <- am+jm
  
  temp.f <- temp[temp$AG=="J" & temp$Sex=="F" & temp$Region>4,]
  temp.m <- temp[temp$AG=="J" & temp$Sex=="M" & temp$Region>4,]
  
  add.f <- temp.f[sample(nrow(temp.f), size=15-f, replace=FALSE),]
  add.m <- temp.m[sample(nrow(temp.m), size=15-m, replace=FALSE),]
  
  # Add to livers
  livers <- rbind(livers, temp2, add.f, add.m)
  
}

lv <- livers %>% group_by(Town) %>% count()
lv <- left_join(lv, ages, by="Town")

livers %>% group_by(Region) %>% count()
ggplot() + 
  geom_polygon(data=nys, aes(x=long, y=lat, group=group), fill="gray40") +
  geom_polygon(data=for.ha2, aes(x=long, y=lat, group=group, fill=Forest), color="gray40") +
  geom_point(data=lv, aes(x=longitude, y=latitude, size=n), color="gray80", pch=21, fill="red") +
  scale_fill_gradient(low="#ffffe5",high="#004529") +
  theme_void() +
  theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))

# write.csv(livers, "data/liver_samples20211207.csv")

livers <- read.csv("data/liver_samples20211207.csv")

# Switch out replaced livers
new <- dat[dat$RegionalID=="2020-5051" | dat$RegionalID=="2020-5143" |
             dat$RegionalID=="2020-6134" | dat$RegionalID=="2020-6247",]

livers <- livers[!(livers$RegionalID=="2020-5066" | livers$RegionalID=="2020-5116" |
                     livers$RegionalID=="2020-6143" | livers$RegionalID=="2020-6231"),]

livers <- rbind(livers, new)
# write.csv(livers, "data/liver_samples20211207-updated.csv")

ggplot() + 
  geom_polygon(data=nys, aes(x=long, y=lat, group=group), fill="gray40") +
  geom_polygon(data=for.ha2, aes(x=long, y=lat, group=group, fill=Forest), color="gray40") +
  geom_point(data=lv18, aes(x=longitude, y=latitude, size=n), color="gray80", pch=21, fill="red") +
  scale_fill_gradient(low="#ffffe5",high="#004529") +
  theme_void() +
  theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))

## Update after going through freezer

reg7 <- subset(dat, Region==7)
reg7 <- subset(reg7, RegionalID!="70001" & 
                 RegionalID!="70073" & 
                 RegionalID!="70085" & 
                 RegionalID!="70097" & 
                 RegionalID!="70100" & 
                 RegionalID!="70170" & 
                 RegionalID!="70179" & 
                 RegionalID!="70180" & 
                 RegionalID!="70183" & 
                 RegionalID!="7972" & 
                 RegionalID!="7988" & ForestGroup==3)
sample(reg7$RegionalID, size=3, replace=FALSE)

revis <- read.csv("data/2020livers_IDs-only.csv")

revis <- left_join(revis, dat, by="RegionalID")
lv20 <- revis %>% group_by(Town) %>% count()
lv20 <- left_join(lv20, centlab, by="Town")

revis %>% group_by(ForestGroup) %>% count()

ggplot() + 
  geom_polygon(data=nys, aes(x=long, y=lat, group=group), fill="gray40") +
  geom_polygon(data=for.ha2, aes(x=long, y=lat, group=group, fill=Forest), color="gray40") +
  geom_point(data=lv20, aes(x=longitude, y=latitude, size=n), color="gray80", pch=21, fill="red") +
  scale_fill_gradient(low="#ffffe5",high="#004529") +
  theme_void() +
  theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))

# write.csv(revis, "data/2020livers_revised.csv")


