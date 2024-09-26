library(tidyverse)
library(sf)
library(terra)

#### Step 1: QA/QC on fisher harvest data ####

# Read in age data
age2016 <- read.csv("data/20211111_2016_age_data.csv")
age2017 <- read.csv("data/20211111_2017_age_data.csv")
age2018 <- read.csv("data/20211111_2018_age_data.csv")
age2019 <- read.csv("data/20211111_2019_age_data.csv")
age2020 <- read.csv("data/20211111_2020_age_data.csv")

# Add empty sample ID column for 2016 & 2017
age2016$Regional.Sample.ID.. <- NA
age2017$Regional.Sample.ID.. <- NA

# Bind together and rename. Fill missing values with NA. Change date type.
ages <- rbind(age2016,age2017,age2018,age2019,age2020)
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
twn <- st_transform(twn, CRS="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
cents <- st_coordinates(twn) %>% as_tibble() %>% select(X, Y) %>% distinct()
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
  dist <- pointDistance(temp[1,2:3], temp[2,2:3], lonlat=FALSE)/1000
  
  new <- data.frame(Town=temp$Town[1], Dist_km=dist)
  dapart <- rbind(dapart, new)
}

age2 <- ages[ages$Town %in% et$Town,]

N.order <- order(age2$Town, decreasing=FALSE)
age2 <- age2[N.order,]


# Join town coordinates
ages <- left_join(ages, centlab, by="Town")

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

write_csv(ages2, "data/analysis-ready/2016-2020_ages_data.csv")