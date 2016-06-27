######## Script for playing with Cape Lowlands 
######## data (von Hase et al.2003) and CAPE NATURE
##########################################
######## Compiled by Nasiphi Ntshanga 2016
######## Last edited: 25 June 2016
##########################################

##########################################
###1) Get libraries and setwd
##########################################

library(vegan)
library(simba)
#for gis analyses
library(sp)
library(raster)
library(rgdal)
library(maptools)
library(rgeos)
library(spatstat)
library(RANN)
library(geosphere)
library(FNN)

if(Sys.getenv("USERNAME")=="jasper") {giswd <- "C:/Users/jasper.SAEON/Documents/Nasiphi's/Data/Von Hase et al 2003 Cape Lowlands Report/Lowland_/";
datwd <-"C:/Users/jasper.SAEON/Documents/Nasiphi's/Data/Von Hase et al 2003 Cape Lowlands Report/Lowland_"}
if(Sys.getenv("USERNAME")=="Receptionist") {giswd <- "C:/Users/Receptionist/Dropbox/Academics/PhD/Data/von Hase et al/Lowland_/";
datwd <- "C:/Users/Receptionist/Dropbox/Academics/PhD/Data/von Hase et al/PRODUCTS/Fieldwork/"}
#if(Sys.getenv("USER")=="jasper") {datwd <- ""; giswd <- ""}
setwd("~/Nasiphi's/Data/Von Hase et al 2003 Cape Lowlands Report/Lowland_")
##########################################
###2) Get and process GIS data
##########################################


#sampled sites
wcpts <- readOGR(paste0(giswd, "Fieldwork/Shapefiles"), layer = "wc_fieldsites_wtm21")
ob_elginpts <- readOGR(paste0(giswd, "Fieldwork/Shapefiles"), layer = "ob&elgin_fieldsites_wtm21")
proj4string(wcpts) <- "+proj=tmerc +lat_0=0 +lon_0=21 +ellps=WGS84 +datum=WGS84 +units=m"
proj4string(ob_elginpts) <- "+proj=tmerc +lat_0=0 +lon_0=21 +ellps=WGS84 +datum=WGS84 +units=m"

#lowland fragments
fragments <- readOGR(paste0(giswd, "Remnants/Remnant_Shape"), layer = "remnants_only_wtm21")
proj4string(fragments) <- "+proj=tmerc +lat_0=0 +lon_0=21 +ellps=WGS84 +datum=WGS84 +units=m"

#Cape Nature data
cape <- readOGR("C:/Users/Receptionist/Dropbox/Academics/PhD/Data/Remnants_CAPE2001_Pence2016beta/CAPE_Natural_Remnants", layer = "cape_untransformed_areas_genf15m_gw")
proj4string(cape) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #EPSG:4326


#project and write out data for google earth visualisation
fragmentsG <- spTransform(fragments, CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
#capeG <- spTransform(cape, CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

#plots
wcptsG <- wcpts@data
coordinates(wcptsG) <- ~ LONG_DD + LAT_DD
proj4string(wcptsG) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
obptsG <- ob_elginpts@data
coordinates(obptsG) <- ~ LONG_DD + LAT_DD
proj4string(obptsG) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

#write out kml (google earth) and shapefile
writeOGR(fragmentsG, dsn = "Nasiphi's/fragments_plotdat.kml", layer= "fragments", driver = "KML")
writeOGR(wcptsG, dsn = "Nasiphi's/wcptsG.kml", layer = "wcptsG", driver="KML")
writeOGR(obptsG, dsn = "Nasiphi's/obptsG.kml", layer = "obptsG", driver="KML")

writeOGR(cape, dsn = "Data/CapeNature.kml", layer = "CapeNature", driver="KML")

### relate plot locaion to fragments

#calculate polygon centroid and calculate Nearest neighbour distance to points
sppcentroids <- getSpPPolygonsLabptSlots(fragmentsG)
centroids = gCentroid(fragmentsG,byid=TRUE)
geocentroid <- centroid(fragmentsG)

#calculate NND between these two objects

#NN objects need objects of equal lenght
#http://stackoverflow.com/questions/16448402/distance-of-point-feature-to-nearest-polygon-in-r
gdist <-gDistance(obptsG, fragmentsG, byid=TRUE, hausdorff=FALSE, densifyFrac = NULL) #eucl. distance
 nn2(fragments, query = ob_elginpts, k = min(5, nrow(data)), treetype = c("kd", "bd"),
    searchtype = c("standard", "priority", "radius"), radius = 0, eps = 0) #

g <- get.knnx(coordinates(fragments), coordinates(wcpts),k=3)
str(g)
plot(fragments, col=2, xlim=c(-1e5,6e5))
plot(wcpts, add=TRUE)
segments(coordinates(wcpts)[,1], coordinates(wcpts)[,2], coordinates(fragments[g$nn.index[,1]])[,1], coordinates(fragments[g$nn.index[,1]])[,2])
#fragdist <- distGeo(centroids,obptsG, a=6378137, f=1/298.257223563)


#get fragments with all plot data
#fragments <- fragments[pts, ]
#subset fragments into dif layers based on data source ('DESC')
#remove buffer 
#renosterveld = subset(fragments, DESC = "Renosterveld in domain")
#renosterveld <- fragments[fragments@data$DESC="Renosterveld in domain",]


##########################################
###2) Get and process plot data
##########################################
#Load csv data

### Boland domain 
bolandsite <- read.csv("Boland-Swartland_Habitat.csv", sep = ";", skip = 3, stringsAsFactors = T))
bolandvegdat <- read.csv("Boland-Swartland_Veg.csv", sep = ";")
bolandaliens <-read.csv("Boland-Swartland_Aliens.csv", sep = ";")
bolandspecial <-read.csv("Boland-Swartland_Special.csv", sep = ";",skip = 3)

##clean and trim
bolandsite <- bolandsite[1:137,1:43] #remove blank cells
rownames(bolandsite) <- bolandsite[,3]
write.csv(bolandsite, "Boland_sitedata.csv")
str(bolandsite)
####add spp data 

## veg dat
bolandveg <- bolandvegdat[1:13,] #create veg description data
#fix dimensions
bolandveg <- t(bolandveg) 
colnames(bolandveg) <- bolandveg[1,]
bolandveg <- as.data.frame(bolandveg[2:138,])
#bolandveg <- bb2num(bolandveg, from = c("r", "+", "1", "2", "3", "4", "5"), to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5))
write.csv(bolandveg,"Boland_veg.csv")

#aliens
bolandaliens <- bolandaliens[1:53,1:138] #remove empty rows and columns
rownames(bolandaliens) <-bolandaliens[,1] 
bolandaliens <- bolandaliens[,2:138] #remove extra row
bolandaliens <- t(bolandaliens) #transpose 
bolandaliens <- as.data.frame(bolandaliens)
#convert BB to numeric
bolandaliens <- bb2num(bolandaliens, from = c("r", "+", "1", "2", "3", "4", "5"), 
                       to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5))
bolandveg$AlienSpecies_Count <- rowSums(decostand(bolandaliens, "pa")) #add alien spp count to veg dataframe
write.csv(bolandaliens,"Boland_aliens.csv")


##spp dat
bolandspp <- bolandvegdat[14:314,]
bolandspp <- t(bolandspp)
colnames(bolandspp) <- bolandspp[1,]
bolandspp <- as.data.frame(bolandspp[2:138,])
bolandspp <- bb2num(bolandspp, from = c("r", "+", "1", "2", "3", "4", "5"), 
                    to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5))
bolandveg$PlantSpeciesCount <- rowSums(decostand(bolandspp, "pa"))
write.csv(bolandspp,"Boland_plant_spp.csv")
rm(bolandvegdat)
## special species??
rm(bolandspecial)

##### Elgin
elginsite <- read.csv("Elgin_Habitat.csv", sep = ";", skip = 4)
elginveg <- read.csv("Elgin_Veg.csv", sep = ";",skip = 1)
elginaliens <-read.csv("Elgin_Basin_Aliens.csv", sep = ",", skip = 1)
elginspp <- read.csv("Elgin_Domspp.csv", sep = ";",skip = 1)
elginspecial <-read.csv("Elgin_Specials.csv", sep = ";",skip = 1,row.names = 1)
#fix stuff 
elginsite <- elginsite[1:4,1:39]
write.csv(elginsite, "Elgin_site.csv")

#aliens
rownames(elginaliens) <- elginaliens[,1]
elginaliens <- t(elginaliens[,2:5])
elginaliens <- as.data.frame(elginaliens)
elginaliens <- bb2num(elginaliens, from = c("r", "+", "1", "2", "3", "4", "5"), 
                      to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5))
elginveg$AlienSpecies_Count <- rowSums(decostand(elginaliens, "pa")) #add alien spp count to veg dataframe
write.csv(elginaliens, "Elgin_aliens.csv")

## Overberg domain
overbergsite <- read.csv("Overberg_Habitat.csv", sep = ";", skip = 3)
overbergveg <- read.csv("Overberg_Veg.csv", sep = ",", skip = 1)
overbergaliens <-read.csv("Overberg_Aliendat.csv", sep = ",", skip = 1)
overbergspp <- read.csv("Overberg_Domspp.csv", sep = ";",skip = 1)
overbergspecial <-read.csv("Overberg_Special.csv", sep = ",",skip = 1)

#fix stuff
overbergsite <- overbergsite[1:155,]
write.csv(overbergsite, "Overberg_site.csv")
#veg
o
overberg <- t(overbergveg)
colnames(overbergveg) <- overbergveg[1,]
overbergveg <- as.data.frame(overbergveg[2:156,])
write.csv(overbergveg, "Overberg_veg.csv")
#aliens
overbergaliens <- overbergaliens[1:35,1:156]
rownames(overbergaliens) <- overbergaliens[,1]
overbergaliens <- t(overbergaliens[,2:156])
overbergaliens <- as.data.frame(overbergaliens)
overbergaliens <- bb2num(overbergaliens, from = c("r", "+", "1", "2", "3", "4", "5"), 
                         to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5))
overbergveg$AlienSpecies_Count <- rowSums(decostand(overbergaliens, "pa"))
write.csv(overbergaliens, "Overberg_aliens.csv")
#spp
overbergspp <- t(overbergspp[,2:157])
colnames(overbergspp) <- overbergspp[1,]
overbergspp <- as.data.frame(overbergspp[2:156,])
overbergspp <- bb2num(overbergspp, from = c("r", "+", "1", "2", "3", "4", "5"), 
                      to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5))
overbergveg$PlantSpecies_Count <- rowSums(decostand(overbergspp, "pa"))
#######################################################



