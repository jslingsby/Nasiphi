######## Script for playing with Cape Lowlands 
######## data (von Hase et al.2003)
##########################################
######## Compiled by Nasiphi Ntshanga 2016
######## Last edited: 10 May 2016
##########################################

##########################################
###1) Get libraries and setwd
##########################################

library(vegan)
library(simba)
library(sp)
library(raster)
library(rgdal)
library(maptools)
library(rgeos)
library(RColorBrewer) #to plot coloutful polygons
library(spatstat)
library(RANN)
library(geosphere)

if(Sys.getenv("USERNAME")=="jasper") {giswd <- "C:/Users/jasper.SAEON/Documents/Nasiphi's/Data/Von Hase et al 2003 Cape Lowlands Report/Lowland_/";
datwd <-"C:/Users/jasper.SAEON/Documents/Nasiphi's/Data/Von Hase et al 2003 Cape Lowlands Report/Lowland_"}
#if(Sys.getenv("USERNAME")=="Receptionist") {giswd <- "C:/Users/Receptionist/Dropbox/Academics/PhD/Data/von Hase et al/Lowland_/";
datwd <- "C:/Users/Receptionist/Dropbox/Academics/PhD/Data/von Hase et al/PRODUCTS/Fieldwork/"}
#if(Sys.getenv("USER")=="jasper") {datwd <- ""; giswd <- ""}

##########################################
###2) Get and process GIS data
##########################################

#Load data and project it

#point data of sampled plots 
wcpts <- readOGR(paste0(giswd, "Fieldwork/Shapefiles"), layer = "wc_fieldsites_wtm21")
ob_elginpts <- readOGR(paste0(giswd, "Fieldwork/Shapefiles"), layer = "ob&elgin_fieldsites_wtm21")
proj4string(wcpts) <- "+proj=tmerc +lat_0=0 +lon_0=21 +ellps=WGS84 +datum=WGS84 +units=m"
proj4string(ob_elginpts) <- "+proj=tmerc +lat_0=0 +lon_0=21 +ellps=WGS84 +datum=WGS84 +units=m"

#fragments
fragments <- readOGR(paste0(giswd, "Remnants/Remnant_Shape"), layer = "remnants_only_wtm21")
proj4string(fragments) <- "+proj=tmerc +lat_0=0 +lon_0=21 +ellps=WGS84 +datum=WGS84 +units=m"

#####merge data into a single dataframe
x <- as.data.frame(coordinates(wcpts)) # Extract the coordinates 
z <- as.data.frame(coordinates(ob_elginpts))
ALLpts <- rbind(x,z)
rm(z)
rm(x)
coordinates(ALLpts) <- ~ coords.x1 + coords.x2 # Set coordinates 
proj4string(ALLpts) <- "+proj=tmerc +lat_0=0 +lon_0=21 +ellps=WGS84 +datum=WGS84 +units=m"


#calculate polygon centroid and calculate Nearest neighbour distance to points
wcptsG <- spTransform(wcpts, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
obptsG <- spTransform(ob_elginpts, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
frags <- spTransform(fragments, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
w <- as.data.frame(coordinates(wcpts)) # Extract the coordinates 
o <- as.data.frame(coordinates(ob_elginpts))
coordinates(w) <- ~ coords.x1 + coords.x2 # Set coordinates
coordinates(o) <- ~ coords.x1 + coords.x2
proj4string(w) <- "+proj=tmerc +lat_0=0 +lon_0=21 +ellps=WGS84 +datum=WGS84 +units=m"
proj4string(o) <- "+proj=tmerc +lat_0=0 +lon_0=21 +ellps=WGS84 +datum=WGS84 +units=m"
o<- spTransform(o, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
w<- spTransform(w, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))


sppcentroids <- getSpPPolygonsLabptSlots(frags)
centroids = gCentroid(frags,byid=TRUE)
geocentroid <- centroid(frags)

#calculate NND between these two objects
#NN objects need objects of equal lenght
gdist <-gDistance(o, frags, byid=FALSE, hausdorff=FALSE, densifyFrac = NULL) #eucl. distance
nn2(frags, query = o, k = min(5, nrow(data)), treetype = c("kd", "bd"),
    searchtype = c("standard", "priority", "radius"), radius = 0, eps = 0) #

#get co ordinates of fragments (after unprojecting ) and of points

fragdist <- distGeo(centroids,o, a=6378137, f=1/298.257223563)


#get fragments with all plot data
#fragments <- fragments[pts, ]
#subset fragments into dif layers based on data source ('DESC')
#remove buffer 
#renosterveld = subset(fragments, DESC = "Renosterveld in domain")
#renosterveld <- fragments[fragments@data$DESC="Renosterveld in domain",]

#project and write out data for google earth visualisation
fragments <- spTransform(fragments, CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
#wc points
wcptsG <- wcpts@data
coordinates(wcptsG) <- ~ LONG_DD + LAT_DD
proj4string(wcptsG) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

obptsG <- ob_elginpts@data
coordinates(obptsG) <- ~ LONG_DD + LAT_DD
proj4string(obptsG) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

#write out kml (google earth) and shapefile
writeOGR(fragments, dsn = "Nasiphi's/fragments_plotdat.kml", layer= "fragments", driver = "KML")
writeOGR(wcptsG, dsn = "Nasiphi's/wcptsG.kml", layer = "wcptsG", driver="KML")
writeOGR(obptsG, dsn = "Nasiphi's/obptsG.kml", layer = "obptsG", driver="KML")
##########################################
###2) Get and process plot data
##########################################
