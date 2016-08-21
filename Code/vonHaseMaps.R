##########################################
######## R code for preparing von Hase field site maps
########
##########################################
######## Compiled by Jasper Slingsby 2016
######## Last edited: 18 Aug 2016
##########################################

##########################################
###1) Get libraries and setwd
##########################################
library(rgdal)
library(dismo)
library(gdata)

##########################################
###2) Get data
##########################################
#Fragments (2003)
#frag <- readOGR(dsn = "/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/fragments.kml", layer = "Remnants")
#Overberg sites
#Obsites <- readOGR(dsn = "/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/obptsG.kml", layer = "obptsG")

#Nasiphi's mapped sites
Poly <- readOGR(dsn = "/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/poly.kml", layer = "My Places")
#Poly <- spTransform(Poly, CRSobj = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "))

##Calculate polygon areas
#sapply(sapply(slot(Poly, "polygons"), function(x) sapply(slot(x,"Polygons"), slot, "area")), sum)/10000


#Get, merge and fix Overberg tabular data
Tab <- read.csv("/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/Overbergcombo.csv", stringsAsFactors = F)
row.names(Tab) <- Tab$Site

T2 <- t(read.xls("/Users/jasper/Documents/Databases/von Hase et al/Lowland_/Fieldwork/Databases_SENSITIVE/Fielddata Lowlands Project, Overberg.xls", header = F, sheet = 2, skip = 1, stringsAsFactors=F, row.names = 1)) #
Tab <- merge(Tab, T2) #Merge 2 datasets by "Site"

Tab$Lat_dd <- as.numeric(gsub(",","\\.",Tab$Lat_dd)) #Fix coordinates
Tab$Long_dd <- as.numeric(gsub(",","\\.",Tab$Long_dd))
coordinates(Tab) <- ~ Long_dd + Lat_dd #Make it a "spatial" object
proj4string(Tab) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" #set projection

#/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/obptsG.kml
#/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/sites.kmz
#/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/wcptsG.kml

##########################################
###3) Veg Map
##########################################

vm <- readOGR(dsn = "/Users/jasper/Documents/GIS/VegMap/Data/vegm2006.shp", layer = "vegm2006")
vm <- spTransform(vm, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

Tab$NatVegType <- as.character(over(Tab, vm)$NAME)

###Trim Tab to von Hase's "ren" (Renosterveld) and the National Veg Map's Western, Central and Eastern Ruens Renosterveld and remove sites that have since been transformed
FTab <- Tab[which(Tab$Veg_type1=="ren"),]
FTab <- FTab[which(FTab$NatVegType %in% c("Western Rûens Shale Renosterveld", "Central Rûens Shale Renosterveld", "Eastern Rûens Shale Renosterveld")),]
FTab <- FTab[-which(FTab$Transformed == "Yes"),]
FTab <- FTab[-which(FTab@data$Geol1 == "Silcrete"),]
FTab <- FTab[-which(FTab@data$Priority < 4),]

#writeOGR(FTab, "/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/OverbergTrimmed19Aug2016.kml", layer = "OverbergTrimmed19Aug2016", driver = "KML")

#x <- FTab[,1]
#names(x)[which(names(x)=="Site")]<-"name"
#writeOGR(x, "/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/OverbergTrimmed19Aug2016.gpx", layer="waypoints", driver="GPX")

ftab <- FTab@data

write.csv(ftab, "/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/OverbergTrimmed19Aug2016.csv")

ftabp <- ftab[,which(colnames(ftab) %in% c("Site", "Altitude..m.", "Nearest.major.locality..dir.from.", "Farm", "Exact.locality", "Slope1", "Aspect", "Geol1", "Notes", "Additional.Comments"))]


1, 24, 25, 27, 29, 30, 37:40, 50:52, 55:58, 65, 66
##########################################
###4) Plot
##########################################

#for(i in 1:nrow(Tab)) {

##Get polygon
#x <- Poly[which(as.character(Poly$Name)==Tab$Polygon[i]),]

##Get point names
#pts <- Tab[which(Tab$Polygon%in%as.character(x$Name)),]
#n <- pts$Site

##Set extent
#bb <- bbox(x)
##b <- (bb - rowMeans(bb)) * 1.5 + rowMeans(bb)

##Get Google Earth image using dismo
#gmap_cdb <- gmap(bb, exp=1.5, type='satellite', filename='', style=NULL, rgb=FALSE, lonlat=T, scale=3)

#pdf(paste("/Users/jasper/Documents/GIS/von Hase sites/", paste(n, collapse="_"),".pdf", sep=""), width=10, height=10*(nrow(gmap_cdb)/ncol(gmap_cdb)))
#plot(gmap_cdb)
#plot(x, add = T, col=NA, border="green")
#points(pts, col="springgreen3", pch=19)
#text(pts, labels = "Site", pos=1)
#scalebar(0.05, col="white", xy = coordinates(pts) + (c(bb[1,2],bb[2,1])-coordinates(pts))/1.1)
#dev.off()
#}
