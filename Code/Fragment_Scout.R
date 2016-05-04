##############################################################################
######## Script to identify lowland fragments of interest
##############################################################################
######## Compiled by Jasper Slingsby 2015
######## Last edited: 6 Dec 2015
######## Data: 
##############################################################################

##############################################################################
###1) Get libraries
##############################################################################

library(rgdal)
library(ggmap)

##############################################################################
###2) Get data
##############################################################################

frag <- readOGR(dsn = "/Users/jasper/Documents/GIS/CFR/Remnants_Pence_2006_on_hodgepodge/vegm09_rems_deg_and_natural_gw/vegm09_rems_deg_and_natural_gw.shp", layer = "vegm09_rems_deg_and_natural_gw")

##############################################################################
###3) Plot map
##############################################################################

LA <- treedat[which(treedat$SIZE=="A" & treedat$STATE=="L"),]
bb <- bbox(trees)
b <- (bb - rowMeans(bb)) * 3 + rowMeans(bb)
mapG <- get_map(location = b, source = "google", maptype = "satellite", crop = TRUE)
MG <- ggmap(mapG) + theme(legend.position="none")
MG <- MG + geom_point(aes(x=LON, y=LAT, colour="green"), data=LA)
MG <- MG + geom_point(aes(x=coords.x1, y=coords.x2, colour="red"), data=ptsdat)
MG <- MG + geom_polygon(data=fSA, aes(x = long, y = lat, group = id, fill="red"), alpha = 0.25)

ggsave("groundtruth_v0.1.png", plot=MG)