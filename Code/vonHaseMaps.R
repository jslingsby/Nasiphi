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
library(raster)

##########################################
###2) Get data
##########################################

frag <- readOGR(dsn = "/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/fragments.kml", layer = "Remnants")
Obsites <- readOGR(dsn = "/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/obptsG.kml", layer = "obptsG")
Poly <- readOGR(dsn = "/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/poly.kml", layer = "My Places")

Obsites@data$Poly <- Obsites %over% Poly

/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/obptsG.kml
/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/sites.kmz
/Users/jasper/Dropbox/Shared/Nasiphi/von Hase/wcptsG.kml