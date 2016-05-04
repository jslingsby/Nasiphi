##########################################
######## R code for preparing climate data for hull 
######## volume calculations 
##########################################
######## Compiled by Nasiphi Ntshanga 2016
######## Last edited: 1 May 2016
##########################################

##########################################
###1) Get libraries and setwd
##########################################
library(rgdal)
library(raster)

if(Sys.getenv("USERNAME")=="Receptionist") {datwd <- "C:/Users/Receptionist/Dropbox/Academics/PhD/Data/ClimateData/"}
if(Sys.getenv("USER")=="jasper") {datwd <- ""}

##########################################
###2) Get and process data
##########################################

c <- brick(paste(datwd,"currentclimate.grd", sep = ""))
f <- brick(paste(datwd,"futureclimate.grd", sep = ""))


#extract layers to use in calculation
#c[[which(names(c)%in%c("map","tmax01"))]]
tmax <- raster(c,layer="tmax01")
tmin <- raster(c,layer="tmin07")
mmp <- raster(c,layer="mmp01")
tmaxf <- raster(f,layer="RCP45_bcccsm11_20812100_tmax01.1")
tminf <- raster(f,layer="RCP45_bcccsm11_20462065_tmin07.1")
mmpf <- raster(f,layer="RCP45_bcccsm11_20812100_mmp01.1")

#write ascii grid for analysis
writeRaster(tmax, filename = "tmax.asc", format="ascii")
writeRaster(tmin, filename = "tmin.asc", format="ascii")
writeRaster(mmp, filename = "mmp.asc", format="ascii")
writeRaster(mmpf, filename = "mmpf.asc", format="ascii")
writeRaster(tmaxf, filename = "tmaxf.asc", format="ascii")
writeRaster(tminf, filename = "tminf.asc", format="ascii")

#need raster with planning units
 
pa <- raster("C:/Users/Receptionist/Dropbox/Academics/PhD/Data/Rasters/PA_npaes_test_web.tif")
proj4string(pa) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
vm <- raster("C:/Users/Receptionist/Dropbox/Academics/PhD/Data/Rasters/VEG12test_web.tif")
proj4string(vm) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

#reproject 
proj4string(c)
proj4string(pa)
#pa and c are in different projections...
paC <- projectRaster(pa, crs = CRS(proj4string(c)))

############################################################
##### 3) Climate stability calculation
############################################################

