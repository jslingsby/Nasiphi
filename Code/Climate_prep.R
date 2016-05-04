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
if(Sys.getenv("USER")=="jasper") {datwd <- "/Users/jasper/Documents/GIS/VegToolsRaw/Rasters/"; climwd <- "/Users/jasper/GIT/Nyasha/Data/Adam/"}

##########################################
###2) Get and process data
##########################################

###Get climate data
c <- brick(paste(climwd,"currentclimate.grd", sep = ""))
f <- brick(paste(climwd,"futureclimate.grd", sep = ""))

###Get a raster with planning units based on the protected areas
pa <- raster(paste(datwd, "PA_npaes_test_web.tif", sep=""))
proj4string(pa) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
vm <- raster(paste(datwd, "VEG12test_web.tif", sep=""))
proj4string(vm) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

###Reproject protected areas and climate data to the same projection, crop, and rasterize to the same grid size
proj4string(c)
proj4string(f)
proj4string(pa)
#pa, f and c are in different projections...

#extract layers to use in calculation
rsh <- c("tmax01", "tmin07", "mmp01") #historical
his <- c[[which(names(c)%in%rsh)]]
rsf <- c("RCP45_bcccsm11_20812100_tmax01.1", "RCP45_bcccsm11_20462065_tmin07.1", "RCP45_bcccsm11_20812100_mmp01.1") #future
fut <- f[[which(names(f)%in%rsf)]]

#trim protected areas to the same extent as the climate data
cE <- extent(projectRaster(his, crs = CRS(proj4string(pa)))) #Get the extent of c when reprojected to the same coord system as pa
pa <- crop(pa, cE) #Crop pa to the extent of c

#reproject and rasterize historical climate to 30m and write out
projectRaster(his, pa, filename = "Data/historical", bylayer = T, suffix = names(his)) 

#reproject and rasterize historical climate to 30m and write out
projectRaster(fut, pa, filename = "Data/future", bylayer = T, suffix = names(fut)) 

x <- raster(list.files("Data/", filename = "Data/future", bylayer = T, suffix = names(his), bylayer=T, format="ascii")

############################################################
##### 3) Climate stability calculation
############################################################

