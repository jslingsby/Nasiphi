##########################################
######## R code for preparing climate data for hull 
######## volume calculations 
##########################################
######## Compiled by Nasiphi Ntshanga 2016
######## Last edited: 5 May 2016
##########################################
########## 1) Get libraries and setwd
########## 2) Prepare planning units data
########## 3) Prepare climate data

##########################################
###1) Get libraries and setwd
##########################################
library(rgdal)
library(raster)

if(Sys.getenv("USERNAME")=="Receptionist") {climwd <- "C:/Users/Receptionist/Dropbox/Academics/PhD/Data/ClimateData/";datwd <- "C:/Users/Receptionist/Dropbox/Academics/PhD/Data/Rasters/"}
if(Sys.getenv("USERNAME")=="jasper") {datwd <- "C:/Users/jasper/Documents/GIS/VegToolsRaw/Rasters/"; climwd <- "/Users/jasper/GIT/Nyasha/Data/Adam/"}
if(Sys.getenv("USER")=="jasper") {datwd <- ""; climwd ""}
##########################################
###2) Prepare planning units data
##########################################

###Get a raster with planning units based on the protected areas, landcover and vegetation types

#vegetation map
vm <- raster(paste(datwd, "VEG12test_web.tif", sep=""))
proj4string(vm) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

#landcover
lc <- raster(paste(datwd,"LC13test_web.tif", sep=""))
proj4string(lc) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

###Set extent
nex <- extent(2030000, 2060000, -4080000, -4010000) #nex <- extent(2020000, 2140000, -4100000, -3950000)

#trim planning units to specified extent
lc <- crop(lc, nex)
vm <- crop(vm, nex)

###convert landcover to 1's and 0's
NAvalue(lc) <- 128
NAvalue(lc) <- 0
tab <-lc@data@attributes[[1]]
tb <- tab[,c("Value","LC3V")]
tb <- rbind(tb, c(128,"NA"))
tb <- rbind(tb, c(0,"NA"))

if(!file.exists("Data/Data/lcs.grd")) {
  subs(lc,tab, by="ID",which="LC3V",subsWithNA=TRUE, filename= "Data/lcs", overwrite=T)
} else {lcsub <- raster("Data/lcs.grd")}

#create a transformed veg map
if(!file.exists("Data/vmt.grd")) {
  transveg <- mask(vm, lcsub, filename="Data/vmt", maskvalue= 0, overwrite=TRUE)
} else {transveg <- raster("Data/vmt.grd")}

###aggregate planning units to 90m and write out the rasters

#landcover
lc <- aggregate(lcsub, 3, fun = "max", na.rm=T) 
writeRaster(lc, "Data/landcover.asc", overwrite=TRUE)

#veg map
vm <- aggregate(vm, 3, fun = "max", na.rm=T) 
writeRaster(vm, "Data/vegmap.asc", overwrite=TRUE)

#transformed veg map
transveg <- aggregate(transveg, 3, fun = "max", na.rm=T) 
writeRaster(transveg, "Data/transveg.asc", overwrite=TRUE)

###############################################
###### 2) Prepare climate data ################
################################################

###Get climate data
his <- brick(paste(climwd,"currentclimate.grd", sep = ""))
f <- brick(paste(climwd,"futureclimate.grd", sep = ""))

###extract climate layers to use in calculation
rsf <- c("RCP45_bcccsm11_20462065_map.1","RCP45_bcccsm11_20462065_mmp01.1","RCP45_bcccsm11_20462065_pptconc.1","RCP45_bcccsm11_20462065_tmax01.1","RCP45_bcccsm11_20462065_tmin07.1") #future
fut <- f[[which(names(f)%in%rsf)]]

rm(f)

### reproject and rasterize climate data to 30m and write out
projectRaster(his, vm, filename = "Data/historical", bylayer = T, suffix = names(his), format="GTiff",overwrite=TRUE) 
projectRaster(fut, vm, filename = "Data/future", bylayer = T, suffix = names(fut), format="GTiff",overwrite=TRUE) 

#read in files and write out as asci rasters
x <- stack(list.files("Data", pattern="historical", full.names=T))
writeRaster(x, filename = "Data/historical", bylayer = T, suffix = names(x), format="ascii")
x <- stack(list.files("Data", pattern="future", full.names=T))
writeRaster(x, filename = "Data/future", bylayer = T, suffix = names(x), format="ascii")
###
