##########################################
######## R code for preparing climate data for hull 
######## volume calculations 
##########################################
######## Compiled by Nasiphi Ntshanga 2016
######## Last edited: 4 May 2016
##########################################

##########################################
###1) Get libraries and setwd
##########################################
library(rgdal)
library(raster)

if(Sys.getenv("USERNAME")=="Receptionist") {climwd <- "C:/Users/Receptionist/Dropbox/Academics/PhD/Data/ClimateData/";datwd <- "C:/Users/Receptionist/Dropbox/Academics/PhD/Data/Rasters/"}
if(Sys.getenv("USER")=="jasper") {datwd <- "/Users/jasper/Documents/GIS/VegToolsRaw/Rasters/"; climwd <- "/Users/jasper/GIT/Nyasha/Data/Adam/"}

##########################################
###2) Get and process data
##########################################

###Set extent
#nex <- extent(2020000, 2140000, -4100000, -3950000)
nex <- extent(2030000, 2060000, -4080000, -4010000)

###Get climate data
his <- brick(paste(climwd,"currentclimate.grd", sep = ""))
f <- brick(paste(climwd,"futureclimate.grd", sep = ""))

###Get a raster with planning units based on the protected areas, landcover and vegetation types

#vegetation map
vm <- raster(paste(datwd, "VEG12test_web.tif", sep=""))
proj4string(vm) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

#landcover
lc <- raster(paste(datwd,"LC13test_web.tif", sep=""))
proj4string(lc) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

#trim planning units to the same extent as the climate data
#nex <- extent(projectRaster(his, crs = CRS(proj4string(pa)))) #Get the extent of c when reprojected to the same coord system as pa

lc <- crop(lc, nex) #Crop pa to chosen extent
vm <- crop(vm, nex) 

#convert landcover to 1's and 0's
NAvalue(lc) <- 128
NAvalue(lc) <- 0
tab <-lc@data@attributes[[1]]
tb <- tab[,c("Value","LC3V")]
tb <- rbind(tb, c(128,"NA"))
tb <- rbind(tb, c(0,"NA"))
subs(lc,tab, by="ID",which="LC3V",subsWithNA=TRUE, filename= "Data/lcs", overwrite=T)
lcsub <- raster(paste(datwd,"Data/lcs.grd"))

###aggregate planning units to 90m and write out the rasters
#landcover
lc <- aggregate(lcsub, 3, fun = "max", na.rm=T) #Aggregate pa into a 90m raster
writeRaster(lc, "Data/landcover.asc")
#veg map
vm <- aggregate(vm, 3, fun = "max", na.rm=T) #Aggregate pa into a 90m raster
writeRaster(vm, "Data/vegmap.asc")

#create transformed veg map
transveg <- mask(vm, lcsub, filename="Data/vmt", maskvalue= 0, overwrite=TRUE)
transveg <- aggregate(transveg, 3, fun = "max", na.rm=T) #Aggregate pa into a 90m raster
writeRaster(transveg, "Data/transveg.asc")

###extract climate layers to use in calculation
#rsh <- c("tmax01", "tmin07", "mmp01", "map","pptconc") #historical
#his <- c
rsf <- c("RCP45_bcccsm11_20462065_map.1","RCP45_bcccsm11_20462065_mmp01.1","RCP45_bcccsm11_20462065_pptconc.1","RCP45_bcccsm11_20462065_tmax01.1","RCP45_bcccsm11_20462065_tmin07.1") #future
fut <- f[[which(names(f)%in%rsf)]]

rm(c)
rm(f)

#reproject and rasterize historical climate to 30m and write out
projectRaster(his, vm, filename = "Data/historical", bylayer = T, suffix = names(his), format="GTiff") 

#reproject and rasterize historical climate to 30m and write out
projectRaster(fut, vm, filename = "Data/future", bylayer = T, suffix = names(fut), format="GTiff") 


#read in files and write out as asci rasters
x <- stack(list.files("Data", pattern="historical", full.names=T))
writeRaster(x, filename = "Data/historical", bylayer = T, suffix = names(x), format="ascii")
x <- stack(list.files("Data", pattern="future", full.names=T))
writeRaster(x, filename = "Data/future", bylayer = T, suffix = names(x), format="ascii")

############################################################
##### 3) Climate stability calculation
############################################################

