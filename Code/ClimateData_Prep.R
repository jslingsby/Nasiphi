######################################################
######## R code for preparing climate data for hull 
######## volume calculations 
#######################################################
######## Compiled by Nasiphi Ntshanga 2016
######## Last edited: 9 June 2016
#######################################################

##########################################
###1) Get libraries and setwd
##########################################
library(rgdal)
library(raster)

if(Sys.getenv("USERNAME")=="nasi") {climwd <- "C:/Users/nasip/Dropbox/Academics/PhD/Data/ClimateData/";datwd <- "C:/Users/nasip/Dropbox/Academics/PhD/Data/Rasters/"}
if(Sys.getenv("USER")=="jasper") {datwd <- "/Users/jasper/Documents/GIS/VegToolsRaw/Rasters/"; climwd <- "/Users/jasper/GIT/Nyasha/Data/Adam/"}

##########################################
###2) Get and process planning unit data
##########################################

##### planning units (natural and transformed)
fragments <- readOGR("C:/Users/nasip/Dropbox/Academics/PhD/Data/von Hase et al/Lowland_/Remnants/Remnant_Shape", layer = "remnants_only_wtm21")
proj4string(fragments) <- "+proj=tmerc +lat_0=0 +lon_0=21 +ellps=WGS84 +datum=WGS84 +units=m"
cape <- readOGR("C:/Users/nasip/Dropbox/Academics/PhD/Data/Remnants_CAPE2001_Pence2016beta/CAPE_Natural_Remnants", layer = "cape_untransformed_areas_genf15m_gw")
proj4string(cape) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #EPSG:4326
cape <- spTransform(cape,lc@crs)

veg <- raster("Data/westcoast.asc")
proj4string(veg) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

veg <- raster("Data/westUT.asc")


#vegetation
vm <- raster(paste(datwd, "VEG12test_web.tif", sep=""))
proj4string(vm) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

#transform cape to landcover
fragments <- spTransform(fragments,lc@crs)

#set extent based on cape extent
lc <- mask(lc, fragments, filename="Data/lcs.grd", inverse=FALSE, 
     updatevalue=NA, updateNA=FALSE, overwrite=T)
#west coast 
westcoast <-raster("Data/wc_reno.asc")
proj4string(westcoast) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
westcoast <- aggregate(westcoast, 3, fun = "max", na.rm=T) #Aggregate pa into a 90m raster
writeRaster(westcoast, "Data/westcoast.asc", overwrite=T)
westcoast <-raster("Data/westcoast.asc")
proj4string(westcoast) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
his <- brick(paste(climwd,"currentclimate.grd", sep = ""))
projectRaster(his, westcoast, filename = "Data/historical", bylayer = T, overwrite = T, suffix = names(his), format="GTiff")
future <- brick(paste(climwd,"climate_future.grd", sep = ""))
projectRaster(future, westcoast, filename = "Data/future", bylayer = T, overwrite = T, suffix = names(future)) #if you write as .tif you lose layer names 

future <-brick("Data/future.grd", sep="")

nex <- extent(1986566, 2866431, -4141382, -3632780) #CAPE extent
lc <- crop(lc, nex) #Crop pa to chosen extent

if(!file.exists("Data/lcs.grd")) {
  NAvalue(lc) <- 128
  NAvalue(lc) <- 0
  tab <-lc@data@attributes[[1]]
  tb <- tab[,c("Value","LC3V")]
  tb <- rbind(tb, c(128,"NA"))
  tb <- rbind(tb, c(0,"NA"))
  subs(lc,tab, by="ID",which="LC3V",subsWithNA=TRUE, filename= "Data/lcs", overwrite=T)
  lc <- raster("Data/lcs.grd")
} else {lc <- raster("Data/lcs.grd")}

### veg map
veg <- raster(paste(datwd, "VEG12test_web.tif", sep=""))
proj4string(veg) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
veg <- crop(veg, nex)

### transformed veg map
if(!file.exists("Data/vegmap.asc")) {
  transveg <- mask(veg, lc, filename="Data/wlc", maskvalue= , overwrite=TRUE)
} else {transveg <- raster("Data/vmt.grd")}

###aggregate planning units to 90m and write out the rasters
#veg map
veg <- aggregate(veg, 3, fun = "max", na.rm=T) #Aggregate pa into a 90m raster
writeRaster(veg, "Data/vegmap.asc", overwrite=T)
transveg <- aggregate(transveg, 3, fun = "max", na.rm=T) #Aggregate pa into a 90m raster
writeRaster(transveg, "Data/transveg.asc", overwrite=T)
rm(lc)
###############################################

#########################################
###3) Get and process climate data 
#########################################

his <- brick(paste(climwd,"currentclimate.grd", sep = ""))
futureSD <- brick(paste(climwd,"futureclimate.grd", sep = ""))
info <- read.csv("C:/Users/nasip/Dropbox/Academics/PhD/Data/ClimateData/futureclimate_info.csv", sep = ",")

## calculate actual future projections per climatic variable
#(load rasters developed by Jasper)....layer names mising climate variable
#futureplus <-brick(paste(climwd,"futureage_plus_Current.grd", sep = ""))
#futuremap <-brick(paste(climwd,"map/futureclimate_plus_Current.grd", sep = ""))
#futuremmp <-brick(paste(climwd,"mmp01/futureclimate_plus_Current.grd", sep = ""))
#futurepptconc <-brick(paste(climwd,"pptconc/futureclimate_plus_Current.grd", sep = ""))
#futuretmax <-brick(paste(climwd,"tmax01/futureclimate_plus_Current.grd", sep = ""))
#futuretmin <-brick(paste(climwd,"tmin07/futureclimate_plus_Current.grd", sep = ""))

if(!file.exists("Data/future.grd")) {
  ###extract climate layers to use in raster for calculation
  #historical
  histmax <- his[[which(names(his)%in%"tmax01")]]
  histmin <- his[[which(names(his)%in%"tmin07")]]
  hismmp <- his[[which(names(his)%in%"mmp01")]]
  hismap <- his[[which(names(his)%in%"map")]]
  hispptconc <- his[[which(names(his)%in%"pptconc")]]
  ########
#future [take all the tmax's from all models]
  ftmax <- c("RCP45_bcccsm11_20462065_tmax01.1","RCP45_bcccsm11_20812100_tmax01.1","RCP45_BNUESM_20462065_tmax01.1","RCP45_BNUESM_20812100_tmax01.1","RCP45_CanESM2_20462065_tmax01.1","RCP45_CanESM2_20812100_tmax01.1","RCP45_CNRMCM5_20462065_tmax01.1","RCP45_CNRMCM5_20812100_tmax01.1","RCP45_FGOALSs2_20462065_tmax01.1","RCP45_FGOALSs2_20812100_tmax01.1","RCP45_GFDLESM2G_20462065_tmax01.1","RCP45_GFDLESM2G_20812100_tmax01.1","RCP45_GFDLESM2M_20462065_tmax01.1","RCP45_GFDLESM2M_20812100_tmax01.1","RCP45_MIROC5_20462065_tmax01.1","RCP45_MIROC5_20812100_tmax01.1","RCP45_MIROCESM_20462065_tmax01.1","RCP45_MIROCESM_20812100_tmax01.1","RCP45_MIROCESMCHEM_20462065_tmax01.1","RCP45_MIROCESMCHEM_20812100_tmax01.1","RCP45_MRICGCM3_20462065_tmax01.1","RCP45_MRICGCM3_20812100_tmax01.1","RCP85_bcccsm11_20462065_tmax01.1","RCP85_bcccsm11_20812100_tmax01.1","RCP85_BNUESM_20462065_tmax01.1","RCP85_BNUESM_20812100_tmax01.1","RCP85_CanESM2_20462065_tmax01.1","RCP85_CanESM2_20812100_tmax01.1","RCP85_CNRMCM5_20462065_tmax01.1","RCP85_CNRMCM5_20812100_tmax01.1","RCP85_FGOALSs2_20462065_tmax01.1","RCP85_FGOALSs2_20812100_tmax01.1","RCP85_GFDLESM2G_20462065_tmax01.1","RCP85_GFDLESM2G_20812100_tmax01.1","RCP85_GFDLESM2M_20462065_tmax01.1","RCP85_GFDLESM2M_20812100_tmax01.1","RCP85_MIROC5_20462065_tmax01.1","RCP85_MIROC5_20812100_tmax01.1","RCP85_MIROCESM_20462065_tmax01.1","RCP85_MIROCESM_20812100_tmax01.1","RCP85_MIROCESMCHEM_20462065_tmax01.1","RCP85_MIROCESMCHEM_20812100_tmax01.1","RCP85_MRICGCM3_20462065_tmax01.1","RCP85_MRICGCM3_20812100_tmax01.1")
  
  ftmin <- c("RCP45_bcccsm11_20462065_tmin07.1","RCP45_bcccsm11_20812100_tmin07.1","RCP45_BNUESM_20462065_tmin07.1","RCP45_BNUESM_20812100_tmin07.1","RCP45_CNRMCM5_20462065_tmin07.1","RCP45_CNRMCM5_20812100_tmin07.1","RCP45_FGOALSs2_20462065_tmin07.1","RCP45_FGOALSs2_20812100_tmin07.1","RCP45_GFDLESM2G_20462065_tmin07.1","RCP45_GFDLESM2G_20812100_tmin07.1","RCP45_GFDLESM2M_20462065_tmin07.1","RCP45_GFDLESM2M_20812100_tmin07.1","RCP45_MIROC5_20462065_tmin07.1","RCP45_MIROC5_20812100_tmin07.1","RCP45_MIROCESM_20462065_tmin07.1","RCP45_MIROCESM_20812100_tmin07.1","RCP45_MIROCESMCHEM_20462065_tmin07.1","RCP45_MIROCESMCHEM_20812100_tmin07.1","RCP45_MRICGCM3_20462065_tmin07.1","RCP45_MRICGCM3_20812100_tmin07.1","RCP85_bcccsm11_20462065_tmin07.1","RCP85_bcccsm11_20812100_tmin07.1","RCP85_BNUESM_20462065_tmin07.1","RCP85_BNUESM_20812100_tmin07.1","RCP85_CanESM2_20462065_tmin07.1","RCP85_CanESM2_20812100_tmin07.1","RCP85_CNRMCM5_20462065_tmin07.1","RCP85_CNRMCM5_20812100_tmin07.1","RCP85_FGOALSs2_20462065_tmin07.1","RCP85_FGOALSs2_20812100_tmin07.1","RCP85_GFDLESM2G_20462065_tmin07.1","RCP85_GFDLESM2G_20812100_tmin07.1","RCP85_GFDLESM2M_20462065_tmin07.1","RCP85_GFDLESM2M_20812100_tmin07.1","RCP85_MIROC5_20462065_tmin07.1","RCP85_MIROC5_20812100_tmin07.1","RCP85_MIROCESM_20462065_tmin07.1","RCP85_MIROCESM_20812100_tmin07.1","RCP85_MIROCESMCHEM_20462065_tmin07.1","RCP85_MIROCESMCHEM_20812100_tmin07.1","RCP85_MRICGCM3_20462065_tmin07.1","RCP85_MRICGCM3_20812100_tmin07.1","RCP45_CanESM2_20462065_tmin07.1","RCP45_CanESM2_20812100_tmin07.1")
  
  fmap <- c("RCP45_bcccsm11_20462065_map.1","RCP45_bcccsm11_20812100_map.1","RCP45_BNUESM_20462065_map.1","RCP45_BNUESM_20812100_map.1","RCP45_CanESM2_20462065_map.1","RCP45_CanESM2_20812100_map.1","RCP45_CNRMCM5_20462065_map.1","RCP45_CNRMCM5_20812100_map.1","RCP45_FGOALSs2_20462065_map.1","RCP45_FGOALSs2_20812100_map.1","RCP45_GFDLESM2G_20462065_map.1","RCP45_GFDLESM2G_20812100_map.1","RCP45_GFDLESM2M_20462065_map.1","RCP45_GFDLESM2M_20812100_map.1","RCP45_MIROC5_20462065_map.1","RCP45_MIROC5_20812100_map.1","RCP45_MIROCESM_20462065_map.1","RCP45_MIROCESM_20812100_map.1","RCP45_MIROCESMCHEM_20462065_map.1","RCP45_MIROCESMCHEM_20812100_map.1","RCP45_MRICGCM3_20462065_map.1","RCP45_MRICGCM3_20812100_map.1","RCP85_bcccsm11_20462065_map.1","RCP85_bcccsm11_20812100_map.1","RCP85_BNUESM_20462065_map.1","RCP85_BNUESM_20812100_map.1","RCP85_CanESM2_20462065_map.1","RCP85_CanESM2_20812100_map.1","RCP85_CNRMCM5_20462065_map.1","RCP85_CNRMCM5_20812100_map.1","RCP85_FGOALSs2_20462065_map.1","RCP85_FGOALSs2_20812100_map.1","RCP85_GFDLESM2G_20462065_map.1","RCP85_GFDLESM2G_20812100_map.1","RCP85_GFDLESM2M_20462065_map.1","RCP85_GFDLESM2M_20812100_map.1","RCP85_MIROC5_20462065_map.1","RCP85_MIROC5_20812100_map.1","RCP85_MIROCESM_20462065_map.1","RCP85_MIROCESM_20812100_map.1","RCP85_MIROCESMCHEM_20462065_map.1","RCP85_MIROCESMCHEM_20812100_map.1","RCP85_MRICGCM3_20462065_map.1","RCP85_MRICGCM3_20812100_map.1")
  
  fmmp <- c("RCP45_bcccsm11_20462065_mmp01.1","RCP45_bcccsm11_20812100_mmp01.1","RCP45_BNUESM_20462065_mmp01.1","RCP45_BNUESM_20812100_mmp01.1","RCP45_CanESM2_20462065_mmp01.1","RCP45_CanESM2_20812100_mmp01.1","RCP45_CNRMCM5_20462065_mmp01.1","RCP45_CNRMCM5_20812100_mmp01.1","RCP45_FGOALSs2_20462065_mmp01.1","RCP45_FGOALSs2_20812100_mmp01.1","RCP45_GFDLESM2G_20462065_mmp01.1","RCP45_GFDLESM2G_20812100_mmp01.1","RCP45_GFDLESM2M_20462065_mmp01.1","RCP45_GFDLESM2M_20812100_mmp01.1","RCP45_MIROC5_20462065_mmp01.1","RCP45_MIROC5_20812100_mmp01.1","RCP45_MIROCESM_20462065_mmp01.1","RCP45_MIROCESM_20812100_mmp01.1","RCP45_MIROCESMCHEM_20462065_mmp01.1","RCP45_MIROCESMCHEM_20812100_mmp01.1","RCP45_MRICGCM3_20462065_mmp01.1","RCP45_MRICGCM3_20812100_mmp01.1","RCP85_bcccsm11_20462065_mmp01.1","RCP85_bcccsm11_20812100_mmp01.1","RCP85_BNUESM_20462065_mmp01.1","RCP85_BNUESM_20812100_mmp01.1","RCP85_CanESM2_20462065_mmp01.1","RCP85_CanESM2_20812100_mmp01.1","RCP85_CNRMCM5_20462065_mmp01.1","RCP85_CNRMCM5_20812100_mmp01.1","RCP85_FGOALSs2_20462065_mmp01.1","RCP85_FGOALSs2_20812100_mmp01.1","RCP85_GFDLESM2G_20462065_mmp01.1","RCP85_GFDLESM2G_20812100_mmp01.1","RCP85_GFDLESM2M_20462065_mmp01.1","RCP85_GFDLESM2M_20812100_mmp01.1","RCP85_MIROC5_20462065_mmp01.1","RCP85_MIROC5_20812100_mmp01.1","RCP85_MIROCESM_20462065_mmp01.1","RCP85_MIROCESM_20812100_mmp01.1","RCP85_MIROCESMCHEM_20462065_mmp01.1","RCP85_MIROCESMCHEM_20812100_mmp01.1","RCP85_MRICGCM3_20462065_mmp01.1","RCP85_MRICGCM3_20812100_mmp01.1")
  
  fpptconc <- c("RCP45_bcccsm11_20462065_pptconc.1","RCP45_bcccsm11_20812100_pptconc.1","RCP45_BNUESM_20462065_pptconc.1","RCP45_BNUESM_20812100_pptconc.1","RCP45_CanESM2_20462065_pptconc.1","RCP45_CanESM2_20812100_pptconc.1","RCP45_CNRMCM5_20462065_pptconc.1","RCP45_CNRMCM5_20812100_pptconc.1","RCP45_FGOALSs2_20462065_pptconc.1","RCP45_FGOALSs2_20812100_pptconc.1", "RCP45_GFDLESM2G_20462065_pptconc.1","RCP45_GFDLESM2G_20812100_pptconc.1","RCP45_GFDLESM2M_20462065_pptconc.1","RCP45_GFDLESM2M_20812100_pptconc.1","RCP45_MIROC5_20462065_pptconc.1","RCP45_MIROC5_20812100_pptconc.1","RCP45_MIROCESM_20462065_pptconc.1","RCP45_MIROCESM_20812100_pptconc.1","RCP45_MIROCESMCHEM_20462065_pptconc.1","RCP45_MIROCESMCHEM_20812100_pptconc.1","RCP45_MRICGCM3_20462065_pptconc.1","RCP45_MRICGCM3_20812100_pptconc.1","RCP85_bcccsm11_20462065_pptconc.1","RCP85_bcccsm11_20812100_pptconc.1","RCP85_BNUESM_20462065_pptconc.1","RCP85_BNUESM_20812100_pptconc.1","RCP85_CanESM2_20462065_pptconc.1","RCP85_CanESM2_20812100_pptconc.1","RCP85_CNRMCM5_20462065_pptconc.1","RCP85_CNRMCM5_20812100_pptconc.1","RCP85_FGOALSs2_20462065_pptconc.1","RCP85_FGOALSs2_20812100_pptconc.1","RCP85_GFDLESM2G_20462065_pptconc.1","RCP85_GFDLESM2G_20812100_pptconc.1","RCP85_GFDLESM2M_20462065_pptconc.1","RCP85_GFDLESM2M_20812100_pptconc.1","RCP85_MIROC5_20462065_pptconc.1","RCP85_MIROC5_20812100_pptconc.1","RCP85_MIROCESM_20462065_pptconc.1","RCP85_MIROCESM_20812100_pptconc.1","RCP85_MIROCESMCHEM_20462065_pptconc.1","RCP85_MIROCESMCHEM_20812100_pptconc.1","RCP85_MRICGCM3_20462065_pptconc.1","RCP85_MRICGCM3_20812100_pptconc.1")
  ######
#create raster stack for each variable 
  futmax <-futureSD[[which(names(futureSD)%in%ftmax)]]   
  futmin <-futureSD[[which(names(futureSD)%in%ftmin)]]
  futppt <-futureSD[[which(names(futureSD)%in%fpptconc)]]
  futmap <-futureSD[[which(names(futureSD)%in%fmap)]]
  futmmp <-futureSD[[which(names(futureSD)%in%fmmp)]]
  x <-nlayers(futmap)
  for (i in 1:x){ 
    plot(futmap[[i]]>0)
  }
  rm(f)
  rm(fmap)
  rm(fmmp)
  rm(fpptconc)
  rm(ftmax)
  rm(ftmin)
  
  ####calculate future value by his+fut
  future_tmax <- futmax + histmax
  future_tmin <- futmin + histmin
  future_map <- futmap + hismap
  future_mmp <- futmmp + hismmp
  future_pptconc <- futppt + hispptconc
  ##name layers
  names(future_tmax)<-ftmax
  names(future_tmin)<-ftmin
  names(future_map)<-fmap
  names(future_mmp)<-fmmp
  names(future_pptconc)<-fpptconc
  
  
  rm(futmax)
  rm(futmin)
  rm(futmap)
  rm(futmmp)
  rm(futppt)
  
  
  #stack or brick into single raster and write out
  future <- addLayer(future_pptconc,future_mmp,future_map,future_tmin,future_tmax)
  writeRaster(future_pptconc, "Data/futurepptconc.grd")
  writeRaster(future_mmp, "Data/futuremmp.grd")
  writeRaster(future_map, "Data/futuremap.grd")
  writeRaster(future_tmin, "Data/futuretmin.grd")
  writeRaster(future_tmax, "Data/futuretmax.grd")
  rm(future_map)
  rm(future_mmp)
  rm(future_pptconc)
  rm(future_tmax)
  rm(future_tmin)
  rm(hismap)
  rm(hismmp)
  rm(hispptconc)
  rm(histmax)
  rm(histmin)
  rm(fmap)
  rm(fmmp)
  rm(fpptconc)
  rm(ftmax)
  rm(ftmin)
} else {future <- brick(paste(climwd,"climate_future.grd", sep = ""))}
rm(futureSD)

#check raster maths
xx <-nlayers(future)
for (i in 1:xx){ 
  plot(future[[i]]>0)
  }

#reproject climate data to planning unit data
########

#reproject and rasterize historical climate to 30m(90m) and write out
projectRaster(his, veg, filename = "Data/historical", bylayer = T, overwrite = T, suffix = names(his), format="GTiff") 

#reproject and rasterize future climate to 30m and write out
projectRaster(future, veg, filename = "Data/future", bylayer = T, overwrite = T, suffix = names(future)) #if you write as .tif you lose layer names 

#read in files and write out as asci rasters
future <- brick("Data/future.grd")

#########################################################
##### 4) Create raster stack for each future projection
########################################################

#extract each projection
###############################

#RCP45
bcc65 <- c("RCP45_bcccsm11_20462065_map.1","RCP45_bcccsm11_20462065_mmp01.1","RCP45_bcccsm11_20462065_pptconc.1","RCP45_bcccsm11_20462065_tmax01.1","RCP45_bcccsm11_20462065_tmin07.1")

bcc100 <- c("RCP45_bcccsm11_20812100_pptconc.1","RCP45_bcccsm11_20812100_mmp01.1","RCP45_bcccsm11_20812100_map.1","RCP45_bcccsm11_20812100_tmin07.1","RCP45_bcccsm11_20812100_tmax01.1")

BNUESM65 <- c("RCP45_BNUESM_20462065_map.1","RCP45_BNUESM_20462065_mmp01.1","RCP45_BNUESM_20462065_pptconc.1","RCP45_BNUESM_20462065_tmax01.1","RCP45_BNUESM_20462065_tmin07.1")

BNUESM100 <- c("RCP45_BNUESM_20812100_map.1","RCP45_BNUESM_20812100_mmp01.1","RCP45_BNUESM_20812100_pptconc.1","RCP45_BNUESM_20812100_tmax01.1","RCP45_BNUESM_20812100_tmin07.1")

CanESM65 <- c("RCP45_CanESM2_20462065_map.1","RCP45_CanESM2_20462065_mmp01.1","RCP45_CanESM2_20462065_pptconc.1","RCP45_CanESM2_20462065_tmax01.1","RCP45_CanESM2_20462065_tmin07.1")

CanESM100 <- c("RCP45_CanESM2_20812100_map.1","RCP45_CanESM2_20812100_mmp01.1","RCP45_CanESM2_20812100_pptconc.1","RCP45_CanESM2_20812100_tmax01.1","RCP45_CanESM2_20812100_tmin07.1")

CNRMCM65 <- c("RCP45_CNRMCM5_20462065_map.1","RCP45_CNRMCM5_20462065_mmp01.1","RCP45_CNRMCM5_20462065_pptconc.1","RCP45_CNRMCM5_20462065_tmax01.1","RCP45_CNRMCM5_20462065_tmin07.1")

CNRMCM100 <- c("RCP45_CNRMCM5_20812100_map.1","RCP45_CNRMCM5_20812100_mmp01.1","RCP45_CNRMCM5_20812100_pptconc.1","RCP45_CNRMCM5_20812100_tmax01.1","RCP45_CNRMCM5_20812100_tmin07.1")

FGOALSs65 <- c("RCP45_FGOALSs2_20462065_map.1","RCP45_FGOALSs2_20462065_mmp01.1","RCP45_FGOALSs2_20462065_pptconc.1","RCP45_FGOALSs2_20462065_tmax01.1","RCP45_FGOALSs2_20462065_tmin07.1")

FGOALSs100 <- c("RCP45_FGOALSs2_20812100_map.1","RCP45_FGOALSs2_20812100_mmp01.1","RCP45_FGOALSs2_20812100_pptconc.1","RCP45_FGOALSs2_20812100_tmax01.1","RCP45_FGOALSs2_20812100_tmin07.1")

GFDLESM65 <- c("RCP45_GFDLESM2G_20462065_map.1", "RCP45_GFDLESM2G_20462065_mmp01.1","RCP45_GFDLESM2G_20462065_pptconc.1","RCP45_GFDLESM2G_20462065_tmax01.1","RCP45_GFDLESM2G_20462065_tmin07.1")

GFDLESM100 <- c("RCP45_GFDLESM2G_20812100_map.1","RCP45_GFDLESM2G_20812100_mmp01.1","RCP45_GFDLESM2G_20812100_pptconc.1","RCP45_GFDLESM2G_20812100_tmax01.1","RCP45_GFDLESM2G_20812100_tmin07.1")

GFDL65 <- c("RCP45_GFDLESM2M_20462065_map.1","RCP45_GFDLESM2M_20462065_mmp01.1","RCP45_GFDLESM2M_20462065_pptconc.1","RCP45_GFDLESM2M_20462065_tmax01.1","RCP45_GFDLESM2M_20462065_tmin07.1")

GFDL100 <- c("RCP45_GFDLESM2M_20812100_map.1","RCP45_GFDLESM2M_20812100_mmp01.1","RCP45_GFDLESM2M_20812100_pptconc.1","RCP45_GFDLESM2M_20812100_tmax01.1","RCP45_GFDLESM2M_20812100_tmin07.1")

MIROC65 <- c("RCP45_MIROC5_20462065_map.1","RCP45_MIROC5_20462065_mmp01.1","RCP45_MIROC5_20462065_pptconc.1","RCP45_MIROC5_20462065_tmax01.1","RCP45_MIROC5_20462065_tmin07.1")

MIROC100 <- c("RCP45_MIROC5_20812100_map.1","RCP45_MIROC5_20812100_mmp01.1","RCP45_MIROC5_20812100_pptconc.1","RCP45_MIROC5_20812100_tmax01.1","RCP45_MIROC5_20812100_tmin07.1")

MIROCES65 <- c("RCP45_MIROCESM_20462065_map.1","RCP45_MIROCESM_20462065_mmp01.1","RCP45_MIROCESM_20462065_pptconc.1","RCP45_MIROCESM_20462065_tmax01.1","RCP45_MIROCESM_20462065_tmin07.1")

MIROCES100 <- c("RCP45_MIROCESM_20812100_map.1","RCP45_MIROCESM_20812100_mmp01.1","RCP45_MIROCESM_20812100_pptconc.1","RCP45_MIROCESM_20812100_tmax01.1","RCP45_MIROCESM_20812100_tmin07.1")

MIROCHEM65 <- c("RCP45_MIROCESMCHEM_20462065_map.1","RCP45_MIROCESMCHEM_20462065_mmp01.1","RCP45_MIROCESMCHEM_20462065_pptconc.1","RCP45_MIROCESMCHEM_20462065_tmax01.1","RCP45_MIROCESMCHEM_20462065_tmin07.1")

MIROCHEM100 <- c("RCP45_MIROCESMCHEM_20812100_map.1","RCP45_MIROCESMCHEM_20812100_mmp01.1","RCP45_MIROCESMCHEM_20812100_pptconc.1","RCP45_MIROCESMCHEM_20812100_tmax01.1","RCP45_MIROCESMCHEM_20812100_tmin07.1")

MRIC65 <- c("RCP45_MRICGCM3_20462065_map.1","RCP45_MRICGCM3_20462065_mmp01.1","RCP45_MRICGCM3_20462065_pptconc.1","RCP45_MRICGCM3_20462065_tmax01.1","RCP45_MRICGCM3_20462065_tmin07.1")

MRIC100 <- c("RCP45_MRICGCM3_20812100_map.1","RCP45_MRICGCM3_20812100_mmp01.1","RCP45_MRICGCM3_20812100_pptconc.1","RCP45_MRICGCM3_20812100_tmax01.1","RCP45_MRICGCM3_20812100_tmin07.1")

bcc100 <- future[[which(names(future)%in%bcc100)]]
bcc65 <-future[[which(names(future)%in%bcc65)]]
BNUESM100 <-future[[which(names(future)%in%BNUESM100)]]
BNUESM65 <-future[[which(names(future)%in%BNUESM65)]]
CanESM100 <-future[[which(names(future)%in%CanESM100)]]
CanESM65 <-future[[which(names(future)%in%CanESM65)]]
CNRMCM100 <-future[[which(names(future)%in%CNRMCM100)]]
CNRMCM65 <-future[[which(names(future)%in%CNRMCM65)]]
FGOALSs100 <-future[[which(names(future)%in%FGOALSs100)]]
FGOALSs65 <-future[[which(names(future)%in%FGOALSs65)]]
GFDL100 <-future[[which(names(future)%in%GFDL100)]]
GFDL65 <-future[[which(names(future)%in%GFDL65)]]
GFDLESM100 <-future[[which(names(future)%in%GFDLESM100)]]
GFDLESM65 <-future[[which(names(future)%in%GFDLESM65)]]
MIROCHEM100 <-future[[which(names(future)%in%MIROCHEM100)]]
MIROCHEM65 <- future[[which(names(future)%in%MIROCHEM65)]]
MIROC65 <- future[[which(names(future)%in%MIROC65)]]
MIROC100 <- future[[which(names(future)%in%MIROC100)]]
MIROCES100 <- future[[which(names(future)%in%MIROCES100)]]
MIROCES65 <- future[[which(names(future)%in%MIROCES65)]]
MRIC65 <- future[[which(names(future)%in%MRIC65)]]
MRIC100 <- future[[which(names(future)%in%MRIC100)]]
#########################################################

#write out raster for each projection

writeRaster(bcc65, overwrite=T,"Data/bcc65.grd")
writeRaster(bcc100, overwrite=T,"Data/bcc100.grd")
writeRaster(BNUESM65, overwrite=T,"Data/BNUESM65.grd")
writeRaster(BNUESM100, overwrite=T,"Data/BNUESM100.grd")
writeRaster(CanESM65, overwrite=T,"Data/CanESM65.grd")
writeRaster(CanESM100, overwrite=T,"Data/CanESM100.grd")
writeRaster(CNRMCM65, overwrite=T,"Data/CNRMCM65.grd")
writeRaster(CNRMCM100, overwrite=T,"Data/CNRMCM100.grd")
writeRaster(FGOALSs65, overwrite=T,"Data/FGOALSs65.grd")
writeRaster(FGOALSs100, overwrite=T,"Data/FGOALSs100.grd")
writeRaster(GFDL100, overwrite=T,"Data/GFDL100.grd")
writeRaster(GFDL65, overwrite=T,"Data/GFDL65.grd")
writeRaster(GFDLESM65, overwrite=T,"Data/GFDLESM65.grd")
writeRaster(GFDLESM100, overwrite=T,"Data/GFDLESM100.grd")
writeRaster(MIROC65, overwrite=T,"Data/MIROC65.grd")
writeRaster(MIROC100, overwrite=T,"Data/MIROC100.grd")
writeRaster(MIROCES65, overwrite=T,"Data/MIROCES65.grd")
writeRaster(MIROCES100, overwrite=T,"Data/MIROCES100.grd")
writeRaster(MIROCHEM65, overwrite=T,"Data/MIROCHEM65.grd")
writeRaster(MIROCHEM100, overwrite=T,"Data/MIROCHEM100.grd")
writeRaster(MRIC65, overwrite=T,"Data/MRIC65.grd")
writeRaster(MRIC100,overwrite=T,"Data/MRIC100.grd")

###################################################################

#RCP85
bcc65 <- c("RCP85_bcccsm11_20462065_map.1","RCP85_bcccsm11_20462065_mmp01.1","RCP85_bcccsm11_20462065_pptconc.1","RCP85_bcccsm11_20462065_tmax01.1","RCP85_bcccsm11_20462065_tmin07.1")

bcc100 <- c("RCP85_bcccsm11_20812100_pptconc.1","RCP85_bcccsm11_20812100_mmp01.1","RCP85_bcccsm11_20812100_map.1","RCP85_bcccsm11_20812100_tmin07.1","RCP85_bcccsm11_20812100_tmax01.1")

BNUESM65 <- c("RCP85_BNUESM_20462065_map.1","RCP85_BNUESM_20462065_mmp01.1","RCP85_BNUESM_20462065_pptconc.1","RCP85_BNUESM_20462065_tmax01.1","RCP85_BNUESM_20462065_tmin07.1")

BNUESM100 <- c("RCP85_BNUESM_20812100_map.1","RCP85_BNUESM_20812100_mmp01.1","RCP85_BNUESM_20812100_pptconc.1","RCP85_BNUESM_20812100_tmax01.1","RCP85_BNUESM_20812100_tmin07.1")

CanESM65 <- c("RCP85_CanESM2_20462065_map.1","RCP85_CanESM2_20462065_mmp01.1","RCP85_CanESM2_20462065_pptconc.1","RCP85_CanESM2_20462065_tmax01.1","RCP85_CanESM2_20462065_tmin07.1")

CanESM100 <- c("RCP85_CanESM2_20812100_map.1","RCP85_CanESM2_20812100_mmp01.1","RCP85_CanESM2_20812100_pptconc.1","RCP85_CanESM2_20812100_tmax01.1","RCP85_CanESM2_20812100_tmin07.1")

CNRMCM65 <- c("RCP85_CNRMCM5_20462065_map.1","RCP85_CNRMCM5_20462065_mmp01.1","RCP85_CNRMCM5_20462065_pptconc.1","RCP85_CNRMCM5_20462065_tmax01.1","RCP85_CNRMCM5_20462065_tmin07.1")

CNRMCM100 <- c("RCP85_CNRMCM5_20812100_map.1","RCP85_CNRMCM5_20812100_mmp01.1","RCP85_CNRMCM5_20812100_pptconc.1","RCP85_CNRMCM5_20812100_tmax01.1","RCP85_CNRMCM5_20812100_tmin07.1")

FGOALSs65 <- c("RCP85_FGOALSs2_20462065_map.1","RCP85_FGOALSs2_20462065_mmp01.1","RCP85_FGOALSs2_20462065_pptconc.1","RCP85_FGOALSs2_20462065_tmax01.1","RCP85_FGOALSs2_20462065_tmin07.1")

FGOALSs100 <- c("RCP85_FGOALSs2_20812100_map.1","RCP85_FGOALSs2_20812100_mmp01.1","RCP85_FGOALSs2_20812100_pptconc.1","RCP85_FGOALSs2_20812100_tmax01.1","RCP85_FGOALSs2_20812100_tmin07.1")

GFDLESM65 <- c("RCP85_GFDLESM2G_20462065_map.1", "RCP85_GFDLESM2G_20462065_mmp01.1","RCP85_GFDLESM2G_20462065_pptconc.1","RCP85_GFDLESM2G_20462065_tmax01.1","RCP85_GFDLESM2G_20462065_tmin07.1")

GFDLESM100 <- c("RCP85_GFDLESM2G_20812100_map.1","RCP85_GFDLESM2G_20812100_mmp01.1","RCP85_GFDLESM2G_20812100_pptconc.1","RCP85_GFDLESM2G_20812100_tmax01.1","RCP85_GFDLESM2G_20812100_tmin07.1")

GFDL65 <- c("RCP85_GFDLESM2M_20462065_map.1","RCP85_GFDLESM2M_20462065_mmp01.1","RCP85_GFDLESM2M_20462065_pptconc.1","RCP85_GFDLESM2M_20462065_tmax01.1","RCP85_GFDLESM2M_20462065_tmin07.1")

GFDL100 <- c("RCP85_GFDLESM2M_20812100_map.1","RCP85_GFDLESM2M_20812100_mmp01.1","RCP85_GFDLESM2M_20812100_pptconc.1","RCP85_GFDLESM2M_20812100_tmax01.1","RCP85_GFDLESM2M_20812100_tmin07.1")

MIROC65 <- c("RCP85_MIROC5_20462065_map.1","RCP85_MIROC5_20462065_mmp01.1","RCP85_MIROC5_20462065_pptconc.1","RCP85_MIROC5_20462065_tmax01.1","RCP85_MIROC5_20462065_tmin07.1")

MIROC100 <- c("RCP85_MIROC5_20812100_map.1","RCP85_MIROC5_20812100_mmp01.1","RCP85_MIROC5_20812100_pptconc.1","RCP85_MIROC5_20812100_tmax01.1","RCP85_MIROC5_20812100_tmin07.1")

MIROCES65 <- c("RCP85_MIROCESM_20462065_map.1","RCP85_MIROCESM_20462065_mmp01.1","RCP85_MIROCESM_20462065_pptconc.1","RCP85_MIROCESM_20462065_tmax01.1","RCP85_MIROCESM_20462065_tmin07.1")

MIROCES100 <- c("RCP85_MIROCESM_20812100_map.1","RCP85_MIROCESM_20812100_mmp01.1","RCP85_MIROCESM_20812100_pptconc.1","RCP85_MIROCESM_20812100_tmax01.1","RCP85_MIROCESM_20812100_tmin07.1")

MIROCHEM65 <- c("RCP85_MIROCESMCHEM_20462065_map.1","RCP85_MIROCESMCHEM_20462065_mmp01.1","RCP85_MIROCESMCHEM_20462065_pptconc.1","RCP85_MIROCESMCHEM_20462065_tmax01.1","RCP85_MIROCESMCHEM_20462065_tmin07.1")

MIROCHEM100 <- c("RCP85_MIROCESMCHEM_20812100_map.1","RCP85_MIROCESMCHEM_20812100_mmp01.1","RCP85_MIROCESMCHEM_20812100_pptconc.1","RCP85_MIROCESMCHEM_20812100_tmax01.1","RCP85_MIROCESMCHEM_20812100_tmin07.1")

MRIC65 <- c("RCP85_MRICGCM3_20462065_map.1","RCP85_MRICGCM3_20462065_mmp01.1","RCP85_MRICGCM3_20462065_pptconc.1","RCP85_MRICGCM3_20462065_tmax01.1","RCP85_MRICGCM3_20462065_tmin07.1")

MRIC100 <- c("RCP85_MRICGCM3_20812100_map.1","RCP85_MRICGCM3_20812100_mmp01.1","RCP85_MRICGCM3_20812100_pptconc.1","RCP85_MRICGCM3_20812100_tmax01.1","RCP85_MRICGCM3_20812100_tmin07.1")

bcc100 <- future[[which(names(future)%in%bcc100)]]
bcc65 <-future[[which(names(future)%in%bcc65)]]
BNUESM100 <-future[[which(names(future)%in%BNUESM100)]]
BNUESM65 <-future[[which(names(future)%in%BNUESM65)]]
CanESM100 <-future[[which(names(future)%in%CanESM100)]]
CanESM65 <-future[[which(names(future)%in%CanESM65)]]
CNRMCM100 <-future[[which(names(future)%in%CNRMCM100)]]
CNRMCM65 <-future[[which(names(future)%in%CNRMCM65)]]
FGOALSs100 <-future[[which(names(future)%in%FGOALSs100)]]
FGOALSs65 <-future[[which(names(future)%in%FGOALSs65)]]
GFDL100 <-future[[which(names(future)%in%GFDL100)]]
GFDL65 <-future[[which(names(future)%in%GFDL65)]]
GFDLESM100 <-future[[which(names(future)%in%GFDLESM100)]]
GFDLESM65 <-future[[which(names(future)%in%GFDLESM65)]]
MIROCHEM100 <-future[[which(names(future)%in%MIROCHEM100)]]
MIROCHEM65 <- future[[which(names(future)%in%MIROCHEM65)]]
MIROC65 <- future[[which(names(future)%in%MIROC65)]]
MIROC100 <- future[[which(names(future)%in%MIROC100)]]
MIROCES100 <- future[[which(names(future)%in%MIROCES100)]]
MIROCES65 <- future[[which(names(future)%in%MIROCES65)]]
MRIC65 <- future[[which(names(future)%in%MRIC65)]]
MRIC100 <- future[[which(names(future)%in%MRIC100)]]
#########################################################

#write out raster for each projection

writeRaster(bcc65, overwrite=T,"Data/rcp85/bcc65.grd")
writeRaster(bcc100, overwrite=T,"Data/rcp85/bcc100.grd")
writeRaster(BNUESM65, overwrite=T,"Data/rcp85/BNUESM65.grd")
writeRaster(BNUESM100, overwrite=T,"Data/rcp85/BNUESM100.grd")
writeRaster(CanESM65, overwrite=T,"Data/rcp85/CanESM65.grd")
writeRaster(CanESM100, overwrite=T,"Data/rcp85/CanESM100.grd")
writeRaster(CNRMCM65, overwrite=T,"Data/rcp85/CNRMCM65.grd")
writeRaster(CNRMCM100, overwrite=T,"Data/rcp85/CNRMCM100.grd")
writeRaster(FGOALSs65, overwrite=T,"Data/rcp85/FGOALSs65.grd")
writeRaster(FGOALSs100, overwrite=T,"Data/rcp85/FGOALSs100.grd")
writeRaster(GFDL100, overwrite=T,"Data/rcp85/GFDL100.grd")
writeRaster(GFDL65, overwrite=T,"Data/rcp85/GFDL65.grd")
writeRaster(GFDLESM65, overwrite=T,"Data/rcp85/GFDLESM65.grd")
writeRaster(GFDLESM100, overwrite=T,"Data/rcp85/GFDLESM100.grd")
writeRaster(MIROC65, overwrite=T,"Data/rcp85/MIROC65.grd")
writeRaster(MIROC100, overwrite=T,"Data/rcp85/MIROC100.grd")
writeRaster(MIROCES65, overwrite=T,"Data/rcp85/MIROCES65.grd")
writeRaster(MIROCES100, overwrite=T,"Data/rcp85/MIROCES100.grd")
writeRaster(MIROCHEM65, overwrite=T,"Data/rcp85/MIROCHEM65.grd")
writeRaster(MIROCHEM100, overwrite=T,"Data/rcp85/MIROCHEM100.grd")
writeRaster(MRIC65, overwrite=T,"Data/rcp85/MRIC65.grd")
writeRaster(MRIC100,overwrite=T,"Data/rcp85/MRIC100.grd")

