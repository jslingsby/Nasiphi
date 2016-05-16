## R scripts for calculating climate stability 
## requires these libraries 
library(rgdal)
library(geometry) 
library(foreign)
library(raster)

#data -- requires ascii rasters describing variables of climate space in present and future, and spatial planning units, all with the same cell size, extent, and geographic projection.

#read in data and set up dataframes for historic and future climate info

#x <- as(x, 'SpatialGridDataFrame') # to convert from raster
#x <- raster("Data/vegmap.asc")
#pu <- readGDAL("Data/vegmap.asc") #raster with numerical planning units
#pu <- readGDAL("Data/landcover.asc")
pu <- readGDAL("Data/transveg.asc")
names(pu@data)<-"pu"
cellsize <-  pu@grid@cellsize
hst <- pu
fut <- pu
planning <- unique(pu@data)
planning <- planning[which(planning$pu!="NA"),]

#read in asciis of historic/future climate data

#hst avg max temp January
hsttmax <- readGDAL("Data/historical_historical_tmax01.asc")
hst@data$tmax<-hsttmax$band1
rm(hsttmax)

#hst avg min temp July
hsttmin <- readGDAL("Data/historical_historical_tmin07.asc")
hst@data$tmin<-hsttmin$band1
rm(hsttmin)

#hst avg climatic water deficit
hstcwd <- readGDAL("Data/historical_historical_mmp01.asc")
hst@data$cwd<-hstcwd$band1
rm(hstcwd)

#hst avg mean annual precipitation
hstmap <- readGDAL("Data/historical_historical_map.asc")
hst@data$map<-hstmap$band1
rm(hstmap)

#hst preciptation concentration
hstpptconc <- readGDAL("Data/historical_historical_pptconc.asc")
hst@data$pptconc<-hstpptconc$band1
rm(hstpptconc)

#fut avg max temp Jan
futtmax <- readGDAL("Data/future_future_RCP45_bcccsm11_20462065_tmax01.1.asc")
fut@data$tmax<-futtmax$band1
rm(futtmax)

#fut avg min temp July
futtmin <- readGDAL("Data/future_future_RCP45_bcccsm11_20462065_tmin07.1.asc")
fut@data$tmin<-futtmin$band1
rm(futtmin)

#fut avg mean annual precipitation
futmap <- readGDAL("Data/future_future_RCP45_bcccsm11_20462065_map.1.asc")
fut@data$map<-futmap$band1
rm(futmap)

#fut preciptation concentration
futpptconc <- readGDAL("Data/future_future_RCP45_bcccsm11_20462065_pptconc.1.asc")
fut@data$pptconc<-futpptconc$band1
rm(futpptconc)

#fut avg climatic water deficit
futcwd <- readGDAL("Data/future_future_RCP45_bcccsm11_20462065_mmp01.1.asc")
fut@data$cwd<-futcwd@data$band1
rm(futcwd)
hstlu<-hst@data
futlu<-fut@data

# normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
both<-rbind(hstlu,futlu)
range01 <- function(x,y){(x-min(y))/(max(y)-min(y))}
hstn<-hstlu
futn<-futlu
hstn$tmax<-range01(hstlu$tmax, both$tmax)
hstn$tmin<-range01(hstlu$tmin, both$tmin)
hstn$cwd<-range01(hstlu$cwd, both$cwd)
hstn$map<-range01(hstlu$map, both$map)
hstn$pptconc<-range01(hstlu$pptconc, both$pptconc)
futn$tmax<-range01(futlu$tmax, both$tmax)
futn$tmin<-range01(futlu$tmin, both$tmin)
futn$cwd<-range01(futlu$cwd, both$cwd)
futn$map<-range01(futlu$map, both$map)
futn$pptconc<-range01(futlu$pptconc, both$pptconc)

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  #convex hull volumes
  hull.hstn <- convhulln(sub.hstn[,2:4], options="FA")
  hull.futn <- convhulln(sub.futn[,2:4], options="FA")
  hull.bothn <- convhulln(sub.bothn[,2:4], options="FA")
  
  #fill in climstab table
  climstab$hullV_hst[i]<-hull.hstn$vol
  climstab$hullV_fut[i]<-hull.futn$vol
  climstab$hullV_both[i]<-hull.bothn$vol
  
  #calculate stability 
  climstab$climstab[i] <- ((hull.hstn$vol + hull.futn$vol) - hull.bothn$vol)/hull.hstn$vol
} #end loop through planning values

#write results to a table
#write.csv(climstab, "Data/climstabtransveg.csv")
#write.csv(climstab, "Data/climstabveg.csv")
write.csv(climstab, "Data/climstabfrag.csv")
