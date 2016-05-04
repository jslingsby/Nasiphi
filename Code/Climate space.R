## R scripts for calculating climate stability 
## requires these libraries 
#library(rgdal)
#library(geometry) 
#library(foreign)

#data -- requires ascii rasters describing variables of climate space in present and future, and spatial planning units, all with the same cell size, extent, and geographic projection.

#read in data and set up dataframes for historic and future climate info

pu <- readGDAL("path") #raster with numerical planning units
names(pu@data)<-"pu"
cellsize <-  pu@grid@cellsize
hst <- pu
fut <- pu
planning <- unique(pu@data)
planning <- planning[which(planning$pu!="NA"),]

#read in asciis of historic/future climate data

#hst avg max temp June, July, August
hsttmax <- readGDAL("path")
hst@data$tmax<-hsttmax$band1
rm(hsttmax)

#hst avg min temp December, January, February
hsttmin <- readGDAL("path")
hst@data$tmin<-hsttmin$band1
rm(hsttmin)

#hst avg climatic water deficit
hstcwd <- readGDAL("path")
hst@data$cwd<-hstcwd$band1
rm(hstcwd)

#fut avg max temp June, July, August
futtmax <- readGDAL("path")
fut@data$tmax<-futtmax@data$band1
rm(futtmax)

#fut avg min temp December, January, February
futtmin <- readGDAL("path")
fut@data$tmin<-futtmin@data$band1
rm(futtmin)

#fut avg climatic water deficit
futcwd <- readGDAL("path")
fut@data$cwd<-futcwd@data$band1
rm(futcwd)
hstlu<-hst@data
futlu<-fut@data

# normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA"),]
both<-rbind(hstlu,futlu)
range01 <- function(x,y){(x-min(y))/(max(y)-min(y))}
hstn<-hstlu
futn<-futlu
hstn$tmax<-range01(hstlu$tmax, both$tmax)
hstn$tmin<-range01(hstlu$tmin, both$tmin)
hstn$cwd<-range01(hstlu$cwd, both$cwd)
futn$tmax<-range01(futlu$tmax, both$tmax)
futn$tmin<-range01(futlu$tmin, both$tmin)
futn$cwd<-range01(futlu$cwd, both$cwd)

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, 
                                                                                    length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==i),]
  sub.futn <- futn[which(futn$pu==i),]
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
write.csv(climstab, "path/climstab.csv")