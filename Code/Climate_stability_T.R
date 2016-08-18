####################################################
############## Climate stability calculation #######
############# using TRANSFORMED vegetation #######
############## Heller et al. 2015 ##################
############## last edited 10 June 2016 ############
####################################################
####### Steps
####### 1) Load libraries and set wd
####### 2) Load and normalize data
####### 3) Hull volume calculation
####### 4) Repeat steps for all climate models
##########################################################

#### 1) load libraries
library(rgdal)
library(geometry) 
library(foreign)
library(raster)

#### 2) Load and normalize data, one scenario at a time

### Planning unit data

#natural veg map
#pu <- readGDAL("Data/vegmap.asc") #raster with numerical planning units
pu <-readGDAL("Data/westUT.asc")
names(pu@data)<-"pu"
cellsize <-  pu@grid@cellsize
hst <- pu
fut <- pu
planning <- unique(pu@data)
planning <- planning[which(planning$pu!="NA"),]


### Historical climate data
hstpptconc <- raster("Data/historical_pptconc.tif")
hstcwd <- raster("Data/historical_mmp01.tif")
hstmap <- raster("Data/historical_map.tif")
hsttmax <- raster("Data/historical_tmax01.tif")
hsttmin <- raster("Data/historical_tmin07.tif")

### Future climate data
#RCP45
###############
### bcc65
#########
bcc65 <- stack("Data/bcc65.grd")
#pptconc
futpptconc <- bcc65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)

#mmp(cwd)

futcwd <- bcc65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd)
hist(futcwd_norm, add=T,col="grey")
hist(hstcwd)
hist(hstcwd_norm, add=T,col="grey")

#map

futmap <- bcc65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/(maxmap - minmap)
futmap_norm <- (futmap - minmap)/(maxmap - minmap)

#tmin

futtmin <- bcc65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax

futtmax <- bcc65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)

### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)

##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  

  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
    {climstab[i,2:5] <- NA}
  else {
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
  }
  } #end loop through planning values

#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegbcc65.csv")
rm(bcc65)
#########################################

#Repeat
#############################
#############
#bcc100
#############
bcc100 <-stack("Data/bcc100.grd")
#pptconc
futpptconc <- bcc100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)

#mmp(cwd)
futcwd <- bcc100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")

#map
futmap <- bcc100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- bcc100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- bcc100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)

### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegbcc100.csv")
rm(bcc100)
#########################################
##############
#BNUESM100
###############
BNUESM100 <-stack ("Data/BNUESM100.grd")
#pptconc
futpptconc <- BNUESM100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- BNUESM100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- BNUESM100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- BNUESM100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- BNUESM100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegBNUESM100.csv")
rm(BNUESM100)
#########################################
##############
#BNUESM65
##########
BNUESM65 <- stack("Data/BNUESM65.grd")
#pptconc
futpptconc <- BNUESM65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- BNUESM65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- BNUESM65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- BNUESM65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- BNUESM65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)

### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegBNUESM65.csv")
rm(BNUESM65)
#########################################
##############
#CanESM100
##########
CanESM100 <- stack("Data/CanESM100.grd")
#pptconc
futpptconc <- CanESM100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- CanESM100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- CanESM100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- CanESM100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- CanESM100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegCanESM100.csv")
rm(CanESM100)
#########################################
##############
#CanESM65
##########
CanESM65 <-stack("Data/CanESM65.grd")
#pptconc
futpptconc <- CanESM65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- CanESM65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- CanESM65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- CanESM65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- CanESM65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegCanESM65.csv")
rm(CanESM65)
#########################################
##############
#CNRMCM100
##########
CNRMCM100 <- stack("Data/CNRMCM100.grd")
#pptconc
futpptconc <- CNRMCM100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- CNRMCM100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- CNRMCM100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- CNRMCM100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- CNRMCM100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegCNRMCM100.csv")
rm(CNRMCM100)
#########################################
##############
#CNRMCM65
##########
CNRMCM65 <- stack("Data/CNRMCM65.grd")
#pptconc
futpptconc <- CNRMCM65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- CNRMCM65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- CNRMCM65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- CNRMCM65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- CNRMCM65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values

#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegCNRMCM65.csv")
rm(CNRMCM65)
#########################################
##############
#FGOALSs100
##########
FGOALSs100 <- stack("Data/FGOALSs100.grd")
#pptconc
futpptconc <- FGOALSs100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- FGOALSs100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- FGOALSs100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- FGOALSs100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- FGOALSs100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values



#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegFGOALSs100.csv")
rm(FGOALSs100)
###################################################
##############
#FGOALSs65
##########
FGOALSs65 <- stack("Data/FGOALSs65.grd")
#pptconc
futpptconc <- FGOALSs65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- FGOALSs65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- FGOALSs65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- FGOALSs65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- FGOALSs65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegFGOALSs65.csv")
rm(FGOALSs65)
#########################################################
##############
#GFDL65
##########
GFDL65 <- stack("Data/GFDL65.grd")
#pptconc
futpptconc <- GFDL65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- GFDL65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- GFDL65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- GFDL65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- GFDL65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegGFDL65.csv")
rm(GFDL65)
#######################################################
##############
#GFDL100
##########
GFDL100 <- stack("Data/GFDL100.grd")
#pptconc
futpptconc <- GFDL100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- GFDL100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- GFDL100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- GFDL100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- GFDL100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegGFDL100.csv")
rm(GFDL100)
##############
#GFDLESM100
##########
GFDLESM100 <- stack("Data/GFDLESM100.grd")
#pptconc
futpptconc <- GFDLESM100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- GFDLESM100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- GFDLESM100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- GFDLESM100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- GFDLESM100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegGFDLESM100.csv")
rm(GFDLESM100)
#######################################################
##############
#GFDLESM65
##########
GFDLESM65 <- stack("Data/GFDLESM65.grd")
#pptconc
futpptconc <- GFDLESM65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- GFDLESM65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- GFDLESM65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- GFDLESM65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- GFDLESM65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegGFDLESM65.csv")
rm(GFDLESM65)
#######################################################
##############
#MIROC65
##########
MIROC65 <- stack("Data/MIROC65.grd")
#pptconc
futpptconc <- MIROC65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MIROC65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MIROC65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MIROC65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MIROC65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegMIROC65.csv")
rm(MIROC65)
#######################################################
##############
#MIROC100
##########
MIROC100 <- stack("Data/MIROC100.grd")
#pptconc
futpptconc <- MIROC100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MIROC100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MIROC100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MIROC100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MIROC100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegMIROC100.csv")
rm(MIROC100)
##############
#MIROCES100
##########
MIROCES100 <- stack("Data/MIROCES100.grd")
#pptconc
futpptconc <- MIROCES100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MIROCES100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MIROCES100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MIROCES100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MIROCES100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegMIROCES100.csv")
rm(MIROCES100)
#######################################################
##############
#MIROCES65
##########
MIROCES65 <- stack("Data/MIROCES65.grd")
#pptconc
futpptconc <- MIROCES65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MIROCES65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MIROCES65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MIROCES65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MIROCES65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegMIROCES65.csv")
rm(MIROCES65)
#######################################################
##############
#MIROCHEM65
##########
MIROCHEM65 <- stack("Data/MIROCHEM65.grd")
#pptconc
futpptconc <- MIROCHEM65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MIROCHEM65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MIROCHEM65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MIROCHEM65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MIROCHEM65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegMIROCHEM65.csv")
rm(MIROCHEM65)
#######################################################
##############
#MIROCHEM100
##########
MIROCHEM100 <-stack("Data/MIROCHEM100.grd")
#pptconc
futpptconc <- MIROCHEM100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MIROCHEM100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MIROCHEM100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MIROCHEM100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MIROCHEM100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegMIROCHEM100.csv")
rm(MIROCHEM100)
#######################################################
##############
#MRIC100
##########
MRIC100 <- stack("Data/MRIC100.grd")
#pptconc
futpptconc <- MRIC100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MRIC100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MRIC100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MRIC100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MRIC100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegMRIC100.csv")
rm(MRIC100)
#######################################################
##############
#MRIC65
##########
MRIC65 <- stack("Data/MRIC65.grd")
#pptconc
futpptconc <- MRIC65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MRIC65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MRIC65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MRIC65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MRIC65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp45/climstabvegMRIC65.csv")
rm(MRIC65)
#######################################################

#RCP85
###############
### bcc65
#########
bcc65 <- stack("Data/rcp85/bcc65.grd")
#pptconc
futpptconc <- bcc65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)

#mmp(cwd)

futcwd <- bcc65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")

#map

futmap <- bcc65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/(maxmap - minmap)
futmap_norm <- (futmap - minmap)/(maxmap - minmap)

#tmin

futtmin <- bcc65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax

futtmax <- bcc65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)

### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)

##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegbcc65.csv")
rm(bcc65)
#########################################

#Repeat
#############################
#############
#bcc100
#############
bcc100 <-stack("Data/rcp85/bcc100.grd")
#pptconc
futpptconc <- bcc100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)

#mmp(cwd)
futcwd <- bcc100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")

#map
futmap <- bcc100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- bcc100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- bcc100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)

### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegbcc100.csv")
rm(bcc100)
#########################################
##############
#BNUESM100
###############
BNUESM100 <-stack ("Data/rcp85/BNUESM100.grd")
#pptconc
futpptconc <- BNUESM100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- BNUESM100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- BNUESM100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- BNUESM100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- BNUESM100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegBNUESM100.csv")
rm(BNUESM100)
#########################################
##############
#BNUESM65
##########
BNUESM65 <- stack("Data/rcp85/BNUESM65.grd")
#pptconc
futpptconc <- BNUESM65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- BNUESM65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- BNUESM65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- BNUESM65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- BNUESM65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)

### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegBNUESM65.csv")
rm(BNUESM65)
#########################################
##############
#CanESM100
##########
CanESM100 <- stack("Data/CanESM100.grd")
#pptconc
futpptconc <- CanESM100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- CanESM100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- CanESM100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- CanESM100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- CanESM100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegCanESM100.csv")
rm(CanESM100)
#########################################
##############
#CanESM65
##########
CanESM65 <-stack("Data/rcp85/CanESM65.grd")
#pptconc
futpptconc <- CanESM65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- CanESM65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- CanESM65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- CanESM65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- CanESM65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegCanESM65.csv")
rm(CanESM65)
#########################################
##############
#CNRMCM100
##########
CNRMCM100 <- stack("Data/rcp85/CNRMCM100.grd")
#pptconc
futpptconc <- CNRMCM100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- CNRMCM100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- CNRMCM100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- CNRMCM100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- CNRMCM100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegCNRMCM100.csv")
rm(CNRMCM100)
#########################################
##############
#CNRMCM65
##########
CNRMCM65 <- stack("Data/rcp85/CNRMCM65.grd")
#pptconc
futpptconc <- CNRMCM65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- CNRMCM65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")

#map
futmap <- CNRMCM65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- CNRMCM65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- CNRMCM65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegCNRMCM65.csv")
rm(CNRMCM65)
#########################################
##############
#FGOALSs100
##########
FGOALSs100 <- stack("Data/rcp85/FGOALSs100.grd")
#pptconc
futpptconc <- FGOALSs100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- FGOALSs100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- FGOALSs100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- FGOALSs100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- FGOALSs100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegFGOALSs100.csv")
rm(FGOALSs100)
###################################################
##############
#FGOALSs65
##########
FGOALSs65 <- stack("Data/rcp85/FGOALSs65.grd")
#pptconc
futpptconc <- FGOALSs65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- FGOALSs65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- FGOALSs65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- FGOALSs65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- FGOALSs65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

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
write.csv(climstab, "Output/T/rcp85/climstabvegFGOALSs65.csv")
rm(FGOALSs65)
#########################################################
##############
#GFDL65
##########
GFDL65 <- stack("Data/rcp85/GFDL65.grd")
#pptconc
futpptconc <- GFDL65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- GFDL65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")

#map
futmap <- GFDL65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- GFDL65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- GFDL65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

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
write.csv(climstab, "Output/T/rcp85/climstabvegGFDL65.csv")
rm(GFDL65)
#######################################################
##############
#GFDL100
##########
GFDL100 <- stack("Data/rcp85/GFDL100.grd")
#pptconc
futpptconc <- GFDL100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- GFDL100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- GFDL100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- GFDL100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- GFDL100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegGFDL100.csv")
rm(GFDL100)
##############
#GFDLESM65
##########
GFDLESM65 <- stack("Data/rcp85/GFDLESM65.grd")
#pptconc
futpptconc <- GFDLESM65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- GFDLESM65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- GFDLESM65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- GFDLESM65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- GFDLESM65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegGFDLESM65.csv")
rm(GFDLESM65)
#######################################################
##############
#GFDLESM100
##########
GFDLESM100 <- stack("Data/rcp85/GFDLESM100.grd")
#pptconc
futpptconc <- GFDLESM100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- GFDLESM100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- GFDLESM100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- GFDLESM100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- GFDLESM100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values

#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegGFDLESM100.csv")
rm(GFDLESM100)
#######################################################
##############
#MIROC65
##########
MIROC65 <- stack("Data/rcp85/MIROC65.grd")
#pptconc
futpptconc <- MIROC65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MIROC65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MIROC65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MIROC65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MIROC65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegMIROC65.csv")
rm(MIROC65)
#######################################################
#######################################################
##############
#MIROC100
##########
MIROC100 <- stack("Data/rcp85/MIROC100.grd")
#pptconc
futpptconc <- MIROC100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MIROC100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MIROC100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MIROC100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MIROC100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegMIROC100.csv")
rm(MIROC100)
##############
#MIROCES100
##########
MIROCES100 <- stack("Data/rcp85/MIROCES100.grd")
#pptconc
futpptconc <- MIROCES100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MIROCES100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MIROCES100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MIROCES100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MIROCES100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegMIROCES100.csv")
rm(MIROCES100)
#######################################################
##############
#MIROCES65
##########
MIROCES65 <- stack("Data/rcp85/MIROCES65.grd")
#pptconc
futpptconc <- MIROCES65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MIROCES65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MIROCES65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MIROCES65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MIROCES65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegMIROCES65.csv")
rm(MIROCES65)
#######################################################
##############
#MIROCHEM65
##########
MIROCHEM65 <- stack("Data/rcp85/MIROCHEM65.grd")
#pptconc
futpptconc <- MIROCHEM65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MIROCHEM65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MIROCHEM65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MIROCHEM65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MIROCHEM65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegMIROCHEM65.csv")
rm(MIROCHEM65)
#######################################################
##############
#MIROCHEM100
##########
MIROCHEM100 <-stack("Data/rcp85/MIROCHEM100.grd")
#pptconc
futpptconc <- MIROCHEM100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MIROCHEM100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MIROCHEM100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MIROCHEM100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MIROCHEM100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegMIROCHEM100.csv")
rm(MIROCHEM100)
#######################################################
##############
#MRIC100
##########
MRIC100 <- stack("Data/rcp85/MRIC100.grd")
#pptconc
futpptconc <- MRIC100@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MRIC100@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MRIC100@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MRIC100@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MRIC100@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegMRIC100.csv")
rm(MRIC100)
#######################################################
##############
#MRIC65
##########
MRIC65 <- stack("Data/rcp85/MRIC65.grd")
#pptconc
futpptconc <- MRIC65@layers[[1]]
maxpptconc <- max(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
minpptconc <- min(c(getValues(hstpptconc), getValues(futpptconc)), na.rm = T)
hstpptconc_norm <- (hstpptconc - minpptconc)/ (maxpptconc - minpptconc)
futpptconc_norm <- (futpptconc - minpptconc)/ (maxpptconc - minpptconc)
#mmp(cwd)
futcwd <- MRIC65@layers[[2]]
maxcwd <- max(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
mincwd <- min(c(getValues(hstcwd), getValues(futcwd)), na.rm = T)
hstcwd_norm <- (hstcwd - mincwd)/ (maxcwd - mincwd)
futcwd_norm <- (futcwd - mincwd)/ (maxcwd - mincwd)
hist(futcwd_norm)
hist(hstcwd_norm, add=T,col="grey")
#map
futmap <- MRIC65@layers[[3]]
maxmap <- max(c(getValues(hstmap), getValues(futmap)), na.rm = T)
minmap <- min(c(getValues(hstmap), getValues(futmap)), na.rm = T)
hstmap_norm <- (hstmap - minmap)/ (maxmap - minmap)
futmap_norm <- (futmap - minmap)/ (maxmap - minmap)

#tmin
futtmin <- MRIC65@layers[[4]]
maxtmn <- max(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
mintmn  <- min(c(getValues(hsttmin), getValues(futtmin)), na.rm = T)
hsttmin_norm <- (hsttmin - mintmn)/(maxtmn - mintmn)
futtmin_norm <- (futtmin - mintmn)/(maxtmn - mintmn)

#tmax
futtmax <- MRIC65@layers[[5]]
maxtmx <- max(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
mintmx <- min(c(getValues(hsttmax), getValues(futtmax)), na.rm = T)
hsttmax_norm <- (hsttmax - mintmx)/(maxtmx - mintmx)
futtmax_norm <- (futtmax - mintmx)/(maxtmx - mintmx)
### add to spatial grid dataframe
hst@data$tmax <- hsttmax_norm@data@values
hst@data$tmin <- hsttmin_norm@data@values
hst@data$cwd <- hstcwd_norm@data@values
hst@data$map <- hstmap_norm@data@values
hst@data$pptconc <- hstpptconc_norm@data@values

fut@data$tmax <- futtmax_norm@data@values
fut@data$tmin <- futtmin_norm@data@values
fut@data$cwd <- futcwd_norm@data@values
fut@data$map <- futmap_norm@data@values
fut@data$pptconc <- futpptconc_norm@data@values

hstlu<-hst@data
futlu<-fut@data

#remove unused objects
##########
#
rm(futcwd_norm)
rm(futcwd)
rm(futmap)
rm(futmap_norm)
rm(futpptconc)
rm(futpptconc_norm)
rm(futtmax)
rm(futtmax_norm)
rm(futtmin_norm)
rm(futtmin)
rm(hstcwd_norm)
rm(hstmap_norm)
rm(hsttmin_norm)
rm(hstpptconc_norm)
rm(hsttmax_norm)
##### 3) Hull volume calculation 
################
#normalizing data for hull creation
hstlu<-hstlu[which(hstlu$tmax!="NA" & hstlu$tmin!="NA" & hstlu$cwd!="NA" & hstlu$map!="NA" & hstlu$pptconc!="NA"),]
futlu<-futlu[which(futlu$tmax!="NA" & futlu$tmin!="NA" & futlu$cwd!="NA" & futlu$map!="NA" & futlu$pptconc!="NA"),]
hstn<-hstlu
futn<-futlu

#create the climstab results table
climstab <- data.frame(planning, hullV_hst=rep(NA, length(planning)), hullV_fut=rep(NA, length(planning)), hullV_both=rep(NA, length(planning)), climstab=rep(NA, length(planning)))

#create function  to check that there is enough variation in each var to run qhull()
checker  <- function(x) {length(unique(signif(x,6)))<6}
#sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0

xx<-length(planning)

#loop through planning values
for (i in 1:xx){ # where XX is a vector of unique planning unit ids
  
  #subset data by each planning unit
  sub.hstn <- hstn[which(hstn$pu==planning[i]),]
  sub.futn <- futn[which(futn$pu==planning[i]),]
  sub.bothn <- rbind(sub.hstn, sub.futn)  
  
  
  if(sum(sapply(sub.hstn[,2:ncol(sub.hstn)], checker))>0) 
  {climstab[i,2:5] <- NA}
  else {
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
  }
} #end loop through planning values


#write results to a table
write.csv(climstab, "Output/T/rcp85/climstabvegMRIC65.csv")
rm(MRIC65)
#######################################################


