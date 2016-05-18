##########################################
######## R code for preparing  data for 
######## transformation & fragmentation 
######## analysis
##########################################
######## Compiled by Nasiphi Ntshanga 2016
######## Last edited: 16 May 2016
##########################################

##########################################
###1) Get libraries and setwd
##########################################
library(raster)
library(gdalUtils)
library(rgdal)
library(dplyr)
library(SDMTools)
library(igraph)

if(Sys.getenv("USERNAME")=="Receptionist") {datwd <- "C:/Users/Receptionist/Dropbox/Academics/PhD/Data/"}
if(Sys.getenv("USER")=="jasper") {datwd <- "/Users/jasper/Documents/GIS/VegToolsRaw/"}
if(Sys.getenv("USERNAME")=="jasper") {datwd <- "~/Nasiphi's/Data/"}
##########################################
###2) Get and process data
##########################################

#set extent
nex <- extent(2030000, 2060000, -4080000, -4010000) #extent(2020000, 2140000, -4100000, -3950000)

#Land cover
lc <- raster(paste(datwd,"Rasters/LC13test_web.tif", sep=""))
proj4string(lc) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs") #EPSG3857 - Web Mercator

#vegetation
vm <- raster(paste(datwd, "Rasters/VEG12test_web.tif", sep=""))
proj4string(vm) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

#Protected areas
pa <- raster(paste(datwd,"/Rasters/PA_npaes_test_web.tif", sep=""))
proj4string(pa) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

#crop data
lc <- crop(lc, nex)
vm <- crop(vm, nex)
pa <- crop(pa, nex)

#insert if/else arguement for lc


## Reclassify landcover data
NAvalue(lc) <- 128
NAvalue(lc) <- 0
tab <-lc@data@attributes[[1]]
tb <- tab[,c("Value","LC3V")]
tb <- rbind(tb, c(128,"NA"))
tb <- rbind(tb, c(0,"NA"))

if(!file.exists("Data/lcs.grd")) {
  subs(lc,tab, by="ID",which="LC3V",subsWithNA=TRUE, filename= "Data/lcs", overwrite=T)
} else {lc <- raster("Data/lcs.grd")}

#create a transformed veg map
if(!file.exists("Data/vmt.grd")) {
  veg <- mask(vm, lc, filename="Data/vmt", maskvalue= 0, overwrite=TRUE)
} else {veg <- raster("Data/vmt.grd")}

### create a raster layer for each veg type
# find a loop function that creates this

##HAVE TO AGGREGATE EACH RASTER
vtab <- vm@data@attributes[[1]]

veg2 <- veg==2 #veg type 2 only
aggregate(veg2, fact=3, fun=mean, expand=TRUE, na.rm=TRUE, filename='veg2ag')
veg2ag <- raster("veg2ag.grd")
veg59 <- veg ==59


#########################################
##### 2) Fragmentation analysis #########
##########################################

#Class statistics for different veg types
classfrag <- ClassStat(veg, cellsize = 30)

#patch stats for each veg type #[need to subset each veg type]

#landcover
conlc = ConnCompLabel(lc) #identifies disjunct patches
lcpatch = PatchStat(conlc,cellsize = 30) #calculate patch statistics
colnames(lcpatch)[1] <- "ID"
levels(conlc) <- lcpatch

#vegtype 2
conv2 = ConnCompLabel(veg2)
v2patch = PatchStat(conv2,cellsize = 30) #calculate patch statistics
colnames(v2patch)[1] <- "ID"
levels(conv2) <- v2patch


#vegtype 2 aggregated to 90m
conv2a = ConnCompLabel(veg2ag)
v2apatch = PatchStat(conv2a,cellsize = 90) #calculate patch statistics
colnames(v2apatch)[1] <- "ID"
levels(conv2a) <- v2apatch

#vegtype 59
conv59 = ConnCompLabel(veg59)
v59patch = PatchStat(conv59,cellsize = 30) #calculate patch statistics
colnames(v59patch)[1] <- "ID"
levels(conv59) <- v59patch

#########################################
##### 3) Transformation analysis #########
##########################################

###calculate rate of transformation (1990 minus 2014)
#calculate min,mean and max natural land cover by fragment and by veg type
#determine min veg type size appropriate for fragstats

###make table

#getvalues
VegType <- getValues(vm)
Fragments <- getValues(fragments)
Landcover <- getValues(lc)
Protected <- getValues(pa)

dat <- data.frame(Landcover, Fragments , Protected, VegType)
rm(VegType)
rm(Fragments)
rm(Landcover)
rm(Protected)

dat <- na.omit(dat)
dtbl <- tbl_df(dat)
rm(dat)

#original <- summarise(group_by(dtbl, VegType), original = (length(Landcover)*900)/10000)
#remaining <- summarise(group_by(dtbl, VegType), remaining = (sum(Landcover)*900)/10000)
#protected <- summarise(group_by(filter(dtbl, Protected>0) , Fragments)) protected = (sum(Landcover)*900)/10000))





