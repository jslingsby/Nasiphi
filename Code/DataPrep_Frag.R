##########################################
######## R code for preparing  data for 
######## transformation & fragmentation 
######## analysis
##########################################
######## Compiled by Nasiphi Ntshanga 2016
######## Last edited: 09 May 2016
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
subs(lc,tab, by="ID",which="LC3V",subsWithNA=TRUE, filename='lcs', overwrite=T) #reclassifies raster to 1s and 0s
lc <- raster("Data/lcs.grd")


## mask veg map based on untransformed landcover to get fragments layer
fragments <- mask(vm, lc, filename="Data/vmt", maskvalue= 0, overwrite=TRUE)

##make map for each veg type and get patch stats within eachveg type
vegtab <- fragments@data@attributes[[1]]

#extract by vegtype names
#use which function  

#########################################
##### 2) Fragmentation analysis #########
##########################################

#calculate fragstats for remnant natural land cover

#fragmentation of all veg types
classfrag = ClassStat(fragments)

#patch stats for each veg type
#[need to subset each veg type]


#landcover
conlc = ConnCompLabel(lc) #identifies disjunct patches
lcpatch = PatchStat(conlc,cellsize = 30) #calculate patch statistics
colnames(lcpatch)[1] <- "ID"
levels(conlc) <- lcpatch
lcclass = ClassStat(conlc) #calculate class statistics


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





