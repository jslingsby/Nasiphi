library(raster)

if(Sys.getenv("USERNAME")=="nasi") {datadir <- "C:/Users/nasip/Dropbox/Nasiphi/Fairfields_Landsat/Daily/"}
if(Sys.getenv("USER")=="jasper") {datadir <- "/Users/jasper/Dropbox/Shared/Nasiphi/Fairfields_Landsat/Daily/"}

x<-stack(paste0(datadir, "20160825_v1_LC8_L1T_TOA_TurtleConservancy_daily__2013-04-15-2016-08-13.tif"))

x <- raster(xmn=385070, xmx=389960, ymn=6190720, ymx=6197660, crs = "+proj=utm +zone=34 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs", resolution = 30, vals=0) #xmn=385970, xmx=389360, ymn=6192720, ymx=6195660
projectRaster(x, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#Get NDVI function
getNDVI=function(file,datefile,prefix)
{ndvi=stack(paste0(datadir,file))
NAvalue(ndvi)=0
offs(ndvi)=-2
gain(ndvi)=.001
dates=as.Date(read.delim(paste0(datadir,datefile),header=F, sep=":")[,2])
names(ndvi)=paste0(prefix,sub("-","",dates))
ndvi=setZ(ndvi,dates)}

#Get NDVI data - need to rename etc
l5=getNDVI(file= "20160826_v1_LT5_L1T_TOA_Fairfields_daily__1984-07-04-2011-03-09.tif", datefile="L5_dates.txt",prefix="L5_")

l7=getNDVI(file="20160826_v1_LE7_L1T_TOA_Fairfields_daily__1999-07-06-2016-08-05.tif", datefile="L7_dates.txt",prefix="L7_")

l8=getNDVI(file="20160826_v1_LC8_L1T_TOA_Fairfields_daily__2013-04-15-2016-08-13.tif", datefile="L8_dates.txt",prefix="L8_")

yr <- as.numeric(substr(names(l5),4,7))
cyr <- cut(yr, breaks=seq(1985, 2015, 5))

mxl5 <- max(l5[[which(cyr == levels(cyr)[1])]], na.rm=T)
mnl5 <- min(l5[[which(cyr == levels(cyr)[1])]], na.rm=T)
mxl7 <- max(l7, na.rm=T)
mxl8 <- max(l8, na.rm=T)

plot(mxl5, zlim=c(0,1))
plot(mxl5 > 0.65, zlim=c(0,1))

plot(mnl5, zlim=c(0,1))
plot(mnl5 > 0.1, zlim=c(0,1))

#plot(mxl5 > 0.65, col= gray(c(0,.25), alpha=.25), add=T)

plot(mxl7, zlim=c(0,1))
plot(mxl7 > 0.7, zlim=c(0,1))
plot(mxl5 > 0.8, col= gray(c(0,.25), alpha=.25), add=T)
