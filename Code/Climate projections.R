##############################################
######### Code to assess climate change model#
######### predictions in the CFR #############
########## Written by Nasiphi ##############
######### last edit 20 June #################
################################################



##########################################
###1) Get libraries and setwd
##########################################
library(rgdal)
library(raster)
library(animation)

if(Sys.getenv("USERNAME")=="nasip") {climwd <- "C:/Users/nasip/Dropbox/Academics/PhD/Data/ClimateData/";datwd <- "C:/Users/nasip/Dropbox/Academics/PhD/Data/Rasters/"}

#########################################
###3) Get and process climate data 
#########################################

#load rasters I created with complete layer names
future_pptconc <-brick("Data/futurepptconc.grd")
future_mmp <- brick("Data/futuremmp.grd")
future_map <-brick("Data/futuremap.grd")
future_tmin <-brick("Data/futuretmin.grd")
future_tmax <-brick("Data/futuretmax.grd")
#################
#load rasters developed by Jasper)....layer names mising climate variable
#futureplus <-brick(paste(climwd,"futureage_plus_Current.grd", sep = ""))
#futuremap <-brick(paste(climwd,"map/futureclimate_plus_Current.grd", sep = ""))
#futuremmp <-brick(paste(climwd,"mmp01/futureclimate_plus_Current.grd", sep = ""))
#futurepptconc <-brick(paste(climwd,"pptconc/futureclimate_plus_Current.grd", sep = ""))
#futuretmax <-brick(paste(climwd,"tmax01/futureclimate_plus_Current.grd", sep = ""))
#futuretmin <-brick(paste(climwd,"tmin07/futureclimate_plus_Current.grd", sep = ""))
#################
#### Plot all models for each variable, plot >0 only
#mmp
oopts = if (.Platform$OS.type == "windows") {
  ani.options(ffmpeg = "C:/Program Files/ImageMagick-7.0.2-Q16/ffmpeg.exe")
}
## usually Linux users do not need to worry about the ffmpeg path as long as
## FFmpeg or avconv has been installed

saveVideo({
  x <-nlayers(future_mmp)
  for (i in 1:x) plot(future_mmp[[i]]>0,
                      main=names(future_mmp[[i]]),
                      legend.lab="Mean precipitation in January")},
  video.name = "Output/mmp.mp4", other.opts = "-pix_fmt yuv420p -b 300k")
# higher bitrate, better quality
ani.options(oopts)
#######################
#map
oopts = if (.Platform$OS.type == "windows") {
  ani.options(ffmpeg = "C:/Program Files/ImageMagick-7.0.2-Q16/ffmpeg.exe")
}

saveVideo({
  x <-nlayers(future_map)
  for (i in 1:x) plot(future_map[[i]]>0,
                      main=names(future_map[[i]]),
                      legend.lab="Mean annual precipitation")},
  video.name = "Output/map.mp4", other.opts = "-pix_fmt yuv420p -b 300k")
# higher bitrate, better quality
ani.options(oopts)
###################
#tmax
oopts = if (.Platform$OS.type == "windows") {
  ani.options(ffmpeg = "C:/Program Files/ImageMagick-7.0.2-Q16/ffmpeg.exe")
}

saveVideo({
  x <-nlayers(future_tmax)
  for (i in 1:x) plot(future_tmax[[i]]>0,
                      main=names(future_tmax[[i]]),
                      legend.lab="Maximum temperature in summer")},
  video.name = "Output/tmax.mp4", other.opts = "-pix_fmt yuv420p -b 300k")
# higher bitrate, better quality
ani.options(oopts)
##############################
#pptconc
oopts = if (.Platform$OS.type == "windows") {
  ani.options(ffmpeg = "C:/Program Files/ImageMagick-7.0.2-Q16/ffmpeg.exe")
}

saveVideo({
  x <-nlayers(future_pptconc)
  for (i in 1:x) plot(future_pptconc[[i]]>0,
                      main=names(future_pptconc[[i]]),
                      legend.lab="Precipitation concerntration")},
  video.name = "Output/pptconc.mp4", other.opts = "-pix_fmt yuv420p -b 300k")
# higher bitrate, better quality
ani.options(oopts)
########################
#tmin
oopts = if (.Platform$OS.type == "windows") {
  ani.options(ffmpeg = "C:/Program Files/ImageMagick-7.0.2-Q16/ffmpeg.exe")
}

saveVideo({
  x <-nlayers(future_tmin)
  for (i in 1:x) plot(future_tmin[[i]]>0,
                      main=names(future_tmin[[i]]),
                      legend.lab="Minimum temperature in Winter")},
  video.name = "Output/tmin.mp4", other.opts = "-pix_fmt yuv420p -b 300k")
# higher bitrate, better quality
ani.options(oopts)

