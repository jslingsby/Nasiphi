##########################################
######## Script for playing with Cape Lowlands 
######## data (von Hase et al.2003)
##########################################
######## Compiled by Nasiphi Ntshanga 2016
######## Last edited: 10 May 2016
##########################################

##########################################
###1) Get libraries and setwd
##########################################
#library(gdata)
#library(xlsx)
#library(XLConnect)
library(vegan)
library(simba)
library(sp)
library(raster)
library(rgdal)
library(maptools)


if(Sys.getenv("USERNAME")=="Receptionist") {giswd <- "C:/Users/Receptionist/Dropbox/Academics/PhD/Data/von Hase et al/Lowland_/";datwd <- "C:/Users/Receptionist/Dropbox/Academics/PhD/Data/von Hase et al/PRODUCTS/Fieldwork/"}

setwd("C:/Users/Receptionist/Dropbox/Academics/PhD/Data/von Hase et al/PRODUCTS/Fieldwork")

if(Sys.getenv("USER")=="jasper") {datwd <- ""; giswd <- ""}

##########################################
###2) Get and process data
##########################################
## load xls worksheets
##########################################################################
#read.xls using gdata
#boland <- read.xls("Fielddata Lowlands Project, Boland-Swartland.xls", sheet = "Habitat", perl = "C:/Strawberry/perl/bin/perl.exe")
#read.xlsx using xlsx package
#boland <- read.xlsx(paste0(datwd,"Fielddata Lowlands Project, Boland-Swartland.xls", sheetIndex = 1, sheetName =  "Habitat", rowIndex = 4, as.data.frame = T, header = T))
#read using XLconnect
#b <- readWorksheetFromFile("Fielddata Lowlands Project, Boland-Swartland.xls", "Habitat", range = 'A4:AP141')
#######3###########################################################
#Load csv data

### Boland domain 
bolandsite <- read.csv("Boland-Swartland_Habitat.csv", sep = ";", skip = 3, stringsAsFactors = T)
bolandvegdat <- read.csv("Boland-Swartland_Veg.csv", sep = ";")
bolandaliens <-read.csv("Boland-Swartland_Aliens.csv", sep = ";")
bolandspecial <-read.csv("Boland-Swartland_Special.csv", sep = ";",skip = 3)

##clean and trim
bolandsite <- bolandsite[1:137,1:43] #remove blank cells
rownames(bolandsite) <- bolandsite[,3]
write.csv(bolandsite, "Boland_sitedata.csv")
str(bolandsite)
####add spp data 

## veg dat
bolandveg <- bolandvegdat[1:13,] #create veg description data
#fix dimensions
bolandveg <- t(bolandveg) 
colnames(bolandveg) <- bolandveg[1,]
bolandveg <- as.data.frame(bolandveg[2:138,])
#bolandveg <- bb2num(bolandveg, from = c("r", "+", "1", "2", "3", "4", "5"), to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5))
write.csv(bolandveg,"Boland_veg.csv")

#aliens
bolandaliens <- bolandaliens[1:53,1:138] #remove empty rows and columns
rownames(bolandaliens) <-bolandaliens[,1] 
bolandaliens <- bolandaliens[,2:138] #remove extra row
bolandaliens <- t(bolandaliens) #transpose 
bolandaliens <- as.data.frame(bolandaliens)
#convert BB to numeric
bolandaliens <- bb2num(bolandaliens, from = c("r", "+", "1", "2", "3", "4", "5"), 
       to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5))
bolandveg$AlienSpecies_Count <- rowSums(decostand(bolandaliens, "pa")) #add alien spp count to veg dataframe
write.csv(bolandaliens,"Boland_aliens.csv")
rm(bolandaliens)

##spp dat
bolandspp <- bolandvegdat[14:314,]
bolandspp <- t(bolandspp)
colnames(bolandspp) <- bolandspp[1,]
bolandspp <- as.data.frame(bolandspp[2:138,])
bolandspp <- bb2num(bolandspp, from = c("r", "+", "1", "2", "3", "4", "5"), 
                       to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5))
bolandveg$PlantSpeciesCount <- rowSums(decostand(bolandspp, "pa"))
write.csv(bolandspp,"Boland_plant_spp.csv")
rm(bolandvegdat)
## special species??

##### Elgin
elginsite <- read.csv("Elgin_Habitat.csv", sep = ";", skip = 4)
elginveg <- read.csv("Elgin_Veg.csv", sep = ";",skip = 1)
elginaliens <-read.csv("Elgin_Aliens.csv", sep = ";", skip = 1)
elginspp <- read.csv("Elgin_Domspp.csv", sep = ";",skip = 1)
elginspecial <-read.csv("Elgin_Specials.csv", sep = ";",skip = 1,row.names = 1)
#fix stuff 
elginsite <- elginsite[1:4,1:39]
write.csv(elginsite, "Elgin_site.csv")

#aliens
rownames(elginaliens) <- elginaliens[,1]
elginaliens <- t(elginaliens[,2:5])
elginaliens <- as.data.frame(elginaliens)
elginaliens <- bb2num(elginaliens, from = c("r", "+", "1", "2", "3", "4", "5"), 
                       to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5))
elginveg$AlienSpecies_Count <- rowSums(decostand(elginaliens, "pa")) #add alien spp count to veg dataframe
write.csv(elginaliens, "Elgin_aliens.csv")

## Overberg domain
overbergsite <- read.csv("Overberg_Habitat.csv", sep = ";", skip = 3)
overbergveg <- read.csv("Overberg_Veg.csv", sep = ";", skip = 1)
overbergaliens <-read.csv("Overberg_Aliens.csv", sep = ";", skip = 1)
overbergspp <- read.csv("Overberg_Domspp.csv", sep = ";",skip = 1)
overbergspecial <-read.csv("Overberg_Special.csv", sep = ";",skip = 1)

#fix stuff
overbergsite <- overbergsite[1:155,]
write.csv(overbergsite, "Overberg_site.csv")
#veg
overbergveg <- t(overbergveg)
colnames(overbergveg) <- overbergveg[1,]
overbergveg <- as.data.frame(overbergveg[2:156,])
write.csv(overbergveg, "Overberg_veg.csv")
#aliens
overbergaliens <- overbergaliens[1:35,1:156]
rownames(overbergaliens) <- overbergaliens[,1]
overbergaliens <- t(overbergaliens[,2:156])
overbergaliens <- as.data.frame(overbergaliens)
overbergaliens <- bb2num(overbergaliens, from = c("r", "+", "1", "2", "3", "4", "5"), 
                      to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5))
overbergveg$AlienSpecies_Count <- rowSums(decostand(overbergaliens, "pa"))
write.csv(overbergaliens, "Overberg_aliens.csv")
#spp
overbergspp <- t(overbergspp[,2:157])
colnames(overbergspp) <- overbergspp[1,]
overbergspp <- as.data.frame(overbergspp[2:156,])
overbergspp <- bb2num(overbergspp, from = c("r", "+", "1", "2", "3", "4", "5"), 
                      to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5))
overbergveg$PlantSpecies_Count <- rowSums(decostand(overbergspp, "pa"))


#########################################################################
###add gis stuff

stdCRS <- "+proj=utm +zone=34 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

fragments <- readOGR(paste0(giswd, "Remnants/Remnant_Shape"), layer = "remnants_only_wtm21")

proj4string(fragments) <- CRS("+proj=utm +zone=34 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

writeOGR(fragments, dsn="lowland_fragments.kml", layer= "remnants_wtm21", driver="KML")

