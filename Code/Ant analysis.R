##########################################
######## R code for OTS Elandsberg ant FFP
##########################################
######## Compiled by Jasper Slingsby 2016
######## edited by Nasiphi Ntshanga 
######## Last edited:4 April 2016
##########################################

##########################################
###1) Get libraries, setwd and get data
##########################################

library(vegan)
library(betapart)
library(ade4)

setwd("~/Nasiphi's/Data")

###########################################################################
#is ant diversity dependant on land cover 

antdat <- read.csv("antsOTS.csv",sep=";", row.names=1)
sitedatAll <- read.csv("ants.csv",sep=";", row.names=1) #coarse site dat inc fallow
#remove plots with zero count
AntIndividualsCount <- rowSums(antdat) #create object calculating row sums
antdat <- subset(antdat, AntIndividualsCount > 0)


#match antdat with site dat
sitedatAll <- sitedatAll[which(rownames(sitedatAll) %in% rownames(antdat)),]

#calculate abundance, estimate species richness
sitedatAll$AntIndividualsCount <- rowSums(antdat)
sitedatAll$AntSpeciesCount <- rowSums(decostand(antdat, "pa")) 
sitedatAll$AntSpeciesRichness <- rarefy(antdat, sample = min(sitedatAll$AntIndividualsCount))

hist(sitedatAll$AntIndividualsCount)
hist(sitedatAll$AntSpeciesCount)
hist(sitedatAll$AntSpeciesRichness)


boxplot(AntSpeciesCount ~ VegClass, data=sitedatAll)
boxplot(AntIndividualsCount ~ VegClass, data=sitedatAll)
boxplot(AntSpeciesRichness ~ VegClass, data=sitedatAll)
###########################################################################
pdf("SppRichness.pdf", width=4, height=4) #Turn on the pdf plotting device, set file name
boxplot(AntSpeciesRichness ~ VegClass, data=sitedatAll)
dev.off() #Turn off the pdf plotting device

pdf("SppCount.pdf", width=4, height=4) #Turn on the pdf plotting device, set file name
boxplot(AntSpeciesCount ~ VegClass, data=sitedatAll)
dev.off() #Turn off the pdf plotting device

pdf("IndividualCount.pdf", width=4, height=4) #Turn on the pdf plotting device, set file name
boxplot(AntIndividualsCount ~ VegClass, data=sitedatAll)
dev.off() #Turn off the pdf plotting device

###########################################################################
#Generalized Linear Model
#load data
sitedat <- read.csv("sitedataOTS.csv")
sitedatsub <-sitedat[which(sitedat$Plot %in% rownames(antdat)),] #subset sitedat based on ant plots (30)

antdatsub <- antdat[-which(sitedatAll$VegClass=="AF"), ]# (31)
table(sitedatAll$VegClass) 

#antdatsub <- antdatsub[rownames(antdatsub)%in%sitedatsub$Plot,] #match site dat with plots

sitedatsub$AntIndividualsCount <- rowSums(antdatsub)
sitedatsub$AntSpeciesCount <- rowSums(decostand(antdatsub, "pa")) 
sitedatsub$AntSpeciesRichness <- rarefy(antdatsub, sample = min(sitedatsub$AntIndividualsCount))
###########################################################################
#GLM for count data
#check for covariates for ant diversity
#[plant diversity],[alien cover], pig disturbance, trampling, veg class,elevation, ndvi, termite mound cc,dung, trampling s,elevation
colnames(sitedatsub)

pairs(scale(sitedatsub[,c("Species_richness", "Species_richness_perennial","Live_veg_cover..S.", "Dung..S.", "Alien_grass_cover", "Alien_all_cover", "Pig_disturbance..C", "Termite_mound..C.", "Dung..S.", "Trampling..S.", "Litter_Cover", "ndvi","Elevation", "Sand")]))

pairs(scale(sitedat[,c("Sand", "Soil_Depth","Live_veg_cover..S.", "Dung..S.", "Bedrock", "Bare_soil", "Gravel.", "Cobble", "Clay", "Loam", "Litter_Cover")]))

fit <- glm(AntSpeciesRichness ~ Dung..S., data = sitedatsub, family = "gaussian")

fit1<-glm(AntSpeciesCount ~ Reserve + VegClass, family = "poisson", data=sitedatsub)
summary(fit1)

fit2<-glm(AntSpeciesCount ~ Reserve, family = "poisson", data=antsitedat)
summary(fit2)


##########################################################################
# calculate ant composition among plots
vegdat <- read.csv("vegdataOTS.csv", row.names=1)
#subset veg dat

vegdatsub <- vegdat[which(rownames(vegdat)%in%rownames(antdatsub)),]

disV = beta.pair (decostand(vegdatsub,"pa"), index.family="sorensen")$beta.sim 
disA = beta.pair (decostand(antdatsub,"pa"), index.family="sorensen")$beta.sim 
disE <- vegdist(scale(sitedatsub[,c("Bare_soil", "Gravel.", "Cobble", "Sand", "Clay", "Loam")]), "euclid")

treeV = hclust(disV, method="average") 
treeA = hclust(disA, method="average")
treeE = hclust(disE, method = "average")

plot(treeV, labels=rownames(vegdatsub), main="", xlab="", sub="", ylab="Sorenson's similarity") 
#looks like we have two main groups of data
# Cut the tree into 2 classes based on similarity among plots
plot(treeA, labels=rownames(antdatsub), main="", xlab="", sub="", ylab="Sorenson's similarity")
plot(treeE, labels=rownames(sitedatsub), main="", xlab="", sub="", ylab="Sorenson's similarity")

which(rownames(sitedatAll)%in%rownames(antdatsub))
sitedatAllsub <- sitedatAll[c(1:16,24:29,31:38,40),]
table(AntTest, sitedatAllsub$VegClass)
#cuts tree (from hclust) into groups of data
#Let's see how that compares to our actual data?

table(AntTest, sitedatAll$VegClass) 

#mantel test
mantel(disE, disV, method="pearson", permutations=999)

###MRPP test of differences in composition between veg types and fallow fields

disAnts = beta.pair (decostand(antdat,"pa"), index.family="sorensen")$beta.sim 
mrpp(disAnts, grouping = sitedatAll$VegClass, permutations = 999)
     
###Indicator power of species

ip <- indpower(antdat)

diag(ip) <- NA
(TIP <- rowMeans(ip, na.rm=TRUE))

## p value calculation for a species
## from Halme et al. 2009
## i is ID for the species
i <- 1
fun <- function(x, i) indpower(x)[i,-i]
## 'c0' randomizes species occurrences
os <- oecosimu(antdat, fun, "c0", i=i, nsimul=999)
## get z values from oecosimu output
z <- os$oecosimu$z
## p-value
(p <- sum(z) / sqrt(length(z)))
