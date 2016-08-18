############################################
###### Code to collate  climate stabilty ##
###### fot TRANSFORMED veg
###### created by Nasiphi Ntshanga ########
###### last edit 9 June 2016 ##############
###########################################

######################################
#### Step 1) Load data and libraries
#### Step 2) Collate results
#### Step 3) Explore results
######################################

######### 1) load data 

#climate models info 
info <- read.csv("C:/Users/nasip/Dropbox/Academics/PhD/Data/ClimateData/futureclimate_info.csv", sep = ",")

#list.files(path = "Data",pattern = "climstabveg")
#RCP45 scenario for untransformed veg map
rcp45bcc65 <- read.csv("Output/T/rcp45/climstabvegbcc65.csv", sep=",")
rcp45bcc100 <- read.csv("Output/T/rcp45/climstabvegbcc100.csv", sep=",")
rcp45BNUESM100 <- read.csv("Output/T/rcp45/climstabvegBNUESM100.csv", sep=",")
rcp45BNUESM65 <- read.csv("Output/T/rcp45/climstabvegBNUESM65.csv", sep=",")
rcp45CanESM100 <- read.csv("Output/T/rcp45/climstabvegCanESM100.csv", sep=",")
rcp45CanESM65 <- read.csv("Output/T/rcp45/climstabvegCanESM65.csv", sep=",")
rcp45CNRMCM100 <- read.csv("Output/T/rcp45/climstabvegCNRMCM100.csv", sep=",")
rcp45CNRMCM65 <- read.csv("Output/T/rcp45/climstabvegCNRMCM65.csv", sep=",")
rcp45FGOALSs100 <- read.csv("Output/T/rcp45/climstabvegFGOALSs100.csv", sep=",")
rcp45FGOALSs65 <- read.csv("Output/T/rcp45/climstabvegFGOALSs65.csv", sep=",")
rcp45GFDL65 <- read.csv("Output/T/rcp45/climstabvegGFDL65.csv", sep=",")
rcp45GFDL100 <- read.csv("Output/T/rcp45/climstabvegGFDL100.csv", sep=",")
rcp45GFDLESML100 <- read.csv("Output/T/rcp45/climstabvegGFDLESM100.csv", sep=",")
rcp45GFDLESML65 <- read.csv("Output/T/rcp45/climstabvegGFDLESM65.csv", sep=",")
rcp45MIROC65 <- read.csv("Output/T/rcp45/climstabvegMIROC65.csv", sep=",")
rcp45MIROC100 <-read.csv("Output/T/rcp45/climstabvegMIROC100.csv", sep=",")
rcp45MIROCES100 <- read.csv("Output/T/rcp45/climstabvegMIROCES100.csv", sep=",")
rcp45MIROCES65 <- read.csv("Output/T/rcp45/climstabvegMIROCES65.csv", sep=",")
rcp45MRIC100 <- read.csv("Output/T/rcp45/climstabvegMRIC100.csv", sep=",")
rcp45MRIC65 <- read.csv("Output/T/rcp45/climstabvegMRIC65.csv", sep=",")
rcp45MIROCHEM65 <-read.csv("Output/T/rcp45/climstabvegMIROCHEM65.csv", sep = ",")
rcp45MIROCHEM100 <-read.csv("Output/T/rcp45/climstabvegMIROCHEM100.csv",sep = ",")

#do the same with rcp85
#list.files(path = "Output/T/rcp85",pattern = "climstabveg")
rcp85bcc65 <- read.csv("Output/T/rcp85/climstabvegbcc65.csv", sep=",")
rcp85bcc100 <- read.csv("Output/T/rcp85/climstabvegbcc100.csv", sep=",")
rcp85BNUESM100 <- read.csv("Output/T/rcp85/climstabvegBNUESM100.csv", sep=",")
rcp85BNUESM65 <- read.csv("Output/T/rcp85/climstabvegBNUESM65.csv", sep=",")
rcp85CanESM100 <- read.csv("Output/T/rcp85/climstabvegCanESM100.csv", sep=",")
rcp85CanESM65 <- read.csv("Output/T/rcp85/climstabvegCanESM65.csv", sep=",")
rcp85CNRMCM100 <- read.csv("Output/T/rcp85/climstabvegCNRMCM100.csv", sep=",")
rcp85CNRMCM65 <- read.csv("Output/T/rcp85/climstabvegCNRMCM65.csv", sep=",")
rcp85FGOALSs100 <- read.csv("Output/T/rcp85/climstabvegFGOALSs100.csv", sep=",")
rcp85FGOALSs65 <- read.csv("Output/T/rcp85/climstabvegFGOALSs65.csv", sep=",")
rcp85GFDL65 <- read.csv("Output/T/rcp85/climstabvegGFDL65.csv", sep=",")
rcp85GFDL100 <- read.csv("Output/T/rcp85/climstabvegGFDL100.csv", sep=",")
rcp85GFDLESML100 <- read.csv("Output/T/rcp85/climstabvegGFDLESM100.csv", sep=",")
rcp85GFDLESML65 <- read.csv("Output/T/rcp85/climstabvegGFDLESM65.csv", sep=",")
rcp85MIROC65 <- read.csv("Output/T/rcp85/climstabvegMIROC65.csv", sep=",")
rcp85MIROC100 <-read.csv("Output/T/rcp85/climstabvegMIROC100.csv", sep=",")
rcp85MIROCES100 <- read.csv("Output/T/rcp85/climstabvegMIROCES100.csv", sep=",")
rcp85MIROCES65 <- read.csv("Output/T/rcp85/climstabvegMIROCES65.csv", sep=",")
rcp85MRIC100 <- read.csv("Output/T/rcp85/climstabvegMRIC100.csv", sep=",")
rcp85MRIC65 <- read.csv("Output/T/rcp85/climstabvegMRIC65.csv", sep=",")
rcp85MIROCHEM65 <-read.csv("Output/T/rcp85/climstabvegMIROCHEM65.csv", sep = ",")
rcp85MIROCHEM100 <-read.csv("Output/T/rcp85/climstabvegMIROCHEM100.csv",sep = ",")

veg <- read.csv("Output/vegtypes.csv")
westcoast <- veg[which(veg$ID %in% rcp45bcc65$planning),]


##### 2) Collate outputs 

### create a dataframe for time frame for 20462065 historical
HullV_hst2065 <- cbind(rcp45MRIC65$planning,rcp45bcc65$hullV_hst,rcp45BNUESM65$hullV_hst,rcp45CanESM65$hullV_hst, rcp45CNRMCM65$hullV_hst,rcp45FGOALSs65$hullV_hst,rcp45GFDL65$hullV_hst,rcp45GFDLESML65$hullV_hst,rcp45MIROC65$hullV_hst,rcp45MIROCES65$hullV_hst,rcp45MIROCHEM65$hullV_hst,rcp45MRIC65$hullV_hst,rcp85bcc65$hullV_hst,rcp85BNUESM65$hullV_hst,rcp85CanESM65$hullV_hst,rcp85CNRMCM65$hullV_hst,rcp85FGOALSs65$hullV_hst,rcp85GFDL65$hullV_hst,rcp85GFDLESML65$hullV_hst,rcp85MIROC65$hullV_hst,rcp85MIROCES65$hullV_hst,rcp85MIROCHEM65$hullV_hst,rcp85MRIC65$hullV_hst)

HullV_hst2065 <- data.frame(HullV_hst2065)
HullV_hst2065<-merge.data.frame(HullV_hst2065,westcoast, by.x = "X1", by.y = "ID")
rownames(HullV_hst2065)<-HullV_hst2065[,28]
colnames(HullV_hst2065)<-(c("vegtype","rcp45bcc65","rcp45BNUESM65", "rcp45CanESM65","rcp45CNRMCM65","rcp45FGOALSs65","rcp45GFDL65", "rcp45GFDLESML65", "rcp45MIROC65","rcp45MIROCES65","rcp45MIROCHEM65", "rcp45MRIC65","rcp85bcc65","rcp85BNUESM65","rcp85CanESM65","rcp85CNRMCM65","rcp85FGOALSs65","rcp85GFDL65",  "rcp85GFDLESML65", "rcp85MIROC65", "rcp85MIROCES65","rcp85MIROCHEM65","rcp85MRIC65"))
HullV_hst2065 <- data.frame(t(HullV_hst2065[,2:23]))

HullV_hst2065 <- cbind(HullV_hst2065,scenario=c("RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45", "RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85"))

### create a dataframe for future

HullV_fut2065 <- cbind(rcp45MRIC65$planning,rcp45bcc65$hullV_fut,rcp45BNUESM65$hullV_fut,rcp45CanESM65$hullV_fut, rcp45CNRMCM65$hullV_fut,rcp45FGOALSs65$hullV_fut,rcp45GFDL65$hullV_fut,rcp45GFDLESML65$hullV_fut,rcp45MIROC65$hullV_fut,rcp45MIROCES65$hullV_fut,rcp45MIROCHEM65$hullV_fut,rcp45MRIC65$hullV_fut,rcp85bcc65$hullV_fut,rcp85BNUESM65$hullV_fut,rcp85CanESM65$hullV_fut,rcp85CNRMCM65$hullV_fut,rcp85FGOALSs65$hullV_fut,rcp85GFDL65$hullV_fut,rcp85GFDLESML65$hullV_fut,rcp85MIROC65$hullV_fut,rcp85MIROCES65$hullV_fut,rcp85MIROCHEM65$hullV_fut,rcp85MRIC65$hullV_fut)

HullV_fut2065 <- data.frame(HullV_fut2065)
HullV_fut2065<-merge.data.frame(HullV_fut2065,westcoast, by.x = "X1", by.y = "ID")
rownames(HullV_fut2065)<-HullV_fut2065[,28]
colnames(HullV_fut2065)<-(c("vegtype","rcp45bcc65","rcp45BNUESM65", "rcp45CanESM65","rcp45CNRMCM65","rcp45FGOALSs65","rcp45GFDL65", "rcp45GFDLESML65", "rcp45MIROC65","rcp45MIROCES65","rcp45MIROCHEM65", "rcp45MRIC65","rcp85bcc65","rcp85BNUESM65","rcp85CanESM65","rcp85CNRMCM65","rcp85FGOALSs65","rcp85GFDL65",  "rcp85GFDLESML65", "rcp85MIROC65", "rcp85MIROCES65","rcp85MIROCHEM65","rcp85MRIC65"))
HullV_fut2065 <- data.frame(t(HullV_fut2065[,2:23]))

HullV_fut2065 <- cbind(HullV_fut2065,scenario=c("RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45", "RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85"))

#### create a dataframe for climstab
climstab2065 <- cbind(rcp45MRIC65$planning,rcp45bcc65$climstab,rcp45BNUESM65$climstab,rcp45CanESM65$climstab, rcp45CNRMCM65$climstab,rcp45FGOALSs65$climstab,rcp45GFDL65$climstab,rcp45GFDLESML65$climstab,rcp45MIROC65$climstab,rcp45MIROCES65$climstab,rcp45MIROCHEM65$climstab,rcp45MRIC65$climstab,rcp85bcc65$climstab,rcp85BNUESM65$climstab,rcp85CanESM65$climstab,rcp85CNRMCM65$climstab,rcp85FGOALSs65$climstab,rcp85GFDL65$climstab,rcp85GFDLESML65$climstab,rcp85MIROC65$climstab,rcp85MIROCES65$climstab,rcp85MIROCHEM65$climstab,rcp85MRIC65$climstab)

climstab2065 <- data.frame(climstab2065)
climstab2065<-merge.data.frame(climstab2065,westcoast, by.x = "X1", by.y = "ID")
rownames(climstab2065)<-climstab2065[,28]
colnames(climstab2065)<-(c("vegtype","rcp45bcc65","rcp45BNUESM65", "rcp45CanESM65","rcp45CNRMCM65","rcp45FGOALSs65","rcp45GFDL65", "rcp45GFDLESML65", "rcp45MIROC65","rcp45MIROCES65","rcp45MIROCHEM65", "rcp45MRIC65","rcp85bcc65","rcp85BNUESM65","rcp85CanESM65","rcp85CNRMCM65","rcp85FGOALSs65","rcp85GFDL65",  "rcp85GFDLESML65", "rcp85MIROC65", "rcp85MIROCES65","rcp85MIROCHEM65","rcp85MRIC65"))
climstab2065 <- data.frame(t(climstab2065[,2:23]))

climstab2065 <- cbind(climstab2065,scenario=c("RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45", "RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85"))

###############################################################################

#reapeat for 20812100 climatic period

#historic hull volume
HullV_hst2100 <- cbind(rcp45MRIC100$planning,rcp45bcc100$hullV_hst,rcp45BNUESM100$hullV_hst,rcp45CanESM100$hullV_hst, rcp45CNRMCM100$hullV_hst,rcp45FGOALSs100$hullV_hst,rcp45GFDL100$hullV_hst,rcp45GFDLESML100$hullV_hst,rcp45MIROC100$hullV_hst,rcp45MIROCES100$hullV_hst,rcp45MIROCHEM100$hullV_hst,rcp45MRIC100$hullV_hst,rcp85bcc100$hullV_hst,rcp85BNUESM100$hullV_hst,rcp85CanESM100$hullV_hst,rcp85CNRMCM100$hullV_hst,rcp85FGOALSs100$hullV_hst,rcp85GFDL100$hullV_hst,rcp85GFDLESML100$hullV_hst,rcp85MIROC100$hullV_hst,rcp85MIROCES100$hullV_hst,rcp85MIROCHEM100$hullV_hst,rcp85MRIC100$hullV_hst)

HullV_hst2100 <- data.frame(HullV_hst2100)
rownames(HullV_hst2100)<-HullV_hst2100[,1]
colnames(HullV_hst2100)<-(c("vegtype","rcp45bcc100","rcp45BNUESM100", "rcp45CanESM100","rcp45CNRMCM100","rcp45FGOALSs100","rcp45GFDL100", "rcp45GFDLESML100", "rcp45MIROC100","rcp45MIROCES100","rcp45MIROCHEM100", "rcp45MRIC100","rcp85bcc100","rcp85BNUESM100","rcp85CanESM100","rcp85CNRMCM100","rcp85FGOALSs100","rcp85GFDL100",  "rcp85GFDLESML100", "rcp85MIROC100", "rcp85MIROCES100","rcp85MIROCHEM100","rcp85MRIC100"))
HullV_hst2100 <- data.frame(t(HullV_hst2100[,2:23]))

HullV_hst2100 <- cbind(HullV_hst2100,scenario=c("RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45", "RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85"))

### create a dataframe for future

HullV_fut2100 <- cbind(rcp45MRIC100$planning,rcp45bcc100$hullV_fut,rcp45BNUESM100$hullV_fut,rcp45CanESM100$hullV_fut, rcp45CNRMCM100$hullV_fut,rcp45FGOALSs100$hullV_fut,rcp45GFDL100$hullV_fut,rcp45GFDLESML100$hullV_fut,rcp45MIROC100$hullV_fut,rcp45MIROCES100$hullV_fut,rcp45MIROCHEM100$hullV_fut,rcp45MRIC100$hullV_fut,rcp85bcc100$hullV_fut,rcp85BNUESM100$hullV_fut,rcp85CanESM100$hullV_fut,rcp85CNRMCM100$hullV_fut,rcp85FGOALSs100$hullV_fut,rcp85GFDL100$hullV_fut,rcp85GFDLESML100$hullV_fut,rcp85MIROC100$hullV_fut,rcp85MIROCES100$hullV_fut,rcp85MIROCHEM100$hullV_fut,rcp85MRIC100$hullV_fut)

HullV_fut2100 <- data.frame(HullV_fut2100)
rownames(HullV_fut2100)<-HullV_fut2100[,1]
colnames(HullV_fut2100)<-(c("vegtype","rcp45bcc100","rcp45BNUESM100", "rcp45CanESM100","rcp45CNRMCM100","rcp45FGOALSs100","rcp45GFDL100", "rcp45GFDLESML100", "rcp45MIROC100","rcp45MIROCES100","rcp45MIROCHEM100", "rcp45MRIC100","rcp85bcc100","rcp85BNUESM100","rcp85CanESM100","rcp85CNRMCM100","rcp85FGOALSs100","rcp85GFDL100",  "rcp85GFDLESML100", "rcp85MIROC100", "rcp85MIROCES100","rcp85MIROCHEM100","rcp85MRIC100"))
HullV_fut2100 <- data.frame(t(HullV_fut2100[,2:23]))

HullV_fut2100 <- cbind(HullV_fut2100,scenario=c("RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45", "RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85"))

#### create a dataframe for climstab
climstab2100 <- cbind(rcp45MRIC100$planning,rcp45bcc100$climstab,rcp45BNUESM100$climstab,rcp45CanESM100$climstab, rcp45CNRMCM100$climstab,rcp45FGOALSs100$climstab,rcp45GFDL100$climstab,rcp45GFDLESML100$climstab,rcp45MIROC100$climstab,rcp45MIROCES100$climstab,rcp45MIROCHEM100$climstab,rcp45MRIC100$climstab,rcp85bcc100$climstab,rcp85BNUESM100$climstab,rcp85CanESM100$climstab,rcp85CNRMCM100$climstab,rcp85FGOALSs100$climstab,rcp85GFDL100$climstab,rcp85GFDLESML100$climstab,rcp85MIROC100$climstab,rcp85MIROCES100$climstab,rcp85MIROCHEM100$climstab,rcp85MRIC100$climstab)

climstab2100 <- data.frame(climstab2100)
rownames(climstab2100)<-climstab2100[,1]
colnames(climstab2100)<-(c("vegtype","rcp45bcc100","rcp45BNUESM100", "rcp45CanESM100","rcp45CNRMCM100","rcp45FGOALSs100","rcp45GFDL100", "rcp45GFDLESML100", "rcp45MIROC100","rcp45MIROCES100","rcp45MIROCHEM100", "rcp45MRIC100","rcp85bcc100","rcp85BNUESM100","rcp85CanESM100","rcp85CNRMCM100","rcp85FGOALSs100","rcp85GFDL100",  "rcp85GFDLESML100", "rcp85MIROC100", "rcp85MIROCES100","rcp85MIROCHEM100","rcp85MRIC100"))
climstab2100 <- data.frame(t(climstab2100[,2:23]))

climstab2100 <- cbind(climstab2100,scenario=c("RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45","RCP45", "RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85","RCP85"))

###########

### write out outputs

write.csv(climstab2100, file ="Output/Climstab2100_T.csv")
write.csv(climstab2065, file = "Output/Climstab2065_T.csv")
write.csv(HullV_fut2100, file = "Output/HullV_fut2100_T.csv")
write.csv(HullV_fut2065, file = "Output/HullV_fut2065_T.csv")
write.csv(HullV_hst2065, file = "Output/HullV_hst2065_T.csv")
write.csv(HullV_hst2100, file = "Output/HullV_hst2100_T.csv")

###########################################################################


