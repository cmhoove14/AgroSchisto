#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Investigating results of Halstead et al 2016 mesocosm experiment
#Data file "R_use.csv" contains results of experiment investigated here
#email choover@berkeley.edu to obtain

require(ggplot2)
require(reshape2)

st.er <- function(x) {
  sd(x)/sqrt(length(x))
} #Function to calculate standard error of the mean

datr<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/R_use.csv")

#Add some variables ##########################
for(i in 1:nrow(datr)){
  datr[i,43]=paste(datr[i,2], datr[i,3], datr[i,4], sep="_")
} #Add unique treatment code to data frame
treats<-unique(datr[,43])
colnames(datr)[43]<-"Treatments"

for(i in 1:nrow(datr)){
  datr[i,44]=100*(datr[i,30]/datr[i,29])
} #Calculate final prevalence of B. glabrata
colnames(datr)[44]<-"bg_prev"
datr$bg_prev[is.na(datr$bg_prev)]<-0

for(i in 1:nrow(datr)){
  datr[i,45]=100*(datr[i,33]/datr[i,32])
} #Calculate final prevalence of B. truncatus
colnames(datr)[45]<-"bt_prev"
datr$bt_prev[is.na(datr$bt_prev)]<-0

#Add atrazine dose data
datr$atra_dose[datr$treat == "A 2x"]<-204
datr$atra_dose[datr$treat != "A 2x" & datr$atra==1]<-102
datr$atra_dose[datr$atra != 1]<-0

#Add chlorpyrifos dose data
datr$chlor_dose[datr$treat == "C 2x"]<-128
datr$chlor_dose[datr$treat != "C 2x" & datr$chlor==1]<-64
datr$chlor_dose[datr$chlor != 1]<-0

#Add phosphorous (fertilizer) dose data
datr$P_dose[datr$treat == "F 2x"]<-880
datr$P_dose[datr$treat != "F 2x" & datr$fert==1]<-440
datr$P_dose[datr$fert != 1]<-0

#Add nitrogen (fertilizer) dose data
datr$N_dose[datr$treat == "F 2x"]<-8800
datr$N_dose[datr$treat != "F 2x" & datr$fert==1]<-4400
datr$N_dose[datr$fert != 1]<-0

varbs<-colnames(datr) #retrieve variable names

#Reshape to long format and merge data set to prepare for plotting #########################
datr2<-reshape(datr, varying = varbs[c(13:42,44:45)],
              v.names="measure",
              timevar="variable",
              times=varbs[c(13:42,44:45)],
              new.row.names=1:1920,
              direction="long")

datr2<-datr2[,-c(1:12,20)] #Get rid of needless variables

aggdata1<-aggregate.data.frame(datr2, by=list(datr2[,1], datr2[,6]), FUN=mean) #calculate means of treatment groups
aggdata1<-aggdata1[,-c(3,4)] #remove unneeded variables 
colnames(aggdata1)<-c("Treatment", "variable", "mean")

aggdata2<-aggregate.data.frame(dat2, by=list(dat2[,1], dat2[,2]), FUN=st.er) #calculate st.error of treatment groups
aggdata2<-aggdata2[,-c(3,4)] #remove unneeded variables 
colnames(aggdata2)<-c("Treatment", "variable", "st.err")  

aggdata<-merge(aggdata1, aggdata2, by=c("Treatment", "variable"))

aggdata$Treatment<-as.character(aggdata$Treatment)

aggdata$Treatment[aggdata$Treatment=="0_0_0"]<- "Control"
aggdata$Treatment[aggdata$Treatment=="1_0_0"]<- "Atrazine"
aggdata$Treatment[aggdata$Treatment=="0_1_0"]<- "ChlorP"
aggdata$Treatment[aggdata$Treatment=="0_0_1"]<- "Fertilizer"
aggdata$Treatment[aggdata$Treatment=="1_1_0"]<- "Atrazine_ChlorP"
aggdata$Treatment[aggdata$Treatment=="1_0_1"]<- "Atrazine_Fertilizer"
aggdata$Treatment[aggdata$Treatment=="0_1_1"]<- "ChlorP_Fertilizer"
aggdata$Treatment[aggdata$Treatment=="1_1_1"]<- "All_Three"

aggdata$Treatment<-as.factor(aggdata$Treatment)
