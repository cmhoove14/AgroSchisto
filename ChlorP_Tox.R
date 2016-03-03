#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

cray<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/Cray.LC50.2012.csv')
  cray.chlor<-subset(cray, Chem=='Chlor')
  cc.agg<-aggregate.data.frame(cray.chlor, by=list(cray.chlor[,2]), FUN=mean)
  cc.agg<-cc.agg[,-c(1,2)]
belo<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/Belo.LC50.2012.csv')
  belo.chlor<-subset(belo, Chem=='Chlor')
  bc.agg<-aggregate.data.frame(belo.chlor, by=list(belo.chlor[,2]), FUN=mean)
  bc.agg<-bc.agg[,-c(1,2)]

