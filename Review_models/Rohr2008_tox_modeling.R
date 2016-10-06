##This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Parameterizing cercarial mortality from agrochemical exposure
require(drc)

#Rohr et al 2008 ########################
rohr<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Rohr2008/rohr2008_all.csv')
  rohrchems<-unique(as.character(rohr$chem))
  rohr$surv = rohr$surv/100
  
  plot(rohr$time_hrs[rohr$chem == rohrchems[1]], rohr$surv[rohr$chem == rohrchems[1]], ylim = c(0,1),
       type = 'l', lwd=2, xlab = 'time (hrs)', ylab = '%survival')
  for(i in 2:length(rohrchems)){
    lines(rohr$time_hrs[rohr$chem == rohrchems[i]], rohr$surv[rohr$chem == rohrchems[i]], col = i+1, lwd=2)
  }
  
cmod <- glm(rohr$surv[rohr$chem == rohrchems[1]] ~ rohr$time_hrs[rohr$chem == rohrchems[1]], family = binomial(link = 'logit'))  