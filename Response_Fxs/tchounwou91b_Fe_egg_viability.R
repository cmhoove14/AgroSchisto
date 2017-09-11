#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(drc)

#miracidial mortality (S. mansoni) from Tchounwou 1991 ####################################
eggv = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/tchounwou91_urea_ammP_schistosome_egg_viability.csv')
  eggv.amm = subset(eggv, chem == 'amm_sulph')
  eggv.ure = subset(eggv, chem == 'urea')

#Amm. sulphate model    
tch91.egv.amm<-drm(hatch/total ~ conc, weights = total, data = eggv.amm, type = 'binomial', 
                   fct = LL.2())
  summary(tch91.egv.amm)
  plot(tch91.egv.amm) 
  
tch91.egv.amm_unc = function(Fe){
  v = sum(rdrm(18, LL.2(), coef(tch91.egv.amm), Fe, yerror = 'rbinom', ypar = eggv.amm$total, onlyY = T)$y) / 
        sum(eggv.amm$total) #Simulated number of expected hatches at concentration /  total hatches with no failures
                              #same as hatches in control group, 100%
  v
}  

points(seq(0,1e7,1e4), sapply(seq(0,1e7,1e4), tch91.egv.amm_unc), pch = 5, col = 2, cex = 0.6)

#Urea model    
tch91.egv.ure<-drm(hatch/total ~ conc, weights = total, data = eggv.ure, type = 'binomial', 
                   fct = LL.2())
summary(tch91.egv.ure)
plot(tch91.egv.ure) 

tch91.egv.ure_unc = function(Fe){
  v = sum(rdrm(27, LL.2(), coef(tch91.egv.ure), Fe, yerror = 'rbinom', ypar = eggv.ure$total, onlyY = T)$y) / 
    sum(eggv.ure$total) #Simulated number of expected hatches at concentration /  total hatches with no failures
  #same as hatches in control group, 100%
  v
}  

points(seq(0,5e7,5e4), sapply(seq(0,5e7,5e4), tch91.egv.ure_unc), pch = 5, col = 2, cex = 0.6)


keep.tch91.egv = c('tch91.egv.amm_unc', 'tch91.egv.amm', 'eggv.amm', 
                   'tch91.egv.ure_unc', 'tch91.egv.ure', 'eggv.ure')