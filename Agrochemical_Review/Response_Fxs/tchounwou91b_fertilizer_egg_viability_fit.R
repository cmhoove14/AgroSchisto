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

#egg viability from Tchounwou 1991 ####################################
eggv = read.csv('Agrochemical_Review/Response_Fxs/Data/tchounwou91_urea_ammP_schistosome_egg_viability.csv')
eggv$conc = eggv$conc/1000
  eggv.amm = subset(eggv, chem == 'amm_sulph')
  eggv.ure = subset(eggv, chem == 'urea')

#Amm. sulphate model ########## 
tch91.egv.amm<-drm(hatch/total ~ conc, weights = total, data = eggv.amm, type = 'binomial', 
                   fct = LL.3(names = c('b', 'd', 'e'),
                              fixed = c(NA, 1, NA)))

  tch91_amm_lc50 <- summary(tch91.egv.amm)$coefficients[2,1]
    tch91_amm_lc50_se <- summary(tch91.egv.amm)$coefficients[2,2]
  tch91_amm_slp <- summary(tch91.egv.amm)$coefficients[1,1]
  
tch91_amm_v_unc<-function(Fe){
  Fer = (Fe/1000)
  lc50 = rnorm(1, tch91_amm_lc50, tch91_amm_lc50_se)
  
  v = 1 / (1+exp(tch91_amm_slp*log(Fer / lc50)))

  return(v)
}

keep.tch91.fe.egg <- c("tch91_amm_lc50", "tch91_amm_lc50_se", "tch91_amm_slp", "tch91_amm_v_unc")

#Urea model #############
tch91.egv.ure<-drm(hatch/total ~ conc, weights = total, data = eggv.ure, type = 'binomial', 
                   fct = LL.3(names = c('b', 'd', 'e'),
                              fixed = c(NA, 1, NA)))

  tch91_ure_lc50 <- summary(tch91.egv.ure)$coefficients[2,1]
    tch91_ure_lc50_se <- summary(tch91.egv.ure)$coefficients[2,2]
  tch91_ure_slp <- summary(tch91.egv.ure)$coefficients[1,1]
  
tch91_ure_v_unc<-function(Fe){
  Fer = (Fe/1000)
  lc50 = rnorm(1, tch91_ure_lc50, tch91_ure_lc50_se)
  
  fN = 1 / (1+exp(tch91_ure_slp*log(Fer / lc50)))

  return(fN)
}

keep.tch91.fe.egg <- c("tch91_ure_lc50", "tch91_ure_lc50_se", "tch91_ure_slp", "tch91_ure_v_unc")
