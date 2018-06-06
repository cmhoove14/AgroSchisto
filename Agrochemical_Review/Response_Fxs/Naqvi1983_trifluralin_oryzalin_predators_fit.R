#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Naqvi 1983 P. clarkii toxicity
require(drc)

#trifluralin ##########
naq83_trif_dat <- data.frame(conc = c(0,1,3,5,6,7,9,12,15,18,20), #ppm
                             mort = c(0,0,0,3,12,17,30,37,57,80,100)/100,
                             tot = c(30,30,30,35,40,30,40,30,40,30,40))

naq83_trif_dat$dead <- round(naq83_trif_dat$mort * naq83_trif_dat$tot)

#DRC model  
naq83_trif_mod <- drm(dead / tot ~ conc, 
                     weights = tot, type = 'binomial',  
                     data = naq83_trif_dat, fct = LL2.2())
  
  naq83_trif_lc50 <- summary(naq83_trif_mod)$coef[2,1]   # Lc50 estimate from model
  naq83_trif_lc50_se <- summary(naq83_trif_mod)$coef[2,2]  #SE of lc50 estimate from model
  naq83_trif_b <- summary(naq83_trif_mod)$coef[1,1]  #slope parameter of log-logistic fitted here
    
muPq_trifluralin_Naqvi83_uncertainty<-function(In){

  e = rnorm(1, naq83_trif_lc50, naq83_trif_lc50_se)
  
  mort <- 1 / (1 + exp(naq83_trif_b*(log(In/1000)-e)))
  
  mort
}

#oryzalin ##########
naq83_oryz_dat <- data.frame(conc = c(0,10,100,200,300,400,500,1000,
                                      1500,2000,5000,8000,10000), #ppm
                             mort = c(0,0,0,0,0,0,0,6,6,19,20,20,20)/100,
                             tot = rep(40,13))

naq83_oryz_dat$dead <- round(naq83_oryz_dat$mort * naq83_oryz_dat$tot)

#DRC model  
naq83_oryz_mod <- drm(dead / tot ~ conc, 
                     weights = tot, type = 'binomial',  
                     data = naq83_oryz_dat, fct = LL2.2())
  
  naq83_oryz_lc50 <- summary(naq83_oryz_mod)$coef[2,1]   # Lc50 estimate from model
  naq83_oryz_lc50_se <- summary(naq83_oryz_mod)$coef[2,2]  #SE of lc50 estimate from model
  naq83_oryz_b <- summary(naq83_oryz_mod)$coef[1,1]  #slope parameter of log-logistic fitted here
    
muPq_oryzalin_Naqvi83_uncertainty<-function(In){
  b = naq83_oryz_b
  lc50 = naq83_oryz_lc50
  lc50_se = naq83_oryz_lc50_se
  
  e = rnorm(1, naq83_oryz_lc50, naq83_oryz_lc50_se)
  
  mort <- 1 / (1 + exp(naq83_oryz_b*(log(In/1000)-e)))
  
  mort
}

keep.naqvi83 <- c("muPq_trifluralin_Naqvi83_uncertainty", "naq83_trif_b", 
                  "naq83_trif_lc50", "naq83_trif_lc50_se",
                  "muPq_oryzalin_Naqvi83_uncertainty", "naq83_oryz_b", 
                  "naq83_oryz_lc50", "naq83_oryz_lc50_se")
