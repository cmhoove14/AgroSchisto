#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Leung 1980 P. clarkii fenitrothion toxicity
require(drc)

#Adults ##########
leung_dat <- data.frame(conc = c(0, 15, 25, 35, 50, 75, 100), #ppm
                        mort = c(0, 0, 0, 1.9, 4.0, 5.6, 17)/100,
                        tot = c(115, 30, 50, 54, 40, 54, 30))

leung_dat$dead <- round(leung_dat$mort * leung_dat$tot)

#DRC model  
leung_fen_mod <- drm(dead / tot ~ conc, 
                     weights = tot, type = 'binomial',  
                     data = leung_dat, fct = LL2.2())
  
  leung_fen_lc50 <- summary(leung_fen_mod)$coef[2,1]   # Lc50 estimate from model
  leung_fen_lc50_se <- summary(leung_fen_mod)$coef[2,2]  #SE of lc50 estimate from model
  leung_fen_b <- summary(leung_fen_mod)$coef[1,1]  #slope parameter of log-logistic fitted here
    
muPq_fenitrothion_Leung_uncertainty<-function(In){
  b = leung_fen_b
  lc50 = leung_fen_lc50
  lc50_se = leung_fen_lc50_se
  
  e = rnorm(1, lc50, lc50_se)
  
  mort <- 1 / (1 + exp(b*(log(In/1000)-e)))
  
  mort
}

#Juveniles ##########
leung_dat_juv <- data.frame(conc = c(0, 0.5, 2, 4, 7, 8), #ppm
                            mort = c(3.7, 1.6, 8.1, 19, 38, 47)/100,
                            tot = c(134, 64, 62, 62, 32, 60))

leung_dat_juv$dead <- round(leung_dat_juv$mort * leung_dat_juv$tot)

#DRC model  
leung_juv_fen_mod <- drm(dead / tot ~ conc, 
                         weights = tot, type = 'binomial',  
                         data = leung_dat_juv, fct = LL2.2())
  
  leung_juv_fen_lc50 <- summary(leung_juv_fen_mod)$coef[2,1]   # Lc50 estimate from model
  leung_juv_fen_lc50_se <- summary(leung_juv_fen_mod)$coef[2,2]  #SE of lc50 estimate from model
  leung_juv_fen_b <- summary(leung_juv_fen_mod)$coef[1,1]  #slope parameter of log-logistic fitted here
    
muPq_juv_fenitrothion_Leung_uncertainty<-function(In){
  b = leung_juv_fen_b
  lc50 = leung_juv_fen_lc50
  lc50_se = leung_juv_fen_lc50_se
  
  e = rnorm(1, lc50, lc50_se)
  
  mort <- 1 / (1 + exp(b*(log(In/1000)-e)))
  
  mort
}


keep.leung13 <- c("muPq_fenitrothion_Leung_uncertainty", "leung_fen_b", 
                  "leung_fen_lc50", "leung_fen_lc50_se", "muPq_juv_fenitrothion_Leung_uncertainty",
                  "leung_juv_fen_b", "leung_juv_fen_lc50", "leung_juv_fen_lc50_se")
