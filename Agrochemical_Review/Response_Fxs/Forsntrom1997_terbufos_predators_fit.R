#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Fornstrom 1997 data from Fig 1 at 24 hrs
require(drc)

forn_dat <- data.frame(conc = c(0, 0.73, 1.5, 3.2, 7.4, 18),
                       mort = c(0, 0, 0, 11, 22, 66)/100,
                       dead = c(0, 0, 0, 1, 2, 6),
                       tot = rep(9, 6))

#DRC model  
forn_terb_mod <- drm(dead / tot ~ conc, 
                     weights = tot, type = 'binomial',  
                     data = forn_dat, fct = LL2.2())
  
  forn_terb_lc50 <- summary(forn_terb_mod)$coef[2,1]   # Lc50 estimate from model
  forn_terb_lc50_se <- summary(forn_terb_mod)$coef[2,2]  #SE of lc50 estimate from model
  forn_terb_b <- summary(forn_terb_mod)$coef[1,1]  #slope parameter of log-logistic fitted here
    
muPq_terb_Fornstrom_uncertainty<-function(In){
  b = forn_terb_b
  lc50 = forn_terb_lc50
  lc50_se = forn_terb_lc50_se
  
  e = rnorm(1, lc50, lc50_se)
  
  mort <- 1 / (1 + exp(b*(log(In)-e)))
  
  mort
}

keep.forn97 <- c("muPq_terb_Fornstrom_uncertainty", "forn_terb_b", 
                 "forn_terb_lc50", "forn_terb_lc50_se")
