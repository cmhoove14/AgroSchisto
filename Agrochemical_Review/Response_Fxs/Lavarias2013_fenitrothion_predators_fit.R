#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Fenitrothion ##########
lav13_dat <- read.csv("Agrochemical_Review/Response_Fxs/Data/Lavarias2013_m_borelli_fenitrothion.csv")
  lav13_dat$probit <- qnorm(lav13_dat$mort/100)
  lav13_dat$log10ppb <- log10(lav13_dat$conc)
  
#y=mx+b  
lav13_mod <- lm(probit ~ log10ppb, data = lav13_dat[-1,])

  lav13_b <- summary(lav13_mod)$coef[1,1]   # Intercept
  lav13_b_se <- summary(lav13_mod)$coef[1,2]  #SE of intercept estimate from model
  lav13_m <- summary(lav13_mod)$coef[2,1]  #slope parameter
    
muPq_fenitrothion_Lavarias13_uncertainty<-function(In){

  b = rnorm(1, lav13_b, lav13_b_se)
  
  mort <- pnorm(lav13_m*log10(In) + b)
  
  mort
}

keep.lav13 <- c("muPq_fenitrothion_Lavarias13_uncertainty", "lav13_b", 
                "lav13_b_se", "lav13_m")
