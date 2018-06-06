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
saro86_dat <- read.csv("Agrochemical_Review/Response_Fxs/Data/Sarojini1986_m_lamerrii.csv")
  saro86_dat$mort <- pnorm(saro86_dat$probit, mean = 5)
  saro86_dat$conc <- 10^saro86_dat$log10ppm
  
#y=mx+b  
saro86_mod <- lm(probit ~ log10ppm, data = saro86_dat)

  saro86_b <- summary(saro86_mod)$coef[1,1]   # Intercept
  saro86_b_se <- summary(saro86_mod)$coef[1,2]  #SE of intercept estimate from model
  saro86_m <- summary(saro86_mod)$coef[2,1]  #slope parameter
    
muPq_fenitrothion_Sarojini86_uncertainty<-function(In){

  b = rnorm(1, saro86_b, saro86_b_se)
  
  mort <- pnorm(saro86_m*log10(In/1000) + b, mean = 5)
  
  mort
}

keep.saro86 <- c("muPq_fenitrothion_Sarojini86_uncertainty", "saro86_b", 
                  "saro86_b_se", "saro86_m")
