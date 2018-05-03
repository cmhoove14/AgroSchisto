#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to shukla and omkar 1984 data found table 1
#endosulfan
lc50.shuk84.endo <- 0.00534
se.lc50.shuk84.endo <- log10(0.00585/0.00534) / qnorm(0.975)

slp.shuk84.endo <- coef(lm(qnorm(c(.25,0.5,0.75)) ~ log10(c(0.00453, 0.00534, 0.00631))))[2]

muPq_endo_shuk84_uncertainty <- function(In){
  lc50 = 10^(rnorm(1, log10(lc50.shuk84.endo), se.lc50.shuk84.endo)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.shuk84.endo * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

keep.shuk84 <- c("muPq_endo_shuk84_uncertainty", "slp.shuk84.endo", "lc50.shuk84.endo", "se.lc50.shuk84.endo")

#methyl-demeton
lc50.shuk84.mede <- 4.156
se.lc50.shuk84.mede <- log10(4.594/4.156) / qnorm(0.975)

slp.shuk84.mede <- coef(lm(qnorm(c(.25,0.5,0.75)) ~ log10(c(3.455, 4.156, 4.998))))[2]

muPq_mede_shuk84_uncertainty <- function(In){
  lc50 = 10^(rnorm(1, log10(lc50.shuk84.mede), se.lc50.shuk84.mede)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.shuk84.mede * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

keep.shuk84 <- c(keep.shuk84, "muPq_mede_shuk84_uncertainty", "slp.shuk84.mede", "lc50.shuk84.mede", "se.lc50.shuk84.mede")

#Carbaryl
lc50.shuk84.carb <- 0.0489
se.lc50.shuk84.carb <- log10(0.0518/0.0489) / qnorm(0.975)

slp.shuk84.carb <- coef(lm(qnorm(c(.25,0.5,0.75)) ~ log10(c(0.0439, 0.0489, 0.0545))))[2]

muPq_carb_shuk84_uncertainty <- function(In){
  lc50 = 10^(rnorm(1, log10(lc50.shuk84.carb), se.lc50.shuk84.carb)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.shuk84.carb * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

keep.shuk84 <- c(keep.shuk84, "muPq_carb_shuk84_uncertainty", "slp.shuk84.carb", "lc50.shuk84.carb", "se.lc50.shuk84.carb")