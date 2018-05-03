#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to omkar and murti 1985 data found table 1
#Endosulfan
lc50.omkar.endo <- 0.0062
se.lc50.omkar.endo <- log10(0.0068/0.0062) / qnorm(0.975)

slp.omkar.endo <- coef(lm(qnorm(c(.25,0.5,0.75)) ~ log10(c(0.0052, 0.0062, 0.0073))))[2]

muPq_endo_omkar85_uncertainty <- function(In){
  lc50 = 10^(rnorm(1, log10(lc50.omkar.endo), se.lc50.omkar.endo)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.omkar.endo * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

keep.omkar <- c("muPq_endo_omkar85_uncertainty", "slp.omkar.endo", "lc50.omkar.endo", "se.lc50.omkar.endo")

#Phosphamidon
lc50.omkar.phos <- 4.825
se.lc50.omkar.phos <- log10(5.429/4.825) / qnorm(0.975)

slp.omkar.phos <- coef(lm(qnorm(c(.25,0.5,0.75)) ~ log10(c(3.783, 4.825, 6.153))))[2]

muPq_phos_omkar85_uncertainty <- function(In){
  lc50 = 10^(rnorm(1, log10(lc50.omkar.phos), se.lc50.omkar.phos)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.omkar.phos * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

keep.omkar <- c(keep.omkar, "muPq_phos_omkar85_uncertainty", "slp.omkar.phos", "lc50.omkar.phos", "se.lc50.omkar.phos")

#Carbaryl
lc50.omkar.carb <- 0.0513
se.lc50.omkar.carb <- log10(0.0541/0.0513) / qnorm(0.975)

slp.omkar.carb <- coef(lm(qnorm(c(.25,0.5,0.75)) ~ log10(c(0.0461, 0.0513, 0.0572))))[2]

muPq_carb_omkar85_uncertainty <- function(In){
  lc50 = 10^(rnorm(1, log10(lc50.omkar.carb), se.lc50.omkar.carb)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.omkar.carb * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

keep.omkar <- c(keep.omkar, "muPq_carb_omkar85_uncertainty", "slp.omkar.carb", "lc50.omkar.carb", "se.lc50.omkar.carb")