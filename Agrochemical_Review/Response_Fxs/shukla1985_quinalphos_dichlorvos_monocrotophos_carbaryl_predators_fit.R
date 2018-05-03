#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to shukla and murti 1985 data found table 1
#quinalphos
lc50.shukla.quin <- 0.796
se.lc50.shukla.quin <- log10(0.884/0.796) / qnorm(0.975)

slp.shukla.quin <- coef(lm(qnorm(c(.25,0.5,0.75)) ~ log10(c(0.655, 0.796, 0.969))))[2]

muPq_quin_shukla85_uncertainty <- function(In){
  lc50 = 10^(rnorm(1, log10(lc50.shukla.quin), se.lc50.shukla.quin)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.shukla.quin * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

keep.shukla <- c("muPq_quin_shukla85_uncertainty", "slp.shukla.quin", "lc50.shukla.quin", "se.lc50.shukla.quin")

#dichlorvos
lc50.shukla.dich <- 1.435
se.lc50.shukla.dich <- log10(1.581/1.435) / qnorm(0.975)

slp.shukla.dich <- coef(lm(qnorm(c(.25,0.5,0.75)) ~ log10(c(1.194, 1.435, 1.726))))[2]

muPq_dich_shukla85_uncertainty <- function(In){
  lc50 = 10^(rnorm(1, log10(lc50.shukla.dich), se.lc50.shukla.dich)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.shukla.dich * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

keep.shukla <- c(keep.shukla, "muPq_dich_shukla85_uncertainty", "slp.shukla.dich", "lc50.shukla.dich", "se.lc50.shukla.dich")

#monocrotophos
lc50.shukla.mono <- 2.107
se.lc50.shukla.mono <- log10(2.186/2.107) / qnorm(0.975)

slp.shukla.mono <- coef(lm(qnorm(c(.25,0.5,0.75)) ~ log10(c(1.963, 2.107, 2.261))))[2]

muPq_mono_shukla85_uncertainty <- function(In){
  lc50 = 10^(rnorm(1, log10(lc50.shukla.mono), se.lc50.shukla.mono)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.shukla.mono * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

keep.shukla <- c(keep.shukla, "muPq_mono_shukla85_uncertainty", "slp.shukla.mono", "lc50.shukla.mono", "se.lc50.shukla.mono")


#Carbaryl
lc50.shukla.carb <- 0.033
se.lc50.shukla.carb <- log10(0.036/0.033) / qnorm(0.975)

slp.shukla.carb <- coef(lm(qnorm(c(.25,0.5,0.75)) ~ log10(c(0.029, 0.033, 0.038))))[2]

muPq_carb_shukla85_uncertainty <- function(In){
  lc50 = 10^(rnorm(1, log10(lc50.shukla.carb), se.lc50.shukla.carb)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.shukla.carb * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

keep.shukla <- c(keep.shukla, "muPq_carb_shukla85_uncertainty", "slp.shukla.carb", "lc50.shukla.carb", "se.lc50.shukla.carb")