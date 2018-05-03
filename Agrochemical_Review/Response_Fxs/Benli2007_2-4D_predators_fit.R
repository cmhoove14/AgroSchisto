#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Benli et al 2007 data found table 1

benli_morts <- c(1,5,10,15,50,85,90,95,99)/100
benli_lcs <- c(3.06, 6.08, 8.78, 11.24, 32.6, 91.06, 116.62, 168.27, 334.69)

#2,4D
lc50.benli.24d <- 32.6
se.lc50.benli.24d <- log10(327.16/32.6) / qnorm(0.975)

slp.benli.24d <- 2.281

muPq_24d_benli07 <- function(In, lc50 = lc50.benli.24d){
  mun = pnorm(slp.benli.24d * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

muPq_24d_benli07_uncertainty <- function(In){
  lc50 = 10^(rnorm(1, log10(lc50.benli.24d), se.lc50.benli.24d)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.benli.24d * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

keep.benli <- c("muPq_24d_benli07_uncertainty", "slp.benli.24d", "lc50.benli.24d", "se.lc50.benli.24d")
