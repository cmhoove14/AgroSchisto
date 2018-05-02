#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Monde et al 2016 data found in tabl2 2 for bulinus globusus 

lc50.monde.endo <- 4160
se.lc50.monde.endo <- log10(12956/4160) / qnorm(0.975)

slp.monde.endo <- coef(lm(qnorm(c(.01,0.5,0.9)) ~ log10(c(807, 4160, 21457))))[2]

muNq_endo_monde16_uncertainty <- function(In){
  lc50 = 10^(rnorm(1, log10(lc50.monde.endo), se.lc50.monde.endo)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.monde.endo * log10(In/lc50)) #Estimate daily mortality (percent)

  return(as.numeric(mun))

}

keep.monde2016 <- c("muNq_endo_monde16_uncertainty", "slp.monde.endo", "lc50.monde.endo", "se.lc50.monde.endo")