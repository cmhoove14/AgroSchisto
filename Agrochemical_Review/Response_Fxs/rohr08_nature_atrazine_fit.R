#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Rohr '08 nature paper comparison of control to atrazine treated mesocosms
#snail reproduction parameter from table 1

rohr08_fNq_uncertainty = function(){ #function based on snail egg masses
  fN = rnorm(1, 30.45, 3.03) / 12.37
  while(fN < 0) fN = rnorm(1, 30.45, 3.03) / 12.37
  fN
}

rohr08_fNq_uncertainty2 = function(){ #function based on snail hatchlings
  fN = rnorm(1, 468.09, 205.85) / 110.64
  while(fN < 0) fN = rnorm(1, 468.09, 205.85) / 110.64
  fN
}

#hist(replicate(1000, rohr08_fN_uncertainty()), breaks = 30)
#hist(replicate(1000, rohr08_fN_uncertainty2()), breaks = 30)