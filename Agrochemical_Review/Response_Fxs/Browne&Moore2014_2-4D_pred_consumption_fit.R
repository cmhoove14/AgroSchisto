#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

# Model with beta regression since data reported is percent consumption
require(betareg)

#Reduced predator consumption rate due to 2,4-D exposure to crayfish Orconectus rusticus
browne_dat <- data.frame("per_con" = c(.144, 0.062, 0.057, -0.024),
                         "se" = c(.177-.144, 0.084-0.062, 0.083-0.057, 0.0049--0.024), 
                         "conc" = c(0, 3.75, 14.07, 32.69))
#Can't deal with a negative percent which they say is basically that the food gained more weight by absorbing water than what the crayfish consumed, so we'll call that 0 and normalize everything to it

browne_dat$per_con <- browne_dat$per_con + abs(min(browne_dat$per_con)) + 1e-6

browne_ref <- browne_dat$per_con[1]
browne_ref_se <- browne_dat$se[1]

#2,4D ###################
  browne_mod= betareg(per_con ~ conc, weights = 1/se, data = browne_dat)

  psiq_24D_browne14<-function(He, ref = browne_ref){
    predict(browne_mod, data.frame(conc = He/1000)) / ref
  }  
    
  psiq_24D_browne14_uncertainty<-function(He){
    predict(browne_mod, newdata = data.frame(conc = He/1000)) / rnorm(1, browne_ref, browne_ref_se)
  }
    
keep.24d.browne14 = c('psiq_24D_browne14_uncertainty', 'browne_mod', 'browne_ref', "browne_ref_se")    

