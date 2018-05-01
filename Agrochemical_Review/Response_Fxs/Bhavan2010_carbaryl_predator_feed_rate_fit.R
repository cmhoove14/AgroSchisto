#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(drc)

#Reduced feeding rate from carbaryl from Bhavan et al paper ########
bhavan_dat <- data.frame("food_con" = c(33.51, 28.47, 19.30, 13.32),
                         "error" = c(1.47, 1.44, 1.59, 1.04),
                         "carbaryl" = c(0, 5.15, 7.73, 15.47))
  
#Zinc ###################
  carb_fr= drm(food_con ~ carbaryl, data = bhavan_dat, type = 'continuous',
               fct = LL.3(names = c('b', 'd', 'e'),
                          fixed = c(NA, max(bhavan_dat$food_con), NA)))

  psi_q_carb_bhavan10<-function(In){
    predict(carb_fr, data.frame(conc = In), interval = 'confidence', level = 0.95) / bhavan_dat$food_con[1]
  }  
    
  par.tricks.carb = c(coef(carb_fr), 'd' = max(bhavan_dat$food_con))[c(1,3,2)]
    
  psi_q_carb_bhavan10_uncertainty<-function(In){
    rdrm(nosim = 1, fct = LL.3(), mpar = par.tricks.carb, yerror = 'rnorm', xerror = In,
         ypar = c(0, predict(carb_fr, data.frame(dose = In), se.fit = T)[2]))$y / bhavan_dat$food_con[1]
  }
    
keep.carb.bhavan10 = c('psi_q_carb_bhavan10_uncertainty', 'par.tricks.carb', 'carb_fr')    
