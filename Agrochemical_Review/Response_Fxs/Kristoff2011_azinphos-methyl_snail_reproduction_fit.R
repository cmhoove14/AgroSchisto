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

#Reduced snail reproduction in Bi. glabrata exposed to azinphos-methyl from Kristoff et al  ########
#Data modeled here is just number of egg masses. Eggs per mass, hatchings per mass, and time to hatching are also reported,
#but no variation in these implies that variation in egg masses is only thing that affects parameter value (in units of hatchlings/snail/day)
kristoff_dat <- data.frame("egg_mass" = c(30, 30, 24, 22),
                           "azmethyl" = c(0, 0.5, 2.5, 5))

kristoff_ref <- kristoff_dat$egg_mass[1]
  
#carbaryl ###################
  azmeth_mod= drm(egg_mass ~ azmethyl, data = kristoff_dat, type = 'continuous',
                  fct = LL.3(names = c('b', 'd', 'e'),
                             fixed = c(NA, max(kristoff_dat$egg_mass), NA)))

  fNq_azmeth_kristoff11<-function(In){
    ins = In/1000
    predict(azmeth_mod, data.frame(azmethyl = ins), interval = 'confidence', level = 0.95) / kristoff_ref
  }  
    
  par.tricks.azmeth = c(coef(azmeth_mod), 'd' = max(kristoff_dat$egg_mass))[c(1,3,2)]
    
  fNq_azmeth_kristoff11_uncertainty<-function(In){
    ins = In/1000
    rdrm(nosim = 1, fct = LL.3(), mpar = par.tricks.azmeth, yerror = 'rnorm', xerror = ins,
         ypar = c(0, predict(azmeth_mod, data.frame(azmethyl = ins), se.fit = T)[2]))$y / kristoff_ref
  }
    
keep.azmeth.kristoff11 = c('fNq_azmeth_kristoff11_uncertainty', 'par.tricks.azmeth', 'azmeth_mod')    
