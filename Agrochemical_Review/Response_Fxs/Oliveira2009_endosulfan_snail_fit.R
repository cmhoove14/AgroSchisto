#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#data and d-r response from Oliveira et al 2009
#eggs per snail over eight weeks
oliv_dat <- data.frame(endo = c(0, 0.001, 0.01, 0.1), 
                       eggs = c(580.1, 580.1, 490, 335.4),
                       eggs_se = c(33.9, 33.9, 31.1, 28.3), #Assume variance proportional to mean 
                       surv = c(0.98,0.93,0.92,0.76), 
                       hatch = c(0.91, 0.99, 0.62, 0.24))   

oliv_dat$net <- oliv_dat$eggs * oliv_dat$surv * oliv_dat$hatch

oliv_dat_ref <- oliv_dat$net[1]
#Estimate d-r function with drc: eggs/snail/week over entire study period as a function of endosulfan concentration
  endo_repro_mod= drm(net ~ endo, data = oliv_dat, type = 'continuous',
                      fct = LL.3(names = c('b', 'd', 'e'),
                                 fixed = c(NA, max(oliv_dat$net), NA)))

  fNq_endo_oliv09<-function(In){
    predict(endo_repro_mod, data.frame(conc = In/1000), interval = 'confidence', level = 0.95)
  }  
    
  par.tricks.endo = c(coef(endo_repro_mod), 'd' = max(oliv_dat$net))[c(1,3,2)]
    
  fNq_endo_oliv09_uncertainty<-function(In){
    fn <- rdrm(nosim = 1, fct = LL.3(), mpar = par.tricks.endo, yerror = 'rnorm', xerror = In/1000,
               ypar = c(0, predict(endo_repro_mod, data.frame(conc = In/1000), se.fit = T)[2]))$y / oliv_dat_ref
    
    if(fn < 0){
      return(0)
    } else {
      return(fn)
    }
  }
  
keep.oliv09 <- c("fNq_endo_oliv09_uncertainty", "par.tricks.endo", "endo_repro_mod", "oliv_dat_ref")