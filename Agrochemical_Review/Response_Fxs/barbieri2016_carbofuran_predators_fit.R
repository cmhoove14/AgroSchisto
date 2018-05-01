#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(LW1949)
#Toxicity of Carbofuran to Macrobrachium olfersii from Barbieri et al 2016

barb_dat <- dataprep(dose = c(0, 0.01, 0.1, 0.5, 1, 2, 3, 4),
                     ntot = rep(45*3, 8),
                     nfx = round(c(0, 0, 0, 6.66, 20, 33.33, 100, 100)/100 * (45*3)))
#DRC analysis
  barb.mupq<-drm(nfx/ntot ~ dose, weights = ntot,  data = barb_dat, type = 'binomial', 
                fct = LL.2(names = c('b', 'e'),
                          fixed = c(NA, NA)))
  
    summary(barb.mupq)

  muPq_carb_barb16_uncertainty<-function(In){
      init = predict(barb.mupq, data.frame(conc=In/1000), se.fit = T)
      mup = rnorm(1, init[1], init[2])

      if(mup > 1) mup = 1

    return(mup)
  }

#LW1949 analysis
barb_slp <- fitLWauto(barb_dat)[2]
barb_lw_pars <- LWestimate(fitLWauto(barb_dat), barb_dat)

  lc50.barb.carb = as.numeric(barb_lw_pars$LWest[1])  
  slp.barb.carb = as.numeric(barb_lw_pars$params[2])   
  se.barb.carb = as.numeric(log10(barb_lw_pars$LWest[3]/lc50.barb.carb) / qnorm(0.975))  
  
  barb_carbofuran_muPq_uncertainty = function(In){
    ins = In/1000
    lc50 = 10^(rnorm(1, log10(lc50.barb.carb), se.barb.carb))
    mun = pnorm(slp.barb.carb * log10(ins/lc50)) 

    return(mun)
  }

  keep.barb2016 <- c("barb_carbofuran_muPq_uncertainty", "slp.barb.carb", "lc50.barb.carb", "se.barb.carb",
                     "muPq_carb_barb16_uncertainty", "barb.mupq")