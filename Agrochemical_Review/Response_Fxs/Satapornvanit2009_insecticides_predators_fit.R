#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(ggplot2)
require(drc)
require(LW1949)

#24-hr mortality functions for predator (macrobrachium) populations from Satapornvanit et al chemosphere paper ########
  sap.mort<-read.csv("Agrochemical_Review/Response_Fxs/Data/satapornvanit2009_m_rosenbergii_mort.csv")
    sap.mort$dead = round(sap.mort$dead) #Make number dead from each trial an integer
  mort.sub = subset(sap.mort, chem !='carbendazim') 
  
  sap.mupq<-drm(dead/total ~ conc, chem, weights = total,  data = mort.sub, type = 'binomial', 
                fct = LL.3(names = c('b', 'd', 'e'),
                          fixed = c(NA, 1, NA)))
#Zinc ########    
  muPq_zinc_satapornvanit09<-function(In){
    predict(sap.mupq, data.frame(conc=In, chem = 'zinc'), interval = 'confidence', level = 0.95)
  }  
    
  muPq_zinc_satapornvanit09_uncertainty<-function(In){
      init = predict(sap.mupq, data.frame(conc=In, chem = 'zinc'), se.fit = T)
      mup = rnorm(1, init[1], init[2])

      if(mup > 1) mup = 1

    return(mup)
  }
    
#LW1949 analysis
zinc.lw1949 = dataprep(dose = sap.mort$conc[sap.mort$chem == 'zinc'], 
                       ntot = sap.mort$total[sap.mort$chem == 'zinc'], 
                       nfx = sap.mort$dead[sap.mort$chem == 'zinc'])
zinc.lwmod = fitLWauto(zinc.lw1949)
fzinc = LWestimate(zinc.lwmod, zinc.lw1949)

  lc50.sat.zinc = as.numeric(fzinc$LWest[1]) #387.8831 
  slp.sat.zinc = as.numeric(fzinc$params[2]) #1.832589  
  se.sat.zinc = as.numeric(log10(fzinc$LWest[3]/lc50.sat.zinc) / qnorm(0.975)) #0.080841 
  
  muPq_zinc_satapornvanit09_uncertainty_lw49 = function(In){
    lc50 = 10^(rnorm(1, log10(lc50.sat.zinc), se.sat.zinc))
    mun = pnorm(slp.sat.zinc * log10(In/lc50)) 

    return(mun)
  }
  
keep.zinc.sat09 = c('muPq_zinc_satapornvanit09_uncertainty_lw49', 'muPq_zinc_satapornvanit09', 'sap.mupq',
                    'muPq_zinc_satapornvanit09_uncertainty', 'zinc.lwmod.sat09') 

#Chlorpyrifos ###########    
  muPq_chlor_satapornvanit09<-function(In){
    predict(sap.mupq, data.frame(conc=In, chem = 'chlorpyrifos'), 
            interval = 'confidence', level = 0.95)
  }  
    
  muPq_chlor_satapornvanit09_uncertainty<-function(In){
      init = predict(sap.mupq, data.frame(conc=In, chem = 'chlorpyrifos'), se.fit = T)
      mup = rnorm(1, init[1], init[2])
      
    if(mup > 1) mup = 1
      
    return(mup)
  }
      
#LW1949 analysis
chlor.lw1949 = dataprep(dose = sap.mort$conc[sap.mort$chem == 'chlorpyrifos'], 
                       ntot = sap.mort$total[sap.mort$chem == 'chlorpyrifos'], 
                       nfx = sap.mort$dead[sap.mort$chem == 'chlorpyrifos'])

chlor.lwmod = fitLWauto(chlor.lw1949)
fchlor = LWestimate(chlor.lwmod, chlor.lw1949)

  lc50.sat.chlor = as.numeric(fchlor$LWest[1])  
  slp.sat.chlor = as.numeric(fchlor$params[2])   
  se.sat.chlor = as.numeric(log10(fchlor$LWest[3]/lc50.sat.chlor) / qnorm(0.975))  
  
  muPq_chlor_satapornvanit09_uncertainty_lw49 = function(In){
    lc50 = 10^(rnorm(1, log10(lc50.sat.chlor), se.sat.chlor))
    mun = pnorm(slp.sat.chlor * log10(In/lc50)) 

    return(mun)
  }

keep.chlor.sat09 = c('muPq_chlor_satapornvanit09_uncertainty_lw49', 'muPq_chlor_satapornvanit09','sap.mupq',
                     'muPq_chlor_satapornvanit09_uncertainty', 'chlor.lwmod.sat09') 

#Dimethoate ################        
  muPq_dim_satapornvanit09<-function(In){
    predict(sap.mupq, data.frame(conc=In, chem = 'dimethoate'), 
            interval = 'confidence', level = 0.95)
  }  
  
  muPq_dim_satapornvanit09_uncertainty<-function(In){
      init = predict(sap.mupq, data.frame(conc=In, chem = 'dimethoate'), se.fit = T)
      mup = rnorm(1, init[1], init[2])
      if(mup > 1) mup = 1
    
    return(mup)
  }
    
#LW1949 analysis
dim.lw1949 = dataprep(dose = sap.mort$conc[sap.mort$chem == 'dimethoate'], 
                        ntot = sap.mort$total[sap.mort$chem == 'dimethoate'], 
                        nfx = sap.mort$dead[sap.mort$chem == 'dimethoate'])

dim.lwmod = fitLWauto(dim.lw1949)
fdim = LWestimate(dim.lwmod, dim.lw1949)

  lc50.sat.dim = as.numeric(fdim$LWest[1])  
  slp.sat.dim = as.numeric(fdim$params[2])   
  se.sat.dim = as.numeric(log10(fdim$LWest[3]/lc50.sat.dim) / qnorm(0.975))  
  
  muPq_dim_satapornvanit09_uncertainty_lw49 = function(In){
    lc50 = 10^(rnorm(1, log10(lc50.sat.dim), se.sat.dim))
    mun = pnorm(slp.sat.dim * log10(In/lc50)) 

    return(mun)
  }

keep.dim.sat09 = c('muPq_dim_satapornvanit09_uncertainty', 'muPq_dim_satapornvanit09','sap.mupq',
                   'muPq_dim_satapornvanit09_uncertainty_lw49', 'dim.lwmod.sat09')

#Profenofos ################      
  muPq_pr_satapornvanit09<-function(In){
    predict(sap.mupq, data.frame(conc=In, chem = 'profenofos'), interval = 'confidence', level = 0.95)
  }  
      
  muPq_prof_satapornvanit09_uncertainty<-function(In){
      init = predict(sap.mupq, data.frame(conc=In, chem = 'profenofos'), se.fit = T)
      mup = rnorm(1, init[1], init[2])
      if(mup > 1) mup = 1
    
    return(mup)
  }
    
#LW1949 analysis
prof.lw1949 = dataprep(dose = sap.mort$conc[sap.mort$chem == 'profenofos'], 
                      ntot = sap.mort$total[sap.mort$chem == 'profenofos'], 
                      nfx = sap.mort$dead[sap.mort$chem == 'profenofos'])

prof.lwmod = fitLWauto(prof.lw1949)
fprof = LWestimate(prof.lwmod, prof.lw1949)

  lc50.sat.prof = as.numeric(fprof$LWest[1])  
  slp.sat.prof = as.numeric(fprof$params[2])   
  se.sat.prof = as.numeric(log10(fprof$LWest[3]/lc50.sat.prof) / qnorm(0.975))  
  
  muPq_prof_satapornvanit09_uncertainty_lw49 = function(In){
    lc50 = 10^(rnorm(1, log10(lc50.sat.prof), se.sat.prof))
    mun = pnorm(slp.sat.prof * log10(In/lc50)) 

    return(mun)
  }
  
keep.prof.sat09 = c('muPq_prof_satapornvanit09_uncertainty', 'muPq_pr_satapornvanit09','sap.mupq',
                    'muPq_prof_satapornvanit09_uncertainty_lw49', 'prof.lwmod.sat09')
  
#Carbendazim ###########       
  muPq_carb_satapornvanit09<-function(In){ #Paper found no consistent effect of carbendazim on mortality, even
      0*In                                 #at levels higher than the solubility of carbendazim in water, but it's a fungicide, so not necessarily surprising
  }  
  
keep.carb.sat09 = c('muPq_carb_satapornvanit09')

#Reduced feeding rate from zinc and Chlorpyrifos from Satapornvanit et al chemosphere paper ########
sap.fr<-read.csv("Agrochemical_Review/Response_Fxs/Data/satapornvanit2009_m_rosenbergii_feed_rate.csv")
  fr.z = subset(sap.fr, chem == 'zinc')
  fr.ch = subset(sap.fr, chem == 'chlorpyrifos')
  
#Zinc ###################
  zinc.fr= drm(feed_rate ~ conc, data = fr.z, type = 'continuous',
               fct = LL.3(names = c('b', 'd', 'e'),
                          fixed = c(NA, max(fr.z$feed_rate), NA)))

  psi_q_zinc_satapornvanit09<-function(In){
    predict(zinc.fr, data.frame(conc = In), interval = 'confidence', level = 0.95) / fr.z$feed_rate[1]
  }  
    
  par.tricksz = c(coef(zinc.fr), 'd' = max(fr.z$feed_rate))[c(1,3,2)]
    
  psiq_zinc_satapornvanit09_uncertainty<-function(In){
    rdrm(nosim = 1, fct = LL.3(), mpar = par.tricksz, yerror = 'rnorm', xerror = In,
         ypar = c(0, predict(zinc.fr, data.frame(dose = In), se.fit = T)[2]))$y / fr.z$feed_rate[1]
  }
    
keep.zinc.sat09 = c(keep.zinc.sat09, 'psiq_zinc_satapornvanit09_uncertainty', 
                    'par.tricksz', 'zinc.fr', 'fr.z')    
#Chlorpyrifos ###################
  chlor.fr= drm(feed_rate ~ conc, data = fr.ch, type = 'continuous',
               fct = LL.3(names = c('b', 'd', 'e'),
                          fixed = c(NA, max(fr.ch$feed_rate), NA)))

  psi_q_chlor_satapornvanit09<-function(In){
    predict(chlor.fr, data.frame(conc = In), interval = 'confidence', level = 0.95) / fr.ch$feed_rate[1]
  }  
  
  par.tricksc = c(coef(chlor.fr), 'd' = max(fr.ch$feed_rate))[c(1,3,2)]
    
    psiq_chlor_satapornvanit09_uncertainty<-function(In){
      rdrm(nosim = 1, fct = LL.3(), mpar = par.tricksc, yerror = 'rnorm', xerror = In,
           ypar = c(0, predict(chlor.fr, data.frame(dose = In), se.fit = T)[2]))$y / fr.ch$feed_rate[1]
    }
    
keep.chlor.sat09 = c(keep.chlor.sat09, 'psiq_chlor_satapornvanit09_uncertainty', 
                    'par.tricksc', 'chlor.fr', 'fr.ch')    

keep.all.sat09 = c(keep.carb.sat09, keep.chlor.sat09, keep.dim.sat09, keep.prof.sat09, keep.zinc.sat09)
#*Note: checked that this produced the desired relationship on predator feeding rate and it does