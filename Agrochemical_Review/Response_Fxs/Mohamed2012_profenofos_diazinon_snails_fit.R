#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Data extraction and model fitting to Mohamed et al 2012 data
source("Agrochemical_Review/Models/litchfield_wilcoxon_get_b1_from_slope.R")

#Basudin (Diazinon) reported LC50 and slope data #####################
  lc50.moh.diaz.report = 14.16
  slp.moh.diaz.report = 1.22
  b1.moh.diaz = get_b1(slp.moh.diaz.report)
  #SE not reported, so borrow from upper range of SEs reported from other studies
    se.lc50.moh.diaz = 0.12
  
#Function to plot estimates without uncertainty
muNq_diaz_mohamed = function(In, lc50 = lc50.moh.diaz.report){
  Ins = (In/1000) #Parameters based on ppm, data input as ppb

    mun = pnorm(b1.moh.diaz * log10(Ins/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }
    
#Create function based on reverse of litchfield and wilcoxon      
muNq_diaz_mohamed_uncertainty = function(In){
  Ins = (In/1000) #Parameters based on ppm, data input as ppb
  lc50 = 10^(rnorm(1, log10(lc50.moh.diaz.report), se.lc50.moh.diaz)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
    mun = pnorm(b1.moh.diaz * log10(Ins/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }
  
#keep vector
  keep.moh.diaz = c('muNq_diaz_mohamed_uncertainty', 'lc50.moh.diaz.report', 'se.lc50.moh.diaz', 'b1.moh.diaz')    

#Selecron (profenofos) reported LC50 and slope data #####################
  lc50.moh.prof.report = 4.02
  slp.moh.prof.report = 1.36
  b1.moh.prof = get_b1(slp.moh.prof.report)
  #SE not reported, so borrow from upper range of SEs reported from other studies
    se.lc50.moh.prof = 0.12
 
#Function to plot estimates without uncertainty
muNq_prof_mohamed = function(In, lc50 = lc50.moh.prof.report){
  Ins = (In/1000) #Parameters based on ppm, data input as ppb

    mun = pnorm(b1.moh.prof * log10(Ins/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }
     
#Create function based on reverse of litchfield and wilcoxon      
muNq_prof_mohamed_uncertainty = function(In){
  Ins = (In/1000) #Parameters based on ppm, data input as ppb
  lc50 = 10^(rnorm(1, log10(lc50.moh.prof.report), se.lc50.moh.prof)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
    mun = pnorm(b1.moh.prof * log10(Ins/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }
  
#keep vector
  keep.moh.prof = c('muNq_prof_mohamed_uncertainty', 'lc50.moh.prof.report', 'se.lc50.moh.prof', 'b1.moh.prof')    

#Reproduction data from tables 2A and B ######
#Diazinon  
#Concentrations tested  
  diaz_conc <- c(0, 1.40, 10.53, 12.25)

#Fecundity data over entire study period from table 2A    
  diaz_fecun <- c(3.50+1.35+2.55+2.70+18.18+14.50+4.30+4.60+8.70+11.20+2.40+3.50+3.57+5.50+5.40+8.20+9.30+11.20+10.60+12.80,
                  3.50+0.70+0.26+1.00+8.60+0+1.1+1.33+0.55+2.1+0+0+0+5+15.3+1+10.5,
                  3.50+1.45+0.75+0.46+17.10+14.7+11+1.25+21+19.3,
                  3.50+0+5.33+1.20)

#proportion of eggs that hatch in each concnetration group    
  diaz_hatch <- c(100,100,51.7,23.3)/100

#Fecundity rate as hatchlings/snail/week  
  diaz_rate <- diaz_fecun*diaz_hatch
  
#Reference value to scale to 0-1 
  moh_repro_ref <- diaz_fecun[1]

#Combine in a data frame  
  diaz_repro_df <- data.frame(conc = diaz_conc,
                              fecun = diaz_fecun,
                              hatch = diaz_hatch,
                              rate = diaz_rate)
  
#Estimate d-r function with drc: eggs/snail/week over entire study period as a function of diazinon concentration
  diaz_repro_mod= drm(rate ~ conc, data = diaz_repro_df, type = 'continuous',
                      fct = LL.3(names = c('b', 'd', 'e'),
                                 fixed = c(NA, max(diaz_repro_df$rate), NA)))

  fNq_moh_diaz_moh12<-function(In){
    predict(diaz_repro_mod, data.frame(conc = In/1000), interval = 'confidence', level = 0.95)
  }  
    
  par.tricks.diaz = c(coef(diaz_repro_mod), 'd' = max(diaz_repro_df$rate))[c(1,3,2)]
    
  fNq_moh_diaz_moh12_uncertainty<-function(In){
    fn <- rdrm(nosim = 1, fct = LL.3(), mpar = par.tricks.diaz, yerror = 'rnorm', xerror = In/1000,
               ypar = c(0, predict(diaz_repro_mod, data.frame(conc = In/1000), se.fit = T)[2]))$y / moh_repro_ref
    
    if(fn < 0){
      return(0)
    } else {
      return(fn)
    }
  }
  
keep.moh.diaz <- c(keep.moh.diaz, "fNq_moh_diaz_moh12_uncertainty", "diaz_repro_mod", "moh_repro_ref", "par.tricks.diaz")  

#Profenofos (selecron)  
#Concentrations tested  
  prof_conc <- c(0, 0.40, 2.44, 3.11)

#Fecundity data over entire study period from table 2A    
  prof_fecun <- c(diaz_fecun[1], #Same control for both chemcials
                  3.50+0.80+0.10+0.77+0.15+0.31+10.3+17.6+8+6,
                  3.50+0.75+0.72+1.80+14.70+6.4,
                  3.50+1.00)

#proportion of eggs that hatch in each concnetration group    
  prof_hatch <- c(100,100,0,0)/100

#Fecundity rate as hatchlings/snail/week  
  prof_rate <- prof_fecun*prof_hatch
  
#Combine in a data frame  
  prof_repro_df <- data.frame(conc = prof_conc,
                              fecun = prof_fecun,
                              hatch = prof_hatch,
                              rate = prof_rate)
  
#Estimate d-r function with drc: eggs/snail/week over entire study period as a function of profenofos concentration  
  prof_repro_mod= drm(rate ~ conc, data = prof_repro_df, type = 'continuous',
                      fct = LL.3(names = c('b', 'd', 'e'),
                                 fixed = c(NA, max(prof_repro_df$rate), NA)))

  fNq_moh_prof_moh12<-function(In){
    predict(prof_repro_mod, data.frame(conc = In/1000), interval = 'confidence', level = 0.95)
  }  
    
  par.tricks.prof = c(coef(prof_repro_mod), 'd' = max(prof_repro_df$rate))[c(1,3,2)]
    
  fNq_moh_prof_moh12_uncertainty<-function(In){
    fn <- rdrm(nosim = 1, fct = LL.3(), mpar = par.tricks.prof, yerror = 'rnorm', xerror = In/1000,
               ypar = c(0, predict(prof_repro_mod, data.frame(conc = In/1000), se.fit = T)[2]))$y / moh_repro_ref
    
    if(fn < 0){
      return(0)
    } else {
      return(fn)
    }
  }

keep.moh.prof <- c(keep.moh.prof, "fNq_moh_prof_moh12_uncertainty", "par.tricks.prof", "moh_repro_ref", "prof_repro_mod")

keep.moh12_all <- c(keep.moh.diaz, keep.moh.prof)

