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
