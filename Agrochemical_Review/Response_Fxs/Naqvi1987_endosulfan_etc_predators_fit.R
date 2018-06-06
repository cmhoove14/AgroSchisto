#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Data extraction and model fitting to Naqvi et al 1987 P. clarkii toxicity
source("Agrochemical_Review/Models/litchfield_wilcoxon_get_b1_from_slope.R")

#Thiodan (endosulfan) #####################
  lc50.naq.endo.report = 423
  slp.naq.endo.report = 1.72
  b1.naq.endo = get_b1(slp.naq.endo.report)
  #get standard error from reported 95% CIs of lc50
    se.lc50.naq.endo = mean(c(log10(503.9/lc50.naq.endo.report), log10(lc50.naq.endo.report/356.3))) / 1.96
  
#Create function based on reverse of litchfield and wilcoxon      
muPq_endo_Naqvi87_uncertainty = function(In){
  lc50 = 10^(rnorm(1, log10(lc50.naq.endo.report), se.lc50.naq.endo)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
    mun = pnorm(b1.naq.endo * log10(In/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }
  
#keep vector
  keep.naq.endo = c('muPq_endo_Naqvi87_uncertainty', 'lc50.naq.endo.report', 'se.lc50.naq.endo', 'b1.naq.endo')    
  
#Treflan (trifluarlin) #############    
  lc50.naq.trif.report = 26
  slp.naq.trif.report = 3.24
  b1.naq.trif = get_b1(slp.naq.trif.report)
  #get standard error from reported 95% CIs of lc50
    se.lc50.naq.trif = mean(c(log10(28.9/lc50.naq.trif.report), log10(lc50.naq.trif.report/23.8))) / 1.96
  
#Create function based on reverse of litchfield and wilcoxon      
muPq_trif_Naqvi87_uncertainty = function(In){
  lc50 = 10^(rnorm(1, log10(lc50.naq.trif.report), se.lc50.naq.trif)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    #Estimate mortality, converting from ppb to ppm
    mun = pnorm(b1.naq.trif * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }
  
#keep vector
  keep.naq.trif = c('muPq_trif_Naqvi87_uncertainty', 'lc50.naq.trif.report', 'se.lc50.naq.trif', 'b1.naq.trif')    

#MSMA (Monosodium methanearsonate) #############    
  lc50.naq.msma.report = 1019
  slp.naq.msma.report = 2.39
  b1.naq.msma = get_b1(slp.naq.msma.report)
  #get standard error from reported 95% CIs of lc50
    se.lc50.naq.msma = mean(c(log10(1123.8/lc50.naq.msma.report), log10(lc50.naq.msma.report/916.8))) / 1.96
  
#Create function based on reverse of litchfield and wilcoxon      
muPq_msma_Naqvi87_uncertainty = function(In){
  lc50 = 10^(rnorm(1, log10(lc50.naq.msma.report), se.lc50.naq.msma)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    #Estimate mortality, converting from ppb to ppm
    mun = pnorm(b1.naq.msma * log10((In/1000)/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }
  
#keep vector
  keep.naq.msma = c('muPq_msma_Naqvi87_uncertainty', 'lc50.naq.msma.report', 'se.lc50.naq.msma', 'b1.naq.msma')    

#Oust (sulfometuron-methyl) ######
#Not found to be tosic even at extremely high conentrations to adults therefore function just returns 0  
  muPq_oust_Naqvi87_uncertainty = function(In){
    return(0)
  }
  