#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
source("Agrochemical_Review/Models/litchfield_wilcoxon_get_b1_from_slope.R")
   
#Toxicity to miracidia from table 5 ###############
#ChlorP miracidia #########
  lc50.hash.pim.ch.report = 0.78
  slp.hash.pim.ch.report = 1.86
  b1.hash.pim.ch.report = get_b1(slp.hash.pim.ch.report)
  
  #get standard error from reported 95% CIs of lc50
    se.lc50.hash.pim.ch = mean(c(log10(1.1 / lc50.hash.pim.ch.report), 
                                 log10(lc50.hash.pim.ch.report / 0.56))) / 1.96
  
#Create function based on reverse of litchfield and wilcoxon      
piM_ch_Hash11_uncertainty = function(In){
  Ins = (In/1000) #Parameters based on ppm, data input as ppb
  lc50 = 10^(rnorm(1, log10(lc50.hash.pim.ch.report), se.lc50.hash.pim.ch)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
    pim = pnorm(-b1.hash.pim.ch.report * log10(Ins/lc50)) #Estimate daily mortality (percent)

    return(pim)
  }
  
#keep vector
  keep.hash.pim.ch = c('piM_ch_Hash11_uncertainty', 'lc50.hash.pim.ch.report', 'se.lc50.hash.pim.ch', 'b1.hash.pim.ch.report')  
  
#Profenofos miracidia #########
  lc50.hash.pim.prof.report = 1.5
  slp.hash.pim.prof.report = 1.64
  b1.hash.pim.prof.report = get_b1(slp.hash.pim.prof.report)
  #Get standard error form reported 95%CI
    se.lc50.hash.pim.prof = mean(c(log10(1.95 / lc50.hash.pim.prof.report), 
                                   log10(lc50.hash.pim.prof.report / 1.15))) / 1.96
  
  piM_pr_Hash11_uncertainty = function(In){
    Ins = In/1000
    lc50 = 10^(rnorm(1, log10(lc50.hash.pim.prof.report), se.lc50.hash.pim.prof))
      piM = pnorm(-b1.hash.pim.prof.report * log10(Ins / lc50))

    return(piM)
  }

#keep vector    
  keep.hash.pim.prof = c('piM_pr_Hash11_uncertainty', 'b1.hash.pim.prof.report', 'lc50.hash.pim.prof.report', 'se.lc50.hash.pim.prof')  
  
#miracidia keep vector #########
keep.hash.pim = c(keep.hash.pim.ch, keep.hash.pim.prof)  

#Toxicity to cercariae from table 5 ###############
#chlorpyrifos cercariae #######
  lc50.hash.pic.ch.report = 0.96
  slp.hash.pic.ch.report = 2.3
  b1.hash.pic.ch.report = get_b1(slp.hash.pic.ch.report)
    se.lc50.hash.pic.ch = mean(c(log10(1.44 / lc50.hash.pic.ch.report), 
                                 log10(lc50.hash.pic.ch.report / 0.62))) / 1.96

  piC_ch_Hash11_uncertainty = function(In){
    Ins = In/1000
    lc50 = 10^(rnorm(1, log10(lc50.hash.pic.ch.report), se.lc50.hash.pic.ch))
    piC = pnorm(-b1.hash.pic.ch.report * log10(Ins / lc50))

    return(piC)
  }
     
  keep.hash.ch.pic = c('piC_ch_Hash11_uncertainty', 'b1.hash.pic.ch.report',
                       'lc50.hash.pic.ch.report', 'se.lc50.hash.pic.ch')  
  

#profenofos cercariae #######
  lc50.hash.pic.prof.report = 1.85
  slp.hash.pic.prof.report = 1.84
  b1.hash.pic.prof.report = get_b1(slp.hash.pic.prof.report)
    se.lc50.hash.pic.prof = mean(c(log10(2.59 / lc50.hash.pic.prof.report), 
                                   log10(lc50.hash.pic.prof.report / 1.32))) / 1.96
  
  piC_pr_Hash11_uncertainty = function(In){
    Ins = In/1000
    lc50 = 10^(rnorm(1, log10(lc50.hash.pic.prof.report), se.lc50.hash.pic.prof))
    piC = pnorm(-b1.hash.pic.prof.report * log10(Ins / lc50))

    return(piC)
  }
  
keep.hash.pic.prof = c('piC_pr_Hash11_uncertainty', 'b1.hash.pic.prof.report', 'lc50.hash.pic.prof.report',
                       'se.lc50.hash.pic.prof')
  
#cercariae keep vector #########
keep.hash.pic = c(keep.hash.ch.pic, keep.hash.pic.prof)  