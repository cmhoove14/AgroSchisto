#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#direct snail toxicity from Atrazine ############  
lc50.bak.atr.report = 1.25
slp.bak.atr.report = 2.48

#Get standard error of reported LC50 based on 95%CI
  se.lc50.bak.atr = mean(c(log10(1.88 / lc50.bak.atr.report), log10(lc50.bak.atr.report / 0.83))) / 1.96
  
#Create function based on reverse of litchfield and wilcoxon      
  muNq_atr_Bakry12_uncertainty = function(He){
    Heu = (He/1000) #Parameters based on ppm, data input as ppb
    lc50 = 10^(rnorm(1, log10(lc50.bak.atr.report), se.lc50.bak.atr)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
    mun = pnorm((slp.bak.atr.report) * log10(Heu/lc50)) #Estimate daily mortality (percent)
    
    return(mun)
  }

#keep vector
keep.bak.atr = c('muNq_atr_Bakry12_uncertainty', 'lc50.bak.atr.report', 'se.lc50.bak.atr', 'slp.bak.atr.report')    

#direct snail toxicity from Glyphosate ############  
lc50.bak.gly.report = 3.15
slp.bak.gly.report = 2.16

#Get standard error of reported LC50 based on 95%CI
  se.lc50.bak.gly = log10(4.82 / lc50.bak.gly.report) / 1.96 
  #Lower limit excluded from this estimate as upper limit leads to more sensible estimate of SE, lower limit it way out of range, assuming a manuscript typo or misprint

#Create function based on reverse of litchfield and wilcoxon      
muNq_gly_Bakry12_uncertainty = function(He){
  Heu = (He/1000) #Parameters based on ppm, data input as ppb
  lc50 = 10^(rnorm(1, log10(lc50.bak.gly.report), se.lc50.bak.gly)) #Estimate lc50 with uncertainty and backtransform from log10 scale
  
  mun = pnorm((slp.bak.gly.report) * log10(Heu/lc50)) #Estimate daily mortality (percent)
  
  return(mun)
}

#keep vector
keep.bak.gly = c('muNq_gly_Bakry12_uncertainty', 'lc50.bak.gly.report', 'se.lc50.bak.gly', 'slp.bak.gly.report')    

#Snail reproduction over time, not analyzed at this point because only control and single dose group for each chemical ######
bakry12<-read.csv('Agrochemical_Review/Response_Fxs/Data/bakry2012.csv')

keep.bak.all = c(keep.bak.atr, keep.bak.gly)