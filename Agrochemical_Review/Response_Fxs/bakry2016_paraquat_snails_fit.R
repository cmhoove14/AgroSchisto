#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Bakry 2012 data
source("Agrochemical_Review/Response_Fxs/bakry2012_atrazine_glyphosate_snails_fit.R")

#direct snail toxicity from paraquat ############  
  lc50.bak.prq.report = 2.3
  slp.bak.prq.report = 1.5

  #NO STANDARD ERROR REPORTED SO BORROW UNCERTAINTY FROM BAKRY 2012 ATRAZINE STUDY
  se.lc50.bak.prq = se.lc50.bak.atr

  muNq_prq_Bakry16_uncertainty = function(He){
    heu = (He/1000)
    lc50 = 10^(rnorm(1, log10(lc50.bak.prq.report), se.lc50.bak.prq))
      mun = pnorm((slp.bak.prq.report) * log10(heu/lc50))
   
    return(mun)
  }
  
#keep vector
keep.bak.prq = c('muNq_prq_Bakry16_uncertainty', 'lc50.bak.prq.report', 'se.lc50.bak.prq', 'slp.bak.prq.report')    
