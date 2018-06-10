#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Ragab 2006 data
source("Agrochemical_Review/Models/litchfield_wilcoxon_get_b1_from_slope.R")

#Ammonium Nitrate Snail toxicity ##########
lc50.amm.rag = 470
slp.amm.rag = get_b1(1.30)
se.lc50.amm.rag = mean(c(log10(507.6 / lc50.amm.rag), log10(lc50.amm.rag / 435.19))) / 1.96 

lc90.amm.rag = 640
  
rag06_muNq_amm = function(Fe, lc = lc50.amm.rag){
  feu = Fe/1000
  lc50 = 10^(rnorm(1, log10(lc50.amm.rag), se.lc50.amm.rag)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
  mun = pnorm(slp.amm.rag * log10(feu/lc50)) #Estimate daily mortality (percent)

  return(mun)
} 

#Pottasium Sulfate Snail toxicity ##########    
lc50.pot.rag = 1900
slp.pot.rag = get_b1(1.27)
se.lc50.pot.rag = mean(c(log10(2280 / lc50.pot.rag), 
                         log10(lc50.pot.rag / 1583.3))) / 1.96 
lc90.pot.rag = 2600

rag06_muNq_pot = function(Fe){
    feu = Fe/1000
    lc50 = 10^(rnorm(1, log10(lc50.pot.rag), se.lc50.pot.rag))
    mun = pnorm(slp.pot.rag * log10(feu / lc50))
  
  return(mun)
}

#Urea Snail toxicity ##########    
lc50.urea.rag = 22000
slp.urea.rag = get_b1(1.28)
se.lc50.urea.rag = mean(c(log10(24860 / lc50.urea.rag), 
                          log10(lc50.urea.rag / 19469))) / 1.96 
lc90.urea.rag = 31000

rag06_muNq_urea = function(Fe){
    feu = Fe/1000
    lc50 = 10^(rnorm(1, log10(lc50.urea.rag), se.lc50.urea.rag))
    mun = pnorm((slp.urea.rag) * log10(feu / lc50))
  
  return(mun)
}

#Keep vector #################
keep.ragab.mun = c('rag06_muNq_urea', 'lc50.urea.rag', 'se.lc50.urea.rag', 'slp.urea.rag',
                   'rag06_muNq_pot', 'lc50.pot.rag', 'se.lc50.pot.rag', 'slp.pot.rag',
                   'rag06_muNq_amm', 'lc50.amm.rag', 'se.lc50.amm.rag', 'slp.amm.rag')
