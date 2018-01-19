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
require(drc)

#direct snail toxicity from Atrazine ############  
mun.bak.prq = data.frame(conc = c(0, 1.1, 2.3, 4.4)*1000,
                         mort = c(0, 0.25, 0.5, 0.9)) 
  mun.bak.prq$ppm = mun.bak.prq$conc/1000
  mun.bak.prq$log10 = log10(mun.bak.prq$ppm)
  mun.bak.prq$probit = qnorm(mun.bak.prq$mort, mean = 5)

  lc50.bak.prq.report = 2.3
  slp.bak.prq.report = 1.5
  #se.lc50.bak.prq = mean(c(log10(1.88 / lc50.bak.prq.report), log10(lc50.bak.prq.report / 0.83))) / 1.96
  
  
  plot(mun.bak.prq$conc, mun.bak.prq$mort, ylim = c(0,1), xlim = c(0,5000),
       pch = 16, xlab = 'atrazine (ppb)', ylab = 'snail mortality',
       main = 'D-R function based on reported values')
    #segments(y0 = 0.5, x0 = 830, y1 = 0.5, x1 = 1880)
  
  muNq_atr_Bakry12_uncertainty = function(He){
    #if(He == 0) mun = 0 else{
    heu = (He/1000)
    lc50 = 10^(rnorm(1, log10(lc50.bak.prq.report), se.lc50.bak.prq))
      mun = pnorm((slp.bak.prq.report) * log10(heu/lc50))
    # }
    # if(mun < 0) mun = 0
    
    return(mun)
  }
  
  points(seq(0,5000,10), sapply(seq(0,5000,10), muNq_atr_Bakry12_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
  
#keep vector
keep.bak.prq = c('muNq_atr_Bakry12_uncertainty', 
                 'lc50.bak.prq.report', 'se.lc50.bak.prq', 'slp.bak.prq.report')    
