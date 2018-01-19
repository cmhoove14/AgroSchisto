#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
    
#Toxicity to miracidia and cercariae from table 5 ###############
#ChlorP miracidia #########
piM.hsh.ch = data.frame(conc = c(0, 0.78, 1.8)*1000,
                        mort = c(0, .50 , .90),
                        surv = 0)
  piM.hsh.ch$surv = 1 - piM.hsh.ch$mort
  piM.hsh.ch$ppm = piM.hsh.ch$conc/1000
  piM.hsh.ch$log10 = log10(piM.hsh.ch$conc)
  piM.hsh.ch$probit = qnorm(piM.hsh.ch$mort, mean = 5)
  
  lc50.hash.pim.ch.report = 0.78
  slp.hash.pim.ch.report = 1.86
  se.lc50.hash.pim.ch = mean(c(log(1.1 / lc50.hash.pim.ch.report), log(lc50.hash.pim.ch.report / 0.56))) / 1.96
  
  piM_ch_Hash11_uncertainty = function(In){
    #if(In == 0) piM = 1 else{
    Ins = In/1000
    lc50 = exp(rnorm(1, log(lc50.hash.pim.ch.report), se.lc50.hash.pim.ch))
    piM = pnorm(-slp.hash.pim.ch.report * log(Ins / lc50))
    #}
   return(piM)
  }

   plot(piM.hsh.ch$conc, piM.hsh.ch$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
       xlab = 'Chlorpyrifos (ppb)', ylab = '8-hr miracidial mortality',
       main = expression(paste(pi[M], ' estimate from Hasheesh 2011')))
    segments(x0 = 560, y0 = 0.5, x1 = 1100, y1 = 0.5)
    
   points(seq(0,4000,10), sapply(seq(0,4000,10), piM_ch_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)
 
keep.hsh.ch.pim = c('piM_ch_Hash11_uncertainty', 'slp.hash.pim.ch.report',
                    'lc50.hash.pim.ch.report', 'se.lc50.hash.pim.ch')  

#Profenofos miracidia #########
piM.hsh.prof = data.frame(conc = c(1.5, 2.51)*1000,
                    mort = c(.50 , .90),
                    surv = 0)
  piM.hsh.prof$surv = 1 - piM.hsh.prof$mort
  piM.hsh.prof$surv = 1 - piM.hsh.prof$mort
  piM.hsh.prof$ppm = piM.hsh.prof$conc/1000
  piM.hsh.prof$log10 = log10(piM.hsh.prof$conc)
  piM.hsh.prof$probit = qnorm(piM.hsh.prof$mort, mean = 5)
  
  lc50.hash.pim.prof.report = 1.5
  slp.hash.pim.prof.report = 1.64
  se.lc50.hash.pim.prof = mean(c(log(1.95 / lc50.hash.pim.prof.report), log(lc50.hash.pim.prof.report / 1.15))) / 1.96
  
  piM_pr_Hash11_uncertainty = function(In){
    #if(In == 0) piM = 1 else{
    Ins = In/1000
    lc50 = exp(rnorm(1, log(lc50.hash.pim.prof.report), se.lc50.hash.pim.prof))
    piM = pnorm(-slp.hash.pim.prof.report * log(Ins / lc50))
    #}
    
    return(piM)
  }
  
  plot(piM.hsh.prof$conc, piM.hsh.prof$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
       xlab = 'Profenofos (ppb)', ylab = '8-hr miracidial mortality',
       main = expression(paste(pi[M], ' estimate from Hasheesh 2011')))
    segments(x0 = 1150, y0 = 0.5, x1 = 1950, y1 = 0.5)
  
  points(seq(0,4000,10), sapply(seq(0,4000,10), piM_pr_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)
  
  keep.hsh.pim.prof = c('piM_pr_Hash11_uncertainty', 'slp.hash.pim.prof.report', 'lc50.hash.pim.prof.report',
                      'se.lc50.hash.pim.prof')  
  
#miracidia keep vector #########
keep.hsh.pim = c(keep.hsh.ch.pim, keep.hsh.pim.prof)  

#profenofos cercariae #######
piC.hsh.prof = data.frame(conc = c(0,1.85, 2.85)*1000,
                          mort = c(0,.50 , .90),
                          surv = 0)
  piC.hsh.prof$surv = 1 - piC.hsh.prof$mort
  piC.hsh.prof$surv = 1 - piC.hsh.prof$mort
  piC.hsh.prof$ppm = piC.hsh.prof$conc/1000
  piC.hsh.prof$log10 = log10(piC.hsh.prof$conc)
  piC.hsh.prof$probit = qnorm(piC.hsh.prof$mort, mean = 5)

  lc50.hash.pic.prof.report = 1.85
  slp.hash.pic.prof.report = 1.84
  se.lc50.hash.pic.prof = mean(c(log(2.59 / lc50.hash.pic.prof.report), log(lc50.hash.pic.prof.report / 1.32))) / 1.96
  
  piC_pr_Hash11_uncertainty = function(In){
    # if(In == 0) piC = 1 else{
    Ins = In/1000
    lc50 = exp(rnorm(1, log(lc50.hash.pic.prof.report), se.lc50.hash.pic.prof))
    piC = pnorm(-slp.hash.pic.prof.report * log(Ins / lc50))
    # }
    
    return(piC)
  }
  
  plot(piC.hsh.prof$conc, piC.hsh.prof$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
       xlab = 'Profenofos (ppb)', ylab = '8-hr cercarial mortality',
       main = expression(paste(pi[C], ' estimate from Hasheesh 2011')))
  segments(x0 = 2590, y0 = 0.5, x1 = 1320, y1 = 0.5)
  
 points(seq(0,4000,10), sapply(seq(0,4000,10), piC_pr_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)  

keep.hsh.pic.prof = c('piC_pr_Hash11_uncertainty', 'slp.hash.pic.prof.report', 'lc50.hash.pic.prof.report',
                      'se.lc50.hash.pic.prof')  

#chlorpyrifos cercariae #######
piC.hsh.ch = data.frame(conc = c(0,0.96, 2.1)*1000,
                        mort = c(0,.50 , .90),
                        surv = 0)
  piC.hsh.ch$surv = 1 - piC.hsh.ch$mort
  piC.hsh.ch$ppm = piC.hsh.ch$conc/1000
  piC.hsh.ch$log10 = log10(piC.hsh.ch$conc)
  piC.hsh.ch$probit = qnorm(piC.hsh.ch$mort, mean = 5)
  
  lc50.hash.pic.ch.report = 0.96
  slp.hash.pic.ch.report = 2.3
  se.lc50.hash.pic.ch = mean(c(log(1.44 / lc50.hash.pic.ch.report), log(lc50.hash.pic.ch.report / 0.62))) / 1.96

     piC_ch_Hash11_uncertainty = function(In){
    #if(In == 0) piC = 1 else{
      Ins = In/1000
      lc50 = exp(rnorm(1, log(lc50.hash.pic.ch.report), se.lc50.hash.pic.ch))
      piC = pnorm(-slp.hash.pic.ch.report * log(Ins / lc50))
    #}
    
    return(piC)
  }
  
plot(piC.hsh.ch$conc, piC.hsh.ch$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
     xlab = 'Chlorpyrifos (ppb)', ylab = '8-hr cercarial mortality',
     main = expression(paste(pi[C], ' estimate from Hasheesh 2011')))
  segments(x0 = 620, y0 = 0.5, x1 = 1440, y1 = 0.5)
  
  points(seq(0,4000,10), sapply(seq(0,4000,10), piC_ch_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)

  
  keep.hsh.ch.pic = c('piC_ch_Hash11_uncertainty', 'slp.hash.pic.ch.report',
                      'lc50.hash.pic.ch.report', 'se.lc50.hash.pic.ch')  
  
#cercariae keep vector #########
keep.hsh.pic = c(keep.hsh.ch.pic, keep.hsh.pic.prof)  