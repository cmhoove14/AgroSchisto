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
piM.hsh.ch = data.frame(conc = c(0.78, 1.8)*1000,
                        mort = c(.50 , .90),
                        surv = 0)
  piM.hsh.ch$surv = 1 - piM.hsh.ch$mort
  piM.hsh.ch$ppm = piM.hsh.ch$conc/1000
  piM.hsh.ch$log10 = log10(piM.hsh.ch$conc)
  piM.hsh.ch$probit = qnorm(piM.hsh.ch$mort, mean = 5)
  
plot(piM.hsh.ch$ppm, piM.hsh.ch$probit, pch = 16)
  lm.hsh.pim.ch = lm(probit ~ ppm, data = piM.hsh.ch)
  
  abline(coef(lm.hsh.pim.ch), lty = 2)
  
  lc50.hash.pim.ch = 0.78
  slp.hash.pim.ch = coef(lm.hsh.pim.ch)[2]
  se.lc50.hash.pim.ch = mean(c(log10(1.1 / lc50.hash.pim.ch), log10(lc50.hash.pim.ch / 0.56))) / 1.96

  plot(piM.hsh.ch$conc, piM.hsh.ch$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
       xlab = 'Chlorpyrifos (ppb)', ylab = '8-hr miracidial mortality',
       main = expression(paste(pi[M], ' estimate from Hasheesh 2011')))
    segments(x0 = 560, y0 = 0.5, x1 = 1100, y1 = 0.5)

fx.piM.hsh.chlor = function(In, lc = lc50.hash.pim.ch){
  Ins = In/1000
  pnorm(-slp.hash.pim.ch * (Ins - lc))
}    

  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piM.hsh.chlor), lty = 2, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piM.hsh.chlor, lc = 1.1), lty = 3, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piM.hsh.chlor, lc = 0.56), lty = 3, col = 2)

piM_ch_Hash11_uncertainty = function(In){
  if(In == 0) piM = 1 else{
    Ins = In/1000
    lc50 = 10^(rnorm(1, log10(lc50.hash.pim.ch), se.lc50.hash.pim.ch))
    piM = pnorm(-slp.hash.pim.ch * (Ins - lc50)) / fx.piM.hsh.chlor(0)
  }
  if(piM > 1) piM = 1
  if(piM < 0) piM = 0
  
  return(piM)
}
  points(seq(0,4000,10), sapply(seq(0,4000,10), piM_ch_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)

keep.hsh.ch.pim = c('piM_ch_Hash11_uncertainty', 'fx.piM.hsh.chlor', 'slp.hash.pim.ch',
                    'lc50.hash.pim.ch', 'se.lc50.hash.pim.ch')  
#Profenofos miracidia #########
piM.hsh.prof = data.frame(conc = c(1.5, 2.51)*1000,
                    mort = c(.50 , .90),
                    surv = 0)
  piM.hsh.prof$surv = 1 - piM.hsh.prof$mort
  piM.hsh.prof$surv = 1 - piM.hsh.prof$mort
  piM.hsh.prof$ppm = piM.hsh.prof$conc/1000
  piM.hsh.prof$log10 = log10(piM.hsh.prof$conc)
  piM.hsh.prof$probit = qnorm(piM.hsh.prof$mort, mean = 5)
  
  plot(piM.hsh.prof$ppm, piM.hsh.prof$probit, pch = 16)
    lm.hsh.pim.prof = lm(probit ~ ppm, data = piM.hsh.prof)
    
  abline(coef(lm.hsh.pim.prof), lty = 2)
  
  lc50.hash.pim.prof = 1.5
  slp.hash.pim.prof = coef(lm.hsh.pim.prof)[2]
  se.lc50.hash.pim.prof = mean(c(log10(1.95 / lc50.hash.pim.prof), log10(lc50.hash.pim.prof / 1.15))) / 1.96
  
plot(piM.hsh.prof$conc, piM.hsh.prof$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
     xlab = 'Profenofos (ppb)', ylab = '8-hr miracidial mortality',
     main = expression(paste(pi[M], ' estimate from Hasheesh 2011')))
  segments(x0 = 1150, y0 = 0.5, x1 = 1950, y1 = 0.5)

fx.piM.hsh.prof = function(In, lc = lc50.hash.pim.prof){
  Ins = In/1000
  pnorm(-slp.hash.pim.prof * (Ins - lc))
}  
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piM.hsh.prof), lty = 2, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piM.hsh.prof, lc = 1.95), lty = 3, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piM.hsh.prof, lc = 1.15), lty = 3, col = 2)
  
piM_pr_Hash11_uncertainty = function(In){
  if(In == 0) piM = 1 else{
    Ins = In/1000
    lc50 = 10^(rnorm(1, log10(lc50.hash.pim.prof), se.lc50.hash.pim.prof))
    piM = pnorm(-slp.hash.pim.prof * (Ins - lc50)) / fx.piM.hsh.prof(0)
  }
  if(piM > 1) piM = 1
  if(piM < 0) piM = 0
  
  return(piM)
}

  points(seq(0,4000,10), sapply(seq(0,4000,10), piM_pr_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)
  
keep.hsh.pim.prof = c('piM_pr_Hash11_uncertainty', 'slp.hash.pim.prof', 'lc50.hash.pim.prof',
                      'se.lc50.hash.pim.prof', 'fx.piM.hsh.prof')  
#miracidia keep vector #########
keep.hsh.pim = c(keep.hsh.ch.pim, keep.hsh.pim.prof)  

#profenofos cercariae #######
piC.hsh.prof = data.frame(conc = c(1.85, 2.85)*1000,
                          mort = c(.50 , .90),
                          surv = 0)
  piC.hsh.prof$surv = 1 - piC.hsh.prof$mort
  piC.hsh.prof$surv = 1 - piC.hsh.prof$mort
  piC.hsh.prof$ppm = piC.hsh.prof$conc/1000
  piC.hsh.prof$log10 = log10(piC.hsh.prof$conc)
  piC.hsh.prof$probit = qnorm(piC.hsh.prof$mort, mean = 5)

plot(piC.hsh.prof$ppm, piC.hsh.prof$probit, pch = 16)
  lm.hsh.pic.prof = lm(probit ~ ppm, data = piC.hsh.prof)
  
  abline(coef(lm.hsh.pic.prof), lty = 2)

  lc50.hash.pic.prof = 1.85
  slp.hash.pic.prof = coef(lm.hsh.pic.prof)[2]
  se.lc50.hash.pic.prof = mean(c(log10(2.59 / lc50.hash.pic.prof), log10(lc50.hash.pic.prof / 1.32))) / 1.96

plot(piC.hsh.prof$conc, piC.hsh.prof$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
     xlab = 'Profenofos (ppb)', ylab = '8-hr cercarial mortality',
     main = expression(paste(pi[C], ' estimate from Hasheesh 2011')))
  segments(x0 = 2590, y0 = 0.5, x1 = 1320, y1 = 0.5)

fx.piC.hsh.prof = function(In, lc = lc50.hash.pic.prof){
  Ins = In/1000
  pnorm(-slp.hash.pic.prof * (Ins - lc))
}  
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piC.hsh.prof), lty = 2, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piC.hsh.prof, lc = 2.59), lty = 3, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piC.hsh.prof, lc = 1.32), lty = 3, col = 2)

  piC_pr_Hash11_uncertainty = function(In){
    if(In == 0) piC = 1 else{
      Ins = In/1000
      lc50 = 10^(rnorm(1, log10(lc50.hash.pic.prof), se.lc50.hash.pic.prof))
      piC = pnorm(-slp.hash.pic.prof * (Ins - lc50)) / fx.piC.hsh.prof(0)
    }
    if(piC > 1) piC = 1
    if(piC < 0) piC = 0
    
    return(piC)
  }
  
    points(seq(0,4000,10), sapply(seq(0,4000,10), piC_pr_Hash11_uncertainty, simplify = T), 
           pch = 5, col = 4, cex = 0.5)

keep.hsh.pic.prof = c('piC_pr_Hash11_uncertainty', 'slp.hash.pic.prof', 'lc50.hash.pic.prof',
                      'se.lc50.hash.pic.prof', 'fx.piC.hsh.prof')  

#chlorpyrifos cercariae #######
piC.hsh.ch = data.frame(conc = c(0.96, 2.1)*1000,
                        mort = c(.50 , .90),
                        surv = 0)
  piC.hsh.ch$surv = 1 - piC.hsh.ch$mort
  piC.hsh.ch$ppm = piC.hsh.ch$conc/1000
  piC.hsh.ch$log10 = log10(piC.hsh.ch$conc)
  piC.hsh.ch$probit = qnorm(piC.hsh.ch$mort, mean = 5)
  
plot(piC.hsh.ch$ppm, piC.hsh.ch$probit, pch = 16)
  lm.hsh.pic.ch = lm(probit ~ ppm, data = piC.hsh.ch)
  
  abline(coef(lm.hsh.pic.ch), lty = 2)
  
  lc50.hash.pic.ch = 0.96
  slp.hash.pic.ch = coef(lm.hsh.pic.ch)[2]
  se.lc50.hash.pic.ch = mean(c(log10(1.44 / lc50.hash.pic.ch), log10(lc50.hash.pic.ch / 0.62))) / 1.96
  
plot(piC.hsh.ch$conc, piC.hsh.ch$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
     xlab = 'Chlorpyrifos (ppb)', ylab = '8-hr cercarial mortality',
     main = expression(paste(pi[C], ' estimate from Hasheesh 2011')))
  segments(x0 = 620, y0 = 0.5, x1 = 1440, y1 = 0.5)
  
  fx.piC.hsh.chlor = function(In, lc = lc50.hash.pic.ch){
    Ins = In/1000
    pnorm(-slp.hash.pic.ch * (Ins - lc))
  }    
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piC.hsh.chlor), lty = 2, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piC.hsh.chlor, lc = 1.44), lty = 3, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piC.hsh.chlor, lc = 0.62), lty = 3, col = 2)
  
  piC_ch_Hash11_uncertainty = function(In){
    if(In == 0) piC = 1 else{
      Ins = In/1000
      lc50 = 10^(rnorm(1, log10(lc50.hash.pic.ch), se.lc50.hash.pic.ch))
      piC = pnorm(-slp.hash.pic.ch * (Ins - lc50)) / fx.piC.hsh.chlor(0)
    }
    if(piC > 1) piC = 1
    if(piC < 0) piC = 0
    
    return(piC)
  }
  points(seq(0,4000,10), sapply(seq(0,4000,10), piC_ch_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)
  
  keep.hsh.ch.pic = c('piC_ch_Hash11_uncertainty', 'fx.piC.hsh.chlor', 'slp.hash.pic.ch',
                      'lc50.hash.pic.ch', 'se.lc50.hash.pic.ch')  
  
#cercariae keep vector #########
keep.hsh.pic = c(keep.hsh.ch.pic, keep.hsh.pic.prof)  