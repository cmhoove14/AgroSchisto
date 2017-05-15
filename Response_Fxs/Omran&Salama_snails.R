#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#herbicide toxicity to Bi. alexandrina from Omran and Salama 2013 ###############

#Load packages, clean/supplement data, fit initial model #############
require(drc)

ons = read.csv('~/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/Omran&Salama2013.csv')
  ons$mort = pnorm(ons$probit, mean = 5)
  ons$surv = 1 - ons$mort

ons.atr = subset(ons, chem == 'atrazine')
ons.gly = subset(ons, chem == 'glyphosate')

#Atrazine mortality ############
plot(ons.atr$log_conc, ons.atr$probit, pch = 16, col = 'gold', #ylim = c(0,7), xlim = c(0, 3),
     xlab = 'atrazine (ppm)', ylab = c('Probit mortality')) 
  lm.ons.atr = lm(probit ~ log_conc, data = ons.atr)
  
  abline(coef(lm.ons.atr), lty = 2, col = 'gold')
  
  lc50.ons.atr = 10^((5 - coef(lm.ons.atr)[1]) / coef(lm.ons.atr)[2])
  slp.ons.atr = coef(lm.ons.atr)[2]

plot(ons.atr$conc, ons.atr$mort, pch = 16, ylim = c(0,1),
     xlab = 'atrazine (ppb)', ylab = c('mortality')) 
    
fx.ons.atr = function(He){
  heu = log10(He/1000)
  pnorm(predict(lm.ons.atr, newdata = data.frame(log_conc = heu), 
                interval = 'confidence', level = 0.95), mean = 5)
} 
  lines(seq(0,5e5,500), sapply(seq(0,5e5,500), fx.ons.atr)[1,], lty = 2, col = 2)
  lines(seq(0,5e5,500), sapply(seq(0,5e5,500), fx.ons.atr)[2,], lty = 3, col = 2)
  lines(seq(0,5e5,500), sapply(seq(0,5e5,500), fx.ons.atr)[3,], lty = 3, col = 2)

ons.munq.atr = function(He){
  if(He == 0) mun = 0 else{
    heu = log10(He/1000)
    init = predict(lm.ons.atr, newdata = data.frame(log_conc = heu), se.fit = T)
    est = rnorm(1, init$fit, init$se.fit)
  while(est < 0){
    est = rnorm(1, init$fit, init$se.fit)
  }
    mun = pnorm(est, mean = 5)
  }
  return(mun)
} 

points(seq(0,5e5,1000), sapply(seq(0,5e5,1000), ons.munq.atr),
       pch = 5, col = 4, cex = 0.5)

keep.ons.atr = c('ons.munq.atr', 'lm.ons.atr')
#zoom just to check out low conc
plot(ons.atr$conc, ons.atr$mort, pch = 16, ylim = c(0,1), xlim = c(0,10000),
     xlab = 'atrazine (ppb)', ylab = c('mortality')) 
  lines(seq(0,10000,20), sapply(seq(0,10000,20), fx.ons.atr)[1,], lty = 2, col = 2)
  lines(seq(0,10000,20), sapply(seq(0,10000,20), fx.ons.atr)[2,], lty = 3, col = 2)
  lines(seq(0,10000,20), sapply(seq(0,10000,20), fx.ons.atr)[3,], lty = 3, col = 2)

  points(seq(0,10000,50), sapply(seq(0,10000,50), ons.munq.atr),
         pch = 5, col = 4, cex = 0.5)

#glyphosate toxicity ############   
plot(ons.gly$log_conc, ons.gly$probit, pch = 16, #xlim = c(0,3), ylim = c(4,6)
     col = 3)  
  lm.ons.gly = lm(probit ~ log_conc, data = ons.gly)  

  abline(coef(lm.ons.gly), lty = 2, col = 3)
  
  lc50.ons.gly = 10^((5 - coef(lm.ons.gly)[1]) / coef(lm.ons.gly)[2])
  slp.ons.gly = coef(lm.ons.gly)[2]
  
fx.ons.gly = function(He){
  heu = log10(He/1000)
  pnorm(predict(lm.ons.gly, newdata = data.frame(log_conc = heu), 
          interval = 'confidence', level = 0.95), mean = 5)
} 

plot(ons.gly$conc, ons.gly$mort, pch = 16, ylim = c(0,1),
     xlab = 'Glyphosate (ppb)', ylab = 'mortality')
  lines(seq(0,5e5,500), sapply(seq(0,5e5,500), fx.ons.gly)[1,], lty = 2, col = 2)
  lines(seq(0,5e5,500), sapply(seq(0,5e5,500), fx.ons.gly)[2,], lty = 3, col = 2)
  lines(seq(0,5e5,500), sapply(seq(0,5e5,500), fx.ons.gly)[3,], lty = 3, col = 2)
  
ons.munq.gly = function(He){
  if(He == 0) mun = 0 else{
    heu = log10(He/1000)
    init = predict(lm.ons.gly, newdata = data.frame(log_conc = heu), se.fit = T)
    est = rnorm(1, init$fit, init$se.fit)
    while(est < 0){
      est = rnorm(1, init$fit, init$se.fit)
    }
    mun = pnorm(est, mean = 5)
  }
  return(mun)
}

points(seq(0,5e5,1000), sapply(seq(0,5e5,1000), ons.munq.gly),
       pch = 5, col = 4, cex = 0.5)

keep.ons.gly = c('ons.munq.gly', 'lm.ons.gly')

#zoom to more environmentally feasible range (still 5000 ppb which is real high)
plot(ons.gly$conc, ons.gly$mort, pch = 16, ylim = c(0,1), xlim = c(0,10000),
     xlab = 'Glyphosate (ppb)', ylab = 'mortality')
  lines(seq(0,10000,100), sapply(seq(0,10000,100), fx.ons.gly)[1,], lty = 2, col = 2)
  lines(seq(0,10000,100), sapply(seq(0,10000,100), fx.ons.gly)[2,], lty = 3, col = 2)
  lines(seq(0,10000,100), sapply(seq(0,10000,100), fx.ons.gly)[3,], lty = 3, col = 2)

  points(seq(0,10000,50), sapply(seq(0,10000,50), ons.munq.gly),
         pch = 5, col = 4, cex = 0.5)
  