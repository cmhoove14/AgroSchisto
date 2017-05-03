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
  ons$mort = pnorm(ons$prob_mort, mean = 5, sd = 1)
  ons$surv = 1 - ons$mort
  ons$dead = round(50*ons$mort)
  ons$live = round(50*ons$surv)
  ons$total = 50
  
ons.atr = subset(ons, chem == 'atrazine')
ons.gly = subset(ons, chem == 'glyphosate')

plot(ons.atr$log_conc, ons.atr$prob_mort, pch = 16, col = 'gold', ylim = c(0,7), xlim = c(0, 3),
     xlab = 'Log10 Concentration', ylab = c('Probit mortality')) 
  points(ons.gly$log_conc, ons.gly$prob_mort, pch = 16, col = 'green')

ons.lm.atr = lm(prob_mort ~ log_conc, data = ons.atr)  
  summary(ons.lm.atr)
  
ons.pred.atr = function(He){
  predict(ons.lm.atr, newdata = data.frame(log_conc = He), interval = 'confidence', level = 0.95)
} 
  lines(seq(0,4,0.01), sapply(seq(0,4,0.01), ons.pred.atr)[1,], lty = 2, col = 'gold')
  lines(seq(0,4,0.01), sapply(seq(0,4,0.01), ons.pred.atr)[2,], lty = 3, col = 'gold')
  lines(seq(0,4,0.01), sapply(seq(0,4,0.01), ons.pred.atr)[3,], lty = 3, col = 'gold')
  
ons.lm.gly = lm(prob_mort ~ log_conc, data = ons.gly)  
  summary(ons.lm.gly)
  
ons.pred.gly = function(He){
  predict(ons.lm.gly, newdata = data.frame(log_conc = He), interval = 'confidence', level = 0.95)
} 
  lines(seq(0,4,0.01), sapply(seq(0,4,0.01), ons.pred.gly)[1,], lty = 2, col = 'green')
  lines(seq(0,4,0.01), sapply(seq(0,4,0.01), ons.pred.gly)[2,], lty = 3, col = 'green')
  lines(seq(0,4,0.01), sapply(seq(0,4,0.01), ons.pred.gly)[3,], lty = 3, col = 'green')
  
title(expression(paste('Omran & Salama herbicide toxicity to ', italic('Bi. alexandrina'))))
legend('bottomright', pch = 16, col = c('gold', 'green'), legend = c('Atrazine', 'Glyphosate'), bty = 'n', cex = 0.8)

#derive function for Atrazine  ############
plot(ons.atr$conc10_ppm, ons.atr$mort, pch = 16,  ylim = c(0,1), xlim = c(0,500),
     xlab = 'atrazine (ppm)', ylab = 'mortality rate')

ons.munq.atr = function(He){
  if(He == 0) munq = 0 else{
    heu = log10(He/1000)  #convert ppb to log10 of ppm
    init = predict(ons.lm.atr, newdata = data.frame(log_conc = heu), se.fit = T) 
    munq = pnorm(rnorm(1, init$fit, init$se.fit), mean = 5, sd = 1)
  }
  return(munq)
}

points(c(0:500), sapply(c(0:500)*1000, ons.munq.atr), pch = 5, cex = 0.5, col = 4)

#zoom to more environmentally feasible range (still 5000 ppb which is real high)
plot(ons.atr$conc10_ppm, ons.atr$mort, pch = 16,  ylim = c(0,1), xlim = c(0,5),
     xlab = 'atrazine (ppm)', ylab = 'mortality rate')

points(seq(0,5,0.01), sapply(seq(0,5,0.01)*1000, ons.munq.atr), pch = 5, cex = 0.5, col = 4)

#derive function for Glyphosate  ############
plot(ons.gly$conc10_ppm, ons.gly$mort, pch = 16,  ylim = c(0,1), xlim = c(0,500),
     xlab = 'glyphosate (ppm)', ylab = 'mortality rate')

ons.munq.gly = function(He){
  if(He == 0) munq = 0 else{
    heu = log10(He/1000)
    init = predict(ons.lm.gly, newdata = data.frame(log_conc = heu), se.fit = T) 
    munq = pnorm(rnorm(1, init$fit, init$se.fit), mean = 5, sd = 1)
  }
  return(munq)
}

points(c(0:500), sapply(c(0:500)*1000, ons.munq.gly), pch = 5, cex = 0.5, col = 4)


#zoom to more environmentally feasible range (still 5000 ppb which is real high)
plot(ons.gly$conc10_ppm, ons.gly$mort, pch = 16,  ylim = c(0,1), xlim = c(0,5),
     xlab = 'glyphosate (ppm)', ylab = 'mortality rate')

points(seq(0,5,0.01), sapply(seq(0,5,0.01)*1000, ons.munq.gly), pch = 5, cex = 0.5, col = 4)