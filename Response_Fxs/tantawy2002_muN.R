#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
#Data extraction and model fitting to Tantawy 2002 data
require(drc)

#Snail toxicity ##########
#Butachlor ##########
  snail.but = data.frame(conc = c(0.65, 1.5, 4.5,6.5,44), # suspected typo in paper @4.5 (reads 45)
                         mort = c(0, .10, .25, .50, .90),
                         dead = c(0, 1, 2.5, 5, 9),
                         surv = 0)
  
  snail.but$surv = 1 - snail.but$mort 
  
  se.but = (log10(10.4) - log10(4.06)) / 1.96 #st. err of lc50 in ppm
  
  snail.but$se = (snail.but$conc / 1.32) * se.but #st. err proportional to concentration
  
  plot(snail.but$conc, snail.but$mort, pch = 16, ylim = c(0,1), xlim = c(0,50),
       xlab = 'Butachlor (ppm)', ylab = 'Snail mortality',
       main = expression(paste('Butachlor toxicity to ', italic('Bi. alexandrina'))))
  for(i in 1:length(unique(snail.but$conc))){
    segments(x0 = snail.but$conc[i] + snail.but$se[i], y0 = snail.but$mort[i],
             x1 = snail.but$conc[i] - snail.but$se[i], y1 = snail.but$mort[i])
  }
  
  but.sn = drm(dead/10 ~ conc, weights = rep(10, nrow(snail.but)), #assume cohort of 10 snails
               data = snail.but, type = 'binomial', 
               fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                          fixed = c(NA, 0, 1, NA)))  
    summary(but.sn)
  
  muN.tant.but = function(He){
    Her = He/1000
    predict(but.sn, data.frame(conc = Her), interval = 'confidence', level = 0.95)
  }  
  
  lines(c(0:50), sapply(c(0:50)*1000, muN.tant.but, simplify = T)[1,], lty = 2, col = 2)
  lines(c(0:50), sapply(c(0:50)*1000, muN.tant.but, simplify = T)[2,], lty = 3, col = 2)
  lines(c(0:50), sapply(c(0:50)*1000, muN.tant.but, simplify = T)[3,], lty = 3, col = 2)
  
  muN.tant.but_uncertainty<-function(He){
    Her = He/1000
    rdrm(1, LL.2(), coef(but.sn), Her, yerror = 'rbinom', ypar = 10)$y / 10  
  }
  
  points(c(0:50), sapply(c(0:50)*1000, muN.tant.but_uncertainty, simplify = T),
         pch = 5, cex = 0.6, col = 4)
#Fluazifop-p-butyl  ##########     
  snail.fpb = data.frame(conc = c(1.76, 4.5, 9,17.6,58), # suspected typo in paper @4.5 (reads 45)
                         mort = c(0, .10, .25, .50, .90),
                         dead = c(0, 1, 2.5, 5, 9),
                         surv = 0)
  
    snail.fpb$surv = 1 - snail.fpb$mort 
  
    se.fpb = (log10(11.73) - log10(26.4)) / 1.96 #st. err of lc50 in ppm
  
      snail.fpb$se = (snail.fpb$conc / 1.32) * se.fpb #st. err proportional to concentration
  
  plot(snail.fpb$conc, snail.fpb$mort, pch = 16, ylim = c(0,1), xlim = c(0,75),
       xlab = 'Fluzifop-p-butyl (ppm)', ylab = 'Snail mortality',
       main = expression(paste('Fluzifop-p-butyl toxicity to ', italic('Bi. alexandrina'))))
  for(i in 1:length(unique(snail.fpb$conc))){
    segments(x0 = snail.fpb$conc[i] + snail.fpb$se[i], y0 = snail.fpb$mort[i],
             x1 = snail.fpb$conc[i] - snail.fpb$se[i], y1 = snail.fpb$mort[i])
  }
  
  fpb.sn = drm(dead/10 ~ conc, weights = rep(10, nrow(snail.fpb)), #assume cohort of 10 snails
               data = snail.fpb, type = 'binomial', 
               fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                          fixed = c(NA, 0, 1, NA)))  
    summary(fpb.sn)
  
  muN.tant.fpb = function(He){
    Her = He/1000
    predict(fpb.sn, data.frame(conc = Her), interval = 'confidence', level = 0.95)
  }  
  
    lines(c(0:75), sapply(c(0:75)*1000, muN.tant.fpb, simplify = T)[1,], lty = 2, col = 2)
    lines(c(0:75), sapply(c(0:75)*1000, muN.tant.fpb, simplify = T)[2,], lty = 3, col = 2)
    lines(c(0:75), sapply(c(0:75)*1000, muN.tant.fpb, simplify = T)[3,], lty = 3, col = 2)
    
  muN.tant.fpb_uncertainty<-function(He){
    Her = He/1000
    rdrm(1, LL.2(), coef(fpb.sn), Her, yerror = 'rbinom', ypar = 10)$y / 10  
  }
  
  points(c(0:75), sapply(c(0:75)*1000, muN.tant.fpb_uncertainty, simplify = T),
         pch = 5, cex = 0.6, col = 4)    
#keep vector #########  
keep.muN.tantawy = c('but.sn', 'fpb.sn', 'muN.tant.but', 'muN.tant.but_uncertainty', 
                     'muN.tant.fpb', 'muN.tant.fpb_uncertainty', 'snail.fpb')  