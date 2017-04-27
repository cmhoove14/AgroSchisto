#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#meta-analysis of atrazine effects on Bi. alexandrina

source('Response_Fxs/bakry2012.R')

atr.meta = data.frame(conc = c(muN.bak$atr*1000, 10100, 101620),
                     mort = c(muN.bak$mort, 0.1, 0.5),
                     tot = rep(50,5))

plot(atr.meta$conc, atr.meta$mort, pch = 16, ylim =c(0,1),
     xlab = 'atrazine (ppb)', ylab = 'snail mortality rate')

  atr.mun.meta = drm(mort ~ conc, weights = tot, data = atr.meta, type = 'binomial',
                     fct = LL.2())
  
    atr.meta.pred = function(He){
      predict(atr.mun.meta, newdata = data.frame(conc = He), interval = 'confidence', level = 0.95)
    }

  lines(seq(0,1.5e5,200), sapply(seq(0,1.5e5,200), atr.meta.pred, simplify = T)[1,], col = 2, lty = 2)
  lines(seq(0,1.5e5,200), sapply(seq(0,1.5e5,200), atr.meta.pred, simplify = T)[2,], col = 2, lty = 3)
  lines(seq(0,1.5e5,200), sapply(seq(0,1.5e5,200), atr.meta.pred, simplify = T)[3,], col = 2, lty = 3)

  mu_Nq_meta_atr_uncertainty<-function(He){
    rdrm(1, LL.2(), coef(atr.mun.meta), He, yerror = 'rbinom', ypar = 50)$y / 50  
  }

    points(seq(0,1.5e5,200), sapply(seq(0,1.5e5,200), mu_Nq_meta_atr_uncertainty, simplify = T),
           col = 4, pch = 5, cex = 0.5)
#Zoom
  plot(atr.meta$conc, atr.meta$mort, pch = 16, ylim =c(0,1), xlim = c(0,500),
       xlab = 'atrazine (ppb)', ylab = 'snail mortality rate')    
  
  lines(c(0:500), sapply(c(0:500), atr.meta.pred, simplify = T)[1,], col = 2, lty = 2)
  lines(c(0:500), sapply(c(0:500), atr.meta.pred, simplify = T)[2,], col = 2, lty = 3)
  lines(c(0:500), sapply(c(0:500), atr.meta.pred, simplify = T)[3,], col = 2, lty = 3)
  
#Compare to function from Bakry 2012 alone
  plot(atr.meta$conc, atr.meta$mort, pch = 16, ylim =c(0,1), xlim = c(0,12000),
       xlab = 'atrazine (ppb)', ylab = 'snail mortality rate')
  
    lines(seq(0,13, 0.1)*1000, sapply(seq(0,13, 0.1), bak12.atr.pred, simplify = T)[1,], col = 2, lty = 2)
    lines(seq(0,13, 0.1)*1000, sapply(seq(0,13, 0.1), bak12.atr.pred, simplify = T)[2,], col = 2, lty = 3)
    lines(seq(0,13, 0.1)*1000, sapply(seq(0,13, 0.1), bak12.atr.pred, simplify = T)[3,], col = 2, lty = 3)
  
    lines(seq(0,13, 0.1)*1000, sapply(seq(0,13, 0.1)*1000, atr.meta.pred, simplify = T)[1,], 
          col = 3, lty = 2)
    lines(seq(0,13, 0.1)*1000, sapply(seq(0,13, 0.1)*1000, atr.meta.pred, simplify = T)[2,], 
          col = 3, lty = 3)
    lines(seq(0,13, 0.1)*1000, sapply(seq(0,13, 0.1)*1000, atr.meta.pred, simplify = T)[3,], 
          col = 3, lty = 3)
    