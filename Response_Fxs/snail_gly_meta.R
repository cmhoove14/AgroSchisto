#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#meta-analysis of glyphosate effects on Bi. alexandrina

source('Response_Fxs/Ghaffar2016_snails.R')
source('Response_Fxs/bakry2012.R')

#mortality effects ########
#combine data from each study (last two data points from Omran & Salama 2013)
gly.mun = data.frame(conc = c(gly.dat$glyphosate, muN.bak$gly*1000, 4200, 41600),
                     mort = c(gly.dat$lcs/100, muN.bak$mort, 0.1, 0.5),
                     tot = c(rep(30,6), rep(50,5)))
plot(gly.mun$conc, gly.mun$mort, pch = 16, ylim = c(0,1), 
     xlab = 'glyphosate (ppb)', ylab = 'mortality')
  
  gly.mun.meta = drm(mort ~ conc, weights = tot, data = gly.mun, type = 'binomial',
                     fct = LL.2())
  gly.meta.pred = function(He){
    predict(gly.mun.meta, newdata = data.frame(conc = He), interval = 'confidence', level = 0.95)
  }
  
  lines(seq(0,45000,100), sapply(seq(0,45000,100), gly.meta.pred, simplify = T)[1,], col = 2, lty = 2)
  lines(seq(0,45000,100), sapply(seq(0,45000,100), gly.meta.pred, simplify = T)[2,], col = 2, lty = 3)
  lines(seq(0,45000,100), sapply(seq(0,45000,100), gly.meta.pred, simplify = T)[3,], col = 2, lty = 3)
  
  mu_Nq_meta_gly_uncertainty<-function(He){
    rdrm(1, LL.2(), coef(gly.mun.meta), He, yerror = 'rbinom', ypar = 50)$y / 50  
  }
  
  points(seq(0,45000,100), sapply(seq(0,45000,100), mu_Nq_meta_gly_uncertainty, simplify = T),
         col = 4, pch = 5, cex = 0.5)
  
#compare meta estimate to individual estimates ########
plot(gly.dat$glyphosate, gly.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,46000),
     xlab = 'Glyphosate (ppb)', ylab = 'snail mortality rate')
  
  points(muN.bak$gly*1000, muN.bak$mort, pch = 17)
  points(c(4200, 41600), c(0.1, 0.5), pch = 15)
  legend('bottomright', legend = c('Abdel-Ghaffar 2016', 'Bakry 2012', 'Omran&Salama 2013'),
         pch = c(16,17,15), cex = 0.7, bty = 'n', title = 'Data')
  
  lines(seq(0,46000,100), sapply(seq(0,46000,100), gly.meta.pred, simplify = T)[1,], 
        lty = 2, lwd = 2)
  lines(seq(0,46000,100), sapply(seq(0,46000,100), mu_Nq_gly_gaf16, simplify = T)[1,],
        lty = 2, col = 2, lwd = 2) 
  lines(seq(0,46000,100), sapply(seq(0,46000,100)/1000, bak12.gly.pred, simplify = T)[1,], 
        col = 3, lty = 2, lwd = 2)
  
  legend('bottom', legend = c('Meta','Abdel-Ghaffar 2016', 'Bakry 2012'),
         lty = 2, lwd = 2, col = c(1,2,3), cex = 0.7, bty = 'n', title = 'Fit D-R fxs')
  
#Reproductive effects ############
plot(gafrep$gly.conc, gafrep$gly.rep/gafrep$but.rep[1], pch = 16, ylim = c(0,1),
     xlab = 'glyphosate (ppb)', ylab = 'relative reproduction rate')
  points(bak.fn$gly.conc, bak.fn$gly.rep/bak.fn$gly.rep[1], pch = 17)
  
  lines(seq(0,10000,50), sapply(seq(0,10000,50), gly.r0.pred, simplify = T)[1,]/gafrep$but.rep[1],
        lty = 2, col = 2) 
  
#No attempt to fit function as raw measurements of reproduction are not directly comparable, but
  #added data point from Bakry12 seems to reinforce the lower limit from the Abdel-Ghaffar '16 study
  #and fit generally with the estimated function