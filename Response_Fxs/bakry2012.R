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

#Snail (B. alexandrina) toxicity; create data frame, clean, etc. ##########
muN.bak = data.frame(atr = c(.330, 1.250, 4.750),
                     gly = c(.840, 3.150, 12.600),
                     mort = c(0.1, 0.5, 0.9))

  muN.bak$surv = 1 - muN.bak$mort

#direct snail toxicity from Atrazine ############  
plot(muN.bak$atr*1000, muN.bak$mort, ylim = c(0,1), xlim = c(0,5000),
       pch = 16, xlab = 'atrazine (ppb)', ylab = 'snail mortality',
     main = 'D-R function based on reported values')
  
    segments(y0 = 0.5, x0 = 830, y1 = 0.5, x1 = 1880)

#function based on reported values    
  lc50.atr.mun = 1.25
    se.lc50.atr.mun = mean(log10(1.88/1.25), log10(1.25/0.83)) / 1.96
  slp.atr.mun = 2.48
  
  fx.mun.atr = function(He, lc = lc50.atr.mun){
    heu = He/1000
    pnorm(slp.atr.mun * log10(heu/lc))
  }

    lines(c(0:5000), sapply(c(0:5000), fx.mun.atr), lty = 2, col = 2)
    lines(c(0:5000), sapply(c(0:5000), fx.mun.atr, lc = 1.88), lty = 3, col = 2)
    lines(c(0:5000), sapply(c(0:5000), fx.mun.atr, lc = 0.83), lty = 3, col = 2)
 
  muNq_atr_Bakry12_uncertainty = function(He){
    heu = He/1000
    lc50 = 10^(rnorm(1, log10(lc50.atr.mun), se.lc50.atr.mun))
    pnorm(slp.atr.mun * log10(heu/lc50))
  }   
    
    points(seq(0,5000,10), sapply(seq(0,5000,10), muNq_atr_Bakry12_uncertainty), 
           pch = 5, col = 4, cex = 0.5)  
    
#Function based on fit to LC values  
 plot(muN.bak$atr*1000, muN.bak$mort, ylim = c(0,1), xlim = c(0,5000),
     pch = 16, xlab = 'atrazine (ppb)', ylab = 'snail mortality',
     main = 'D-R function based on fit to LC values')    
  segments(y0 = 0.5, x0 = 830, y1 = 0.5, x1 = 1880)
 
  bak12.atr.drm = drm(mort ~ atr, weights = rep(50,3), data = muN.bak, type = 'binomial',
                      fct = LL.2())
  
    bak12.atr.pred = function(He){
      predict(bak12.atr.drm, newdata = data.frame(atr = He), interval = 'confidence', level = 0.95)
    }
    
  lines(seq(0,13, 0.1)*1000, sapply(seq(0,13, 0.1), bak12.atr.pred, simplify = T)[1,], col = 2, lty = 2)
  lines(seq(0,13, 0.1)*1000, sapply(seq(0,13, 0.1), bak12.atr.pred, simplify = T)[2,], col = 2, lty = 3)
  lines(seq(0,13, 0.1)*1000, sapply(seq(0,13, 0.1), bak12.atr.pred, simplify = T)[3,], col = 2, lty = 3)
  
  mu_Nq_atr_bak12_uncertainty<-function(He){
    He.u = He/1000
    rdrm(1, LL.2(), coef(bak12.atr.drm), He.u, yerror = 'rbinom', ypar = 50)$y / 50  
  }
  
  points(seq(0,13, 0.01)*1000, sapply(seq(0,13, 0.01)*1000, mu_Nq_atr_bak12_uncertainty, simplify = T),
         col = 4, pch = 5, cex = 0.5)
  
#direct snail toxicity from Glyphosate ############  
  plot(muN.bak$gly*1000, muN.bak$mort, ylim = c(0,1), xlim = c(0,13000),
       pch = 16, xlab = 'glyphosate (ppb)', ylab = 'snail mortality',
       main = 'D-R function based on reported values')  
  
      segments(y0 = 0.5, x0 = 890, y1 = 0.5, x1 = 4820)
#function based on provided values
  lc50.gly.mun = 3.15
    se.lc50.gly.mun = mean(log10(4.82/3.15), log10(3.15/0.89)) / 1.96
  slp.gly.mun = 2.16
  
  fx.mun.gly = function(He, lc = lc50.gly.mun){
    heu = He/1000
    pnorm(slp.gly.mun * log10(heu/lc))
  }
  
    lines(seq(0,13000, 13), sapply(seq(0,13000, 13), fx.mun.gly), lty = 2, col = 2)
    lines(seq(0,13000, 13), sapply(seq(0,13000, 13), fx.mun.gly, lc = 4.82), lty = 3, col = 2)
    lines(seq(0,13000, 13), sapply(seq(0,13000, 13), fx.mun.gly, lc = 0.89), lty = 3, col = 2)
  
  muNq_gly_Bakry12_uncertainty = function(He){
    heu = He/1000
    lc50 = 10^(rnorm(1, log10(lc50.gly.mun), se.lc50.gly.mun))
    pnorm(slp.gly.mun * log10(heu/lc50))
  }   
  
    points(seq(0,13000, 20), sapply(seq(0,13000, 20), muNq_gly_Bakry12_uncertainty), 
           pch = 5, col = 4, cex = 0.5)  

#function based on fit to LC values  
plot(muN.bak$gly*1000, muN.bak$mort, ylim = c(0,1), xlim = c(0,13000),
     pch = 16, xlab = 'glyphosate (ppb)', ylab = 'snail mortality',
     main = 'D-R function fit to LC values')  
    
    segments(y0 = 0.5, x0 = 890, y1 = 0.5, x1 = 4820)
    
  bak12.gly.drm = drm(mort ~ gly, weights = rep(50,3), data = muN.bak, type = 'binomial',
                      fct = LL.2())
  
    bak12.gly.pred = function(He){
      predict(bak12.gly.drm, newdata = data.frame(gly = He), interval = 'confidence', level = 0.95)
    }

    lines(seq(0,13, 0.1)*1000, sapply(seq(0,13, 0.1), bak12.gly.pred, simplify = T)[1,], col = 2, lty = 2)
    lines(seq(0,13, 0.1)*1000, sapply(seq(0,13, 0.1), bak12.gly.pred, simplify = T)[2,], col = 2, lty = 3)
    lines(seq(0,13, 0.1)*1000, sapply(seq(0,13, 0.1), bak12.gly.pred, simplify = T)[3,], col = 2, lty = 3)
    
    mu_Nq_gly_bak12_uncertainty<-function(He){
      He.u = He/1000
      rdrm(1, LL.2(), coef(bak12.gly.drm), He.u, yerror = 'rbinom', ypar = 50)$y / 50  
    }
    
    points(seq(0,13, 0.01)*1000, sapply(seq(0,13, 0.01)*1000, mu_Nq_gly_bak12_uncertainty, simplify = T),
           col = 4, pch = 5, cex = 0.5)
    
  
  keep.bak12 = c('muN.bak', 'mu_Nq_atr_bak12_uncertainty', 'mu_Nq_gly_bak12_uncertainty', 
                 'bak12.atr.drm', 'bak12.gly.drm')

#Now lets look at reproduction over time ###########
bakry12<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/bakry2012.csv')
  bak12 = subset(bakry12, time != 0)

#get estimates of hatchlings/snail/day across all weeks for each treatment   
fn.ctrl = mean(bak12$hatch_per[bak12$chem == 'control' & bak12$surv != 0] / 7)  
  fn.ctrl.sd = sd(bak12$hatch_per[bak12$chem == 'control' & bak12$surv != 0] / 7)
  
fn.atr = mean(bak12$hatch_per[bak12$chem == 'atrazine' & bak12$surv != 0] / 7)  
  fn.atr.sd = sd(bak12$hatch_per[bak12$chem == 'atrazine' & bak12$surv != 0] / 7)

fn.gly = mean(bak12$hatch_per[bak12$chem == 'glyphosate' & bak12$surv != 0] / 7)  
  fn.gly.sd = sd(bak12$hatch_per[bak12$chem == 'glyphosate' & bak12$surv != 0] / 7)
  

#assume lc90 halts all reproduction to provide third data point
bak.fn = data.frame(gly = c(0, .840, 12.600),
                    atr = c(0, .330, 4.750),
                    gly.r = c(fn.ctrl, fn.atr, 0),
                    atr.r = c(fn.ctrl, fn.gly, 0))  

#glyphosate influence on reproduction ############
fn.bak.gly = drm(gly.r ~ gly, data = bak.fn, type = 'continuous',
                 fct = L.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                           fixed = c(NA, 0, bak.fn$gly.r[1], NA)))

  summary(fn.bak.gly)

fn.bak.gly.pred = function(He){
  heu = He/1000
  predict(fn.bak.gly, newdata = data.frame(gly = heu), interval = 'confidence', level = 0.95)
}  

plot(bak.fn$gly*1000, bak.fn$gly.r / bak.fn$gly.r[1] , ylim = c(0,1), pch = 16,
     xlab = 'glyphosate (ppb)', ylab = 'relative hatchlings/snail/week')

  lines(seq(0,13000, 13), sapply(seq(0,13000, 13), fn.bak.gly.pred)[1,] / bak.fn$gly.r[1], col = 2, lty=2)
  lines(seq(0,13000, 13), sapply(seq(0,13000, 13), fn.bak.gly.pred)[2,] / bak.fn$gly.r[1], col = 2, lty=3)
  lines(seq(0,13000, 13), sapply(seq(0,13000, 13), fn.bak.gly.pred)[3,] / bak.fn$gly.r[1], col = 2, lty=3)

fN.gly.fx.uncertainty = function(He){
  if(He == 0) fn = 1 else{
    heu = He/1000
  init = predict(fn.bak.gly, newdata = data.frame(gly = heu), se.fit = T)
  fn = rnorm(1, init[1], init[2]) / bak.fn$gly.r[1]
  while(fn < 0 && fn > 1.00000){
    fn = rnorm(1, init[1], init[2]) / bak.fn$gly.r[1]
  }
}
  return(fn)
} #normalized to 1, upper limit at 1, lower limit at 0

  points(seq(0,13000, 13), sapply(seq(0,13000, 13), fN.gly.fx.uncertainty), pch = 5, cex = 0.5, col = 4)

#atrazine influence on reproduction ############
fn.bak.atr = drm(atr.r ~ atr, data = bak.fn, type = 'continuous',
                 fct = L.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                           fixed = c(NA, 0, bak.fn$atr.r[1], NA)))

summary(fn.bak.atr)

fn.bak.atr.pred = function(He){
  heu = He/1000
  predict(fn.bak.atr, newdata = data.frame(atr = heu), interval = 'confidence', level = 0.95)
}  

plot(bak.fn$atr*1000, bak.fn$atr.r / bak.fn$atr.r[1] , ylim = c(0,1), pch = 16,
     xlab = 'atrazine (ppb)', ylab = 'relative hatchlings/snail/week')

  lines(seq(0,5000,5), sapply(seq(0,5000,5), fn.bak.atr.pred)[1,] / bak.fn$atr.r[1], col = 2, lty=2)
  lines(seq(0,5000,5), sapply(seq(0,5000,5), fn.bak.atr.pred)[2,] / bak.fn$atr.r[1], col = 2, lty=3)
  lines(seq(0,5000,5), sapply(seq(0,5000,5), fn.bak.atr.pred)[3,] / bak.fn$atr.r[1], col = 2, lty=3)

fN.atr.fx.uncertainty = function(He){
  if(He == 0) fn = 1 else{
    heu = He/1000
  init = predict(fn.bak.atr, newdata = data.frame(atr = heu), se.fit = T)
  fn = rnorm(1, init[1], init[2]) / bak.fn$atr.r[1]
  while(fn < 0 && fn > 1.00000){
    fn = rnorm(1, init[1], init[2]) / bak.fn$atr.r[1]
  }
}
  return(fn)
} #normalized to 1, upper limit at 1, lower limit at 0

points(seq(0,5000, 10), sapply(seq(0,5000, 10), fN.atr.fx.uncertainty), pch = 5, cex = 0.5, col = 4)
