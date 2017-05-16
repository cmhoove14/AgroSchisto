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
mun.bak.atr = data.frame(conc = c(.330, 1.250, 4.750)*1000,
                         mort = c(0.1, 0.5, 0.9)) 
  mun.bak.atr$ppm = mun.bak.atr$conc/1000
  mun.bak.atr$log10 = log10(mun.bak.atr$ppm)
  mun.bak.atr$probit = qnorm(mun.bak.atr$mort, mean = 5)

plot(mun.bak.atr$log10, mun.bak.atr$probit, pch = 16)    
  lm.bak.atr = lm(probit ~ log10, data = mun.bak.atr)
  
  abline(coef(lm.bak.atr), lty = 2)
  
  lc50.bak.atr = 10^((5 - coef(lm.bak.atr)[1]) / coef(lm.bak.atr)[2])
  slp.bak.atr = coef(lm.bak.atr)[2]
  #get standard error from reported 95% CIs of lc50
  se.lc50.bak.atr = mean(c(log10(1.88 / lc50.bak.atr), log10(lc50.bak.atr / 0.83))) / 1.96
  
plot(mun.bak.atr$conc, mun.bak.atr$mort, ylim = c(0,1), xlim = c(0,5000),
       pch = 16, xlab = 'atrazine (ppb)', ylab = 'snail mortality',
     main = 'D-R function based on reported values')
    segments(y0 = 0.5, x0 = 830, y1 = 0.5, x1 = 1880)
    
  fx.bak.atr = function(He, lc = lc50.bak.atr){
    heu = (He/1000)
    pnorm(slp.bak.atr * log10(heu/(lc)))
  }  
    lines(seq(0,5000,10), sapply(seq(0,5000,10), fx.bak.atr), lty = 2, col = 2)
    lines(seq(0,5000,10), sapply(seq(0,5000,10), fx.bak.atr, lc = 0.83), lty = 2, col = 2)
    lines(seq(0,5000,10), sapply(seq(0,5000,10), fx.bak.atr, lc = 1.88), lty = 2, col = 2)
    
muNq_atr_Bakry12_uncertainty = function(He){
  if(He == 0) mun = 0 else{
    heu = (He/1000)
    lc50 = 10^(rnorm(1, log10(lc50.bak.atr), se.lc50.bak.atr))
  while(lc50 < 0){
    lc50 = 10^(rnorm(1, log10(lc50.bak.atr), se.lc50.bak.atr))
  }  
    mun = pnorm((slp.bak.atr) * log10(heu/lc50)) - fx.bak.atr(0)
  }
  while(mun < 0){
    lc50 = 10^(rnorm(1, log10(lc50.bak.atr), se.lc50.bak.atr))
    mun = pnorm((slp.bak.atr) * log10(heu/lc50)) - fx.bak.atr(0)
  } 
  return(mun)
}
points(seq(0,5000,10), sapply(seq(0,5000,10), muNq_atr_Bakry12_uncertainty), 
       pch = 5, col = 4, cex = 0.5)

#keep vector
keep.bak.atr = c('muNq_atr_Bakry12_uncertainty', 'fx.bak.atr',
                 'lc50.bak.atr', 'se.lc50.bak.atr', 'slp.bak.atr')    

#direct snail toxicity from Glyphosate ############  
mun.bak.gly = data.frame(conc = c(.840, 3.150, 12.600)*1000,
                         mort = c(0.1, 0.5, 0.9))
  
  mun.bak.gly$ppm = mun.bak.gly$conc/1000
  mun.bak.gly$log10 = log10(mun.bak.gly$ppm)
  mun.bak.gly$probit = qnorm(mun.bak.gly$mort, mean = 5)

plot(mun.bak.gly$log10, mun.bak.gly$probit, pch = 16)
  lm.bak.gly = lm(probit ~ log10, data = mun.bak.gly)
  
  abline(coef(lm.bak.gly), lty = 2)
  
  lc50.bak.gly = 10^((5 - coef(lm.bak.gly)[1]) / coef(lm.bak.gly)[2])
  slp.bak.gly = coef(lm.bak.gly)[2]
  #get standard error from reported 95% CIs of lc50
  se.lc50.bak.gly = mean(c(log10(4.82 / lc50.bak.gly), log10(lc50.bak.gly / 0.89))) / 1.96
  
      
  plot(mun.bak.gly$conc, mun.bak.gly$mort, ylim = c(0,1), xlim = c(0,13000),
       pch = 16, xlab = 'glyphosate (ppb)', ylab = 'snail mortality',
       main = 'D-R function based on reported values')  
  
      segments(y0 = 0.5, x0 = 890, y1 = 0.5, x1 = 4820)
      
fx.bak.gly = function(He, lc = lc50.bak.gly){
  heu = (He/1000)
  pnorm(slp.bak.gly * log10(heu/(lc)))
}  
  lines(seq(0,13000,13), sapply(seq(0,13000,13), fx.bak.gly), lty = 2, col = 2)
  lines(seq(0,13000,13), sapply(seq(0,13000,13), fx.bak.gly, lc = 0.89), lty = 2, col = 2)
  lines(seq(0,13000,13), sapply(seq(0,13000,13), fx.bak.gly, lc = 4.82), lty = 2, col = 2)

muNq_gly_Bakry12_uncertainty = function(He){
  if(He == 0) mun = 0 else{
    heu = (He/1000)
    lc50 = 10^(rnorm(1, log10(lc50.bak.gly), se.lc50.bak.gly))
  while(lc50 < 0){
    lc50 = 10^(rnorm(1, log10(lc50.bak.gly), se.lc50.bak.gly))
  }  
    mun = pnorm((slp.bak.gly) * log10(heu/lc50)) - fx.bak.gly(0)
  }
  while(mun < 0){
    lc50 = 10^(rnorm(1, log10(lc50.bak.gly), se.lc50.bak.gly))
    mun = pnorm((slp.bak.gly) * log10(heu/lc50)) - fx.bak.gly(0)
  } 
  return(mun)
}
points(seq(0,13000,25), sapply(seq(0,13000,25), muNq_gly_Bakry12_uncertainty), 
       pch = 5, col = 4, cex = 0.5)

#keep vector
keep.bak.gly = c('muNq_gly_Bakry12_uncertainty', 'fx.bak.gly',
                 'lc50.bak.gly', 'se.lc50.bak.gly', 'slp.bak.gly')    

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

  lines(seq(0,13000, 25), sapply(seq(0,13000, 25), fn.bak.gly.pred)[1,] / bak.fn$gly.r[1], col = 2, lty=2)
  lines(seq(0,13000, 25), sapply(seq(0,13000, 25), fn.bak.gly.pred)[2,] / bak.fn$gly.r[1], col = 2, lty=3)
  lines(seq(0,13000, 25), sapply(seq(0,13000, 25), fn.bak.gly.pred)[3,] / bak.fn$gly.r[1], col = 2, lty=3)

fN.gly.fx.uncertainty = function(He){
  if(He == 0) fn = 1 else{
    heu = He/1000
  init = predict(fn.bak.gly, newdata = data.frame(gly = heu), se.fit = T)
  fn = rnorm(1, init[1], init[2]) / bak.fn$gly.r[1]
  while(fn < 0 || fn > 1.00000){
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

  lines(seq(0,5000,25), sapply(seq(0,5000,25), fn.bak.atr.pred)[1,] / bak.fn$atr.r[1], col = 2, lty=2)
  lines(seq(0,5000,25), sapply(seq(0,5000,25), fn.bak.atr.pred)[2,] / bak.fn$atr.r[1], col = 2, lty=3)
  lines(seq(0,5000,25), sapply(seq(0,5000,25), fn.bak.atr.pred)[3,] / bak.fn$atr.r[1], col = 2, lty=3)

fN.atr.fx.uncertainty = function(He){
  if(He == 0) fn = 1 else{
    heu = He/1000
  init = predict(fn.bak.atr, newdata = data.frame(atr = heu), se.fit = T)
  fn = rnorm(1, init[1], init[2]) / bak.fn$atr.r[1]
  while(fn < 0 || fn > 1.00000){
    fn = rnorm(1, init[1], init[2]) / bak.fn$atr.r[1]
  }
}
  return(fn)
} #normalized to 1, upper limit at 1, lower limit at 0

points(seq(0,5000, 10), sapply(seq(0,5000, 10), fN.atr.fx.uncertainty), pch = 5, cex = 0.5, col = 4)
