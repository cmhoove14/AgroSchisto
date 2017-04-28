#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Ghaffar 2016 LARVAL (S.mansoni) data
require(drc)

  but.vals = c(0, 556, 2417, 3906, 5560, 8703)
  gly.vals = c(0, 1506, 3875, 9174, 15062, 26249)
  pen.vals = c(0, 214.8, 535, 1299, 2148, 3762)
  
  L.3.fx = function(t, lc50 = lc50, slp = slp){
    1 / (1+exp(slp*(t - lc50)))
  }

#Cercarial (S. mansoni) toxicity ###############
  cer<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/ghaffar2016_cercariae.csv')

#butralin cercariae ###############
  plot(cer$time_hrs[cer$conc==0], cer$surv[cer$conc==0], 
       pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 1:length(unique(cer$conc[cer$chem == 'butralin']))){
      points(cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[i] & cer$chem == 'butralin'], 
             cer$surv[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[i] & cer$chem == 'butralin'], pch=17,
             col = i+1)
  }
  
#fit to control points *****************************************************************
gaf.cer.cont<-drm(cer$surv[cer$conc==0] ~ cer$time_hrs[cer$conc==0],
                  type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.cer.cont, newdata = data.frame(conc = seq(0,24,0.1))), lty=2)

  auc.cer.cont=integrate(f = L.3.fx, lc50 = coef(gaf.cer.cont)[2], slp = coef(gaf.cer.cont)[1],
                         lower=0, upper=24)[1]$value    

#556 ppb *********************************************************************************
gaf.cer.but556<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[1] & cer$chem == 'butralin'] ~ 
                    cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[1] & cer$chem == 'butralin'],
                    type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.but556, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 2)

auc.cer.but556=integrate(f = L.3.fx, lc50 = coef(gaf.cer.but556)[2], slp = coef(gaf.cer.but556)[1],
                       lower=0, upper=24)[1]$value   
    
#2417 ppb *********************************************************************************
gaf.cer.but2417<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[2] & cer$chem == 'butralin'] ~ 
                     cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[2] & cer$chem == 'butralin'],
                     type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.but2417, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 3)

auc.cer.but2417=integrate(f = L.3.fx, lc50 = coef(gaf.cer.but2417)[2], slp = coef(gaf.cer.but2417)[1],
                         lower=0, upper=24)[1]$value   

#3906 ppb ***********************************************************************************
gaf.cer.but3906<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[3] & cer$chem == 'butralin'] ~ 
                     cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[3] & cer$chem == 'butralin'],
                     type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.but3906, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 4)

auc.cer.but3906=integrate(f = L.3.fx, lc50 = coef(gaf.cer.but3906)[2], slp = coef(gaf.cer.but3906)[1],
                          lower=0, upper=24)[1]$value   

#5560 ppb **********************************************************************************
gaf.cer.but5560<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[4] & cer$chem == 'butralin'] ~ 
                     cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[4] & cer$chem == 'butralin'],
                     type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.but5560, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 5)

auc.cer.but5560=integrate(f = L.3.fx, lc50 = coef(gaf.cer.but5560)[2], slp = coef(gaf.cer.but5560)[1],
                          lower=0, upper=24)[1]$value   

#8703 ppb **********************************************************************************
gaf.cer.but8703<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[5] & cer$chem == 'butralin'] ~ 
                     cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[5] & cer$chem == 'butralin'],
                     type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.but8703, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 6)

auc.cer.but8703=integrate(f = L.3.fx, lc50 = coef(gaf.cer.but8703)[2], slp = coef(gaf.cer.but8703)[1],
                          lower=0, upper=24)[1]$value   

    title('Ghaffar2016 butralin toxicity to cercariae')
    legend('topright', legend = c('control', 556,2417,3906,5560,8703), pch = c(16,17,17,17,17,17), 
           col = c(1:6), cex=0.8, bty = 'n')
    
  gaf.but.aucs.cer = c(auc.cer.cont,auc.cer.but556,auc.cer.but2417,
                       auc.cer.but3906,auc.cer.but5560,auc.cer.but8703)/auc.cer.cont
    
    

#Compile L.2 data for functional responses #############
but.cer = data.frame(e = c(coef(gaf.cer.cont)[2], coef(gaf.cer.but556)[2], coef(gaf.cer.but2417)[2],
                           coef(gaf.cer.but3906)[2], coef(gaf.cer.but5560)[2], coef(gaf.cer.but8703)[2]),
                     e.se = c(summary(gaf.cer.cont)$coefficients[2,2], 
                              summary(gaf.cer.but556)$coefficients[2,2],
                              summary(gaf.cer.but2417)$coefficients[2,2], 
                              summary(gaf.cer.but3906)$coefficients[2,2],
                              summary(gaf.cer.but5560)$coefficients[2,2],
                              summary(gaf.cer.but8703)$coefficients[2,2]),
                     b = c(coef(gaf.cer.cont)[1], coef(gaf.cer.but556)[1], coef(gaf.cer.but2417)[1],
                           coef(gaf.cer.but3906)[1], coef(gaf.cer.but5560)[1], coef(gaf.cer.but8703)[1]),
                     b.se = c(summary(gaf.cer.cont)$coefficients[1,2], 
                              summary(gaf.cer.but556)$coefficients[1,2],
                              summary(gaf.cer.but2417)$coefficients[1,2], 
                              summary(gaf.cer.but3906)$coefficients[1,2],
                              summary(gaf.cer.but5560)$coefficients[1,2],
                              summary(gaf.cer.but8703)$coefficients[1,2]),
                     butralin = but.vals,
                     logbutralin = log(but.vals+1))

plot(but.cer$butralin, but.cer$e, pch = 16, xlab = 'butralin (ppb)', 
     ylab = 'L.2 Parameters (cercariae)', ylim = c(0,10))

  points(but.cer$butralin+200, but.cer$b, pch = 17, col=2)

for(i in 1:length(but.cer$butralin)){
  segments(x0 = but.cer$butralin[i], y0 = but.cer$e[i] + but.cer$e.se[i],
           x1 = but.cer$butralin[i], y1 = but.cer$e[i] - but.cer$e.se[i])
  segments(x0 = but.cer$butralin[i]+200, y0 = but.cer$b[i] + but.cer$b.se[i],
           x1 = but.cer$butralin[i]+200, y1 = but.cer$b[i] - but.cer$b.se[i], col=2)
} 

#functions fit to L.2 parameters ############
but.cer.lm.e = lm(e ~ butralin, weights = e.se^-1, data = but.cer) #linear response of LC50
  but.cer.pred = function(but){
    predict(but.cer.lm.e, newdata = data.frame(butralin = but), 
            interval = 'confidence', level = 0.95)
  }
    
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred, simplify = T)[3,], lty = 3)

but.cer.lm.e2 = lm(e ~ logbutralin, weights = e.se^-1, data = but.cer) #log-linear response of LC50
  but.cer.pred2 = function(but){
    predict(but.cer.lm.e2, newdata = data.frame(logbutralin = log(but+1)), 
            interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred2, simplify = T)[3,], lty = 3, col=3)

  AIC(but.cer.lm.e, but.cer.lm.e2)  #Linear is a better fit    

but.cer.lm.b = lm(b ~ butralin, weights = b.se^-1, data = but.cer)   
  but.cer.pred.b = function(but){
    predict(but.cer.lm.b, newdata = data.frame(butralin = but), interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred.b, simplify = T)[3,], lty = 3, col = 2)

legend('topleft', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
       legend = c('LC50 - Linear',
                  'LC50 - Exponential',
                  'Slp - linear',
                  '95% CI'), cex = 0.7)  
legend('top', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
       cex = 0.8, bty = 'n', title = 'Observed points')

but.piC.pred = function(He){
  e = as.numeric(predict(but.cer.lm.e, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
  b = as.numeric(predict(but.cer.lm.b, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  
  while(e.use < 0){
    e.use = rnorm(1, e[1], e[2])
  }      #resample if negative
  
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piC.ghaf_butr.lin_unc = function(He){
  if(He == 0) piC = 1
  else(piC = but.piC.pred(He) / but.piC.pred(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

but.piC.pred2 = function(He){
  e = as.numeric(predict(but.cer.lm.e2, newdata = data.frame(logbutralin = log(He+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(but.cer.lm.b, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  
  while(e.use < 0){
    e.use = rnorm(1, e[1], e[2])
  }      #resample if negative
  
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piC.ghaf_butr.exp_unc = function(He){
  if(He == 0) piC = 1
  else(piC = but.piC.pred2(He) / but.piC.pred2(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

#plot to test function ########
plot(but.cer$butralin, gaf.but.aucs.cer,
     xlab = 'Butralin (ppb)', ylab = 'relative AUC (cercariae-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0, 18000,50)/2, sapply(seq(0, 18000,50)/2, piC.ghaf_butr.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0, 18000,50)/2, sapply(seq(0, 18000,50)/2, piC.ghaf_butr.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)

#glyphosate cercariae ###############
plot(cer$time_hrs[cer$conc==0], cer$surv[cer$conc==0], 
      pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    for(i in 1:length(unique(cer$conc[cer$chem == 'glyphosate']))){
      points(cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[i] & cer$chem == 'glyphosate'], 
             cer$surv[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[i] & cer$chem == 'glyphosate'], pch=17,
             col = i+1)
    }
  
  lines(x=seq(0,24,0.1), y=predict(gaf.cer.cont, newdata = data.frame(conc = seq(0,24,0.1))), lty=2)
  
#1506 ppb *********************************************************************************
gaf.cer.gly1506<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[1] & cer$chem == 'glyphosate'] ~ 
                     cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[1] & cer$chem == 'glyphosate'],
                     type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.gly1506, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 2)

auc.cer.gly1506=integrate(f = L.3.fx, lc50 = coef(gaf.cer.gly1506)[2], slp = coef(gaf.cer.gly1506)[1],
                         lower=0, upper=24)[1]$value   

#3875 ppb *********************************************************************************
gaf.cer.gly3875<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[2] & cer$chem == 'glyphosate'] ~ 
                     cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[2] & cer$chem == 'glyphosate'],
                     type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.gly3875, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 3)

auc.cer.gly3875=integrate(f = L.3.fx, lc50 = coef(gaf.cer.gly3875)[2], slp = coef(gaf.cer.gly3875)[1],
                          lower=0, upper=24)[1]$value   

#9174 ppb ***********************************************************************************
gaf.cer.gly9174<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[3] & cer$chem == 'glyphosate'] ~ 
                     cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[3] & cer$chem == 'glyphosate'],
                     type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.gly9174, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 4)

auc.cer.gly9174=integrate(f = L.3.fx, lc50 = coef(gaf.cer.gly9174)[2], slp = coef(gaf.cer.gly9174)[1],
                          lower=0, upper=24)[1]$value   

#15062 ppb **********************************************************************************
gaf.cer.gly15062<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[4] & cer$chem == 'glyphosate'] ~ 
                      cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[4] & cer$chem == 'glyphosate'],
                      type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.gly15062, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 5)

auc.cer.gly15062=integrate(f = L.3.fx, lc50 = coef(gaf.cer.gly15062)[2], slp = coef(gaf.cer.gly15062)[1],
                          lower=0, upper=24)[1]$value   

#26249 ppb **********************************************************************************
gaf.cer.gly26249<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[5] & cer$chem == 'glyphosate'] ~ 
                      cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[5] & cer$chem == 'glyphosate'],
                      type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.gly26249, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 6)

auc.cer.gly26249=integrate(f = L.3.fx, lc50 = coef(gaf.cer.gly26249)[2], slp = coef(gaf.cer.gly26249)[1],
                           lower=0, upper=24)[1]$value   

    title('Ghaffar2016 glyphosate toxicity to cercariae')
    legend('topright', legend = c('control', 1506,3875,9174,15062,26249), pch = c(16,17,17,17,17,17), 
           col = c(1:6), cex=0.8, bty = 'n') 
    
    gaf.gly.aucs.cer = c(auc.cer.cont, auc.cer.gly1506, auc.cer.gly3875, 
                         auc.cer.gly9174, auc.cer.gly15062, auc.cer.gly26249)/auc.cer.cont

#Compile L.2 parameters ######## 
gly.cer = data.frame(e = c(coef(gaf.cer.cont)[2], coef(gaf.cer.gly1506)[2], coef(gaf.cer.gly3875)[2],
                           coef(gaf.cer.gly9174)[2], coef(gaf.cer.gly15062)[2], coef(gaf.cer.gly26249)[2]),
                     e.se = c(summary(gaf.cer.cont)$coefficients[2,2], 
                              summary(gaf.cer.gly1506)$coefficients[2,2],
                              summary(gaf.cer.gly3875)$coefficients[2,2], 
                              summary(gaf.cer.gly9174)$coefficients[2,2],
                              summary(gaf.cer.gly15062)$coefficients[2,2],
                              summary(gaf.cer.gly26249)$coefficients[2,2]),
                     b = c(coef(gaf.cer.cont)[1], coef(gaf.cer.gly1506)[1], coef(gaf.cer.gly3875)[1],
                           coef(gaf.cer.gly9174)[1], coef(gaf.cer.gly15062)[1], coef(gaf.cer.gly26249)[1]),
                     b.se = c(summary(gaf.cer.cont)$coefficients[1,2], 
                              summary(gaf.cer.gly1506)$coefficients[1,2],
                              summary(gaf.cer.gly3875)$coefficients[1,2], 
                              summary(gaf.cer.gly9174)$coefficients[1,2],
                              summary(gaf.cer.gly15062)$coefficients[1,2],
                              summary(gaf.cer.gly26249)$coefficients[1,2]),
                     glyphosate = gly.vals,
                     logglyphosate = log(gly.vals+1))

plot(gly.cer$glyphosate, gly.cer$e, pch = 16, xlab = 'glyphosate (ppb)', 
     ylab = 'L.2 Parameters (cercariae)', ylim = c(0,10))

  points(gly.cer$glyphosate+200, gly.cer$b, pch = 17, col=2)

for(i in 1:length(gly.cer$glyphosate)){
  segments(x0 = gly.cer$glyphosate[i], y0 = gly.cer$e[i] + gly.cer$e.se[i],
           x1 = gly.cer$glyphosate[i], y1 = gly.cer$e[i] - gly.cer$e.se[i])
  segments(x0 = gly.cer$glyphosate[i]+200, y0 = gly.cer$b[i] + gly.cer$b.se[i],
           x1 = gly.cer$glyphosate[i]+200, y1 = gly.cer$b[i] - gly.cer$b.se[i], col=2)
} 

#fit functions to LL.2 parameters ##########  
gly.cer.lm.e = lm(e ~ glyphosate, weights = e.se^-1, data = gly.cer) #linear response of LC50
  gly.cer.pred = function(gly){
    predict(gly.cer.lm.e, newdata = data.frame(glyphosate = gly), 
            interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred, simplify = T)[3,], lty = 3)

gly.cer.lm.e2 = lm(e ~ logglyphosate, weights = e.se^-1, data = gly.cer) #log-linear response of LC50
  gly.cer.pred2 = function(gly){
    predict(gly.cer.lm.e2, newdata = data.frame(logglyphosate = log(gly+1)), 
            interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred2, simplify = T)[3,], lty = 3, col=3)

  AIC(gly.cer.lm.e, gly.cer.lm.e2)  #Exponential is a better fit    

gly.cer.lm.b = lm(b ~ glyphosate, weights = b.se^-1, data = gly.cer)   
  gly.cer.pred.b = function(gly){
    predict(gly.cer.lm.b, newdata = data.frame(glyphosate = gly), interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

gly.piC.pred = function(He){
  e = as.numeric(predict(gly.cer.lm.e, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
  b = as.numeric(predict(gly.cer.lm.b, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  
  while(e.use < 0){
    e.use = rnorm(1, e[1], e[2])
  }      #resample if negative
  
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piC.ghaf_gly.lin_unc = function(He){
  if(He == 0) piC = 1
  else(piC = gly.piC.pred(He) / gly.piC.pred(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

gly.piC.pred2 = function(He){
  e = as.numeric(predict(gly.cer.lm.e2, newdata = data.frame(logglyphosate = log(He+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(gly.cer.lm.b, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  
  while(e.use < 0){
    e.use = rnorm(1, e[1], e[2])
  }      #resample if negative
  
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piC.ghaf_gly.exp_unc = function(He){
  if(He == 0) piC = 1
  else(piC = gly.piC.pred2(He) / gly.piC.pred2(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

#plot to test function ########
plot(gly.cer$glyphosate, gaf.gly.aucs.cer,
     xlab = 'Glyphosate (ppb)', ylab = 'relative AUC (cercariae-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0,9000,10)*3, sapply(seq(0,9000,10)*3, piC.ghaf_gly.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,9000,10)*3, sapply(seq(0,9000,10)*3, piC.ghaf_gly.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)



#pendimethalin cercariae ###############
plot(cer$time_hrs[cer$conc==0], cer$surv[cer$conc==0], 
      pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    for(i in 1:length(unique(cer$conc[cer$chem == 'pendimethalin']))){
      points(cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[i] & cer$chem == 'pendimethalin'], 
             cer$surv[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[i] & cer$chem == 'pendimethalin'], pch=17,
             col = i+1)
    }
  lines(x=seq(0,24,0.1), y=predict(gaf.cer.cont, newdata = data.frame(conc = seq(0,24,0.1))), lty=2)
  
#215 ppb *********************************************************************************
gaf.cer.pen215<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[1] & cer$chem == 'pendimethalin'] ~ 
                    cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[1] & cer$chem == 'pendimethalin'],
                    type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.pen215, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 2)

auc.cer.pen215=integrate(f = L.3.fx, lc50 = coef(gaf.cer.pen215)[2], slp = coef(gaf.cer.pen215)[1],
                          lower=0, upper=24)[1]$value   

#535 ppb *********************************************************************************
gaf.cer.pen535<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[2] & cer$chem == 'pendimethalin'] ~ 
                    cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[2] & cer$chem == 'pendimethalin'],
                    type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.pen535, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 3)

auc.cer.pen535=integrate(f = L.3.fx, lc50 = coef(gaf.cer.pen535)[2], slp = coef(gaf.cer.pen535)[1],
                         lower=0, upper=24)[1]$value   

#1299 ppb ***********************************************************************************
gaf.cer.pen1299<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[3] & cer$chem == 'pendimethalin'] ~ 
                     cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[3] & cer$chem == 'pendimethalin'],
                     type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.pen1299, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 4)

auc.cer.pen1299=integrate(f = L.3.fx, lc50 = coef(gaf.cer.pen1299)[2], slp = coef(gaf.cer.pen1299)[1],
                         lower=0, upper=24)[1]$value   

#2148 ppb **********************************************************************************
gaf.cer.pen2148<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[4] & cer$chem == 'pendimethalin'] ~ 
                     cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[4] & cer$chem == 'pendimethalin'],
                     type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.pen2148, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 5)

auc.cer.pen2148=integrate(f = L.3.fx, lc50 = coef(gaf.cer.pen2148)[2], slp = coef(gaf.cer.pen2148)[1],
                          lower=0, upper=24)[1]$value   

#3762 ppb **********************************************************************************
gaf.cer.pen3762<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[5] & cer$chem == 'pendimethalin'] ~ 
                     cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[5] & cer$chem == 'pendimethalin'],
                     type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
    
lines(x=seq(0,24,0.1), y=predict(gaf.cer.pen3762, newdata = data.frame(conc = seq(0,24,0.1))), 
      lty=2, col = 6)

auc.cer.pen3762=integrate(f = L.3.fx, lc50 = coef(gaf.cer.pen3762)[2], slp = coef(gaf.cer.pen3762)[1],
                          lower=0, upper=24)[1]$value   

    title('Ghaffar2016 pendimethalin toxicity to cercariae')
    legend('topright', legend = c('control', 215,535,1299,2148,3762), pch = c(16,17,17,17,17,17),
           col = c(1:6), cex=0.8, bty = 'n')  
    
    gaf.pen.aucs.cer = c(auc.cer.cont, auc.cer.pen215, auc.cer.pen535, auc.cer.pen1299, 
                         auc.cer.pen2148, auc.cer.pen3762)/auc.cer.cont

#Compile L.2 parameters ######## 
pen.cer = data.frame(e = c(coef(gaf.cer.cont)[2], coef(gaf.cer.pen215)[2], coef(gaf.cer.pen535)[2],
                           coef(gaf.cer.pen1299)[2], coef(gaf.cer.pen2148)[2], coef(gaf.cer.pen3762)[2]),
                     e.se = c(summary(gaf.cer.cont)$coefficients[2,2], 
                              summary(gaf.cer.pen215)$coefficients[2,2],
                              summary(gaf.cer.pen535)$coefficients[2,2], 
                              summary(gaf.cer.pen1299)$coefficients[2,2],
                              summary(gaf.cer.pen2148)$coefficients[2,2],
                              summary(gaf.cer.pen3762)$coefficients[2,2]),
                     b = c(coef(gaf.cer.cont)[1], coef(gaf.cer.pen215)[1], coef(gaf.cer.pen535)[1],
                           coef(gaf.cer.pen1299)[1], coef(gaf.cer.pen2148)[1], coef(gaf.cer.pen3762)[1]),
                     b.se = c(summary(gaf.cer.cont)$coefficients[1,2], 
                              summary(gaf.cer.pen215)$coefficients[1,2],
                              summary(gaf.cer.pen535)$coefficients[1,2], 
                              summary(gaf.cer.pen1299)$coefficients[1,2],
                              summary(gaf.cer.pen2148)$coefficients[1,2],
                              summary(gaf.cer.pen3762)$coefficients[1,2]),
                     pendimethalin = pen.vals,
                     logpendimethalin = log(pen.vals+1))

plot(pen.cer$pendimethalin, pen.cer$e, pch = 16, xlab = 'Pendimethalin (ppb)', 
     ylab = 'L.2 Parameters (cercariae)', ylim = c(0,12))

  points(pen.cer$pendimethalin+50, pen.cer$b, pch = 17, col=2)

  for(i in 1:length(pen.cer$pendimethalin)){
    segments(x0 = pen.cer$pendimethalin[i], y0 = pen.cer$e[i] + pen.cer$e.se[i],
             x1 = pen.cer$pendimethalin[i], y1 = pen.cer$e[i] - pen.cer$e.se[i])
    segments(x0 = pen.cer$pendimethalin[i]+50, y0 = pen.cer$b[i] + pen.cer$b.se[i],
             x1 = pen.cer$pendimethalin[i]+50, y1 = pen.cer$b[i] - pen.cer$b.se[i], col=2)
  } 

#fit functions to LL.2 parameters ##########  
pen.cer.lm.e = lm(e ~ pendimethalin, weights = e.se^-1, data = pen.cer) #linear response of LC50
  pen.cer.pred = function(pen){
    predict(pen.cer.lm.e, newdata = data.frame(pendimethalin = pen), 
            interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred, simplify = T)[3,], lty = 3)

pen.cer.lm.e2 = lm(e ~ logpendimethalin, weights = e.se^-1, data = pen.cer) #log-linear response of LC50
  pen.cer.pred2 = function(pen){
    predict(pen.cer.lm.e2, newdata = data.frame(logpendimethalin = log(pen+1)), 
            interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred2, simplify = T)[3,], lty = 3, col=3)
  
  AIC(pen.cer.lm.e, pen.cer.lm.e2)  #Exponential is a better fit    

pen.cer.lm.b = lm(b ~ pendimethalin, weights = b.se^-1, data = pen.cer)   
  pen.cer.pred.b = function(pen){
    predict(pen.cer.lm.b, newdata = data.frame(pendimethalin = pen), interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

pen.piC.pred = function(He){
  e = as.numeric(predict(pen.cer.lm.e, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
  b = as.numeric(predict(pen.cer.lm.b, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  
  while(e.use < 0){
    e.use = rnorm(1, e[1], e[2])
  }      #resample if negative
  
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piC.ghaf_pen.lin_unc = function(He){
  if(He == 0) piC = 1
  else(piC = pen.piC.pred(He) / pen.piC.pred(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

pen.piC.pred2 = function(He){
  e = as.numeric(predict(pen.cer.lm.e2, newdata = data.frame(logpendimethalin = log(He+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(pen.cer.lm.b, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  
  while(e.use < 0){
    e.use = rnorm(1, e[1], e[2])
  }      #resample if negative
  
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piC.ghaf_pen.exp_unc = function(He){
  if(He == 0) piC = 1
  else(piC = pen.piC.pred2(He) / pen.piC.pred2(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

#plot to test function ########
plot(pen.cer$pendimethalin, gaf.pen.aucs.cer,
     xlab = 'pendimethalin (ppb)', ylab = 'relative AUC (cercariae-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0,4000,25), sapply(seq(0,4000,25), piC.ghaf_pen.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,4000,25), sapply(seq(0,4000,25), piC.ghaf_pen.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)


