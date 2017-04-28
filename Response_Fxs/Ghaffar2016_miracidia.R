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

#Miracidial (S. mansoni) toxicity ###############
mir<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Miracidial Mortality/ghaffar2016_miracidia.csv')

#butralin ###############
plot(mir$time_hrs[mir$conc==0], mir$surv[mir$conc==0], 
     pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24),
     main = 'Abdel-Ghaffar 2016: Butralin toxicity to miracidia')
  for(i in 1:length(unique(mir$conc[mir$chem == 'butralin']))){
    points(mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[i] & mir$chem == 'butralin'], 
           mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[i] & mir$chem == 'butralin'], pch=17,
           col = i+1)
  }

#fit to control points *****************************************************************
  gaf.mir.cont<-drm(mir$surv[mir$conc==0] ~ mir$time_hrs[mir$conc==0],
                    type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.cont, newdata = data.frame(conc = seq(0,24,0.1))), lty=2)
  
  auc.mir.cont=integrate(f = L.3.fx, lc50 = coef(gaf.mir.cont)[2], slp = coef(gaf.mir.cont)[1],
                         lower=0, upper=24)[1]$value  

#556 ppb *********************************************************************************
  gaf.mir.but556<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[1] & mir$chem == 'butralin'] ~ 
                        mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[1] & mir$chem == 'butralin'],
                      type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.but556, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 2)
  
  auc.mir.but556=integrate(f = L.3.fx, lc50 = coef(gaf.mir.but556)[2], slp = coef(gaf.mir.but556)[1],
                         lower=0, upper=24)[1]$value  
  
#2417 ppb *********************************************************************************
  gaf.mir.but2417<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[2] & mir$chem == 'butralin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[2] & mir$chem == 'butralin'],
                       type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.but2417, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 3)
  
  auc.mir.but2417=integrate(f = L.3.fx, lc50 = coef(gaf.mir.but2417)[2], slp = coef(gaf.mir.but2417)[1],
                           lower=0, upper=24)[1]$value    

#3906 ppb ***********************************************************************************
  gaf.mir.but3906<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[3] & mir$chem == 'butralin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[3] & mir$chem == 'butralin'],
                       type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.but3906, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 4)
  
  auc.mir.but3906=integrate(f = L.3.fx, lc50 = coef(gaf.mir.but3906)[2], slp = coef(gaf.mir.but3906)[1],
                            lower=0, upper=24)[1]$value  

#5560 ppb **********************************************************************************
  gaf.mir.but5560<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[4] & mir$chem == 'butralin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[4] & mir$chem == 'butralin'],
                       type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.but5560, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 5)
  
  auc.mir.but5560=integrate(f = L.3.fx, lc50 = coef(gaf.mir.but5560)[2], slp = coef(gaf.mir.but5560)[1],
                            lower=0, upper=24)[1]$value      

#8703 ppb **********************************************************************************
  gaf.mir.but8703<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[5] & mir$chem == 'butralin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[5] & mir$chem == 'butralin'],
                       type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.but8703, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 6)
  
  auc.mir.but8703=integrate(f = L.3.fx, lc50 = coef(gaf.mir.but8703)[2], slp = coef(gaf.mir.but8703)[1],
                            lower=0, upper=24)[1]$value      
  
  legend('topright', legend = c('control', 556,2417,3906,5560,8703), pch = c(16,17,17,17,17,17), 
         col = c(1:6), cex=0.8, bty = 'n')
  
  gaf.but.aucs.mir = c(auc.mir.cont,auc.mir.but556,auc.mir.but2417,auc.mir.but3906,auc.mir.but5560,auc.mir.but8703)/auc.mir.cont

#Compile L.2 data for functional responses #############
but.mir = data.frame(e = c(coef(gaf.mir.cont)[2], coef(gaf.mir.but556)[2], coef(gaf.mir.but2417)[2],
                           coef(gaf.mir.but3906)[2], coef(gaf.mir.but5560)[2], coef(gaf.mir.but8703)[2]),
                     e.se = c(summary(gaf.mir.cont)$coefficients[2,2], 
                              summary(gaf.mir.but556)$coefficients[2,2],
                              summary(gaf.mir.but2417)$coefficients[2,2], 
                              summary(gaf.mir.but3906)$coefficients[2,2],
                              summary(gaf.mir.but5560)$coefficients[2,2],
                              summary(gaf.mir.but8703)$coefficients[2,2]),
                     b = c(coef(gaf.mir.cont)[1], coef(gaf.mir.but556)[1], coef(gaf.mir.but2417)[1],
                           coef(gaf.mir.but3906)[1], coef(gaf.mir.but5560)[1], coef(gaf.mir.but8703)[1]),
                     b.se = c(summary(gaf.mir.cont)$coefficients[1,2], 
                              summary(gaf.mir.but556)$coefficients[1,2],
                              summary(gaf.mir.but2417)$coefficients[1,2], 
                              summary(gaf.mir.but3906)$coefficients[1,2],
                              summary(gaf.mir.but5560)$coefficients[1,2],
                              summary(gaf.mir.but8703)$coefficients[1,2]),
                     butralin = but.vals,
                     logbutralin = log(but.vals+1))

plot(but.mir$butralin, but.mir$e, pch = 16, xlab = 'butralin (ppb)', 
     ylab = 'L.2 Parameters (miracidia)', ylim = c(0,8))

  points(but.mir$butralin+100, but.mir$b, pch = 17, col=2)

  for(i in 1:length(but.mir$butralin)){
    segments(x0 = but.mir$butralin[i], y0 = but.mir$e[i] + but.mir$e.se[i],
             x1 = but.mir$butralin[i], y1 = but.mir$e[i] - but.mir$e.se[i])
    segments(x0 = but.mir$butralin[i]+100, y0 = but.mir$b[i] + but.mir$b.se[i],
             x1 = but.mir$butralin[i]+100, y1 = but.mir$b[i] - but.mir$b.se[i], col=2)
  } 
  
#functions fit to L.2 parameters ############
but.mir.lm.e = lm(e ~ butralin, weights = e.se^-1, data = but.mir) #linear response of LC50
  but.mir.pred = function(but){
    predict(but.mir.lm.e, newdata = data.frame(butralin = but), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.mir.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.mir.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.mir.pred, simplify = T)[3,], lty = 3)

but.mir.lm.e2 = lm(e ~ logbutralin, weights = e.se^-1, data = but.mir) #log-linear response of LC50
  but.mir.pred2 = function(but){
    predict(but.mir.lm.e2, newdata = data.frame(logbutralin = log(but+1)), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.mir.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.mir.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.mir.pred2, simplify = T)[3,], lty = 3, col=3)

AIC(but.mir.lm.e, but.mir.lm.e2)  #Linear is a better fit    

but.mir.lm.b = lm(b ~ butralin, weights = b.se^-1, data = but.mir)   
  but.mir.pred.b = function(but){
    predict(but.mir.lm.b, newdata = data.frame(butralin = but), interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.mir.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.mir.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.mir.pred.b, simplify = T)[3,], lty = 3, col = 2)

legend('topleft', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
       legend = c('LC50 - Linear',
                  'LC50 - Exponential',
                  'Slp - linear',
                  '95% CI'), cex = 0.7)  
legend('top', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
       cex = 0.8, bty = 'n', title = 'Observed points')

but.piM.pred = function(He){
  e = as.numeric(predict(but.mir.lm.e, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
  b = as.numeric(predict(but.mir.lm.b, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  
  while(e.use < 0){
    e.use = rnorm(1, e[1], e[2])
  }      #resample if negative
  
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use,
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piM.ghaf_butr.lin_unc = function(He){
  if(He == 0) piC = 1
  else(piC = but.piM.pred(He) / but.piM.pred(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

but.piM.pred2 = function(He){
  e = as.numeric(predict(but.mir.lm.e2, newdata = data.frame(logbutralin = log(He+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(but.mir.lm.b, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  
  while(e.use < 0){
    e.use = rnorm(1, e[1], e[2])
  }      #resample if negative
  
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use,
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piM.ghaf_butr.exp_unc = function(He){
  if(He == 0) piC = 1
  else(piC = but.piM.pred2(He) / but.piM.pred2(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

#plot to test function ########
plot(but.mir$butralin, c(auc.mir.cont, auc.mir.but556, auc.mir.but2417, 
                          auc.mir.but3906, auc.mir.but5560, auc.mir.but8703)/auc.mir.cont,
     xlab = 'Butralin (ppb)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0, 18000,50)/2, sapply(seq(0, 18000,50)/2, piM.ghaf_butr.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0, 18000,50)/2, sapply(seq(0, 18000,50)/2, piM.ghaf_butr.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)

#glyphosate ###############
  plot(mir$time_hrs[mir$conc==0], mir$surv[mir$conc==0], 
       pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 1:length(unique(mir$conc[mir$chem == 'glyphosate']))){
    points(mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[i] & mir$chem == 'glyphosate'], 
           mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[i] & mir$chem == 'glyphosate'], pch=17,
           col = i+1)
  }
    
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.cont, newdata = data.frame(conc = seq(0,24,0.1))), lty=2)
    
#1506 ppb *********************************************************************************
  gaf.mir.gly1506<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[1] & mir$chem == 'glyphosate'] ~ 
                        mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[1] & mir$chem == 'glyphosate'],
                      type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.gly1506, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 2)
  
  auc.mir.gly1506=integrate(f = L.3.fx, lc50 = coef(gaf.mir.gly1506)[2], slp = coef(gaf.mir.gly1506)[1],
                           lower=0, upper=24)[1]$value
  
#3875 ppb *********************************************************************************
  gaf.mir.gly3875<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[2] & mir$chem == 'glyphosate'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[2] & mir$chem == 'glyphosate'],
                       type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.gly3875, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 3)
  
  auc.mir.gly3875=integrate(f = L.3.fx, lc50 = coef(gaf.mir.gly3875)[2], slp = coef(gaf.mir.gly3875)[1],
                            lower=0, upper=24)[1]$value
  
#9174 ppb ***********************************************************************************
  gaf.mir.gly9174<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[3] & mir$chem == 'glyphosate'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[3] & mir$chem == 'glyphosate'],
                       type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.gly9174, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 4)
  
  auc.mir.gly9174=integrate(f = L.3.fx, lc50 = coef(gaf.mir.gly9174)[2], slp = coef(gaf.mir.gly9174)[1],
                            lower=0, upper=24)[1]$value
  
#15062 ppb **********************************************************************************
  gaf.mir.gly15062<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[4] & mir$chem == 'glyphosate'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[4] & mir$chem == 'glyphosate'],
                       type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.gly15062, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 5)
  
  auc.mir.gly15062=integrate(f = L.3.fx, lc50 = coef(gaf.mir.gly15062)[2], slp = coef(gaf.mir.gly15062)[1],
                            lower=0, upper=24)[1]$value    
  
#26249 ppb **********************************************************************************
  gaf.mir.gly26249<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[5] & mir$chem == 'glyphosate'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[5] & mir$chem == 'glyphosate'],
                       type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.gly26249, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 6)
  
  auc.mir.gly26249=integrate(f = L.3.fx, lc50 = coef(gaf.mir.gly26249)[2], slp = coef(gaf.mir.gly26249)[1],
                             lower=0, upper=24)[1]$value   
  
  title('Ghaffar2016 glyphosate toxicity to miracidia')
  legend('topright', legend = c('control', 1506,3875,9174,15062,26249), pch = c(16,17,17,17,17,17), 
         col = c(1:6), cex=0.8, bty = 'n') 
  
  gaf.gly.aucs.mir = c(auc.mir.cont, auc.mir.gly1506, auc.mir.gly3875, auc.mir.gly9174, auc.mir.gly15062, auc.mir.gly26249)/auc.mir.cont

#Compile L.2 parameters ######## 
gly.mir = data.frame(e = c(coef(gaf.mir.cont)[2], coef(gaf.mir.gly1506)[2], coef(gaf.mir.gly3875)[2],
                           coef(gaf.mir.gly9174)[2], coef(gaf.mir.gly15062)[2], coef(gaf.mir.gly26249)[2]),
                     e.se = c(summary(gaf.mir.cont)$coefficients[2,2], 
                              summary(gaf.mir.gly1506)$coefficients[2,2],
                              summary(gaf.mir.gly3875)$coefficients[2,2], 
                              summary(gaf.mir.gly9174)$coefficients[2,2],
                              summary(gaf.mir.gly15062)$coefficients[2,2],
                              summary(gaf.mir.gly26249)$coefficients[2,2]),
                     b = c(coef(gaf.mir.cont)[1], coef(gaf.mir.gly1506)[1], coef(gaf.mir.gly3875)[1],
                           coef(gaf.mir.gly9174)[1], coef(gaf.mir.gly15062)[1], coef(gaf.mir.gly26249)[1]),
                     b.se = c(summary(gaf.mir.cont)$coefficients[1,2], 
                              summary(gaf.mir.gly1506)$coefficients[1,2],
                              summary(gaf.mir.gly3875)$coefficients[1,2], 
                              summary(gaf.mir.gly9174)$coefficients[1,2],
                              summary(gaf.mir.gly15062)$coefficients[1,2],
                              summary(gaf.mir.gly26249)$coefficients[1,2]),
                     glyphosate = gly.vals,
                     logglyphosate = log(gly.vals+1))

plot(gly.mir$glyphosate, gly.mir$e, pch = 16, xlab = 'glyphosate (ppb)', 
     ylab = 'LL.2 Parameters (miracidia)', ylim = c(0,10))

  points(gly.mir$glyphosate+200, gly.mir$b, pch = 17, col=2)

for(i in 1:length(gly.mir$glyphosate)){
  segments(x0 = gly.mir$glyphosate[i], y0 = gly.mir$e[i] + gly.mir$e.se[i],
           x1 = gly.mir$glyphosate[i], y1 = gly.mir$e[i] - gly.mir$e.se[i])
  segments(x0 = gly.mir$glyphosate[i]+200, y0 = gly.mir$b[i] + gly.mir$b.se[i],
           x1 = gly.mir$glyphosate[i]+200, y1 = gly.mir$b[i] - gly.mir$b.se[i], col=2)
} 

#fit functions to L.2 parameters ##########  
gly.mir.lm.e = lm(e ~ glyphosate, weights = e.se^-1, data = gly.mir) #linear response of LC50
  gly.mir.pred = function(gly){
    predict(gly.mir.lm.e, newdata = data.frame(glyphosate = gly), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.mir.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.mir.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.mir.pred, simplify = T)[3,], lty = 3)

gly.mir.lm.e2 = lm(e ~ logglyphosate, weights = e.se^-1, data = gly.mir) #log-linear response of LC50
  gly.mir.pred2 = function(gly){
    predict(gly.mir.lm.e2, newdata = data.frame(logglyphosate = log(gly+1)), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.mir.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.mir.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.mir.pred2, simplify = T)[3,], lty = 3, col=3)

AIC(gly.mir.lm.e, gly.mir.lm.e2)  #Exponential is a better fit    

gly.mir.lm.b = lm(b ~ glyphosate, weights = b.se^-1, data = gly.mir)   
  gly.mir.pred.b = function(gly){
    predict(gly.mir.lm.b, newdata = data.frame(glyphosate = gly), interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.mir.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.mir.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.mir.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

gly.piM.pred = function(He){
  e = as.numeric(predict(gly.mir.lm.e, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
  b = as.numeric(predict(gly.mir.lm.b, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  
  while(e.use < 0){
    e.use = rnorm(1, e[1], e[2])
  }      #resample if negative
  
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piM.ghaf_gly.lin_unc = function(He){
  if(He == 0) piC = 1
  else(piC = gly.piM.pred(He) / gly.piM.pred(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

gly.piM.pred2 = function(He){
  e = as.numeric(predict(gly.mir.lm.e2, newdata = data.frame(logglyphosate = log(He+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(gly.mir.lm.b, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  
  while(e.use < 0){
    e.use = rnorm(1, e[1], e[2])
  }      #resample if negative
  
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piM.ghaf_gly.exp_unc = function(He){
  if(He == 0) piC = 1
  else(piC = gly.piM.pred2(He) / gly.piM.pred2(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

#plot to test function ########
plot(gly.mir$glyphosate, gaf.gly.aucs.mir,
     xlab = 'Glyphosate (ppb)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0,9000,50)*3, sapply(seq(0,9000,50)*3, piM.ghaf_gly.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,9000,50)*3, sapply(seq(0,9000,50)*3, piM.ghaf_gly.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)

#pendimethalin ###############
  plot(mir$time_hrs[mir$conc==0], mir$surv[mir$conc==0], 
       pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 1:length(unique(mir$conc[mir$chem == 'pendimethalin']))){
    points(mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[i] & mir$chem == 'pendimethalin'], 
           mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[i] & mir$chem == 'pendimethalin'], pch=17,
           col = i+1)
  }
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.cont, newdata = data.frame(conc = seq(0,24,0.1))), lty=2)
  
#215 ppb *********************************************************************************
  gaf.mir.pen215<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[1] & mir$chem == 'pendimethalin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[1] & mir$chem == 'pendimethalin'],
                       type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.pen215, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 2)
  
  auc.mir.pen215=integrate(f = L.3.fx, lc50 = coef(gaf.mir.pen215)[2], slp = coef(gaf.mir.pen215)[1],
                            lower=0, upper=24)[1]$value  
  
#535 ppb *********************************************************************************
  gaf.mir.pen535<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[2] & mir$chem == 'pendimethalin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[2] & mir$chem == 'pendimethalin'],
                       type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.pen535, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 3)
  
  auc.mir.pen535=integrate(f = L.3.fx, lc50 = coef(gaf.mir.pen535)[2], slp = coef(gaf.mir.pen535)[1],
                           lower=0, upper=24)[1]$value 
  
#1299 ppb ***********************************************************************************
  gaf.mir.pen1299<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[3] & mir$chem == 'pendimethalin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[3] & mir$chem == 'pendimethalin'],
                       type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.pen1299, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 4)
  
  auc.mir.pen1299=integrate(f = L.3.fx, lc50 = coef(gaf.mir.pen1299)[2], slp = coef(gaf.mir.pen1299)[1],
                           lower=0, upper=24)[1]$value 
  
#2148 ppb **********************************************************************************
  gaf.mir.pen2148<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[4] & mir$chem == 'pendimethalin'] ~ 
                          mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[4] & mir$chem == 'pendimethalin'],
                        type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.pen2148, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 5)
  
  auc.mir.pen2148=integrate(f = L.3.fx, lc50 = coef(gaf.mir.pen2148)[2], slp = coef(gaf.mir.pen2148)[1],
                            lower=0, upper=24)[1]$value   
  
#3762 ppb **********************************************************************************
  gaf.mir.pen3762<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[5] & mir$chem == 'pendimethalin'] ~ 
                          mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[5] & mir$chem == 'pendimethalin'],
                        type = 'binomial', fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  
  lines(x=seq(0,24,0.1), y=predict(gaf.mir.pen3762, newdata = data.frame(conc = seq(0,24,0.1))), 
        lty=2, col = 6)
  
  auc.mir.pen3762=integrate(f = L.3.fx, lc50 = coef(gaf.mir.pen3762)[2], slp = coef(gaf.mir.pen3762)[1],
                            lower=0, upper=24)[1]$value   
  
  title('Ghaffar2016 pendimethalin toxicity to miracidia')
  legend('topright', legend = c('control', 215,535,1299,2148,3762), pch = c(16,17,17,17,17,17), 
         col = c(1:6), cex=0.8, bty='n')  
  
  gaf.pen.aucs.mir = c(auc.mir.cont, auc.mir.pen215, auc.mir.pen535, auc.mir.pen1299, auc.mir.pen2148, auc.mir.pen3762)/auc.mir.cont
  
#Compile L.2 parameters ######## 
pen.mir = data.frame(e = c(coef(gaf.mir.cont)[2], coef(gaf.mir.pen215)[2], coef(gaf.mir.pen535)[2],
                           coef(gaf.mir.pen1299)[2], coef(gaf.mir.pen2148)[2], coef(gaf.mir.pen3762)[2]),
                     e.se = c(summary(gaf.mir.cont)$coefficients[2,2], 
                              summary(gaf.mir.pen215)$coefficients[2,2],
                              summary(gaf.mir.pen535)$coefficients[2,2], 
                              summary(gaf.mir.pen1299)$coefficients[2,2],
                              summary(gaf.mir.pen2148)$coefficients[2,2],
                              summary(gaf.mir.gly26249)$coefficients[2,2]),
                     b = c(coef(gaf.mir.cont)[1], coef(gaf.mir.pen215)[1], coef(gaf.mir.pen535)[1],
                           coef(gaf.mir.pen1299)[1], coef(gaf.mir.pen2148)[1], coef(gaf.mir.pen3762)[1]),
                     b.se = c(summary(gaf.mir.cont)$coefficients[1,2], 
                              summary(gaf.mir.pen215)$coefficients[1,2],
                              summary(gaf.mir.pen535)$coefficients[1,2], 
                              summary(gaf.mir.pen1299)$coefficients[1,2],
                              summary(gaf.mir.pen2148)$coefficients[1,2],
                              summary(gaf.mir.pen3762)$coefficients[1,2]),
                     pendimethalin = pen.vals,
                     logpendimethalin = log(pen.vals+1))

plot(pen.mir$pendimethalin, pen.mir$e, pch = 16, xlab = 'Pendimethalin (ppb)', 
     ylab = 'L.2 Parameters (miracidia)', ylim = c(0,12))

  points(pen.mir$pendimethalin+50, pen.mir$b, pch = 17, col=2)

  for(i in 1:length(pen.mir$pendimethalin)){
    segments(x0 = pen.mir$pendimethalin[i], y0 = pen.mir$e[i] + pen.mir$e.se[i],
             x1 = pen.mir$pendimethalin[i], y1 = pen.mir$e[i] - pen.mir$e.se[i])
    segments(x0 = pen.mir$pendimethalin[i]+50, y0 = pen.mir$b[i] + pen.mir$b.se[i],
             x1 = pen.mir$pendimethalin[i]+50, y1 = pen.mir$b[i] - pen.mir$b.se[i], col=2)
  } 

#fit functions to L.2 parameters ##########  
pen.mir.lm.e = lm(e ~ pendimethalin, weights = e.se^-1, data = pen.mir) #linear response of LC50
  pen.mir.pred = function(pen){
    predict(pen.mir.lm.e, newdata = data.frame(pendimethalin = pen), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,4000,25), sapply(seq(0,4000,25), pen.mir.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,4000,25), sapply(seq(0,4000,25), pen.mir.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,4000,25), sapply(seq(0,4000,25), pen.mir.pred, simplify = T)[3,], lty = 3)

pen.mir.lm.e2 = lm(e ~ logpendimethalin, weights = e.se^-1, data = pen.mir) #log-linear response of LC50
  pen.mir.pred2 = function(pen){
    predict(pen.mir.lm.e2, newdata = data.frame(logpendimethalin = log(pen+1)), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,4000,25), sapply(seq(0,4000,25), pen.mir.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,4000,25), sapply(seq(0,4000,25), pen.mir.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,4000,25), sapply(seq(0,4000,25), pen.mir.pred2, simplify = T)[3,], lty = 3, col=3)

AIC(pen.mir.lm.e, pen.mir.lm.e2)  #Exponential is a better fit    

pen.mir.lm.b = lm(b ~ pendimethalin, weights = b.se^-1, data = pen.mir)   
  pen.mir.pred.b = function(pen){
    predict(pen.mir.lm.b, newdata = data.frame(pendimethalin = pen), interval = 'confidence', level = 0.95)
  }

    lines(seq(0,4000,25), sapply(seq(0,4000,25), pen.mir.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,4000,25), sapply(seq(0,4000,25), pen.mir.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,4000,25), sapply(seq(0,4000,25), pen.mir.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

pen.piM.pred = function(He){
  e = as.numeric(predict(pen.mir.lm.e, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
  b = as.numeric(predict(pen.mir.lm.b, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  
  while(e.use < 0){
    e.use = rnorm(1, e[1], e[2])
  }      #resample if negative
  
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piM.ghaf_pen.lin_unc = function(He){
  if(He == 0) piC = 1
  else(piC = pen.piM.pred(He) / pen.piM.pred(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

pen.piM.pred2 = function(He){
  e = as.numeric(predict(pen.mir.lm.e2, newdata = data.frame(logpendimethalin = log(He+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(pen.mir.lm.b, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  
  while(e.use < 0){
    e.use = rnorm(1, e[1], e[2])
  }      #resample if negative
  
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piM.ghaf_pen.exp_unc = function(He){
  if(He == 0) piC = 1
  else(piC = pen.piM.pred2(He) / pen.piM.pred2(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate


#plot to test function ########
plot(pen.mir$pendimethalin, gaf.pen.aucs.mir,
     xlab = 'pendimethalin (ppb)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0,4000,10), sapply(seq(0,4000,10), piM.ghaf_pen.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,4000,10), sapply(seq(0,4000,10), piM.ghaf_pen.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)
