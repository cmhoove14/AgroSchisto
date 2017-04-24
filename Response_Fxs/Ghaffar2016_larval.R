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
                    type = 'binomial', fct = LL.2())
  
  gaf.mir.cont.surv = function(t){
    (1/(1+exp(gaf.mir.cont$coefficients[1]*(log(t/gaf.mir.cont$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.cont.surv(seq(0,24,0.1)), lty=2)
  
  auc.mir.cont=integrate(f = gaf.mir.cont.surv, lower=0, upper=24)[1]$value  

#556 ppb *********************************************************************************
  gaf.mir.but556<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[1] & mir$chem == 'butralin'] ~ 
                        mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[1] & mir$chem == 'butralin'],
                      type = 'binomial', fct = LL.2())
  
  gaf.mir.but556.surv = function(t){
    (1/(1+exp(gaf.mir.but556$coefficients[1]*(log(t/gaf.mir.but556$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.but556.surv(seq(0,24,0.1)), lty=2, col=2)
  
  auc.mir.but556=integrate(f = gaf.mir.but556.surv, lower=0, upper=24)[1]$value  

#2417 ppb *********************************************************************************
  gaf.mir.but2417<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[2] & mir$chem == 'butralin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[2] & mir$chem == 'butralin'],
                       type = 'binomial', fct = LL.2())
  
  gaf.mir.but2417.surv = function(t){
    (1/(1+exp(gaf.mir.but2417$coefficients[1]*(log(t/gaf.mir.but2417$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.but2417.surv(seq(0,24,0.1)), lty=2, col=3)
  
  auc.mir.but2417=integrate(f = gaf.mir.but2417.surv, lower=0, upper=24)[1]$value  

#3906 ppb ***********************************************************************************
  gaf.mir.but3906<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[3] & mir$chem == 'butralin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[3] & mir$chem == 'butralin'],
                       type = 'binomial', fct = LL.2())
  
  gaf.mir.but3906.surv = function(t){
    (1/(1+exp(gaf.mir.but3906$coefficients[1]*(log(t/gaf.mir.but3906$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.but3906.surv(seq(0,24,0.1)), lty=2, col=4)
  
  auc.mir.but3906=integrate(f = gaf.mir.but3906.surv, lower=0, upper=24)[1]$value  

#5560 ppb **********************************************************************************
  gaf.mir.but5560<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[4] & mir$chem == 'butralin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[4] & mir$chem == 'butralin'],
                       type = 'binomial', fct = LL.2())
  
  gaf.mir.but5560.surv = function(t){
    (1/(1+exp(gaf.mir.but5560$coefficients[1]*(log(t/gaf.mir.but5560$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.but5560.surv(seq(0,24,0.1)), lty=2, col=5)
  
  auc.mir.but5560=integrate(f = gaf.mir.but5560.surv, lower=0, upper=24)[1]$value     

#8703 ppb **********************************************************************************
  gaf.mir.but8703<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[5] & mir$chem == 'butralin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[5] & mir$chem == 'butralin'],
                       type = 'binomial', fct = LL.2())
  
  gaf.mir.but8703.surv = function(t){
    (1/(1+exp(gaf.mir.but8703$coefficients[1]*(log(t/gaf.mir.but8703$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.but8703.surv(seq(0,24,0.1)), lty=2, col=6)
  
  auc.mir.but8703=integrate(f = gaf.mir.but8703.surv, lower=0, upper=24)[1]$value   
  
  legend('topright', legend = c('control', 556,2417,3906,5560,8703), pch = c(16,17,17,17,17,17), col = c(1:6), cex=0.8)
  
  gaf.but.aucs.mir = c(auc.mir.cont,auc.mir.but556,auc.mir.but2417,auc.mir.but3906,auc.mir.but5560,auc.mir.but8703)/auc.mir.cont

#Compile LL.2 data for functional responses #############
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
     ylab = 'LL.2 Parameters (miracidia)', ylim = c(0,13))

  points(but.mir$butralin+100, but.mir$b, pch = 17, col=2)

  for(i in 1:length(but.mir$butralin)){
    segments(x0 = but.mir$butralin[i], y0 = but.mir$e[i] + but.mir$e.se[i],
             x1 = but.mir$butralin[i], y1 = but.mir$e[i] - but.mir$e.se[i])
    segments(x0 = but.mir$butralin[i]+100, y0 = but.mir$b[i] + but.mir$b.se[i],
             x1 = but.mir$butralin[i]+100, y1 = but.mir$b[i] - but.mir$b.se[i], col=2)
  } 
  
#functions fit to LL.2 parameters ############
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
    predict(but.mir.lm.e2, newdata = data.frame(logbutralin = log(but)+1), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.mir.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.mir.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.mir.pred2, simplify = T)[3,], lty = 3, col=3)

AIC(but.mir.lm.e, but.mir.lm.e2)  #Exponential is a better fit    

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
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
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
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
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
    lines(x=seq(0,24,0.1), y=gaf.mir.cont.surv(seq(0,24,0.1)), lty=2)

#1506 ppb *********************************************************************************
  gaf.mir.gly1506<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[1] & mir$chem == 'glyphosate'] ~ 
                        mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[1] & mir$chem == 'glyphosate'],
                      type = 'binomial', fct = LL.2())
  
  gaf.mir.gly1506.surv = function(t){
    (1/(1+exp(gaf.mir.gly1506$coefficients[1]*(log(t/gaf.mir.gly1506$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.gly1506.surv(seq(0,24,0.1)), lty=2, col=2)
  
  auc.mir.gly1506=integrate(f = gaf.mir.gly1506.surv, lower=0, upper=24)[1]$value  
  
#3875 ppb *********************************************************************************
  gaf.mir.gly3875<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[2] & mir$chem == 'glyphosate'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[2] & mir$chem == 'glyphosate'],
                       type = 'binomial', fct = LL.2())
  
  gaf.mir.gly3875.surv = function(t){
    (1/(1+exp(gaf.mir.gly3875$coefficients[1]*(log(t/gaf.mir.gly3875$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.gly3875.surv(seq(0,24,0.1)), lty=2, col=3)
  
  auc.mir.gly3875=integrate(f = gaf.mir.gly3875.surv, lower=0, upper=24)[1]$value  
  
#9174 ppb ***********************************************************************************
  gaf.mir.gly9174<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[3] & mir$chem == 'glyphosate'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[3] & mir$chem == 'glyphosate'],
                       type = 'binomial', fct = LL.2())
  
  gaf.mir.gly9174.surv = function(t){
    (1/(1+exp(gaf.mir.gly9174$coefficients[1]*(log(t/gaf.mir.gly9174$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.gly9174.surv(seq(0,24,0.1)), lty=2, col=4)
  
  auc.mir.gly9174=integrate(f = gaf.mir.gly9174.surv, lower=0, upper=24)[1]$value  
  
#15062 ppb **********************************************************************************
  gaf.mir.gly15062<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[4] & mir$chem == 'glyphosate'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[4] & mir$chem == 'glyphosate'],
                       type = 'binomial', fct = LL.2())
  
  gaf.mir.gly15062.surv = function(t){
    (1/(1+exp(gaf.mir.gly15062$coefficients[1]*(log(t/gaf.mir.gly15062$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.gly15062.surv(seq(0,24,0.1)), lty=2, col=5)
  
  auc.mir.gly15062=integrate(f = gaf.mir.gly15062.surv, lower=0, upper=24)[1]$value     
  
#26249 ppb **********************************************************************************
  gaf.mir.gly26249<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[5] & mir$chem == 'glyphosate'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[5] & mir$chem == 'glyphosate'],
                       type = 'binomial', fct = LL.2())
  
  gaf.mir.gly26249.surv = function(t){
    (1/(1+exp(gaf.mir.gly26249$coefficients[1]*(log(t/gaf.mir.gly26249$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.gly26249.surv(seq(0,24,0.1)), lty=2, col=6)
  
  auc.mir.gly26249=integrate(f = gaf.mir.gly26249.surv, lower=0, upper=24)[1]$value   
  
  title('Ghaffar2016 glyphosate toxicity to miracidia')
  legend('topright', legend = c('control', 1506,3875,9174,15062,26249), pch = c(16,17,17,17,17,17), col = c(1:6), cex=0.8) 
  
  gaf.gly.aucs.mir = c(auc.mir.cont, auc.mir.gly1506, auc.mir.gly3875, auc.mir.gly9174, auc.mir.gly15062, auc.mir.gly26249)/auc.mir.cont

#Compile LL.2 parameters ######## 
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
     ylab = 'LL.2 Parameters (miracidia)', ylim = c(0,8))

  points(gly.mir$glyphosate+200, gly.mir$b, pch = 17, col=2)

for(i in 1:length(gly.mir$glyphosate)){
  segments(x0 = gly.mir$glyphosate[i], y0 = gly.mir$e[i] + gly.mir$e.se[i],
           x1 = gly.mir$glyphosate[i], y1 = gly.mir$e[i] - gly.mir$e.se[i])
  segments(x0 = gly.mir$glyphosate[i]+200, y0 = gly.mir$b[i] + gly.mir$b.se[i],
           x1 = gly.mir$glyphosate[i]+200, y1 = gly.mir$b[i] - gly.mir$b.se[i], col=2)
} 

#fit functions to LL.2 parameters ##########  
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
    predict(gly.mir.lm.e2, newdata = data.frame(logglyphosate = log(gly)+1), 
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
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
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
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
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
  points(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, piM.ghaf_gly.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, piM.ghaf_gly.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)

#pendimethalin ###############
  plot(mir$time_hrs[mir$conc==0], mir$surv[mir$conc==0], 
       pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 1:length(unique(mir$conc[mir$chem == 'pendimethalin']))){
    points(mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[i] & mir$chem == 'pendimethalin'], 
           mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[i] & mir$chem == 'pendimethalin'], pch=17,
           col = i+1)
  }
  lines(x=seq(0,24,0.1), y=gaf.mir.cont.surv(seq(0,24,0.1)), lty=2)
  
#215 ppb *********************************************************************************
  gaf.mir.pen215<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[1] & mir$chem == 'pendimethalin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[1] & mir$chem == 'pendimethalin'],
                       type = 'binomial', fct = LL.2())
  
  gaf.mir.pen215.surv = function(t){
    (1/(1+exp(gaf.mir.pen215$coefficients[1]*(log(t/gaf.mir.pen215$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.pen215.surv(seq(0,24,0.1)), lty=2, col=2)
  
  auc.mir.pen215=integrate(f = gaf.mir.pen215.surv, lower=0, upper=24)[1]$value  
  
#535 ppb *********************************************************************************
  gaf.mir.pen535<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[2] & mir$chem == 'pendimethalin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[2] & mir$chem == 'pendimethalin'],
                       type = 'binomial', fct = LL.2())
  
  gaf.mir.pen535.surv = function(t){
    (1/(1+exp(gaf.mir.pen535$coefficients[1]*(log(t/gaf.mir.pen535$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.pen535.surv(seq(0,24,0.1)), lty=2, col=3)
  
  auc.mir.pen535=integrate(f = gaf.mir.pen535.surv, lower=0, upper=24)[1]$value  
  
#1299 ppb ***********************************************************************************
  gaf.mir.pen1299<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[3] & mir$chem == 'pendimethalin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[3] & mir$chem == 'pendimethalin'],
                       type = 'binomial', fct = LL.2())
  
  gaf.mir.pen1299.surv = function(t){
    (1/(1+exp(gaf.mir.pen1299$coefficients[1]*(log(t/gaf.mir.pen1299$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.pen1299.surv(seq(0,24,0.1)), lty=2, col=4)
  
  auc.mir.pen1299=integrate(f = gaf.mir.pen1299.surv, lower=0, upper=24)[1]$value  
  
#2148 ppb **********************************************************************************
  gaf.mir.pen2148<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[4] & mir$chem == 'pendimethalin'] ~ 
                          mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[4] & mir$chem == 'pendimethalin'],
                        type = 'binomial', fct = LL.2())
  
  gaf.mir.pen2148.surv = function(t){
    (1/(1+exp(gaf.mir.pen2148$coefficients[1]*(log(t/gaf.mir.pen2148$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.pen2148.surv(seq(0,24,0.1)), lty=2, col=5)
  
  auc.mir.pen2148=integrate(f = gaf.mir.pen2148.surv, lower=0, upper=24)[1]$value     
  
#3762 ppb **********************************************************************************
  gaf.mir.pen3762<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[5] & mir$chem == 'pendimethalin'] ~ 
                          mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[5] & mir$chem == 'pendimethalin'],
                        type = 'binomial', fct = LL.2())
  
  gaf.mir.pen3762.surv = function(t){
    (1/(1+exp(gaf.mir.pen3762$coefficients[1]*(log(t/gaf.mir.pen3762$coefficients[2])))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.pen3762.surv(seq(0,24,0.1)), lty=2, col=6)
  
  auc.mir.pen3762=integrate(f = gaf.mir.pen3762.surv, lower=0, upper=24)[1]$value   
  
  title('Ghaffar2016 pendimethalin toxicity to miracidia')
  legend('topright', legend = c('control', 215,535,1299,2148,3762), pch = c(16,17,17,17,17,17), col = c(1:6), cex=0.8)  
  
  gaf.pen.aucs.mir = c(auc.mir.cont, auc.mir.pen215, auc.mir.pen535, auc.mir.pen1299, auc.mir.pen2148, auc.mir.pen3762)/auc.mir.cont
  
#Compile LL.2 parameters ######## 
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
     ylab = 'LL.2 Parameters (miracidia)', ylim = c(0,12))

points(pen.mir$pendimethalin+100, pen.mir$b, pch = 17, col=2)

for(i in 1:length(pen.mir$pendimethalin)){
  segments(x0 = pen.mir$pendimethalin[i], y0 = pen.mir$e[i] + pen.mir$e.se[i],
           x1 = pen.mir$pendimethalin[i], y1 = pen.mir$e[i] - pen.mir$e.se[i])
  segments(x0 = pen.mir$pendimethalin[i]+100, y0 = pen.mir$b[i] + pen.mir$b.se[i],
           x1 = pen.mir$pendimethalin[i]+100, y1 = pen.mir$b[i] - pen.mir$b.se[i], col=2)
} 

#fit functions to LL.2 parameters ##########  
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
    predict(gly.mir.lm.e2, newdata = data.frame(logglyphosate = log(gly)+1), 
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
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
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
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piM.ghaf_gly.exp_unc = function(He){
  if(He == 0) piC = 1
  else(piC = gly.piM.pred2(He) / gly.piM.pred2(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate



#Fit pi_M functions ################  
  mirgaf.df = data.frame('but.conc'=but.vals,
                         'gly.conc'=gly.vals,
                         'pen.conc'=pen.vals,
                         'but.auc'=gaf.but.aucs.mir,
                         'gly.auc'=gaf.gly.aucs.mir,
                         'pen.auc'=gaf.pen.aucs.mir)
  
  plot(mirgaf.df$gly.conc, mirgaf.df$gly.auc, pch = 16, col = 'firebrick', xlim = c(0,27000), ylim = c(0,1),
       xlab = 'Herbicide Concentration', ylab = 'relative auc')
    points(mirgaf.df$but.conc, mirgaf.df$but.auc, pch = 16, col = 'darkorange')
    points(mirgaf.df$pen.conc, mirgaf.df$pen.auc, pch = 16, col = 'indianred2')
#Glyphosate   
  gly.piM.pred = nls(gly.auc ~ exp(-b*gly.conc), data=mirgaf.df, start = list(b=0.0001))
    summary(gly.piM.pred)
    
  pi_M_gly_gaf16 = function(He){
      exp(-summary(gly.piM.pred)$parameters[1]*He) 
  }
  
  gly.piM.pred.l2 = drm(mirgaf.df$gly.auc  ~ mirgaf.df$gly.conc, type = 'binomial', fct = LL2.2())
  
  pi_M_gly_gaf16.l2 = function(He){
    (1/(1+exp(gly.piM.pred.l2$coefficients[1]*(log(He)-gly.piM.pred.l2$coefficients[2]))))
  }
    
    lines(c(0:27000), pi_M_gly_gaf16(c(0:27000)), lty=2, col='firebrick')  
    lines(c(0:27000), pi_M_gly_gaf16.l2(c(0:27000)), lty=3, col='firebrick')  
    
#butralin   
  but.piM.pred = nls(but.auc ~ exp(-b*but.conc), data=mirgaf.df, start = list(b=0.001))
    summary(but.piM.pred)
    
  pi_M_but_gaf16 = function(He){
      exp(-summary(but.piM.pred)$parameters[1]*He) 
  }
  
  but.piM.pred.l2 = drm(mirgaf.df$but.auc  ~ mirgaf.df$but.conc, type = 'binomial', fct = LL2.2())
  
  pi_M_but_gaf16.l2 = function(He){
    (1/(1+exp(but.piM.pred.l2$coefficients[1]*(log(He)-but.piM.pred.l2$coefficients[2]))))
  }
    
    lines(c(0:27000), pi_M_but_gaf16(c(0:27000)), lty=2, col='darkorange')  
    lines(c(0:27000), pi_M_but_gaf16.l2(c(0:27000)), lty=3, col='darkorange')  
    
#pendimethalin   
  pen.piM.pred = nls(pen.auc ~ exp(-b*pen.conc), data=mirgaf.df, start = list(b=0.005))
    summary(pen.piM.pred)
    
  pi_M_pen_gaf16 = function(He){
      exp(-summary(pen.piM.pred)$parameters[1]*He) 
  }
  
  pen.piM.pred.l2 = drm(mirgaf.df$pen.auc  ~ mirgaf.df$pen.conc, type = 'binomial', fct = LL2.2())
  
  pi_M_pen_gaf16.l2 = function(He){
    (1/(1+exp(pen.piM.pred.l2$coefficients[1]*(log(He)-pen.piM.pred.l2$coefficients[2]))))
  }
    
    lines(c(0:27000), pi_M_pen_gaf16(c(0:27000)), lty=2, col='indianred2')  
    lines(c(0:27000), pi_M_pen_gaf16.l2(c(0:27000)), lty=3, col='indianred2')  

#Cercarial (S. mansoni) toxicity ###############
  cer<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/ghaffar2016_cercariae.csv')
    
#butralin ###############
  plot(cer$time_hrs[cer$conc==0], cer$surv[cer$conc==0], 
       pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 1:length(unique(cer$conc[cer$chem == 'butralin']))){
      points(cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[i] & cer$chem == 'butralin'], 
             cer$surv[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[i] & cer$chem == 'butralin'], pch=17,
             col = i+1)
    }
    
#fit to control points *****************************************************************
  gaf.cer.cont<-drm(cer$surv[cer$conc==0] ~ cer$time_hrs[cer$conc==0],
                      type = 'binomial', fct = LL2.2())
    
    gaf.cer.cont.surv = function(t){
      (1/(1+exp(gaf.cer.cont$coefficients[1]*(log(t)-gaf.cer.cont$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.cont.surv(seq(0,24,0.1)), lty=2)
    
    auc.cer.cont=integrate(f = gaf.cer.cont.surv, lower=0, upper=24)[1]$value  
    
#556 ppb *********************************************************************************
    gaf.cer.but556<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[1] & cer$chem == 'butralin'] ~ 
                          cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[1] & cer$chem == 'butralin'],
                        type = 'binomial', fct = LL2.2())
    
    gaf.cer.but556.surv = function(t){
      (1/(1+exp(gaf.cer.but556$coefficients[1]*(log(t)-gaf.cer.but556$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.but556.surv(seq(0,24,0.1)), lty=2, col=2)
    
    auc.cer.but556=integrate(f = gaf.cer.but556.surv, lower=0, upper=24)[1]$value  
    
#2417 ppb *********************************************************************************
    gaf.cer.but2417<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[2] & cer$chem == 'butralin'] ~ 
                           cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[2] & cer$chem == 'butralin'],
                         type = 'binomial', fct = LL2.2())
    
    gaf.cer.but2417.surv = function(t){
      (1/(1+exp(gaf.cer.but2417$coefficients[1]*(log(t)-gaf.cer.but2417$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.but2417.surv(seq(0,24,0.1)), lty=2, col=3)
    
    auc.cer.but2417=integrate(f = gaf.cer.but2417.surv, lower=0, upper=24)[1]$value  
    
#3906 ppb ***********************************************************************************
    gaf.cer.but3906<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[3] & cer$chem == 'butralin'] ~ 
                           cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[3] & cer$chem == 'butralin'],
                         type = 'binomial', fct = LL2.2())
    
    gaf.cer.but3906.surv = function(t){
      (1/(1+exp(gaf.cer.but3906$coefficients[1]*(log(t)-gaf.cer.but3906$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.but3906.surv(seq(0,24,0.1)), lty=2, col=4)
    
    auc.cer.but3906=integrate(f = gaf.cer.but3906.surv, lower=0, upper=24)[1]$value  
    
#5560 ppb **********************************************************************************
    gaf.cer.but5560<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[4] & cer$chem == 'butralin'] ~ 
                           cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[4] & cer$chem == 'butralin'],
                         type = 'binomial', fct = LL2.2())
    
    gaf.cer.but5560.surv = function(t){
      (1/(1+exp(gaf.cer.but5560$coefficients[1]*(log(t)-gaf.cer.but5560$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.but5560.surv(seq(0,24,0.1)), lty=2, col=5)
    
    auc.cer.but5560=integrate(f = gaf.cer.but5560.surv, lower=0, upper=24)[1]$value     
    
#8703 ppb **********************************************************************************
    gaf.cer.but8703<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[5] & cer$chem == 'butralin'] ~ 
                           cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'butralin'])[5] & cer$chem == 'butralin'],
                         type = 'binomial', fct = LL2.2())
    
    gaf.cer.but8703.surv = function(t){
      (1/(1+exp(gaf.cer.but8703$coefficients[1]*(log(t)-gaf.cer.but8703$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.but8703.surv(seq(0,24,0.1)), lty=2, col=6)
    
    auc.cer.but8703=integrate(f = gaf.cer.but8703.surv, lower=0, upper=24)[1]$value   
    
    title('Ghaffar2016 butralin toxicity to cercariae')
    legend('topright', legend = c('control', 556,2417,3906,5560,8703), pch = c(16,17,17,17,17,17), col = c(1:6), cex=0.8)
    
    gaf.but.aucs.cer = c(auc.cer.cont,auc.cer.but556,auc.cer.but2417,auc.cer.but3906,auc.cer.but5560,auc.cer.but8703)/auc.cer.cont
#glyphosate ###############
    plot(cer$time_hrs[cer$conc==0], cer$surv[cer$conc==0], 
         pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    for(i in 1:length(unique(cer$conc[cer$chem == 'glyphosate']))){
      points(cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[i] & cer$chem == 'glyphosate'], 
             cer$surv[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[i] & cer$chem == 'glyphosate'], pch=17,
             col = i+1)
    }
    lines(x=seq(0,24,0.1), y=gaf.cer.cont.surv(seq(0,24,0.1)), lty=2)
    
#1506 ppb *********************************************************************************
    gaf.cer.gly1506<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[1] & cer$chem == 'glyphosate'] ~ 
                           cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[1] & cer$chem == 'glyphosate'],
                         type = 'binomial', fct = LL2.2())
    
    gaf.cer.gly1506.surv = function(t){
      (1/(1+exp(gaf.cer.gly1506$coefficients[1]*(log(t)-gaf.cer.gly1506$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.gly1506.surv(seq(0,24,0.1)), lty=2, col=2)
    
    auc.cer.gly1506=integrate(f = gaf.cer.gly1506.surv, lower=0, upper=24)[1]$value  
    
#3875 ppb *********************************************************************************
    gaf.cer.gly3875<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[2] & cer$chem == 'glyphosate'] ~ 
                           cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[2] & cer$chem == 'glyphosate'],
                         type = 'binomial', fct = LL2.2())
    
    gaf.cer.gly3875.surv = function(t){
      (1/(1+exp(gaf.cer.gly3875$coefficients[1]*(log(t)-gaf.cer.gly3875$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.gly3875.surv(seq(0,24,0.1)), lty=2, col=3)
    
    auc.cer.gly3875=integrate(f = gaf.cer.gly3875.surv, lower=0, upper=24)[1]$value  
    
#9174 ppb ***********************************************************************************
    gaf.cer.gly9174<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[3] & cer$chem == 'glyphosate'] ~ 
                           cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[3] & cer$chem == 'glyphosate'],
                         type = 'binomial', fct = LL2.2())
    
    gaf.cer.gly9174.surv = function(t){
      (1/(1+exp(gaf.cer.gly9174$coefficients[1]*(log(t)-gaf.cer.gly9174$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.gly9174.surv(seq(0,24,0.1)), lty=2, col=4)
    
    auc.cer.gly9174=integrate(f = gaf.cer.gly9174.surv, lower=0, upper=24)[1]$value  
    
#15062 ppb **********************************************************************************
    gaf.cer.gly15062<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[4] & cer$chem == 'glyphosate'] ~ 
                            cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[4] & cer$chem == 'glyphosate'],
                          type = 'binomial', fct = LL2.2())
    
    gaf.cer.gly15062.surv = function(t){
      (1/(1+exp(gaf.cer.gly15062$coefficients[1]*(log(t)-gaf.cer.gly15062$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.gly15062.surv(seq(0,24,0.1)), lty=2, col=5)
    
    auc.cer.gly15062=integrate(f = gaf.cer.gly15062.surv, lower=0, upper=24)[1]$value     
    
#26249 ppb **********************************************************************************
    gaf.cer.gly26249<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[5] & cer$chem == 'glyphosate'] ~ 
                            cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'glyphosate'])[5] & cer$chem == 'glyphosate'],
                          type = 'binomial', fct = LL2.2())
    
    gaf.cer.gly26249.surv = function(t){
      (1/(1+exp(gaf.cer.gly26249$coefficients[1]*(log(t)-gaf.cer.gly26249$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.gly26249.surv(seq(0,24,0.1)), lty=2, col=6)
    
    auc.cer.gly26249=integrate(f = gaf.cer.gly26249.surv, lower=0, upper=24)[1]$value   
    
    title('Ghaffar2016 glyphosate toxicity to cercariae')
    legend('topright', legend = c('control', 1506,3875,9174,15062,26249), pch = c(16,17,17,17,17,17), col = c(1:6), cex=0.8) 
    
    gaf.gly.aucs.cer = c(auc.cer.cont, auc.cer.gly1506, auc.cer.gly3875, auc.cer.gly9174, auc.cer.gly15062, auc.cer.gly26249)/auc.cer.cont
#pendimethalin ###############
    plot(cer$time_hrs[cer$conc==0], cer$surv[cer$conc==0], 
         pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    for(i in 1:length(unique(cer$conc[cer$chem == 'pendimethalin']))){
      points(cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[i] & cer$chem == 'pendimethalin'], 
             cer$surv[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[i] & cer$chem == 'pendimethalin'], pch=17,
             col = i+1)
    }
    lines(x=seq(0,24,0.1), y=gaf.cer.cont.surv(seq(0,24,0.1)), lty=2)
    
#215 ppb *********************************************************************************
    gaf.cer.pen215<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[1] & cer$chem == 'pendimethalin'] ~ 
                          cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[1] & cer$chem == 'pendimethalin'],
                        type = 'binomial', fct = LL2.2())
    
    gaf.cer.pen215.surv = function(t){
      (1/(1+exp(gaf.cer.pen215$coefficients[1]*(log(t)-gaf.cer.pen215$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.pen215.surv(seq(0,24,0.1)), lty=2, col=2)
    
    auc.cer.pen215=integrate(f = gaf.cer.pen215.surv, lower=0, upper=24)[1]$value  
    
#535 ppb *********************************************************************************
    gaf.cer.pen535<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[2] & cer$chem == 'pendimethalin'] ~ 
                          cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[2] & cer$chem == 'pendimethalin'],
                        type = 'binomial', fct = LL2.2())
    
    gaf.cer.pen535.surv = function(t){
      (1/(1+exp(gaf.cer.pen535$coefficients[1]*(log(t)-gaf.cer.pen535$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.pen535.surv(seq(0,24,0.1)), lty=2, col=3)
    
    auc.cer.pen535=integrate(f = gaf.cer.pen535.surv, lower=0, upper=24)[1]$value  
    
#1299 ppb ***********************************************************************************
    gaf.cer.pen1299<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[3] & cer$chem == 'pendimethalin'] ~ 
                           cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[3] & cer$chem == 'pendimethalin'],
                         type = 'binomial', fct = LL2.2())
    
    gaf.cer.pen1299.surv = function(t){
      (1/(1+exp(gaf.cer.pen1299$coefficients[1]*(log(t)-gaf.cer.pen1299$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.pen1299.surv(seq(0,24,0.1)), lty=2, col=4)
    
    auc.cer.pen1299=integrate(f = gaf.cer.pen1299.surv, lower=0, upper=24)[1]$value  
    
#2148 ppb **********************************************************************************
    gaf.cer.pen2148<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[4] & cer$chem == 'pendimethalin'] ~ 
                           cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[4] & cer$chem == 'pendimethalin'],
                         type = 'binomial', fct = LL2.2())
    
    gaf.cer.pen2148.surv = function(t){
      (1/(1+exp(gaf.cer.pen2148$coefficients[1]*(log(t)-gaf.cer.pen2148$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.pen2148.surv(seq(0,24,0.1)), lty=2, col=5)
    
    auc.cer.pen2148=integrate(f = gaf.cer.pen2148.surv, lower=0, upper=24)[1]$value     
    
#3762 ppb **********************************************************************************
    gaf.cer.pen3762<-drm(cer$surv[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[5] & cer$chem == 'pendimethalin'] ~ 
                           cer$time_hrs[cer$conc==unique(cer$conc[cer$chem == 'pendimethalin'])[5] & cer$chem == 'pendimethalin'],
                         type = 'binomial', fct = LL2.2())
    
    gaf.cer.pen3762.surv = function(t){
      (1/(1+exp(gaf.cer.pen3762$coefficients[1]*(log(t)-gaf.cer.pen3762$coefficients[2]))))
    } 
    
    lines(x=seq(0,24,0.1), y=gaf.cer.pen3762.surv(seq(0,24,0.1)), lty=2, col=6)
    
    auc.cer.pen3762=integrate(f = gaf.cer.pen3762.surv, lower=0, upper=24)[1]$value   
    
    title('Ghaffar2016 pendimethalin toxicity to cercariae')
    legend('topright', legend = c('control', 215,535,1299,2148,3762), pch = c(16,17,17,17,17,17), col = c(1:6), cex=0.8)  
    
    gaf.pen.aucs.cer = c(auc.cer.cont, auc.cer.pen215, auc.cer.pen535, auc.cer.pen1299, auc.cer.pen2148, auc.cer.pen3762)/auc.cer.cont
#Fit pi_C functions ################  
    cergaf.df = data.frame('but.conc'=but.vals,
                           'gly.conc'=gly.vals,
                           'pen.conc'=pen.vals,
                           'but.auc'=gaf.but.aucs.cer,
                           'gly.auc'=gaf.gly.aucs.cer,
                           'pen.auc'=gaf.pen.aucs.cer)
    
    plot(cergaf.df$gly.conc, cergaf.df$gly.auc, pch = 16, col = 'firebrick', xlim = c(0,27000), ylim = c(0,1),
         xlab = 'Herbicide Concentration', ylab = 'relative auc')
    points(cergaf.df$but.conc, cergaf.df$but.auc, pch = 16, col = 'darkorange')
    points(cergaf.df$pen.conc, cergaf.df$pen.auc, pch = 16, col = 'indianred2')
#Glyphosate   
  gly.piC.pred = nls(gly.auc ~ exp(-b*gly.conc), data=cergaf.df, start = list(b=0.0001))
    summary(gly.piC.pred)
    
  pi_C_gly_gaf16 = function(He){
      exp(-summary(gly.piC.pred)$parameters[1]*He) 
  }
  
  gly.piC.pred.l2 = drm(cergaf.df$gly.auc  ~ cergaf.df$gly.conc, type = 'binomial', fct = LL2.2())
  
  pi_C_gly_gaf16.l2 = function(He){
      (1/(1+exp(gly.piC.pred.l2$coefficients[1]*(log(He)-gly.piC.pred.l2$coefficients[2]))))
    }
    
    lines(c(0:27000), pi_C_gly_gaf16(c(0:27000)), lty=2, col='firebrick')  
    lines(c(0:27000), pi_C_gly_gaf16.l2(c(0:27000)), lty=3, col='firebrick')  
    
#butralin   
  but.piC.pred = nls(but.auc ~ exp(-b*but.conc), data=cergaf.df, start = list(b=0.001))
    summary(but.piC.pred)
    
  pi_C_but_gaf16 = function(He){
      exp(-summary(but.piC.pred)$parameters[1]*He) 
  }
  
  but.piC.pred.l2 = drm(cergaf.df$but.auc  ~ cergaf.df$but.conc, type = 'binomial', fct = LL2.2())
  
  pi_C_but_gaf16.l2 = function(He){
    (1/(1+exp(but.piC.pred.l2$coefficients[1]*(log(He)-but.piC.pred.l2$coefficients[2]))))
  }
  
    lines(c(0:27000), pi_C_but_gaf16(c(0:27000)), lty=2, col='darkorange')  
    lines(c(0:27000), pi_C_but_gaf16.l2(c(0:27000)), lty=3, col='darkorange')   
#pendimethalin   
  pen.piC.pred = nls(pen.auc ~ exp(-b*pen.conc), data=cergaf.df, start = list(b=0.005))
    summary(pen.piC.pred)
    
  pi_C_pen_gaf16 = function(He){
      exp(-summary(pen.piC.pred)$parameters[1]*He) 
  }
  
  pen.piC.pred.l2 = drm(cergaf.df$pen.auc  ~ cergaf.df$pen.conc, type = 'binomial', fct = LL2.2())
  
  pi_C_pen_gaf16.l2 = function(He){
    (1/(1+exp(pen.piC.pred.l2$coefficients[1]*(log(He)-pen.piC.pred.l2$coefficients[2]))))
  }
  
    lines(c(0:27000), pi_C_pen_gaf16(c(0:27000)), lty=2, col='indianred2')  
    lines(c(0:27000), pi_C_pen_gaf16.l2(c(0:27000)), lty=3, col='indianred2')  