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

#miracidia toxicity ############
tant<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/Tantawy2002.csv')
  mir<-subset(tant, larv == 'miracidia')
    
#butachlor toxicity to miracidia  ###############
    plot(mir$time_hrs[mir$conc==0 & mir$chem == 'butachlor'], 
         mir$surv[mir$conc==0 & mir$chem == 'butachlor']/100, 
         pch=17, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    for(i in 2:length(unique(mir$conc[mir$chem == 'butachlor']))){
      points(mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[i] & mir$chem == 'butachlor'], 
             mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[i] & mir$chem == 'butachlor']/100, pch=16,
             col = i)
    }
    
#fit to control points
    tant.mir.but0<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[1] & mir$chem == 'butachlor']/100 ~ 
                          mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[1] & mir$chem == 'butachlor'],
                        type = 'binomial', fct = LL.2())
    
    tant.mir.but0.surv = function(t){
      (1/(1+exp(tant.mir.but0$coefficients[1]*(log(t/tant.mir.but0$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.mir.but0.surv(c(0:24)), lty = 2)
    
    auc.mir.but0=integrate(f = tant.mir.but0.surv, lower=0, upper=24)[1]$value  
    
    
#650 ppb
    tant.mir.but650<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[2] & mir$chem == 'butachlor']/100 ~ 
                            mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[2] & mir$chem == 'butachlor'],
                          type = 'binomial', fct = LL.2())
    
    tant.mir.but650.surv = function(t){
      (1/(1+exp(tant.mir.but650$coefficients[1]*(log(t/tant.mir.but650$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.mir.but650.surv(c(0:24)), lty=2, col=2)
    
    auc.mir.but650=integrate(f = tant.mir.but650.surv, lower=0, upper=24)[1]$value  
    
    
#1500 ppb
    tant.mir.but1500<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[3] & mir$chem == 'butachlor']/100 ~ 
                             mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[3] & mir$chem == 'butachlor'],
                           type = 'binomial', fct = LL.2())
    
    tant.mir.but1500.surv = function(t){
      (1/(1+exp(tant.mir.but1500$coefficients[1]*(log(t/tant.mir.but1500$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.mir.but1500.surv(c(0:24)), lty=2, col=3)
    
    auc.mir.but1500=integrate(f = tant.mir.but1500.surv, lower=0, upper=24)[1]$value  
    
#4500 ppb
    tant.mir.but4500<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[4] & mir$chem == 'butachlor']/100 ~ 
                             mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[4] & mir$chem == 'butachlor'],
                           type = 'binomial', fct = LL.2())
    
    tant.mir.but4500.surv = function(t){
      (1/(1+exp(tant.mir.but4500$coefficients[1]*(log(t/tant.mir.but4500$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.mir.but4500.surv(c(0:24)), lty=2, col=4)
    
    auc.mir.but4500=integrate(f = tant.mir.but4500.surv, lower=0, upper=24)[1]$value  
    
#6500 ppb
    tant.mir.but6500<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[5] & mir$chem == 'butachlor']/100 ~ 
                             mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[5] & mir$chem == 'butachlor'],
                           type = 'binomial', fct = LL.2())
    
    tant.mir.but6500.surv = function(t){
      (1/(1+exp(tant.mir.but6500$coefficients[1]*(log(t/tant.mir.but6500$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.mir.but6500.surv(c(0:24)), lty=2, col=5)
    title('butachlor toxicity to miracidia')
    legend('topright', legend = c('control', 650,1500,4500,6500), pch = c(17,16,16,16,16), col = c(1:5), cex=0.8)
    
    auc.mir.but6500=integrate(f = tant.mir.but6500.surv, lower=0, upper=24)[1]$value     
    
#Compile LL.2 data for functional responses #############
mir.but = data.frame(e = c(coef(tant.mir.but0)[2], coef(tant.mir.but650)[2], coef(tant.mir.but1500)[2],
                            coef(tant.mir.but4500)[2], coef(tant.mir.but6500)[2]),
                      e.se = c(summary(tant.mir.but0)$coefficients[2,2], 
                               summary(tant.mir.but650)$coefficients[2,2],
                               summary(tant.mir.but1500)$coefficients[2,2], 
                               summary(tant.mir.but4500)$coefficients[2,2],
                               summary(tant.mir.but6500)$coefficients[2,2]),
                      b = c(coef(tant.mir.but0)[1], coef(tant.mir.but650)[1], coef(tant.mir.but1500)[1],
                            coef(tant.mir.but4500)[1], coef(tant.mir.but6500)[1]),
                      b.se = c(summary(tant.mir.but0)$coefficients[1,2], 
                               summary(tant.mir.but650)$coefficients[1,2],
                               summary(tant.mir.but1500)$coefficients[1,2], 
                               summary(tant.mir.but4500)$coefficients[1,2],
                               summary(tant.mir.but6500)$coefficients[1,2]),
                      but = c(0,650,1500,4500,6500),
                      logbut = log(c(0,650,1500,4500,6500)+1))

plot(mir.but$but, mir.but$e, pch = 16, xlab = 'butachlor (ppm)', ylab = 'LL.2 Parameters (miracidia',
     ylim = c(0,20))

  points(mir.but$but+50, mir.but$b, pch = 17, col=2)

  for(i in 1:length(mir.but$but)){
    segments(x0 = mir.but$but[i], y0 = mir.but$e[i] + mir.but$e.se[i],
             x1 = mir.but$but[i], y1 = mir.but$e[i] - mir.but$e.se[i])
    segments(x0 = mir.but$but[i]+50, y0 = mir.but$b[i] + mir.but$b.se[i],
             x1 = mir.but$but[i]+50, y1 = mir.but$b[i] - mir.but$b.se[i], col=2)
  } 
  
#fit models to LL.2 parameters across concentration ########
el.but.mir = lm(e ~ but, weights = e.se^-1, data = mir.but) #linear response of LC50
  el.pred.mir = function(but){
    predict(el.but.mir, newdata = data.frame(but = but), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred.mir, simplify = T)[1,], lty = 2)
    lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred.mir, simplify = T)[2,], lty = 3)
    lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred.mir, simplify = T)[3,], lty = 3)

el.but2.mir = lm(e ~ logbut, weights = e.se^-1, data = mir.but) #log-linear response of LC50
  el.pred2.mir = function(but){
    predict(el.but2.mir, newdata = data.frame(logbut = log(but)+1), 
            interval = 'confidence', level = 0.95)
  }

  lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred2.mir, simplify = T)[1,], lty = 2, col=3)
  lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred2.mir, simplify = T)[2,], lty = 3, col=3)
  lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred2.mir, simplify = T)[3,], lty = 3, col=3)

AIC(el.but.mir, el.but2.mir)  #exponential is a better fit    

bl.but.mir = lm(b ~ but, weights = b.se^-1, data = mir.but)   
  bl.pred.mir = function(but){
    predict(bl.but.mir, newdata = data.frame(but = but), interval = 'confidence', level = 0.95)
  }

  lines(seq(0,7000,7), sapply(seq(0,7000,7), bl.pred.mir, simplify = T)[1,], lty = 2, col = 2)
  lines(seq(0,7000,7), sapply(seq(0,7000,7), bl.pred.mir, simplify = T)[2,], lty = 3, col = 2)
  lines(seq(0,7000,7), sapply(seq(0,7000,7), bl.pred.mir, simplify = T)[3,], lty = 3, col = 2)
  
  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

but.fx.lin.mir = function(He){
  e = as.numeric(predict(el.but.mir, newdata = data.frame(but = He), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bl.but.mir, newdata = data.frame(but = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piM.tant02_but.lin_unc = function(He){
  if(He == 0) piC = 1
  else(piC = but.fx.lin.mir(He) / but.fx.lin.mir(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

but.fx.exp.mir = function(He){
  e = as.numeric(predict(el.but2.mir, newdata = data.frame(logbut = log(He+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bl.but.mir, newdata = data.frame(but = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piM.tant02_but.exp_unc = function(He){
  if(He == 0) piC = 1
  else(piC = but.fx.exp.mir(He) / but.fx.exp.mir(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

#plot sample outputs compared to observed points ##########
plot(mir.but$but, c(auc.mir.but0, auc.mir.but650, auc.mir.but1500, 
                     auc.mir.but4500, auc.mir.but6500)/auc.mir.but0,
     xlab = 'Butachlor (ppb)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0, 6500,10), sapply(seq(0, 6500,10), piM.tant02_but.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0, 6500,10), sapply(seq(0, 6500,10), piM.tant02_but.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)
    
#fluazifop-p-butyl toxicity to cercariae ###########
    plot(mir$time_hrs[mir$conc==0 & mir$chem == 'fluazifop-p-butyl'], 
         mir$surv[mir$conc==0 & mir$chem == 'fluazifop-p-butyl']/100, 
         pch=17, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    for(i in 2:length(unique(mir$conc[mir$chem == 'fluazifop-p-butyl']))){
      points(mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[i] & mir$chem == 'fluazifop-p-butyl'], 
             mir$surv[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[i] & mir$chem == 'fluazifop-p-butyl']/100, pch=16,
             col = i)
    }
    
#fit to control points
    tant.mir.fpb0<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[1] & mir$chem == 'fluazifop-p-butyl']/100 ~ 
                          mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[1] & mir$chem == 'fluazifop-p-butyl'],
                        type = 'binomial', fct = LL.2())
    
    tant.mir.fpb0.surv = function(t){
      (1/(1+exp(tant.mir.fpb0$coefficients[1]*(log(t/tant.mir.fpb0$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.mir.fpb0.surv(c(0:24)))
    
    auc.mir.fpb0=integrate(f = tant.mir.fpb0.surv, lower=0, upper=24)[1]$value  
    
    
#1760 ppb
    tant.mir.fpb1760<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[2] & mir$chem == 'fluazifop-p-butyl']/100 ~ 
                             mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[2] & mir$chem == 'fluazifop-p-butyl'],
                           type = 'binomial', fct = LL.2())
    
    tant.mir.fpb1760.surv = function(t){
      (1/(1+exp(tant.mir.fpb1760$coefficients[1]*(log(t/tant.mir.fpb1760$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.mir.fpb1760.surv(c(0:24)), lty=2, col=2)
    
    auc.mir.fpb1760=integrate(f = tant.mir.fpb1760.surv, lower=0, upper=24)[1]$value  
    
    
#4500 ppb
    tant.mir.fpb4500<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[3] & mir$chem == 'fluazifop-p-butyl']/100 ~ 
                             mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[3] & mir$chem == 'fluazifop-p-butyl'],
                           type = 'binomial', fct = LL.2())
    
    tant.mir.fpb4500.surv = function(t){
      (1/(1+exp(tant.mir.fpb4500$coefficients[1]*(log(t/tant.mir.fpb4500$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.mir.fpb4500.surv(c(0:24)), lty=2, col=3)
    
    auc.mir.fpb4500=integrate(f = tant.mir.fpb4500.surv, lower=0, upper=24)[1]$value  
    
#9000 ppb
    tant.mir.fpb9000<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[4] & mir$chem == 'fluazifop-p-butyl']/100 ~ 
                             mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[4] & mir$chem == 'fluazifop-p-butyl'],
                           type = 'binomial', fct = LL.2())
    
    tant.mir.fpb9000.surv = function(t){
      (1/(1+exp(tant.mir.fpb9000$coefficients[1]*(log(t/tant.mir.fpb9000$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.mir.fpb9000.surv(c(0:24)), lty=2, col=4)
    
    auc.mir.fpb9000=integrate(f = tant.mir.fpb9000.surv, lower=0, upper=24)[1]$value  
    
#17600 ppb
    tant.mir.fpb17600<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[5] & mir$chem == 'fluazifop-p-butyl']/100 ~ 
                              mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[5] & mir$chem == 'fluazifop-p-butyl'],
                            type = 'binomial', fct = LL.2())
    
    tant.mir.fpb17600.surv = function(t){
      (1/(1+exp(tant.mir.fpb17600$coefficients[1]*(log(t/tant.mir.fpb17600$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.mir.fpb17600.surv(c(0:24)), lty=2, col=5)
    title('fluazifop-p-butyl toxicity to miracidia')
    legend('topright', legend = c('control', 1760,4500,9000,17600), pch = c(17,16,16,16,16), col = c(1:5), cex=0.8)
    
    auc.mir.fpb17600=integrate(f = tant.mir.fpb17600.surv, lower=0, upper=24)[1]$value  
    
    
#Compile LL.2 data for functional responses #############
mir.fpb = data.frame(e = c(coef(tant.mir.fpb0)[2], coef(tant.mir.fpb1760)[2], coef(tant.mir.fpb4500)[2],
                           coef(tant.mir.fpb9000)[2], coef(tant.mir.fpb17600)[2]),
                     e.se = c(summary(tant.mir.fpb0)$coefficients[2,2], 
                              summary(tant.mir.fpb1760)$coefficients[2,2],
                              summary(tant.mir.fpb4500)$coefficients[2,2], 
                              summary(tant.mir.fpb9000)$coefficients[2,2],
                              summary(tant.mir.fpb17600)$coefficients[2,2]),
                     b = c(coef(tant.mir.fpb0)[1], coef(tant.mir.fpb1760)[1], coef(tant.mir.fpb4500)[1],
                           coef(tant.mir.fpb9000)[1], coef(tant.mir.fpb17600)[1]),
                     b.se = c(summary(tant.mir.fpb0)$coefficients[1,2], 
                              summary(tant.mir.fpb1760)$coefficients[1,2],
                              summary(tant.mir.fpb4500)$coefficients[1,2], 
                              summary(tant.mir.fpb9000)$coefficients[1,2],
                              summary(tant.mir.fpb17600)$coefficients[1,2]),
                     fpb = c(0,1760,4500,9000,17600),
                     logfpb = log(c(0,1760,4500,9000,17600)+1))

  plot(mir.fpb$fpb, mir.fpb$e, pch = 16, xlab = 'butachlor (ppm)', ylab = 'LL.2 Parameters (miracidia',
       ylim = c(0,20))
  
  points(mir.fpb$fpb+150, mir.fpb$b, pch = 17, col=2)
  
  for(i in 1:length(mir.fpb$fpb)){
    segments(x0 = mir.fpb$fpb[i], y0 = mir.fpb$e[i] + mir.fpb$e.se[i],
             x1 = mir.fpb$fpb[i], y1 = mir.fpb$e[i] - mir.fpb$e.se[i])
    segments(x0 = mir.fpb$fpb[i]+150, y0 = mir.fpb$b[i] + mir.fpb$b.se[i],
             x1 = mir.fpb$fpb[i]+150, y1 = mir.fpb$b[i] - mir.fpb$b.se[i], col=2)
  } 

#fit models to LL.2 parameters across concentration ########
el.fpb.mir = lm(e ~ fpb, weights = e.se^-1, data = mir.fpb) #linear response of LC50
  el.pred.mir.fpb = function(fpb){
    predict(el.fpb.mir, newdata = data.frame(fpb = fpb), 
            interval = 'confidence', level = 0.95)
  }

  lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred.mir.fpb, simplify = T)[1,], lty = 2)
  lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred.mir.fpb, simplify = T)[2,], lty = 3)
  lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred.mir.fpb, simplify = T)[3,], lty = 3)

el.fpb2.mir = lm(e ~ logfpb, weights = e.se^-1, data = mir.fpb) #log-linear response of LC50
  el.pred2.mir.fpb = function(fpb){
    predict(el.fpb2.mir, newdata = data.frame(logfpb = log(fpb)+1), 
            interval = 'confidence', level = 0.95)
  }

  lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred2.mir.fpb, simplify = T)[1,], lty = 2, col=3)
  lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred2.mir.fpb, simplify = T)[2,], lty = 3, col=3)
  lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred2.mir.fpb, simplify = T)[3,], lty = 3, col=3)

    AIC(el.fpb.mir, el.fpb2.mir)  #Linear is a better fit    

bl.fpb.mir = lm(b ~ fpb, weights = b.se^-1, data = mir.fpb)   
  bl.pred.mir.fpb = function(fpb){
    predict(bl.fpb.mir, newdata = data.frame(fpb = fpb), interval = 'confidence', level = 0.95)
  }

  lines(seq(0,18000,100), sapply(seq(0,18000,100), bl.pred.mir.fpb, simplify = T)[1,], lty = 2, col = 2)
  lines(seq(0,18000,100), sapply(seq(0,18000,100), bl.pred.mir.fpb, simplify = T)[2,], lty = 3, col = 2)
  lines(seq(0,18000,100), sapply(seq(0,18000,100), bl.pred.mir.fpb, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

fpb.fx.lin.mir = function(He){
  e = as.numeric(predict(el.fpb.mir, newdata = data.frame(fpb = He), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bl.fpb.mir, newdata = data.frame(fpb = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piM.tant02_fpb.lin_unc = function(He){
  if(He == 0) piC = 1
  else(piC = fpb.fx.lin.mir(He) / fpb.fx.lin.mir(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

fpb.fx.exp.mir = function(He){
  e = as.numeric(predict(el.fpb2.mir, newdata = data.frame(logfpb = log(He+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bl.fpb.mir, newdata = data.frame(fpb = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piM.tant02_fpb.exp_unc = function(He){
  if(He == 0) piC = 1
  else(piC = fpb.fx.exp.mir(He) / fpb.fx.exp.mir(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

#plot sample output compared to observed AUCs ##########
plot(mir.fpb$fpb, c(auc.mir.fpb0, auc.mir.fpb1760, auc.mir.fpb4500, 
                     auc.mir.fpb9000, auc.mir.fpb17600)/auc.mir.fpb0,
     xlab = 'Fluazifop-p-butyl (ppb)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0, 18000,50), sapply(seq(0, 18000,50), piM.tant02_fpb.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0, 18000,50), sapply(seq(0, 18000,50), piM.tant02_fpb.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)

#keep vector ##########
  keep.tantawy.piM = c('bl.but.mir', 'el.but.mir', 'el.but2.mir', 'bl.fpb.mir', 'el.fpb.mir', 'el.fpb2.mir',
                       'but.fx.exp.mir', 'but.fx.lin.mir', 'fpb.fx.exp.mir', 'fpb.fx.lin.mir',
                       'piM.tant02_but.exp_unc', 'piM.tant02_but.lin_unc',
                       'piM.tant02_fpb.exp_unc', 'piM.tant02_fpb.lin_unc', 
                       'auc.mir.fpb0', 'auc.mir.fpb1760', 'auc.mir.fpb4500', 
                       'auc.mir.fpb9000', 'auc.mir.fpb17600', 'auc.mir.but0', 'auc.mir.but650', 
                       'auc.mir.but1500', 'auc.mir.but4500', 'auc.mir.but6500')    
  
  
