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

#cercarial toxicity ############
  tant<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/Tantawy2002.csv')
  cerc<-subset(tant, larv == 'cercariae')
  
#butachlor toxicity to cercariae  ###############
  plot(cerc$time_hrs[cerc$conc==0 & cerc$chem == 'butachlor'], 
       cerc$surv[cerc$conc==0 & cerc$chem == 'butachlor']/100, 
       pch=17, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 2:length(unique(cerc$conc[cerc$chem == 'butachlor']))){
    points(cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[i] & cerc$chem == 'butachlor'], 
           cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[i] & cerc$chem == 'butachlor']/100, pch=16,
           col = i)
  }
  
#fit to control points
  tant.cerc.but0<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[1] & cerc$chem == 'butachlor']/100 ~ 
                      cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[1] & cerc$chem == 'butachlor'],
                      type = 'binomial', fct = LL.2())
  
  tant.cerc.but0.surv = function(t){
    (1/(1+exp(tant.cerc.but0$coefficients[1]*(log(t/tant.cerc.but0$coefficients[2])))))
  } 
  
    lines(x=c(0:24), y=tant.cerc.but0.surv(c(0:24)), lty = 2)
    
    auc.cerc.but0=integrate(f = tant.cerc.but0.surv, lower=0, upper=24)[1]$value  
    
    
#650 ppb
    tant.cerc.but650<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[2] & cerc$chem == 'butachlor']/100 ~ 
                          cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[2] & cerc$chem == 'butachlor'],
                        type = 'binomial', fct = LL.2())
    
    tant.cerc.but650.surv = function(t){
      (1/(1+exp(tant.cerc.but650$coefficients[1]*(log(t/tant.cerc.but650$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.but650.surv(c(0:24)), lty=2, col=2)
    
    auc.cerc.but650=integrate(f = tant.cerc.but650.surv, lower=0, upper=24)[1]$value  
    
    
#1500 ppb
    tant.cerc.but1500<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[3] & cerc$chem == 'butachlor']/100 ~ 
                            cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[3] & cerc$chem == 'butachlor'],
                          type = 'binomial', fct = LL.2())
    
    tant.cerc.but1500.surv = function(t){
      (1/(1+exp(tant.cerc.but1500$coefficients[1]*(log(t/tant.cerc.but1500$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.but1500.surv(c(0:24)), lty=2, col=3)
    
    auc.cerc.but1500=integrate(f = tant.cerc.but1500.surv, lower=0, upper=24)[1]$value  
    
#4500 ppb
    tant.cerc.but4500<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[4] & cerc$chem == 'butachlor']/100 ~ 
                             cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[4] & cerc$chem == 'butachlor'],
                           type = 'binomial', fct = LL.2())
    
    tant.cerc.but4500.surv = function(t){
      (1/(1+exp(tant.cerc.but4500$coefficients[1]*(log(t/tant.cerc.but4500$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.but4500.surv(c(0:24)), lty=2, col=4)
    
    auc.cerc.but4500=integrate(f = tant.cerc.but4500.surv, lower=0, upper=24)[1]$value  
    
#6500 ppb
    tant.cerc.but6500<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[5] & cerc$chem == 'butachlor']/100 ~ 
                             cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[5] & cerc$chem == 'butachlor'],
                           type = 'binomial', fct = LL.2())
    
    tant.cerc.but6500.surv = function(t){
      (1/(1+exp(tant.cerc.but6500$coefficients[1]*(log(t/tant.cerc.but6500$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.but6500.surv(c(0:24)), lty=2, col=5)
    title('butachlor toxicity to cercariae')
    legend('topright', legend = c('control', 650,1500,4500,6500), pch = c(17,16,16,16,16), col = c(1:5), cex=0.8)
    
    auc.cerc.but6500=integrate(f = tant.cerc.but6500.surv, lower=0, upper=24)[1]$value  

#compile butachlor data for function ###############
cerc.but = data.frame(e = c(coef(tant.cerc.but0)[2], coef(tant.cerc.but650)[2], coef(tant.cerc.but1500)[2],
                            coef(tant.cerc.but4500)[2], coef(tant.cerc.but6500)[2]),
                      e.se = c(summary(tant.cerc.but0)$coefficients[2,2], 
                               summary(tant.cerc.but650)$coefficients[2,2],
                               summary(tant.cerc.but1500)$coefficients[2,2], 
                               summary(tant.cerc.but4500)$coefficients[2,2],
                               summary(tant.cerc.but6500)$coefficients[2,2]),
                      b = c(coef(tant.cerc.but0)[1], coef(tant.cerc.but650)[1], coef(tant.cerc.but1500)[1],
                            coef(tant.cerc.but4500)[1], coef(tant.cerc.but6500)[1]),
                      b.se = c(summary(tant.cerc.but0)$coefficients[1,2], 
                               summary(tant.cerc.but650)$coefficients[1,2],
                               summary(tant.cerc.but1500)$coefficients[1,2], 
                               summary(tant.cerc.but4500)$coefficients[1,2],
                               summary(tant.cerc.but6500)$coefficients[1,2]),
                      but = c(0,650,1500,4500,6500),
                      logbut = log(c(0,650,1500,4500,6500)+1))
    
plot(cerc.but$but, cerc.but$e, pch = 16, xlab = 'butachlor (ppm)', ylab = 'LL.2 Parameters',
     ylim = c(0,20))

  points(cerc.but$but+50, cerc.but$b, pch = 17, col=2)
  
    for(i in 1:length(cerc.but$but)){
      segments(x0 = cerc.but$but[i], y0 = cerc.but$e[i] + cerc.but$e.se[i],
               x1 = cerc.but$but[i], y1 = cerc.but$e[i] - cerc.but$e.se[i])
      segments(x0 = cerc.but$but[i]+50, y0 = cerc.but$b[i] + cerc.but$b.se[i],
               x1 = cerc.but$but[i]+50, y1 = cerc.but$b[i] - cerc.but$b.se[i], col=2)
    }    
  
#fit models to LL.2 parameters across concentration ########
  el.but = lm(e ~ but, weights = e.se^-1, data = cerc.but) #linear response of LC50
    el.pred = function(but){
      predict(el.but, newdata = data.frame(but = but), 
              interval = 'confidence', level = 0.95)
    }
    
    lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred, simplify = T)[3,], lty = 3)

  el.but2 = lm(e ~ logbut, weights = e.se^-1, data = cerc.but) #log-linear response of LC50
    el.pred2 = function(but){
      predict(el.but2, newdata = data.frame(logbut = log(but)+1), 
              interval = 'confidence', level = 0.95)
    }
    
  lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred2, simplify = T)[1,], lty = 2, col=3)
  lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred2, simplify = T)[2,], lty = 3, col=3)
  lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred2, simplify = T)[3,], lty = 3, col=3)
    
  AIC(el.but, el.but2)  #Linear is a better fit    
  
  bl.but = lm(b ~ but, weights = b.se^-1, data = cerc.but)   
    bl.pred = function(but){
      predict(bl.but, newdata = data.frame(but = but), interval = 'confidence', level = 0.95)
    }
  
    lines(seq(0,7000,7), sapply(seq(0,7000,7), bl.pred, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,7000,7), sapply(seq(0,7000,7), bl.pred, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,7000,7), sapply(seq(0,7000,7), bl.pred, simplify = T)[3,], lty = 3, col = 2)
    legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
           legend = c('LC50 - Linear',
                      'LC50 - Exponential',
                      'Slp - linear',
                      '95% CI'), cex = 0.7)  
    legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
           cex = 0.8, bty = 'n', title = 'Observed points')
      
but.fx.lin = function(He){
  e = as.numeric(predict(el.but, newdata = data.frame(but = He), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bl.but, newdata = data.frame(but = He), se.fit = TRUE)[1:2])
      
  e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
      
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
                      lower=0, upper=24)[1]$value
      auc
}  #function to estimate AUC 
 
piC.tant02_but.lin_unc = function(He){
  if(He == 0) piC = 1
  else(piC = but.fx.lin(He) / but.fx.lin(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

but.fx.exp = function(He){
  e = as.numeric(predict(el.but2, newdata = data.frame(logbut = log(He+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bl.but, newdata = data.frame(but = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
                  lower=0, upper=24)[1]$value
  auc
}  #function to estimate AUC 

piC.tant02_but.exp_unc = function(He){
  if(He == 0) piC = 1
  else(piC = but.fx.exp(He) / but.fx.exp(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  #Parameter estimate

#plot sample output compared to observed AUCs ##########
  plot(cerc.but$but, c(auc.cerc.but0, auc.cerc.but650, auc.cerc.but1500, 
         auc.cerc.but4500, auc.cerc.but6500)/auc.cerc.but0,
        xlab = 'Butachlor (ppb)', ylab = 'relative AUC (cercariae-hours)',
       pch = 16, ylim = c(0,1))
    points(seq(0, 6500,10), sapply(seq(0, 6500,10), piC.tant02_but.lin_unc, simplify = T),
           pch = 5, col = 4, cex = 0.5)
    points(seq(0, 6500,10), sapply(seq(0, 6500,10), piC.tant02_but.exp_unc, simplify = T),
           pch = 5, col = 2, cex = 0.5)
     
#fluazifop-p-butyl toxicity to cercariae ###########
    plot(cerc$time_hrs[cerc$conc==0 & cerc$chem == 'fluazifop-p-butyl'], 
         cerc$surv[cerc$conc==0 & cerc$chem == 'fluazifop-p-butyl']/100, 
         pch=17, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    for(i in 2:length(unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl']))){
      points(cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[i] & cerc$chem == 'fluazifop-p-butyl'], 
             cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[i] & cerc$chem == 'fluazifop-p-butyl']/100, pch=16,
             col = i)
    }
    
  #fit to control points
    tant.cerc.fpb0<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[1] & cerc$chem == 'fluazifop-p-butyl']/100 ~ 
                          cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[1] & cerc$chem == 'fluazifop-p-butyl'],
                        type = 'binomial', fct = LL.2())
    
    tant.cerc.fpb0.surv = function(t){
      (1/(1+exp(tant.cerc.fpb0$coefficients[1]*(log(t/tant.cerc.fpb0$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.fpb0.surv(c(0:24)), lty = 2)
    
    auc.cerc.fpb0=integrate(f = tant.cerc.fpb0.surv, lower=0, upper=24)[1]$value  
    
    
  #1760 ppb
    tant.cerc.fpb1760<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[2] & cerc$chem == 'fluazifop-p-butyl']/100 ~ 
                            cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[2] & cerc$chem == 'fluazifop-p-butyl'],
                          type = 'binomial', fct = LL.2())
    
    tant.cerc.fpb1760.surv = function(t){
      (1/(1+exp(tant.cerc.fpb1760$coefficients[1]*(log(t/tant.cerc.fpb1760$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.fpb1760.surv(c(0:24)), lty=2, col=2)
    
    auc.cerc.fpb1760=integrate(f = tant.cerc.fpb1760.surv, lower=0, upper=24)[1]$value  
    
    
  #4500 ppb
    tant.cerc.fpb4500<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[3] & cerc$chem == 'fluazifop-p-butyl']/100 ~ 
                             cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[3] & cerc$chem == 'fluazifop-p-butyl'],
                           type = 'binomial', fct = LL.2())
    
    tant.cerc.fpb4500.surv = function(t){
      (1/(1+exp(tant.cerc.fpb4500$coefficients[1]*(log(t/tant.cerc.fpb4500$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.fpb4500.surv(c(0:24)), lty=2, col=3)
    
    auc.cerc.fpb4500=integrate(f = tant.cerc.fpb4500.surv, lower=0, upper=24)[1]$value  
    
  #9000 ppb
    tant.cerc.fpb9000<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[4] & cerc$chem == 'fluazifop-p-butyl']/100 ~ 
                             cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[4] & cerc$chem == 'fluazifop-p-butyl'],
                           type = 'binomial', fct = LL.2())
    
    tant.cerc.fpb9000.surv = function(t){
      (1/(1+exp(tant.cerc.fpb9000$coefficients[1]*(log(t/tant.cerc.fpb9000$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.fpb9000.surv(c(0:24)), lty=2, col=4)
    
    auc.cerc.fpb9000=integrate(f = tant.cerc.fpb9000.surv, lower=0, upper=24)[1]$value  
    
  #17600 ppb
    tant.cerc.fpb17600<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[5] & cerc$chem == 'fluazifop-p-butyl']/100 ~ 
                             cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[5] & cerc$chem == 'fluazifop-p-butyl'],
                           type = 'binomial', fct = LL.2())
    
    tant.cerc.fpb17600.surv = function(t){
      (1/(1+exp(tant.cerc.fpb17600$coefficients[1]*(log(t/tant.cerc.fpb17600$coefficients[2])))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.fpb17600.surv(c(0:24)), lty=2, col=5)
    title('fluazifop-p-butyl toxicity to cercariae')
    legend('topright', legend = c('control', 1760,4500,9000,17600), pch = c(17,16,16,16,16), col = c(1:5), cex=0.8)
    
    auc.cerc.fpb17600=integrate(f = tant.cerc.fpb17600.surv, lower=0, upper=24)[1]$value  
    
#compile fluazifop-p-butyl data for function ###############
cerc.fpb = data.frame(e = c(coef(tant.cerc.fpb0)[2], coef(tant.cerc.fpb1760)[2], coef(tant.cerc.fpb4500)[2],
                                coef(tant.cerc.fpb9000)[2], coef(tant.cerc.fpb17600)[2]),
                      e.se = c(summary(tant.cerc.fpb0)$coefficients[2,2], 
                                   summary(tant.cerc.fpb1760)$coefficients[2,2],
                                   summary(tant.cerc.fpb4500)$coefficients[2,2], 
                                   summary(tant.cerc.fpb9000)$coefficients[2,2],
                                   summary(tant.cerc.fpb17600)$coefficients[2,2]),
                      b = c(coef(tant.cerc.fpb0)[1], coef(tant.cerc.fpb1760)[1], coef(tant.cerc.fpb4500)[1],
                                coef(tant.cerc.fpb9000)[1], coef(tant.cerc.fpb17600)[1]),
                      b.se = c(summary(tant.cerc.fpb0)$coefficients[1,2], 
                                   summary(tant.cerc.fpb1760)$coefficients[1,2],
                                   summary(tant.cerc.fpb4500)$coefficients[1,2], 
                                   summary(tant.cerc.fpb9000)$coefficients[1,2],
                                   summary(tant.cerc.fpb17600)$coefficients[1,2]),
                      fpb = c(0,1760,4500,9000,17600),
                      logfpb = log(c(0,1760,4500,9000,17600)+1))
    
    plot(cerc.fpb$fpb, cerc.fpb$e, pch = 16, xlab = 'Fluazifop-p-butyl (ppb)', ylab = 'LL.2 Parameters',
         ylim = c(0,20))
    
    points(cerc.fpb$fpb+100, cerc.fpb$b, pch = 17, col=2)
    
    for(i in 1:length(cerc.fpb$fpb)){
      segments(x0 = cerc.fpb$fpb[i], y0 = cerc.fpb$e[i] + cerc.fpb$e.se[i],
               x1 = cerc.fpb$fpb[i], y1 = cerc.fpb$e[i] - cerc.fpb$e.se[i])
      segments(x0 = cerc.fpb$fpb[i]+100, y0 = cerc.fpb$b[i] + cerc.fpb$b.se[i],
               x1 = cerc.fpb$fpb[i]+100, y1 = cerc.fpb$b[i] - cerc.fpb$b.se[i], col=2)
    }    
    
#fit models to LL.2 parameters across concentration ########
  el.fpb = lm(e ~ fpb, weights = e.se^-1, data = cerc.fpb) #linear response of LC50
    el.pred.fpb = function(fpb){
      predict(el.fpb, newdata = data.frame(fpb = fpb), 
              interval = 'confidence', level = 0.95)
    }
    
    lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred.fpb, simplify = T)[1,], lty = 2)
    lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred.fpb, simplify = T)[2,], lty = 3)
    lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred.fpb, simplify = T)[3,], lty = 3)
    
  el.fpb2 = lm(e ~ logfpb, weights = e.se^-1, data = cerc.fpb) #log-linear response of LC50
    el.pred.fpb2 = function(fpb){
      predict(el.fpb2, newdata = data.frame(logfpb = log(fpb)+1), 
              interval = 'confidence', level = 0.95)
    }
    
    lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred.fpb2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred.fpb2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred.fpb2, simplify = T)[3,], lty = 3, col=3)
    
      AIC(el.fpb, el.fpb2)  #exponential is a better fit    
    
  bl.fpb = lm(b ~ fpb, weights = b.se^-1, data = cerc.fpb)   
    bl.pred.fpb = function(fpb){
      predict(bl.fpb, newdata = data.frame(fpb = fpb), interval = 'confidence', level = 0.95)
    }
    
    lines(seq(0,18000,100), sapply(seq(0,18000,100), bl.pred.fpb, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,18000,100), sapply(seq(0,18000,100), bl.pred.fpb, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,18000,100), sapply(seq(0,18000,100), bl.pred.fpb, simplify = T)[3,], lty = 3, col = 2)
    
    legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
           legend = c('LC50 - Linear',
                      'LC50 - Exponential',
                      'Slp - linear',
                      '95% CI'), cex = 0.7)  
    legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
           cex = 0.8, bty = 'n', title = 'Observed points')
    
    fpb.fx.exp = function(He){
      e = as.numeric(predict(el.fpb2, newdata = data.frame(logfpb = log(He+1)), se.fit = TRUE)[1:2])
      b = as.numeric(predict(bl.fpb, newdata = data.frame(fpb = He), se.fit = TRUE)[1:2])
      
      e.use = rnorm(1, e[1], e[2])
      b.use = rnorm(1, b[1], b[2])
      
      auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
                      lower=0, upper=24)[1]$value
      auc
    }  #function to estimate AUC with exponential fit to LC50
    
    piC.tant02_fpb.exp_unc = function(He){
      if(He == 0) piC = 1
      else(piC = fpb.fx.exp(He) / fpb.fx.exp(0))
      if(piC > 1) piC = 1
      else(return(piC))
    }  #Parameter estimate with exponential function
    
    fpb.fx.lin = function(He){
      e = as.numeric(predict(el.fpb, newdata = data.frame(fpb = He), se.fit = TRUE)[1:2])
      b = as.numeric(predict(bl.fpb, newdata = data.frame(fpb = He), se.fit = TRUE)[1:2])
      
      e.use = rnorm(1, e[1], e[2])
      b.use = rnorm(1, b[1], b[2])
      
      auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
                      lower=0, upper=24)[1]$value
      auc
    }  #function to estimate AUC with linear fit to LC50
    
    piC.tant02_fpb.lin_unc = function(He){
      if(He == 0) piC = 1
      else(piC = fpb.fx.lin(He) / fpb.fx.lin(0))
      if(piC > 1) piC = 1
      else(return(piC))
    }  #Parameter estimate with linear function
    
#plot sample output compared to observed AUCs ##########
  plot(cerc.fpb$fpb, c(auc.cerc.fpb0, auc.cerc.fpb1760, auc.cerc.fpb4500, 
                        auc.cerc.fpb9000, auc.cerc.fpb17600)/auc.cerc.fpb0,
        xlab = 'Fluazifop-p-butyl (ppb)', ylab = 'relative AUC (cercariae-hours)',
        pch = 16, ylim = c(0,1))
    points(seq(0, 18000,100), sapply(seq(0, 18000,100), piC.tant02_fpb.lin_unc, simplify = T),
           pch = 5, col = 4, cex = 0.5)
    points(seq(0, 18000,100), sapply(seq(0, 18000,100), piC.tant02_fpb.exp_unc, simplify = T),
           pch = 5, col = 2, cex = 0.5)
    
    #looks like linear response performs better at low concentrations, exponential performs better
      #at high concentrations; so function used will be determined by EEC / concentration range tested
#keep vector ########    
keep.tantawy.piC = c('bl.but', 'el.but', 'el.but2', 'bl.fpb', 'el.fpb', 'el.fpb2',
                     'but.fx.exp', 'but.fx.lin', 'fpb.fx.exp', 'fpb.fx.lin',
                     'piC.tant02_but.exp_unc', 'piC.tant02_but.lin_unc',
                     'piC.tant02_fpb.exp_unc', 'piC.tant02_fpb.lin_unc', 
                     'auc.cerc.fpb0', 'auc.cerc.fpb1760', 'auc.cerc.fpb4500', 
                     'auc.cerc.fpb9000', 'auc.cerc.fpb17600', 'auc.cerc.but0', 'auc.cerc.but650', 
                     'auc.cerc.but1500', 'auc.cerc.but4500', 'auc.cerc.but6500')    
    