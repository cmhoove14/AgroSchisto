#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(drc)

ll4 = function(hi,lo,slp,lc,x){
  lo + ((hi-lo)/(1+exp(slp*(log(x)-lc))))
}

cerc.g = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/Griggs2008.csv')
  time = c(0:25)

#Plot control survival curve, estimate time-dependent die-off function, and get auc from 0-24 hours ########
cerc.g0 = subset(cerc.g, chem == 'solvent')
  plot(x = cerc.g0$time_hrs,y = cerc.g0$surv/100,
       xlab = 'time (hrs)', ylab = 'prop alive', pch = 16, xlim = c(0,25), ylim = c(0,1))

  grg.ctrl = drm(alive/total ~ time_hrs, total, data = cerc.g0, 
                 type = 'binomial', fct = LL2.2())
  
  grg.cerc.ctrl.surv = function(t){
    (1/(1+exp(grg.ctrl$coefficients[1]*(log(t)-grg.ctrl$coefficients[2]))))
  } 
    
  lines(time, grg.cerc.ctrl.surv(time), lty=2)

  auc.grg.ctrl=integrate(f = grg.cerc.ctrl.surv, lower=0, upper=24)[1]$value
  
  ctrl.grg.aucs<-as.numeric()
  
  for(i in 1:10000){
    s.use = rnorm(1, summary(grg.ctrl)$coefficients[1,1], summary(grg.ctrl)$coefficients[1,2])
    lc.use = rnorm(1, summary(grg.ctrl)$coefficients[2,1], summary(grg.ctrl)$coefficients[2,2])
    
    ctrl.grg.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                                 0, 24)[1]$value
  }
  
#Function for low dose treatment (atrazine 15 ppb & metolachlor 10 ppb combined)  ##########
  cerc.g15 = subset(cerc.g, conc == 15)
  points(cerc.g15$time_hrs, cerc.g15$surv/100, pch = 17)
  
  grg.15 = drm(alive/total ~ time_hrs, total, data = cerc.g15, 
               type = 'binomial', fct = LL2.2())
  
  grg.cerc15 = function(t){
    (1/(1+exp(grg.15$coefficients[1]*(log(t)-grg.15$coefficients[2]))))
  } 
  
  lines(x=time, y=grg.cerc15(time), lty=3)
  
  auc.grg.atr15=integrate(f = grg.cerc15, lower=0, upper=24)[1]$value
  
  grg15.aucs<-as.numeric()
  
  for(i in 1:10000){
    s.use = rnorm(1, summary(grg.15)$coefficients[1,1], summary(grg.15)$coefficients[1,2])
    lc.use = rnorm(1, summary(grg.15)$coefficients[2,1], summary(grg.15)$coefficients[2,2])
    
    grg15.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                                 0, 24)[1]$value
  }
  
#Function for low dose treatment (atrazine 100 ppb & metolachlor 85 ppb combined)   ###############
cerc.g100 = subset(cerc.g, conc == 100)
  points(cerc.g100$time_hrs, cerc.g100$surv/100, pch = 15)
  
  grg.100 = drm(alive/total ~ time_hrs, total, data = cerc.g100, 
                type = 'binomial', fct = LL2.2())
  
  grg.cerc100 = function(t){
    (1/(1+exp(grg.100$coefficients[1]*(log(t)-grg.100$coefficients[2]))))
  } 
  
  lines(x=time, y=grg.cerc100(time), lty=4)
  
  auc.grg.atr100=integrate(f = grg.cerc100, lower=0, upper=24)[1]$value
  
  grg100.aucs<-as.numeric()
  
  for(i in 1:10000){
    s.use = rnorm(1, summary(grg.100)$coefficients[1,1], summary(grg.100)$coefficients[1,2])
    lc.use = rnorm(1, summary(grg.100)$coefficients[2,1], summary(grg.100)$coefficients[2,2])
    
    grg100.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                              0, 24)[1]$value
  }  
  title(main='Griggs2008 Atrazine-Cercarial mortality (E.trivolvis)')
    legend('bottomleft', legend = c('control', '15ppb', '100ppb'), pch = c(16,17,15), cex=0.8)
#Derive functional response of pi_C to atrazine concentration #############
grg.df = data.frame(atr = c(0,15,100),
                    cent.auc = c(mean(ctrl.grg.aucs), mean(grg15.aucs), mean(grg100.aucs)),
                    var.auc = c(var(ctrl.grg.aucs), var(grg15.aucs), var(grg100.aucs)),
                    pt.auc = c(auc.grg.ctrl, auc.grg.atr15, auc.grg.atr100),
                    drc.pred = 0,
                    drc.se = 0)  
    
  plot(grg.df$atr, grg.df$cent.auc, pch = 16, xlab = 'Atrazine (ppb)', ylab = 'AUC', ylim = c(0,20))
    for(i in 1:length(unique(grg.df$atr))){
      segments(x0 = grg.df$atr[i], y0 = (grg.df$cent.auc[i] - sqrt(grg.df$var.auc[i])),
               x1 = grg.df$atr[i], y1 = (grg.df$cent.auc[i] + sqrt(grg.df$var.auc[i])))
    }
    
    points(grg.df$atr, grg.df$pt.auc, pch = 17, col=2, cex=0.9)
    
  piC.grg.atr = drm(cent.auc ~ atr, data = grg.df, weights = var.auc^-1, type = 'continuous',
                    fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                               fixed = c(NA, 0, mean(ctrl.grg.aucs), NA)))
#Function doesn't converge, likely because of lack of data... Have to revisit later.    
  grg.df[,5:6] = predict(piC.grg.atr, newdata = grg.df, se.fit = TRUE)
    points(grg.df$atr, grg.df$drc.pred, pch = 17, col = 3, cex = 0.9)
    
    legend('bottomleft', legend = c('mean +/- SD', 'observed', 'drc-predict'), 
           pch = c(16,17,17,17), col = c(1:4), cex = 0.7)
    
#Check what the function looks like and generate function to use in subsequent models ##########################
new.df = data.frame(atr = c(0:200),
                    pred = 0,
                    pred.lo = 0,
                    pred.hi = 0)
    
  new.df[,2:4] = predict(piC.grg.atr, newdata = new.df, interval = 'confidence', level = 0.95)
    
    lines(new.df$atr, new.df$pred, lty=2, col=2)
    lines(new.df$atr, new.df$pred.lo, lty=3, col=2)
    lines(new.df$atr, new.df$pred.hi, lty=3, col=2)