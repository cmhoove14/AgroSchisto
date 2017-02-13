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

#mirarial mortality (S. mansoni) from Tchounwou 1992 ##################
mir = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Miracidial Mortality/tchounwou08_miracidia.csv')
  mir.mal = subset(mir, chem == 'mal')
  time = seq(0,25,0.1)

#Tchounwou Data plotted ############
plot(mir.mal$time_hrs[mir.mal$conc==0], mir.mal$surv[mir.mal$conc==0], pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 2:length(unique(mir.mal$conc))){
    points(mir.mal$time_hrs[mir.mal$conc==unique(mir.mal$conc)[i]], 
           mir.mal$surv[mir.mal$conc==unique(mir.mal$conc)[i]], pch=16,
           col = i)
  }
  legend('topright', legend = c(0, 30, 60, 90, 120, 150), title = 'Malathion (ppm)',
         pch = c(17,rep(16,5)), col = c(1,2:6), cex = 0.8)
  
#Fit to tchounwou control ########
#Malathion experiment control points
  mir.ctrl = subset(mir.mal, conc == 0)
  ctrl.mod = drm(alive/total ~ time_hrs, total, data = mir.ctrl, 
                 type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,ctrl.mod$coefficients[1], ctrl.mod$coefficients[2], time), lty=2)
  
  fx.ctrl = function(t){
    (1/(1+exp(ctrl.mod$coefficients[1]*(log(t)-ctrl.mod$coefficients[2]))))
  } 
  
  auc.ctrl=integrate(f = fx.ctrl, lower=0, upper=24)[1]$value
  
  ctrl.aucs<-as.numeric()
  
  for(i in 1:10000){
    s.use = rnorm(1, summary(ctrl.mod)$coefficients[1,1], summary(ctrl.mod)$coefficients[1,2])
    lc.use = rnorm(1, summary(ctrl.mod)$coefficients[2,1], summary(ctrl.mod)$coefficients[2,2])
    
    ctrl.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                             0, 24)[1]$value
  }

#Fit log logistic model to data #######

#Malathion = 30000 ppb ********************************************************************************************  
  mir.mal30 = subset(mir.mal, conc == 30000)
  mal30.mod = drm(alive/total ~ time_hrs, total, data = mir.mal30, 
                  type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,mal30.mod$coefficients[1], mal30.mod$coefficients[2], time), lty=2, col=2)
  
  fx.mal30 = function(t){
    (1/(1+exp(mal30.mod$coefficients[1]*(log(t)-mal30.mod$coefficients[2]))))
  } 
  
  auc.mal30=integrate(f = fx.mal30, lower=0, upper=24)[1]$value
  
  mal30.aucs<-as.numeric()
  
  for(i in 1:10000){
    s.use = rnorm(1, summary(mal30.mod)$coefficients[1,1], summary(mal30.mod)$coefficients[1,2])
    lc.use = rnorm(1, summary(mal30.mod)$coefficients[2,1], summary(mal30.mod)$coefficients[2,2])
    
    mal30.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                              0, 24)[1]$value
  }

#Malathion = 60000 ppb ********************************************************************************************  
  mir.mal60 = subset(mir.mal, conc == 60000)
  mal60.mod = drm(alive/total ~ time_hrs, total, data = mir.mal60, 
                  type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,mal60.mod$coefficients[1], mal60.mod$coefficients[2], time), lty=2, col=3)
  
  fx.mal60 = function(t){
    (1/(1+exp(mal60.mod$coefficients[1]*(log(t)-mal60.mod$coefficients[2]))))
  } 
  
  auc.mal60=integrate(f = fx.mal60, lower=0, upper=24)[1]$value
  
  mal60.aucs<-as.numeric()
  
  for(i in 1:10000){
    s.use = rnorm(1, summary(mal60.mod)$coefficients[1,1], summary(mal60.mod)$coefficients[1,2])
    lc.use = rnorm(1, summary(mal60.mod)$coefficients[2,1], summary(mal60.mod)$coefficients[2,2])
    
    mal60.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                              0, 24)[1]$value
  }

#Malathion = 90000 ********************************************************************************************  
  mir.mal90 = subset(mir.mal, conc == 90000)
  mal90.mod = drm(alive/total ~ time_hrs, total, data = mir.mal90, 
                  type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,mal90.mod$coefficients[1], mal90.mod$coefficients[2], time), lty=2, col=4)
  
  fx.mal90 = function(t){
    (1/(1+exp(mal90.mod$coefficients[1]*(log(t)-mal90.mod$coefficients[2]))))
  } 
  
  auc.mal90=integrate(f = fx.mal90, lower=0, upper=24)[1]$value
  
  mal90.aucs<-as.numeric()
  
  for(i in 1:10000){
    s.use = rnorm(1, summary(mal90.mod)$coefficients[1,1], summary(mal90.mod)$coefficients[1,2])
    lc.use = rnorm(1, summary(mal90.mod)$coefficients[2,1], summary(mal90.mod)$coefficients[2,2])
    
    mal90.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                              0, 24)[1]$value
  }

#Malathion = 120000 ********************************************************************************************  
  mir.mal120 = subset(mir.mal, conc == 120000)
  mal120.mod = drm(alive/total ~ time_hrs, total, data = mir.mal120, 
                  type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,mal120.mod$coefficients[1], mal120.mod$coefficients[2], time), lty=2, col=5)
  
  fx.mal120 = function(t){
    (1/(1+exp(mal120.mod$coefficients[1]*(log(t)-mal120.mod$coefficients[2]))))
  } 
  
  auc.mal120=integrate(f = fx.mal120, lower=0, upper=24)[1]$value
  
  mal120.aucs<-as.numeric()
  
  for(i in 1:10000){
    s.use = rnorm(1, summary(mal120.mod)$coefficients[1,1], summary(mal120.mod)$coefficients[1,2])
    lc.use = rnorm(1, summary(mal120.mod)$coefficients[2,1], summary(mal120.mod)$coefficients[2,2])
    
    mal120.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                              0, 24)[1]$value
  }

#Malathion = 150000 ********************************************************************************************  
  mir.mal150 = subset(mir.mal, conc == 150000)
  mal150.mod = drm(alive/total ~ time_hrs, total, data = mir.mal150, 
                  type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,mal150.mod$coefficients[1], mal150.mod$coefficients[2], time), lty=2, col=6)
  
  fx.mal150 = function(t){
    (1/(1+exp(mal150.mod$coefficients[1]*(log(t)-mal150.mod$coefficients[2]))))
  } 
  
  auc.mal150=integrate(f = fx.mal150, lower=0, upper=24)[1]$value
  
  mal150.aucs<-as.numeric()
  
  for(i in 1:10000){
    s.use = rnorm(1, summary(mal150.mod)$coefficients[1,1], summary(mal150.mod)$coefficients[1,2])
    lc.use = rnorm(1, summary(mal150.mod)$coefficients[2,1], summary(mal150.mod)$coefficients[2,2])
    
    mal150.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                              0, 24)[1]$value
  }
  
#Derive functional response of pi_M to malathion concentration #############
#Check decrease in AUC across malathion concentration treating AUC as continuous variable
  mir.df = data.frame(mal = unique(mir.mal$conc),
                      cent.auc = c(median(ctrl.aucs), median(mal30.aucs), median(mal60.aucs),
                                   median(mal90.aucs), median(mal120.aucs), median(mal150.aucs)),
                      var.auc = c(var(ctrl.aucs), var(mal30.aucs), var(mal60.aucs),
                                  var(mal90.aucs), var(mal120.aucs), var(mal150.aucs)),
                      pt.auc = c(auc.ctrl, auc.mal30, auc.mal60, auc.mal90, auc.mal120, auc.mal150),
                      drc.pred = 0,
                      drc.se = 0)
  
  plot(mir.df$mal, mir.df$cent.auc, pch = 16, cex = 1.2, xlab = 'Malathion (ppb)', ylab = 'AUC', ylim = c(0,12))
    for(i in 1:length(unique(mir.df$mal))){
      segments(x0 = mir.df$mal[i], y0 = (mir.df$cent.auc[i] - sqrt(mir.df$var.auc[i])),
               x1 = mir.df$mal[i], y1 = (mir.df$cent.auc[i] + sqrt(mir.df$var.auc[i])))
    }
  points(mir.df$mal, mir.df$pt.auc, pch = 17, col=2, cex=0.9)
  
  piM.tch.mal.init = drm(cent.auc ~ mal, data = mir.df, weights = var.auc^-1, type = 'continuous',
                         fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                                    fixed = c(NA, 0, median(ctrl.aucs), NA)))
  
  mir.df[,5:6] = predict(piM.tch.mal.init, newdata = mir.df, se.fit = TRUE)
    points(mir.df$mal, mir.df$drc.pred, pch = 17, col = 3, cex = 0.9)
    
    legend('bottomleft', legend = c('median +/- SD', 'observed', 'drc-predict'), 
           pch = c(16,17,17), col = c(1:3), cex = 0.7)
    
  
#Check what the function looks like and generate function to use in subsequent models ##########################
  new.df = data.frame(mal = c(0:150000),
                      pred = 0,
                      pred.lo = 0,
                      pred.hi = 0)
    
  new.df[,2:4] = predict(piM.tch.mal.init, newdata = new.df, interval = 'confidence', level = 0.95)
    
    lines(new.df$mal, new.df$pred, lty=2, col=2)
    lines(new.df$mal, new.df$pred.lo, lty=3, col=2)
    lines(new.df$mal, new.df$pred.hi, lty=3, col=2)
    
  
#Create new plot for pi_C estimates which are scaled from 0 to 1
  plot(new.df$mal, new.df$pred / new.df$pred[1], type = 'l', lty = 2, xlab = 'Malathion concentration', ylab = expression(pi[M]), ylim = c(0,1.5))
    lines(new.df$mal, new.df$pred.hi / new.df$pred[1], lty = 3)
    lines(new.df$mal, new.df$pred.lo / new.df$pred[1], lty = 3)
    
#Function to return point estimate   
    pi_M_malathion_tch92 = function(In){
      predict(piM.tch.mal.init, newdata = data.frame(mal = In)) / predict(piM.tch.mal.init, newdata = data.frame(mal = 0))
    }
#Check it's returning the right value
  lines(c(0:150000), pi_M_malathion_tch92(c(0:150000)), lty = 2, col = 2)
    
#Function for sampling uncertainty space
    pi_M_malathion_tch92_uncertainty = function(In){
      c1 = predict(piM.tch.mal.init, newdata = data.frame(mal = 0))
      
      ts = predict(piM.tch.mal.init, newdata = data.frame(mal = In), se.fit = TRUE)

      piM = rnorm(1,ts[1],ts[2]) / c1
      
      piM
    }
    
    #for(i in seq(0, 150000, length.out = 10000)){
    #  points(i, pi_M_malathion_tch92_uncertainty(i), pch=17, col=4, cex=0.5)
    #}
    
#Character vector of objects to keep from this script ###############
  keep.tch92.piM = c('piM.tch.mal.init','pi_M_malathion_tch92_uncertainty')    