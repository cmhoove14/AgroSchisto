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

#Cercarial mortality (S. mansoni) from Tchounwou 1992 ##################
cerc = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/tchounwou92.csv')
  cerc.t = subset(cerc, chem == 'malathion')
  cerc.t$dead = round(cerc.t$dead)
time = seq(0,25,0.1)

#Tchounwou Data plotted ############
  plot(cerc.t$time_hrs[cerc.t$conc==0], cerc.t$surv[cerc.t$conc==0]/100, pch=17, xlab = 'time(hrs)',
       ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    for(i in 2:length(unique(cerc.t$conc[cerc.t$chem == 'malathion']))){
      points(cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[i]], 
             cerc.t$surv[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[i]]/100, pch=16,
             col = i)
    }
  legend('topright', title = 'Mal(ppm)', legend = c(0,50,100,150,200,250), pch = c(17,rep(16,4)),
         col = c(1:6), cex=0.7)
  
#Fit models to survival curves and Add lines of best fit to plot ######################################
  #Control ********************************************************************************************
  cerc.ctrl = subset(cerc.t, conc == 0)
  ctrl.mod = drm(alive/total ~ time_hrs, total, data = cerc.ctrl, 
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
  
#50 ppm ********************************************************************************************  
  cerc.mal50 = subset(cerc.t, conc == 50000)
  mal50.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal50, 
                 type = 'binomial', fct = LL2.2())
  
    lines(time, ll4(1,0,mal50.mod$coefficients[1], mal50.mod$coefficients[2], time), lty=2, col=2)
    
    fx.mal50 = function(t){
      (1/(1+exp(mal50.mod$coefficients[1]*(log(t)-mal50.mod$coefficients[2]))))
    } 
    
    auc.mal50=integrate(f = fx.mal50, lower=0, upper=24)[1]$value
    
    mal50.aucs<-as.numeric()
    
    for(i in 1:10000){
      s.use = rnorm(1, summary(mal50.mod)$coefficients[1,1], summary(mal50.mod)$coefficients[1,2])
      lc.use = rnorm(1, summary(mal50.mod)$coefficients[2,1], summary(mal50.mod)$coefficients[2,2])
      
      mal50.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                                0, 24)[1]$value
    }
    
#100ppm ********************************************************************************************
  cerc.mal100 = subset(cerc.t, conc == 100000)
  mal100.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal100, 
                  type = 'binomial', fct = LL2.2())
  
    lines(time, ll4(1,0,mal100.mod$coefficients[1], mal100.mod$coefficients[2], time), lty=2, col=3)
    
    fx.mal100 = function(t){
      (1/(1+exp(mal100.mod$coefficients[1]*(log(t)-mal100.mod$coefficients[2]))))
    } 
    
    auc.mal100=integrate(f = fx.mal100, lower=0, upper=24)[1]$value
    
  mal100.aucs<-as.numeric()
    
    for(i in 1:10000){
      s.use = rnorm(1, summary(mal100.mod)$coefficients[1,1], summary(mal100.mod)$coefficients[1,2])
      lc.use = rnorm(1, summary(mal100.mod)$coefficients[2,1], summary(mal100.mod)$coefficients[2,2])
      
      mal100.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                                 0, 24)[1]$value
    }  
  
#150 ppm ********************************************************************************************
  cerc.mal150 = subset(cerc.t, conc == 150000)
  mal150.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal150, 
                   type = 'binomial', fct = LL2.2())
  
    lines(time, ll4(1,0,mal150.mod$coefficients[1], mal150.mod$coefficients[2], time), lty=2, col=4)
    
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
  
#200 ppm ********************************************************************************************
  cerc.mal200 = subset(cerc.t, conc == 200000)
  mal200.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal200, 
                   type = 'binomial', fct = LL2.2())
  
    lines(time, ll4(1,0,mal200.mod$coefficients[1], mal200.mod$coefficients[2], time), lty=2, col=5)
    
    fx.mal200 = function(t){
      (1/(1+exp(mal200.mod$coefficients[1]*(log(t)-mal200.mod$coefficients[2]))))
    } 
    
    auc.mal200=integrate(f = fx.mal200, lower=0, upper=24)[1]$value
    
  mal200.aucs<-as.numeric()
    
    for(i in 1:10000){
      s.use = rnorm(1, summary(mal200.mod)$coefficients[1,1], summary(mal200.mod)$coefficients[1,2])
      lc.use = rnorm(1, summary(mal200.mod)$coefficients[2,1], summary(mal200.mod)$coefficients[2,2])
      
      mal200.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                                 0, 24)[1]$value
    }  
  
#250ppm ********************************************************************************************
  cerc.mal250 = subset(cerc.t, conc == 250000)
  mal250.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal250, 
                   type = 'binomial', fct = LL2.2())
  
    lines(time, ll4(1,0,mal250.mod$coefficients[1], mal250.mod$coefficients[2], time), lty=2, col=6)
  
    fx.mal250 = function(t){
      (1/(1+exp(mal250.mod$coefficients[1]*(log(t)-mal250.mod$coefficients[2]))))
    } 
    
    auc.mal250=integrate(f = fx.mal250, lower=0, upper=24)[1]$value
    
  mal250.aucs<-as.numeric()

    for(i in 1:10000){
      s.use = rnorm(1, summary(mal250.mod)$coefficients[1,1], summary(mal250.mod)$coefficients[1,2])
      lc.use = rnorm(1, summary(mal250.mod)$coefficients[2,1], summary(mal250.mod)$coefficients[2,2])
      
      mal250.aucs[i] = integrate(function(t) {1/(1+exp(s.use*(log(t)-lc.use)))},
                                 0, 24, rel.tol = 24)[1]$value #relative tolerance added to prevent divergent integeral error
    }
  
#Plot AUC estimates and fit d-r function to the data ####################################################
  tch.df = data.frame(mal = c(0,50,100,150,200,250)*1000,
                      cent.auc = c(median(ctrl.aucs), median(mal50.aucs), median(mal100.aucs),
                                   median(mal150.aucs), median(mal200.aucs), median(mal250.aucs)),
                      var.auc = c(var(ctrl.aucs), var(mal50.aucs), var(mal100.aucs),
                                  var(mal150.aucs), var(mal200.aucs), var(mal250.aucs)),
                      pt.auc = c(auc.ctrl, auc.mal50, auc.mal100, auc.mal150, auc.mal200, auc.mal250),
                      drc.pred = 0,
                      drc.se = 0,
                      nls.pred = 0)
  
  plot(tch.df$mal, tch.df$cent.auc, pch = 16, cex = 1.2, xlab = 'Malathion (ppb)', ylab = 'AUC', ylim = c(0,20))
    for(i in 1:length(unique(tch.df$mal))){
      segments(x0 = tch.df$mal[i], y0 = (tch.df$cent.auc[i] - sqrt(tch.df$var.auc[i])),
               x1 = tch.df$mal[i], y1 = (tch.df$cent.auc[i] + sqrt(tch.df$var.auc[i])))
    }
  points(tch.df$mal, tch.df$pt.auc, pch = 17, col=2, cex=0.9)
  
  piC.tch.mal.init = drm(cent.auc ~ mal, data = tch.df, weights = var.auc^-1, type = 'continuous',
                         fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                                    fixed = c(NA, 0, median(ctrl.aucs), NA)))
  
  piC.tch.mal = nls(cent.auc ~ (median(ctrl.aucs)/(1+exp(slp*(log(mal)-log(lc))))), data = tch.df, weights = var.auc^-1, 
                    start = list(slp = piC.tch.mal.init$coefficients[1], lc = piC.tch.mal.init$coefficients[2]))
    
  tch.df[,5:6] = predict(piC.tch.mal.init, newdata = tch.df, se.fit = TRUE)
    points(tch.df$mal, tch.df$drc.pred, pch = 17, col = 3, cex = 0.9)
  
  tch.df[,7] = predict(piC.tch.mal, newdata = tch.df, se.fit = TRUE)
    points(tch.df$mal, tch.df$nls.pred, pch = 17, col = 4, cex = 0.9)

  legend('bottomleft', legend = c('median +/- SD', 'observed', 'drc-predict', 'nls-predict'), 
         pch = c(16,17,17,17), col = c(1:4), cex = 0.7)

#Check what the function looks like and generate function to use in subsequent models ##########################
  new.df = data.frame(mal = c(0:250000),
                      pred = 0,
                      pred.lo = 0,
                      pred.hi = 0)
  
  new.df[,2:4] = predict(piC.tch.mal.init, newdata = new.df, interval = 'confidence', level = 0.95)
  
  lines(new.df$mal, new.df$pred, lty=2, col=2)
  lines(new.df$mal, new.df$pred.lo, lty=3, col=2)
  lines(new.df$mal, new.df$pred.hi, lty=3, col=2)
 
#Create new plot for pi_C estimates which are scaled from 0 to 1
  plot(new.df$mal, new.df$pred / new.df$pred[1], type = 'l', lty = 2, xlab = 'Malathion concentration', ylab = expression(pi[C]), ylim = c(0,1.5))
    lines(new.df$mal, new.df$pred.hi / new.df$pred[1], lty = 3)
    lines(new.df$mal, new.df$pred.lo / new.df$pred[1], lty = 3)
    
#Function to return point estimate   
  pi_C_malathion_tch92 = function(In){
   predict(piC.tch.mal.init, newdata = data.frame(mal = In)) / predict(piC.tch.mal.init, newdata = data.frame(mal = 0))
  }
  #Check it's returning the right value
    lines(c(0:250000), pi_C_malathion_tch92(c(0:250000)), lty = 2, col = 2)

#Function for sampling uncertainty space
  pi_C_malathion_tch92_uncertainty = function(In){
    c1 = predict(piC.tch.mal.init, newdata = data.frame(mal = 0))
    
    ts = predict(piC.tch.mal.init, newdata = data.frame(mal = In), se.fit = TRUE)

    piC = rnorm(1,ts[1],ts[2]) / c1
    
    piC
  }

#Character vector of objects to keep from this script ###############
  keep.tch92.piC = c('piC.tch.mal.init','pi_C_malathion_tch92_uncertainty')