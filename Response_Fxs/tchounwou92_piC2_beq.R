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
  lo + ((hi-lo)/(1+exp(slp*(log(x/lc)))))
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
                 type = 'binomial', fct = LL.2())

  lines(time, ll4(1,0,ctrl.mod$coefficients[1], ctrl.mod$coefficients[2], time), lty=2)
  
  cerc.fx0 = function(t){
    1/(1+exp(summary(ctrl.mod)$coefficients[1]*(log(t/summary(ctrl.mod)$coefficients[2]))))
  }
  
  auc.cerc0 = integrate(cerc.fx0, lower = 0, upper = 24)$value
  
#50 ppm ********************************************************************************************  
  cerc.mal50 = subset(cerc.t, conc == 50000)
  mal50.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal50, 
                  type = 'binomial', fct = LL.2())
  
  lines(time, ll4(1,0,mal50.mod$coefficients[1], mal50.mod$coefficients[2], time), lty=2, col=2)
  
  cerc.fx50 = function(t){
    1/(1+exp(summary(mal50.mod)$coefficients[1]*(log(t/summary(mal50.mod)$coefficients[2]))))
  }
  
  auc.cerc50 = integrate(cerc.fx50, lower = 0, upper = 24)$value
  
#100ppm ********************************************************************************************
  cerc.mal100 = subset(cerc.t, conc == 100000)
  mal100.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal100, 
                   type = 'binomial', fct = LL.2())
  
  lines(time, ll4(1,0,mal100.mod$coefficients[1], mal100.mod$coefficients[2], time), lty=2, col=3)
  
  cerc.fx100 = function(t){
    1/(1+exp(summary(mal100.mod)$coefficients[1]*(log(t/summary(mal100.mod)$coefficients[2]))))
  }
  
  auc.cerc100 = integrate(cerc.fx100, lower = 0, upper = 24)$value
  
  
#150 ppm ********************************************************************************************
  cerc.mal150 = subset(cerc.t, conc == 150000)
  mal150.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal150, 
                   type = 'binomial', fct = LL.2())
  
  lines(time, ll4(1,0,mal150.mod$coefficients[1], mal150.mod$coefficients[2], time), lty=2, col=4)
  
  cerc.fx150 = function(t){
    1/(1+exp(summary(mal150.mod)$coefficients[1]*(log(t/summary(mal150.mod)$coefficients[2]))))
  }
  
  auc.cerc150 = integrate(cerc.fx150, lower = 0, upper = 24)$value
  
#200 ppm ********************************************************************************************
  cerc.mal200 = subset(cerc.t, conc == 200000)
  mal200.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal200, 
                   type = 'binomial', fct = LL.2())
  
  lines(time, ll4(1,0,mal200.mod$coefficients[1], mal200.mod$coefficients[2], time), lty=2, col=5)
  
  cerc.fx200 = function(t){
    1/(1+exp(summary(mal200.mod)$coefficients[1]*(log(t/summary(mal200.mod)$coefficients[2]))))
  }
  
  auc.cerc200 = integrate(cerc.fx200, lower = 0, upper = 24)$value
  
#250ppm ********************************************************************************************
  cerc.mal250 = subset(cerc.t, conc == 250000)
  mal250.mod = drm(alive/total ~ time_hrs, total, data = cerc.mal250, 
                   type = 'binomial', fct = LL.2())
  
  lines(time, ll4(1,0,mal250.mod$coefficients[1], mal250.mod$coefficients[2], time), lty=2, col=6)
  
  cerc.fx250 = function(t){
    1/(1+exp(summary(mal250.mod)$coefficients[1]*(log(t/summary(mal250.mod)$coefficients[2]))))
  }
  
  auc.cerc250 = integrate(cerc.fx250, lower = 0, upper = 24)$value

#Store relative auc values as observed points ###############    
  cerc.auc.mal = data.frame(mal = unique(cerc$conc[cerc$chem == 'malathion']),
                            auc = c(auc.cerc0, auc.cerc50, auc.cerc100, 
                                    auc.cerc150, auc.cerc200, auc.cerc250),
                            piC = c(auc.cerc0, auc.cerc50, auc.cerc100, 
                                    auc.cerc150, auc.cerc200, auc.cerc250) / auc.cerc0)

#Data frame of LL.2 parameters across maltahion concentrations ##################    
parms.df = data.frame(mal = c(0,50,100,150,200,250),
                      logmal = log(c(0,50,100,150,200,250)+1),
                      e = c(coef(ctrl.mod)[2], coef(mal50.mod)[2], coef(mal100.mod)[2],
                            coef(mal150.mod)[2], coef(mal200.mod)[2], coef(mal250.mod)[2]),
                      e.se = c(summary(ctrl.mod)$coefficients[2,2], summary(mal50.mod)$coefficients[2,2],
                               summary(mal100.mod)$coefficients[2,2], summary(mal150.mod)$coefficients[2,2],
                               summary(mal200.mod)$coefficients[2,2], summary(mal250.mod)$coefficients[2,2]),
                      b = c(coef(ctrl.mod)[1], coef(mal50.mod)[1], coef(mal100.mod)[1],
                            coef(mal150.mod)[1], coef(mal200.mod)[1], coef(mal250.mod)[1]),
                      b.se = c(summary(ctrl.mod)$coefficients[1,2], summary(mal50.mod)$coefficients[1,2],
                               summary(mal100.mod)$coefficients[1,2], summary(mal150.mod)$coefficients[1,2],
                               summary(mal200.mod)$coefficients[1,2], summary(mal250.mod)$coefficients[1,2]))

  plot(parms.df$mal, parms.df$e, pch = 16, xlab = 'malathion (ppm)', ylab = 'LL.2 Parameters',
       ylim = c(0, 15))
    points(parms.df$mal, parms.df$b, pch = 17, col=2)
    for(i in 1:length(parms.df$mal)){
      segments(x0 = parms.df$mal[i], y0 = parms.df$e[i] + parms.df$e.se[i],
               x1 = parms.df$mal[i], y1 = parms.df$e[i] - parms.df$e.se[i])
      segments(x0 = parms.df$mal[i], y0 = parms.df$b[i] + parms.df$b.se[i],
               x1 = parms.df$mal[i], y1 = parms.df$b[i] - parms.df$b.se[i], col=2)
    }
    legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
           cex = 0.8, bty = 'n', title = 'Observed points')
    
  e.mod = lm(e ~ mal, weights = e.se^-1, data = parms.df) 
  e.mod2 = lm(e ~ logmal, weights = e.se^-1, data = parms.df)
    AIC(e.mod, e.mod2)
  
  #Log-linear model fits better and asymptotes at 0 so we'll use it  
  b.mod = lm(b ~ mal, weights = b.se^-1, data = parms.df) 


  mod.df= data.frame(mal = c(0:250),
                     logmal = log(c(0:250)+1),
                     pred.e = 0,
                     pred.e.se = 0,
                     pred.e2 = 0,
                     pred.e2.se = 0,
                     pred.b = 0,
                     pred.b.se = 0)
  
  mod.df[,3:4] = predict(e.mod, newdata = mod.df, se.fit = TRUE)[1:2]
  mod.df[,5:6] = predict(e.mod2, newdata = mod.df, se.fit = TRUE)[1:2]
  mod.df[,7:8] = predict(b.mod, newdata = mod.df, se.fit = TRUE)[1:2]
    
  lines(mod.df$mal, mod.df$pred.e, lty = 2)
    lines(mod.df$mal, mod.df$pred.e + 1.96*mod.df$pred.e.se, lty = 3)
    lines(mod.df$mal, mod.df$pred.e - 1.96*mod.df$pred.e.se, lty = 3)
    
  lines(mod.df$mal, mod.df$pred.e2, lty = 2, col=4)
    lines(mod.df$mal, mod.df$pred.e2 + 1.96*mod.df$pred.e2.se, lty = 3, col=4)
    lines(mod.df$mal, mod.df$pred.e2 - 1.96*mod.df$pred.e2.se, lty = 3, col=4)  
    
  lines(mod.df$mal, mod.df$pred.b, lty = 2, col=2)
    lines(mod.df$mal, mod.df$pred.b + 1.96*mod.df$pred.b.se, lty = 3, col=2)
    lines(mod.df$mal, mod.df$pred.b - 1.96*mod.df$pred.b.se, lty = 3, col=2)
    
  legend('top', lty = c(2,2,2,3), col = c(1,4,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  

#Function to estimate survival curve as function of malathion conc #####################################
pred.fx = function(In){
  e = as.numeric(predict(e.mod2, newdata = data.frame(logmal = log(In/1000+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(b.mod, newdata = data.frame(mal = In/1000), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    if(e.use <= 0) e.use = mal250.mod$coefficients[2]
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
                  lower=0, upper=24)[1]$value
  auc
}  

#Final function and keep vector ###############
piC.tch92_mal_unc = function(In){
  if(In == 0) piC = 1
  else(piC = pred.fx(In) / pred.fx(0))
  if(piC > 1) piC = 1
  else(return(piC))
}  

keep.tch92.beq = c('cerc.auc.mal', 'piC.tch92_mal_unc', 'mal250.mod', 'pred.fx', 'e.mod2', 'b.mod')  

#Qualitative model validation ###############
#Regenerate plot of observed data
plot(cerc.t$time_hrs[cerc.t$conc==0], cerc.t$surv[cerc.t$conc==0]/100, pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24), 
     main = expression(paste(pi[C], ' Model sim-val')))
  for(i in 2:length(unique(cerc.t$conc[cerc.t$chem == 'malathion']))){
    points(cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[i]], 
           cerc.t$surv[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[i]]/100, pch=16,
           col = i)
  }
legend('topright', title = 'Mal(ppm)', legend = c(0,50,100,150,200,250), pch = c(17,rep(16,4)),
       col = c(1:6), cex=0.7)

#function to plot model predictions
pred.fx.plot = function(In, clr){
  e = as.numeric(predict(e.mod2, newdata = data.frame(logmal = log(In/1000+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(b.mod, newdata = data.frame(mal = In/1000), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    if(e.use <= 0) e.use = mal250.mod$coefficients[2]
  b.use = rnorm(1, b[1], b[2])
  
  lines(time, ll4(1,0, b.use, e.use, time), lty=2, col = clr)
}   

#plot model predictions
for(i in c(0,50,100,150,200,250)*1000){
  c = i/50000 + 1
  print(c)
    replicate(10, pred.fx.plot(In = i, clr = c))
  }
  
#plot model output compared to observed points
set.seed = 0
plot(seq(0,250000, 1000), sapply(seq(0,250000, 1000), piC.tch92_mal_unc), pch = 17, cex = 0.5,
     xlab = 'Malathion (ppb)', ylab = expression(paste(pi[C], 'estimate', sep = ' ')))
  points(cerc.auc.mal$mal, cerc.auc.mal$piC, pch = 16, col=2)