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

#miracidial mortality (S. mansoni) from Tchounwou 1991 ####################################
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
#Fit LL.2 function to each conentration time series  ####################################
mir.ctrl = subset(mir.mal, conc == 0)
ctrl.mod = drm(alive/total ~ time_hrs, total, data = mir.ctrl, 
               type = 'binomial', fct = LL.2())

  lines(time, ll4(1,0,ctrl.mod$coefficients[1], ctrl.mod$coefficients[2], time), lty=2)
  
  mir.fx0 = function(t){
    1/(1+exp(summary(ctrl.mod)$coefficients[1]*(log(t/summary(ctrl.mod)$coefficients[2]))))
  }
  
  auc.mir0 = integrate(mir.fx0, lower = 0, upper = 24)$value

#Malathion = 30000 ppb ********************************************************************************************  
mir.mal30 = subset(mir.mal, conc == 30000)
mal30.mod = drm(alive/total ~ time_hrs, total, data = mir.mal30, 
                type = 'binomial', fct = LL.2())

  lines(time, ll4(1,0,mal30.mod$coefficients[1], mal30.mod$coefficients[2], time), lty=2, col=2)
  
  mir.fx30 = function(t){
    1/(1+exp(summary(mal30.mod)$coefficients[1]*(log(t/summary(mal30.mod)$coefficients[2]))))
  }
  
  auc.mir30 = integrate(mir.fx30, lower = 0, upper = 24)$value

#Malathion = 60000 ppb ********************************************************************************************  
mir.mal60 = subset(mir.mal, conc == 60000)
mal60.mod = drm(alive/total ~ time_hrs, total, data = mir.mal60, 
                type = 'binomial', fct = LL.2())

  lines(time, ll4(1,0,mal60.mod$coefficients[1], mal60.mod$coefficients[2], time), lty=2, col=3)
  
  mir.fx60 = function(t){
    1/(1+exp(summary(mal60.mod)$coefficients[1]*(log(t/summary(mal60.mod)$coefficients[2]))))
  }
  
  auc.mir60 = integrate(mir.fx60, lower = 0, upper = 24)$value

#Malathion = 90000 ********************************************************************************************  
mir.mal90 = subset(mir.mal, conc == 90000)
mal90.mod = drm(alive/total ~ time_hrs, total, data = mir.mal90, 
                type = 'binomial', fct = LL.2())

  lines(time, ll4(1,0,mal90.mod$coefficients[1], mal90.mod$coefficients[2], time), lty=2, col=4)
  
  mir.fx90 = function(t){
    1/(1+exp(summary(mal90.mod)$coefficients[1]*(log(t/summary(mal90.mod)$coefficients[2]))))
  }
  
  auc.mir90 = integrate(mir.fx90, lower = 0, upper = 24)$value

#Malathion = 120000 ********************************************************************************************  
mir.mal120 = subset(mir.mal, conc == 120000)
mal120.mod = drm(alive/total ~ time_hrs, total, data = mir.mal120, 
                 type = 'binomial', fct = LL.2())

  lines(time, ll4(1,0,mal120.mod$coefficients[1], mal120.mod$coefficients[2], time), lty=2, col=5)
  
  mir.fx120 = function(t){
    1/(1+exp(summary(mal120.mod)$coefficients[1]*(log(t/summary(mal120.mod)$coefficients[2]))))
  }
  
  auc.mir120 = integrate(mir.fx120, lower = 0, upper = 24)$value

#Malathion = 150000 ********************************************************************************************  
mir.mal150 = subset(mir.mal, conc == 150000)
mal150.mod = drm(alive/total ~ time_hrs, total, data = mir.mal150, 
                 type = 'binomial', fct = LL.2())

  lines(time, ll4(1,0,mal150.mod$coefficients[1], mal150.mod$coefficients[2], time), lty=2, col=6)
  
  mir.fx150 = function(t){
    1/(1+exp(summary(mal150.mod)$coefficients[1]*(log(t/summary(mal150.mod)$coefficients[2]))))
  }
  
  auc.mir150 = integrate(mir.fx150, lower = 0, upper = 24)$value
  
#Data frame of AUC vals ######## ************************************************************************  
  mir.auc.mal = data.frame(mal = unique(mir$conc[mir$chem == 'mal']),
                           auc = c(auc.mir0, auc.mir30, auc.mir60, 
                                   auc.mir90, auc.mir120, auc.mir150),
                           piM = c(auc.mir0, auc.mir30, auc.mir60, 
                                   auc.mir90, auc.mir120, auc.mir150) / auc.mir0)

#Create data frame with parameter values and malathion concentrations #######################
mirp.df = data.frame(mal = c(0,30,60,90,120,150),
                      e = c(coef(ctrl.mod)[2], coef(mal30.mod)[2], coef(mal60.mod)[2],
                            coef(mal90.mod)[2], coef(mal120.mod)[2], coef(mal150.mod)[2]),
                      e.se = c(summary(ctrl.mod)$coefficients[2,2], summary(mal30.mod)$coefficients[2,2],
                               summary(mal60.mod)$coefficients[2,2], summary(mal90.mod)$coefficients[2,2],
                               summary(mal120.mod)$coefficients[2,2], summary(mal150.mod)$coefficients[2,2]),
                      b = c(coef(ctrl.mod)[1], coef(mal30.mod)[1], coef(mal60.mod)[1],
                            coef(mal90.mod)[1], coef(mal120.mod)[1], coef(mal150.mod)[1]),
                      b.se = c(summary(ctrl.mod)$coefficients[1,2], summary(mal30.mod)$coefficients[1,2],
                               summary(mal60.mod)$coefficients[1,2], summary(mal90.mod)$coefficients[1,2],
                               summary(mal120.mod)$coefficients[1,2], summary(mal150.mod)$coefficients[1,2]))

plot(mirp.df$mal, mirp.df$e, pch = 16, xlab = 'malathion (ppb)', ylab = 'LL.2 Parameters',
     ylim = c(0, 10))
points(mirp.df$mal, mirp.df$b, pch = 17, col=2)
  for(i in 1:length(mirp.df$mal)){
    segments(x0 = mirp.df$mal[i], y0 = mirp.df$e[i] + mirp.df$e.se[i],
             x1 = mirp.df$mal[i], y1 = mirp.df$e[i] - mirp.df$e.se[i])
    segments(x0 = mirp.df$mal[i], y0 = mirp.df$b[i] + mirp.df$b.se[i],
             x1 = mirp.df$mal[i], y1 = mirp.df$b[i] - mirp.df$b.se[i], col=2)
  }

  em.mod = lm(e ~ mal, weights = e.se^-1, data = mirp.df) 
  bm.mod = lm(b ~ mal, weights = b.se^-1, data = mirp.df) 

mod2df= data.frame(mal = c(0:150),
                   pred.e = 0,
                   pred.e.se = 0,
                   pred.b = 0,
                   pred.b.se = 0)

mod2df[,2:3] = predict(em.mod, newdata = mod2df, se.fit = TRUE)[1:2]
mod2df[,4:5] = predict(bm.mod, newdata = mod2df, se.fit = TRUE)[1:2]

lines(mod2df$mal, mod2df$pred.e, lty = 2)
  lines(mod2df$mal, mod2df$pred.e + 1.96*mod2df$pred.e.se, lty = 3)
  lines(mod2df$mal, mod2df$pred.e - 1.96*mod2df$pred.e.se, lty = 3)

lines(mod2df$mal, mod2df$pred.b, lty = 2, col=2)
  lines(mod2df$mal, mod2df$pred.b + 1.96*mod2df$pred.b.se, lty = 3, col=2)
  lines(mod2df$mal, mod2df$pred.b - 1.96*mod2df$pred.b.se, lty = 3, col=2)

#Create function to generate d-r function #####################
predm.fx = function(In){
  e = as.numeric(predict(em.mod, newdata = data.frame(mal = In/1000), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bm.mod, newdata = data.frame(mal = In/1000), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    if(e.use <= 0) e.use = mal150.mod$coefficients[2]
  b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t/e.use)))))}, 
                  lower=0, upper=24)[1]$value
  auc
}  
  
#Final:generate relative miracidia-hrs function ***NOTE:Doesn't work well at very high concentrations as***
piM.tch91_mal_unc = function(In){
  piM = predm.fx(In) / predm.fx(0)
  if(piM > 1) piM = 1
  else(return(piM))
}

keep.tch91.beq = c('mir.auc.mal', 'piM.tch91_mal_unc', 'mal150.mod', 'predm.fx', 'em.mod', 'bm.mod')

#Qualitative model validation ###############
#Regenerate plot of observed data
plot(mir.mal$time_hrs[mir.mal$conc==0], mir.mal$surv[mir.mal$conc==0], pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 2:length(unique(mir.mal$conc))){
    points(mir.mal$time_hrs[mir.mal$conc==unique(mir.mal$conc)[i]], 
           mir.mal$surv[mir.mal$conc==unique(mir.mal$conc)[i]], pch=16,
           col = i)
  }
  legend('topright', legend = c(0, 30, 60, 90, 120, 150), title = 'Malathion (ppm)',
         pch = c(17,rep(16,5)), col = c(1,2:6), cex = 0.8)

#function to plot model predictions
predm.fx.plot = function(In, clr){
  e = as.numeric(predict(em.mod, newdata = data.frame(mal = In/1000), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bm.mod, newdata = data.frame(mal = In/1000), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    if(e.use <= 0) e.use = mal150.mod$coefficients[2]
  b.use = rnorm(1, b[1], b[2])
  
  lines(time, ll4(1,0, b.use, e.use, time), lty=2, col = clr)
}   

#plot model predictions
for(i in c(0,30,60,90,120,150)*1000){
  c = i/30000 + 1
  print(c)
  replicate(10, predm.fx.plot(In = i, clr = c))
}

#plot model output compared to observed points
plot(seq(0, 150000, 1001)/1000, sapply(seq(0, 150000, 1001), piM.tch91_mal_unc), pch = 17, cex = 0.5,
     xlab = 'Malathion (ppm)', ylab = expression(paste(pi[M])),
     main = 'Sample Output of miracidial mortality function')
  points(mir.auc.mal$mal/1000, mir.auc.mal$piM, pch = 16, col=2)