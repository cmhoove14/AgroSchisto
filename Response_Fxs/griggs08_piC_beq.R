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

L.3.fx = function(t, lc50 = lc50, slp = slp){
  1 / (1+exp(slp*(t - lc50)))
}

cerc.g = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/Griggs2008.csv')
  cerc.g$prop_surv = cerc.g$alive / cerc.g$total
  time = seq(0,25,0.1)
cerc.g0 = subset(cerc.g, chem != 'control')

#Plot control survival curve, estimate time-dependent die-off function, and get auc from 0-24 hours ########
grg.mod = drm(prop_surv ~ time_hrs, conc, weights = total, data = cerc.g0, type = 'binomial', 
                 fct = L.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))
  summary(grg.mod)
  plot(grg.mod)

  plot(x = cerc.g0$time_hrs[cerc.g0$conc == 0], y = cerc.g0$prop_surv[cerc.g0$conc == 0],
       xlab = 'time (hrs)', ylab = 'prop alive', pch = 16, cex = 1.2, xlim = c(0,25), ylim = c(0,1))
    lines(time, predict(grg.mod, 
                        data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in unique(cerc.g0$conc)[c(2:3)]){
    points(cerc.g0$time_hrs[cerc.g0$conc == i], cerc.g0$prop_surv[cerc.g0$conc == i], 
           pch = 17, col = i+1)
    lines(time, predict(grg.mod, 
                        data.frame(time_hrs=time, conc = i)), lty = 2, col = i+1)
  } 

  title(main='Griggs2008 Atrazine-Cercarial mortality (E.trivolvis)')
    legend('bottomleft', legend = c('control', '15ppb', '100ppb'), 
           pch = c(16,17,17), col = c(1,8,5), cex=0.8, bty = 'n')
    
#Get estimate of cercariae-hrs for each concentration    
grg08_atr_aucs = as.numeric()
    
  for(j in 1:length(unique(cerc.g0$conc))){
    fx = function(t){
      predict(grg.mod, newdata = data.frame(time_hrs = t, conc = unique(cerc.g0$conc)[j]))
    }
    grg08_atr_aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }

#Create data frame with parameter values and atrazine concentrations #######################
grgc.df = data.frame(atr = c(0,15,100),
                     logatr = log(c(0,15,100)+1),
                     e = c(summary(grg.mod)$coefficients[c(4:6),1]),
                     e.se = c(summary(grg.mod)$coefficients[c(4:6),2]),
                     b = c(summary(grg.mod)$coefficients[c(1:3),1]),
                     b.se = c(summary(grg.mod)$coefficients[c(1:3),2]))
    
plot(grgc.df$atr, grgc.df$e, pch = 16, xlab = 'atrazine (ppb)', ylab = 'LL.2 Parameters',
     ylim = c(0,17), xlim = c(0,300))
  points(grgc.df$atr+3, grgc.df$b, pch = 17, col=2)
    for(i in 1:length(grgc.df$atr)){
      segments(x0 = grgc.df$atr[i], y0 = grgc.df$e[i] + grgc.df$e.se[i],
               x1 = grgc.df$atr[i], y1 = grgc.df$e[i] - grgc.df$e.se[i])
      segments(x0 = grgc.df$atr[i]+3, y0 = grgc.df$b[i] + grgc.df$b.se[i],
               x1 = grgc.df$atr[i]+3, y1 = grgc.df$b[i] - grgc.df$b.se[i], col=2)
    }
#parameters as function of atrazine  
  eg.mod = lm(e ~ atr, weights = e.se^-1, data = grgc.df) 
  eg.mod2 = lm(e ~ logatr, weights = e.se^-1, data = grgc.df)  
    AIC(eg.mod, eg.mod2) #exponential fits way better
  bg.mod = lm(b ~ atr, weights = b.se^-1, data = grgc.df) 
  
  modgdf= data.frame(atr = c(0:300),
                     logatr = log(c(0:300)+1),
                     pred.e = 0,
                     pred.e.se = 0,
                     pred.e2 = 0,
                     pred.e2.se = 0,
                     pred.b = 0,
                     pred.b.se = 0)
  
  modgdf[,3:4] = predict(eg.mod, newdata = modgdf, se.fit = TRUE)[1:2]
  modgdf[,5:6] = predict(eg.mod2, newdata = modgdf, se.fit = TRUE)[1:2]
  modgdf[,7:8] = predict(bg.mod, newdata = modgdf, se.fit = TRUE)[1:2]
  
lines(modgdf$atr, modgdf$pred.e, lty = 2)
  lines(modgdf$atr, modgdf$pred.e + 1.96*modgdf$pred.e.se, lty = 3)
  lines(modgdf$atr, modgdf$pred.e - 1.96*modgdf$pred.e.se, lty = 3)

lines(modgdf$atr, modgdf$pred.e2, lty = 2, col=4)
  lines(modgdf$atr, modgdf$pred.e2 + 1.96*modgdf$pred.e2.se, lty = 3, col=4)
  lines(modgdf$atr, modgdf$pred.e2 - 1.96*modgdf$pred.e2.se, lty = 3, col=4)  
  
lines(modgdf$atr, modgdf$pred.b, lty = 2, col=2)
  lines(modgdf$atr, modgdf$pred.b + 1.96*modgdf$pred.b.se, lty = 3, col=2)
  lines(modgdf$atr, modgdf$pred.b - 1.96*modgdf$pred.b.se, lty = 3, col=2)
  
legend('topright', legend = c('slp', 'lc50'), pch = c(17, 16), col=c(2,1), cex = 0.7, bty='n')  
title('Griggs 08 cercarial survival parameters')  
#Create function to generate d-r function with linear fit to lc50 parameter#####################
piC.grg08_atr_unc = function(He){
  if(He == 0) piC = 1 else{
    e = as.numeric(predict(eg.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bg.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, e[1], e[2])
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, stop.on.error = FALSE)[1]$value
    piC = auc/grg08_atr_aucs[1]
    if(piC > 1) piC = 1
  }
  return(piC)
  }  
  
#Create function to generate d-r function with exponential fit to lc50 parameter#####################
piC.grg08_atr_unc2 = function(He){
    if(He == 0) piC = 1 else{
      e = as.numeric(predict(eg.mod2, newdata = data.frame(logatr = log(He+1)), se.fit = TRUE)[1:2])
      b = as.numeric(predict(bg.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
      
      e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
      b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, e[1], e[2])
      auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, stop.on.error = FALSE)[1]$value
      piC = auc/grg08_atr_aucs[1]
      if(piC > 1) piC = 1
    }
    return(piC)
    
  }  
  
#Compare the linear and exponential functions and store key items ##############
plot(c(0,15,100), grg08_atr_aucs/grg08_atr_aucs[1], ylim = c(0,1),pch = 16, cex = 1.2,
     xlab = 'atrazine (ppb)', ylab = expression(paste(pi[C], 'estimate')))
  
  points(c(0:500), sapply(c(0:500), piC.grg08_atr_unc), pch = 5, col = 2, cex = 0.5)
  points(c(0:500), sapply(c(0:500), piC.grg08_atr_unc2), pch = 5, cex = 0.5, col = 4)
  legend('bottomleft', legend = c('linear', 'exponential'), pch = 5, col = c(2,4),
         title = 'Function fit to lc50 parameter', cex = 0.7, bty='n')
  
keep.grg08 = c('L.3.fx', 'grg08_atr_aucs', 'piC.grg08_atr_unc', 'piC.grg08_atr_unc2', 'grgc.df',
               'eg.mod', 'eg.mod2', 'bg.mod')

#Plot sample function output across atrazine conc
plot(c(0:500), sapply(c(0:500), piC.grg08_atr_unc2), pch = 17, cex = 0.5,
     xlab = 'Atrazine (ppb)', ylab = expression(paste(pi[C])),
     main = 'Sample Output of cercarial mortality function')