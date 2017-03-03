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

cerc.g = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/Griggs2008.csv')
  time = c(0:25)

#Plot control survival curve, estimate time-dependent die-off function, and get auc from 0-24 hours ########
cerc.g0 = subset(cerc.g, chem == 'solvent')
  plot(x = cerc.g0$time_hrs,y = cerc.g0$surv/100,
       xlab = 'time (hrs)', ylab = 'prop alive', pch = 16, xlim = c(0,25), ylim = c(0,1))

  grg.ctrl = drm(alive/total ~ time_hrs, total, data = cerc.g0, 
                 type = 'binomial', fct = LL.2())
    
  lines(time, ll4(1,0,coef(grg.ctrl)[1], coef(grg.ctrl)[2], time), lty=2)
  
  grg.fx0 = function(t){
    1/(1+exp(summary(grg.ctrl)$coefficients[1]*(log(t/summary(grg.ctrl)$coefficients[2]))))
  }
  
  auc.grg0 = integrate(grg.fx0, lower = 0, upper = 24)$value
  
#Function for low dose treatment (atrazine 15 ppb & metolachlor 10 ppb combined)  ##########
cerc.g15 = subset(cerc.g, conc == 15)
  points(cerc.g15$time_hrs, cerc.g15$surv/100, pch = 17)
  
  grg.15 = drm(alive/total ~ time_hrs, total, data = cerc.g15, 
               type = 'binomial', fct = LL.2())
  
  lines(time, ll4(1,0,coef(grg.15)[1], coef(grg.15)[2], time), lty=3)
  
  grg.fx15 = function(t){
    1/(1+exp(summary(grg.15)$coefficients[1]*(log(t/summary(grg.15)$coefficients[2]))))
  }
  
  auc.grg15 = integrate(grg.fx15, lower = 0, upper = 24)$value
  
#Function for low dose treatment (atrazine 100 ppb & metolachlor 85 ppb combined)   ###############
cerc.g100 = subset(cerc.g, conc == 100)
  points(cerc.g100$time_hrs, cerc.g100$surv/100, pch = 15)
  
  grg.100 = drm(alive/total ~ time_hrs, total, data = cerc.g100, 
                type = 'binomial', fct = LL.2())
  
  lines(time, ll4(1,0,coef(grg.100)[1], coef(grg.100)[2], time), lty=4)
  
  grg.fx100 = function(t){
    1/(1+exp(summary(grg.100)$coefficients[1]*(log(t/summary(grg.100)$coefficients[2]))))
  }
  
  auc.grg100 = integrate(grg.fx100, lower = 0, upper = 24)$value
  
  title(main='Griggs2008 Atrazine-Cercarial mortality (E.trivolvis)')
    legend('bottomleft', legend = c('control', '15ppb', '100ppb'), pch = c(16,17,15), cex=0.8)
    
    grgatr.auc = data.frame(atr = unique(cerc.g$conc),
                            auc = c(auc.grg0, auc.grg15, auc.grg100),
                            piC = c(auc.grg0, auc.grg15, auc.grg100)/auc.grg0)

#Create data frame with parameter values and atrazine concentrations #######################
grgc.df = data.frame(atr = c(0,15,100),
                     logatr = log(c(0,15,100)+1),
                     e = c(coef(grg.ctrl)[2], coef(grg.15)[2], coef(grg.100)[2]),
                     e.se = c(summary(grg.ctrl)$coefficients[2,2], summary(grg.15)$coefficients[2,2],
                                  summary(grg.100)$coefficients[2,2]),
                     b = c(coef(grg.ctrl)[1], coef(grg.15)[1], coef(grg.100)[1]),
                     b.se = c(summary(grg.ctrl)$coefficients[1,2], summary(grg.15)$coefficients[1,2],
                              summary(grg.100)$coefficients[1,2]))
    
  plot(grgc.df$atr, grgc.df$e, pch = 16, xlab = 'atrazine (ppb)', ylab = 'LL.2 Parameters',
       ylim = c(2,18), xlim = c(0,300))
  points(grgc.df$atr, grgc.df$b, pch = 17, col=2)
    for(i in 1:length(grgc.df$atr)){
      segments(x0 = grgc.df$atr[i], y0 = grgc.df$e[i] + grgc.df$e.se[i],
               x1 = grgc.df$atr[i], y1 = grgc.df$e[i] - grgc.df$e.se[i])
      segments(x0 = grgc.df$atr[i], y0 = grgc.df$b[i] + grgc.df$b.se[i],
               x1 = grgc.df$atr[i], y1 = grgc.df$b[i] - grgc.df$b.se[i], col=2)
    }
#parameters as function of atrazine  
  eg.mod = lm(e ~ atr, weights = e.se^-1, data = grgc.df) 
  eg.mod2 = lm(e ~ logatr, weights = e.se^-1, data = grgc.df)  
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
  
legend('topright', legend = c('slp', 'lc50'), pch = c(17, 16), col=c(2,1), cex = 0.7)  
title('Griggs 08 cercarial survival parameters')  
#Create function to generate d-r function with linear fit to lc50 parameter#####################
  grg.fx = function(He){
    e = as.numeric(predict(eg.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bg.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    
    if(e.use <= 0){auc=0} else {
      auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t / e.use)))))}, 
                      lower=0, upper=24, stop.on.error = FALSE)[1]$value
    }
    
    auc
  }  
  
#Final:generate relative cercariae-hrs function  
  piC.grg08_atr_unc = function(He){
    piC = grg.fx(He) / grg.fx(0)
    if(piC > 1) {piC = 1} else {
      (return(piC))
    }
  }
  
#Create function to generate d-r function with linear fit to lc50 parameter#####################
  grg.fx2 = function(He){
    e = as.numeric(predict(eg.mod2, newdata = data.frame(logatr = log(He+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bg.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    
    if(e.use <= 0){auc=0} else {
      auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t / e.use)))))}, 
                      lower=0, upper=24, stop.on.error = FALSE)[1]$value
    }
    
    auc
  }  
  
  #Final:generate relative cercariae-hrs function  
  piC.grg08_atr_unc2 = function(He){
    piC = grg.fx(He) / grg.fx(0)
    if(piC > 1) {piC = 1} else {
      (return(piC))
    }
  }
  
#Compare the linear and exponential functions and store key items ##############
plot(c(0:500), sapply(c(0:500), piC.grg08_atr_unc), pch = 1, ylim = c(0,1),
     xlab = 'atrazine (ppb)', ylab = expression(paste(pi[C], 'estimate')))
  points(c(0:500), sapply(c(0:500), piC.grg08_atr_unc2), pch = 1, col = 4)
  legend('bottomleft', legend = c('linear', 'exponential'), pch = 1, col = c(1,4),
         title = 'Function fit to lc50 parameter', cex = 0.7)
  
keep.grg08 = c('grgatr.auc', 'piC.grg08_atr_unc', 'grg.fx', 'piC.grg08_atr_unc2', 'grg.fx2',
               'eg.mod', 'eg.mod2', 'bg.mod')

#Plot sample function output across atrazine conc
plot(c(0:500), sapply(c(0:500), piC.grg08_atr_unc2), pch = 17, cex = 0.5,
     xlab = 'Atrazine (ppb)', ylab = expression(paste(pi[C])),
     main = 'Sample Output of cercarial mortality function')