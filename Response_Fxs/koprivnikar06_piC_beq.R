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

#Cercarial mortality (E. trivolvis) from Koprivnikar 2006 ##################
kop.c = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/Koprivnikar2006.csv')
  kop.c[,5:8] = kop.c[,5:8]/100 #Convert survival measures to proportions
  time = c(0:25)
  
kop.cc = subset(kop.c, chem == 'control')  

plot(x = kop.cc$time_hrs, y = kop.cc$surv,
     xlab = 'time (hrs)', ylab = 'prop dead', pch = 16)

#We'll use nls for this study becuse of its implementation of weights
#Since the authors only report standard error of the cercarial survival experiments,
#using the inverse of the standard error in nls appropriately incoporates uncertainty into
#model fitting. Using drm function requires knowledge of the sample size to properly
#weight observations

  kop.ctrl = nls(surv ~ 1/(1+exp(slp*(log(time_hrs / lc50)))), data = kop.cc,
                 start = list(slp = 5, lc50 = 2.5), weights = surv.se^-1)
  
  lines(time, ll4(1,0,summary(kop.ctrl)$coefficients[1],
                  summary(kop.ctrl)$coefficients[2], time), lty=2)
  
  kop.fx0 = function(t){
    1/(1+exp(summary(kop.ctrl)$coefficients[1]*(log(t/summary(kop.ctrl)$coefficients[2]))))
  }
  
  auc.kop0 = integrate(kop.fx0, lower = 0, upper = 24)$value

#20 ppb atrazine ##################
kop.20 = subset(kop.c, conc == 20)
  
  points(x = kop.20$time_hrs, y = kop.20$surv, pch = 17)

  kop.20mod = nls(surv ~ 1/(1+exp(slp*(log(time_hrs / lc50)))), data = kop.20,
                  start = list(slp = 5, lc50 = 2.5), weights = surv.se^-1)
  
  lines(time, ll4(1,0,summary(kop.20mod)$coefficients[1],
                  summary(kop.20mod)$coefficients[2], time), lty=3)
  
  kop.fx20 = function(t){
    1/(1+exp(summary(kop.20mod)$coefficients[1]*(log(t/summary(kop.20mod)$coefficients[2]))))
  }
  
  auc.kop20 = integrate(kop.fx20, lower = 0, upper = 24)$value
  
#200 ppb atrazine ##################
kop.200 = subset(kop.c, conc == 200)
  
  points(x = kop.200$time_hrs, y = kop.200$surv, pch = 15)
  
  kop.200mod = nls(surv ~ 1/(1+exp(slp*(log(time_hrs / lc50)))), data = kop.200,
                   start = list(slp = 5, lc50 = 2.5), weights = surv.se^-1)
  
  lines(time, ll4(1,0,summary(kop.200mod)$coefficients[1],
                  summary(kop.200mod)$coefficients[2], time), lty=4)
  
  kop.fx200 = function(t){
    1/(1+exp(summary(kop.200mod)$coefficients[1]*(log(t/summary(kop.200mod)$coefficients[2]))))
  }
  
  auc.kop200 = integrate(kop.fx200, lower = 0, upper = 24)$value
  
  title(main='Koprivnikar2006 Cercarial mortality (E.trivolvis)')
  legend('topright', legend = c('control', '20ppb', '200ppb'), pch = c(16,17,15), cex=0.8)
  
  kopatr.auc = data.frame(atr = unique(kop.c$conc),
                          auc = c(auc.kop0, auc.kop20, auc.kop200),
                          piC = c(auc.kop0, auc.kop20, auc.kop200)/auc.kop0)

#Create data frame with parameter values and atrazine concentrations #######################
kopatr.df = data.frame(atr = c(0,20,200),
                       logatr = log(c(0,20,200)+1),
                       e = c(coef(kop.ctrl)[2], coef(kop.20mod)[2], coef(kop.200mod)[2]),
                       e.se = c(summary(kop.ctrl)$coefficients[2,2], summary(kop.20mod)$coefficients[2,2],
                                summary(kop.200mod)$coefficients[2,2]),
                       b = c(coef(kop.ctrl)[1], coef(kop.20mod)[1], coef(kop.200mod)[1]),
                       b.se = c(summary(kop.ctrl)$coefficients[1,2], summary(kop.20mod)$coefficients[1,2],
                                summary(kop.200mod)$coefficients[1,2]))
  
plot(kopatr.df$atr, kopatr.df$e, pch = 16, xlab = 'atrazine (ppb)', ylab = 'LL.2 Parameters',
     ylim = c(0, 20), xlim = c(0,400))
  points(kopatr.df$atr, kopatr.df$b, pch = 17, col=2)
  for(i in 1:length(kopatr.df$atr)){
    segments(x0 = kopatr.df$atr[i], y0 = kopatr.df$e[i] + kopatr.df$e.se[i],
             x1 = kopatr.df$atr[i], y1 = kopatr.df$e[i] - kopatr.df$e.se[i])
    segments(x0 = kopatr.df$atr[i], y0 = kopatr.df$b[i] + kopatr.df$b.se[i],
             x1 = kopatr.df$atr[i], y1 = kopatr.df$b[i] - kopatr.df$b.se[i], col=2)
  }
  
  ek.mod = lm(e ~ atr, weights = e.se^-1, data = kopatr.df) 
  bk.mod = lm(b ~ atr, weights = b.se^-1, data = kopatr.df) 
  bk.mod2 = lm(b~logatr, weights = b.se^-1, data = kopatr.df)
  
modkopdf= data.frame(atr = c(0:400),
                     logatr = log(c(0:400)+1),
                     pred.e = 0,
                     pred.e.se = 0,
                     pred.b = 0,
                     pred.b.se = 0,
                     pred.b2 = 0,
                     pred.b2.se = 0)
  
  modkopdf[,3:4] = predict(ek.mod, newdata = modkopdf, se.fit = TRUE)[1:2]
  modkopdf[,5:6] = predict(bk.mod, newdata = modkopdf, se.fit = TRUE)[1:2]
  modkopdf[,7:8] = predict(bk.mod2, newdata = modkopdf, se.fit = TRUE)[1:2]
  
  
lines(modkopdf$atr, modkopdf$pred.e, lty = 2)
  lines(modkopdf$atr, modkopdf$pred.e + 1.96*modkopdf$pred.e.se, lty = 3)
  lines(modkopdf$atr, modkopdf$pred.e - 1.96*modkopdf$pred.e.se, lty = 3)
  
lines(modkopdf$atr, modkopdf$pred.b, lty = 2, col=2)
  lines(modkopdf$atr, modkopdf$pred.b + 1.96*modkopdf$pred.b.se, lty = 3, col=2)
  lines(modkopdf$atr, modkopdf$pred.b - 1.96*modkopdf$pred.b.se, lty = 3, col=2)
  
lines(modkopdf$atr, modkopdf$pred.b2, lty = 2, col=3)
  lines(modkopdf$atr, modkopdf$pred.b2 + 1.96*modkopdf$pred.b2.se, lty = 3, col=3)
  lines(modkopdf$atr, modkopdf$pred.b2 - 1.96*modkopdf$pred.b2.se, lty = 3, col=3)
  
legend('topright', legend = c('slp', 'lc50'), pch = c(17, 16), col=c(2,1), cex = 0.7)  
title('Koprivnikar 06 cercarial survival parameters')  
#Create function to generate d-r function with linear b function#####################
  predk.fx = function(He){
    e = as.numeric(predict(ek.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bk.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    
    #lines(time, ll4(1,0,b.use, log(e.use), time), lty=2, col=3)
    
    if(e.use <= 0){auc=0} else{
      auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t / e.use)))))}, 
                      lower=0, upper=24, stop.on.error = FALSE)[1]$value
    }
    
    auc
  }  
  
#Final:generate relative cercariae-hrs function  
  piC.kop_atr_unc = function(He){
    piC = predk.fx(He) / predk.fx(0)
    if(piC > 1) {piC = 1} else {
      return(piC)
    }
  }
  
#Create function to generate d-r function with exponential b function#####################
  predk.fx2 = function(He){
    e = as.numeric(predict(ek.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bk.mod2, newdata = data.frame(logatr = log(He+1)), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    
    #lines(time, ll4(1,0,b.use, log(e.use), time), lty=2, col=3)
    
    if(e.use <= 0){auc=0} else{
      auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t / e.use)))))}, 
                      lower=0, upper=24, stop.on.error = FALSE)[1]$value
    }
    
    auc
  }  
  
  #Final:generate relative cercariae-hrs function  
  piC.kop_atr_unc2 = function(He){
    piC = predk.fx2(He) / predk.fx2(0)
    if(piC > 1) {piC = 1} else {
      return(piC)
    }
  }

#Compare the linear and exponential functions and store key items ##############
  plot(c(0:500), sapply(c(0:500), piC.kop_atr_unc), pch = 1, col = 2, ylim = c(0,1),
       xlab = 'atrazine (ppb)', ylab = expression(paste(pi[C], 'estimate')))
    points(c(0:500), sapply(c(0:500), piC.kop_atr_unc2), pch = 1, col = 3)
    legend('bottomleft', legend = c('linear', 'exponential'), pch = 1, col = c(2,3),
           title = 'Function fit to slp parameter', cex = 0.7)
  
  
  keep.kop06.beq = c('kopatr.auc', 'piC.kop_atr_unc', 'predk.fx', 'piC.kop_atr_unc2', 'predk.fx2',
                     'ek.mod', 'bk.mod', 'bk.mod2')