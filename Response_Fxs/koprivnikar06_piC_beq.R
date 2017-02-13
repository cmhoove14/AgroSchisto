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

#20 ppb atrazine ##################
kop.20 = subset(kop.c, conc == 20)
  
  points(x = kop.20$time_hrs, y = kop.20$surv, pch = 17)

  kop.20mod = nls(surv ~ 1/(1+exp(slp*(log(time_hrs / lc50)))), data = kop.20,
                  start = list(slp = 5, lc50 = 2.5), weights = surv.se^-1)
  
  lines(time, ll4(1,0,summary(kop.20mod)$coefficients[1],
                  summary(kop.20mod)$coefficients[2], time), lty=3)
  
#200 ppb atrazine ##################
kop.200 = subset(kop.c, conc == 200)
  
  points(x = kop.200$time_hrs, y = kop.200$surv, pch = 15)
  
  kop.200mod = nls(surv ~ 1/(1+exp(slp*(log(time_hrs / lc50)))), data = kop.200,
                   start = list(slp = 5, lc50 = 2.5), weights = surv.se^-1)
  
  lines(time, ll4(1,0,summary(kop.200mod)$coefficients[1],
                  summary(kop.200mod)$coefficients[2], time), lty=4)
  
  title(main='Koprivnikar2006 Cercarial mortality (E.trivolvis)')
  legend('topright', legend = c('control', '20ppb', '200ppb'), pch = c(16,17,15), cex=0.8)

#Create data frame with parameter values and atrazine concentrations #######################
kopatr.df = data.frame(atr = c(0,20,200),
                       e = c(coef(kop.ctrl)[2], coef(kop.20mod)[2], coef(kop.200mod)[2]),
                       e.se = c(summary(kop.ctrl)$coefficients[2,2], summary(kop.20mod)$coefficients[2,2],
                                summary(kop.200mod)$coefficients[2,2]),
                       b = c(coef(kop.ctrl)[1], coef(kop.20mod)[1], coef(kop.200mod)[1]),
                       b.se = c(summary(kop.ctrl)$coefficients[1,2], summary(kop.20mod)$coefficients[1,2],
                                summary(kop.200mod)$coefficients[1,2]))
  
plot(kopatr.df$atr, kopatr.df$e, pch = 16, xlab = 'atrazine (ppb)', ylab = 'LL.2 Parameters',
     ylim = c(0, 20))
  points(kopatr.df$atr, kopatr.df$b, pch = 17, col=2)
  for(i in 1:length(kopatr.df$atr)){
    segments(x0 = kopatr.df$atr[i], y0 = kopatr.df$e[i] + kopatr.df$e.se[i],
             x1 = kopatr.df$atr[i], y1 = kopatr.df$e[i] - kopatr.df$e.se[i])
    segments(x0 = kopatr.df$atr[i], y0 = kopatr.df$b[i] + kopatr.df$b.se[i],
             x1 = kopatr.df$atr[i], y1 = kopatr.df$b[i] - kopatr.df$b.se[i], col=2)
  }
  
  ek.mod = lm(e ~ atr, weights = e.se^-1, data = kopatr.df) 
  bk.mod = lm(b ~ atr, weights = b.se^-1, data = kopatr.df) 
  
modkopdf= data.frame(atr = c(0:250),
                     pred.e = 0,
                     pred.e.se = 0,
                     pred.b = 0,
                     pred.b.se = 0)
  
  modkopdf[,2:3] = predict(ek.mod, newdata = modkopdf, se.fit = TRUE)[1:2]
  modkopdf[,4:5] = predict(bk.mod, newdata = modkopdf, se.fit = TRUE)[1:2]
  
lines(modkopdf$atr, modkopdf$pred.e, lty = 2)
  lines(modkopdf$atr, modkopdf$pred.e + 1.96*modkopdf$pred.e.se, lty = 3)
  lines(modkopdf$atr, modkopdf$pred.e - 1.96*modkopdf$pred.e.se, lty = 3)
  
lines(modkopdf$atr, modkopdf$pred.b, lty = 2, col=2)
  lines(modkopdf$atr, modkopdf$pred.b + 1.96*modkopdf$pred.b.se, lty = 3, col=2)
  lines(modkopdf$atr, modkopdf$pred.b - 1.96*modkopdf$pred.b.se, lty = 3, col=2)
  
legend('topright', legend = c('slp', 'lc50'), pch = c(17, 16), col=c(2,1), cex = 0.7)  
title('Koprivnikar 06 cercarial survival parameters')  
#Create function to generate d-r function #####################
  predk.fx = function(He){
    e = as.numeric(predict(ek.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bk.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    
    lines(time, ll4(1,0,b.use, log(e.use), time), lty=2, col=3)
    
    auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t)-e.use))))}, 
                    lower=0, upper=24)[1]$value
    auc
  }  
  
#Final:generate relative cercariae-hrs function  
  piC.kop_atr_unc = function(He){
    piC = predk.fx(In) / predk.fx(0)
    if(piC > 1) piC = 1
    else(return(piC))
  }
  
  keep.kop06.beq = c('piC.kop_atr_unc', 'predk.fx', 'ek.mod', 'bk.mod')