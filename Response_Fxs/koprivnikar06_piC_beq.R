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

#Cercarial mortality (E. trivolvis) from Koprivnikar 2006 ##################
kop.c = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/Koprivnikar2006.csv')
  kop.c[,5:8] = kop.c[,5:8]/100 #Convert survival measures to proportions
  time = c(0:25)
  
kop.cc = subset(kop.c, chem == 'control')  

plot(x = kop.cc$time_hrs, y = kop.cc$surv,
     xlab = 'time (hrs)', ylab = 'prop dead', pch = 16)

#We'll use nls for this study becuse of its implementation of weights
#Since the authors only report standard error of the cercarial survival experiments,
#using the inverse of the standard error in nls appropriately incorporates uncertainty into
#model fitting. Using drm function requires knowledge of the sample size to properly
#weight observations
#Control experiment ##########
  kop.ctrl = nls(surv ~ 1/(1+exp(slp*(time_hrs - lc50))), data = kop.cc,
                 start = list(slp = 5, lc50 = 2.5), weights = surv.se^-1)
  
  lines(time, L.3.fx(time, lc50 = summary(kop.ctrl)$coefficients[2],
                     slp = summary(kop.ctrl)$coefficients[1]), lty=2)
  
  auc.kop0 = integrate(L.3.fx, lc50 = summary(kop.ctrl)$coefficients[2], slp = summary(kop.ctrl)$coefficients[1], 
                       lower = 0, upper = 24)$value

#20 ppb atrazine ##################
kop.20 = subset(kop.c, conc == 20)
  
  points(x = kop.20$time_hrs, y = kop.20$surv, pch = 16, col = 2)

  kop.20mod = nls(surv ~ 1/(1+exp(slp*(time_hrs - lc50))), data = kop.20,
                  start = list(slp = 5, lc50 = 2.5), weights = surv.se^-1)
  
  lines(time, L.3.fx(time, lc50 = summary(kop.20mod)$coefficients[2],
                     slp = summary(kop.20mod)$coefficients[1]), lty=2, col = 2)
  
  auc.kop20 = integrate(L.3.fx, lc50 = summary(kop.20mod)$coefficients[2], slp = summary(kop.20mod)$coefficients[1], 
                        lower = 0, upper = 24)$value
  
#200 ppb atrazine ##################
kop.200 = subset(kop.c, conc == 200)
  
  points(x = kop.200$time_hrs, y = kop.200$surv, pch = 16, col = 3)
  
  kop.200mod = nls(surv ~ 1/(1+exp(slp*(time_hrs - lc50))), data = kop.200,
                   start = list(slp = 5, lc50 = 2.5), weights = surv.se^-1)
  
  lines(time, L.3.fx(time, lc50 = summary(kop.200mod)$coefficients[2],
                     slp = summary(kop.200mod)$coefficients[1]), lty=2, col = 3)
  
  auc.kop200 = integrate(L.3.fx, lc50 = summary(kop.200mod)$coefficients[2], slp = summary(kop.200mod)$coefficients[1], 
                         lower = 0, upper = 24)$value
  
  title(main= expression(paste('Koprivnikar 2006 Cercarial mortality ', italic('(E.trivolvis)'))))
  legend('topright', legend = c('control', '20ppb', '200ppb'), pch = c(16), col = c(1,2,3), cex=0.8, bty = 'n')
  
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
     ylim = c(0, 15), xlim = c(0,400))
  points(kopatr.df$atr, kopatr.df$b, pch = 17, col=2)
  for(i in 1:length(kopatr.df$atr)){
    segments(x0 = kopatr.df$atr[i], y0 = kopatr.df$e[i] + kopatr.df$e.se[i],
             x1 = kopatr.df$atr[i], y1 = kopatr.df$e[i] - kopatr.df$e.se[i])
    segments(x0 = kopatr.df$atr[i], y0 = kopatr.df$b[i] + kopatr.df$b.se[i],
             x1 = kopatr.df$atr[i], y1 = kopatr.df$b[i] - kopatr.df$b.se[i], col=2)
  }
  
  bk.mod = lm(b ~ atr, weights = b.se^-1, data = kopatr.df) 
    bk.mod.pred = function(He){
      predict(bk.mod, newdata = data.frame(atr = He), interval = 'confidence', level = 0.95)
    }
    
  ek.mod = lm(e ~ atr, weights = e.se^-1, data = kopatr.df) 
    ek.mod.pred = function(He){
      predict(ek.mod, newdata = data.frame(atr = He), interval = 'confidence', level = 0.95)
    }
    
  ek.mod2 = lm(e ~ logatr, weights = e.se^-1, data = kopatr.df)
    ek.mod.pred2 = function(He){
      predict(ek.mod2, newdata = data.frame(logatr = log(He+1)), interval = 'confidence', level = 0.95)
    }
    
  AIC(ek.mod, ek.mod2)  #linear is a better fit

lines(c(0:400), sapply(c(0:400), bk.mod.pred)[1,], lty = 2, col = 2)
  lines(c(0:400), sapply(c(0:400), bk.mod.pred)[2,], lty = 3, col = 2)
  lines(c(0:400), sapply(c(0:400), bk.mod.pred)[3,], lty = 3, col = 2)
  
lines(c(0:400), sapply(c(0:400), ek.mod.pred)[1,], lty = 2)
  lines(c(0:400), sapply(c(0:400), ek.mod.pred)[2,], lty = 3)
  lines(c(0:400), sapply(c(0:400), ek.mod.pred)[3,], lty = 3)
  
lines(c(0:400), sapply(c(0:400), ek.mod.pred2)[1,], lty = 2, col = 4)
  lines(c(0:400), sapply(c(0:400), ek.mod.pred2)[2,], lty = 3, col = 4)
  lines(c(0:400), sapply(c(0:400), ek.mod.pred2)[3,], lty = 3, col = 4)
  

  
legend('topright', legend = c('slp', 'lc50'), pch = c(17, 16), col=c(2,1), cex = 0.7, bty = 'n')  
title('Koprivnikar 06 cercarial survival parameters') 

#Create function to generate d-r function with linear b function#####################
  predk.fx = function(He){
    e = as.numeric(predict(ek.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bk.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    
    while(e.use < 0){
      e.use = rnorm(1, e[1], e[2])
    }      #resample if negative
    
    auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use, 
                    lower=0, upper=24)[1]$value
    
    auc
  }  
  
#Final: generate relative cercariae-hrs function  
  piC.kop_atr_unc = function(He){
    if(He == 0) piC = 1 
    else (piC = predk.fx(He) / auc.kop0)
    if(piC > 1) piC = 1
    return(piC)
  }
  
#Create function to generate d-r function with exponential b function#####################
  predk.fx2 = function(He){
    e = as.numeric(predict(ek.mod2, newdata = data.frame(logatr = log(He+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bk.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    
    while(e.use < 0){
      e.use = rnorm(1, e[1], e[2])
    }      #resample if negative
    
    auc = integrate(f = L.3.fx, lc50 = e.use, slp = b.use, 
                    lower=0, upper=24)[1]$value
    
    auc
  }  
  
  #Final:generate relative cercariae-hrs function  
  piC.kop_atr_unc2 = function(He){
    if(He == 0) piC = 1 
    else (piC = predk.fx2(He) / auc.kop0)
    if(piC > 1) piC = 1
    return(piC)
  }

#Compare the linear and exponential functions and store key items ##############
  plot(c(0,20,200), c(auc.kop0, auc.kop20, auc.kop200)/auc.kop0, pch = 16, xlim = c(0,500), ylim = c(0,1),
       xlab = 'Atrazine (ppb)', ylab = expression(paste(pi[C], 'estimate')))
    points(seq(0,500,5), sapply(seq(0,500,5), piC.kop_atr_unc), pch = 5, col = 4, cex = 0.5)
    points(seq(0,500,5), sapply(seq(0,500,5), piC.kop_atr_unc2), pch = 5, col = 2, cex = 0.5)
    legend('bottomleft', legend = c('linear', 'exponential'), pch = 5, col = c(4,2), bty = 'n',
           title = 'Function fit to lc50 parameter', cex = 0.7)
  
  
  keep.kop06.beq = c('kopatr.auc', 'piC.kop_atr_unc', 'predk.fx', 'piC.kop_atr_unc2', 'predk.fx2',
                     'ek.mod', 'bk.mod', 'bk.mod2')