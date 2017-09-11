#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Ghaffar 2016 LARVAL (S.mansoni) data
require(drc)

  but.vals = c(0, 556, 2417, 3906, 5560, 8703)
  gly.vals = c(0, 1506, 3875, 9174, 15062, 26249)
  pen.vals = c(0, 214.8, 535, 1299, 2148, 3762)
  
  ll3 = function(slp,lc,x){
    1/(1+exp(slp*(log(x/lc))))
  }
  
  time = seq(0,25,0.1)

#Cercarial (S. mansoni) toxicity ###############
  cer<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/ghaffar2016_cercariae.csv')
    cer$alive = 100*cer$surv
    cer$dead = 100-cer$alive
    cer$total = 100
    
  cer.ctrl = subset(cer, chem == 'control')  
    
    gaf.ctrl.drc.cer = drm(alive/total ~ time_hrs, weights = total, data = cer.ctrl, type = 'binomial', 
                      fct = LL.4(names = c('b', 'c', 'd', 'e'),
                                fixed = c(NA, 0, 1, NA)))
    summary(gaf.ctrl.drc.cer)

#butralin cercariae ###############
  cer.but = subset(cer, chem == 'butralin')
  
  gaf.but.drc.cer = drm(alive/total ~ time_hrs, conc, weights = total, data = cer.but, type = 'binomial', 
                    fct = LL.4(names = c('b', 'c', 'd', 'e'),
                               fixed = c(NA, 0, 1, NA)))
  summary(gaf.but.drc.cer)
  plot(gaf.but.drc.cer) 
  
  plot(cer.ctrl$time_hrs, cer.ctrl$surv, ylim = c(0,1), xlim = c(0,24), pch=16, 
       xlab = 'time(hrs)', ylab = 'prop surviving')
    lines(time, ll3(lc = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], x = time), lty = 2)
    
    for(i in 1:length(unique(cer.but$conc))){
      points(cer.but$time_hrs[cer.but$conc==unique(cer.but$conc)[i]], 
             cer.but$surv[cer.but$conc==unique(cer.but$conc)[i]], pch=17,
             col = i+1)
      lines(time, ll3(lc = gaf.but.drc.cer$coefficients[i+5], slp = gaf.but.drc.cer$coefficients[i], x = time), lty = 2, col = i+1)
    }
  legend('topright', title = 'Butralin (ppb)', legend = but.vals, pch = c(16,rep(17,5)),
         col = c(1:6), cex=0.7, bty = 'n')
  title('Ghaffar2016 butralin toxicity to cercariae')

#Get estimate of cercariae-hrs for each concentration    
gaf.but.aucs.cer = as.numeric()
  gaf.but.aucs.cer[1] = integrate(f = ll3, lc = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], 
                                  lower=0, upper=24)[1]$value  
  
  for(j in 1:length(unique(cer.but$conc))){
    gaf.but.aucs.cer[j+1] = integrate(f = ll3, lc = gaf.but.drc.cer$coefficients[j+5], slp = gaf.but.drc.cer$coefficients[j],
                                      lower=0, upper=24)[1]$value  
  }
  
#Compile L.2 data for functional responses #############
but.cer = data.frame(butralin = but.vals,
                     logbutralin = log(but.vals+1),
                     e = c(summary(gaf.ctrl.drc.cer)$coefficients[2,1], summary(gaf.but.drc.cer)$coefficients[c(6:10),1]),
                     e.se = c(summary(gaf.ctrl.drc.cer)$coefficients[2,2], summary(gaf.but.drc.cer)$coefficients[c(6:10),2]),
                     b = c(summary(gaf.ctrl.drc.cer)$coefficients[1,1], summary(gaf.but.drc.cer)$coefficients[c(1:5),1]),
                     b.se = c(summary(gaf.ctrl.drc.cer)$coefficients[1,2], summary(gaf.but.drc.cer)$coefficients[c(1:5),2]))

plot(but.cer$butralin, but.cer$e, pch = 16, xlab = 'butralin (ppb)', 
     ylab = 'L.2 Parameters (cercariae)', ylim = c(0,10))

  points(but.cer$butralin+200, but.cer$b, pch = 17, col=2)

for(i in 1:length(but.cer$butralin)){
  segments(x0 = but.cer$butralin[i], y0 = but.cer$e[i] + but.cer$e.se[i],
           x1 = but.cer$butralin[i], y1 = but.cer$e[i] - but.cer$e.se[i])
  segments(x0 = but.cer$butralin[i]+200, y0 = but.cer$b[i] + but.cer$b.se[i],
           x1 = but.cer$butralin[i]+200, y1 = but.cer$b[i] - but.cer$b.se[i], col=2)
} 

#functions fit to L.2 parameters ############
but.cer.lm.e = lm(e ~ butralin, weights = e.se^-1, data = but.cer) #linear response of LC50
  but.cer.pred = function(but){
    predict(but.cer.lm.e, newdata = data.frame(butralin = but), 
            interval = 'confidence', level = 0.95)
  }
    
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred, simplify = T)[3,], lty = 3)

but.cer.lm.e2 = lm(e ~ logbutralin, weights = e.se^-1, data = but.cer) #log-linear response of LC50
  but.cer.pred2 = function(but){
    predict(but.cer.lm.e2, newdata = data.frame(logbutralin = log(but+1)), 
            interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred2, simplify = T)[3,], lty = 3, col=3)

  AIC(but.cer.lm.e, but.cer.lm.e2)  #exponential is a better fit    

but.cer.lm.b = lm(b ~ butralin, weights = b.se^-1, data = but.cer)   
  but.cer.pred.b = function(but){
    predict(but.cer.lm.b, newdata = data.frame(butralin = but), interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,9000,100), sapply(seq(0,9000,100), but.cer.pred.b, simplify = T)[3,], lty = 3, col = 2)

legend('topleft', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
       legend = c('LC50 - Linear',
                  'LC50 - Exponential',
                  'Slp - linear',
                  '95% CI'), cex = 0.7)  
legend('top', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
       cex = 0.8, bty = 'n', title = 'Observed points')

piC.ghaf_butr.lin_unc = function(He){
  if(He == 0) piC = 1 else{
    e = as.numeric(predict(but.cer.lm.e, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(but.cer.lm.b, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, e[1], e[2])
    auc = integrate(ll3, lc = e.use, slp = b.use, lower=0, upper=24)[1]$value
    piC = auc/gaf.but.aucs.cer[1]
    if(piC > 1) piC = 1
  }
  return(piC)
}  #function to estimate AUC 

piC.ghaf_butr.exp_unc = function(He){
  if(He == 0) piC = 1 else{
    e = as.numeric(predict(but.cer.lm.e2, newdata = data.frame(logbutralin = log(He+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(but.cer.lm.b, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, e[1], e[2])
    auc = integrate(ll3, lc = e.use, slp = b.use, lower=0, upper=24)[1]$value
    piC = auc/gaf.but.aucs.cer[1]
    if(piC > 1) piC = 1
  }
  return(piC)
}  #function to estimate AUC 

keep.gaf.but.cer = c('piC.ghaf_butr.lin_unc', 'but.cer.lm.e', 'but.cer.lm.b', 'll3', 'gaf.but.aucs.cer',
                     'piC.ghaf_butr.exp_unc', 'but.cer.lm.e2')
#plot to test function ########
plot(but.cer$butralin, gaf.but.aucs.cer/gaf.but.aucs.cer[1],
     xlab = 'Butralin (ppb)', ylab = 'relative AUC (cercariae-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0, 18000,50)/2, sapply(seq(0, 18000,50)/2, piC.ghaf_butr.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0, 18000,50)/2, sapply(seq(0, 18000,50)/2, piC.ghaf_butr.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)

#glyphosate cercariae ###############
cer.gly = subset(cer, chem == 'glyphosate')
  
  gaf.gly.drc.cer = drm(alive/total ~ time_hrs, conc, weights = total, data = cer.gly, type = 'binomial', 
                    fct = LL.4(names = c('b', 'c', 'd', 'e'),
                              fixed = c(NA, 0, 1, NA)))
  summary(gaf.gly.drc.cer)
  plot(gaf.gly.drc.cer) 
  
  plot(cer$time_hrs[cer$conc==0], cer$surv[cer$conc==0], 
       pch=16, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    lines(time, ll3(lc = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], x = time), lty = 2)
  
  for(i in 1:length(unique(cer.gly$conc))){
    points(cer.gly$time_hrs[cer.gly$conc==unique(cer.gly$conc)[i]], 
           cer.gly$surv[cer.gly$conc==unique(cer.gly$conc)[i]], pch=17,
           col = i+1)
    lines(time, ll3(lc = gaf.gly.drc.cer$coefficients[i+5], slp = gaf.gly.drc.cer$coefficients[i], x = time), lty = 2, col = i+1)
  }
  legend('topright', title = 'Glyphosate (ppb)', legend = gly.vals, pch = c(16,rep(17,5)),
         col = c(1:6), cex=0.7, bty = 'n')
  title('Ghaffar2016 Glyphosate toxicity to cercariae')
  
#Get estimate of cercariae-hrs for each concentration    
gaf.gly.aucs.cer = as.numeric()
  gaf.gly.aucs.cer[1] = integrate(f = ll3, lc = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], 
                                  lower=0, upper=24)[1]$value  
  
  for(j in 1:length(unique(cer.gly$conc))){
    gaf.gly.aucs.cer[j+1] = integrate(f = ll3, lc = gaf.gly.drc.cer$coefficients[j+5], slp = gaf.gly.drc.cer$coefficients[j],
                                      lower=0, upper=24)[1]$value  
  }
  

#Compile L.2 parameters ######## 
gly.cer = data.frame(glyphosate = gly.vals,
                     logglyphosate = log(gly.vals+1),
                     e = c(summary(gaf.ctrl.drc.cer)$coefficients[2,1], summary(gaf.gly.drc.cer)$coefficients[c(6:10),1]),
                     e.se = c(summary(gaf.ctrl.drc.cer)$coefficients[2,2], summary(gaf.gly.drc.cer)$coefficients[c(6:10),2]),
                     b = c(summary(gaf.ctrl.drc.cer)$coefficients[1,1], summary(gaf.gly.drc.cer)$coefficients[c(1:5),1]),
                     b.se = c(summary(gaf.ctrl.drc.cer)$coefficients[1,2], summary(gaf.gly.drc.cer)$coefficients[c(1:5),2]))

plot(gly.cer$glyphosate, gly.cer$e, pch = 16, xlab = 'glyphosate (ppb)', 
     ylab = 'L.2 Parameters (cercariae)', ylim = c(0,6))

  points(gly.cer$glyphosate+200, gly.cer$b, pch = 17, col=2)

for(i in 1:length(gly.cer$glyphosate)){
  segments(x0 = gly.cer$glyphosate[i], y0 = gly.cer$e[i] + gly.cer$e.se[i],
           x1 = gly.cer$glyphosate[i], y1 = gly.cer$e[i] - gly.cer$e.se[i])
  segments(x0 = gly.cer$glyphosate[i]+200, y0 = gly.cer$b[i] + gly.cer$b.se[i],
           x1 = gly.cer$glyphosate[i]+200, y1 = gly.cer$b[i] - gly.cer$b.se[i], col=2)
} 

#fit functions to LL.2 parameters ##########  
gly.cer.lm.e = lm(e ~ glyphosate, weights = e.se^-1, data = gly.cer) #linear response of LC50
  gly.cer.pred = function(gly){
    predict(gly.cer.lm.e, newdata = data.frame(glyphosate = gly), 
            interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred, simplify = T)[3,], lty = 3)

gly.cer.lm.e2 = lm(e ~ logglyphosate, weights = e.se^-1, data = gly.cer) #log-linear response of LC50
  gly.cer.pred2 = function(gly){
    predict(gly.cer.lm.e2, newdata = data.frame(logglyphosate = log(gly+1)), 
            interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred2, simplify = T)[3,], lty = 3, col=3)

  AIC(gly.cer.lm.e, gly.cer.lm.e2)  #Exponential is a better fit    

gly.cer.lm.b = lm(b ~ glyphosate, weights = b.se^-1, data = gly.cer)   
  gly.cer.pred.b = function(gly){
    predict(gly.cer.lm.b, newdata = data.frame(glyphosate = gly), interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,9000,100)*3, sapply(seq(0,9000,100)*3, gly.cer.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

piC.ghaf_gly.lin_unc = function(He){
  if(He == 0) piC = 1 else{
    e = as.numeric(predict(gly.cer.lm.e, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(gly.cer.lm.b, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, e[1], e[2])
    auc = integrate(ll3, lc = e.use, slp = b.use, lower=0, upper=24)[1]$value
    piC = auc/gaf.gly.aucs.cer[1]
    if(piC > 1) piC = 1
  }
  return(piC)
}  #function to estimate AUC 

piC.ghaf_gly.exp_unc = function(He){
  if(He == 0) piC = 1 else{
    e = as.numeric(predict(gly.cer.lm.e2, newdata = data.frame(logglyphosate = log(He+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(gly.cer.lm.b, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, e[1], e[2])
    auc = integrate(ll3, lc = e.use, slp = b.use, lower=0, upper=24)[1]$value
    piC = auc/gaf.gly.aucs.cer[1]
    if(piC > 1) piC = 1
  }
  return(piC)
}  #function to estimate AUC 

keep.gaf.gly.cer = c('piC.ghaf_gly.lin_unc', 'gly.cer.lm.e', 'gly.cer.lm.b', 'll3', 'gaf.gly.aucs.cer',
                     'piC.ghaf_gly.exp_unc', 'gly.cer.lm.e2')

#plot to test function ########
plot(gly.cer$glyphosate, gaf.gly.aucs.cer/gaf.gly.aucs.cer[1],
     xlab = 'Glyphosate (ppb)', ylab = 'relative AUC (cercariae-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0,9000,50)*3, sapply(seq(0,9000,50)*3, piC.ghaf_gly.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,9000,50)*3, sapply(seq(0,9000,50)*3, piC.ghaf_gly.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)



#pendimethalin cercariae ###############
cer.pen = subset(cer, chem == 'pendimethalin')

  gaf.pen.drc.cer = drm(alive/total ~ time_hrs, conc, weights = total, data = cer.pen, type = 'binomial', 
                    fct = LL.4(names = c('b', 'c', 'd', 'e'),
                               fixed = c(NA, 0, 1, NA)))
  summary(gaf.pen.drc.cer)
  plot(gaf.pen.drc.cer) 

plot(cer$time_hrs[cer$conc==0], cer$surv[cer$conc==0], 
     pch=16, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, ll3(lc = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], x = time), lty = 2)
  
  for(i in 1:length(unique(cer.pen$conc))){
    points(cer.pen$time_hrs[cer.pen$conc==unique(cer.pen$conc)[i]], 
           cer.pen$surv[cer.pen$conc==unique(cer.pen$conc)[i]], pch=17,
           col = i+1)
    lines(time, ll3(lc = gaf.pen.drc.cer$coefficients[i+5], slp = gaf.pen.drc.cer$coefficients[i], x = time), lty = 2, col = i+1)
  }
    legend('topright', title = 'Pendimethalin (ppb)', legend = pen.vals, pch = c(16,rep(17,5)),
           col = c(1:6), cex=0.7, bty = 'n')
    title('Ghaffar2016 Pendimethalin toxicity to cercariae')

#Get estimate of cercariae-hrs for each concentration    
gaf.pen.aucs.cer = as.numeric()
gaf.pen.aucs.cer[1] = integrate(f = ll3, lc = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], 
                                lower=0, upper=24)[1]$value  

  for(j in 1:length(unique(cer.pen$conc))){
    gaf.pen.aucs.cer[j+1] = integrate(f = ll3, lc = gaf.pen.drc.cer$coefficients[j+5], slp = gaf.pen.drc.cer$coefficients[j],
                                      lower=0, upper=24)[1]$value  
  }

#Compile L.2 parameters ######## 
pen.cer = data.frame(pendimethalin = pen.vals,
                     logpendimethalin = log(pen.vals+1),
                     e = c(summary(gaf.ctrl.drc.cer)$coefficients[2,1], summary(gaf.pen.drc.cer)$coefficients[c(6:10),1]),
                     e.se = c(summary(gaf.ctrl.drc.cer)$coefficients[2,2], summary(gaf.pen.drc.cer)$coefficients[c(6:10),2]),
                     b = c(summary(gaf.ctrl.drc.cer)$coefficients[1,1], summary(gaf.pen.drc.cer)$coefficients[c(1:5),1]),
                     b.se = c(summary(gaf.ctrl.drc.cer)$coefficients[1,2], summary(gaf.pen.drc.cer)$coefficients[c(1:5),2]))

plot(pen.cer$pendimethalin, pen.cer$e, pch = 16, xlab = 'Pendimethalin (ppb)', 
     ylab = 'L.2 Parameters (cercariae)', ylim = c(0,6))

  points(pen.cer$pendimethalin+50, pen.cer$b, pch = 17, col=2)

  for(i in 1:length(pen.cer$pendimethalin)){
    segments(x0 = pen.cer$pendimethalin[i], y0 = pen.cer$e[i] + pen.cer$e.se[i],
             x1 = pen.cer$pendimethalin[i], y1 = pen.cer$e[i] - pen.cer$e.se[i])
    segments(x0 = pen.cer$pendimethalin[i]+50, y0 = pen.cer$b[i] + pen.cer$b.se[i],
             x1 = pen.cer$pendimethalin[i]+50, y1 = pen.cer$b[i] - pen.cer$b.se[i], col=2)
  } 

#fit functions to LL.2 parameters ##########  
pen.cer.lm.e = lm(e ~ pendimethalin, weights = e.se^-1, data = pen.cer) #linear response of LC50
  pen.cer.pred = function(pen){
    predict(pen.cer.lm.e, newdata = data.frame(pendimethalin = pen), 
            interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred, simplify = T)[3,], lty = 3)

pen.cer.lm.e2 = lm(e ~ logpendimethalin, weights = e.se^-1, data = pen.cer) #log-linear response of LC50
  pen.cer.pred2 = function(pen){
    predict(pen.cer.lm.e2, newdata = data.frame(logpendimethalin = log(pen+1)), 
            interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred2, simplify = T)[3,], lty = 3, col=3)
  
  AIC(pen.cer.lm.e, pen.cer.lm.e2)  #Exponential is a better fit    

pen.cer.lm.b = lm(b ~ pendimethalin, weights = b.se^-1, data = pen.cer)   
  pen.cer.pred.b = function(pen){
    predict(pen.cer.lm.b, newdata = data.frame(pendimethalin = pen), interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,4000,100), sapply(seq(0,4000,100), pen.cer.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

piC.ghaf_pen.lin_unc = function(He){
  if(He == 0) piC = 1 else{
    e = as.numeric(predict(pen.cer.lm.e, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(pen.cer.lm.b, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, e[1], e[2])
    auc = integrate(ll3, lc = e.use, slp = b.use, lower=0, upper=24)[1]$value
    piC = auc/gaf.pen.aucs.cer[1]
    if(piC > 1) piC = 1
  }
  return(piC)
}  #function to estimate AUC 
  
piC.ghaf_pen.exp_unc = function(He){
  if(He == 0) piC = 1 else{
    e = as.numeric(predict(pen.cer.lm.e2, newdata = data.frame(logpendimethalin = log(He+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(pen.cer.lm.b, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, e[1], e[2])
    auc = integrate(ll3, lc = e.use, slp = b.use, lower=0, upper=24)[1]$value
    piC = auc/gaf.pen.aucs.cer[1]
    if(piC > 1) piC = 1
  }
  return(piC)
}  #function to estimate AUC 

keep.gaf.pen.cer = c('piC.ghaf_pen.lin_unc', 'pen.cer.lm.e', 'pen.cer.lm.b', 'll3', 'gaf.pen.aucs.cer',
                     'piC.ghaf_pen.exp_unc', 'pen.cer.lm.e2')

#plot to test function ########
plot(pen.cer$pendimethalin, gaf.pen.aucs.cer/gaf.pen.aucs.cer[1],
     xlab = 'pendimethalin (ppb)', ylab = 'relative AUC (cercariae-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0,4000,25), sapply(seq(0,4000,25), piC.ghaf_pen.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,4000,25), sapply(seq(0,4000,25), piC.ghaf_pen.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)



#Final keep vector
keep.gaf.piC = c(keep.gaf.pen.cer, keep.gaf.gly.cer, keep.gaf.but.cer)  