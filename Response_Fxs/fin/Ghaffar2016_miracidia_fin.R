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

  but.vals = c(0, 556, 2417, 3906, 5560, 8703)/1000
  gly.vals = c(0, 1506, 3875, 9174, 15062, 26249)/1000
  pen.vals = c(0, 214.8, 535, 1299, 2148, 3762)/1000
  
L.3.fx = function(t, lc50 = lc50, slp = slp){
    1 / (1+exp(slp*log(t / lc50)))
  }  
  time = seq(0,25,0.1)

#Miracidial (S. mansoni) toxicity ###############
mir<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Miracidial Mortality/ghaffar2016_miracidia.csv')
  mir$conc = mir$conc/1000
  mir$alive = 100*mir$surv
  mir$dead = 100-mir$alive
  mir$total = 100
  
  mir.ctrl = subset(mir, chem == 'control')  
  
  gaf.ctrl.drc.mir = drm(alive/total ~ time_hrs, weights = total, data = mir.ctrl, type = 'binomial', 
                     fct = LL.3(names = c('b', 'd', 'e'),
                                fixed = c(NA, 1, NA)))
  summary(gaf.ctrl.drc.mir)
#butralin ###############
mir.but = subset(mir, chem == 'butralin')

gaf.but.drc.mir = drm(alive/total ~ time_hrs, conc, weights = total, data = mir.but, type = 'binomial', 
                  fct = LL.3(names = c('b', 'd', 'e'),
                             fixed = c(NA, 1, NA)))
  summary(gaf.but.drc.mir)
  plot(gaf.but.drc.mir) 

plot(mir.ctrl$time_hrs, mir.ctrl$surv, ylim = c(0,1), xlim = c(0,24), pch=16, 
     xlab = 'time(hrs)', ylab = 'prop surviving',
     main = 'Abdel-Ghaffar 2016: Butralin toxicity to miracidia')
  lines(time, L.3.fx(lc50 = gaf.ctrl.drc.mir$coefficients[2], slp = gaf.ctrl.drc.mir$coefficients[1], t =time), lty = 2)
  
  for(i in 1:length(unique(mir.but$conc))){
    points(mir.but$time_hrs[mir.but$conc==unique(mir.but$conc)[i]], 
           mir.but$surv[mir.but$conc==unique(mir.but$conc)[i]], pch=17,
           col = i+1)
    lines(time, L.3.fx(lc50 = gaf.but.drc.mir$coefficients[i+5], slp = gaf.but.drc.mir$coefficients[i], t =time), 
          lty = 2, col = i+1)
  }
  
  legend('topright', title = 'Butralin (ppb)', legend = but.vals, pch = c(16,rep(17,5)),
         col = c(1:6), cex=0.7, bty = 'n')

#Get estimate of miracidia-hrs for each concentration    
  gaf.but.aucs.mir = as.numeric()
  gaf.but.aucs.mir[1] = integrate(f = L.3.fx, lc50 = gaf.ctrl.drc.mir$coefficients[2], slp = gaf.ctrl.drc.mir$coefficients[1], 
                                  lower=0, upper=24)[1]$value  
  
  for(j in 1:length(unique(mir.but$conc))){
    gaf.but.aucs.mir[j+1] = integrate(f = L.3.fx, lc50 = gaf.but.drc.mir$coefficients[j+5], slp = gaf.but.drc.mir$coefficients[j],
                                      lower=0, upper=24)[1]$value  
  }
  
#Compile L.2 data for functional responses #############
but.mir = data.frame(butralin = but.vals,
                     logbutralin = log(but.vals+1),
                     e = c(summary(gaf.ctrl.drc.mir)$coefficients[2,1], summary(gaf.but.drc.mir)$coefficients[c(6:10),1]),
                     e.se = c(summary(gaf.ctrl.drc.mir)$coefficients[2,2], summary(gaf.but.drc.mir)$coefficients[c(6:10),2]),
                     b = c(summary(gaf.ctrl.drc.mir)$coefficients[1,1], summary(gaf.but.drc.mir)$coefficients[c(1:5),1]),
                     b.se = c(summary(gaf.ctrl.drc.mir)$coefficients[1,2], summary(gaf.but.drc.mir)$coefficients[c(1:5),2]))
  
plot(but.mir$butralin, but.mir$e, pch = 16, xlab = 'butralin (ppb)', 
     ylab = 'L.2 Parameters (miracidia)', ylim = c(0,6))

  points(but.mir$butralin, but.mir$b, pch = 17, col=2)

  for(i in 1:length(but.mir$butralin)){
    segments(x0 = but.mir$butralin[i], y0 = but.mir$e[i] + but.mir$e.se[i],
             x1 = but.mir$butralin[i], y1 = but.mir$e[i] - but.mir$e.se[i])
    segments(x0 = but.mir$butralin[i], y0 = but.mir$b[i] + but.mir$b.se[i],
             x1 = but.mir$butralin[i], y1 = but.mir$b[i] - but.mir$b.se[i], col=2)
  } 
  
#functions fit to L.2 parameters ############
but.mir.lm.e = lm(e ~ butralin, weights = e.se^-1, data = but.mir) #linear response of LC50
  but.mir.pred = function(but){
    predict(but.mir.lm.e, newdata = data.frame(butralin = but), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred, simplify = T)[3,], lty = 3)

but.mir.lm.e2 = lm(e ~ logbutralin, weights = e.se^-1, data = but.mir) #log-linear response of LC50
  but.mir.pred2 = function(but){
    predict(but.mir.lm.e2, newdata = data.frame(logbutralin = log(but+1)), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred2, simplify = T)[3,], lty = 3, col=3)

AIC(but.mir.lm.e, but.mir.lm.e2)  #Exponential is a better fit    

but.mir.lm.b = lm(b ~ butralin, weights = b.se^-1, data = but.mir)   
  but.mir.pred.b = function(but){
    predict(but.mir.lm.b, newdata = data.frame(butralin = but), interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred.b, simplify = T)[3,], lty = 3, col = 2)

legend('topright', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
       legend = c('LC50 - Linear',
                  'LC50 - Exponential',
                  'Slp - linear',
                  '95% CI'), cex = 0.7)  
legend('top', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
       cex = 0.8, bty = 'n', title = 'Observed points')

auc.but.mir.lin0 = function(He){
  e0 = as.numeric(predict(but.mir.lm.e, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(but.mir.lm.b, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
  auc0
}

piM.ghaf_butr.lin_unc = function(He){
  #if(He == 0) piM = 1 else{
    Heu = He/1000
    e = as.numeric(predict(but.mir.lm.e, newdata = data.frame(butralin = Heu), se.fit = TRUE)[1:2])
    b = as.numeric(predict(but.mir.lm.b, newdata = data.frame(butralin = Heu), se.fit = TRUE)[1:2])
  
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piM = auc/auc.but.mir.lin0(0) #gaf.but.aucs.mir[1]
  #  if(piM > 1) piM = 1
  #}
  return(piM)
}  #function to estimate AUC 

auc.but.mir.exp0 = function(He){
  e0 = as.numeric(predict(but.mir.lm.e2, newdata = data.frame(logbutralin = log(He+1)), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(but.mir.lm.b, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
  auc0
}


piM.ghaf_butr.exp_unc = function(He){
  #if(He == 0) piM = 1 else{
    Heu = He/1000
    e = as.numeric(predict(but.mir.lm.e2, newdata = data.frame(logbutralin = log(Heu+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(but.mir.lm.b, newdata = data.frame(butralin = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, e[1], e[2])
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piM = auc/auc.but.mir.exp0(0) #gaf.but.aucs.mir[1]
  #  if(piM > 1) piM = 1
  #}
  return(piM)
}  #function to estimate AUC 

keep.gaf.but.mir = c('piM.ghaf_butr.lin_unc', 'but.mir.lm.e', 'but.mir.lm.b', 'L.3.fx', 'gaf.but.aucs.mir',
                     'piM.ghaf_butr.exp_unc', 'but.mir.lm.e2', 'auc.but.mir.lin0', 'auc.but.mir.exp0')

#plot to test function ########
plot(but.mir$butralin, gaf.but.aucs.mir/gaf.but.aucs.mir[1],
     xlab = 'Butralin (ppm)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0, 9000,100)/1000, sapply(seq(0, 9000, 100), piM.ghaf_butr.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0, 9000, 100)/1000, sapply(seq(0, 9000, 100), piM.ghaf_butr.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)

#glyphosate ###############
mir.gly = subset(mir, chem == 'glyphosate')

gaf.gly.drc.mir = drm(alive/total ~ time_hrs, conc, weights = total, data = mir.gly, type = 'binomial', 
                      fct = LL.3(names = c('b', 'd', 'e'),
                                 fixed = c(NA, 1, NA)))
  summary(gaf.gly.drc.mir)
  plot(gaf.gly.drc.mir) 

plot(mir.ctrl$time_hrs, mir.ctrl$surv, ylim = c(0,1), xlim = c(0,24), pch=16, 
     xlab = 'time(hrs)', ylab = 'prop surviving',
     main = 'Abdel-Ghaffar 2016: Glyphosate toxicity to miracidia')
  lines(time, L.3.fx(lc50 = gaf.ctrl.drc.mir$coefficients[2], slp = gaf.ctrl.drc.mir$coefficients[1], t =time), lty = 2)
  
  for(i in 1:length(unique(mir.gly$conc))){
    points(mir.gly$time_hrs[mir.gly$conc==unique(mir.gly$conc)[i]], 
           mir.gly$surv[mir.gly$conc==unique(mir.gly$conc)[i]], pch=17,
           col = i+1)
    lines(time, L.3.fx(lc50 = gaf.gly.drc.mir$coefficients[i+5], slp = gaf.gly.drc.mir$coefficients[i], t =time), 
          lty = 2, col = i+1)
  }

  legend('topright', title = 'Glyphosate (ppm)', legend = gly.vals, pch = c(16,rep(17,5)),
         col = c(1:6), cex=0.7, bty = 'n')

#Get estimate of miracidia-hrs for each concentration    
  gaf.gly.aucs.mir = as.numeric()
  gaf.gly.aucs.mir[1] = integrate(f = L.3.fx, lc50 = gaf.ctrl.drc.mir$coefficients[2], slp = gaf.ctrl.drc.mir$coefficients[1], 
                                  lower=0, upper=24)[1]$value  
  
  for(j in 1:length(unique(mir.gly$conc))){
    gaf.gly.aucs.mir[j+1] = integrate(f = L.3.fx, lc50 = gaf.gly.drc.mir$coefficients[j+5], slp = gaf.gly.drc.mir$coefficients[j],
                                      lower=0, upper=24)[1]$value  
  }

#Compile L.2 parameters ######## 
gly.mir = data.frame(glyphosate = gly.vals,
                     logglyphosate = log(gly.vals+1),
                     e = c(summary(gaf.ctrl.drc.mir)$coefficients[2,1], summary(gaf.gly.drc.mir)$coefficients[c(6:10),1]),
                     e.se = c(summary(gaf.ctrl.drc.mir)$coefficients[2,2], summary(gaf.gly.drc.mir)$coefficients[c(6:10),2]),
                     b = c(summary(gaf.ctrl.drc.mir)$coefficients[1,1], summary(gaf.gly.drc.mir)$coefficients[c(1:5),1]),
                     b.se = c(summary(gaf.ctrl.drc.mir)$coefficients[1,2], summary(gaf.gly.drc.mir)$coefficients[c(1:5),2]))

plot(gly.mir$glyphosate, gly.mir$e, pch = 16, xlab = 'glyphosate (ppb)', 
     ylab = 'LL.2 Parameters (miracidia)', ylim = c(0,6))

  points(gly.mir$glyphosate, gly.mir$b, pch = 17, col=2)

for(i in 1:length(gly.mir$glyphosate)){
  segments(x0 = gly.mir$glyphosate[i], y0 = gly.mir$e[i] + gly.mir$e.se[i],
           x1 = gly.mir$glyphosate[i], y1 = gly.mir$e[i] - gly.mir$e.se[i])
  segments(x0 = gly.mir$glyphosate[i], y0 = gly.mir$b[i] + gly.mir$b.se[i],
           x1 = gly.mir$glyphosate[i], y1 = gly.mir$b[i] - gly.mir$b.se[i], col=2)
} 

#fit functions to L.2 parameters ##########  
gly.mir.lm.e = lm(e ~ glyphosate, weights = e.se^-1, data = gly.mir) #linear response of LC50
  gly.mir.pred = function(gly){
    predict(gly.mir.lm.e, newdata = data.frame(glyphosate = gly), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred, simplify = T)[3,], lty = 3)

gly.mir.lm.e2 = lm(e ~ logglyphosate, weights = e.se^-1, data = gly.mir) #log-linear response of LC50
  gly.mir.pred2 = function(gly){
    predict(gly.mir.lm.e2, newdata = data.frame(logglyphosate = log(gly+1)), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred2, simplify = T)[3,], lty = 3, col=3)

AIC(gly.mir.lm.e, gly.mir.lm.e2)  #Exponential is a better fit    

gly.mir.lm.b = lm(b ~ glyphosate, weights = b.se^-1, data = gly.mir)   
  gly.mir.pred.b = function(gly){
    predict(gly.mir.lm.b, newdata = data.frame(glyphosate = gly), interval = 'confidence', level = 0.95)
  }

    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topleft', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

auc.gly.mir.lin0 = function(He){
    e0 = as.numeric(predict(gly.mir.lm.e, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(gly.mir.lm.b, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
    
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
    auc0
  }
  
piM.ghaf_gly.lin_unc = function(He){
  #if(He == 0) piM = 1 else{
    Heu = He/1000
    e = as.numeric(predict(gly.mir.lm.e, newdata = data.frame(glyphosate = Heu), se.fit = TRUE)[1:2])
    b = as.numeric(predict(gly.mir.lm.b, newdata = data.frame(glyphosate = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    #print(e.use)
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    #print(b.use)
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piM = auc/auc.gly.mir.lin0(0) #gaf.gly.aucs.mir[1]
  #  if(piM > 1) piM = 1
  #}
  return(piM)
}  #function to estimate AUC 

auc.gly.mir.exp0 = function(He){
  e0 = as.numeric(predict(gly.mir.lm.e2, newdata = data.frame(logglyphosate = log(He+1)), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(gly.mir.lm.b, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
  auc0
}


piM.ghaf_gly.exp_unc = function(He){
  #if(He == 0) piM = 1 else{
    Heu = He/1000
    e = as.numeric(predict(gly.mir.lm.e2, newdata = data.frame(logglyphosate = log(Heu+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(gly.mir.lm.b, newdata = data.frame(glyphosate = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    #print(e.use)
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, e[1], e[2])
    #print(b.use)
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piM = auc/auc.gly.mir.exp0(0) #gaf.gly.aucs.mir[1]
  #  if(piM > 1) piM = 1
  #}
  return(piM)
}  #function to estimate AUC 

keep.gaf.gly.mir = c('piM.ghaf_gly.lin_unc', 'gly.mir.lm.e', 'gly.mir.lm.b', 'L.3.fx', 'gaf.gly.aucs.mir',
                     'piM.ghaf_gly.exp_unc', 'gly.mir.lm.e2', 'auc.gly.mir.lin0', 'auc.gly.mir.exp0')

#plot to test function ########
plot(gly.mir$glyphosate, gaf.gly.aucs.mir/gaf.gly.aucs.mir[1],
     xlab = 'Glyphosate (ppb)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0,25500,100)/1000, sapply(seq(0,25500,100), piM.ghaf_gly.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,25500,100)/1000, sapply(seq(0,25500,100), piM.ghaf_gly.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)

#pendimethalin ###############
mir.pen = subset(mir, chem == 'pendimethalin')

gaf.pen.drc.mir = drm(alive/total ~ time_hrs, conc, weights = total, data = mir.pen, type = 'binomial', 
                      fct = LL.3(names = c('b', 'd', 'e'),
                                 fixed = c(NA, 1, NA)))
  summary(gaf.pen.drc.mir)
  plot(gaf.pen.drc.mir) 

plot(mir.ctrl$time_hrs, mir.ctrl$surv, ylim = c(0,1), xlim = c(0,24), pch=16, 
     xlab = 'time(hrs)', ylab = 'prop surviving',
     main = 'Abdel-Ghaffar 2016: Pendimethalin toxicity to miracidia')
  lines(time, L.3.fx(lc50 = gaf.ctrl.drc.mir$coefficients[2], slp = gaf.ctrl.drc.mir$coefficients[1], t =time), lty = 2)

  for(i in 1:length(unique(mir.pen$conc))){
    points(mir.pen$time_hrs[mir.pen$conc==unique(mir.pen$conc)[i]], 
           mir.pen$surv[mir.pen$conc==unique(mir.pen$conc)[i]], pch=17,
           col = i+1)
    lines(time, L.3.fx(lc50 = gaf.pen.drc.mir$coefficients[i+5], slp = gaf.pen.drc.mir$coefficients[i], t =time), 
          lty = 2, col = i+1)
  }

  legend('topright', title = 'pendimethalin (ppm)', legend = pen.vals, pch = c(16,rep(17,5)),
         col = c(1:6), cex=0.7, bty = 'n')

#Get estimate of miracidia-hrs for each concentration    
gaf.pen.aucs.mir = as.numeric()
  gaf.pen.aucs.mir[1] = integrate(f = L.3.fx, lc50 = gaf.ctrl.drc.mir$coefficients[2], slp = gaf.ctrl.drc.mir$coefficients[1], 
                                  lower=0, upper=24)[1]$value  
  
  for(j in 1:length(unique(mir.pen$conc))){
    gaf.pen.aucs.mir[j+1] = integrate(f = L.3.fx, lc50 = gaf.pen.drc.mir$coefficients[j+5], slp = gaf.pen.drc.mir$coefficients[j],
                                      lower=0, upper=24)[1]$value  
  }
  
#Compile L.2 parameters ######## 
pen.mir = data.frame(pendimethalin = pen.vals,
                     logpendimethalin = log(pen.vals+1),
                     e = c(summary(gaf.ctrl.drc.mir)$coefficients[2,1], summary(gaf.pen.drc.mir)$coefficients[c(6:10),1]),
                     e.se = c(summary(gaf.ctrl.drc.mir)$coefficients[2,2], summary(gaf.pen.drc.mir)$coefficients[c(6:10),2]),
                     b = c(summary(gaf.ctrl.drc.mir)$coefficients[1,1], summary(gaf.pen.drc.mir)$coefficients[c(1:5),1]),
                     b.se = c(summary(gaf.ctrl.drc.mir)$coefficients[1,2], summary(gaf.pen.drc.mir)$coefficients[c(1:5),2]))

plot(pen.mir$pendimethalin, pen.mir$e, pch = 16, xlab = 'Pendimethalin (ppb)', 
     ylab = 'L.2 Parameters (miracidia)', ylim = c(0,5))

  points(pen.mir$pendimethalin, pen.mir$b, pch = 17, col=2)

  for(i in 1:length(pen.mir$pendimethalin)){
    segments(x0 = pen.mir$pendimethalin[i], y0 = pen.mir$e[i] + pen.mir$e.se[i],
             x1 = pen.mir$pendimethalin[i], y1 = pen.mir$e[i] - pen.mir$e.se[i])
    segments(x0 = pen.mir$pendimethalin[i], y0 = pen.mir$b[i] + pen.mir$b.se[i],
             x1 = pen.mir$pendimethalin[i], y1 = pen.mir$b[i] - pen.mir$b.se[i], col=2)
  } 

#fit functions to L.2 parameters ##########  
pen.mir.lm.e = lm(e ~ pendimethalin, weights = e.se^-1, data = pen.mir) #linear response of LC50
  pen.mir.pred = function(pen){
    predict(pen.mir.lm.e, newdata = data.frame(pendimethalin = pen), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred, simplify = T)[3,], lty = 3)

pen.mir.lm.e2 = lm(e ~ logpendimethalin, weights = e.se^-1, data = pen.mir) #log-linear response of LC50
  pen.mir.pred2 = function(pen){
    predict(pen.mir.lm.e2, newdata = data.frame(logpendimethalin = log(pen+1)), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred2, simplify = T)[3,], lty = 3, col=3)

AIC(pen.mir.lm.e, pen.mir.lm.e2)  #Exponential is a better fit    

pen.mir.lm.b = lm(b ~ pendimethalin, weights = b.se^-1, data = pen.mir)   
  pen.mir.pred.b = function(pen){
    predict(pen.mir.lm.b, newdata = data.frame(pendimethalin = pen), interval = 'confidence', level = 0.95)
  }

    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topleft', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

auc.pen.mir.lin0 = function(He){
    e0 = as.numeric(predict(pen.mir.lm.e, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(pen.mir.lm.b, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
    
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
    auc0
  }
  
  
piM.ghaf_pen.lin_unc = function(He){
  #if(He == 0) piM = 1 else{
    Heu = He/1000
    e = as.numeric(predict(pen.mir.lm.e, newdata = data.frame(pendimethalin = Heu), se.fit = TRUE)[1:2])
    b = as.numeric(predict(pen.mir.lm.b, newdata = data.frame(pendimethalin = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piM = auc/auc.pen.mir.lin0(0) #gaf.pen.aucs.mir[1]
  #  if(piM > 1) piM = 1
  #}
  return(piM)
}  #function to estimate AUC 

auc.pen.mir.exp0 = function(He){
  e0 = as.numeric(predict(pen.mir.lm.e2, newdata = data.frame(logpendimethalin = log(He+1)), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(pen.mir.lm.b, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
  auc0
}

piM.ghaf_pen.exp_unc = function(He){
  #if(He == 0) piM = 1 else{
    Heu = He/1000
    e = as.numeric(predict(pen.mir.lm.e2, newdata = data.frame(logpendimethalin = log(Heu+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(pen.mir.lm.b, newdata = data.frame(pendimethalin = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    #print(e.use)
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    #print(b.use)
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piM = auc/auc.pen.mir.exp0(0) #gaf.pen.aucs.mir[1]
  # if(piM > 1) piM = 1
  #}
  return(piM)
}  #function to estimate AUC 

keep.gaf.pen.mir = c('piM.ghaf_pen.lin_unc', 'pen.mir.lm.e', 'pen.mir.lm.b', 'L.3.fx', 'gaf.pen.aucs.mir',
                     'piM.ghaf_pen.exp_unc', 'pen.mir.lm.e2', 'auc.pen.mir.lin0', 'auc.pen.mir.exp0')

#plot to test function ########
plot(pen.mir$pendimethalin, gaf.pen.aucs.mir/gaf.pen.aucs.mir[1],
     xlab = 'pendimethalin (ppm)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0,4000,10)/1000, sapply(seq(0,4000,10), piM.ghaf_pen.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,4000,10)/1000, sapply(seq(0,4000,10), piM.ghaf_pen.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)

#final keep vector #############
keep.gaf.piM = c(keep.gaf.but.mir, keep.gaf.gly.mir, keep.gaf.pen.mir)  