#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
#Data extraction and model fitting to Tantawy 2002 data
require(drc)
L.3.fx = function(t, lc50 = lc50, slp = slp){
  1 / (1+exp(slp*log(t / lc50)))
}
time = c(0:25)


#miracidia toxicity ############
tant<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/Tantawy2002.csv')
  tant$total = 100
    mir<-subset(tant, larv == 'miracidia')
    mir.but = subset(mir, chem == 'butachlor')
    mir.fpb = subset(mir, chem == 'fluazifop-p-butyl')
  
#butachlor toxicity to miracidia  ###############
tant.piM.but<-drm(surv/total ~ time_hrs, conc, weights = total,  
                  data = mir.but, type = 'binomial', 
                  fct = LL.4(names = c('b', 'c', 'd', 'e'),
                             fixed = c(NA, 0, 1, NA)))
  summary(tant.piM.but)
  plot(tant.piM.but) 

plot(mir.but$time_hrs[mir.but$conc==0], mir.but$surv[mir.but$conc==0]/100, 
     pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tant.piM.but, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in 2:length(unique(mir.but$conc))){
    points(mir.but$time_hrs[mir.but$conc==unique(mir.but$conc)[i]], 
           mir.but$surv[mir.but$conc==unique(mir.but$conc)[i]]/100, pch=16,
           col = i)
    lines(time, predict(tant.piM.but, 
                        data.frame(time_hrs=time, conc = unique(mir.but$conc)[i])),
          lty = 2, col = i)
  }

title('butachlor toxicity to miracidia')
legend('topright', legend = c('control', 650,1500,4500,6500), 
       pch = c(17,16,16,16,16), col = c(1:5), cex=0.8, bty = 'n')  

#Get estimate of cercariae-hrs for each concentration    
tant.but.piM.aucs = as.numeric()

for(j in 1:length(unique(mir.but$conc))){
  fx = function(t){
    predict(tant.piM.but, newdata = data.frame(time_hrs = t, conc = unique(mir.but$conc)[j]))
  }
  tant.but.piM.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
}

#Compile LL.2 data for functional responses #############
mir.but = data.frame(e = summary(tant.piM.but)$coefficients[c(6:10),1],
                     e.se = summary(tant.piM.but)$coefficients[c(6:10),2],
                     b = summary(tant.piM.but)$coefficients[c(1:5),1],
                     b.se = summary(tant.piM.but)$coefficients[c(1:5),2],
                     but = c(0,650,1500,4500,6500),
                     logbut = log(c(0,650,1500,4500,6500)+1))

plot(mir.but$but, mir.but$e, pch = 16, xlab = 'butachlor (ppm)', ylab = 'LL.2 Parameters (miracidia',
     ylim = c(0,20))

  points(mir.but$but+50, mir.but$b, pch = 17, col=2)

  for(i in 1:length(mir.but$but)){
    segments(x0 = mir.but$but[i], y0 = mir.but$e[i] + mir.but$e.se[i],
             x1 = mir.but$but[i], y1 = mir.but$e[i] - mir.but$e.se[i])
    segments(x0 = mir.but$but[i]+50, y0 = mir.but$b[i] + mir.but$b.se[i],
             x1 = mir.but$but[i]+50, y1 = mir.but$b[i] - mir.but$b.se[i], col=2)
  } 
  
#fit models to LL.2 parameters across concentration ########
el.but.mir = lm(e ~ but, weights = e.se^-1, data = mir.but) #linear response of LC50
  el.pred.mir = function(but){
    predict(el.but.mir, newdata = data.frame(but = but), 
            interval = 'confidence', level = 0.95)
  }

    lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred.mir, simplify = T)[1,], lty = 2)
    lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred.mir, simplify = T)[2,], lty = 3)
    lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred.mir, simplify = T)[3,], lty = 3)

el.but2.mir = lm(e ~ logbut, weights = e.se^-1, data = mir.but) #log-linear response of LC50
  el.pred2.mir = function(but){
    predict(el.but2.mir, newdata = data.frame(logbut = log(but+1)), 
            interval = 'confidence', level = 0.95)
  }

  lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred2.mir, simplify = T)[1,], lty = 2, col=3)
  lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred2.mir, simplify = T)[2,], lty = 3, col=3)
  lines(seq(0,7000,7), sapply(seq(0,7000,7), el.pred2.mir, simplify = T)[3,], lty = 3, col=3)

  AIC(el.but.mir, el.but2.mir)  #exponential is a better fit    

bl.but.mir = lm(b ~ but, weights = b.se^-1, data = mir.but)   
  bl.pred.mir = function(but){
    predict(bl.but.mir, newdata = data.frame(but = but), interval = 'confidence', level = 0.95)
  }

  lines(seq(0,7000,7), sapply(seq(0,7000,7), bl.pred.mir, simplify = T)[1,], lty = 2, col = 2)
  lines(seq(0,7000,7), sapply(seq(0,7000,7), bl.pred.mir, simplify = T)[2,], lty = 3, col = 2)
  lines(seq(0,7000,7), sapply(seq(0,7000,7), bl.pred.mir, simplify = T)[3,], lty = 3, col = 2)
  
  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

piM.tant02_but.lin_unc = function(He){
  if(He == 0) piM = 1 else{
    e = as.numeric(predict(el.but.mir, newdata = data.frame(but = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bl.but.mir, newdata = data.frame(but = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, e[1], e[2])
    
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24)[1]$value
    piM = auc/tant.but.piM.aucs[1]
    if(piM > 1) piM = 1
  }
  return(piM)
}  #function to estimate AUC 

piM.tant02_but.exp_unc = function(He){
  if(He == 0) piM = 1 else{
  e = as.numeric(predict(el.but2.mir, newdata = data.frame(logbut = log(He+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bl.but.mir, newdata = data.frame(but = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
  while(b.use < 0) b.use = rnorm(1, e[1], e[2])
  
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24)[1]$value
  piM = auc/tant.but.piM.aucs[1]
  if(piM > 1) piM = 1
  }
  return(piM)
}  #function to estimate AUC 

#plot sample outputs compared to observed points ##########
plot(mir.but$but, tant.but.piM.aucs/tant.but.piM.aucs[1],
     xlab = 'Butachlor (ppb)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0, 6500,25), sapply(seq(0, 6500,25), piM.tant02_but.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0, 6500,25), sapply(seq(0, 6500,25), piM.tant02_but.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)
  
#keep vector for butachlor toxicity to miracidia
  keep.tant.but.piM = c('piM.tant02_but.lin_unc', 'piM.tant02_but.exp_unc', 
                        'el.but.mir', 'el.but2.mir', 'bl.but.mir', 'L.3.fx', 'tant.but.piM.aucs')
    
#fluazifop-p-butyl toxicity to miracidia ###########
tant.piM.fpb<-drm(surv/total ~ time_hrs, conc, weights = total,  
                  data = mir.fpb, type = 'binomial', 
                  fct = LL.4(names = c('b', 'c', 'd', 'e'),
                             fixed = c(NA, 0, 1, NA)))
  summary(tant.piM.fpb)
  plot(tant.piM.fpb) 

plot(mir.fpb$time_hrs[mir.fpb$conc==0], mir.fpb$surv[mir.fpb$conc==0]/100, 
     pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tant.piM.fpb, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in 2:length(unique(mir.fpb$conc))){
    points(mir.fpb$time_hrs[mir.fpb$conc==unique(mir.fpb$conc)[i]], 
           mir.fpb$surv[mir.fpb$conc==unique(mir.fpb$conc)[i]]/100, pch=16,
           col = i)
    lines(time, predict(tant.piM.fpb, 
                        data.frame(time_hrs=time, conc = unique(mir.fpb$conc)[i])),
          lty = 2, col = i)
  }

  title('fluazifop-p-butyl toxicity to miracidia')
  legend('topright', legend = c('control', 1760,4500,9000,17600), 
         pch = c(17,16,16,16,16), col = c(1:5), cex=0.8, bty = 'n')  

#Get estimate of cercariae-hrs for each concentration    
  tant.fpb.piM.aucs = as.numeric()

  for(j in 1:length(unique(mir.fpb$conc))){
    fx = function(t){
      predict(tant.piM.fpb, newdata = data.frame(time_hrs = t, conc = unique(mir.fpb$conc)[j]))
    }
    tant.fpb.piM.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }

#Compile LL.2 data for functional responses #############
mir.fpb = data.frame(e = summary(tant.piM.fpb)$coefficients[c(6:10),1],
                     e.se = summary(tant.piM.fpb)$coefficients[c(6:10),2],
                     b = summary(tant.piM.fpb)$coefficients[c(1:5),1],
                     b.se = summary(tant.piM.fpb)$coefficients[c(1:5),2],
                     fpb = c(0, 1760,4500,9000,17600),
                     logfpb = log(c(0, 1760,4500,9000,17600)+1))
  
  plot(mir.fpb$fpb, mir.fpb$e, pch = 16, xlab = 'butachlor (ppm)', ylab = 'LL.2 Parameters (miracidia',
       ylim = c(0,20))
  
  points(mir.fpb$fpb+50, mir.fpb$b, pch = 17, col=2)
  
  for(i in 1:length(mir.fpb$fpb)){
    segments(x0 = mir.fpb$fpb[i], y0 = mir.fpb$e[i] + mir.fpb$e.se[i],
             x1 = mir.fpb$fpb[i], y1 = mir.fpb$e[i] - mir.fpb$e.se[i])
    segments(x0 = mir.fpb$fpb[i]+50, y0 = mir.fpb$b[i] + mir.fpb$b.se[i],
             x1 = mir.fpb$fpb[i]+50, y1 = mir.fpb$b[i] - mir.fpb$b.se[i], col=2)
  } 
  
#fit models to LL.2 parameters across concentration ########
el.fpb.mir = lm(e ~ fpb, weights = e.se^-1, data = mir.fpb) #linear response of LC50
  el.pred.mir.fpb = function(fpb){
    predict(el.fpb.mir, newdata = data.frame(fpb = fpb), 
            interval = 'confidence', level = 0.95)
  }

  lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred.mir.fpb, simplify = T)[1,], lty = 2)
  lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred.mir.fpb, simplify = T)[2,], lty = 3)
  lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred.mir.fpb, simplify = T)[3,], lty = 3)

el.fpb2.mir = lm(e ~ logfpb, weights = e.se^-1, data = mir.fpb) #log-linear response of LC50
  el.pred2.mir.fpb = function(fpb){
    predict(el.fpb2.mir, newdata = data.frame(logfpb = log(fpb+1)), 
            interval = 'confidence', level = 0.95)
  }

  lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred2.mir.fpb, simplify = T)[1,], lty = 2, col=3)
  lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred2.mir.fpb, simplify = T)[2,], lty = 3, col=3)
  lines(seq(0,18000,100), sapply(seq(0,18000,100), el.pred2.mir.fpb, simplify = T)[3,], lty = 3, col=3)

    AIC(el.fpb.mir, el.fpb2.mir)  #Linear is a slightly better fit    

bl.fpb.mir = lm(b ~ fpb, weights = b.se^-1, data = mir.fpb)   
  bl.pred.mir.fpb = function(fpb){
    predict(bl.fpb.mir, newdata = data.frame(fpb = fpb), interval = 'confidence', level = 0.95)
  }

  lines(seq(0,18000,100), sapply(seq(0,18000,100), bl.pred.mir.fpb, simplify = T)[1,], lty = 2, col = 2)
  lines(seq(0,18000,100), sapply(seq(0,18000,100), bl.pred.mir.fpb, simplify = T)[2,], lty = 3, col = 2)
  lines(seq(0,18000,100), sapply(seq(0,18000,100), bl.pred.mir.fpb, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

piM.tant02_fpb.lin_unc = function(He){
  if(He == 0) piM = 1 else{
  e = as.numeric(predict(el.fpb.mir, newdata = data.frame(fpb = He), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bl.fpb.mir, newdata = data.frame(fpb = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
  while(b.use < 0) b.use = rnorm(1, e[1], e[2])
  
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24)[1]$value
  piM = auc/tant.fpb.piM.aucs[1]
  if(piM > 1) piM = 1
  }
  return(piM)
}  #function to estimate AUC 

piM.tant02_fpb.exp_unc = function(He){
  if(He == 0) piM = 1 else{
  e = as.numeric(predict(el.fpb2.mir, newdata = data.frame(logfpb = log(He+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bl.fpb.mir, newdata = data.frame(fpb = He), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
  while(b.use < 0) b.use = rnorm(1, e[1], e[2])
  
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24)[1]$value
  piM = auc/tant.fpb.piM.aucs[1]
  if(piM > 1) piM = 1
  }
  return(piM)
}  #function to estimate AUC 

#plot sample output compared to observed AUCs ##########
plot(mir.fpb$fpb, tant.fpb.piM.aucs/tant.fpb.piM.aucs[1],
     xlab = 'Fluazifop-p-butyl (ppb)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))
  points(seq(0, 18000,50), sapply(seq(0, 18000,50), piM.tant02_fpb.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0, 18000,50), sapply(seq(0, 18000,50), piM.tant02_fpb.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)
  
#keep vector for fluazifop-p-butyl toxicity to miracidia
  keep.tant.fpb.piM = c('piM.tant02_fpb.lin_unc', 'piM.tant02_fpb.exp_unc', 
                        'el.fpb.mir', 'el.fpb2.mir', 'bl.fpb.mir', 'L.3.fx', 'tant.fpb.piM.aucs')
  
#keep vector ##########
  keep.tantawy.piM = c(keep.tant.but.piM, keep.tant.fpb.piM)    
  
  
