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

#cercarial toxicity ############
tant<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/Tantawy2002.csv')
tant$total = 100
tant$conc = tant$conc/1000
cerc<-subset(tant, larv == 'cercariae')
  cerc.but = subset(cerc, chem == 'butachlor')
  cerc.fpb = subset(cerc, chem == 'fluazifop-p-butyl')

#butachlor toxicity to cercariae  ###############
tant.piC.but<-drm(surv/total ~ time_hrs, conc, weights = total,  data = cerc.but, type = 'binomial', 
                  fct = LL.3(names = c('b', 'd', 'e'),
                            fixed = c(NA, 1, NA)))
  summary(tant.piC.but)
  plot(tant.piC.but) 
  
  plot(cerc.but$time_hrs[cerc.but$conc==0], cerc.but$surv[cerc.but$conc==0]/100, 
       pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tant.piC.but, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in 2:length(unique(cerc.but$conc))){
    points(cerc.but$time_hrs[cerc.but$conc==unique(cerc.but$conc)[i]], 
           cerc.but$surv[cerc.but$conc==unique(cerc.but$conc)[i]]/100, pch=16,
           col = i)
    lines(time, predict(tant.piC.but, 
                        data.frame(time_hrs=time, conc = unique(cerc.but$conc)[i])),
          lty = 2, col = i)
  }
  
  title('butachlor toxicity to cercariae')
  legend('topright', legend = c('control', .650,1.500,4.500,6.500), title = 'butachlor (ppm)',
         pch = c(17,16,16,16,16), col = c(1:5), cex=0.8, bty = 'n')  

#Get estimate of cercariae-hrs for each concentration    
tant.but.pic.aucs = as.numeric()

  for(j in 1:length(unique(cerc.but$conc))){
    fx = function(t){
      predict(tant.piC.but, newdata = data.frame(time_hrs = t, conc = unique(cerc.but$conc)[j]))
    }
    tant.but.pic.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }

#compile butachlor data for function ###############
cerc.but = data.frame(e = summary(tant.piC.but)$coefficients[c(6:10),1],
                      e.se = summary(tant.piC.but)$coefficients[c(6:10),2],
                      b = summary(tant.piC.but)$coefficients[c(1:5),1],
                      b.se = summary(tant.piC.but)$coefficients[c(1:5),2],
                      but = unique(cerc.but$conc),
                      logbut = log(unique(cerc.but$conc)+1))
    
plot(cerc.but$but, cerc.but$e, pch = 16, xlab = 'butachlor (ppm)', ylab = 'LL.2 Parameters',
     ylim = c(0,20))

  points(cerc.but$but, cerc.but$b, pch = 17, col=2)
  
    for(i in 1:length(cerc.but$but)){
      segments(x0 = cerc.but$but[i], y0 = cerc.but$e[i] + cerc.but$e.se[i],
               x1 = cerc.but$but[i], y1 = cerc.but$e[i] - cerc.but$e.se[i])
      segments(x0 = cerc.but$but[i], y0 = cerc.but$b[i] + cerc.but$b.se[i],
               x1 = cerc.but$but[i], y1 = cerc.but$b[i] - cerc.but$b.se[i], col=2)
    }    
  
#fit models to L.2 parameters across concentration ########
  el.but = lm(e ~ but, weights = e.se^-1, data = cerc.but) #linear response of LC50
    el.pred = function(but){
      predict(el.but, newdata = data.frame(but = but), 
              interval = 'confidence', level = 0.95)
    }
    
    lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred)[1,], lty = 2)
    lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred)[2,], lty = 3)
    lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred)[3,], lty = 3)

  el.but2 = lm(e ~ logbut, weights = e.se^-1, data = cerc.but) #log-linear response of LC50
    el.pred2 = function(but){
      predict(el.but2, newdata = data.frame(logbut = log(but+1)), 
              interval = 'confidence', level = 0.95)
    }
    
  lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred2)[1,], lty = 2, col=3)
  lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred2)[2,], lty = 3, col=3)
  lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred2)[3,], lty = 3, col=3)
    
  AIC(el.but, el.but2)  #exponential is a better fit    
  
  bl.but = lm(b ~ but, weights = b.se^-1, data = cerc.but)   
    bl.pred = function(but){
      predict(bl.but, newdata = data.frame(but = but), interval = 'confidence', level = 0.95)
    }
  
    lines(seq(0,7,.07), sapply(seq(0,7,.07), bl.pred)[1,], lty = 2, col = 2)
    lines(seq(0,7,.07), sapply(seq(0,7,.07), bl.pred)[2,], lty = 3, col = 2)
    lines(seq(0,7,.07), sapply(seq(0,7,.07), bl.pred)[3,], lty = 3, col = 2)
    
    legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
           legend = c('LC50 - Linear',
                      'LC50 - Exponential',
                      'Slp - linear',
                      '95% CI'), cex = 0.7)  
    legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
           cex = 0.7, bty = 'n', title = 'Observed points')
     
piC.tant02_but.lin_auc0 = function(He){
  Heu = He/1000
  e0 = as.numeric(predict(el.but, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(bl.but, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                  stop.on.error = FALSE)[1]$value
  return(auc0)
}

piC.tant02_but.lin_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(el.but, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bl.but, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
        
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                    stop.on.error = FALSE)[1]$value
    piC = auc/piC.tant02_but.lin_auc0(0)  #tant.but.pic.aucs[1]
  
    return(piC)
} 

piC.tant02_but.exp_auc0 = function(He){
  Heu = He/1000
  e0 = as.numeric(predict(el.but2, newdata = data.frame(logbut = log(Heu+1)), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(bl.but, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                   stop.on.error = FALSE)[1]$value
  return(auc0)
}

piC.tant02_but.exp_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(el.but2, newdata = data.frame(logbut = log(Heu+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bl.but, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                    stop.on.error = FALSE)[1]$value
    piC = auc/piC.tant02_but.exp_auc0(0)  #tant.but.pic.aucs[1]
  
  return(piC)
}  #Parameter estimate

#plot sample output compared to observed AUCs ##########
  plot(cerc.but$but, tant.but.pic.aucs / tant.but.pic.aucs[1],
        xlab = 'Butachlor (ppb)', ylab = 'relative AUC (cercariae-hours)',
       pch = 16, ylim = c(0,1))
    points(seq(0, 6500,50)/1000, sapply(seq(0, 6500,50), piC.tant02_but.lin_unc),
           pch = 5, col = 4, cex = 0.5)
    points(seq(0, 6500,50)/1000, sapply(seq(0, 6500,50), piC.tant02_but.exp_unc),
           pch = 5, col = 2, cex = 0.5)
    
  keep.tant.but.pic = c('piC.tant02_but.lin_unc', 'piC.tant02_but.lin_auc0', 'piC.tant02_but.exp_unc', 'piC.tant02_but.exp_auc0', 'bl.but',
                        'el.but2', 'el.but', 'L.3.fx')      
    
#fluazifop-p-butyl toxicity to cercariae ###########
  tant.piC.fpb<-drm(surv/total ~ time_hrs, conc, weights = total,  data = cerc.fpb, type = 'binomial', 
                    fct = LL.3(names = c('b', 'd', 'e'),
                              fixed = c(NA, 1, NA)))
  summary(tant.piC.fpb)
  plot(tant.piC.fpb) 
  
  plot(cerc.fpb$time_hrs[cerc.fpb$conc==0], cerc.fpb$surv[cerc.fpb$conc==0]/100, 
       pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tant.piC.fpb, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in 2:length(unique(cerc.fpb$conc))){
    points(cerc.fpb$time_hrs[cerc.fpb$conc==unique(cerc.fpb$conc)[i]], 
           cerc.fpb$surv[cerc.fpb$conc==unique(cerc.fpb$conc)[i]]/100, pch=16,
           col = i)
    lines(time, predict(tant.piC.fpb, 
                        data.frame(time_hrs=time, conc = unique(cerc.fpb$conc)[i])),
          lty = 2, col = i)
  }
  
  title('fluazifop-p-butyl toxicity to cercariae')
  legend('topright', legend = c('control', 1.760,4.500,9.000,17.600), 
         pch = c(17,16,16,16,16), col = c(1:5), cex=0.8, bty = 'n')  
  
#Get estimate of cercariae-hrs for each concentration    
  tant.fpb.pic.aucs = as.numeric()
  
  for(j in 1:length(unique(cerc.fpb$conc))){
    fx = function(t){
      predict(tant.piC.fpb, newdata = data.frame(time_hrs = t, conc = unique(cerc.fpb$conc)[j]))
    }
    tant.fpb.pic.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }
  
  
#compile fluazifop-p-butyl data for function ###############
  cerc.fpb = data.frame(e = summary(tant.piC.fpb)$coefficients[c(6:10),1],
                        e.se = summary(tant.piC.fpb)$coefficients[c(6:10),2],
                        b = summary(tant.piC.fpb)$coefficients[c(1:5),1],
                        b.se = summary(tant.piC.fpb)$coefficients[c(1:5),2],
                        fpb = unique(cerc.fpb$conc),
                        logfpb = log(unique(cerc.fpb$conc)+1))
  
  plot(cerc.fpb$fpb, cerc.fpb$e, pch = 16, xlab = 'fluazifop-p-butyl (ppm)', 
       ylab = 'L.2 Parameters', ylim = c(0,20))
  
  points(cerc.fpb$fpb, cerc.fpb$b, pch = 17, col=2)
  
  for(i in 1:length(cerc.fpb$fpb)){
    segments(x0 = cerc.fpb$fpb[i], y0 = cerc.fpb$e[i] + cerc.fpb$e.se[i],
             x1 = cerc.fpb$fpb[i], y1 = cerc.fpb$e[i] - cerc.fpb$e.se[i])
    segments(x0 = cerc.fpb$fpb[i]+50, y0 = cerc.fpb$b[i] + cerc.fpb$b.se[i],
             x1 = cerc.fpb$fpb[i]+50, y1 = cerc.fpb$b[i] - cerc.fpb$b.se[i], col=2)
  }    
  
#fit models to L.2 parameters across concentration ########
el.fpb = lm(e ~ fpb, weights = e.se^-1, data = cerc.fpb) #linear response of LC50
  el.pred.fpb = function(fpb){
    predict(el.fpb, newdata = data.frame(fpb = fpb), 
            interval = 'confidence', level = 0.95)
  }
    
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.fpb)[1,], lty = 2)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.fpb)[2,], lty = 3)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.fpb)[3,], lty = 3)
    
el.fpb2 = lm(e ~ logfpb, weights = e.se^-1, data = cerc.fpb) #log-linear response of LC50
  el.pred.fpb2 = function(fpb){
    predict(el.fpb2, newdata = data.frame(logfpb = log(fpb+1)), 
            interval = 'confidence', level = 0.95)
  }
    
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.fpb2)[1,], lty = 2, col=3)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.fpb2)[2,], lty = 3, col=3)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.fpb2)[3,], lty = 3, col=3)
    
    AIC(el.fpb, el.fpb2)  #exponential is a better fit    
    
bl.fpb = lm(b ~ fpb, weights = b.se^-1, data = cerc.fpb)   
  bl.pred.fpb = function(fpb){
    predict(bl.fpb, newdata = data.frame(fpb = fpb), interval = 'confidence', level = 0.95)
  }
    
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), bl.pred.fpb)[1,], lty = 2, col = 2)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), bl.pred.fpb)[2,], lty = 3, col = 2)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), bl.pred.fpb)[3,], lty = 3, col = 2)
    
    legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
           legend = c('LC50 - Linear',
                      'LC50 - Exponential',
                      'Slp - linear',
                      '95% CI'), cex = 0.7)  
    legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
           cex = 0.8, bty = 'n', title = 'Observed points')
 

piC.tant02_fpb.exp_auc0 = function(He){
    Heu = He/1000
    e0 = as.numeric(predict(el.fpb2, newdata = data.frame(logfpb = log(Heu+1)), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(bl.fpb, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], bo[2])
    
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value
    
  return(auc0)
}  #Parameter estimate with exponential function

piC.tant02_fpb.exp_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(el.fpb2, newdata = data.frame(logfpb = log(Heu+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bl.fpb, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24)[1]$value
    piC = auc/piC.tant02_fpb.exp_auc0(0)  #tant.fpb.pic.aucs[1]

  return(piC)
}  #Parameter estimate with exponential function
  
piC.tant02_fpb.lin_auc0 = function(He){
  Heu = He/1000
  e0 = as.numeric(predict(el.fpb, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(bl.fpb, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], bo[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                   stop.on.error = FALSE)[1]$value
  
  return(auc0)
}  #Parameter estimate with exponential function

piC.tant02_fpb.lin_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(el.fpb, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bl.fpb, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                    stop.on.error = FALSE)[1]$value
    piC = auc/piC.tant02_fpb.lin_auc0(0)  #tant.fpb.pic.aucs[1]
  
   return(piC)
}  #Parameter estimate with linear function
    
#plot sample output compared to observed AUCs ##########
  plot(cerc.fpb$fpb, tant.fpb.pic.aucs/tant.fpb.pic.aucs[1],
        xlab = 'Fluazifop-p-butyl (ppb)', ylab = 'relative AUC (cercariae-hours)',
        pch = 16, ylim = c(0,1))
    points(seq(0, 18000,100)/1000, sapply(seq(0, 18000,100), piC.tant02_fpb.lin_unc),
           pch = 5, col = 4, cex = 0.5)
    points(seq(0, 18000,100)/1000, sapply(seq(0, 18000,100), piC.tant02_fpb.exp_unc),
           pch = 5, col = 2, cex = 0.5)
    
#keep vector ########    
keep.tant.fpb.pic = c('piC.tant02_fpb.lin_unc', 'piC.tant02_fpb.lin_auc0', 'piC.tant02_fpb.exp_unc', 'piC.tant02_fpb.exp_auc0', 'bl.fpb',
                      'el.fpb2', 'el.fpb', 'L.3.fx') 
    
keep.tantawy.piC = c(keep.tant.but.pic, keep.tant.fpb.pic)    
    