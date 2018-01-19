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
  1 / (1+exp(slp*log(t / lc50)))
}
  time = seq(0,25,0.1)

#miracidial mortality (S. mansoni) from Tchounwou 1991 ####################################
mir2 = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Miracidial Mortality/tchounwou91_urea_ammP_miracidial_mortality.csv')
  mir2$conc = mir2$conc/1000
  mir.amm = subset(mir2, chem == 'amm_sulph')
  mir.ure = subset(mir2, chem == 'urea')

#Ammonium Sulphate concentration ############
tch91.piM.amm<-drm(alive/total ~ time_hrs, conc, weights = total, data = mir.amm, type = 'binomial', 
                   fct = LL.3(names = c('b', 'd', 'e'),
                              fixed = c(NA, 1, NA)))
  summary(tch91.piM.amm)
  plot(tch91.piM.amm) 

plot(mir.amm$time_hrs[mir.amm$conc==0], mir.amm$alive[mir.amm$conc==0]/mir.amm$total[mir.amm$conc==0], 
     pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tch91.piM.amm, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in 2:length(unique(mir.amm$conc))){
    points(mir.amm$time_hrs[mir.amm$conc==unique(mir.amm$conc)[i]], 
           mir.amm$alive[mir.amm$conc==unique(mir.amm$conc)[i]] /
           mir.amm$total[mir.amm$conc==unique(mir.amm$conc)[i]], pch=16,
           col = i)
    lines(time, predict(tch91.piM.amm, 
                        data.frame(time_hrs=time, conc = unique(mir.amm$conc)[i])),
          lty = 2, col = i)
  }

title('ammonium phosphate toxicity to miracidia')
legend('topright', legend = c('control', unique(mir.amm$conc)[-1]), 
       pch = c(17,16,16,16,16,16), col = c(1:6), cex=0.8, bty = 'n')  

#Get estimate of miracidia-hrs for each concentration    
tch91.amm.pim.aucs = as.numeric()
  
  for(j in 1:length(unique(mir.amm$conc))){
    fx = function(t){
      predict(tch91.piM.amm, newdata = data.frame(time_hrs = t, conc = unique(mir.amm$conc)[j]))
    }
    tch91.amm.pim.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }


#Create data frame with parameter values and ammonium sulphate concentrations #######################
miramm.df = data.frame(amm = unique(mir.amm$conc),
                     logamm = log(unique(mir.amm$conc)+1),
                     e = summary(tch91.piM.amm)$coefficients[c(7:12), 1],
                     e.se = summary(tch91.piM.amm)$coefficients[c(7:12), 2],
                     b = summary(tch91.piM.amm)$coefficients[c(1:6), 1],
                     b.se = summary(tch91.piM.amm)$coefficients[c(1:6), 2])

plot(miramm.df$amm, miramm.df$e, pch = 16, xlab = 'amm. sulphate (ppm)', ylab = 'LL.2 Parameters',
     ylim = c(0, 7))
  points(miramm.df$amm, miramm.df$b, pch = 17, col=2)
    for(i in 1:length(miramm.df$amm)){
      segments(x0 = miramm.df$amm[i], y0 = miramm.df$e[i] + miramm.df$e.se[i],
               x1 = miramm.df$amm[i], y1 = miramm.df$e[i] - miramm.df$e.se[i])
      segments(x0 = miramm.df$amm[i], y0 = miramm.df$b[i] + miramm.df$b.se[i],
               x1 = miramm.df$amm[i], y1 = miramm.df$b[i] - miramm.df$b.se[i], col=2)
    }
#linear fit
tch91.e.amm.lin = lm(e ~ amm, weights = e.se^-1, data = miramm.df) 
  tch91.e.amm.lin.pred = function(amm){
    predict(tch91.e.amm.lin, newdata = data.frame(amm = amm), interval = 'confidence', level = 0.95)
  }
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.e.amm.lin.pred)[1,],
          lty = 2)
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.e.amm.lin.pred)[2,],
          lty = 3)
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.e.amm.lin.pred)[3,],
          lty = 3)
    
#exponential fit
tch91.e.amm.exp = lm(e ~ logamm, weights = e.se^-1, data = miramm.df) 
  tch91.e.amm.exp.pred = function(amm){
      predict(tch91.e.amm.exp, newdata = data.frame(logamm = log(amm+1)), interval = 'confidence', level = 0.95)
    }
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.e.amm.exp.pred)[1,],
          lty = 2, col = 3)
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.e.amm.exp.pred)[2,],
          lty = 3, col = 3)
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.e.amm.exp.pred)[3,],
          lty = 3, col = 3)    
    
AIC(tch91.e.amm.lin, tch91.e.amm.exp) #exponential is a better fit
    
tch91.b.amm = lm(b ~ amm, weights = b.se^-1, data = miramm.df) 
  tch91.b.amm.pred = function(amm){
    predict(tch91.b.amm, newdata = data.frame(amm = amm), interval = 'confidence', level = 0.95)
  }
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.b.amm.pred)[1,],
          lty = 2, col = 2)
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.b.amm.pred)[2,],
          lty = 3, col = 2)
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.b.amm.pred)[3,],
          lty = 3, col = 2)

mir.auc.amm = data.frame(amm = miramm.df$amm,
                         piM = c(tch91.amm.pim.aucs)/tch91.amm.pim.aucs[1])

#Create function to generate d-r function #####################
piM.tch91_amm_exp_auc0 = function(In){
    Ins = In/1000
    e0 = as.numeric(predict(tch91.e.amm.exp, newdata = data.frame(logamm = log(Ins+1)), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(tch91.b.amm, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
    
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value

  return(auc0)
}  

piM.tch91_amm_unc = function(In){
    Ins = In/1000
  e = as.numeric(predict(tch91.e.amm.exp, newdata = data.frame(logamm = log(Ins+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch91.b.amm, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, e[1], e[2])
  
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                  stop.on.error = FALSE)[1]$value
  piM = auc/piM.tch91_amm_exp_auc0(0) #tch91.amm.pim.aucs[1]

  return(piM)
} 

piM.tch91_amm_lin_auc0 = function(In){
  Ins = In/1000
  e0 = as.numeric(predict(tch91.e.amm.lin, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(tch91.b.amm, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                   stop.on.error = FALSE)[1]$value
  
  return(auc0)
}  

piM.tch91_amm_lin = function(In){
    Ins = In/1000
    e = as.numeric(predict(tch91.e.amm.lin, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
    b = as.numeric(predict(tch91.b.amm, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                    stop.on.error = FALSE)[1]$value
    piM = auc/piM.tch91_amm_lin_auc0(0) #tch91.amm.pim.aucs[1]
  
  return(piM)
}  
  
keep.tch91.amm = c('tch91.amm.pim.aucs', 'piM.tch91_amm_lin', 'piM.tch91_amm_exp_auc0', 'piM.tch91_amm_unc', 'piM.tch91_amm_lin_auc0',
                   'L.3.fx', 'miramm.df', 'tch91.e.amm.exp', 'tch91.e.amm.lin', 'tch91.b.amm', 'mir.auc.amm')

#Qualitative model validation ###############
#Regenerate plot of observed data
plot(mir.amm$time_hrs[mir.amm$conc==0], mir.amm$surv[mir.amm$conc==0], pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 2:length(unique(mir.amm$conc))){
    points(mir.amm$time_hrs[mir.amm$conc==unique(mir.amm$conc)[i]], 
           mir.amm$surv[mir.amm$conc==unique(mir.amm$conc)[i]], pch=16,
           col = i)
  }
  title('ammonium phosphate toxicity to miracidia')
  legend('topright', legend = c('control', unique(mir.amm$conc)[-1]), 
         pch = c(17,16,16,16,16,16), col = c(1:6), cex=0.8, bty = 'n')

#function to plot model predictions
predm.amm.plot = function(In, clr){
  Ins = In/1000
  e = as.numeric(predict(tch91.e.amm.lin, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch91.b.amm, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  lines(time, L.3.fx(time, slp = b.use, lc50 = e.use), lty=2, col = clr)
}   

#plot model predictions
for(i in c(unique(mir.amm$conc)*1000)){
  c = i/1000000 + 1
  print(c)
  replicate(10, predm.amm.plot(In = i, clr = c))
}

#plot AUC predictions
plot(mir.auc.amm[,1], mir.auc.amm[,2], pch = 16, cex = 1.2)
  points(seq(0,max(miramm.df$amm)+1,25), sapply(seq(0,max(miramm.df$amm)*1000,25000) , piM.tch91_amm_lin),
         pch = 17, cex = 0.5, col=4)
  points(seq(0,max(miramm.df$amm)+1,25), sapply(seq(0,max(miramm.df$amm)*1000,25000) , piM.tch91_amm_unc),
         pch = 17, cex = 0.5, col=2)
  
  
#Urea concentration ############
tch91.piM.ure<-drm(alive/total ~ time_hrs, conc, weights = total, data = mir.ure, type = 'binomial', 
                   fct = LL.3(names = c('b', 'd', 'e'),
                              fixed = c(NA, 1, NA)))
summary(tch91.piM.ure)
plot(tch91.piM.ure) 

plot(mir.ure$time_hrs[mir.ure$conc==0], mir.ure$alive[mir.ure$conc==0]/mir.ure$total[mir.ure$conc==0], 
     pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tch91.piM.ure, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
    for(i in 2:length(unique(mir.ure$conc))){
      points(mir.ure$time_hrs[mir.ure$conc==unique(mir.ure$conc)[i]], 
             mir.ure$alive[mir.ure$conc==unique(mir.ure$conc)[i]] /
               mir.ure$total[mir.ure$conc==unique(mir.ure$conc)[i]], pch=16,
             col = i)
      lines(time, predict(tch91.piM.ure, 
                          data.frame(time_hrs=time, conc = unique(mir.ure$conc)[i])),
            lty = 2, col = i)
    }

title('urea toxicity to miracidia')
legend('topright', legend = c('control', unique(mir.ure$conc)[-1]), 
       pch = c(17,16,16,16,16,16), col = c(1:6), cex=0.8, bty = 'n')  

#Get estimate of miracidia-hrs for each concentration    
tch91.ure.pim.aucs = as.numeric()

  for(j in 1:length(unique(mir.ure$conc))){
    fx = function(t){
      predict(tch91.piM.ure, newdata = data.frame(time_hrs = t, conc = unique(mir.ure$conc)[j]))
    }
    tch91.ure.pim.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }


#Create data frame with parameter values and urea concentrations #######################
mirure.df = data.frame(ure = unique(mir.ure$conc),
                       logure = log(unique(mir.ure$conc)+1),
                       e = summary(tch91.piM.ure)$coefficients[c(7:12), 1],
                       e.se = summary(tch91.piM.ure)$coefficients[c(7:12), 2],
                       b = summary(tch91.piM.ure)$coefficients[c(1:6), 1],
                       b.se = summary(tch91.piM.ure)$coefficients[c(1:6), 2])

plot(mirure.df$ure, mirure.df$e, pch = 16, xlab = 'urea (ppb)', ylab = 'LL.2 Parameters',
     ylim = c(0, 7))
  points(mirure.df$ure, mirure.df$b, pch = 17, col=2)
  for(i in 1:length(mirure.df$ure)){
    segments(x0 = mirure.df$ure[i], y0 = mirure.df$e[i] + mirure.df$e.se[i],
             x1 = mirure.df$ure[i], y1 = mirure.df$e[i] - mirure.df$e.se[i])
    segments(x0 = mirure.df$ure[i], y0 = mirure.df$b[i] + mirure.df$b.se[i],
             x1 = mirure.df$ure[i], y1 = mirure.df$b[i] - mirure.df$b.se[i], col=2)
  }
#linear fit
tch91.e.ure.lin = lm(e ~ ure, weights = e.se^-1, data = mirure.df) 
  tch91.e.ure.lin.pred = function(ure){
    predict(tch91.e.ure.lin, newdata = data.frame(ure = ure), interval = 'confidence', level = 0.95)
  }
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.e.ure.lin.pred)[1,],
          lty = 2)
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.e.ure.lin.pred)[2,],
          lty = 3)
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.e.ure.lin.pred)[3,],
          lty = 3)
  
#exponential fit
tch91.e.ure.exp = lm(e ~ logure, weights = e.se^-1, data = mirure.df) 
  tch91.e.ure.exp.pred = function(ure){
    predict(tch91.e.ure.exp, newdata = data.frame(logure = log(ure+1)), interval = 'confidence', level = 0.95)
  }
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.e.ure.exp.pred)[1,],
          lty = 2, col = 3)
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.e.ure.exp.pred)[2,],
          lty = 3, col = 3)
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.e.ure.exp.pred)[3,],
          lty = 3, col = 3)    

AIC(tch91.e.ure.lin, tch91.e.ure.exp) #Linear is a better fit

  tch91.b.ure = lm(b ~ ure, weights = b.se^-1, data = mirure.df) 
    tch91.b.ure.pred = function(ure){
      predict(tch91.b.ure, newdata = data.frame(ure = ure), interval = 'confidence', level = 0.95)
    }
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.b.ure.pred)[1,],
          lty = 2, col = 2)
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.b.ure.pred)[2,],
          lty = 3, col = 2)
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.b.ure.pred)[3,],
          lty = 3, col = 2)

mir.auc.ure = data.frame(ure = mirure.df$ure,
                         piM = c(tch91.ure.pim.aucs)/tch91.ure.pim.aucs[1])

#Create function to generate d-r function #####################

piM.tch91_ure_lin_auc0 = function(In){
  Ins = In/1000
  e0 = as.numeric(predict(tch91.e.ure.lin, newdata = data.frame(ure = Ins), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(tch91.b.ure, newdata = data.frame(ure = Ins), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                   stop.on.error = FALSE)[1]$value
  
  return(auc0)
}  

piM.tch91_ure_unc = function(In){
  Ins = In/1000
  e = as.numeric(predict(tch91.e.ure.lin, newdata = data.frame(ure = Ins), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch91.b.ure, newdata = data.frame(ure = Ins), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                  stop.on.error = FALSE)[1]$value
  piM = auc/piM.tch91_ure_lin_auc0(0) #tch91.amm.pim.aucs[1]
  
  return(piM)
}  

keep.tch91.ure = c('tch91.ure.pim.aucs', 'tch91.e.ure.exp', 'piM.tch91_ure_unc', 'piM.tch91_ure_lin_auc0', 'L.3.fx', 'mirure.df',
                   'tch91.e.ure.lin', 'tch91.b.ure', 'mir.auc.ure')

#Qualitative model validation ###############
#Regenerate plot of observed data
plot(mir.ure$time_hrs[mir.ure$conc==0], mir.ure$surv[mir.ure$conc==0], pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
for(i in 2:length(unique(mir.ure$conc))){
  points(mir.ure$time_hrs[mir.ure$conc==unique(mir.ure$conc)[i]], 
         mir.ure$surv[mir.ure$conc==unique(mir.ure$conc)[i]], pch=16,
         col = i)
}
title('urea toxicity to miracidia')
legend('topright', legend = c('control', unique(mir.ure$conc)[-1]), 
       pch = c(17,16,16,16,16,16), col = c(1:6), cex=0.8, bty = 'n')

#function to plot model predictions
predm.ure.plot = function(In, clr){
  Ins = In/1000
  e = as.numeric(predict(tch91.e.ure.lin, newdata = data.frame(ure = Ins), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch91.b.ure, newdata = data.frame(ure = Ins), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  lines(time, L.3.fx(time, slp = b.use, lc50 = e.use), lty=2, col = clr)
}   

#plot model predictions
for(i in c(unique(mir.ure$conc)*1000)){
  c = i/2000000 + 1
  print(c)
  replicate(10, predm.ure.plot(In = i, clr = c))
}

#plot AUC predictions
plot(mir.auc.ure[,1], mir.auc.ure[,2], pch = 16, cex = 1.2)
  points(seq(0,max(mirure.df$ure)+1,25), sapply(seq(0,max(mirure.df$ure)*1000,25000) , piM.tch91_ure_unc),
         pch = 17, cex = 0.5, col=4)


#Final keep vector #########
keep.tch91.Fe = c(keep.tch91.amm, keep.tch91.ure)