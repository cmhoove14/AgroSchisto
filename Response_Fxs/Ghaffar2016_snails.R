#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Ghaffar 2016 SNAIL (B. alexandrina) data
require(drc)

#Snail (30 B. alexandrina 6-8mm shell width) toxicity ##########
  #Reported lc50/slope combinations did not fit provided data well therefore functions fit to the 
    #reported lc0, lc10, lc25, lc50, and lc90 values were used to model excess daily per capita death rate
#Butralin ##############
but.dat = data.frame(lcs = c(0, 0, 10, 25, 50, 90),
                     butralin = c(0, 556, 2417, 3906, 5560, 8703))
  

#Fit function based on provided lc values
  gaf.muNq.but<-drm(lcs/100 ~ butralin, data = but.dat, weights = rep(30, 6),
                    type = 'binomial', fct = LL.2())
    summary(gaf.muNq.but)
    
  plot(but.dat$butralin, but.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,13000),
       xlab = 'Butralin (ppb)', ylab = 'snail mortality rate')
    segments(x0 = 3700, x1 = 8340, y0 = 0.5, y1 = 0.5, lty = 1, col = 1)
    segments(x0 = 6220, x1 = 12800, y0 = 0.9, y1 = 0.9, lty = 1, col = 1)
  
    mu_Nq_butr_gaf16<-function(He){
      predict(gaf.muNq.but, data.frame(butralin = He), 
              interval = 'confidence', level = 0.95)
    }  
  
  lines(seq(0,13000,100), sapply(seq(0,13000,100), mu_Nq_butr_gaf16, simplify = T)[1,],
        lty = 2, col = 2) 
  lines(seq(0,13000,100), sapply(seq(0,13000,100), mu_Nq_butr_gaf16, simplify = T)[2,],
        lty = 3, col = 2) 
  lines(seq(0,13000,100), sapply(seq(0,13000,100), mu_Nq_butr_gaf16, simplify = T)[3,],
        lty = 3, col = 2) 
    
    mu_Nq_butr_gaf16_uncertainty<-function(He){
      rdrm(1, LL.2(), coef(gaf.muNq.but), He, yerror = 'rbinom', ypar = 30)$y / 30  
    }

    points(seq(0,13000,100), sapply(seq(0,13000,100), 
                                   mu_Nq_butr_gaf16_uncertainty, 
                                   simplify = T),
           pch = 5, col = 4, cex = 0.5)

#Glyphosate #############
gly.dat = data.frame(lcs = c(0, 0, 10, 25, 50, 90),
                     glyphosate = c(0, 1506, 3875, 9174, 15062, 26249))
  
#Fit function based on provided lc values
  gaf.muNq.gly<-drm(lcs/100 ~ glyphosate, data = gly.dat, weights = rep(30, 6),
                    type = 'binomial', fct = LL.2())
    summary(gaf.muNq.gly)
    
  plot(gly.dat$glyphosate, gly.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,33000),
       xlab = 'Glyphosate (ppb)', ylab = 'snail mortality rate')
    segments(x0 = 9130, x1 = 16570, y0 = 0.5, y1 = 0.5, lty = 1, col = 1)
    segments(x0 = 23870, x1 = 28900, y0 = 0.9, y1 = 0.9, lty = 1, col = 1)
    
  
  mu_Nq_gly_gaf16<-function(He){
    predict(gaf.muNq.gly, data.frame(butralin = He), 
            interval = 'confidence', level = 0.95)
  }  
  
  lines(seq(0,33000,150), sapply(seq(0,33000,150), mu_Nq_gly_gaf16, simplify = T)[1,],
        lty = 2, col = 2) 
  lines(seq(0,33000,150), sapply(seq(0,33000,150), mu_Nq_gly_gaf16, simplify = T)[2,],
        lty = 3, col = 2) 
  lines(seq(0,33000,150), sapply(seq(0,33000,150), mu_Nq_gly_gaf16, simplify = T)[3,],
        lty = 3, col = 2) 
  
    mu_Nq_gly_gaf16_uncertainty<-function(He){
      rdrm(1, LL.2(), coef(gaf.muNq.gly), He, yerror = 'rbinom', ypar = 30)$y / 30  
    }
    
    points(seq(0,33000,200), sapply(seq(0,33000,200), 
                                    mu_Nq_gly_gaf16_uncertainty, 
                                    simplify = T),
           pch = 5, col = 4, cex = 0.5)
    
#Pendimethalin ##############
pen.dat = data.frame(lcs = c(0, 0, 10, 25, 50, 90),
                     pendimethalin = c(0, 214.8, 535, 1299, 2148, 3762))

  
#Fit function based on provided lc values
  gaf.muNq.pen<-drm(lcs/100 ~ pendimethalin, data = pen.dat, weights = rep(30, 6),
                    type = 'binomial', fct = LL.2())
    summary(gaf.muNq.pen)
    
  plot(pen.dat$pendimethalin, pen.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,7000),
       xlab = 'pendimethalin (ppb)', ylab = 'snail mortality rate')
    segments(x0 = 1430, x1 = 3220, y0 = 0.5, y1 = 0.5, lty = 1, col = 1)
    segments(x0 = 2400, x1 = 6420, y0 = 0.9, y1 = 0.9, lty = 1, col = 1)
    
  mu_Nq_pen_gaf16<-function(He){
    predict(gaf.muNq.pen, data.frame(glyphosate = He), 
            interval = 'confidence', level = 0.95)  
  }
  
  lines(seq(0,7000,50), sapply(seq(0,7000,50), mu_Nq_pen_gaf16, simplify = T)[1,],
        lty = 2, col = 2) 
  lines(seq(0,7000,50), sapply(seq(0,7000,50), mu_Nq_pen_gaf16, simplify = T)[2,],
        lty = 3, col = 2) 
  lines(seq(0,7000,50), sapply(seq(0,7000,50), mu_Nq_pen_gaf16, simplify = T)[3,],
        lty = 3, col = 2) 
  
  mu_Nq_pen_gaf16_uncertainty<-function(He){
    rdrm(1, LL.2(), coef(gaf.muNq.pen), He, yerror = 'rbinom', ypar = 30)$y / 30  
  }
  
  points(seq(0,7000,50), sapply(seq(0,7000,50), 
                                  mu_Nq_pen_gaf16_uncertainty, 
                                  simplify = T),
         pch = 5, col = 4, cex = 0.5)

#Longitudinal juvenile snail (B. alexandrina <3.6mm shell width) toxicity ###########
gaf16.juv<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/ghaffar2016_juv.csv')
  cont.surv = gaf16.juv$prop_surv[gaf16.juv$chem == 'control']
  days = gaf16.juv$time[gaf16.juv$chem == 'control']
  
  plot(days, cont.surv, pch = 16, xlab = 'time (days)', ylab = 'Prop survival', ylim = c(0,1))

#Function to generate mortality estimates and compare to observed    
  pred = function(N0 = 100, t, rate){
    Nt = N0*exp(t*rate)
    return(Nt)
  }

#Get background mortality rate    
  bakground<-nls(cont.surv ~ 1 * exp(-r * days), start = list(r = 0.01))
    summary(bakground)
    
  lines(c(0:60), pred(t = c(0:60), rate = -summary(bakground)$parameters[1])/100, lty=2)

#Fill for butralin data  
  for(i in 1:length(unique(gaf16.juv$conc[gaf16.juv$chem == 'butralin']))){
    points(gaf16.juv$time[gaf16.juv$chem == 'butralin' & gaf16.juv$conc == unique(gaf16.juv$conc[gaf16.juv$chem == 'butralin'])[i]],
           gaf16.juv$prop_surv[gaf16.juv$chem == 'butralin' & gaf16.juv$conc == unique(gaf16.juv$conc[gaf16.juv$chem == 'butralin'])[i]],
           pch = 17, col = i+1)
  }
  
  gaf16.juv$pred = 0
  gaf16.juv$pred[gaf16.juv$chem == 'butralin'] = 
    pred(t = gaf16.juv$time[gaf16.juv$chem == 'butralin'],
         rate = -mu_Nq_butr_gaf16(gaf16.juv$conc[gaf16.juv$chem == 'butralin'])[,1] - 
                summary(bakground)$parameters[1])/100
  
  for(i in 1:length(unique(gaf16.juv$conc[gaf16.juv$chem == 'butralin']))){
    lines(gaf16.juv$time[gaf16.juv$chem == 'butralin' & 
                           gaf16.juv$conc == unique(gaf16.juv$conc[gaf16.juv$chem == 'butralin'])[i]],
          gaf16.juv$pred[gaf16.juv$chem == 'butralin' & 
                           gaf16.juv$conc == unique(gaf16.juv$conc[gaf16.juv$chem == 'butralin'])[i]],
           lty=2, col = i+1)
  }
  
  legend('right', title = 'Butralin', legend = c('control', '556 ppb', '2417 ppb', '3906 ppb'), 
         pch = c(16,17,17,17), col = c(1,2,3,4), cex = 0.7)
  legend('topright',  legend = 'predicted', lty = 2, cex = 0.7)
  title(main = 'Ghaffar 2016 long. survival to Butralin',
        sub = 'predicted vs observed (juveniles)')
  
#Fill for glyphosate data  
  plot(days, cont.surv, pch = 16, xlab = 'time (days)', ylab = 'Prop survival', ylim = c(0,1))
  lines(c(0:60), pred(t = c(0:60), rate = -summary(bakground)$parameters[1])/100, lty=2)
  
  for(i in 1:length(unique(gaf16.juv$conc[gaf16.juv$chem == 'glyphosate']))){
    points(gaf16.juv$time[gaf16.juv$chem == 'glyphosate' & gaf16.juv$conc == unique(gaf16.juv$conc[gaf16.juv$chem == 'glyphosate'])[i]],
           gaf16.juv$prop_surv[gaf16.juv$chem == 'glyphosate' & gaf16.juv$conc == unique(gaf16.juv$conc[gaf16.juv$chem == 'glyphosate'])[i]],
           pch = 17, col = i+1)
  }
  
  gaf16.juv$pred[gaf16.juv$chem == 'glyphosate'] = 
    pred(t = gaf16.juv$time[gaf16.juv$chem == 'glyphosate'],
         rate = -mu_Nq_gly_gaf16(gaf16.juv$conc[gaf16.juv$chem == 'glyphosate'])[,1] - 
                summary(bakground)$parameters[1])/100
  
  for(i in 1:length(unique(gaf16.juv$conc[gaf16.juv$chem == 'glyphosate']))){
    lines(gaf16.juv$time[gaf16.juv$chem == 'glyphosate' & 
                           gaf16.juv$conc == unique(gaf16.juv$conc[gaf16.juv$chem == 'glyphosate'])[i]],
          gaf16.juv$pred[gaf16.juv$chem == 'glyphosate' & 
                           gaf16.juv$conc == unique(gaf16.juv$conc[gaf16.juv$chem == 'glyphosate'])[i]],
          lty=2, col = i+1)
  }
  
  legend('right', title = 'glyphosate', legend = c('control', '1506 ppb', '3875 ppb', '9174 ppb'), 
         pch = c(16,17,17,17), col = c(1,2,3,4), cex = 0.7)
  legend('topright',  legend = 'predicted', lty = 2, cex = 0.7)
  title(main = 'Ghaffar 2016 long. survival to Glyphosate',
        sub = 'predicted vs observed (juveniles)') 
  
#Fill for pendimethalin data  
  plot(days, cont.surv, pch = 16, xlab = 'time (days)', ylab = 'Prop survival', ylim = c(0,1))
  lines(c(0:60), pred(t = c(0:60), rate = -summary(bakground)$parameters[1])/100, lty=2)
  
  for(i in 1:length(unique(gaf16.juv$conc[gaf16.juv$chem == 'pendimethalin']))){
    points(gaf16.juv$time[gaf16.juv$chem == 'pendimethalin' & gaf16.juv$conc == unique(gaf16.juv$conc[gaf16.juv$chem == 'pendimethalin'])[i]],
           gaf16.juv$prop_surv[gaf16.juv$chem == 'pendimethalin' & gaf16.juv$conc == unique(gaf16.juv$conc[gaf16.juv$chem == 'pendimethalin'])[i]],
           pch = 17, col = i+1)
  }
  
  gaf16.juv$pred[gaf16.juv$chem == 'pendimethalin'] = 
    pred(t = gaf16.juv$time[gaf16.juv$chem == 'pendimethalin'],
         rate = -mu_Nq_pen_gaf16(gaf16.juv$conc[gaf16.juv$chem == 'pendimethalin'])[,1] - 
                summary(bakground)$parameters[1])/100
  
  for(i in 1:length(unique(gaf16.juv$conc[gaf16.juv$chem == 'pendimethalin']))){
    lines(gaf16.juv$time[gaf16.juv$chem == 'pendimethalin' & 
                           gaf16.juv$conc == unique(gaf16.juv$conc[gaf16.juv$chem == 'pendimethalin'])[i]],
          gaf16.juv$pred[gaf16.juv$chem == 'pendimethalin' & 
                           gaf16.juv$conc == unique(gaf16.juv$conc[gaf16.juv$chem == 'pendimethalin'])[i]],
          lty=2, col = i+1)
  }
  
  legend('right', title = 'pendimethalin', legend = c('control', '214 ppb', '535 ppb', '1299 ppb'), 
         pch = c(16,17,17,17), col = c(1,2,3,4), cex = 0.7)
  legend('topright',  legend = 'predicted', lty = 2, cex = 0.7)
  title(main = 'Ghaffar 2016 long. survival to Pendimethalin',
        sub = 'predicted vs observed (juveniles)')  
  
#longitudinal adult survival (B. alexandrina 6-8mm shell width) toxicity #############
  gaf16.ad<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/ghaffar2016_adult.csv')
    cont.surv2 = gaf16.ad$prop_surv[gaf16.ad$chem == 'control']
    days2 = gaf16.ad$time[gaf16.ad$chem == 'control']
  
  plot(days2, cont.surv2, pch = 16, xlab = 'time (days)', ylab = 'Prop survival', ylim = c(0,1))
  
#Get background mortality rate    
  bakground2<-nls(cont.surv2 ~ 1 * exp(-r * days2), start = list(r = 0.01))
    summary(bakground2)
  
  lines(c(0:60), pred(t = c(0:60), rate = -summary(bakground2)$parameters[1])/100, lty=2)
  
#Fill for butralin data  
  for(i in 1:length(unique(gaf16.ad$conc[gaf16.ad$chem == 'butralin']))){
    points(gaf16.ad$time[gaf16.ad$chem == 'butralin' & gaf16.ad$conc == unique(gaf16.ad$conc[gaf16.ad$chem == 'butralin'])[i]],
           gaf16.ad$prop_surv[gaf16.ad$chem == 'butralin' & gaf16.ad$conc == unique(gaf16.ad$conc[gaf16.ad$chem == 'butralin'])[i]],
           pch = 17, col = i+1)
  }
  
  gaf16.ad$pred = 0
  gaf16.ad$pred[gaf16.ad$chem == 'butralin'] = 
    pred(t = gaf16.ad$time[gaf16.ad$chem == 'butralin'],
         rate = -mu_Nq_butr_gaf16(gaf16.ad$conc[gaf16.ad$chem == 'butralin'])[,1] - 
           summary(bakground2)$parameters[1])/100
  
  for(i in 1:length(unique(gaf16.ad$conc[gaf16.ad$chem == 'butralin']))){
    lines(gaf16.ad$time[gaf16.ad$chem == 'butralin' & 
                           gaf16.ad$conc == unique(gaf16.ad$conc[gaf16.ad$chem == 'butralin'])[i]],
          gaf16.ad$pred[gaf16.ad$chem == 'butralin' & 
                           gaf16.ad$conc == unique(gaf16.ad$conc[gaf16.ad$chem == 'butralin'])[i]],
          lty=2, col = i+1)
  }
  
  legend('right', title = 'Butralin', legend = c('control', '556 ppb', '2417 ppb', '3906 ppb'), 
         pch = c(16,17,17,17), col = c(1,2,3,4), cex = 0.7)
  legend('topright',  legend = 'predicted', lty = 2, cex = 0.7)
  title(main = 'Ghaffar 2016 long. survival to Butralin',
        sub = 'predicted vs observed (adults)')
  
#Fill for glyphosate data  
  plot(days2, cont.surv2, pch = 16, xlab = 'time (days)', ylab = 'Prop survival', ylim = c(0,1))
  lines(c(0:60), pred(t = c(0:60), rate = -summary(bakground2)$parameters[1])/100, lty=2)
  
  for(i in 1:length(unique(gaf16.ad$conc[gaf16.ad$chem == 'glyphosate']))){
    points(gaf16.ad$time[gaf16.ad$chem == 'glyphosate' & gaf16.ad$conc == unique(gaf16.ad$conc[gaf16.ad$chem == 'glyphosate'])[i]],
           gaf16.ad$prop_surv[gaf16.ad$chem == 'glyphosate' & gaf16.ad$conc == unique(gaf16.ad$conc[gaf16.ad$chem == 'glyphosate'])[i]],
           pch = 17, col = i+1)
  }
  
  gaf16.ad$pred[gaf16.ad$chem == 'glyphosate'] = 
    pred(t = gaf16.ad$time[gaf16.ad$chem == 'glyphosate'],
         rate = -mu_Nq_gly_gaf16(gaf16.ad$conc[gaf16.ad$chem == 'glyphosate'])[,1] - 
           summary(bakground2)$parameters[1])/100
  
  for(i in 1:length(unique(gaf16.ad$conc[gaf16.ad$chem == 'glyphosate']))){
    lines(gaf16.ad$time[gaf16.ad$chem == 'glyphosate' & 
                           gaf16.ad$conc == unique(gaf16.ad$conc[gaf16.ad$chem == 'glyphosate'])[i]],
          gaf16.ad$pred[gaf16.ad$chem == 'glyphosate' & 
                           gaf16.ad$conc == unique(gaf16.ad$conc[gaf16.ad$chem == 'glyphosate'])[i]],
          lty=2, col = i+1)
  }
  
  legend('right', title = 'glyphosate', legend = c('control', '1506 ppb', '3875 ppb', '9174 ppb'), 
         pch = c(16,17,17,17), col = c(1,2,3,4), cex = 0.7)
  legend('topright',  legend = 'predicted', lty = 2, cex = 0.7)
  title(main = 'Ghaffar 2016 long. survival to Glyphosate',
        sub = 'predicted vs observed (adults)') 
  
#Fill for pendimethalin data  
  plot(days2, cont.surv2, pch = 16, xlab = 'time (days)', ylab = 'Prop survival', ylim = c(0,1))
  lines(c(0:60), pred(t = c(0:60), rate = -summary(bakground2)$parameters[1])/100, lty=2)
  
  for(i in 1:length(unique(gaf16.ad$conc[gaf16.ad$chem == 'pendimethalin']))){
    points(gaf16.ad$time[gaf16.ad$chem == 'pendimethalin' & gaf16.ad$conc == unique(gaf16.ad$conc[gaf16.ad$chem == 'pendimethalin'])[i]],
           gaf16.ad$prop_surv[gaf16.ad$chem == 'pendimethalin' & gaf16.ad$conc == unique(gaf16.ad$conc[gaf16.ad$chem == 'pendimethalin'])[i]],
           pch = 17, col = i+1)
  }
  
  gaf16.ad$pred[gaf16.ad$chem == 'pendimethalin'] = 
    pred(t = gaf16.ad$time[gaf16.ad$chem == 'pendimethalin'],
         rate = -mu_Nq_pen_gaf16(gaf16.ad$conc[gaf16.ad$chem == 'pendimethalin'])[,1] - 
           summary(bakground2)$parameters[1])/100
  
  for(i in 1:length(unique(gaf16.ad$conc[gaf16.ad$chem == 'pendimethalin']))){
    lines(gaf16.ad$time[gaf16.ad$chem == 'pendimethalin' & 
                           gaf16.ad$conc == unique(gaf16.ad$conc[gaf16.ad$chem == 'pendimethalin'])[i]],
          gaf16.ad$pred[gaf16.ad$chem == 'pendimethalin' & 
                           gaf16.ad$conc == unique(gaf16.ad$conc[gaf16.ad$chem == 'pendimethalin'])[i]],
          lty=2, col = i+1)
  }
  
  legend('right', title = 'pendimethalin', legend = c('control', '214 ppb', '535 ppb', '1299 ppb'), 
         pch = c(16,17,17,17), col = c(1,2,3,4), cex = 0.7)
  legend('topright',  legend = 'predicted', lty = 2, cex = 0.7)
  title(main = 'Ghaffar 2016 long. survival to Pendimethalin',
        sub = 'predicted vs observed (adults)')
  
#Reductions in adult snail (B. alexandrina) reproduction ###########  
  #Paper reports reproduction as R0: the summed product of live snails * eggs produced
  #we just want reduction in eggs produced as the model already takes additional mortality 
  #into account; therefore we normalize reported R0s by relative survival between treatment groups
  #to get relative estimate of eggs/snail/week
sn.wk = as.numeric()  #snails/week estimates for each treatment
for(i in 1:length(unique(gaf16.ad$conc))){
  sn.wk[i] = sum(gaf16.ad$prop_surv[gaf16.ad$conc == unique(gaf16.ad$conc)[i]])
}  

htch.wk = c(0.99, 0.92, 0.82, 0.24,
            0.89, 0.58, 0.46,
            0.9, 0.3, 0.05)   #vector of hatching proportions (hatches/egg)

#vector of concnetrations and reproduction measured as approximate hatchlings/snail/week
  gafrep = data.frame('but.conc'=but.dat$butralin[1:4],
                      'gly.conc'=gly.dat$glyphosate[1:4],
                      'pen.conc'=pen.dat$pendimethalin[1:4],
                      'but.rep'=c(44.231/sn.wk[1]*htch.wk[1], 7.05/sn.wk[2]*htch.wk[2], 
                                  4.52/sn.wk[3]*htch.wk[3], 4.04/sn.wk[4]*htch.wk[4]),
                      'gly.rep'=c(44.231/sn.wk[1]*htch.wk[1], 4.87/sn.wk[5]*htch.wk[5], 
                                  4.20/sn.wk[6]*htch.wk[6], 4.23/sn.wk[7]*htch.wk[7]),
                      'pen.rep'=c(44.231/sn.wk[1]*htch.wk[1], 5.18/sn.wk[8]*htch.wk[8], 
                                  4.75/sn.wk[9]*htch.wk[9], 4.27/sn.wk[10]*htch.wk[10]))
  
#Reductions from butralin ##########
plot(gafrep$but.conc, gafrep$but.rep/gafrep$but.rep[1], pch = 16, ylim = c(0,1),
     xlab = 'Butralin (ppb)', ylab = 'relative reproduction rate')
    
  but.r0 = drm(but.rep ~ but.conc, data = gafrep, type = 'continuous',
               fct = LL.4(names = c("b", "c", "d", "e"),
                          fixed = c(NA, gafrep$but.rep[4], gafrep$but.rep[1], NA))) 
    summary(but.r0)
    
  but.r0.pred = function(He){
    predict(but.r0, newdata = data.frame(but.conc = He), interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,4500,10), sapply(seq(0,4500,10), but.r0.pred, simplify = T)[1,] / gafrep$but.rep[1],
          lty = 2, col = 2) 
    lines(seq(0,4500,10), sapply(seq(0,4500,10), but.r0.pred, simplify = T)[2,] / gafrep$but.rep[1],
          lty = 3, col = 2) 
    lines(seq(0,4500,10), sapply(seq(0,4500,10), but.r0.pred, simplify = T)[3,] / gafrep$but.rep[1],
          lty = 3, col = 2) 
  
 fN.butr.fx.uncertainty = function(He){
   init = predict(but.r0, newdata = data.frame(but.conc = He), se.fit = T)
   fn = rnorm(1, init[1], init[2]) / gafrep$but.rep[1]
   while(fn < 0 && fn > 1.00000){
     fn = rnorm(1, init[1], init[2]) / gafrep$but.rep[1]
   }
   return(fn)
 } #normalized to 1
 
 points(seq(0,4500,10), 
        sapply(seq(0,4500,10), fN.butr.fx.uncertainty, simplify = T),
        pch = 5, col = 4, cex = 0.5) 

#Reductions from glyphosate ##########
 plot(gafrep$gly.conc, gafrep$gly.rep/gafrep$but.rep[1], pch = 16, ylim = c(0,1),
      xlab = 'glyphosate (ppb)', ylab = 'relative reproduction rate')
 
   gly.r0 = drm(gly.rep ~ gly.conc, data = gafrep, type = 'continuous',
                fct = LL.4(names = c("b", "c", "d", "e"),
                          fixed = c(NA, gafrep$gly.rep[4], gafrep$but.rep[1], NA)))  
   summary(gly.r0)
   
 gly.r0.pred = function(He){
   predict(gly.r0, newdata = data.frame(gly.conc = He), interval = 'confidence', level = 0.95)
 }
 
   lines(c(0:1000, seq(1001, 10000,100)), sapply(c(0:1000, seq(1001, 10000,100)), 
                                                 gly.r0.pred, simplify = T)[1,]/gafrep$but.rep[1],
         lty = 2, col = 2) 
   lines(c(0:1000, seq(1001, 10000,100)), sapply(c(0:1000, seq(1001, 10000,100)), 
                                                 gly.r0.pred, simplify = T)[2,]/gafrep$but.rep[1],
         lty = 3, col = 2) 
   lines(c(0:1000, seq(1001, 10000,100)), sapply(c(0:1000, seq(1001, 10000,100)), 
                                                 gly.r0.pred, simplify = T)[3,]/gafrep$but.rep[1],
         lty = 3, col = 2) 
 
   fN.gly.fx.uncertainty = function(He){
     init = predict(gly.r0, newdata = data.frame(gly.conc = He), se.fit = T)
     fn = rnorm(1, init[1], init[2]) / gafrep$but.rep[1]
     while(fn < 0 && fn > 1){
       fn = rnorm(1, init[1], init[2]) / gafrep$but.rep[1]
     }
     return(fn)
   } #normalized to 1
 
     points(seq(1, 10000,100), 
            sapply(seq(1, 10000,100), fN.gly.fx.uncertainty, simplify = T),
            pch = 5, col = 4, cex = 0.5) 
     
#Reductions from pendimethalin ##########
plot(gafrep$pen.conc, gafrep$pen.rep/gafrep$but.rep[1], pch = 16, ylim = c(0,1),
     xlab = 'pendimethalin (ppb)', ylab = 'relative reproduction rate')
     
  pen.r0 = drm(pen.rep ~ pen.conc, data = gafrep, type = 'continuous',
               fct = LL.4(names = c("b", "c", "d", "e"),
                         fixed = c(NA, gafrep$pen.rep[4], gafrep$but.rep[1], NA)))  
    summary(pen.r0)
     
  pen.r0.pred = function(He){
      predict(pen.r0, newdata = data.frame(pen.conc = He), interval = 'confidence', level = 0.95)
  }
     
     lines(c(0:1300), sapply(c(0:1300), pen.r0.pred, simplify = T)[1,]/gafrep$but.rep[1],
           lty = 2, col = 2) 
     lines(c(0:1300), sapply(c(0:1300), pen.r0.pred, simplify = T)[2,]/gafrep$but.rep[1],
           lty = 3, col = 2) 
     lines(c(0:1300), sapply(c(0:1300), pen.r0.pred, simplify = T)[3,]/gafrep$but.rep[1],
           lty = 3, col = 2) 
     
  fN.pen.fx.uncertainty = function(He){
    init = predict(pen.r0, newdata = data.frame(pen.conc = He), se.fit = T)
    fn = rnorm(1, init[1], init[2]) / gafrep$but.rep[1]
    while(fn < 0 && fn > 1){
      fn = rnorm(1, init[1], init[2]) / gafrep$but.rep[1]
    }
    return(fn)  
  } #normalized to 1
     
     points(seq(0,1300, 10), 
            sapply(seq(0,1300, 10), fN.pen.fx.uncertainty, simplify = T),
            pch = 5, col = 4, cex = 0.5) #plot with reference back to raw control value