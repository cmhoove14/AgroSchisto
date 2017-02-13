#Data extraction and model fitting to Ghaffar 2016 SNAIL (B. alexandrina) data
require(drc)

#Snail (B. alexandrina 6-8mm shell width) toxicity ##########
  #Reported lc50/slope combinations did not fit provided data well therefore functions fit to the 
    #reported lc0, lc10, lc25, lc50, and lc90 values were used to model excess daily per capita death rate
#Butralin
  lcs = c(0, 0, 10, 25, 50, 90)
  lc.but.vals = c(0, 556, 2417, 3906, 5560, 8703)

#Fit function based on provided lc values
  gaf.muNq.but<-drm(lcs/100 ~ lc.but.vals, type = 'binomial', fct = LL.2())
    summary(gaf.muNq.but)
    
  mu_Nq_butr_gaf16<-function(In){
    1/(1+exp(gaf.muNq.but$coefficients[1]*(log(In)-log(gaf.muNq.but$coefficients[2]))))
  }

plot(c(0:30000), mu_Nq_butr_gaf16(c(0:30000)), lwd = 2, lty=2, type = 'l', col = 'darkorange', 
     xlab = 'herbicide concentration (ppb)', ylab = 'mu_Nq', ylim = c(0,1), main = 'herbicide toxicity to snails, Ghaffar2016')
  points(lc.but.vals, lcs/100, pch = 17, col = 'darkorange')

#Glyphosate
  lc.gly.vals = c(0, 1506, 3875, 9174, 15062, 26249)
  
#Fit function based on provided lc values
  gaf.muNq.gly<-drm(lcs/100 ~ lc.gly.vals, type = 'binomial', fct = LL.2())
    summary(gaf.muNq.gly)
  
  mu_Nq_gly_gaf16<-function(In){
    1/(1+exp(gaf.muNq.gly$coefficients[1]*(log(In)-log(gaf.muNq.gly$coefficients[2]))))
  }
  
  points(lc.gly.vals, lcs/100, pch = 17, col = 'firebrick')
  lines(c(0:30000), mu_Nq_gly_gaf16(c(0:30000)), lwd = 2, lty=2, col = 'firebrick')

#Pendimethalin
  lc.pen.vals = c(0, 214.8, 535, 1299, 2148, 3762)
  
#Fit function based on provided lc values
  gaf.muNq.pen<-drm(lcs/100 ~ lc.pen.vals, type = 'binomial', fct = LL.2())
    summary(gaf.muNq.pen)
  
  mu_Nq_pen_gaf16<-function(In){
    1/(1+exp(gaf.muNq.pen$coefficients[1]*(log(In)-log(gaf.muNq.pen$coefficients[2]))))
  }
  
  points(lc.pen.vals, lcs/100, pch = 17, col = 'indianred2')
  lines(c(0:30000), mu_Nq_pen_gaf16(c(0:30000)), lwd = 2, lty=2, col = 'indianred2')
  legend('bottomright', legend = c('Butralin', 'Glyphosate', 'Pendimethalin'), lwd=2, 
         col = c('darkorange','firebrick','indianred2'), cex=0.75)
  
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
         rate = -mu_Nq_butr_gaf16(gaf16.juv$conc[gaf16.juv$chem == 'butralin']) - 
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
         rate = -mu_Nq_gly_gaf16(gaf16.juv$conc[gaf16.juv$chem == 'glyphosate']) - 
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
         rate = -mu_Nq_pen_gaf16(gaf16.juv$conc[gaf16.juv$chem == 'pendimethalin']) - 
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
         rate = -mu_Nq_butr_gaf16(gaf16.ad$conc[gaf16.ad$chem == 'butralin']) - 
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
         rate = -mu_Nq_gly_gaf16(gaf16.ad$conc[gaf16.ad$chem == 'glyphosate']) - 
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
         rate = -mu_Nq_pen_gaf16(gaf16.ad$conc[gaf16.ad$chem == 'pendimethalin']) - 
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
  gafrep = data.frame('but.conc'=lc.but.vals[1:4],
                      'gly.conc'=lc.gly.vals[1:4],
                      'pen.conc'=lc.pen.vals[1:4],
                      'but.rep'=c(44.231, 7.05, 4.52, 4.04)/44.231,
                      'gly.rep'=c(44.231, 4.87, 4.20, 4.23)/44.231,
                      'pen.rep'=c(44.231, 5.18, 4.75, 4.27)/44.231)
  
  plot(gafrep$gly.conc, gafrep$gly.rep, pch = 16, col = 'firebrick', xlim = c(0,10000), ylim = c(0,1),
       xlab = 'Herbicide Concentration', ylab = 'relative reproduction (R0)')
    points(gafrep$but.conc, gafrep$but.rep, pch = 16, col = 'darkorange')
    points(gafrep$pen.conc, gafrep$pen.rep, pch = 16, col = 'indianred2')
    
#Fit function for butralin
  but.fN.pred = nls(but.rep ~ exp(-b*but.conc), data=gafrep, start = list(b=0.001))
    summary(but.fN.pred)
    
    f_N_but_gaf16 = function(He){
      exp(-summary(but.fN.pred)$parameters[1]*He) 
    }
    
    but.fN.pred.l2 = drm(gafrep$but.rep  ~ gafrep$but.conc, type = 'binomial', fct = LL2.2())
    
    f_N_but_gaf16.l2 = function(He){
      (1/(1+exp(but.fN.pred.l2$coefficients[1]*(log(He)-but.fN.pred.l2$coefficients[2]))))
    }
    
    lines(c(0:10000), f_N_but_gaf16(c(0:10000)), lty=2, col='darkorange')  
    lines(c(0:10000), f_N_but_gaf16.l2(c(0:10000)), lty=3, col='darkorange')   
    
#Fit function for glyphosate
  gly.fN.pred = nls(gly.rep ~ exp(-b*gly.conc), data=gafrep, start = list(b=0.001))
    summary(gly.fN.pred)
    
    f_N_gly_gaf16 = function(He){
      exp(-summary(gly.fN.pred)$parameters[1]*He) 
    }
    
    gly.fN.pred.l2 = drm(gafrep$gly.rep  ~ gafrep$gly.conc, type = 'binomial', fct = LL2.2())
    
    f_N_gly_gaf16.l2 = function(He){
      (1/(1+exp(gly.fN.pred.l2$coefficients[1]*(log(He)-gly.fN.pred.l2$coefficients[2]))))
    }
    
    lines(c(0:10000), f_N_gly_gaf16(c(0:10000)), lty=2, col='firebrick')  
    lines(c(0:10000), f_N_gly_gaf16.l2(c(0:10000)), lty=3, col='firebrick')   
    
#Fit function for pendimethalin
  pen.fN.pred = nls(pen.rep ~ exp(-b*pen.conc), data=gafrep, start = list(b=0.001))
    summary(pen.fN.pred)
    
    f_N_pen_gaf16 = function(He){
      exp(-summary(pen.fN.pred)$parameters[1]*He) 
    }
    
    pen.fN.pred.l2 = drm(gafrep$pen.rep  ~ gafrep$pen.conc, type = 'binomial', fct = LL2.2())
    
    f_N_pen_gaf16.l2 = function(He){
      (1/(1+exp(pen.fN.pred.l2$coefficients[1]*(log(He)-pen.fN.pred.l2$coefficients[2]))))
    }
    
    lines(c(0:10000), f_N_pen_gaf16(c(0:10000)), lty=2, col='firebrick')  
    lines(c(0:10000), f_N_pen_gaf16.l2(c(0:10000)), lty=3, col='firebrick')   