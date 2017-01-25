#Data extraction and model fitting to Ghaffar 2016 LARVAL (S.mansoni) data
require(drc)

  but.vals = c(0, 556, 2417, 3906, 5560, 8703)
  gly.vals = c(0, 1506, 3875, 9174, 15062, 26249)
  pen.vals = c(0, 214.8, 535, 1299, 2148, 3762)

#Miracidial (S. mansoni) toxicity ###############
  mir<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Miracidial Mortality/ghaffar2016_miracidia.csv')

#butralin ###############
  plot(mir$time_hrs[mir$conc==0], mir$surv[mir$conc==0], 
       pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 1:length(unique(mir$conc[mir$chem == 'butralin']))){
    points(mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[i] & mir$chem == 'butralin'], 
           mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[i] & mir$chem == 'butralin'], pch=17,
           col = i+1)
  }

#fit to control points *****************************************************************
  gaf.mir.cont<-drm(mir$surv[mir$conc==0] ~ mir$time_hrs[mir$conc==0],
                    type = 'binomial', fct = LL2.2())
  
  gaf.mir.cont.surv = function(t){
    (1/(1+exp(gaf.mir.cont$coefficients[1]*(log(t)-gaf.mir.cont$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.cont.surv(seq(0,24,0.1)), lty=2)
  
  auc.mir.cont=integrate(f = gaf.mir.cont.surv, lower=0, upper=24)[1]$value  

#556 ppb *********************************************************************************
  gaf.mir.but556<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[1] & mir$chem == 'butralin'] ~ 
                        mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[1] & mir$chem == 'butralin'],
                      type = 'binomial', fct = LL2.2())
  
  gaf.mir.but556.surv = function(t){
    (1/(1+exp(gaf.mir.but556$coefficients[1]*(log(t)-gaf.mir.but556$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.but556.surv(seq(0,24,0.1)), lty=2, col=2)
  
  auc.mir.but556=integrate(f = gaf.mir.but556.surv, lower=0, upper=24)[1]$value  

#2417 ppb *********************************************************************************
  gaf.mir.but2417<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[2] & mir$chem == 'butralin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[2] & mir$chem == 'butralin'],
                       type = 'binomial', fct = LL2.2())
  
  gaf.mir.but2417.surv = function(t){
    (1/(1+exp(gaf.mir.but2417$coefficients[1]*(log(t)-gaf.mir.but2417$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.but2417.surv(seq(0,24,0.1)), lty=2, col=3)
  
  auc.mir.but2417=integrate(f = gaf.mir.but2417.surv, lower=0, upper=24)[1]$value  

#3906 ppb ***********************************************************************************
  gaf.mir.but3906<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[3] & mir$chem == 'butralin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[3] & mir$chem == 'butralin'],
                       type = 'binomial', fct = LL2.2())
  
  gaf.mir.but3906.surv = function(t){
    (1/(1+exp(gaf.mir.but3906$coefficients[1]*(log(t)-gaf.mir.but3906$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.but3906.surv(seq(0,24,0.1)), lty=2, col=4)
  
  auc.mir.but3906=integrate(f = gaf.mir.but3906.surv, lower=0, upper=24)[1]$value  

#5560 ppb **********************************************************************************
  gaf.mir.but5560<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[4] & mir$chem == 'butralin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[4] & mir$chem == 'butralin'],
                       type = 'binomial', fct = LL2.2())
  
  gaf.mir.but5560.surv = function(t){
    (1/(1+exp(gaf.mir.but5560$coefficients[1]*(log(t)-gaf.mir.but5560$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.but5560.surv(seq(0,24,0.1)), lty=2, col=5)
  
  auc.mir.but5560=integrate(f = gaf.mir.but5560.surv, lower=0, upper=24)[1]$value     

#8703 ppb **********************************************************************************
  gaf.mir.but8703<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[5] & mir$chem == 'butralin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butralin'])[5] & mir$chem == 'butralin'],
                       type = 'binomial', fct = LL2.2())
  
  gaf.mir.but8703.surv = function(t){
    (1/(1+exp(gaf.mir.but8703$coefficients[1]*(log(t)-gaf.mir.but8703$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.but8703.surv(seq(0,24,0.1)), lty=2, col=6)
  
  auc.mir.but8703=integrate(f = gaf.mir.but8703.surv, lower=0, upper=24)[1]$value   
  
  title('Ghaffar2016 butralin toxicity to miracidia')
  legend('topright', legend = c('control', 556,2417,3906,5560,8703), pch = c(16,17,17,17,17,17), col = c(1:6), cex=0.8)
  
  gaf.but.aucs = c(auc.mir.cont,auc.mir.but556,auc.mir.but2417,auc.mir.but3906,auc.mir.but5560,auc.mir.but8703)/auc.mir.cont
#glyphosate ###############
  plot(mir$time_hrs[mir$conc==0], mir$surv[mir$conc==0], 
       pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 1:length(unique(mir$conc[mir$chem == 'glyphosate']))){
    points(mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[i] & mir$chem == 'glyphosate'], 
           mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[i] & mir$chem == 'glyphosate'], pch=17,
           col = i+1)
  }
    lines(x=seq(0,24,0.1), y=gaf.mir.cont.surv(seq(0,24,0.1)), lty=2)

#1506 ppb *********************************************************************************
  gaf.mir.gly1506<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[1] & mir$chem == 'glyphosate'] ~ 
                        mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[1] & mir$chem == 'glyphosate'],
                      type = 'binomial', fct = LL2.2())
  
  gaf.mir.gly1506.surv = function(t){
    (1/(1+exp(gaf.mir.gly1506$coefficients[1]*(log(t)-gaf.mir.gly1506$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.gly1506.surv(seq(0,24,0.1)), lty=2, col=2)
  
  auc.mir.gly1506=integrate(f = gaf.mir.gly1506.surv, lower=0, upper=24)[1]$value  
  
#3875 ppb *********************************************************************************
  gaf.mir.gly3875<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[2] & mir$chem == 'glyphosate'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[2] & mir$chem == 'glyphosate'],
                       type = 'binomial', fct = LL2.2())
  
  gaf.mir.gly3875.surv = function(t){
    (1/(1+exp(gaf.mir.gly3875$coefficients[1]*(log(t)-gaf.mir.gly3875$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.gly3875.surv(seq(0,24,0.1)), lty=2, col=3)
  
  auc.mir.gly3875=integrate(f = gaf.mir.gly3875.surv, lower=0, upper=24)[1]$value  
  
#9174 ppb ***********************************************************************************
  gaf.mir.gly9174<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[3] & mir$chem == 'glyphosate'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[3] & mir$chem == 'glyphosate'],
                       type = 'binomial', fct = LL2.2())
  
  gaf.mir.gly9174.surv = function(t){
    (1/(1+exp(gaf.mir.gly9174$coefficients[1]*(log(t)-gaf.mir.gly9174$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.gly9174.surv(seq(0,24,0.1)), lty=2, col=4)
  
  auc.mir.gly9174=integrate(f = gaf.mir.gly9174.surv, lower=0, upper=24)[1]$value  
  
#15062 ppb **********************************************************************************
  gaf.mir.gly15062<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[4] & mir$chem == 'glyphosate'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[4] & mir$chem == 'glyphosate'],
                       type = 'binomial', fct = LL2.2())
  
  gaf.mir.gly15062.surv = function(t){
    (1/(1+exp(gaf.mir.gly15062$coefficients[1]*(log(t)-gaf.mir.gly15062$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.gly15062.surv(seq(0,24,0.1)), lty=2, col=5)
  
  auc.mir.gly15062=integrate(f = gaf.mir.gly15062.surv, lower=0, upper=24)[1]$value     
  
#26249 ppb **********************************************************************************
  gaf.mir.gly26249<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[5] & mir$chem == 'glyphosate'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'glyphosate'])[5] & mir$chem == 'glyphosate'],
                       type = 'binomial', fct = LL2.2())
  
  gaf.mir.gly26249.surv = function(t){
    (1/(1+exp(gaf.mir.gly26249$coefficients[1]*(log(t)-gaf.mir.gly26249$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.gly26249.surv(seq(0,24,0.1)), lty=2, col=6)
  
  auc.mir.gly26249=integrate(f = gaf.mir.gly26249.surv, lower=0, upper=24)[1]$value   
  
  title('Ghaffar2016 glyphosate toxicity to miracidia')
  legend('topright', legend = c('control', 1506,3875,9174,15062,26249), pch = c(16,17,17,17,17,17), col = c(1:6), cex=0.8) 
  
  gaf.gly.aucs = c(auc.mir.cont, auc.mir.gly1506, auc.mir.gly3875, auc.mir.gly9174, auc.mir.gly15062, auc.mir.gly26249)/auc.mir.cont
#pendimethalin ###############
  plot(mir$time_hrs[mir$conc==0], mir$surv[mir$conc==0], 
       pch=16, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 1:length(unique(mir$conc[mir$chem == 'pendimethalin']))){
    points(mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[i] & mir$chem == 'pendimethalin'], 
           mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[i] & mir$chem == 'pendimethalin'], pch=17,
           col = i+1)
  }
  lines(x=seq(0,24,0.1), y=gaf.mir.cont.surv(seq(0,24,0.1)), lty=2)
  
#215 ppb *********************************************************************************
  gaf.mir.pen215<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[1] & mir$chem == 'pendimethalin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[1] & mir$chem == 'pendimethalin'],
                       type = 'binomial', fct = LL2.2())
  
  gaf.mir.pen215.surv = function(t){
    (1/(1+exp(gaf.mir.pen215$coefficients[1]*(log(t)-gaf.mir.pen215$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.pen215.surv(seq(0,24,0.1)), lty=2, col=2)
  
  auc.mir.pen215=integrate(f = gaf.mir.pen215.surv, lower=0, upper=24)[1]$value  
  
#535 ppb *********************************************************************************
  gaf.mir.pen535<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[2] & mir$chem == 'pendimethalin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[2] & mir$chem == 'pendimethalin'],
                       type = 'binomial', fct = LL2.2())
  
  gaf.mir.pen535.surv = function(t){
    (1/(1+exp(gaf.mir.pen535$coefficients[1]*(log(t)-gaf.mir.pen535$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.pen535.surv(seq(0,24,0.1)), lty=2, col=3)
  
  auc.mir.pen535=integrate(f = gaf.mir.pen535.surv, lower=0, upper=24)[1]$value  
  
#1299 ppb ***********************************************************************************
  gaf.mir.pen1299<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[3] & mir$chem == 'pendimethalin'] ~ 
                         mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[3] & mir$chem == 'pendimethalin'],
                       type = 'binomial', fct = LL2.2())
  
  gaf.mir.pen1299.surv = function(t){
    (1/(1+exp(gaf.mir.pen1299$coefficients[1]*(log(t)-gaf.mir.pen1299$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.pen1299.surv(seq(0,24,0.1)), lty=2, col=4)
  
  auc.mir.pen1299=integrate(f = gaf.mir.pen1299.surv, lower=0, upper=24)[1]$value  
  
#2148 ppb **********************************************************************************
  gaf.mir.pen2148<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[4] & mir$chem == 'pendimethalin'] ~ 
                          mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[4] & mir$chem == 'pendimethalin'],
                        type = 'binomial', fct = LL2.2())
  
  gaf.mir.pen2148.surv = function(t){
    (1/(1+exp(gaf.mir.pen2148$coefficients[1]*(log(t)-gaf.mir.pen2148$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.pen2148.surv(seq(0,24,0.1)), lty=2, col=5)
  
  auc.mir.pen2148=integrate(f = gaf.mir.pen2148.surv, lower=0, upper=24)[1]$value     
  
#3762 ppb **********************************************************************************
  gaf.mir.pen3762<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[5] & mir$chem == 'pendimethalin'] ~ 
                          mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'pendimethalin'])[5] & mir$chem == 'pendimethalin'],
                        type = 'binomial', fct = LL2.2())
  
  gaf.mir.pen3762.surv = function(t){
    (1/(1+exp(gaf.mir.pen3762$coefficients[1]*(log(t)-gaf.mir.pen3762$coefficients[2]))))
  } 
  
  lines(x=seq(0,24,0.1), y=gaf.mir.pen3762.surv(seq(0,24,0.1)), lty=2, col=6)
  
  auc.mir.pen3762=integrate(f = gaf.mir.pen3762.surv, lower=0, upper=24)[1]$value   
  
  title('Ghaffar2016 pendimethalin toxicity to miracidia')
  legend('topright', legend = c('control', 215,535,1299,2148,3762), pch = c(16,17,17,17,17,17), col = c(1:6), cex=0.8)  
  
  gaf.pen.aucs = c(auc.mir.cont, auc.mir.pen215, auc.mir.pen535, auc.mir.pen1299, auc.mir.pen2148, auc.mir.pen3762)/auc.mir.cont
#Fit pi_M functions ################  
  mirgaf.df = data.frame('but.conc'=but.vals,
                         'gly.conc'=gly.vals,
                         'pen.conc'=pen.vals,
                         'but.auc'=gaf.but.aucs,
                         'gly.auc'=gaf.gly.aucs,
                         'pen.auc'=gaf.pen.aucs)
  
  plot(mirgaf.df$gly.conc, mirgaf.df$gly.auc, pch = 16, col = 'green', xlim = c(0,27000), ylim = c(0,1),
       xlab = 'Herbicide Concentration', ylab = 'relative auc')
    points(mirgaf.df$but.conc, mirgaf.df$but.auc, pch = 16, col = 'darkolivegreen')
    points(mirgaf.df$pen.conc, mirgaf.df$pen.auc, pch = 16, col = 'cyan4')
#Glyphosate   
  gly.piM.pred = nls(gly.auc ~ exp(-b*gly.conc), data=mirgaf.df, start = list(b=0.0001))
    summary(gly.piM.pred)
    
  pi_M_gly_gaf16 = function(He){
      exp(-summary(gly.piM.pred)$parameters[1]*He) 
    }
    
    lines(c(0:27000), pi_M_gly_gaf16(c(0:27000)), lty=2, col='green')  
    
#butralin   
  but.piM.pred = nls(but.auc ~ exp(-b*but.conc), data=mirgaf.df, start = list(b=0.001))
    summary(but.piM.pred)
    
  pi_M_but_gaf16 = function(He){
      exp(-summary(but.piM.pred)$parameters[1]*He) 
    }
    
    lines(c(0:27000), pi_M_but_gaf16(c(0:27000)), lty=2, col='darkolivegreen')  
    
#pendimethalin   
  pen.piM.pred = nls(pen.auc ~ exp(-b*pen.conc), data=mirgaf.df, start = list(b=0.005))
    summary(pen.piM.pred)
    
  pi_M_pen_gaf16 = function(He){
      exp(-summary(pen.piM.pred)$parameters[1]*He) 
    }
    
    lines(c(0:27000), pi_M_pen_gaf16(c(0:27000)), lty=2, col='cyan4')  