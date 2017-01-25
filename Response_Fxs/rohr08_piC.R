require(drc)

cerc = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/all_cercarial_mortality.csv')
  time = c(0:25)

#Plot control survival curve, estimate time-dependent die-off function, and get auc from 0-24 hours ########
plot(x = cerc$time_hrs[cerc$chem == 'control' & cerc$Study == 'rohr2008'], 
     y = cerc$surv[cerc$chem == 'control' & cerc$Study == 'rohr2008']/100,
     xlab = 'time (hrs)', ylab = 'prop dead', pch = 16)

cerc.ctrl = drm(cerc$mort[cerc$chem == 'control' & cerc$Study == 'rohr2008']/100 ~ 
                  cerc$time_hrs[cerc$chem == 'control' & cerc$Study == 'rohr2008'], 
                type = 'binomial', fct = LL2.2())

rohr.cerc.ctrl.surv = function(t){
  1-(1/(1+exp(cerc.ctrl$coefficients[1]*(log(t)-cerc.ctrl$coefficients[2]))))
} 

  lines(x=time, y=rohr.cerc.ctrl.surv(time), lty=2)
  
  auc.rohr.ctrl=integrate(f = rohr.cerc.ctrl.surv, lower=0, upper=24)[1]$value
  
#Atrazine ########
points(x = cerc$time_hrs[cerc$chem == 'atrazine' & cerc$Study == 'rohr2008'], 
       y = cerc$surv[cerc$chem == 'atrazine' & cerc$Study == 'rohr2008']/100,
       col = 'gold', pch = 16)

  cerc.atr = drm(cerc$mort[cerc$chem == 'atrazine' & cerc$Study == 'rohr2008']/100 ~ 
                    cerc$time_hrs[cerc$chem == 'atrazine' & cerc$Study == 'rohr2008'], 
                  type = 'binomial', fct = LL2.2())
  rohr.cerc.atr = function(t){
    1 - (1/(1+exp(cerc.atr$coefficients[1]*(log(t)-cerc.atr$coefficients[2]))))
  }  
  
  lines(x=time, y=rohr.cerc.atr(time), lty=2, col = 'gold')  
  
  auc.rohr.atr201=integrate(f = rohr.cerc.atr, lower=0, upper=24)[1]$value
  
#Malathion ########
points(x = cerc$time_hrs[cerc$chem == 'malathion' & cerc$Study == 'rohr2008'], 
       y = cerc$surv[cerc$chem == 'malathion' & cerc$Study == 'rohr2008']/100,
       col = 'red', pch = 16)
  
  cerc.mal = drm(cerc$mort[cerc$chem == 'malathion' & cerc$Study == 'rohr2008']/100 ~ 
                   cerc$time_hrs[cerc$chem == 'malathion' & cerc$Study == 'rohr2008'], 
                 type = 'binomial', fct = LL2.2())
  rohr.cerc.mal = function(t){
    1 - (1/(1+exp(cerc.mal$coefficients[1]*(log(t)-cerc.mal$coefficients[2]))))
  }  
  
  lines(x=time, y=rohr.cerc.mal(time), lty=2, col = 'red')  
  
  auc.rohr.mal9.6=integrate(f = rohr.cerc.mal, lower=0, upper=24)[1]$value
#carbaryl ########
points(x = cerc$time_hrs[cerc$chem == 'carbaryl' & cerc$Study == 'rohr2008'], 
       y = cerc$surv[cerc$chem == 'carbaryl' & cerc$Study == 'rohr2008']/100,
       col = 'purple', pch = 16)
  
  cerc.car = drm(cerc$mort[cerc$chem == 'carbaryl' & cerc$Study == 'rohr2008']/100 ~ 
                   cerc$time_hrs[cerc$chem == 'carbaryl' & cerc$Study == 'rohr2008'], 
                 type = 'binomial', fct = LL2.2())
  rohr.cerc.car = function(t){
    1 - (1/(1+exp(cerc.car$coefficients[1]*(log(t)-cerc.car$coefficients[2]))))
  }  
  
  lines(x=time, y=rohr.cerc.car(time), lty=2, col = 'purple')  
  
  auc.rohr.car33.5=integrate(f = rohr.cerc.car, lower=0, upper=24)[1]$value  
  
#glyphosate ########
points(x = cerc$time_hrs[cerc$chem == 'glyphosate' & cerc$Study == 'rohr2008'], 
       y = cerc$surv[cerc$chem == 'glyphosate' & cerc$Study == 'rohr2008']/100,
       col = 'green', pch = 16)
  
  cerc.gly = drm(cerc$mort[cerc$chem == 'glyphosate' & cerc$Study == 'rohr2008']/100 ~ 
                   cerc$time_hrs[cerc$chem == 'glyphosate' & cerc$Study == 'rohr2008'], 
                 type = 'binomial', fct = LL2.2())
  rohr.cerc.gly = function(t){
    1 - (1/(1+exp(cerc.gly$coefficients[1]*(log(t)-cerc.gly$coefficients[2]))))
  }  
  
  lines(x=time, y=rohr.cerc.gly(time), lty=2, col = 'green')  
  
  auc.rohr.gly3700=integrate(f = rohr.cerc.gly, lower=0, upper=24)[1]$value  
  
#Summaries ##########
  legend('bottomleft', legend = c('control', 'atrazine', 'carbaryl', 'glyphosate', 'malathion'),
         col = c('black', 'gold', 'purple', 'green', 'red'), pch = 16, cex=0.8)
  title(main = 'Rohr 2008 Cercarial mortality')
  legend('topright', title = 'auc', 
         legend = c(paste('control= ', round(auc.rohr.ctrl,2), sep=''), 
                    paste('atrazine= ', round(auc.rohr.atr201,2), sep=''), 
                    paste('carbaryl= ', round(auc.rohr.car33.5,2), sep=''), 
                    paste('glyphosate= ', round(auc.rohr.gly3700,2), sep=''), 
                    paste('malathion= ', round(auc.rohr.mal9.6,2), sep='')), cex=0.7)
  
  rel.auc.rohr = c('ctrl' = auc.rohr.ctrl/auc.rohr.ctrl,
                   'atr201' = auc.rohr.atr201/auc.rohr.ctrl,
                   'carb33.5' = auc.rohr.car33.5/auc.rohr.ctrl,
                   'gly3700' = auc.rohr.gly3700/auc.rohr.ctrl,
                   'mal9.6' = auc.rohr.mal9.6/auc.rohr.ctrl)
#Atrazine "model" (yes, I know there are only two data points)  ###########
  atr.dat = data.frame('conc' = c(0,201),
                       'auc24' = c(rel.auc.rohr[1], rel.auc.rohr[2]))
  
    atr.m<-nls(auc24 ~ exp(-b*(conc+1e-6)), data = atr.dat,
               start = list(b=0.01))
      summary(atr.m)
    
    atr.con = c(0:201)
    
    pi_C_atr_rohr08 = function(He){
      exp(-summary(atr.m)$parameters[1]*(He)) 
    } 
    
    plot(atr.dat$conc, atr.dat$auc24, pch = 16, xlab = 'Atrazine concentration', ylab = '24-hr AUC', ylim = c(0,1))
      lines(atr.con, pi_C_atr_rohr08(atr.con), lty=2, col='red')
      
#Carbaryl "model" (yes, I know there are only two data points)  ###########
  carb.dat = data.frame('conc' = c(0,33.5),
                       'auc24' = c(rel.auc.rohr[1], rel.auc.rohr[3]))
  
    carb.m<-nls(auc24 ~ exp(-b*(conc+1e-6)), data = carb.dat,
               start = list(b=0.01))
      summary(carb.m)
    
    carb.con = c(0:40)
    
    pi_C_carb_rohr08 = function(In){
      exp(-summary(carb.m)$parameters[1]*(In)) 
    } 
    
    plot(carb.dat$conc, carb.dat$auc24, pch = 16, xlab = 'Carbaryl concentration', ylab = '24-hr AUC', ylim = c(0,1))
      lines(carb.con, pi_C_carb_rohr08(carb.con), lty=2, col='red')
    
#glyphosate "model" (yes, I know there are only two data points)  ###########
  gly.dat = data.frame('conc' = c(0,3700),
                       'auc24' = c(rel.auc.rohr[1], rel.auc.rohr[4]))
  
  gly.m<-nls(auc24 ~ exp(-b*(conc+1e-6)), data = gly.dat,
              start = list(b=0.001))
    summary(gly.m)
  
  gly.con = c(0:4000)
  
  pi_C_gly_rohr08 = function(He){
    exp(-summary(gly.m)$parameters[1]*(He)) 
  } 
  
  plot(gly.dat$conc, gly.dat$auc24, pch = 16, xlab = 'glyphosate concentration', ylab = '24-hr AUC', ylim = c(0,1))
    lines(gly.con, pi_C_gly_rohr08(gly.con), lty=2, col='red')
    
#Malathion "model" (yes, I know there are only two data points)  ###########
    
  mal.dat = data.frame('conc' = c(0,9.6),
                       'auc24' = c(rel.auc.rohr[1], rel.auc.rohr[5]))
  
  mal.m<-nls(auc24 ~ exp(-b*(conc+1e-6)), data = mal.dat,
             start = list(b=0.001))
    summary(mal.m)
  
  mal.con = seq(0,10,0.1)
  
  pi_C_mal_rohr08 = function(In){
    exp(-summary(mal.m)$parameters[1]*(In)) 
  } 
  
  plot(mal.dat$conc, mal.dat$auc24, pch = 16, xlab = 'Malathion concentration', ylab = '24-hr AUC', ylim = c(0,1))
    lines(mal.con, pi_C_mal_rohr08(mal.con), lty=2, col='red')