require(drc)

#Cercarial mortality (S. mansoni) from Tchounwou 1992 ##################
cerc = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/all_cercarial_mortality.csv')
  cerc.t = subset(cerc, Study == 'tchounwou92')
time = seq(0,25,0.1)

#Reference from Rohr 2008 ################
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
auc.rohr.ctrl.7hr=integrate(f = rohr.cerc.ctrl.surv, lower=0, upper=7)[1]$value

#Tchounwou Data plotted ############
plot(cerc.t$time_hrs[cerc.t$conc==0], cerc.t$surv[cerc.t$conc==0]/100, pch=17, xlab = 'time(hrs',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 1:length(unique(cerc.t$conc[cerc.t$chem == 'malathion']))){
    points(cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[i]], 
           cerc.t$surv[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[i]]/100, pch=16,
           col = 1+i)
  }
#Fit to tchounwou control ########
#Malathion experiment control points
  cerc.mal.ctrl = drm(cerc.t$mort[cerc.t$id == '5_7']/100 ~ 
                      cerc.t$time_hrs[cerc.t$id == '5_7'], 
                   type = 'binomial', fct = LL2.2())
  
  tch.cerc.mal0.surv = function(t){
    1-(1/(1+exp(cerc.mal.ctrl$coefficients[1]*(log(t)-cerc.mal.ctrl$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.cerc.mal0.surv(time), lty=2)
  
  auc.tch.mal0=integrate(f = tch.cerc.mal0.surv, lower=0, upper=24)[1]$value
  auc.tch.mal0.7hr = integrate(f = tch.cerc.mal0.surv, lower=0, upper=7)[1]$value
  
#Bayluscide experiment control points  
  cerc.bay.ctrl = drm(cerc.t$mort[cerc.t$id == '5_1']/100 ~ 
                      cerc.t$time_hrs[cerc.t$id == '5_1'], 
                      type = 'binomial', fct = LL2.2())
  
  tch.cerc.bay0.surv = function(t){
    1-(1/(1+exp(cerc.bay.ctrl$coefficients[1]*(log(t)-cerc.bay.ctrl$coefficients[2]))))
  } #Extremeley low mortality leads to poor model fit
  
  lines(x=time, y=tch.cerc.bay0.surv(time), lty=2)
  
  auc.tch.bay0=integrate(f = tch.cerc.bay0.surv, lower=0, upper=24)[1]$value
  auc.tch.bay0.7hr = integrate(f = tch.cerc.bay0.surv, lower=0, upper=7)[1]$value

#Fit log logistic model to data #######

#Malathion = 50000 ppb
  cerc.mal50 = drm(cerc.t$mort[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[1]]/100 ~ 
                      cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[1]], 
                    type = 'binomial', fct = LL2.2())
  
  tch.cerc.mal50.surv = function(t){
    1-(1/(1+exp(cerc.mal50$coefficients[1]*(log(t)-cerc.mal50$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.cerc.mal50.surv(time), lty=2, col=2)
  
  auc.tch.mal50=integrate(f = tch.cerc.mal50.surv, lower=0, upper=24)[1]$value
  auc.tch.mal50.7hr = integrate(f = tch.cerc.mal50.surv, lower=0, upper=7)[1]$value
  
#Malathion = 100000 ppb
  cerc.mal100 = drm(cerc.t$mort[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[2]]/100 ~ 
                   cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[2]], 
                  type = 'binomial', fct = LL2.2())
  
  tch.cerc.mal100.surv = function(t){
    1-(1/(1+exp(cerc.mal100$coefficients[1]*(log(t)-cerc.mal100$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.cerc.mal100.surv(time), lty=2, col=3)
  
  auc.tch.mal100=integrate(f = tch.cerc.mal100.surv, lower=0, upper=24)[1]$value
  auc.tch.mal100.7hr = integrate(f = tch.cerc.mal100.surv, lower=0, upper=7)[1]$value
  
#Malathion = 150000
  cerc.mal150 = drm(cerc.t$mort[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[3]]/100 ~ 
                      cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[3]], 
                    type = 'binomial', fct = LL2.2())
  
  tch.cerc.mal150.surv = function(t){
    1-(1/(1+exp(cerc.mal150$coefficients[1]*(log(t)-cerc.mal150$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.cerc.mal150.surv(time), lty=2, col=4)
  
  auc.tch.mal150=integrate(f = tch.cerc.mal150.surv, lower=0, upper=24)[1]$value
  auc.tch.mal150.7hr = integrate(f = tch.cerc.mal150.surv, lower=0, upper=7)[1]$value
  
#Malathion = 200000
  cerc.mal200 = drm(cerc.t$mort[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[4]]/100 ~ 
                      cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[4]], 
                    type = 'binomial', fct = LL2.2())
  
  tch.cerc.mal200.surv = function(t){
    1-(1/(1+exp(cerc.mal200$coefficients[1]*(log(t)-cerc.mal200$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.cerc.mal200.surv(time), lty=2, col=5)
  
  auc.tch.mal200=integrate(f = tch.cerc.mal200.surv, lower=0, upper=24)[1]$value
  auc.tch.mal200.7hr = integrate(f = tch.cerc.mal200.surv, lower=0, upper=7)[1]$value
  
#Malathion = 250000
  cerc.mal250 = drm(cerc.t$mort[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[5]]/100 ~ 
                      cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[5]], 
                    type = 'binomial', fct = LL2.2())
  
  tch.cerc.mal250.surv = function(t){
    1-(1/(1+exp(cerc.mal250$coefficients[1]*(log(t)-cerc.mal250$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.cerc.mal250.surv(time), lty=2, col=6)
  
  auc.tch.mal250=integrate(f = tch.cerc.mal250.surv, lower=0, upper=24)[1]$value  
  auc.tch.mal250.7hr = integrate(f = tch.cerc.mal250.surv, lower=0, upper=7)[1]$value
  
#Derive functional response of pi_C to malathion concentration #############
  lines(x=time, y=rohr.cerc.ctrl.surv(time), lty=3)

#Check relative decrease in AUC across malathion concentration
  rel.mal.auc = c(auc.tch.mal0, auc.tch.mal50, auc.tch.mal100, auc.tch.mal150,
                  auc.tch.mal200, auc.tch.mal250)/auc.tch.mal0
  rel.mal.auc.7hr = c(auc.tch.mal0.7hr, auc.tch.mal50.7hr, auc.tch.mal100.7hr, 
                      auc.tch.mal150.7hr, auc.tch.mal200.7hr, auc.tch.mal250.7hr) / 
                    auc.tch.mal0.7hr
  mal.con = c(0,50,100,150,200,250)
  mal.df = data.frame('conc' = mal.con,
                      'auc24' = rel.mal.auc,
                      'auc7' = rel.mal.auc.7hr)
  
  plot(mal.con, rel.mal.auc, pch = 16, ylim=c(0,1), ylab = 'relative auc',
       xlab = 'malathion conentration (ppm)')
  
  mal.piC.predict = nls(auc24 ~ exp(-b*(conc+1e-6)), data = mal.df,
                        start = list(b=0.1))
    summary(mal.piC.predict)

  dose = c(0:250)*1000
  
  pi_C_mal_tch92 = function(q){
    exp(-summary(mal.piC.predict)$parameters[1]*(q/1000)) 
  } #Model derived in ppm; convert concentration to ppb for response
  
  lines(dose/1000, pi_C_mal_tch92(dose), lty=2, col='red')
  title(main = expression(paste(pi[C], '(malathion) from Tchounwou92 data')))
  legend('topright', cex = 0.7,
         legend = paste('b= ', round(summary(mal.piC.predict)$parameters[1],3)))