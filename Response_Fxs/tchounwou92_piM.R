require(drc)

#mirarial mortality (S. mansoni) from Tchounwou 1992 ##################
mir = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Miracidial Mortality/tchounwou08_miracidia.csv')
  mir.mal = subset(mir, chem == 'mal')
  mir.mal$mort = 1 - mir.mal$prop_surv
  time = seq(0,25,0.1)

#Tchounwou Data plotted ############
plot(mir.mal$time[mir.mal$conc_ppm==0], mir.mal$prop_surv[mir.mal$conc_ppm==0], pch=17, xlab = 'time(hrs',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 2:length(unique(mir.mal$conc_ppm))){
    points(mir.mal$time[mir.mal$conc==unique(mir.mal$conc_ppm)[i]], 
           mir.mal$prop_surv[mir.mal$conc==unique(mir.mal$conc_ppm)[i]], pch=16,
           col = 1+i)
  }
#Fit to tchounwou control ########
#Malathion experiment control points
  mir.mal.ctrl = drm(mir.mal$mort[mir.mal$conc_ppm==0] ~ 
                      mir.mal$time[mir.mal$conc_ppm==0], 
                   type = 'binomial', fct = LL2.2())
  
  tch.mir.mal0.surv = function(t){
    1-(1/(1+exp(mir.mal.ctrl$coefficients[1]*(log(t)-mir.mal.ctrl$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.mir.mal0.surv(time), lty=2)
  
  auc.tch.mir.mal0=integrate(f = tch.mir.mal0.surv, lower=0, upper=24)[1]$value
  auc.tch.mir.mal0.6hr = integrate(f = tch.mir.mal0.surv, lower=0, upper=6)[1]$value
  
#Fit log logistic model to data #######

#Malathion = 30000 ppb
  mir.mal30 = drm(mir.mal$mort[mir.mal$conc_ppm == 30] ~ mir.mal$time[mir.mal$conc_ppm == 30], 
                  type = 'binomial', fct = LL2.2())
  
  tch.mir.mal30.surv = function(t){
    1-(1/(1+exp(mir.mal30$coefficients[1]*(log(t)-mir.mal30$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.mir.mal30.surv(time), lty=2, col=3)
  
  auc.tch.mir.mal30=integrate(f = tch.mir.mal30.surv, lower=0, upper=24)[1]$value
  auc.tch.mir.mal30.6hr = integrate(f = tch.mir.mal30.surv, lower=0, upper=6)[1]$value
  
#Malathion = 60000 ppb
  mir.mal60 = drm(mir.mal$mort[mir.mal$conc_ppm==60] ~ mir.mal$time[mir.mal$conc_ppm==60], 
                  type = 'binomial', fct = LL2.2())
  
  tch.mir.mal60.surv = function(t){
    1-(1/(1+exp(mir.mal60$coefficients[1]*(log(t)-mir.mal60$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.mir.mal60.surv(time), lty=2, col=4)
  
  auc.tch.mir.mal60=integrate(f = tch.mir.mal60.surv, lower=0, upper=24)[1]$value
  auc.tch.mir.mal60.6hr = integrate(f = tch.mir.mal60.surv, lower=0, upper=6)[1]$value
  
#Malathion = 90000
  mir.mal90 = drm(mir.mal$mort[mir.mal$conc_ppm==90] ~ mir.mal$time[mir.mal$conc_ppm==90], 
                    type = 'binomial', fct = LL2.2())
  
  tch.mir.mal90.surv = function(t){
    1-(1/(1+exp(mir.mal90$coefficients[1]*(log(t)-mir.mal90$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.mir.mal90.surv(time), lty=2, col=5)
  
  auc.tch.mir.mal90=integrate(f = tch.mir.mal90.surv, lower=0, upper=24)[1]$value
  auc.tch.mir.mal90.6hr = integrate(f = tch.mir.mal90.surv, lower=0, upper=6)[1]$value
  
#Malathion = 120000
  mir.mal120 = drm(mir.mal$mort[mir.mal$conc_ppm == 120] ~ mir.mal$time[mir.mal$conc_ppm == 120], 
                    type = 'binomial', fct = LL2.2())
  
  tch.mir.mal120.surv = function(t){
    1-(1/(1+exp(mir.mal120$coefficients[1]*(log(t)-mir.mal120$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.mir.mal120.surv(time), lty=2, col=6)
  
  auc.tch.mir.mal120=integrate(f = tch.mir.mal120.surv, lower=0, upper=24)[1]$value
  auc.tch.mir.mal120.6hr = integrate(f = tch.mir.mal120.surv, lower=0, upper=6)[1]$value
  
#Malathion = 150000
  mir.mal150 = drm(mir.mal$mort[mir.mal$conc_ppm == 150] ~ mir.mal$time[mir.mal$conc_ppm == 150], 
                    type = 'binomial', fct = LL2.2())
  
  tch.mir.mal150.surv = function(t){
    1-(1/(1+exp(mir.mal150$coefficients[1]*(log(t)-mir.mal150$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.mir.mal150.surv(time), lty=2, col=7)
  
  auc.tch.mir.mal150=integrate(f = tch.mir.mal150.surv, lower=0, upper=24)[1]$value  
  auc.tch.mir.mal150.6hr = integrate(f = tch.mir.mal150.surv, lower=0, upper=6)[1]$value
  
#Check relative decrease in AUC across malathion concentration
  rel.mal.mir.auc = c(auc.tch.mir.mal0, auc.tch.mir.mal30, auc.tch.mir.mal60, auc.tch.mir.mal90,
                  auc.tch.mir.mal120, auc.tch.mir.mal150)/auc.tch.mir.mal0
  rel.mal.mir.auc.6hr = c(auc.tch.mir.mal0.6hr, auc.tch.mir.mal30.6hr, auc.tch.mir.mal60.6hr, auc.tch.mir.mal90.6hr,
                      auc.tch.mir.mal120.6hr, auc.tch.mir.mal150.6hr) / auc.tch.mir.mal0.6hr
  mal.con = c(0,30,60,90,120,150)
  mal.df = data.frame('conc' = mal.con,
                      'auc24' = rel.mal.mir.auc,
                      'auc6' = rel.mal.mir.auc.6hr)
  
  plot(mal.con, rel.mal.mir.auc, pch = 16, ylim=c(0,1), ylab = 'relative auc',
       xlab = 'malathion conentration (ppm)')
  
  mal.piM.predict = lm(auc24 ~ conc + 0, data = mal.df)
  
  mal.piM.predict = nls(auc24 ~ -m*(conc)+1, data = mal.df,
                        start = list(m=0.1))
    summary(mal.piM.predict)

  dose = c(0:150)*1000
  
  pi_M_mal_tch92 = function(q){
    -summary(mal.piM.predict)$parameters[1]*(q/1000)+1 
  } #Model derived in ppm; convert concentration to ppb for response
  
  lines(dose/1000, pi_M_mal_tch92(dose), lty=2, col='red')
  title(main = expression(paste(pi[M], '(malathion) from Tchounwou92 data')))
  