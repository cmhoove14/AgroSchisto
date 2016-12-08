require(drc)

#Cercarial mortality (E. trivolvis) from Koprivnikar 2006 ##################
cerc = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/all_cercarial_mortality.csv')
  cerc.k = subset(cerc, Study == 'koprivnikar2006')
time = c(0:25)

plot(x = cerc.k$time_hrs[cerc.k$chem == 'control'], 
     y = cerc.k$surv[cerc.k$chem == 'control']/100,
     xlab = 'time (hrs)', ylab = 'prop dead', pch = 16)

cerc.k.ctrl = drm(cerc.k$mort[cerc.k$chem == 'control']/100 ~ 
                    cerc.k$time_hrs[cerc.k$chem == 'control'], 
                type = 'binomial', fct = LL2.2())

kop.cerc.ctrl.surv = function(t){
  1-(1/(1+exp(cerc.k.ctrl$coefficients[1]*(log(t)-cerc.k.ctrl$coefficients[2]))))
} 

lines(x=time, y=kop.cerc.ctrl.surv(time), col='grey40')

auc.kop.ctrl=integrate(f = kop.cerc.ctrl.surv, lower=0, upper=24)[1]$value

#20 ppb atrazine ##################
points(x = cerc.k$time_hrs[cerc.k$conc == 20], 
     y = cerc.k$surv[cerc.k$conc == 20]/100, pch = 17)

cerc.k.20 = drm(cerc.k$mort[cerc.k$conc == 20]/100 ~ 
                    cerc.k$time_hrs[cerc.k$conc == 20], 
                  type = 'binomial', fct = LL2.2())

kop.cerc.20 = function(t){
  1-(1/(1+exp(cerc.k.20$coefficients[1]*(log(t)-cerc.k.20$coefficients[2]))))
} 

lines(x=time, y=kop.cerc.20(time), lty=2)

auc.kop.atr20=integrate(f = kop.cerc.20, lower=0, upper=24)[1]$value

#200 ppb atrazine ##################
points(x = cerc.k$time_hrs[cerc.k$conc == 200], 
       y = cerc.k$surv[cerc.k$conc == 200]/100, pch = 15)

cerc.k.200 = drm(cerc.k$mort[cerc.k$conc == 200]/100 ~ 
                  cerc.k$time_hrs[cerc.k$conc == 200], 
                type = 'binomial', fct = LL2.2())

kop.cerc.200 = function(t){
  1-(1/(1+exp(cerc.k.200$coefficients[1]*(log(t)-cerc.k.200$coefficients[2]))))
} 

lines(x=time, y=kop.cerc.200(time), lty=3)

auc.kop.atr200=integrate(f = kop.cerc.200, lower=0, upper=24)[1]$value

#Summaries #############
title(main='Koprivnikar2006 Cercarial mortality (E.trivolvis)')
legend('topright', legend = c('control', '20ppb', '200ppb'), pch = c(16,17,15), cex=0.8)

rel.auc.kop = c(auc.kop.ctrl, auc.kop.atr20, auc.kop.atr200)/auc.kop.ctrl

plot(c(0,20,200), rel.auc.kop, pch = 16, ylim=c(0,1), ylab = 'relative auc',
     xlab = 'atrazine conentration (ppb)')

kop.df = data.frame('auc24' = rel.auc.kop,
                    'conc' = c(0,20,200))

atr.piC.k.predict = nls(auc24 ~ exp(-b*(conc+1e-6)), data = kop.df,
                      start = list(b=0.1))
  summary(atr.piC.k.predict)

atr.k.con = c(0:201)

pi_C_atr_kop06 = function(q){
  exp(-summary(atr.piC.k.predict)$parameters[1]*(q)) 
} #Model derived in ppm; convert concentration to ppb for response

lines(atr.k.con, pi_C_atr_kop06(atr.k.con), lty=2, col='red')
title(main = expression(paste(pi[C], '(atrazine) from Koprivnikar06')))
legend('topright', cex = 0.7,
       legend = paste('b= ', round(summary(atr.piC.k.predict)$parameters[1],3)))