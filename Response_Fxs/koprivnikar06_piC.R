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

ll4 = function(hi,lo,slp,lc,x){
  lo + ((hi-lo)/(1+exp(slp*(log(x)-log(lc)))))
}


require(drc)

#Cercarial mortality (E. trivolvis) from Koprivnikar 2006 ##################
kop.c = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/Koprivnikar2006.csv')
  kop.c[,5:8] = kop.c[,5:8]/100 #Convert survival measures to proportions
  time = c(0:25)
  
kop.cc = subset(kop.c, chem == 'control')  

plot(x = kop.cc$time_hrs, y = kop.cc$surv,
     xlab = 'time (hrs)', ylab = 'prop dead', pch = 16)

  kop.ctrl = nls(surv ~ 1/(1+exp(slp*(log(time_hrs)-log(lc50)))), data = kop.cc,
                 start = list(slp = 5, lc50 = 2.5), weights = surv.se^-1,
                 lower = c(0,0), algorithm = 'port')
  
  lines(time, ll4(1,0,summary(kop.ctrl)$coefficients[1],
                  summary(kop.ctrl)$coefficients[2], time), lty=2)

kop.c.ctrl = drm(kop.cc$surv ~ kop.cc$time_hrs, 
                type = 'binomial', fct = LL2.2())

lines(time, ll4(1,0,coef(kop.c.ctrl)[1], exp(coef(kop.c.ctrl)[2]), time), lty=3)

kop.cerc.ctrl.surv = function(t){
  1-(1/(1+exp(kop.c.ctrl$coefficients[1]*(log(t)-kop.c.ctrl$coefficients[2]))))
} 

lines(x=time, y=kop.cerc.ctrl.surv(time))

auc.kop.ctrl=integrate(f = kop.cerc.ctrl.surv, lower=0, upper=24)[1]$value

#20 ppb atrazine ##################
points(x = kop.c$time_hrs[kop.c$conc == 20], 
       y = kop.c$surv[kop.c$conc == 20]/100, pch = 17)

kop.c.20 = drm(kop.c$mort[kop.c$conc == 20]/100 ~ 
                    kop.c$time_hrs[kop.c$conc == 20], 
                  type = 'binomial', fct = LL2.2())

kop.cerc.20 = function(t){
  1-(1/(1+exp(kop.c.20$coefficients[1]*(log(t)-kop.c.20$coefficients[2]))))
} 

lines(x=time, y=kop.cerc.20(time), lty=2)

auc.kop.atr20=integrate(f = kop.cerc.20, lower=0, upper=24)[1]$value

#200 ppb atrazine ##################
points(x = kop.c$time_hrs[kop.c$conc == 200], 
       y = kop.c$surv[kop.c$conc == 200]/100, pch = 15)

kop.c.200 = drm(kop.c$mort[kop.c$conc == 200]/100 ~ 
                  kop.c$time_hrs[kop.c$conc == 200], 
                type = 'binomial', fct = LL2.2())

kop.cerc.200 = function(t){
  1-(1/(1+exp(kop.c.200$coefficients[1]*(log(t)-kop.c.200$coefficients[2]))))
} 

lines(x=time, y=kop.cerc.200(time), lty=3)

auc.kop.atr200=integrate(f = kop.cerc.200, lower=0, upper=24)[1]$value

  title(main='Koprivnikar2006 Cercarial mortality (E.trivolvis)')
  legend('topright', legend = c('control', '20ppb', '200ppb'), pch = c(16,17,15), cex=0.8)

#Derive functional response of pi_C to atrazine/metolachlor concentration #############
#Check decrease in AUC across atrazine/metolachlor concentration treating AUC as continuous variable
  auc.kop = c(auc.kop.ctrl, auc.kop.atr20, auc.kop.atr200)
  
  kop.df = data.frame('auc24' = auc.kop,
                      'conc' = c(0,20,200))
#Atrazine function  
  atr.piC.k.predict = drm(auc24 ~ conc, data = kop.df, type = 'continuous',
                          fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                                     fixed = c(NA, auc.kop.atr200, auc.kop.ctrl, NA)))
  summary(atr.piC.k.predict)
  plot(atr.piC.k.predict, col = 'goldenrod')
  
  pi_C_atr_kop06 = function(He){
    predict(atr.piC.k.predict, data.frame(conc = He)) / predict(atr.piC.k.predict, data.frame(conc = 0)) 
  } 

#Plot and check fits    
  plot(c(0,20,200), auc.kop / auc.kop.ctrl, pch = 16, ylim=c(0,1), ylab = 'relative auc', col = 'goldenrod',
       xlab = 'atrazine conentration (ppb)')
  
  piC.kop06.df = data.frame(conc = c(0:200),
                            Prediction = 0,
                            Lower = 0,
                            Upper = 0)
  
  piC.kop06.df[,2:4] <- predict(atr.piC.k.predict, newdata = piC.kop06.df, 
                                interval = 'confidence', level = 0.95)
  
  lines(piC.kop06.df$conc, piC.kop06.df$Prediction / piC.kop06.df$Prediction[1], col='goldenrod', lty=2)
  lines(piC.kop06.df$conc, piC.kop06.df$Lower / piC.kop06.df$Prediction[1], col='goldenrod', lty=3)
  lines(piC.kop06.df$conc, piC.kop06.df$Upper / piC.kop06.df$Prediction[1], col='goldenrod', lty=3)
  
  title(main = expression(paste(pi[C], '(atrazine) from Koprivnikar06')))