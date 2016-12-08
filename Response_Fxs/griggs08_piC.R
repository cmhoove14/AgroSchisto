require(drc)

cerc = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/all_cercarial_mortality.csv')
  cerc.g = subset(cerc, Study == 'griggs2008')
time = c(0:25)

#Plot control survival curve, estimate time-dependent die-off function, and get auc from 0-24 hours ########
plot(x = cerc.g$time_hrs[cerc.g$chem == 'control'], 
     y = cerc.g$surv[cerc.g$chem == 'control']/100,
     xlab = 'time (hrs)', ylab = 'prop alive', pch = 16, xlim = c(0,25), ylim = c(0,1))

  grg.ctrl = drm(cerc.g$mort[cerc.g$chem == 'control']/100 ~ 
                   cerc.g$time_hrs[cerc.g$chem == 'control'], 
                  type = 'binomial', fct = LL2.2())
  
  grg.cerc.ctrl.surv = function(t){
    1-(1/(1+exp(grg.ctrl$coefficients[1]*(log(t)-grg.ctrl$coefficients[2]))))
  } 
  
  lines(x=time, y=grg.cerc.ctrl.surv(time), lty=2)
  
  auc.grg.ctrl=integrate(f = grg.cerc.ctrl.surv, lower=0, upper=24)[1]$value
  
#Function for low dose treatment (atrazine 15 ppb & metolachlor 10 ppb combined)  ##########
  points(cerc.g$time_hrs[cerc.g$conc == 15],
         cerc.g$surv[cerc.g$conc == 15]/100, pch = 17)
  
  grg.low = drm(cerc$mort[cerc$conc == 15 & cerc$Study == 'griggs2008']/100 ~ 
                   cerc$time_hrs[cerc$conc == 15 & cerc$Study == 'griggs2008'], 
                 type = 'binomial', fct = LL2.2())
  
  grg.cerc.low = function(t){
    1-(1/(1+exp(grg.low$coefficients[1]*(log(t)-grg.low$coefficients[2]))))
  } 
  
  lines(x=time, y=grg.cerc.low(time), lty=2)
  
  auc.grg.atr15=integrate(f = grg.cerc.low, lower=0, upper=24)[1]$value
  
#Function for low dose treatment (atrazine 100 ppb & metolachlor 85 ppb combined)   ###############
  points(cerc.g$time_hrs[cerc.g$conc == 100],
         cerc.g$surv[cerc.g$conc == 100]/100, pch = 15)
  
  grg.hi = drm(cerc$mort[cerc$conc == 100 & cerc$Study == 'griggs2008']/100 ~ 
                  cerc$time_hrs[cerc$conc == 100 & cerc$Study == 'griggs2008'], 
                type = 'binomial', fct = LL2.2())
  
  grg.cerc.hi = function(t){
    1-(1/(1+exp(grg.hi$coefficients[1]*(log(t)-grg.hi$coefficients[2]))))
  } 
  
  lines(x=time, y=grg.cerc.hi(time), lty=3)
  
  auc.grg.atr100=integrate(f = grg.cerc.hi, lower=0, upper=24)[1]$value  
#Summarize ########
  title(main='Griggs2008 Cercarial mortality (E.trivolvis)')
  legend('bottomleft', legend = c('control', '15ppb', '100ppb'), pch = c(16,17,15), cex=0.8)
  
  rel.auc.grg = c(auc.grg.ctrl,auc.grg.atr15, auc.grg.atr100)/auc.grg.ctrl
  
  plot(c(0,15,100), rel.auc.grg, pch = 16, ylim=c(0,1), ylab = 'relative auc',
       xlab = 'atrazine conentration (ppb)')
  
  grg.df = data.frame('auc24' = rel.auc.grg,
                      'conc' = c(0,15,100))
  
  atr.piC.g.predict = nls(auc24 ~ exp(-b*(conc+1e-6)), data = grg.df,
                          start = list(b=0.1))
  summary(atr.piC.g.predict)
  
  atr.g.con = c(0:201)
  
  pi_C_atr_grg08 = function(q){
    exp(-summary(atr.piC.g.predict)$parameters[1]*(q)) 
  } #Model derived in ppm; convert concentration to ppb for response
  
  lines(atr.g.con, pi_C_atr_grg08(atr.g.con), lty=2, col='red')
  title(main = expression(paste(pi[C], '(atrazine) from Griggs08')))
  legend('bottomright', cex = 0.7,
         legend = paste('b= ', round(summary(atr.piC.g.predict)$parameters[1],3)))