#Data extraction and model fitting to Tantawy 2002 data
require(drc)

#Snail toxicity ##########
#Butachlor
  lc50.n.but<-6.5
  slp.n.but<-2.1
  
    mu_Nq_but_tant02<-function(In){
      Ins = In/1000
      1 - 1/(1+exp(slp.n.but*(log(Ins)-log(lc50.n.but))))
    } 
    
    mu_Nq_but_tant02(lc50.n.but)
    
    plot(c(0:50000), mu_Nq_but_tant02(c(0:50000)), lwd = 2, type = 'l', xlab = 'butachlor concentration (ppb)',
         ylab = 'mu_Nq', ylim = c(0,1), main = 'butachlor toxicity to snails, tantawy2002')

#Fluazifop-p-butyl      
  lc50.n.fpb<-17.6
  slp.n.fpb<-1.8
    
  mu_Nq_fpb_tant02<-function(In){
    Ins = In/1000
    1 - 1/(1+exp(slp.n.fpb*(log(Ins)-log(lc50.n.fpb))))
  } 
    
  mu_Nq_fpb_tant02(lc50.n.fpb)
    
  plot(c(0:60000), mu_Nq_fpb_tant02(c(0:60000)), lwd = 2, type = 'l', xlab = 'fluazifop-p-butyl concentration (ppb)',
       ylab = 'mu_Nq', ylim = c(0,1), main = 'fluazifop-p-butyl toxicity to snails, tantawy2002')
    
#cercarial toxicity ############
  tant<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/Tantawy2002.csv')
  cerc<-subset(tant, larv == 'cercariae')
  
#butachlor toxicity to cercariae  ###############
  plot(cerc$time_hrs[cerc$conc==0 & cerc$chem == 'butachlor'], 
       cerc$surv[cerc$conc==0 & cerc$chem == 'butachlor']/100, 
       pch=17, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 2:length(unique(cerc$conc[cerc$chem == 'butachlor']))){
    points(cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[i] & cerc$chem == 'butachlor'], 
           cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[i] & cerc$chem == 'butachlor']/100, pch=16,
           col = i)
  }
  
#fit to control points
  tant.cerc.but0<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[1] & cerc$chem == 'butachlor']/100 ~ 
                      cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[1] & cerc$chem == 'butachlor'],
                      type = 'binomial', fct = LL2.2())
  
  tant.cerc.but0.surv = function(t){
    (1/(1+exp(tant.cerc.but0$coefficients[1]*(log(t)-tant.cerc.but0$coefficients[2]))))
  } 
  
    lines(x=c(0:24), y=tant.cerc.but0.surv(c(0:24)))
    
    auc.cerc.but0=integrate(f = tant.cerc.but0.surv, lower=0, upper=24)[1]$value  
    
    
#650 ppb
    tant.cerc.but650<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[2] & cerc$chem == 'butachlor']/100 ~ 
                          cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[2] & cerc$chem == 'butachlor'],
                        type = 'binomial', fct = LL2.2())
    
    tant.cerc.but650.surv = function(t){
      (1/(1+exp(tant.cerc.but650$coefficients[1]*(log(t)-tant.cerc.but650$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.but650.surv(c(0:24)), lty=2, col=2)
    
    auc.cerc.but650=integrate(f = tant.cerc.but650.surv, lower=0, upper=24)[1]$value  
    
    
#1500 ppb
    tant.cerc.but1500<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[3] & cerc$chem == 'butachlor']/100 ~ 
                            cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[3] & cerc$chem == 'butachlor'],
                          type = 'binomial', fct = LL2.2())
    
    tant.cerc.but1500.surv = function(t){
      (1/(1+exp(tant.cerc.but1500$coefficients[1]*(log(t)-tant.cerc.but1500$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.but1500.surv(c(0:24)), lty=2, col=3)
    
    auc.cerc.but1500=integrate(f = tant.cerc.but1500.surv, lower=0, upper=24)[1]$value  
    
#4500 ppb
    tant.cerc.but4500<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[4] & cerc$chem == 'butachlor']/100 ~ 
                             cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[4] & cerc$chem == 'butachlor'],
                           type = 'binomial', fct = LL2.2())
    
    tant.cerc.but4500.surv = function(t){
      (1/(1+exp(tant.cerc.but4500$coefficients[1]*(log(t)-tant.cerc.but4500$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.but4500.surv(c(0:24)), lty=2, col=4)
    
    auc.cerc.but4500=integrate(f = tant.cerc.but4500.surv, lower=0, upper=24)[1]$value  
    
#6500 ppb
    tant.cerc.but6500<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[5] & cerc$chem == 'butachlor']/100 ~ 
                             cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'butachlor'])[5] & cerc$chem == 'butachlor'],
                           type = 'binomial', fct = LL2.2())
    
    tant.cerc.but6500.surv = function(t){
      (1/(1+exp(tant.cerc.but6500$coefficients[1]*(log(t)-tant.cerc.but6500$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.but6500.surv(c(0:24)), lty=2, col=5)
    title('butachlor toxicity to cercariae')
    legend('topright', legend = c('control', 650,1500,4500,6500), pch = c(17,16,16,16,16), col = c(1:5), cex=0.8)
    
    auc.cerc.but6500=integrate(f = tant.cerc.but6500.surv, lower=0, upper=24)[1]$value  
  
    
#fluazifop-p-butyl toxicity to cercariae ###########
    plot(cerc$time_hrs[cerc$conc==0 & cerc$chem == 'fluazifop-p-butyl'], 
         cerc$surv[cerc$conc==0 & cerc$chem == 'fluazifop-p-butyl']/100, 
         pch=17, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    for(i in 2:length(unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl']))){
      points(cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[i] & cerc$chem == 'fluazifop-p-butyl'], 
             cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[i] & cerc$chem == 'fluazifop-p-butyl']/100, pch=16,
             col = i)
    }
    
    #fit to control points
    tant.cerc.fpb0<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[1] & cerc$chem == 'fluazifop-p-butyl']/100 ~ 
                          cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[1] & cerc$chem == 'fluazifop-p-butyl'],
                        type = 'binomial', fct = LL2.2())
    
    tant.cerc.fpb0.surv = function(t){
      (1/(1+exp(tant.cerc.fpb0$coefficients[1]*(log(t)-tant.cerc.fpb0$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.fpb0.surv(c(0:24)))
    
    auc.cerc.fpb0=integrate(f = tant.cerc.fpb0.surv, lower=0, upper=24)[1]$value  
    
    
    #1760 ppb
    tant.cerc.fpb1760<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[2] & cerc$chem == 'fluazifop-p-butyl']/100 ~ 
                            cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[2] & cerc$chem == 'fluazifop-p-butyl'],
                          type = 'binomial', fct = LL2.2())
    
    tant.cerc.fpb1760.surv = function(t){
      (1/(1+exp(tant.cerc.fpb1760$coefficients[1]*(log(t)-tant.cerc.fpb1760$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.fpb1760.surv(c(0:24)), lty=2, col=2)
    
    auc.cerc.fpb1760=integrate(f = tant.cerc.fpb1760.surv, lower=0, upper=24)[1]$value  
    
    
    #4500 ppb
    tant.cerc.fpb4500<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[3] & cerc$chem == 'fluazifop-p-butyl']/100 ~ 
                             cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[3] & cerc$chem == 'fluazifop-p-butyl'],
                           type = 'binomial', fct = LL2.2())
    
    tant.cerc.fpb4500.surv = function(t){
      (1/(1+exp(tant.cerc.fpb4500$coefficients[1]*(log(t)-tant.cerc.fpb4500$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.fpb4500.surv(c(0:24)), lty=2, col=3)
    
    auc.cerc.fpb4500=integrate(f = tant.cerc.fpb4500.surv, lower=0, upper=24)[1]$value  
    
    #9000 ppb
    tant.cerc.fpb9000<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[4] & cerc$chem == 'fluazifop-p-butyl']/100 ~ 
                             cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[4] & cerc$chem == 'fluazifop-p-butyl'],
                           type = 'binomial', fct = LL2.2())
    
    tant.cerc.fpb9000.surv = function(t){
      (1/(1+exp(tant.cerc.fpb9000$coefficients[1]*(log(t)-tant.cerc.fpb9000$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.fpb9000.surv(c(0:24)), lty=2, col=4)
    
    auc.cerc.fpb9000=integrate(f = tant.cerc.fpb9000.surv, lower=0, upper=24)[1]$value  
    
    #17600 ppb
    tant.cerc.fpb17600<-drm(cerc$surv[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[5] & cerc$chem == 'fluazifop-p-butyl']/100 ~ 
                             cerc$time_hrs[cerc$conc==unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])[5] & cerc$chem == 'fluazifop-p-butyl'],
                           type = 'binomial', fct = LL2.2())
    
    tant.cerc.fpb17600.surv = function(t){
      (1/(1+exp(tant.cerc.fpb17600$coefficients[1]*(log(t)-tant.cerc.fpb17600$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.cerc.fpb17600.surv(c(0:24)), lty=2, col=5)
    title('fluazifop-p-butyl toxicity to cercariae')
    legend('topright', legend = c('control', 1760,4500,9000,17600), pch = c(17,16,16,16,16), col = c(1:5), cex=0.8)
    
    auc.cerc.fpb17600=integrate(f = tant.cerc.fpb17600.surv, lower=0, upper=24)[1]$value  
    
#Fit functional response to auc data across concentration of both herbicides ###############
  but.conc = unique(cerc$conc[cerc$chem == 'butachlor'])
  fpb.conc = unique(cerc$conc[cerc$chem == 'fluazifop-p-butyl'])
  but.auc.cerc = c(auc.cerc.but0, auc.cerc.but650, auc.cerc.but1500, auc.cerc.but4500, auc.cerc.but6500) / auc.cerc.but0
  fpb.auc.cerc = c(auc.cerc.fpb0, auc.cerc.fpb1760, auc.cerc.fpb4500, auc.cerc.fpb9000, auc.cerc.fpb17600) / auc.cerc.fpb0
  
#Plot and function for butachlor  
  plot(but.conc, but.auc.cerc, pch = 16, xlab = 'butachlor concentration (ppb)', ylab = 'relative 24-hr auc', ylim = c(0,1))
  
  tant.but.piC.predict = nls(but.auc.cerc ~ exp(-b*(but.conc)), start = list(b=0.0001))
    summary(tant.but.piC.predict)
  
  pi_C_but_tant02 = function(In){
    exp(-summary(tant.but.piC.predict)$parameters[1]*In) 
  }
  
  lines(c(0:7000), pi_C_but_tant02(c(0:7000)), lty=2, col='red')
    title(main = expression(paste(pi[C], '(butachlor) from Tantawy02 data')))
    
#Plot and function for fluazifop-p-butyl
  plot(fpb.conc, fpb.auc.cerc, pch = 16, xlab = 'fluazifop-p-butyl concentration (ppb)', ylab = 'relative 24-hr auc', ylim = c(0,1))
    
  tant.fpb.piC.predict = nls(fpb.auc.cerc ~ exp(-b*(fpb.conc)), start = list(b=0.0001))
    summary(tant.fpb.piC.predict)
    
  pi_C_fpb_tant02 = function(In){
    exp(-summary(tant.fpb.piC.predict)$parameters[1]*In) 
  }
  
  lines(c(0:20000), pi_C_fpb_tant02(c(0:20000)), lty=2, col='red')
    title(main = expression(paste(pi[C], '(fluazifop-p-butyl) from Tantawy02 data')))
    
#miracidia toxicity ############
  mir<-subset(tant, larv == 'miracidia')
    
#butachlor toxicity to miracidia  ###############
    plot(mir$time_hrs[mir$conc==0 & mir$chem == 'butachlor'], 
         mir$surv[mir$conc==0 & mir$chem == 'butachlor']/100, 
         pch=17, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    for(i in 2:length(unique(mir$conc[mir$chem == 'butachlor']))){
      points(mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[i] & mir$chem == 'butachlor'], 
             mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[i] & mir$chem == 'butachlor']/100, pch=16,
             col = i)
    }
    
#fit to control points
    tant.mir.but0<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[1] & mir$chem == 'butachlor']/100 ~ 
                          mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[1] & mir$chem == 'butachlor'],
                        type = 'binomial', fct = LL2.2())
    
    tant.mir.but0.surv = function(t){
      (1/(1+exp(tant.mir.but0$coefficients[1]*(log(t)-tant.mir.but0$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.mir.but0.surv(c(0:24)))
    
    auc.mir.but0=integrate(f = tant.mir.but0.surv, lower=0, upper=24)[1]$value  
    
    
#650 ppb
    tant.mir.but650<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[2] & mir$chem == 'butachlor']/100 ~ 
                            mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[2] & mir$chem == 'butachlor'],
                          type = 'binomial', fct = LL2.2())
    
    tant.mir.but650.surv = function(t){
      (1/(1+exp(tant.mir.but650$coefficients[1]*(log(t)-tant.mir.but650$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.mir.but650.surv(c(0:24)), lty=2, col=2)
    
    auc.mir.but650=integrate(f = tant.mir.but650.surv, lower=0, upper=24)[1]$value  
    
    
#1500 ppb
    tant.mir.but1500<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[3] & mir$chem == 'butachlor']/100 ~ 
                             mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[3] & mir$chem == 'butachlor'],
                           type = 'binomial', fct = LL2.2())
    
    tant.mir.but1500.surv = function(t){
      (1/(1+exp(tant.mir.but1500$coefficients[1]*(log(t)-tant.mir.but1500$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.mir.but1500.surv(c(0:24)), lty=2, col=3)
    
    auc.mir.but1500=integrate(f = tant.mir.but1500.surv, lower=0, upper=24)[1]$value  
    
#4500 ppb
    tant.mir.but4500<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[4] & mir$chem == 'butachlor']/100 ~ 
                             mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[4] & mir$chem == 'butachlor'],
                           type = 'binomial', fct = LL2.2())
    
    tant.mir.but4500.surv = function(t){
      (1/(1+exp(tant.mir.but4500$coefficients[1]*(log(t)-tant.mir.but4500$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.mir.but4500.surv(c(0:24)), lty=2, col=4)
    
    auc.mir.but4500=integrate(f = tant.mir.but4500.surv, lower=0, upper=24)[1]$value  
    
#6500 ppb
    tant.mir.but6500<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[5] & mir$chem == 'butachlor']/100 ~ 
                             mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'butachlor'])[5] & mir$chem == 'butachlor'],
                           type = 'binomial', fct = LL2.2())
    
    tant.mir.but6500.surv = function(t){
      (1/(1+exp(tant.mir.but6500$coefficients[1]*(log(t)-tant.mir.but6500$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.mir.but6500.surv(c(0:24)), lty=2, col=5)
    title('butachlor toxicity to miracidia')
    legend('topright', legend = c('control', 650,1500,4500,6500), pch = c(17,16,16,16,16), col = c(1:5), cex=0.8)
    
    auc.mir.but6500=integrate(f = tant.mir.but6500.surv, lower=0, upper=24)[1]$value     
    
#fluazifop-p-butyl toxicity to cercariae ###########
    plot(mir$time_hrs[mir$conc==0 & mir$chem == 'fluazifop-p-butyl'], 
         mir$surv[mir$conc==0 & mir$chem == 'fluazifop-p-butyl']/100, 
         pch=17, xlab = 'time(hrs', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    for(i in 2:length(unique(mir$conc[mir$chem == 'fluazifop-p-butyl']))){
      points(mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[i] & mir$chem == 'fluazifop-p-butyl'], 
             mir$surv[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[i] & mir$chem == 'fluazifop-p-butyl']/100, pch=16,
             col = i)
    }
    
#fit to control points
    tant.mir.fpb0<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[1] & mir$chem == 'fluazifop-p-butyl']/100 ~ 
                          mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[1] & mir$chem == 'fluazifop-p-butyl'],
                        type = 'binomial', fct = LL2.2())
    
    tant.mir.fpb0.surv = function(t){
      (1/(1+exp(tant.mir.fpb0$coefficients[1]*(log(t)-tant.mir.fpb0$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.mir.fpb0.surv(c(0:24)))
    
    auc.mir.fpb0=integrate(f = tant.mir.fpb0.surv, lower=0, upper=24)[1]$value  
    
    
#1760 ppb
    tant.mir.fpb1760<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[2] & mir$chem == 'fluazifop-p-butyl']/100 ~ 
                             mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[2] & mir$chem == 'fluazifop-p-butyl'],
                           type = 'binomial', fct = LL2.2())
    
    tant.mir.fpb1760.surv = function(t){
      (1/(1+exp(tant.mir.fpb1760$coefficients[1]*(log(t)-tant.mir.fpb1760$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.mir.fpb1760.surv(c(0:24)), lty=2, col=2)
    
    auc.mir.fpb1760=integrate(f = tant.mir.fpb1760.surv, lower=0, upper=24)[1]$value  
    
    
#4500 ppb
    tant.mir.fpb4500<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[3] & mir$chem == 'fluazifop-p-butyl']/100 ~ 
                             mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[3] & mir$chem == 'fluazifop-p-butyl'],
                           type = 'binomial', fct = LL2.2())
    
    tant.mir.fpb4500.surv = function(t){
      (1/(1+exp(tant.mir.fpb4500$coefficients[1]*(log(t)-tant.mir.fpb4500$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.mir.fpb4500.surv(c(0:24)), lty=2, col=3)
    
    auc.mir.fpb4500=integrate(f = tant.mir.fpb4500.surv, lower=0, upper=24)[1]$value  
    
#9000 ppb
    tant.mir.fpb9000<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[4] & mir$chem == 'fluazifop-p-butyl']/100 ~ 
                             mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[4] & mir$chem == 'fluazifop-p-butyl'],
                           type = 'binomial', fct = LL2.2())
    
    tant.mir.fpb9000.surv = function(t){
      (1/(1+exp(tant.mir.fpb9000$coefficients[1]*(log(t)-tant.mir.fpb9000$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.mir.fpb9000.surv(c(0:24)), lty=2, col=4)
    
    auc.mir.fpb9000=integrate(f = tant.mir.fpb9000.surv, lower=0, upper=24)[1]$value  
    
#17600 ppb
    tant.mir.fpb17600<-drm(mir$surv[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[5] & mir$chem == 'fluazifop-p-butyl']/100 ~ 
                              mir$time_hrs[mir$conc==unique(mir$conc[mir$chem == 'fluazifop-p-butyl'])[5] & mir$chem == 'fluazifop-p-butyl'],
                            type = 'binomial', fct = LL2.2())
    
    tant.mir.fpb17600.surv = function(t){
      (1/(1+exp(tant.mir.fpb17600$coefficients[1]*(log(t)-tant.mir.fpb17600$coefficients[2]))))
    } 
    
    lines(x=c(0:24), y=tant.mir.fpb17600.surv(c(0:24)), lty=2, col=5)
    title('fluazifop-p-butyl toxicity to miracidia')
    legend('topright', legend = c('control', 1760,4500,9000,17600), pch = c(17,16,16,16,16), col = c(1:5), cex=0.8)
    
    auc.mir.fpb17600=integrate(f = tant.mir.fpb17600.surv, lower=0, upper=24)[1]$value  
    
    
#Fit functional response to auc data across concentration of both herbicides ###############
    but.auc.mir = c(auc.mir.but0, auc.mir.but650, auc.mir.but1500, auc.mir.but4500, auc.mir.but6500) / auc.mir.but0
    fpb.auc.mir = c(auc.mir.fpb0, auc.mir.fpb1760, auc.mir.fpb4500, auc.mir.fpb9000, auc.mir.fpb17600) / auc.mir.fpb0
    
  #Plot and function for butachlor  
    plot(but.conc, but.auc.mir, pch = 16, xlab = 'butachlor concentration (ppb)', ylab = 'relative 24-hr auc', ylim = c(0,1))
    
    tant.but.piM.predict = nls(but.auc.mir ~ exp(-b*(but.conc)), start = list(b=0.0001))
      summary(tant.but.piM.predict)
    
    pi_M_but_tant02 = function(In){
      exp(-summary(tant.but.piM.predict)$parameters[1]*In) 
    }
    
    lines(c(0:7000), pi_M_but_tant02(c(0:7000)), lty=2, col='red')
    title(main = expression(paste(pi[M], '(butachlor) from Tantawy02 data')))
    
  #Plot and function for fluazifop-p-butyl
    plot(fpb.conc, fpb.auc.mir, pch = 16, xlab = 'fluazifop-p-butyl concentration (ppb)', ylab = 'relative 24-hr auc', ylim = c(0,1))
    
    tant.fpb.piM.predict = nls(fpb.auc.mir ~ exp(-b*(fpb.conc)), start = list(b=0.0001))
      summary(tant.fpb.piM.predict)
    
    pi_M_fpb_tant02 = function(In){
      exp(-summary(tant.fpb.piM.predict)$parameters[1]*In) 
    }
    
    lines(c(0:20000), pi_M_fpb_tant02(c(0:20000)), lty=2, col='red')
    title(main = expression(paste(pi[M], '(fluazifop-p-butyl) from Tantawy02 data')))
    