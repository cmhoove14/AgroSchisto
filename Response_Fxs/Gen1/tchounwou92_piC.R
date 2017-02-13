require(drc)

#Cercarial mortality (S. mansoni) from Tchounwou 1992 ##################
cerc = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/tchounwou92.csv')
  cerc.t = subset(cerc, chem == 'malathion')
  cerc.t$dead = round(cerc.t$dead)
time = seq(0,25,0.1)

#Tchounwou Data plotted ############
plot(cerc.t$time_hrs[cerc.t$conc==0], cerc.t$surv[cerc.t$conc==0]/100, pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 2:length(unique(cerc.t$conc[cerc.t$chem == 'malathion']))){
    points(cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[i]], 
           cerc.t$surv[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[i]]/100, pch=16,
           col = i)
  }
legend('topright', title = 'Mal(ppm)', legend = c(0,50,100,150,200,250), pch = c(17,rep(16,4)),
       col = c(1:6), cex=0.7)
#Fit to tchounwou control ########
#Malathion experiment control points
  cerc.ctrl = subset(cerc.t, conc == 0)
  cerc.mal.ctrl = drm(alive/total ~ time_hrs, total, data = cerc.ctrl, 
                      type = 'binomial', fct = LL2.2())
    
    summary(cerc.mal.ctrl)

  cerc.ctrl.df = data.frame(time_hrs = seq(0,24,0.1),
                            Prediction = 0,
                            SE = 0,
                            Lower = 0,
                            Upper = 0)
    
    cerc.ctrl.df[2:3] = predict(cerc.mal.ctrl, newdata = cerc.ctrl.df, 
                                se.fit = TRUE)
      cerc.ctrl.df$Lower = cerc.ctrl.df$Prediction - cerc.ctrl.df$SE
        cerc.ctrl.df$Lower[which(cerc.ctrl.df$Lower<0)] <- 0
      cerc.ctrl.df$Upper = cerc.ctrl.df$Prediction + cerc.ctrl.df$SE
        cerc.ctrl.df$Upper[which(cerc.ctrl.df$Upper>1)] <- 1
      
  lines(cerc.ctrl.df$time_hrs, cerc.ctrl.df$Lower, lty=3)
  lines(cerc.ctrl.df$time_hrs, cerc.ctrl.df$Upper, lty=3)
  
  tch.cerc.mal0.surv = function(t){
    (1/(1+exp(cerc.mal.ctrl$coefficients[1]*(log(t)-cerc.mal.ctrl$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.cerc.mal0.surv(time), lty=2)
  
  auc.tch.mal0=integrate(f = tch.cerc.mal0.surv, lower=0, upper=24)[1]$value

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

#Malathion = 100000 ppb
  cerc.mal100 = drm(cerc.t$mort[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[2]]/100 ~ 
                   cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[2]], 
                  type = 'binomial', fct = LL2.2())
  
  tch.cerc.mal100.surv = function(t){
    1-(1/(1+exp(cerc.mal100$coefficients[1]*(log(t)-cerc.mal100$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.cerc.mal100.surv(time), lty=2, col=3)
  
  auc.tch.mal100=integrate(f = tch.cerc.mal100.surv, lower=0, upper=24)[1]$value

#Malathion = 150000
  cerc.mal150 = drm(cerc.t$mort[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[3]]/100 ~ 
                      cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[3]], 
                    type = 'binomial', fct = LL2.2())
  
  tch.cerc.mal150.surv = function(t){
    1-(1/(1+exp(cerc.mal150$coefficients[1]*(log(t)-cerc.mal150$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.cerc.mal150.surv(time), lty=2, col=4)
  
  auc.tch.mal150=integrate(f = tch.cerc.mal150.surv, lower=0, upper=24)[1]$value

#Malathion = 200000
  cerc.mal200 = drm(cerc.t$mort[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[4]]/100 ~ 
                      cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[4]], 
                    type = 'binomial', fct = LL2.2())
  
  tch.cerc.mal200.surv = function(t){
    1-(1/(1+exp(cerc.mal200$coefficients[1]*(log(t)-cerc.mal200$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.cerc.mal200.surv(time), lty=2, col=5)
  
  auc.tch.mal200=integrate(f = tch.cerc.mal200.surv, lower=0, upper=24)[1]$value

#Malathion = 250000
  cerc.mal250 = drm(cerc.t$mort[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[5]]/100 ~ 
                      cerc.t$time_hrs[cerc.t$conc==unique(cerc.t$conc[cerc.t$chem == 'malathion'])[5]], 
                    type = 'binomial', fct = LL2.2())
  
  tch.cerc.mal250.surv = function(t){
    1-(1/(1+exp(cerc.mal250$coefficients[1]*(log(t)-cerc.mal250$coefficients[2]))))
  } 
  
  lines(x=time, y=tch.cerc.mal250.surv(time), lty=2, col=6)
  
  auc.tch.mal250=integrate(f = tch.cerc.mal250.surv, lower=0, upper=24)[1]$value  

#Derive functional response of pi_C to malathion concentration #############
#Check decrease in AUC across malathion concentration treating AUC as continuous variable
  mal.auc = c(auc.tch.mal0, auc.tch.mal50, auc.tch.mal100, auc.tch.mal150,
              auc.tch.mal200, auc.tch.mal250)

  mal.con = c(0,50,100,150,200,250)
  mal.df = data.frame('conc' = mal.con,
                      'auc24' = mal.auc)
#Fit 3-parameter log-logistic model (constrained with lower limit=0 and fit slope, upper limit, and ED50)
  mal.piC.predict = drm(mal.auc ~ mal.con, data = mal.df, type = 'continuous',
                        fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                                   fixed = c(NA, 0, auc.tch.mal0, NA)))
    summary(mal.piC.predict)
    plot(mal.piC.predict)
  
  pi_C_mal_tch92 = function(In){
    Ins = In/1000   #Model derived in ppm; convert concentration to ppb for response
    predict(mal.piC.predict, data.frame(mal.con = Ins)) / predict(mal.piC.predict, data.frame(mal.con = 0))
  }
  
#Function to normalize results to 1; check fit to data
  piC.mal.tch92.df = data.frame(mal.con = c(0:250),
                                Prediction = 0,
                                Lower = 0,
                                Upper = 0)
  
  piC.mal.tch92.df[,2:4] <- predict(mal.piC.predict, newdata = piC.mal.tch92.df, 
                                    interval = 'confidence', level = 0.95)
    
  plot(mal.df$conc, mal.df$auc24 / mal.df$auc24[1], pch = 16, xlab = 'malathion conc (ppb)', ylab = 'relative auc')
  
    lines(piC.mal.tch92.df$mal.con, piC.mal.tch92.df$Prediction / piC.mal.tch92.df$Prediction[1], col = 2, lty=2)
    lines(piC.mal.tch92.df$mal.con, piC.mal.tch92.df$Lower / piC.mal.tch92.df$Prediction[1], col = 2, lty=3)
    lines(piC.mal.tch92.df$mal.con, piC.mal.tch92.df$Upper / piC.mal.tch92.df$Prediction[1], col = 2, lty=3)
  
  title(main = expression(paste(pi[C], '(malathion) from Tchounwou92 data')))
  