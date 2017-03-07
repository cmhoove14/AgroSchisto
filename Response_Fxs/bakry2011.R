#Data extraction and model fitting to Bakry 2011 data
require(drc)

#Snail (H. duryi) toxicity ##########
#Malathion
  mun.mal = data.frame(conc = c(.176, .480, .830, 1.760, 3.360),
                       mort = c(0,   .10, .25, .50 , .90),
                       surv = 0)
  mun.mal$surv = 1 - mun.mal$mort
  
  se = (log10(3.12) - log10(1.76)) / 1.96 #st. err of lc50 in ppm
  
  mun.mal$se = (mun.mal$conc / 1.760) * se #st. err proportional to concentration

  plot(mun.mal$conc, mun.mal$mort, pch = 16, cex = 1.2, 
       xlab = 'Malathion (ppm)', ylab = 'prop dead', ylim = c(0,1), xlim = c(0,5))
    for(i in 1:length(unique(mun.mal$conc))){
      segments(x0 = mun.mal$conc[i] + mun.mal$se[i], y0 = mun.mal$mort[i],
               x1 = mun.mal$conc[i] - mun.mal$se[i], y1 = mun.mal$mort[i])
    }
  
  bak11.mod = drm(mort ~ conc, data = mun.mal, weights = se^-1, type = 'binomial',
                  fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                  fixed = c(NA, 0, 1, NA)))
#    summary(bak11.mod)
    
  muNq_mal_Bakry11<-function(In){
    Ins = In/1000
    predict(bak11.mod, data.frame(conc = Ins))
  }  
  
  bak.mal.df = data.frame(conc = seq(0,5,0.001),
                           Prediction = 0,
                           Lower = 0,
                           Upper = 0)
  
  bak.mal.df[,2:4] <- predict(bak11.mod, newdata = bak.mal.df, 
                               interval = 'confidence', level = 0.95)

    lines(bak.mal.df$conc, bak.mal.df$Prediction, col = 2, lty=2)
    lines(bak.mal.df$conc, bak.mal.df$Lower, col = 2, lty=3)
    lines(bak.mal.df$conc, bak.mal.df$Upper, col = 2, lty=3)
    
    muNq_mal_Bakry11_uncertainty<-function(In){
      Ins = In/1000
      rdrm(1, LL.2(), coef(bak11.mod), Ins, yerror = 'rbinom', ypar = 100)$y / 100 #estimate deaths / live 
    }
#Plotting function to check estimates vs data and model  
    points(seq(0,5,0.01), sapply(seq(0, 5000, 10), muNq_mal_Bakry11_uncertainty), 
           col=2, pch=1, cex=0.5)
    
    
#Deltamethrin     
mun.del = data.frame(conc = c(.482, 1.21, 2.034, 4.82, 7.26),
                     mort = c(0,   .10, .25, .50 , .90),
                     surv = 0)
  mun.del$surv = 1 - mun.del$mort

se2 = (log10(7.7) - log10(4.82)) / 1.96 #st. err of lc50 in ppm

  mun.del$se2 = (mun.del$conc / 4.82) * se2 #st. err proportional to concentration

plot(mun.del$conc, mun.del$mort, pch = 16, cex = 1.2, 
     xlab = 'Deltamethrin (ppm)', ylab = 'prop dead', ylim = c(0,1), xlim = c(0,10))
  for(i in 1:length(unique(mun.del$conc))){
    segments(x0 = mun.del$conc[i] + mun.del$se2[i], y0 = mun.del$mort[i],
             x1 = mun.del$conc[i] - mun.del$se2[i], y1 = mun.del$mort[i])
  }

bak11.mod2 = drm(mort ~ conc, data = mun.del, weights = se2^-1, type = 'binomial',
                fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                           fixed = c(NA, 0, 1, NA)))
  summary(bak11.mod2)

muNq_del_Bakry11<-function(In){
  Ins = In/1000
  predict(bak11.mod2, data.frame(conc = Ins))
}  

bak.del.df = data.frame(conc = seq(0,10,0.01),
                        Prediction = 0,
                        Lower = 0,
                        Upper = 0)

  bak.del.df[,2:4] <- predict(bak11.mod2, newdata = bak.del.df, 
                              interval = 'confidence', level = 0.95)
  
    lines(bak.del.df$conc, bak.del.df$Prediction, col = 4, lty=2)
    lines(bak.del.df$conc, bak.del.df$Lower, col = 4, lty=3)
    lines(bak.del.df$conc, bak.del.df$Upper, col = 4, lty=3)

  muNq_del_Bakry11_uncertainty<-function(In){
      Ins = In/1000
      rdrm(1, LL.2(), coef(bak11.mod2), Ins, yerror = 'rbinom', ypar = 100)$y / 100 #estimate deaths / live 
    }
#Plotting function to check estimates vs data and model  
points(seq(0,10,0.01), sapply(seq(0, 10000, 10), muNq_del_Bakry11_uncertainty), 
       col=4, pch=1, cex=0.5)

#The paper also provides info on longitudinal survival of snail cohorts exposed to LC10 of each insecticide ################
  #So let's compare that data with the expected long. survival from the model at the same concentration
  bakry11<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/bakry2011.csv')

#Plot data from paper  
  plot(bakry11$time[bakry11$chem == 'control'], bakry11$surv[bakry11$chem == 'control'], ylim = c(0,50), pch = 16, cex = 1.2,
       xlab = 'time (days)', ylab = 'n alive')
    points(bakry11$time[bakry11$chem == 'malathion'], bakry11$surv[bakry11$chem == 'malathion'], col = 2, pch = 16, cex = 1.2)
    points(bakry11$time[bakry11$chem == 'deltamethrin'], bakry11$surv[bakry11$chem == 'deltamethrin'], col = 6, pch = 16, cex = 1.2)
    title('long. snail survival exposed to LC10s')
    legend('bottomleft', legend = c('control', 'deltamethrin', 'malathion'), col = c(1,2,6), pch = 16, cex=0.7, 
           title = 'observed')
    
#Generate model predictions    
  #Background mortality seems to be pretty constant (constant slope in control group), so let's use control group to check that value
    m = (bakry11$surv[bakry11$chem == 'control'][7] - bakry11$surv[bakry11$chem == 'control'][1]) /
        (bakry11$time[bakry11$chem == 'control'][7] - bakry11$time[bakry11$chem == 'control'][1])

    m.df = data.frame('days' = c(0:28),
                      'control' = 0,
                      'mal' = 0,
                      'delt' = 0)
    
    m.df[1,c(2:4)] = 50
    
    for(i in 1:28){
      m.df[i+1, 2] = m.df[i,2] + m  #subtract deaths per day (from control slope derived above)
      m.df[i+1, 3] = m.df[i,3] + m - m.df[i,3] * muNq_mal_Bakry11(480)  #for agrochemical toxicity, additional mortality as per capita deaths
      m.df[i+1, 4] = m.df[i,4] + m - m.df[i,4] * muNq_del_Bakry11(1210) #for agrochemical toxicity, additional mortality as per capita deaths
    }
    
    lines(m.df$days, m.df$control, lty = 2, lwd=2)
    lines(m.df$days, m.df$mal, lty = 2, col=2, lwd=2)
    lines(m.df$days, m.df$delt, lty = 2, col=6, lwd=2)
      legend('left', legend = c('control', 'deltamethrin', 'malathion'), col = c(1,2,6),  lty=2, cex=0.7, 
             title = 'modeled')
      
#Now lets look at reproduction over time ###########
  #The paper only provides mean eggs / snail over 4 weeks exposure to lc10 of each insecticide
  #This isn't enough data to generate d-r function without some assumptions
  #So we'll assume there's no reproduction at lc90 and fit 2-parameter ll function with 1 degree freedom
  #Values for controls are equal to the sum of mean eggs per day produced in the control up to the day where
  #reproduction stops in the lc10 test groups (day4 for malathion, day 14 for deltamethrin)    
      
  fn.bak = data.frame(mal = c(0, .480, 3.360), #0 ppm = control experiment; assume 0 reproduction at lc90
                      del = c(0, 1.21, 7.26),  #0 ppm = control experiment; assume 0 reproduction at lc90
                      mal.r = c(5.2 + 11.4, 3.9, 0),  
                      del.r = c(5.2 + 11.4 + 10 + 8.8, 14.33, 0))  
      
  fn.bak.mal = drm(mal.r ~ mal, data = fn.bak, type = 'continuous',
                   fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                              fixed = c(NA, 0, fn.bak$mal.r[1], NA)))
  
    summary(fn.bak.mal)
    
  fn.bak.del = drm(del.r ~ del, data = fn.bak, type = 'continuous',
                   fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                              fixed = c(NA, 0, fn.bak$del.r[1], NA)))
    summary(fn.bak.del)
      
  plot(fn.bak$mal, fn.bak$mal.r / fn.bak$mal.r[1] , ylim = c(0,1), xlim = c(0,8), pch = 16, col = 'red',
       xlab = 'conc (ppm)', ylab = 'relative mean eggs/snail')
    points(fn.bak$del, fn.bak$del.r / fn.bak$del.r[1], pch = 16, col = 'purple')

    bak.mal.df2 = data.frame(mal = seq(0,10,0.001),
                             Prediction = 0,
                             Lower = 0,
                             Upper = 0)
    
    bak.mal.df2[,2:4] <- predict(fn.bak.mal, newdata = bak.mal.df2, 
                                interval = 'confidence', level = 0.95)
    
      lines(bak.mal.df2$mal, bak.mal.df2$Prediction / fn.bak$mal.r[1], col = 2, lty=2)
      lines(bak.mal.df2$mal, bak.mal.df2$Lower / fn.bak$mal.r[1], col = 2, lty=3)
      lines(bak.mal.df2$mal, bak.mal.df2$Upper / fn.bak$mal.r[1], col = 2, lty=3)
      
    bak.del.df2 = data.frame(del = seq(0,10,0.001),
                             Prediction = 0,
                             Lower = 0,
                             Upper = 0)
    
    bak.del.df2[,2:4] <- predict(fn.bak.del, newdata = bak.del.df2, 
                                 interval = 'confidence', level = 0.95)
    
    lines(bak.del.df2$del, bak.del.df2$Prediction / fn.bak$del.r[1], col = 'purple', lty=2)
    lines(bak.del.df2$del, bak.del.df2$Lower / fn.bak$del.r[1], col = 'purple', lty=3)
    lines(bak.del.df2$del, bak.del.df2$Upper / fn.bak$del.r[1], col = 'purple', lty=3)
    
  fNq_mal_Bakry11_uncertainty = function(In){
    Ins = In/1000 #Input in parts per billion, model in ppm
    c1 = predict(fn.bak.mal, newdata = data.frame(mal = 0))
    
    ts = predict(fn.bak.mal, newdata = data.frame(mal = Ins), se.fit = TRUE)
    
    fN = rnorm(1,ts[1],ts[2]) / c1
    
    fN
  } 
#Plotting function to check estimates vs data and model  
  #for(i in c(0:10000)){
  #  points(i/1000, fNq_mal_Bakry11_uncertainty(i), pch = 17, cex = 0.4, col = 4)
  #}
  
  fNq_del_Bakry11_uncertainty = function(In){
    Ins = In/1000 #Input in parts per billion, model in ppm
    c1 = predict(fn.bak.del, newdata = data.frame(del = 0))
    
    ts = predict(fn.bak.del, newdata = data.frame(del = Ins), se.fit = TRUE)

    fN = rnorm(1,ts[1],ts[2]) / c1
    
    fN
  } 
#Plotting function to check estimates vs data and model  
  #for(i in c(0:10000)){
  #  points(i/1000, fNq_del_Bakry11_uncertainty(i), pch = 17, cex = 0.4, col = 3)
  #}
  
  keep.bak11.N = c('mun.mal', 'muNq_mal_Bakry11_uncertainty', 'bak11.mod', 'muNq_del_Bakry11_uncertainty', 'bak11.mod2',
                   'fNq_mal_Bakry11_uncertainty', 'fn.bak.mal', 'fNq_del_Bakry11_uncertainty', 'fn.bak.del')