#Data extraction and model fitting to Bakry 2011 data
require(drc)

#Snail (H. duryi) toxicity ##########
#Malathion #####################
  mun.mal = data.frame(conc = c(.176, .480, .830, 1.760, 3.360),
                       mort = c(0,   .10, .25, .50 , .90),
                       surv = 0)
  mun.mal$surv = 1 - mun.mal$mort
  
    plot(log10(mun.mal$conc[c(2:5)]), qnorm(mun.mal$mort[c(2:5)], mean = 5), pch = 16)
  
  plot(mun.mal$conc*1000, mun.mal$mort, pch = 16, 
       xlab = 'Malathion (ppb)', ylab = 'prop dead', ylim = c(0,1), xlim = c(0,5000),
       main = 'D-R function based on reported values')
      
    segments(x0 = 3120, y0 = 0.5, x1 = 990, y1 = 0.5)

  
#function based on provided parameters
  lc50.mal.mun = 1.76
    se.lc50.mal.mun = mean(log10(3.12/lc50.mal.mun), log10(lc50.mal.mun/0.99)) / 1.96
  slp.mal.mun = 2.68  
  
  fx.mun.mal = function(In, lc = lc50.mal.mun){
    Ins = In/1000
    pnorm(slp.mal.mun * log10(Ins/lc))
  }
  
    lines(c(0:5000), sapply(c(0:5000), fx.mun.mal), lty = 2, col = 2)
    lines(c(0:5000), sapply(c(0:5000), fx.mun.mal, lc = 3.12), lty = 3, col = 2)
    lines(c(0:5000), sapply(c(0:5000), fx.mun.mal, lc = 0.99), lty = 3, col = 2)
    
  muNq_mal_Bakry11_uncertainty = function(In){
    Ins = In/1000
    lc50 = 10^(rnorm(1, log10(lc50.mal.mun), se.lc50.mal.mun))
    pnorm(slp.mal.mun * log10(Ins/lc50))
  }   
  
  points(seq(0,5000,5), sapply(seq(0,5000,5), muNq_mal_Bakry11_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
  
#function based on fit to lc values
  plot(mun.mal$conc, mun.mal$mort, pch = 16, cex = 1.2, 
       xlab = 'Malathion (ppm)', ylab = 'prop dead', ylim = c(0,1), xlim = c(0,5),
       main = 'D-R function based on reported values')
    segments(x0 = 3.120, y0 = 0.5, x1 = .990, y1 = 0.5)
  
  bak11.mod = drm(mort ~ conc, data = mun.mal, type = 'binomial',
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
    
    muNq_mal_Bakry11_uncertainty2<-function(In){
      Ins = In/1000
      rdrm(1, LL.2(), coef(bak11.mod), Ins, yerror = 'rbinom', ypar = 100)$y / 100 #estimate deaths / live 
    }
#Plotting function to check estimates vs data and model  
    points(seq(0,5,0.01), sapply(seq(0, 5000, 10), muNq_mal_Bakry11_uncertainty2), 
           col=4, pch=5, cex=0.5)
    
    
#Deltamethrin  #############    
mun.del = data.frame(conc = c(.482, 1.21, 2.034, 4.82, 7.26),
                     mort = c(0,   .10, .25, .50 , .90),
                     surv = 0)
  mun.del$surv = 1 - mun.del$mort

plot(mun.del$conc*1000, mun.del$mort, pch = 16, cex = 1.2, ylim = c(0,1), xlim = c(0,10000),
     xlab = 'Deltamethrin (ppb)', ylab = 'prop dead', 
     main = 'D-R function based on reported values')
    segments(x0 = 3100, y0 = 0.5, x1 = 7700, y1 = 0.5)

#function based on provided parameters

lc50.del.mun = 4.82
#Pretty sure lower confidence limit is supposed to be 3.1, not 0.31, 
  #don't want to assume typos, but makes sense and function ultimately isn't affected too much
  se.lc50.del.mun = mean(log10(7.7/lc50.del.mun), log10(lc50.del.mun/3.1)) / 1.96
slp.del.mun = 2.74  

fx.mun.del = function(In, lc = lc50.del.mun){
  Ins = In/1000
  pnorm(slp.del.mun * log10(Ins/lc))
}

  lines(seq(0,10000,10), sapply(seq(0,10000,10), fx.mun.del), lty = 2, col = 2)
  lines(seq(0,10000,10), sapply(seq(0,10000,10), fx.mun.del, lc = 7.7), lty = 3, col = 2)
  lines(seq(0,10000,10), sapply(seq(0,10000,10), fx.mun.del, lc = 3.1), lty = 3, col = 2)

muNq_del_Bakry11_uncertainty = function(In){
  Ins = In/1000
  lc50 = 10^(rnorm(1, log10(lc50.del.mun), se.lc50.del.mun))
  pnorm(slp.del.mun * log10(Ins/lc50))
}   

points(seq(0,10000,10), sapply(seq(0,10000,10), muNq_del_Bakry11_uncertainty), 
       pch = 5, col = 4, cex = 0.5)

#Function based on fit to lc values
plot(mun.del$conc, mun.del$mort, pch = 16, cex = 1.2, ylim = c(0,1), xlim = c(0,10),
     xlab = 'Deltamethrin (ppm)', ylab = 'prop dead', 
     main = 'D-R function based on reported values')
  segments(x0 = .310, y0 = 0.5, x1 = 7.700, y1 = 0.5)

bak11.mod2 = drm(mort ~ conc, data = mun.del, type = 'binomial',
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
  
    lines(bak.del.df$conc, bak.del.df$Prediction, col = 2, lty=2)
    lines(bak.del.df$conc, bak.del.df$Lower, col = 2, lty=3)
    lines(bak.del.df$conc, bak.del.df$Upper, col = 2, lty=3)

  muNq_del_Bakry11_uncertainty2<-function(In){
      Ins = In/1000
      rdrm(1, LL.2(), coef(bak11.mod2), Ins, yerror = 'rbinom', ypar = 100)$y / 100 #estimate deaths / live 
    }
#Plotting function to check estimates vs data and model  
points(seq(0,10,0.01), sapply(seq(0, 10000, 10), muNq_del_Bakry11_uncertainty2), 
       col=4, pch=5, cex=0.5)

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
  #So we'll assume there's no reproduction at lc90 and fit 2-parameter logistic function with 1 degree freedom
  #Values for controls are equal to the sum of mean eggs per day produced in the control up to the day where
  #reproduction stops in the lc10 test groups (day4 for malathion, day 14 for deltamethrin)    
bakry11<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/bakry2011.csv')
  
#get mean eggs/snail/day across all time points for each treatment   
  fn.ctrl = mean(bakry11$eggs_snail_day[bakry11$chem == 'control' & bakry11$surv != 0])  
    fn.ctrl.sd = sd(bakry11$eggs_snail_day[bakry11$chem == 'control' & bakry11$surv != 0])
  
  fn.mal = mean(bakry11$eggs_snail_day[bakry11$chem == 'malathion' & bakry11$surv != 0])  
    fn.mal.sd = sd(bakry11$eggs_snail_day[bakry11$chem == 'malathion' & bakry11$surv != 0])
      
  fn.del = mean(bakry11$eggs_snail_day[bakry11$chem == 'deltamethrin' & bakry11$surv != 0])  
    fn.del.sd = sd(bakry11$eggs_snail_day[bakry11$chem == 'deltamethrin' & bakry11$surv != 0])
      
      
#assume lc90 halts all reproduction to provide third data point
  fn.bak = data.frame(mal = c(0, 0.48, 3.36),
                      del = c(0, 1.21, 7.26),
                      del.r = c(fn.ctrl, fn.mal, 0),
                      mal.r = c(fn.ctrl, fn.del, 0))  
      
#Malathion reproduction ##########      
  fn.bak.mal = drm(mal.r ~ mal, data = fn.bak, type = 'continuous',
                   fct = L.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                              fixed = c(NA, 0, fn.bak$mal.r[1], NA)))
  
    summary(fn.bak.mal)
    
    fn.bak.mal.pred = function(In){
      Ins = In/1000
      predict(fn.bak.mal, newdata = data.frame(mal = Ins), interval = 'confidence', level = 0.95)
    }  
    
plot(fn.bak$mal*1000, fn.bak$mal.r / fn.bak$mal.r[1] , ylim = c(0,1), pch = 16,
       xlab = 'malathion (ppm)', ylab = 'relative mean eggs/snail')
    
  lines(seq(0,4000,4), sapply(seq(0,4000,4), fn.bak.mal.pred)[1,] / fn.bak$mal.r[1], col = 2, lty=2)
  lines(seq(0,4000,4), sapply(seq(0,4000,4), fn.bak.mal.pred)[2,] / fn.bak$mal.r[1], col = 2, lty=3)
  lines(seq(0,4000,4), sapply(seq(0,4000,4), fn.bak.mal.pred)[3,] / fn.bak$mal.r[1], col = 2, lty=3)
    
  fN.mal.fx.uncertainty = function(In){
    if(In == 0) fn = 1 else{
      Ins = In/1000
    init = predict(fn.bak.mal, newdata = data.frame(mal = Ins), se.fit = T)
    fn = rnorm(1, init[1], init[2]) / fn.bak$mal.r[1]
    while(fn < 0 && fn > 1.00000){
      fn = rnorm(1, init[1], init[2]) / fn.bak$mal.r[1]
    }
  }
     return(fn)
} #normalized to 1, upper limit at 1, lower limit at 0
  
  points(seq(0,4000,10), sapply(seq(0,4000,10), fN.mal.fx.uncertainty), pch = 5, cex = 0.5, col = 4)

#deltamethrin reproduction ##########      
  fn.bak.del = drm(del.r ~ del, data = fn.bak, type = 'continuous',
                   fct = L.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                              fixed = c(NA, 0, fn.bak$del.r[1], NA)))
    summary(fn.bak.del)
  
  fn.bak.del.pred = function(In){
    Ins = In/1000
    predict(fn.bak.del, newdata = data.frame(del = Ins), interval = 'confidence', level = 0.95)
  }   
  
  plot(fn.bak$del*1000, fn.bak$del.r / fn.bak$del.r[1] , ylim = c(0,1), pch = 16,
       xlab = 'deltamethrin (ppm)', ylab = 'relative mean eggs/snail')
  
    lines(seq(0,7000,7), sapply(seq(0,7000,7), fn.bak.del.pred)[1,] / fn.bak$del.r[1], col = 2, lty=2)
    lines(seq(0,7000,7), sapply(seq(0,7000,7), fn.bak.del.pred)[2,] / fn.bak$del.r[1], col = 2, lty=3)
    lines(seq(0,7000,7), sapply(seq(0,7000,7), fn.bak.del.pred)[3,] / fn.bak$del.r[1], col = 2, lty=3)
  
  fNq_del_Bakry11_uncertainty = function(In){
    if(In == 0) fn = 1 else{
      Ins = In/1000
    init = predict(fn.bak.del, newdata = data.frame(del = Ins), se.fit = T)
    fn = rnorm(1, init[1], init[2]) / fn.bak$del.r[1]
    while(fn < 0 && fn > 1.00000){
      fn = rnorm(1, init[1], init[2]) / fn.bak$del.r[1]
    }
  }
    return(fn)
} 
  
  points(seq(0,7000,10), sapply(seq(0,7000,10), fNq_del_Bakry11_uncertainty), pch = 5, cex = 0.5, col = 4)
  
  keep.bak11.N = c('mun.mal', 'muNq_mal_Bakry11_uncertainty', 'bak11.mod', 'muNq_del_Bakry11_uncertainty', 'bak11.mod2',
                   'fNq_mal_Bakry11_uncertainty', 'fn.bak.mal', 'fNq_del_Bakry11_uncertainty', 'fn.bak.del')