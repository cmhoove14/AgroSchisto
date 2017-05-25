#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Data extraction and model fitting to Bakry 2011 data
require(drc)

#Snail (H. duryi) toxicity ##########
#Malathion #####################
  mun.mal = data.frame(conc = c(.176, .480, .830, 1.760, 3.360)*1000,
                       mort = c(0,   .10, .25, .50 , .90),
                       surv = 0)
  mun.mal$surv = 1 - mun.mal$mort
  mun.mal$ppm = mun.mal$conc/1000
  mun.mal$probit = qnorm(mun.mal$mort, mean = 5)
  
  plot(mun.mal$ppm, mun.mal$probit, pch = 16)
  #looks like LC10 was estimated differently, so only fit to parameters derived from L&W model
  mun.mal.sub = subset(mun.mal, probit >= 0) 
    lm.bak.mal = lm(probit ~ ppm, data = mun.mal.sub)
    
    abline(coef(lm.bak.mal), lty = 2)
    
    lc50.bak.mal = (5 - coef(lm.bak.mal)[1]) / coef(lm.bak.mal)[2]
    slp.bak.mal = coef(lm.bak.mal)[2]
    #get standard error from reported 95% CIs of lc50
    se.lc50.bak.mal = mean(c(log10(3.12/lc50.bak.mal), log10(lc50.bak.mal/0.99))) / 1.96
    
  plot(mun.mal$conc, mun.mal$mort, pch = 16, 
       xlab = 'Malathion (ppb)', ylab = 'prop dead', ylim = c(0,1), xlim = c(0,5000),
       main = 'D-R function based on reported values')
    segments(x0 = 3120, y0 = 0.5, x1 = 990, y1 = 0.5)

  fx.bak.mal = function(In, lc = lc50.bak.mal){
    Ins = In/1000
    pnorm(slp.bak.mal * (Ins - lc))
  }
    lines(seq(0,5000,50), sapply(seq(0,5000,50), fx.bak.mal), lty = 2, col = 2)
    lines(seq(0,5000,50), sapply(seq(0,5000,50), fx.bak.mal, lc = 0.99), lty = 2, col = 2)
    lines(seq(0,5000,50), sapply(seq(0,5000,50), fx.bak.mal, lc = 3.12), lty = 2, col = 2)
    
  muNq_mal_Bakry11_uncertainty = function(In){
    if(In == 0) mun = 0 else{
      Ins = (In/1000)
      lc50 = 10^(rnorm(1, log10(lc50.bak.mal), se.lc50.bak.mal))
      mun = pnorm((slp.bak.mal) * (Ins-lc50)) - fx.bak.mal(0)
    }
    if(mun < 0) mun = 0
 
    return(mun)
  }
    points(seq(0,5000,25), sapply(seq(0,5000,25), muNq_mal_Bakry11_uncertainty), 
           pch = 5, col = 4, cex = 0.5)
    
  #keep vector
  keep.bak.mal = c('muNq_mal_Bakry11_uncertainty', 'fx.bak.mal', 'mun.mal',
                   'lc50.bak.mal', 'se.lc50.bak.mal', 'slp.bak.mal')    
    
#Deltamethrin  #############    
mun.del = data.frame(conc = c(.482, 1.21, 2.034, 4.82, 7.26)*1000,
                     mort = c(0,   .10, .25, .50 , .90),
                     surv = 0)
  mun.del$surv = 1 - mun.del$mort
  mun.del$ppm = mun.del$conc/1000
  mun.del$probit = qnorm(mun.del$mort, mean = 5)

plot(mun.del$ppm, mun.del$probit, pch = 16)
  mun.del.sub = subset(mun.del, probit >= 0)
  lm.bak.del = lm(probit ~ ppm, data = mun.del.sub)
  
  abline(coef(lm.bak.del), lty = 2)  
  
  lc50.bak.del = (5 - coef(lm.bak.del)[1]) / coef(lm.bak.del)[2]
  slp.bak.del = coef(lm.bak.del)[2]
  #get standard error from reported 95% CIs of lc50; assuming 0.31 is typo meaning 3.1
  se.lc50.bak.del = mean(c(log10(7.7/lc50.bak.del), log10(lc50.bak.del/3.1))) / 1.96
  
  fx.bak.del = function(In, lc = lc50.bak.del){
    Ins = In/1000
    pnorm(slp.bak.del * (Ins - lc))
  }
  
plot(mun.del$conc, mun.del$mort, pch = 16, ylim = c(0,1), xlim = c(0,10000),
     xlab = 'Deltamethrin (ppb)', ylab = 'prop dead', 
     main = 'D-R function based on reported values')
    segments(x0 = 3100, y0 = 0.5, x1 = 7700, y1 = 0.5)
    
  lines(seq(0,10000,100), sapply(seq(0,10000,100), fx.bak.del), lty = 2, col = 2)
  lines(seq(0,10000,100), sapply(seq(0,10000,100), fx.bak.del, lc = 7.7), lty = 2, col = 2)
  lines(seq(0,10000,100), sapply(seq(0,10000,100), fx.bak.del, lc = 3.1), lty = 2, col = 2)

  muNq_del_Bakry11_uncertainty = function(In){
    if(In == 0) mun = 0 else{
      Ins = (In/1000)
      lc50 = 10^(rnorm(1, log10(lc50.bak.del), se.lc50.bak.del))
      mun = pnorm((slp.bak.del) * (Ins-lc50)) - fx.bak.del(0)
    }
    if(mun < 0) mun = 0
    
    return(mun)
  }
  points(seq(0,10000,50), sapply(seq(0,10000,50), muNq_del_Bakry11_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
  
  #keep vector
  keep.bak.del = c('muNq_del_Bakry11_uncertainty', 'fx.bak.del', 'mun.del',
                   'lc50.bak.del', 'se.lc50.bak.del', 'slp.bak.del')    
  
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
    m = (bakry11$surv[bakry11$chem == 'control'][7] - 
           bakry11$surv[bakry11$chem == 'control'][1]) /
        (bakry11$time[bakry11$chem == 'control'][7] - 
           bakry11$time[bakry11$chem == 'control'][1])

    m.df = data.frame('days' = c(0:28),
                      'control' = 0,
                      'mal' = 0,
                      'delt' = 0)
    
    m.df[1,c(2:4)] = 50
    
    for(i in 1:28){
      m.df[i+1, 2] = m.df[i,2] + m  #subtract deaths per day (from control slope derived above)
      m.df[i+1, 3] = m.df[i,3] - m.df[i,3] * fx.bak.mal(480)  #for agrochemical toxicity, additional mortality as per capita deaths
      m.df[i+1, 4] = m.df[i,4] - m.df[i,4] * fx.bak.del(1210) #for agrochemical toxicity, additional mortality as per capita deaths
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
    
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fn.bak.mal.pred)[1,] / fn.bak$mal.r[1], col = 2, lty=2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fn.bak.mal.pred)[2,] / fn.bak$mal.r[1], col = 2, lty=3)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fn.bak.mal.pred)[3,] / fn.bak$mal.r[1], col = 2, lty=3)
    
  fN.mal.fx.uncertainty = function(In){
    if(In == 0) fn = 1 else{
      Ins = In/1000
      init = predict(fn.bak.mal, newdata = data.frame(mal = Ins), se.fit = T)
    if(init[1] == 0) fn = 0 else{
      fn = rnorm(1, init[1], init[2]) / fn.bak$mal.r[1]
    while(fn < 0 || fn > 1.00000){
        fn = rnorm(1, init[1], init[2]) / fn.bak$mal.r[1]
      }
    }
  }
     return(fn)
} #normalized to 1, upper limit at 1, lower limit at 0
  
  points(seq(0,4000,4), sapply(seq(0,4000,4), fN.mal.fx.uncertainty), pch = 5, cex = 0.5, col = 4)
keep.bak.mal = c(keep.bak.mal, 'fn.bak', 'fn.bak.mal', 'fN.mal.fx.uncertainty')
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
  
    lines(seq(0,7000,10), sapply(seq(0,7000,10), fn.bak.del.pred)[1,] / fn.bak$del.r[1], col = 2, lty=2)
    lines(seq(0,7000,10), sapply(seq(0,7000,10), fn.bak.del.pred)[2,] / fn.bak$del.r[1], col = 2, lty=3)
    lines(seq(0,7000,10), sapply(seq(0,7000,10), fn.bak.del.pred)[3,] / fn.bak$del.r[1], col = 2, lty=3)
  
  fNq_del_Bakry11_uncertainty = function(In){
    if(In == 0) fn = 1 else{
      Ins = In/1000
    init = predict(fn.bak.del, newdata = data.frame(del = Ins), se.fit = T)
    if(init[1] == 0) fn = 0 else{
      fn = rnorm(1, init[1], init[2]) / fn.bak$del.r[1]
    while(fn < 0 || fn > 1.00000){
        fn = rnorm(1, init[1], init[2]) / fn.bak$del.r[1]
      }
    }
  }
    return(fn)
} 
  
  points(seq(0,7000,7), sapply(seq(0,7000,7), fNq_del_Bakry11_uncertainty), pch = 5, cex = 0.5, col = 4)
keep.bak.del = c(keep.bak.del, 'fn.bak', 'fn.bak.del', 'fNq_del_Bakry11_uncertainty')  
#full keep vector ###############  
  keep.bak11.N = c('mun.mal', 'muNq_mal_Bakry11_uncertainty', 'muNq_del_Bakry11_uncertainty', 'fn.bak',
                   'fN.mal.fx.uncertainty', 'fn.bak.mal', 'fNq_del_Bakry11_uncertainty', 'fn.bak.del')