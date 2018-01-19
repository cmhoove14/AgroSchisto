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
  mun.mal = data.frame(conc = c(0, .176, .480, .830, 1.760, 3.360)*1000,
                       mort = c(0, 0,   .10, .25, .50 , .90),
                       surv = 0)
  mun.mal$surv = 1 - mun.mal$mort
  mun.mal$ppm = mun.mal$conc/1000
  mun.mal$ppmlog10 = log10(mun.mal$ppm)
  mun.mal$probit = qnorm(mun.mal$mort, mean = 5)
  
plot(mun.mal$conc, mun.mal$mort, pch = 16, 
     xlab = 'Malathion (ppb)', ylab = 'prop dead', ylim = c(0,1), xlim = c(0,5000),
     main = 'D-R function based on reported values')
  segments(x0 = 3120, y0 = 0.5, x1 = 990, y1 = 0.5)
  
  lc50.bak.mal.report = 1.760
  slp.bak.mal.report = 2.74
  #get standard error from reported 95% CIs of lc50
  se.lc50.bak.mal = mean(c(log10(3.12/lc50.bak.mal.report), log10(lc50.bak.mal.report/0.99))) / 1.96
  
muNq_mal_Bakry11_uncertainty = function(In){
  #if(In == 0) mun = 0 else{
  Ins = (In/1000)
  lc50 = 10^(rnorm(1, log10(lc50.bak.mal.report), se.lc50.bak.mal))
    mun = pnorm((slp.bak.mal.report) * log10(Ins/lc50))
    #}
    #if(mun < 0) mun = 0
    
    return(mun)
  }
  
  points(seq(0,5000,25), sapply(seq(0,5000,25), muNq_mal_Bakry11_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
  
  #keep vector
  keep.bak.mal = c('muNq_mal_Bakry11_uncertainty', 'mun.mal',
                   'lc50.bak.mal.report', 'se.lc50.bak.mal', 'slp.bak.mal.report')    
  
#Deltamethrin  #############    
mun.del = data.frame(conc = c(0, .482, 1.21, 2.034, 4.82, 7.26)*1000,
                     mort = c(0, 0,   .10, .25, .50 , .90),
                     surv = 0)
  mun.del$surv = 1 - mun.del$mort
  mun.del$ppm = mun.del$conc/1000
  mun.del$probit = qnorm(mun.del$mort, mean = 5)

  plot(mun.del$conc, mun.del$mort, pch = 16, ylim = c(0,1), xlim = c(0,10000),
       xlab = 'Deltamethrin (ppb)', ylab = 'prop dead', 
       main = 'D-R function based on reported values')
  segments(x0 = 3100, y0 = 0.5, x1 = 7700, y1 = 0.5)
  
  
  lc50.bak.del.report = 4.82
  slp.bak.del.report = 2.74
  se.lc50.bak.del = mean(c(log10(7.7/4.82), log10(4.82/3.1))) / 1.96
  
  muNq_del_Bakry11_uncertainty = function(In){
    #if(In == 0) mun = 0 else{
    Ins = (In/1000)
    lc50 = 10^(rnorm(1, log10(lc50.bak.del.report), se.lc50.bak.mal))
    mun = pnorm((slp.bak.del.report) * log10(Ins/lc50))
    #}
    #if(mun < 0) mun = 0
    
    return(mun)
  }
  
  points(seq(0,10000,50), sapply(seq(0,10000,50), muNq_del_Bakry11_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
  
  #keep vector
  keep.bak.del = c('muNq_del_Bakry11_uncertainty', 'mun.del',
                   'lc50.bak.del.report', 'se.lc50.bak.del', 'slp.bak.del.report')    
    
#Now lets look at reproduction over time ###########
  #The paper only provides mean eggs / snail over 4 weeks exposure to lc10 of each insecticide
  #This isn't enough data to generate d-r function without some assumptions
  #So we'll assume there's no reproduction at lc90 and fit 2-parameter logistic function with 1 degree freedom
  #Values for controls are equal to the sum of mean eggs per day produced in the control up to the day where
  #reproduction stops in the lc10 test groups (day4 for malathion, day 14 for deltamethrin)    
bakry11<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/bakry2011.csv')
  
#get mean eggs/snail/day across all time points for each treatment   
  fn.ctrl11 = mean(bakry11$eggs_snail_day[bakry11$chem == 'control' & bakry11$surv != 0])  
    fn.ctrl.sd11 = sd(bakry11$eggs_snail_day[bakry11$chem == 'control' & bakry11$surv != 0])
  
  fn.mal = mean(bakry11$eggs_snail_day[bakry11$chem == 'malathion' & bakry11$surv != 0])  
    fn.mal.sd = sd(bakry11$eggs_snail_day[bakry11$chem == 'malathion' & bakry11$surv != 0])
      
  fn.del = mean(bakry11$eggs_snail_day[bakry11$chem == 'deltamethrin' & bakry11$surv != 0])  
    fn.del.sd = sd(bakry11$eggs_snail_day[bakry11$chem == 'deltamethrin' & bakry11$surv != 0])
      
#Functions that estimate parameter value in single concentration group #########
    fNq_mal_Bakry11_uncertainty = function(waste){ #function based on snail egg masses
      wst = waste
      fN = rnorm(1, fn.mal, fn.mal.sd) / fn.ctrl11
      while(fN < 0) fN = rnorm(1, fn.mal, fn.mal.sd) / fn.ctrl11
      fN
    }
    
    #hist(sapply(rep(1,1000), fNq_mal_Bakry11_uncertainty, simplify = T))
    
    fNq_del_Bakry11_uncertainty = function(waste){ #function based on snail hatchlings
      wst = waste
      fN = rnorm(1, fn.del, fn.del.sd) / fn.ctrl11
      while(fN < 0) fN = rnorm(1, fn.del, fn.del.sd) / fn.ctrl11
      fN
    }  
    
    #hist(sapply(rep(1,1000), fNq_del_Bakry11_uncertainty, simplify = T))
    
  
keep.bak.mal = c(keep.bak.mal, 'fNq_mal_Bakry11_uncertainty',
                 'fn.mal', 'fn.mal.sd', 'fn.ctrl11')
  
  keep.bak.del = c(keep.bak.del, 'fNq_del_Bakry11_uncertainty',
                   'fn.del', 'fn.del.sd', 'fn.ctrl11')  
#full keep vector ###############  
  keep.bak11.N = c(keep.bak.mal, keep.bak.del)