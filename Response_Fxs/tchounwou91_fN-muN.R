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

#direct mortality to snails #############
sn = read.csv('~/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/tchounwou1991.csv')
  sn = subset(sn, conc > 0)
  sn$conc_ppb = sn$conc*1000
  sn$mort = 20 - sn$surv
  sn$log10ppm = log10(sn$conc)
  sn$probit = qnorm(sn$prop_dead, mean = 5)
  
plot(sn$log10ppm, sn$probit, pch = 16, ylim = c(2.5,7.5), xlim = c(0, log10(max(sn$conc)+100)),
     xlab = 'Malathion log10(ppm)', ylab = 'probit mortality', 
     main = 'Tchounwou 91 B. havenensis mortality data')

  lm.tch91.mal = lm(probit ~ log10ppm, data = sn)
    abline(coef(lm.tch91.mal), lty = 2)
  #add uncertainty of lc5, lc50, and lc95 from paper
    segments(x0 = log10(c(33.97, 158.91, 442.71)),
             x1 = log10(c(102.02, 244.16, 767.641)), #767 inferred b/ of typo in manuscript
             y0 = qnorm(c(0.05, 0.5, 0.95), mean = 5),
             y1 = qnorm(c(0.05, 0.5, 0.95), mean = 5))
  
  fx.tch91.mal = function(In){
    ins = log10(In/1000)
    predict(lm.tch91.mal, newdata = data.frame(log10ppm = ins), interval = 'confidence', level = 0.95)
  } 
    lines(seq(0,3.5, 0.01), sapply(10^(seq(0,3.5, 0.01))*1000, fx.tch91.mal)[2,], lty = 3)
    lines(seq(0,3.5, 0.01), sapply(10^(seq(0,3.5, 0.01))*1000, fx.tch91.mal)[3,], lty = 3)
    
  lc50.tch91.mal = 10^((5 - coef(lm.tch91.mal)[1]) / coef(lm.tch91.mal)[2])
    slp.tch91.mal = coef(lm.tch91.mal)[2]
    #get standard error from reported 95% CIs of lc50
    se.lc50.tch91.mal = mean(c(log10(244.16 / lc50.tch91.mal), log10(lc50.tch91.mal/158.91))) / 1.96 #st. err of lc50 in ppm

#Model of adult daily mortality ###########
muNq_mal_tch91_uncertainty = function(In){
  if(In == 0) mun = 0 else{
    ins = (In/1000)
    lc50 = 10^(rnorm(1, log10(lc50.tch91.mal), se.lc50.tch91.mal))
  while(lc50 < 0){
      lc50 = 10^(rnorm(1, log10(lc50.tch91.mal), se.lc50.tch91.mal))
    }  
    mun = pnorm((slp.tch91.mal) * log10(ins/lc50))
  }
  if(mun < 0) mun = 0
     
  return(mun)
}

plot(sn$conc_ppb, sn$prop_dead, pch = 16, xlim = c(0, 1e6), ylim = c(0,1),
     xlab = 'malathion (ppb)', ylab = 'mortality')
  points(seq(0, 1e6, 5e3), sapply(seq(0, 1e6, 5e3), muNq_mal_tch91_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
  
#egg viability ###########
eg = read.csv('~/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/tchounwou1991.csv')
  eg = subset(eg, conc > 0)
  eg$conc_ppb = eg$conc*1000
  eg$log10ppm = log10(eg$conc)
  eg$probit = qnorm(eg$prop_dead, mean = 5)
  eg$conc_ppb = eg$conc*1000
 
plot(eg$log10ppm, eg$probit, pch = 16, ylim = c(2.5,7.5), xlim = c(0, log10(max(sn$conc)+100)),
     xlab = 'Malathion log10(ppm)', ylab = 'probit mortality', 
     main = 'Tchounwou 91 B. havenensis egg viability data')

  lm.tch91.mal.eg = lm(probit ~ log10ppm, data = eg)
  abline(coef(lm.tch91.mal.eg), lty = 2)
  
#add uncertainty of lc5, lc50, and lc95 from paper
  segments(x0 = log10(c(1.23, 35.49, 139.99)),
           x1 = log10(c(71.72, 143.97, 2119.52)), #767 inferred b/ of typo in manuscript
           y0 = qnorm(c(0.05, 0.5, 0.95), mean = 5),
           y1 = qnorm(c(0.05, 0.5, 0.95), mean = 5))
  
  fx.tch91.mal.eg = function(In){
    ins = log10(In/1000)
    predict(lm.tch91.mal.eg, newdata = data.frame(log10ppm = ins), interval = 'confidence', level = 0.95)
  } 
  lines(seq(0,3.5, 0.01), sapply(10^(seq(0,3.5, 0.01))*1000, fx.tch91.mal.eg)[2,], lty = 3)
  lines(seq(0,3.5, 0.01), sapply(10^(seq(0,3.5, 0.01))*1000, fx.tch91.mal.eg)[3,], lty = 3)
  
  lc50.tch91.mal.eg = 10^((5 - coef(lm.tch91.mal.eg)[1]) / coef(lm.tch91.mal.eg)[2])
  slp.tch91.mal.eg = coef(lm.tch91.mal.eg)[2]
  #get standard error from reported 95% CIs of lc50
  se.lc50.tch91.mal.eg = mean(c(log10(143.97 / lc50.tch91.mal.eg), log10(lc50.tch91.mal.eg/35.49))) / 1.96 #st. err of lc50 in ppm
  
#Model of egg viability###########
  fNq_mal_tch91_uncertainty = function(In){
    if(In == 0) fn = 1 else{
      ins = (In/1000)
      lc50 = 10^(rnorm(1, log10(lc50.tch91.mal.eg), se.lc50.tch91.mal.eg))
      while(lc50 < 0){
        lc50 = 10^(rnorm(1, log10(lc50.tch91.mal.eg), se.lc50.tch91.mal.eg))
      }  
      fn = pnorm((-slp.tch91.mal.eg) * log10(ins/lc50))#neg. slope - normailize to 
    }
    if(fn < 0) fn = 0
      
    return(fn)
  }
  
  plot(eg$conc_ppb, eg$prop_surv, pch = 16, xlim = c(0, 4e5), ylim = c(0,1),
       xlab = 'malathion (ppb)', ylab = 'egg viability')
  points(seq(0, 4e5, 1e3), sapply(seq(0, 4e5, 1e3), fNq_mal_tch91_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
  
#Keep vector
  keep.tch91.snail = c('sn', 'lc50.tch91.mal', 'se.lc50.tch91.mal', 'slp.tch91.mal',
                       'eg', 'lc50.tch91.mal.eg', 'se.lc50.tch91.mal.eg', 'slp.tch91.mal.eg',
                       'fNq_mal_tch91_uncertainty', 'muNq_mal_tch91_uncertainty')