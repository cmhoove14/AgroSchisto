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
  sn$log10ppm = log10(sn$conc+1)
  sn$probit = qnorm(sn$prop_dead, mean = 5)
  
  lc50.tch91.mal = 202.93
  slp.tch91.mal = 3.59
  se.lc50.tch91.mal = mean(c(log10(244.16 / lc50.tch91.mal), log10(lc50.tch91.mal/158.91))) / 1.96
  
  muNq_mal_tch91_uncertainty = function(In){
      ins = (In/1000)
      lc50 = 10^(rnorm(1, log10(lc50.tch91.mal), se.lc50.tch91.mal))
      mun = pnorm((slp.tch91.mal) * log10(ins/lc50))

    return(mun)
  }
  
  plot(sn$conc_ppb, sn$prop_dead, pch = 16, xlim = c(0, 1e6), ylim = c(0,1),
       xlab = 'malathion (ppb)', ylab = 'mortality')
    segments(x0 = 158.91*1000, x1 = 244.16*1000, y0 = 0.5, y1 = 0.5)
    points(seq(0, 1e6, 5e3), sapply(seq(0, 1e6, 5e3), muNq_mal_tch91_uncertainty), 
           pch = 5, col = 4, cex = 0.5)
  
#egg viability ###########
eg = data.frame(time = 1,
                chem = 'malathion',
                conc = c(0, 42.63, 94.78, 210.73),
                prop_dead = c(0, .05, .5,.95))
  #eg = subset(eg, conc > 0)
  eg$conc_ppb = eg$conc*1000
  eg$log10ppm = log10(eg$conc)
  eg$probit = qnorm(eg$prop_dead, mean = 5)
  eg$prop_surv = 1-eg$prop_dead
  
  lc50.tch91.mal.eg = 94.78
  slp.tch91.mal.eg = 4.74
  #get standard error from reported 95% CIs of lc50
  se.lc50.tch91.mal.eg = mean(c(log10(143.97 / lc50.tch91.mal.eg), log10(lc50.tch91.mal.eg/35.49))) / 1.96 #st. err of lc50 in ppm
  
  fNq_mal_tch91_uncertainty = function(In){
      ins = (In/1000)
      lc50 = 10^(rnorm(1, log10(lc50.tch91.mal.eg), se.lc50.tch91.mal.eg))
      fn = pnorm((-slp.tch91.mal.eg) * log10(ins/lc50))#neg. slope - normailize to 
    
    return(fn)
  }
  
  plot(eg$conc_ppb, eg$prop_surv, pch = 16, xlim = c(0, 4e5), ylim = c(0,1),
       xlab = 'malathion (ppb)', ylab = 'egg viability')
    segments(x0 = 35.49*1000, x1 = 143.97*1000, y0 = 0.5, y1 = 0.5)
  
  points(seq(0, 4e5, 1e3), sapply(seq(0, 4e5, 1e3), fNq_mal_tch91_uncertainty), 
         pch = 5, col = 4, cex = 0.5)

#Keep vector
  keep.tch91.snail = c('sn', 'lc50.tch91.mal', 'se.lc50.tch91.mal', 'slp.tch91.mal',
                       'eg', 'lc50.tch91.mal.eg', 'se.lc50.tch91.mal.eg', 'slp.tch91.mal.eg',
                       'fNq_mal_tch91_uncertainty', 'muNq_mal_tch91_uncertainty')