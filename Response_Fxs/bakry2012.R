#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Bakry 2012 data
require(drc)

#Snail (B. alexandrina) toxicity; create data frame, clean, etc. ##########
muN.bak = data.frame(atr = c(.330, 1.250, 4.750),
                     atr.se = 0,
                     gly = c(.840, 3.150, 12.600),
                     gly.se = 0,
                     mort = c(0.1, 0.5, 0.9))

  muN.bak$surv = 1 - muN.bak$mort

#Standard errors imputed assuming proportional error to concentration
  
  se.atr = (log10(1.88) - log10(1.25)) / 1.96 #st. err of lc50 in ppm
  se.gly = (log10(4.82) - log10(3.15)) / 1.96 #error bars are assymetrical in this data even after log transformation
                                              #so not entirely sure how to handle that...

  muN.bak$atr.se = (muN.bak$atr / 1.25) * se.atr #st. err proportional to concentration
  muN.bak$gly.se = (muN.bak$gly / 3.15) * se.gly #st. err proportional to concentration

#direct snail toxicity from Atrazine ############  
plot(muN.bak$atr, muN.bak$mort, ylim = c(0,1), xlim = c(0,13),
       pch = 16, xlab = 'atrazine (ppm)', ylab = 'mortality')
  for(i in 1:length(muN.bak$atr)){
    segments(y0 = muN.bak$mort[i], x0 = muN.bak$atr[i] + muN.bak$atr.se[i],
             y1 = muN.bak$mort[i], x1 = muN.bak$atr[i] - muN.bak$atr.se[i])

  }
  
  bak12.atr.drm = drm(mort ~ atr, weights = rep(50,3), data = muN.bak, type = 'binomial',
                      fct = LL.2())
  
    bak12.atr.pred = function(He){
      predict(bak12.atr.drm, newdata = data.frame(atr = He), interval = 'confidence', level = 0.95)
    }
    
  lines(seq(0,13, 0.1), sapply(seq(0,13, 0.1), bak12.atr.pred, simplify = T)[1,], col = 2, lty = 2)
  lines(seq(0,13, 0.1), sapply(seq(0,13, 0.1), bak12.atr.pred, simplify = T)[2,], col = 2, lty = 3)
  lines(seq(0,13, 0.1), sapply(seq(0,13, 0.1), bak12.atr.pred, simplify = T)[3,], col = 2, lty = 3)
  
  mu_Nq_atr_bak12_uncertainty<-function(He){
    He.u = He/1000
    rdrm(1, LL.2(), coef(bak12.atr.drm), He.u, yerror = 'rbinom', ypar = 50)$y / 50  
  }
  
  points(seq(0,13, 0.1), sapply(seq(0,13, 0.1)*1000, mu_Nq_atr_bak12_uncertainty, simplify = T),
         col = 4, pch = 5, cex = 0.5)
  
#direct snail toxicity from Glyphosate ############  
  plot(muN.bak$gly, muN.bak$mort, ylim = c(0,1), xlim = c(0,13),
       pch = 16, xlab = 'glyphosate (ppm)', ylab = 'mortality')  
    for(i in 1:length(muN.bak$atr)){
      segments(y0 = muN.bak$mort[i], x0 = muN.bak$gly[i] + muN.bak$gly.se[i],
               y1 = muN.bak$mort[i], x1 = muN.bak$gly[i] - muN.bak$gly.se[i])
    }
    
  bak12.gly.drm = drm(mort ~ gly, weights = rep(50,3), data = muN.bak, type = 'binomial',
                      fct = LL.2())
  
    bak12.gly.pred = function(He){
      predict(bak12.gly.drm, newdata = data.frame(gly = He), interval = 'confidence', level = 0.95)
    }

    lines(seq(0,13, 0.1), sapply(seq(0,13, 0.1), bak12.gly.pred, simplify = T)[1,], col = 2, lty = 2)
    lines(seq(0,13, 0.1), sapply(seq(0,13, 0.1), bak12.gly.pred, simplify = T)[2,], col = 2, lty = 3)
    lines(seq(0,13, 0.1), sapply(seq(0,13, 0.1), bak12.gly.pred, simplify = T)[3,], col = 2, lty = 3)
    
    mu_Nq_gly_bak12_uncertainty<-function(He){
      He.u = He/1000
      rdrm(1, LL.2(), coef(bak12.gly.drm), He.u, yerror = 'rbinom', ypar = 50)$y / 50  
    }
    
    points(seq(0,13, 0.1), sapply(seq(0,13, 0.1)*1000, mu_Nq_gly_bak12_uncertainty, simplify = T),
           col = 4, pch = 5, cex = 0.5)
    
  
  keep.bak12 = c('muN.bak', 'mu_Nq_atr_bak12_uncertainty', 'mu_Nq_gly_bak12_uncertainty', 
                 'bak12.atr.drm', 'bak12.gly.drm')

#Now lets look at reproduction over time ###########
bakry12<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/bakry2012.csv')
  bak12 = subset(bakry12, time != 0)
  
sn.wk = as.numeric()  #snails/week estimates for each treatment
  for(i in 1:length(unique(bak12$conc))){
    sn.wk[i] = sum(bak12$prop_surv[bak12$conc == unique(bak12$conc)[i]])
  }  
  
hatches = as.numeric()
  for(i in 1:length(unique(bak12$conc))){
    hatches[i] = sum(bak12$hatch[bak12$conc == unique(bak12$conc)[i]])
  } 

fn = hatches / sn.wk

bak.fn = data.frame(gly.conc = c(0, 840),
                    atr.conc = c(0, 330),
                    gly.rep = c(fn[1], fn[3]),
                    atr.rep = c(fn[1], fn[2]))  
#No function attempted to fit because only control and single concentration data point for each agroC