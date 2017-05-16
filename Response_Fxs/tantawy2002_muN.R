#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
#Data extraction and model fitting to Tantawy 2002 data
require(drc)

#Snail toxicity ##########
#Butachlor ##########
snail.but = data.frame(conc = c(0.65, 1.5, 4.5,6.5,44)*1000, # suspected typo in paper @4.5 (reads 45)
                       mort = c(0, .10, .25, .50, .90),
                       dead = c(0, 1, 2.5, 5, 9),
                       surv = 0)
  snail.but$ppm = snail.but$conc / 1000
  snail.but$probit = qnorm(snail.but$mort, mean = 5)
  snail.but$log10ppm = log10(snail.but$ppm)
  snail.but$surv = 1 - snail.but$mort 
  
plot(snail.but$log10ppm, snail.but$probit, pch = 16,
     xlab = 'butralin log10(ppb)', ylab = 'probit mortality')
snail.but2 = subset(snail.but, probit > 0)
  tant.but.lm = lm(probit ~ log10ppm, data = snail.but2)
  abline(coef(tant.but.lm), lty = 2)
  
  lc50.tant.but = 10^((5 - coef(tant.but.lm)[1]) / coef(tant.but.lm)[2])
  slp.tant.but = coef(tant.but.lm)[2]
  #get standard error from reported 95% CIs of lc50
  se.lc50.tant.but = mean(c(log10(10.4 / 6.5), log10(6.5/4.06))) / 1.96 #st. err of lc50 in ppm
  
plot(snail.but$conc, snail.but$mort, pch = 16, ylim = c(0,1), xlim = c(0,max(snail.but$conc)+100),
     xlab = 'Butachlor (ppb)', ylab = 'Snail mortality',
     main = expression(paste('Butachlor toxicity to ', italic('Bi. alexandrina'))))
  segments(x0 = 4060, x1 = 10400, y0 = 0.5, y1 = 0.5)
  
  fx.tant.but = function(He, lc = lc50.tant.but){
    heu = (He/1000)
    pnorm(slp.tant.but * log10(heu/(lc)))
  }  
  
  lines(seq(0, 50000, 50), sapply(seq(0, 50000, 50), fx.tant.but), lty = 2, col = 2)
  lines(seq(0, 50000, 50), sapply(seq(0, 50000, 50), fx.tant.but, lc = 4.06), lty = 3, col = 2)
  lines(seq(0, 50000, 50), sapply(seq(0, 50000, 50), fx.tant.but, lc = 10.4), lty = 3, col = 2)
  
  muN.tant.but_uncertainty<-function(He){
    if(He == 0) mun = 0 else{
      heu = (He/1000)
      lc50 = 10^(rnorm(1, log10(lc50.tant.but), se.lc50.tant.but))
      while(lc50 < 0){
        lc50 = 10^(rnorm(1, log10(lc50.tant.but), se.lc50.tant.but))
      }  
      mun = pnorm((slp.tant.but) * log10(heu/lc50))
    }
    while(mun < 0){
      lc50 = 10^(rnorm(1, log10(lc50.tant.but), se.lc50.tant.but))
      mun = pnorm((slp.tant.but) * log10(heu/lc50))
    } 
    return(mun)
  }
  
  points(seq(0, 50000, 100), sapply(seq(0, 50000, 100), muN.tant.but_uncertainty),
         pch = 5, cex = 0.5, col = 4)
keep.tant.but = c('muN.tant.but_uncertainty', 'lc50.tant.but', 'se.lc50.tant.but', 'slp.tant.but')  
#Fluazifop-p-butyl  ##########     
snail.fpb = data.frame(conc = c(1.76, 4.5, 9,17.6,58)*1000, 
                       mort = c(0, .10, .25, .50, .90),
                       dead = c(0, 1, 2.5, 5, 9),
                       surv = 0)
  snail.fpb$ppm = snail.fpb$conc / 1000
  snail.fpb$probit = qnorm(snail.fpb$mort, mean = 5)
  snail.fpb$log10ppm = log10(snail.fpb$ppm)
  snail.fpb$surv = 1 - snail.fpb$mort 

plot(snail.fpb$log10ppm, snail.fpb$probit, pch = 16,
     xlab = 'fluazifop-p-butyl log10(ppm)', ylab = 'probit mortality')
  snail.fpb2 = subset(snail.fpb, probit > 0)
  tant.fpb.lm = lm(probit ~ log10ppm, data = snail.fpb2)
    abline(coef(tant.fpb.lm), lty = 2)

  lc50.tant.fpb = 10^((5 - coef(tant.fpb.lm)[1]) / coef(tant.fpb.lm)[2])
  slp.tant.fpb = coef(tant.fpb.lm)[2]
  #get standard error from reported 95% CIs of lc50
  se.lc50.tant.fpb = mean(c(log10(26.4 / 17.6), log10(17.6/11.73))) / 1.96 #st. err of lc50 in ppm

plot(snail.fpb$conc, snail.fpb$mort, pch = 16, ylim = c(0,1), xlim = c(0,max(snail.fpb$conc)+100),
     xlab = 'fluazifop-p-butyl (ppb)', ylab = 'Snail mortality',
     main = expression(paste('fluazifop-p-butyl toxicity to ', italic('Bi. alexandrina'))))
  segments(x0 = 11730, x1 = 26400, y0 = 0.5, y1 = 0.5)

fx.tant.fpb = function(He, lc = lc50.tant.fpb){
  heu = (He/1000)
  pnorm(slp.tant.fpb * log10(heu/(lc)))
}  
  
  lines(seq(0, 60000, 50), sapply(seq(0, 60000, 50), fx.tant.fpb), lty = 2, col = 2)
  lines(seq(0, 60000, 50), sapply(seq(0, 60000, 50), fx.tant.fpb, lc = 11.73), lty = 3, col = 2)
  lines(seq(0, 60000, 50), sapply(seq(0, 60000, 50), fx.tant.fpb, lc = 26.4), lty = 3, col = 2)

muN.tant.fpb_uncertainty<-function(He){
  if(He == 0) mun = 0 else{
    heu = (He/1000)
    lc50 = 10^(rnorm(1, log10(lc50.tant.fpb), se.lc50.tant.fpb))
    while(lc50 < 0){
      lc50 = 10^(rnorm(1, log10(lc50.tant.fpb), se.lc50.tant.fpb))
    }  
    mun = pnorm((slp.tant.fpb) * log10(heu/lc50))
  }
  while(mun < 0){
    lc50 = 10^(rnorm(1, log10(lc50.tant.fpb), se.lc50.tant.fpb))
    mun = pnorm((slp.tant.fpb) * log10(heu/lc50))
  } 
  return(mun)
}

  points(seq(0, 60000, 100), sapply(seq(0, 60000, 100), muN.tant.fpb_uncertainty),
         pch = 5, cex = 0.5, col = 4)
  
keep.tant.fpb = c('muN.tant.fpb_uncertainty', 'fx.tant.fpb', 'lc50.tant.fpb', 
                  'se.lc50.tant.fpb', 'slp.tant.fpb')  
#keep vector #########  
keep.muN.tantawy = c('muN.tant.but_uncertainty', 'fx.tant.but', 'lc50.tant.but', 
                     'se.lc50.tant.but', 'slp.tant.but',
                     'muN.tant.fpb_uncertainty', 'fx.tant.fpb', 'lc50.tant.fpb', 
                     'se.lc50.tant.fpb', 'slp.tant.fpb')  