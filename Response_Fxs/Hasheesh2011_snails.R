#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Hasheesh and Mohamed 2011 data and analysis assessing toxicity of Chlorpyrifos and Profenofos ###############
# Direct toxicity to snails (Bu. truncatus); daily mortality rate ########

#Chlorpyrifos #########
mun.hsh.ch = data.frame(conc = c(.72, 1.32, 2.82)*1000,
                        mort = c(.25, .50 , .90),
                        surv = 0)
  mun.hsh.ch$surv = 1 - mun.hsh.ch$mort
  mun.hsh.ch$probit = qnorm(mun.hsh.ch$mort, mean = 5)
  mun.hsh.ch$ppm = mun.hsh.ch$conc/1000
  mun.hsh.ch$log10 = log10(mun.hsh.ch$ppm)
  
plot(mun.hsh.ch$ppm, mun.hsh.ch$probit, pch = 16)  
  hash.ch.lm = lm(probit ~ ppm, data = mun.hsh.ch)
  
  abline(coef(hash.ch.lm), lty = 2)
  
  lc50.ch.hash = (5 - coef(hash.ch.lm)[1]) / coef(hash.ch.lm)[2]
  slp.ch.hash = coef(hash.ch.lm)[2]
  se.lc50.ch.hash = mean(c(1.98-lc50.ch.hash , lc50.ch.hash-0.88)) / 1.96
 
plot(mun.hsh.ch$conc, mun.hsh.ch$mort, pch = 16, ylim = c(0,1), xlim = c(0,3500),
     xlab = 'Chlorpyrifos (ppb)', ylab = 'prop dead')
  segments(x0 = 880, y0 = 0.5, x1 = 1980, y1 = 0.5)
  
fx.hsh.ch = function(In, lc = lc50.ch.hash){
  Ins = In/1000
  pnorm(slp.ch.hash * (Ins - lc))
}
  lines(seq(0,3500,10), sapply(seq(0,3500,10), fx.hsh.ch), lty = 2, col = 2)
  lines(seq(0,3500,10), sapply(seq(0,3500,10), fx.hsh.ch, lc = 1.98), lty = 2, col = 2)
  lines(seq(0,3500,10), sapply(seq(0,3500,10), fx.hsh.ch, lc = 0.88), lty = 2, col = 2)
  
muNq_ch_hash11_uncertainty = function(In){
  if(In == 0) mun = 0 else{
    Ins = (In/1000)
    lc50 = (rnorm(1, lc50.ch.hash, se.lc50.ch.hash))
    mun = pnorm((slp.ch.hash) * (Ins-lc50)) - fx.hsh.ch(0)
  }
  while(mun < 0){
    lc50 = (rnorm(1, lc50.ch.hash, se.lc50.ch.hash))
    mun = pnorm((slp.ch.hash) * (Ins-lc50)) - fx.hsh.ch(0)
  } 
  return(mun)
}
points(seq(0,3500,10), sapply(seq(0,3500,10), muNq_ch_hash11_uncertainty), 
       pch = 5, col = 4, cex = 0.5)

#keep vector
keep.hsh.ch = c('muNq_ch_hash11_uncertainty', 'fx.hsh.ch',
                 'lc50.ch.hash', 'se.lc50.ch.hash', 'slp.ch.hash')    

#Profenofos #########
mun.hsh.prof = data.frame(conc = c(1.4, 2.5, 3.72)*1000,
                    mort = c(.25, .50 , .90),
                    surv = 0)
  mun.hsh.prof$surv = 1 - mun.hsh.prof$mort
  mun.hsh.prof$ppm = mun.hsh.prof$conc/1000
  mun.hsh.prof$log10 = log10(mun.hsh.prof$ppm)
  mun.hsh.prof$probit = qnorm(mun.hsh.prof$mort, mean = 5)
    
plot(mun.hsh.prof$ppm, mun.hsh.prof$probit, pch = 16)
  hash.prof.lm = lm(probit ~ ppm, data = mun.hsh.prof)
  
  abline(coef(hash.prof.lm), lty = 2)
  
  lc50.prof.hash = (5 - coef(hash.prof.lm)[1]) / coef(hash.prof.lm)[2]
  slp.prof.hash = coef(hash.prof.lm)[2]
  se.lc50.prof.hash = mean(c(3.33-lc50.prof.hash , lc50.prof.hash-1.88)) / 1.96
  
fx.hsh.prof = function(In, lc = lc50.prof.hash){
  Ins = In/1000
  pnorm(slp.prof.hash * (Ins - lc))
}  
    
  plot(mun.hsh.prof$conc, mun.hsh.prof$mort, pch = 16, ylim = c(0,1), xlim = c(0,5000),
       xlab = 'Profenofos (ppm)', ylab = 'prop dead')
      segments(x0 = 1880, y0 = 0.5, x1 = 3330, y1 = 0.5)
      
    lines(seq(0,5000,50), sapply(seq(0,5000,50), fx.hsh.prof), lty = 2, col = 2)
    lines(seq(0,5000,50), sapply(seq(0,5000,50), fx.hsh.prof, lc = 3.33), lty = 2, col = 2)
    lines(seq(0,5000,50), sapply(seq(0,5000,50), fx.hsh.prof, lc = 1.88), lty = 2, col = 2)

muNq_prof_hash11_uncertainty = function(In){
  if(In == 0) mun = 0 else{
    Ins = (In/1000)
    lc50 = (rnorm(1, lc50.prof.hash, se.lc50.prof.hash))
    mun = pnorm((slp.prof.hash) * (Ins-lc50)) - fx.hsh.prof(0)
  }
  while(mun < 0){
    lc50 = (rnorm(1, lc50.prof.hash, se.lc50.prof.hash))
    mun = pnorm((slp.prof.hash) * (Ins-lc50)) - fx.hsh.prof(0)
  } 
  return(mun)
}
points(seq(0,5000,10), sapply(seq(0,5000,10), muNq_prof_hash11_uncertainty), 
       pch = 5, col = 4, cex = 0.5)

#keep vector
keep.hsh.ch = c('muNq_prof_hash11_uncertainty', 'fx.hsh.prof',
                'lc50.prof.hash', 'se.lc50.prof.hash', 'slp.prof.hash')    

#Direct toxicity to snails affecting reproduction (Table 2) #######
fn<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail reproduction/Hasheesh2011_snail_mort_repro_weekly.csv')

#longitudinal measurement of mean eggs/snail          
  plot(fn$time[fn$chem == 'control'], fn$eggs_snail_day[fn$chem == 'control'], type = 'l', lwd=2,
        xlab = 'time (weeks)', ylab = 'eggs/snail/day', ylim = c(0,max(fn$eggs_snail_day)+2),
        main = 'Snail reproduction over time Hasheesh 2011')
    lines(fn$time[fn$chem == 'chlorpyrifos'], fn$eggs_snail_day[fn$chem == 'chlorpyrifos'], col = 'red', lwd=2)
    lines(fn$time[fn$chem == 'profenofos'], fn$eggs_snail_day[fn$chem == 'profenofos'], col = 'orange', lwd=2) 
    legend('topleft', legend = c('control', 'chlorpyrifos', 'profenofos'),
           lwd = 2, col = c(1,2, 'orange'), cex = 0.8)
    
#get estimate of mean eggs/snail/day for each treatment
  ctrl.eggs = mean(fn$eggs_snail_day[fn$chem == 'control' & fn$surv != 0])
  chlor.eggs = mean(fn$eggs_snail_day[fn$chem == 'chlorpyrifos' & fn$surv != 0])
  prof.eggs = mean(fn$eggs_snail_day[fn$chem == 'profenofos' & fn$surv != 0])

   
fn.hash = data.frame(chlor.conc = c(0, 720, 2820),
                     prof.conc = c(0, 1400, 3720),
                     chlor.rep = c(ctrl.eggs, chlor.eggs, 0),
                     prof.rep = c(ctrl.eggs, prof.eggs, 0))

#chlorpyrifos reductions measured in mean eggs/snail/day  ###########
  plot(fn.hash$chlor.conc, fn.hash$chlor.rep / fn.hash$chlor.rep[1], ylim = c(0,1), 
       xlab = 'chlorpyrifos concentration (ppb)', ylab = 'relative decrease in fecundity', pch = 16)
    
ch.fn.red.hash = drm(chlor.rep ~ chlor.conc, data = fn.hash, type = 'continuous',
                 fct = L.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                           fixed = c(NA, 0, fn.hash$chlor.rep[1], NA)))

  summary(ch.fn.red.hash)
  
  fn.hash.chlor.pred = function(In){
    predict(ch.fn.red.hash, newdata = data.frame(chlor.conc = In), interval = 'confidence', level = 0.95)
  }  

  lines(seq(0,3000, 3), sapply(seq(0,3000, 3), fn.hash.chlor.pred)[1,] / fn.hash$chlor.rep[1], col = 2, lty=2)
  lines(seq(0,3000, 3), sapply(seq(0,3000, 3), fn.hash.chlor.pred)[2,] / fn.hash$chlor.rep[1], col = 2, lty=3)
  lines(seq(0,3000, 3), sapply(seq(0,3000, 3), fn.hash.chlor.pred)[3,] / fn.hash$chlor.rep[1], col = 2, lty=3)

fN.hash.chlor.uncertainty = function(In){
  if(In == 0) fn = 1 else{
    init = predict(ch.fn.red.hash, newdata = data.frame(chlor.conc = In), se.fit = T)
  fn = rnorm(1, init[1], init[2]) / fn.hash$chlor.rep[1]
  while(fn < 0 || fn > 1.00000){
    fn = rnorm(1, init[1], init[2]) / fn.hash$chlor.rep[1]
  }
}
  return(fn)
} #normalized to 1, upper limit at 1, lower limit at 0

points(seq(0,3000, 3), sapply(seq(0,3000, 3), fN.hash.chlor.uncertainty), pch = 5, cex = 0.5, col = 4)

#profenofos reductions measured in mean eggs/snail/day  ###########
plot(fn.hash$prof.conc, fn.hash$prof.rep / fn.hash$prof.rep[1], ylim = c(0,1), 
     xlab = 'profenofos concentration (ppb)', ylab = 'relative decrease in fecundity', pch = 16)

pr.fn.red.hash = drm(prof.rep ~ prof.conc, data = fn.hash, type = 'continuous',
                     fct = L.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                               fixed = c(NA, 0, fn.hash$prof.rep[1], NA)))

  summary(pr.fn.red.hash)

  fn.hash.prof.pred = function(In){
    predict(pr.fn.red.hash, newdata = data.frame(prof.conc = In), interval = 'confidence', level = 0.95)
  }  

lines(seq(0,4000, 4), sapply(seq(0,4000, 4), fn.hash.prof.pred)[1,] / fn.hash$prof.rep[1], col = 2, lty=2)
lines(seq(0,4000, 4), sapply(seq(0,4000, 4), fn.hash.prof.pred)[2,] / fn.hash$prof.rep[1], col = 2, lty=3)
lines(seq(0,4000, 4), sapply(seq(0,4000, 4), fn.hash.prof.pred)[3,] / fn.hash$prof.rep[1], col = 2, lty=3)

fN.hash.prof.uncertainty = function(In){
  if(In == 0) fn = 1 else{
    init = predict(pr.fn.red.hash, newdata = data.frame(prof.conc = In), se.fit = T)
    fn = rnorm(1, init[1], init[2]) / fn.hash$prof.rep[1]
    while(fn < 0 || fn > 1.00000){
      fn = rnorm(1, init[1], init[2]) / fn.hash$prof.rep[1]
    }
  }
  return(fn)
} #normalized to 1, upper limit at 1, lower limit at 0

points(seq(0,4000, 4), sapply(seq(0,4000, 4), fN.hash.prof.uncertainty), pch = 5, cex = 0.5, col = 4)

