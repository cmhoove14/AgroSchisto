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
mun.ch = data.frame(conc = c(.72, 1.32, 2.82),
                     mort = c(.25, .50 , .90),
                     surv = 0)
mun.ch$surv = 1 - mun.ch$mort
 
plot(mun.ch$conc*1000, mun.ch$mort, pch = 16, cex = 1.2, ylim = c(0,1), xlim = c(0,3500),
     xlab = 'Chlorpyrifos (ppm)', ylab = 'prop dead')
  segments(x0 = 880, y0 = 0.5, x1 = 1980, y1 = 0.5)
  
lc50.chlor.hash = 1.32
  se.lc50.chlor.hash = mean(log(1.98/lc50.chlor.hash), log(lc50.chlor.hash/0.88)) / 1.96
slp.chlor.hash = 2.5

#function based on provided values
  fx.mun.chlor = function(In, lc = lc50.chlor.hash){
    ins = In/1000
    pnorm(slp.chlor.hash * log(ins/lc))
  }
  
    lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.mun.chlor), lty = 2, col = 2)
    lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.mun.chlor, lc = 1.98), lty = 3, col = 2)
    lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.mun.chlor, lc = 0.88), lty = 3, col = 2)

mu_Nq_ch_Hash11_uncertainty = function(In){
  ins = In/1000
  lc50 = exp(rnorm(1, log(lc50.chlor.hash), se.lc50.chlor.hash))
  pnorm(slp.chlor.hash * log(ins/lc50))
}
points(seq(0,4000,10), sapply(seq(0,4000,10), mu_Nq_ch_Hash11_uncertainty, simplify = T), 
       pch = 5, col = 4, cex = 0.5)

#Profenofos #########
mun.prof = data.frame(conc = c(1.4, 2.5, 3.72),
                    mort = c(.25, .50 , .90),
                    surv = 0)
  mun.prof$surv = 1 - mun.prof$mort
    
  plot(mun.prof$conc*1000, mun.prof$mort, pch = 16, cex = 1.2, ylim = c(0,1), xlim = c(0,5000),
       xlab = 'Profenofos (ppm)', ylab = 'prop dead')
      segments(x0 = 1880, y0 = 0.5, x1 = 3330, y1 = 0.5)

lc50.prof.hash = 2.5
  se.lc50.prof.hash = mean(log(3.33/lc50.prof.hash), log(lc50.prof.hash/1.88)) / 1.96
slp.prof.hash = 1.6

#function based on provided values
fx.mun.prof = function(In, lc = lc50.prof.hash){
  ins = In/1000
  pnorm(slp.prof.hash * log(ins/lc))
}

  lines(seq(0,5000,10), sapply(seq(0,5000,10), fx.mun.prof), lty = 2, col = 2)
  lines(seq(0,5000,10), sapply(seq(0,5000,10), fx.mun.prof, lc = 3.33), lty = 3, col = 2)
  lines(seq(0,5000,10), sapply(seq(0,5000,10), fx.mun.prof, lc = 1.88), lty = 3, col = 2)

mu_Nq_prof_Hash11_uncertainty = function(In){
  ins = In/1000
  lc50 = exp(rnorm(1, log(lc50.prof.hash), se.lc50.prof.hash))
  pnorm(slp.prof.hash * log(ins/lc50))
}
points(seq(0,5000,10), sapply(seq(0,5000,10), mu_Nq_prof_Hash11_uncertainty, simplify = T), 
       pch = 5, col = 4, cex = 0.5)


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
  while(fn < 0 && fn > 1.00000){
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
    while(fn < 0 && fn > 1.00000){
      fn = rnorm(1, init[1], init[2]) / fn.hash$prof.rep[1]
    }
  }
  return(fn)
} #normalized to 1, upper limit at 1, lower limit at 0

points(seq(0,4000, 4), sapply(seq(0,4000, 4), fN.hash.prof.uncertainty), pch = 5, cex = 0.5, col = 4)

