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
mun.hsh.ch = data.frame(conc = c(0,.72, 1.32, 2.82)*1000,
                        mort = c(0,.25, .50 , .90),
                        surv = 0)
  mun.hsh.ch$surv = 1 - mun.hsh.ch$mort
  mun.hsh.ch$probit = qnorm(mun.hsh.ch$mort, mean = 5)
  mun.hsh.ch$ppm = mun.hsh.ch$conc/1000
  mun.hsh.ch$log = log(mun.hsh.ch$ppm)
  
  lc50.ch.hash.report = 1.32
  slp.ch.hash.report = 2.5
  se.lc50.ch.hash = mean(c(log(1.98/lc50.ch.hash.report) , log(lc50.ch.hash.report/0.88))) / 1.96

  muNq_ch_hash11_uncertainty = function(In){
    #if(In == 0) mun = 0 else{
    Ins = (In/1000)
    lc50 = exp(rnorm(1, log(lc50.ch.hash.report), se.lc50.ch.hash))
    mun = pnorm((slp.ch.hash.report) * log(Ins/lc50)) 
    #}
    
    return(mun)
  }
  
  plot(mun.hsh.ch$conc, mun.hsh.ch$mort, pch = 16, ylim = c(-0.1,1), xlim = c(0,3500),
       xlab = 'Chlorpyrifos (ppb)', ylab = 'prop dead')
    segments(x0 = 880, y0 = 0.5, x1 = 1980, y1 = 0.5)
  
  points(seq(0,3500,10), sapply(seq(0,3500,10), muNq_ch_hash11_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
 
#keep vector
keep.hsh.ch = c('muNq_ch_hash11_uncertainty', 
                 'lc50.ch.hash.report', 'se.lc50.ch.hash', 'slp.ch.hash.report')    

#Profenofos #########
mun.hsh.prof = data.frame(conc = c(0,1.4, 2.5, 3.72)*1000,
                    mort = c(0,.25, .50 , .90),
                    surv = 0)
  mun.hsh.prof$surv = 1 - mun.hsh.prof$mort
  mun.hsh.prof$ppm = mun.hsh.prof$conc/1000
  mun.hsh.prof$log10 = log10(mun.hsh.prof$ppm)
  mun.hsh.prof$probit = qnorm(mun.hsh.prof$mort, mean = 5)
    
  lc50.prof.hash.report = 2.5
  slp.prof.hash.report = 1.6
  se.lc50.prof.hash = mean(c(log(3.33/lc50.prof.hash.report) , log(lc50.prof.hash.report/1.88))) / 1.96
  
  muNq_prof_hash11_uncertainty = function(In){
    #if(In == 0) mun = 0 else{
    Ins = (In/1000)
    lc50 = exp(rnorm(1, log(lc50.prof.hash.report), se.lc50.prof.hash))
    mun = pnorm((slp.prof.hash.report) * log(Ins/lc50))
    #}

    return(mun)
  }
  
  plot(mun.hsh.prof$conc, mun.hsh.prof$mort, pch = 16, ylim = c(0,1), xlim = c(0,5000),
       xlab = 'Profenofos (ppm)', ylab = 'prop dead')
    segments(x0 = 1880, y0 = 0.5, x1 = 3330, y1 = 0.5)
  
  points(seq(0,5000,10), sapply(seq(0,5000,10), muNq_prof_hash11_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
  
#keep vector
keep.hsh.prof = c('muNq_prof_hash11_uncertainty',
                'lc50.prof.hash.report', 'se.lc50.prof.hash', 'slp.prof.hash.report')    

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
    ctrl.eggs.sd = sd(fn$eggs_snail_day[fn$chem == 'control' & fn$surv != 0])
  
  chlor.eggs = mean(fn$eggs_snail_day[fn$chem == 'chlorpyrifos' & fn$surv != 0])
    chlor.eggs.sd = sd(fn$eggs_snail_day[fn$chem == 'chlorpyrifos' & fn$surv != 0])
  
  prof.eggs = mean(fn$eggs_snail_day[fn$chem == 'profenofos' & fn$surv != 0])
    prof.eggs.sd = sd(fn$eggs_snail_day[fn$chem == 'profenofos' & fn$surv != 0])
  
#Functions that estimate parameter value in single concentration group #########
    fN.hash.chlor.uncertainty = function(waste){ #function based on snail egg masses
      wst = waste
      fN = rnorm(1, chlor.eggs, chlor.eggs.sd) / ctrl.eggs
      while(fN < 0) fN = rnorm(1, chlor.eggs, chlor.eggs.sd) / ctrl.eggs
      fN
    }
    
    #hist(sapply(rep(1,10000), fN.hash.chlor.uncertainty, simplify = T))
    
    fN.hash.prof.uncertainty = function(waste){ #function based on snail hatchlings
      wst = waste
      fN = rnorm(1, prof.eggs, prof.eggs.sd) / ctrl.eggs
      while(fN < 0) fN = rnorm(1, prof.eggs, prof.eggs.sd) / ctrl.eggs
      fN
    }  
    
    #hist(sapply(rep(1,10000), fN.hash.prof.uncertainty, simplify = T))
    
    
#Assume lc90 halts reproduction for third data point ###########   
fn.hash = data.frame(chlor.conc = c(0, 720, 2820),
                     prof.conc = c(0, 1400, 3720),
                     chlor.rep = c(ctrl.eggs, chlor.eggs, 0),
                     prof.rep = c(ctrl.eggs, prof.eggs, 0))

#chlorpyrifos reductions measured in mean eggs/snail/day  ###########
  plot(fn.hash$chlor.conc, fn.hash$chlor.rep / fn.hash$chlor.rep[1], ylim = c(0,1), 
       xlab = 'chlorpyrifos concentration (ppb)', ylab = 'relative decrease in fecundity', pch = 16)
    
ch.fn.red.hash = drm(chlor.rep ~ chlor.conc, data = fn.hash, type = 'continuous',
                 fct = LL.3(names = c("b", "d", "e"),
                           fixed = c(NA, fn.hash$chlor.rep[1], NA)))

  summary(ch.fn.red.hash)
  
  fn.hash.chlor.pred = function(In){
    predict(ch.fn.red.hash, newdata = data.frame(chlor.conc = In), interval = 'confidence', level = 0.95)
  }  

  lines(seq(0,3000, 3), sapply(seq(0,3000, 3), fn.hash.chlor.pred)[1,] / fn.hash$chlor.rep[1], col = 2, lty=2)
  lines(seq(0,3000, 3), sapply(seq(0,3000, 3), fn.hash.chlor.pred)[2,] / fn.hash$chlor.rep[1], col = 2, lty=3)
  lines(seq(0,3000, 3), sapply(seq(0,3000, 3), fn.hash.chlor.pred)[3,] / fn.hash$chlor.rep[1], col = 2, lty=3)

fN.hash.chlor.uncertainty3pts = function(In){
  if(In == 0) fn = 1 else{
    init = predict(ch.fn.red.hash, newdata = data.frame(chlor.conc = In), se.fit = T)
  fn = rnorm(1, init[1], init[2]) / fn.hash$chlor.rep[1]
}
  if(fn < 0) fn = 0 
  if(fn > 1) fn = 1
  
  return(fn)
} #normalized to 1, upper limit at 1, lower limit at 0

points(seq(0,3000, 3), sapply(seq(0,3000, 3), fN.hash.chlor.uncertainty3pts), pch = 5, cex = 0.5, col = 4)

keep.hsh.ch = c(keep.hsh.ch, 'fN.hash.chlor.uncertainty3pts', 'ch.fn.red.hash', 'fn.hash', 'fN.hash.chlor.uncertainty',
                'chlor.eggs', 'chlor.eggs.sd', 'ctrl.eggs')
#profenofos reductions measured in mean eggs/snail/day  ###########
plot(fn.hash$prof.conc, fn.hash$prof.rep / fn.hash$prof.rep[1], ylim = c(0,1), 
     xlab = 'profenofos concentration (ppb)', ylab = 'relative decrease in fecundity', pch = 16)

pr.fn.red.hash = drm(prof.rep ~ prof.conc, data = fn.hash, type = 'continuous',
                     fct = LL.3(names = c("b", "d", "e"),
                               fixed = c(NA, fn.hash$prof.rep[1], NA)))

  summary(pr.fn.red.hash)

  fn.hash.prof.pred = function(In){
    predict(pr.fn.red.hash, newdata = data.frame(prof.conc = In), interval = 'confidence', level = 0.95)
  }  

lines(seq(0,4000, 4), sapply(seq(0,4000, 4), fn.hash.prof.pred)[1,] / fn.hash$prof.rep[1], col = 2, lty=2)
lines(seq(0,4000, 4), sapply(seq(0,4000, 4), fn.hash.prof.pred)[2,] / fn.hash$prof.rep[1], col = 2, lty=3)
lines(seq(0,4000, 4), sapply(seq(0,4000, 4), fn.hash.prof.pred)[3,] / fn.hash$prof.rep[1], col = 2, lty=3)

fN.hash.prof.uncertainty3pts = function(In){
  if(In == 0) fn = 1 else{
    init = predict(pr.fn.red.hash, newdata = data.frame(prof.conc = In), se.fit = T)
    fn = rnorm(1, init[1], init[2]) / fn.hash$prof.rep[1]
  }
    if(fn < 0) fn = 0 
    if(fn > 1) fn = 1
  
  return(fn)
} #normalized to 1, upper limit at 1, lower limit at 0

points(seq(0,4000, 4), sapply(seq(0,4000, 4), fN.hash.prof.uncertainty3pts), pch = 5, cex = 0.5, col = 4)

keep.hsh.prof = c(keep.hsh.prof, 'fN.hash.prof.uncertainty3pts', 'pr.fn.red.hash', 'fn.hash', 'fN.hash.prof.uncertainty',
                  'prof.eggs', 'prof.eggs.sd', 'ctrl.eggs')
keep.hsh.all = c(keep.hsh.ch, keep.hsh.prof)