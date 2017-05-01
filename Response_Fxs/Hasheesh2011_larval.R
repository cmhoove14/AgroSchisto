#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
    
#Toxicity to miracidia and cercariae from table 5 ###############
#ChlorP miracidia #########
piM.ch = data.frame(conc = c(0.78, 1.8),
                    mort = c(.50 , .90),
                    surv = 0)
piM.ch$surv = 1 - piM.ch$mort

  plot(piM.ch$conc*1000, piM.ch$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
       xlab = 'Chlorpyrifos (ppb)', ylab = '8-hr miracidial mortality',
       main = expression(paste(pi[M], ' estimate from Hasheesh 2011')))
    segments(x0 = 560, y0 = 0.5, x1 = 1100, y1 = 0.5)

lc50.chlor.hash.mir = 0.78
  se.lc50.chlor.hash.mir = mean(log(1.1/lc50.chlor.hash.mir), log(lc50.chlor.hash.mir/0.56)) / 1.96
slp.chlor.hash.mir = 1.86

#function based on provided values
fx.piM.chlor = function(In, lc = lc50.chlor.hash.mir){
  ins = In/1000
  1-pnorm(slp.chlor.hash.mir * log(ins/lc))
}

  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piM.chlor), lty = 2, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piM.chlor, lc = 1.1), lty = 3, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piM.chlor, lc = 0.56), lty = 3, col = 2)

piM_ch_Hash11_uncertainty = function(In){
  ins = In/1000
  lc50 = exp(rnorm(1, log(lc50.chlor.hash.mir), se.lc50.chlor.hash.mir))
  1 - pnorm(slp.chlor.hash.mir * log(ins/lc50))
}
  points(seq(0,4000,4), sapply(seq(0,4000,4), piM_ch_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)

#Profenofos miracidia #########
piM.pr = data.frame(conc = c(1.5, 2.51),
                    mort = c(.50 , .90),
                    surv = 0)
  piM.pr$surv = 1 - piM.pr$mort
  
plot(piM.pr$conc*1000, piM.pr$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
     xlab = 'Profenofos (ppb)', ylab = '8-hr miracidial mortality',
     main = expression(paste(pi[M], ' estimate from Hasheesh 2011')))
  segments(x0 = 1150, y0 = 0.5, x1 = 1950, y1 = 0.5)
  
lc50.prof.hash.mir = 1.5
  se.lc50.prof.hash.mir = mean(log(1.95/lc50.prof.hash.mir), log(lc50.prof.hash.mir/1.15)) / 1.96
slp.prof.hash.mir = 1.64
  
#function based on provided values
fx.piM.prof = function(In, lc = lc50.prof.hash.mir){
  ins = In/1000
  1-pnorm(slp.prof.hash.mir * log(ins/lc))
}
  
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piM.prof), lty = 2, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piM.prof, lc = 1.95), lty = 3, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piM.prof, lc = 1.15), lty = 3, col = 2)
  
piM_pr_Hash11_uncertainty = function(In){
  ins = In/1000
  lc50 = exp(rnorm(1, log(lc50.prof.hash.mir), se.lc50.prof.hash.mir))
  1 - pnorm(slp.prof.hash.mir * log(ins/lc50))
}

  points(seq(0,4000,4), sapply(seq(0,4000,4), piM_pr_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)
  
#profenofos cercariae #######
piC.pr = data.frame(conc = c(1.85, 2.85),
                    mort = c(.50 , .90),
                    surv = 0)
  piC.pr$surv = 1 - piC.pr$mort
  
plot(piC.pr$conc*1000, piC.pr$surv, pch = 16, ylim = c(0,1), xlim = c(0,5000),
     xlab = 'Profenofos (ppb)', ylab = '8-hr cercarial mortality',
     main = expression(paste(pi[C], ' estimate from Hasheesh 2011')))
  segments(x0 = 1320, y0 = 0.5, x1 = 2590, y1 = 0.5)
  
lc50.prof.hash.cerc = 1.85
  se.lc50.prof.hash.cerc = mean(log(2.59/lc50.prof.hash.cerc), log(lc50.prof.hash.cerc/1.32)) / 1.96
slp.prof.hash.cerc = 1.84
  
#function based on provided values
fx.piC.prof = function(In, lc = lc50.prof.hash.cerc){
  ins = In/1000
  1-pnorm(slp.prof.hash.cerc * log(ins/lc))
}
  
  lines(seq(0,5000,10), sapply(seq(0,5000,10), fx.piC.prof), lty = 2, col = 2)
  lines(seq(0,5000,10), sapply(seq(0,5000,10), fx.piC.prof, lc = 2.59), lty = 3, col = 2)
  lines(seq(0,5000,10), sapply(seq(0,5000,10), fx.piC.prof, lc = 1.32), lty = 3, col = 2)
  
  piC_pr_Hash11_uncertainty = function(In){
    ins = In/1000
    lc50 = exp(rnorm(1, log(lc50.prof.hash.cerc), se.lc50.prof.hash.cerc))
    1 - pnorm(slp.prof.hash.cerc * log(ins/lc50))
  }
  
  points(seq(0,5000,5), sapply(seq(0,5000,5), piC_pr_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)
  
#chlorpyrifos cercariae #######
piC.ch = data.frame(conc = c(0.96, 2.1),
                    mort = c(.50 , .90),
                    surv = 0)
  piC.ch$surv = 1 - piC.ch$mort

plot(piC.ch$conc*1000, piC.ch$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
     xlab = 'Chlorpyrifos (ppb)', ylab = '8-hr cercarial mortality',
     main = expression(paste(pi[C], ' estimate from Hasheesh 2011')))
  segments(x0 = 620, y0 = 0.5, x1 = 1440, y1 = 0.5)

lc50.chlor.hash.cerc = 0.96
  se.lc50.chlor.hash.cerc = mean(log(1.44/lc50.chlor.hash.cerc), log(lc50.chlor.hash.cerc/0.62)) / 1.96
slp.chlor.hash.cerc = 2.3

#function based on provided values
fx.piC.chlor = function(In, lc = lc50.chlor.hash.cerc){
  ins = In/1000
  1-pnorm(slp.chlor.hash.cerc * log(ins/lc))
}

  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piC.chlor), lty = 2, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piC.chlor, lc = 1.44), lty = 3, col = 2)
  lines(seq(0,4000,10), sapply(seq(0,4000,10), fx.piC.chlor, lc = 0.62), lty = 3, col = 2)

piC_ch_Hash11_uncertainty = function(In){
  ins = In/1000
  lc50 = exp(rnorm(1, log(lc50.chlor.hash.cerc), se.lc50.chlor.hash.cerc))
  1 - pnorm(slp.chlor.hash.cerc * log(ins/lc50))
}

  points(seq(0,4000,4), sapply(seq(0,4000,4), piC_ch_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)
