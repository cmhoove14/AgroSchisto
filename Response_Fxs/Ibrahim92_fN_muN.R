#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Toxicity to Biomphalaria snails from Ibrahim 1992 ###################
require(drc)
require(LW1949)

# reduction in snail fecundity: from table 1, relative change in number of juveniles ####
#produced over the 8 week experiment
snail.repro = data.frame(dose = c(0,125,250,500), #chlorpyrifos in ppb
                         juvs = c(4225, 3005, 1749, 1313),
                        #convert 12 week juvenile cohort to juveniles per snail per day
                         juvs.sn.day = c(4225, 3005, 1749, 1313) / 30 / (7*12)) 
  snail.repro$rel.juvs.sn.day = snail.repro$juvs.sn.day / snail.repro$juvs.sn.day[1]

#model with juveniles/snail/day as outcome  
chlor.fN.predict = drm(juvs.sn.day ~ dose, data = snail.repro, type = 'continuous',
                       fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                                  fixed = c(NA, 0, max(snail.repro$juvs.sn.day), NA)))
  summary(chlor.fN.predict)
  
plot(snail.repro$dose, snail.repro$juvs.sn.day/snail.repro$juvs.sn.day[1], pch = 16, ylim = c(0,1),
     ylab = 'relative juveniles/snail', xlab = 'ChlorP ppb')
  
f_N_chlor_ibr92 = function(In){
    predict(chlor.fN.predict, data.frame(dose = In), interval = 'confidence', level = 0.95)
}

    lines(seq(0, 1000, 10), sapply(seq(0, 1000, 10), f_N_chlor_ibr92, simplify = T)[1,]/snail.repro$juvs.sn.day[1],
            lty = 2, col = 2)
    lines(seq(0, 1000, 10), sapply(seq(0, 1000, 10), f_N_chlor_ibr92, simplify = T)[2,]/snail.repro$juvs.sn.day[1],
            lty = 3, col = 2)
    lines(seq(0, 1000, 10), sapply(seq(0, 1000, 10), f_N_chlor_ibr92, simplify = T)[3,]/snail.repro$juvs.sn.day[1],
            lty = 3, col = 2)
    
  f_N_chlor_ibr92_uncertainty<-function(In){
    if(In == 0) fn = 1 else{
      init = predict(chlor.fN.predict, newdata = data.frame(mal = In), se.fit = T)
      fn = rnorm(1, init[1], init[2]) / snail.repro$juvs.sn.day[1]
      if(fn < 0) fn = 0
      if(fn > 1) fn = 1
    }
    return(fn)
  }
  
    points(seq(0, 1000, 2), sapply(seq(0, 1000, 2), f_N_chlor_ibr92_uncertainty),
               pch = 5, cex = 0.6, col = 4)
    title(main = expression(paste('Ibrahim1992 - Chlorpyrifos reproductive toxicity to ',
                                  italic('Bi. alexandrina', sep = ''))))
      legend('topright', pch = c(16, 5), legend = c('Obs. points', 'Est. points'), col = c(1,4),
             cex = 0.7)

keep.ibr.fn.ch = c('snail.repro', 'chlor.fN.predict', 'f_N_chlor_ibr92_uncertainty')    
#Snail mortality #############
  #-- relative change in survival over entire experiment period (12 weeks)
  #Not advised to use this as mortality estimate as outcomes were only assessed at 5 weeks.
      #though authors note that all snails died within the first week @ 500ppb
  snail.mort = data.frame(dose = c(0,125,250,500),
                          total = rep(30,4),
                          dead = c(5,7,8,30))
  
  snail.mort$mort = snail.mort$dead / snail.mort$total
  snail.mort$mean.daily.rate = snail.mort$mort / 84
  
  
  ibr_muNq<-drm(dead/total ~ dose, weights = total, data = snail.mort,
                type = 'binomial', fct = L.4(names = c('b', 'd', 'c', 'e'),
                                              fixed = c(NA, snail.mort$mort[1], 1, NA)))
  
    summary(ibr_muNq)

  plot(snail.mort$dose, snail.mort$mort - snail.mort$mort[1], pch = 16, ylim = c(0,1),
       xlab = 'Chlorpyrifos (ppb)', ylab = 'relative mortality')  
  
  mu_N_chlor_ibr92 = function(In){
      predict(ibr_muNq, data.frame(dose = In), interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0, 1000, 10), 
          sapply(seq(0, 1000, 10), mu_N_chlor_ibr92, simplify = T)[1,] - snail.mort$mort[1],
          lty = 2, col = 2)
    lines(seq(0, 1000, 10), 
          sapply(seq(0, 1000, 10), mu_N_chlor_ibr92, simplify = T)[2,] - snail.mort$mort[1],
          lty = 3, col = 2)
    lines(seq(0, 1000, 10), 
          sapply(seq(0, 1000, 10), mu_N_chlor_ibr92, simplify = T)[3,] - snail.mort$mort[1],
          lty = 3, col = 2)
  
par.tricks.ibr.muN = c(coef(ibr_muNq), 
                       'Upper Limit:(Intercept)' = 1, 
                       'Lower Limit:(Intercept)' = snail.mort$mort[1])[c(1,4,3,2)]  
    
  mu_N_chlor_ibr92_uncertainty<-function(In){
    if(In == 0) mun = 0 else{
      mun = rdrm(nosim = 1, fct = L.4(), mpar = par.tricks.ibr.muN, yerror = 'rbinom', xerror = In,
               ypar = 30)$y / 30 - snail.mort$mort[1]
      if(mun < 0) mun = 0
    }
    return(mun)
  }
  
    points(seq(0, 1000, 2), sapply(seq(0, 1000, 2), mu_N_chlor_ibr92_uncertainty, simplify = T),
               pch = 5, cex = 0.5, col = 4)
    title(main = expression(paste('Ibrahim1992 - Chlorpyrifos direct snail toxicity ',
                                  italic('Bi. alexandrina', sep = ''))))
      legend('bottomright', pch = c(16, 5), legend = c('Obs. points', 'Est. points'), col = c(1,4),
             cex = 0.7)