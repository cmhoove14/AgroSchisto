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
# reduction in snail fecundity: from table 1, relative change in number of juveniles produced over the 8 week experiment
snail.repro = data.frame(dose = c(0,125,250,500)*1000, #chlorpyrifos in ppb
                         juvs = c(4225, 3005, 1749, 1313),
                        #convert 8 week juvenile cohort to juveniles per snail per day
                         juvs.sn.day = c(4225, 3005, 1749, 1313) / 30 / (7*8)) 
  snail.repro$rel.juvs.sn.day = snail.repro$juvs.sn.day / snail.repro$juvs.sn.day[1]

chlor.fN.predict = drm(juvs ~ dose, data = snail.repro, type = 'continuous',
                       fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                                  fixed = c(NA, 0, max(snail.repro$juvs), NA)))
  summary(chlor.fN.predict)
  
plot(snail.repro$dose, snail.repro$juvs/snail.repro$juvs[1], pch = 16, ylim = c(0,1),
     ylab = 'relative juveniles /day over 8 weeks', xlab = 'ChlorP ppb')
  
  f_N_chlor_ibr92 = function(In){
      predict(chlor.fN.predict, data.frame(dose = In), interval = 'confidence', level = 0.95)
  }

    lines(seq(0, 5e5, 100), sapply(seq(0, 5e5, 100), f_N_chlor_ibr92, simplify = T)[1,]/snail.repro$juvs[1],
            lty = 2, col = 2)
    lines(seq(0, 5e5, 100), sapply(seq(0, 5e5, 100), f_N_chlor_ibr92, simplify = T)[2,]/snail.repro$juvs[1],
            lty = 3, col = 2)
    lines(seq(0, 5e5, 100), sapply(seq(0, 5e5, 100), f_N_chlor_ibr92, simplify = T)[3,]/snail.repro$juvs[1],
            lty = 3, col = 2)
    
par.tricks = c(coef(chlor.fN.predict), 'Upper Limit:(Intercept)' = 4.245811e+03)[c(1,3,2)]

  f_N_chlor_ibr92_uncertainty<-function(In){
    rdrm(nosim = 1, fct = LL.3(), mpar = par.tricks, yerror = 'rnorm', xerror = In,
         ypar = c(0, predict(chlor.fN.predict, data.frame(dose = In), se.fit = T)[2]))$y / snail.repro$juvs[1]
  }
  
    points(seq(0, 5e5, 1000), sapply(seq(0, 5e5, 1000), f_N_chlor_ibr92_uncertainty, simplify = T),
               pch = 5, cex = 0.6, col = 4)
    title(main = expression(paste('Ibrahim1992 - Chlorpyrifos reproductive toxicity to ',
                                  italic('Bi. alexandrina', sep = ''))))
      legend('topright', pch = c(16, 5), legend = c('Obs. points', 'Est. points'), col = c(1,4),
             cex = 0.7)

#zoom to more environmentally feasible range
  plot(snail.repro$dose, snail.repro$juvs/snail.repro$juvs[1], pch = 16, ylim = c(0,1), xlim = c(0,200),
     ylab = 'relative juveniles /day over 8 weeks', xlab = 'ChlorP ppb')
  
  lines(c(0:200), sapply(c(0:200), f_N_chlor_ibr92, simplify = T)[1,]/snail.repro$juvs[1],
            lty = 2, col = 2)
    lines(c(0:200), sapply(c(0:200), f_N_chlor_ibr92, simplify = T)[2,]/snail.repro$juvs[1],
            lty = 3, col = 2)
    lines(c(0:200), sapply(c(0:200), f_N_chlor_ibr92, simplify = T)[3,]/snail.repro$juvs[1],
            lty = 3, col = 2)
    
    points(c(0:200), sapply(c(0:200), f_N_chlor_ibr92_uncertainty, simplify = T),
               pch = 5, cex = 0.6, col = 4)
    title(main = expression(paste('Ibrahim1992 - Chlorpyrifos reproductive toxicity to ',
                                  italic('Bi. alexandrina', sep = ''))))
      legend('topright', pch = c(16, 5), legend = c('Obs. points', 'Est. points'), col = c(1,4),
             cex = 0.7)
    
    
#Snail mortality #############
  #-- relative change in survival over entire experiment period (12 weeks)
  snail.mort = data.frame(dose = c(0,125,250,500)*1000,
                          total = rep(30,4),
                          dead = c(5,7,8,30))
  
  snail.mort$mort = snail.mort$dead / snail.mort$total
  snail.mort$mean.daily.rate = snail.mort$mort / 84
  
  
  
  ibr_muNq<-drm(dead/total ~ dose, weights = total, data = snail.mort,
                type = 'binomial', fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                                              fixed = c(NA, 0, 1, NA)))
  
    summary(ibr_muNq)

  plot(snail.mort$dose, snail.mort$mort, pch = 16, ylim = c(0,1))  
  
  mu_N_chlor_ibr92 = function(In){
      predict(ibr_muNq, data.frame(dose = In), interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0, 5e5, 100), sapply(seq(0, 5e5, 100), mu_N_chlor_ibr92, simplify = T)[1,],
            lty = 2, col = 2)
    lines(seq(0, 5e5, 100), sapply(seq(0, 5e5, 100), mu_N_chlor_ibr92, simplify = T)[2,],
            lty = 3, col = 2)
    lines(seq(0, 5e5, 100), sapply(seq(0, 5e5, 100), mu_N_chlor_ibr92, simplify = T)[3,],
            lty = 3, col = 2)
  
  mu_N_chlor_ibr92_uncertainty<-function(In){
    rdrm(nosim = 1, fct = LL.2(), mpar = coef(ibr_muNq), yerror = 'rbinom', xerror = In,
         ypar = 30)$y / 30
  }
  
    points(seq(0, 5e5, 1000), sapply(seq(0, 5e5, 1000), mu_N_chlor_ibr92_uncertainty, simplify = T),
               pch = 5, cex = 0.6, col = 4)
    title(main = expression(paste('Ibrahim1992 - Chlorpyrifos direct snail toxicity ',
                                  italic('Bi. alexandrina', sep = ''))))
      legend('bottomright', pch = c(16, 5), legend = c('Obs. points', 'Est. points'), col = c(1,4),
             cex = 0.7)