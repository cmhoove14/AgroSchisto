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
source("Agrochemical_Review/Response_Fxs/Ibrahim1992_chlorpyrifos_snails_fit.R")

f_N_chlor_ibr92 = function(In){
    predict(chlor.fN.predict, data.frame(dose = In), interval = 'confidence', level = 0.95)
}

png("Agrochemical_Review/Response_Fxs/Plots/Ibrahim1992/ibrahim1992_chlorpyrifos_snail_reproduction.png")

plot(snail.repro$dose, snail.repro$juvs.sn.day/snail.repro$juvs.sn.day[1], pch = 16, ylim = c(0,1),
     ylab = 'relative juveniles/snail', xlab = 'ChlorP ppb')
  
    lines(seq(0, 1000, 10), sapply(seq(0, 1000, 10), f_N_chlor_ibr92, simplify = T)[1,]/snail.repro$juvs.sn.day[1],
            lty = 2, col = 2)
    lines(seq(0, 1000, 10), sapply(seq(0, 1000, 10), f_N_chlor_ibr92, simplify = T)[2,]/snail.repro$juvs.sn.day[1],
            lty = 3, col = 2)
    lines(seq(0, 1000, 10), sapply(seq(0, 1000, 10), f_N_chlor_ibr92, simplify = T)[3,]/snail.repro$juvs.sn.day[1],
            lty = 3, col = 2)
    
  set.seed(43093)
  
    points(seq(0, 1000, 2), sapply(seq(0, 1000, 2), f_N_chlor_ibr92_uncertainty),
               pch = 5, cex = 0.5, col = 4)
    title(main = expression(paste('Ibrahim1992 - Chlorpyrifos reproductive toxicity to ',
                                  italic('Bi. alexandrina', sep = ''))))
      legend('topright', pch = c(16, 5), legend = c('Observed', 'Sampled'), col = c(1,4),
             cex = 0.7, bty = 'n')

dev.off()      
#Snail mortality #############
  #-- relative change in survival over entire experiment period (12 weeks)
  #Not advised to use this as mortality estimate as outcomes were only assessed at 5 weeks.
      #though authors note that all snails died within the first week @ 500ppb
  mu_N_chlor_ibr92 = function(In){
      predict(ibr_muNq, data.frame(dose = In), interval = 'confidence', level = 0.95)
  }

png("Agrochemical_Review/Response_Fxs/Plots/Ibrahim1992/ibrahim1992_chlorpyrifos_snail_mortality.png")
 
  plot(snail.mort$dose, snail.mort$mort, pch = 16, ylim = c(0,1),
       xlab = 'Chlorpyrifos (ppb)', ylab = 'mortality rate')  
  
    lines(seq(0, 1000, 10), 
          sapply(seq(0, 1000, 10), mu_N_chlor_ibr92, simplify = T)[1,],
          lty = 2, col = 2)
    lines(seq(0, 1000, 10), 
          sapply(seq(0, 1000, 10), mu_N_chlor_ibr92, simplify = T)[2,],
          lty = 3, col = 2)
    lines(seq(0, 1000, 10), 
          sapply(seq(0, 1000, 10), mu_N_chlor_ibr92, simplify = T)[3,],
          lty = 3, col = 2)
  
    points(seq(0, 1000, 2), sapply(seq(0, 1000, 2), mu_N_chlor_ibr92_uncertainty, simplify = T),
               pch = 5, cex = 0.5, col = 4)
    title(main = expression(paste('Ibrahim1992 - Chlorpyrifos direct snail toxicity ',
                                  italic('Bi. alexandrina', sep = ''))))
      legend('bottomright', pch = c(16, 5), legend = c('Observed', 'Sampled'), col = c(1,4), bty = 'n',
             cex = 0.7)
      
dev.off()