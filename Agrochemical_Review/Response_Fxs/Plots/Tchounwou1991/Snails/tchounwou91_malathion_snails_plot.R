#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/tchounwou91_malathion_snails_fit.R")
#direct mortality to snails #############
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Snails/tchounwou1991_malathion_snail_mortality.png")
  plot(sn$conc, sn$prop_dead, pch = 16, xlim = c(0, max(sn$conc)), ylim = c(0,1),
       xlab = 'malathion (ppm)', ylab = 'mortality')
    segments(x0 = 158.91, x1 = 244.16, y0 = 0.5, y1 = 0.5)
    
  set.seed(43093)  
    
    points(seq(0, 1e3, 5), sapply(seq(0, 1e6, 5e3), muNq_mal_tch91_uncertainty), 
           pch = 5, col = 4, cex = 0.5)
    
  legend('bottomright', pch = c(16,5), col = c(1,4), legend = c("observed", "simulated"),
         cex = 0.8, bty = 'n')
  
dev.off()  
#egg viability ###########
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Snails/tchounwou1991_malathion_snail_fecundity.png")
  
  plot(eg$conc, 1 - eg$prop_dead, pch = 16, xlim = c(0, max(eg$conc)), ylim = c(0,1),
       xlab = 'malathion (ppm)', ylab = 'egg viability')
    segments(x0 = 35.49, x1 = 143.97, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)  

  points(seq(0, 4e2, 1), sapply(seq(0, 4e5, 1e3), fNq_mal_tch91_uncertainty), 
         pch = 5, col = 4, cex = 0.5)

  legend('bottomleft', pch = c(16,5), col = c(1,4), legend = c("observed", "simulated"),
         cex = 0.8, bty = 'n')
  
dev.off()  