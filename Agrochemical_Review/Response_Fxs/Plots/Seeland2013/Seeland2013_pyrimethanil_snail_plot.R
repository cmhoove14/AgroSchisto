#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Seeland et al 2013 pyrimethanil toxicity to snails#########

source("Agrochemical_Review/Response_Fxs/Seeland2013_pyrimethanil_snail_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Seeland2013/Seeland2013_pyrimethanil_snail_repro.png")

plot(Seel13_repro_dat$conc, 1-Seel13_repro_dat$mort, 
     xlim = c(0, 2), ylim = c(0,1), pch = 16,
     xlab = "Pyrimethanil (ppm)", ylab = "mortality",
     main = "Seeland et al 2013 Physella acuta reproduction")

  set.seed(43093)
  
  points(seq(0,2,0.01), sapply(seq(0,2,0.01)*1000, fNq_pyrimethanil_Seeland13_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("topright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  
  
png("Agrochemical_Review/Response_Fxs/Plots/Seeland2013/Seeland2013_pyrimethanil_snail_mort.png")

plot(Seel13_mort_dat$conc, Seel13_mort_dat$mort, 
     xlim = c(0, 2), ylim = c(0,1), pch = 16,
     xlab = "Pyrimethanil (ppm)", ylab = "mortality",
     main = "Seeland et al 2013 Physella acuta mortality")

  set.seed(43093)
  
  points(seq(0,2,0.01), sapply(seq(0,2,0.01)*1000, muNq_pyrimethanil_Seeland13_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("topleft", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  