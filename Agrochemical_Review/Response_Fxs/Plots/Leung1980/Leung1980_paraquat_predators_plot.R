#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Fornstrom 1997 data from Fig 1 at 24 hrs
source("Agrochemical_Review/Response_Fxs/Leung1980_paraquat_predators_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Leung1980/Leung1980_pred_mortality_data&sim.png")

plot(leung_dat$conc, leung_dat$mort, pch = 16, ylim = c(0,1), xlim = c(0,110),
     xlab = "Paraquat (ppm)", ylab = "mortality",
     main = "Leung 1980 Paraquat toxicity")

  set.seed(43093)
  
  points(c(0:120), sapply(c(0:120)*1000, muPq_paraquat_Leung_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("topleft", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)
  
dev.off()

png("Agrochemical_Review/Response_Fxs/Plots/Leung1980/Leung1980_juv_pred_mortality_data&sim.png")

plot(leung_dat_juv$conc, leung_dat_juv$mort, pch = 16, ylim = c(0,1), xlim = c(0,10),
     xlab = "Paraquat (ppm)", ylab = "mortality",
     main = "Leung 1980 Paraquat toxicity to juveniles")

  set.seed(43093)
  
  points(seq(0,10,0.1), sapply(seq(0,10,0.1)*1000, muPq_juv_paraquat_Leung_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("topleft", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)
  
dev.off()
