#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Naqvi 1983 P. clarkii toxicity
source("Agrochemical_Review/Response_Fxs/Naqvi1983_trifluralin_oryzalin_predators_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Naqvi/Naqvi1983_trifluralin_p_clarkii.png")

plot(naq83_trif_dat$conc, naq83_trif_dat$mort, pch = 16, ylim = c(0,1), xlim = c(0,25),
     xlab = "trifluralin (ppm)", ylab = "mortality", 
     main = "Naqvi 1983 trifluralin toxicity to P.clarkii")
  
  set.seed(43093)
  
  points(seq(0,25,0.25), sapply(seq(0,25,0.25)*1000, muPq_trifluralin_Naqvi83_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("topleft", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  
#oryzalin ##########
png("Agrochemical_Review/Response_Fxs/Plots/Naqvi/Naqvi1983_oryzalin_p_clarkii.png")

plot(naq83_oryz_dat$conc, naq83_oryz_dat$mort, pch = 16, ylim = c(0,1), xlim = c(0,25000),
     xlab = "oryzalin (ppm)", ylab = "mortality", 
     main = "Naqvi 1983 oryzalin toxicity to P.clarkii")
  
  set.seed(43093)
  
  points(seq(0,25000,250), sapply(seq(0,25000,250)*1000, muPq_oryzalin_Naqvi83_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("topleft", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)


dev.off()