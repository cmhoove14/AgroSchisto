#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Data extraction and model fitting to Naqvi et al 1987 P. clarkii toxicity
source("Agrochemical_Review/Response_Fxs/Naqvi1987_endosulfan_etc_predators_fit.R")

#Thiodan (endosulfan) #####################
png("Agrochemical_Review/Response_Fxs/Plots/Naqvi/Naqvi1987_endosulfan_adult_p_clarkii.png")

plot(c(0, lc50.naq.endo.report), c(0,0.5), pch = 16, xlim = c(0, lc50.naq.endo.report*2), ylim = c(0,1),
     xlab = "endosulfan (ppb)", ylab = "mortality", main = "Endosulfan toxicity to P. clarkii, Naqvi 1987")
  segments(x0 = 356.3, x1 = 503.9, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,1000,5), sapply(seq(0,1000,5), muPq_endo_Naqvi87_uncertainty), pch = 5, col = 4, cex = 0.5)
  
  legend("topleft", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()

#Treflan (trifluaralin) #############    
png("Agrochemical_Review/Response_Fxs/Plots/Naqvi/Naqvi1987_treflan_adult_p_clarkii.png")

plot(c(0, lc50.naq.trif.report), c(0,0.5), pch = 16, xlim = c(0, lc50.naq.trif.report*2), ylim = c(0,1),
     xlab = "treflan (ppm)", ylab = "mortality", main = "Treflan toxicity to P. clarkii, Naqvi 1987")
  segments(x0 = 28.9, x1 = 23.8, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,50,0.5), sapply(seq(0,50,0.5)*1000, muPq_trif_Naqvi87_uncertainty), pch = 5, col = 4, cex = 0.5)
  
  legend("topleft", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()

#MSMA (Monosodium methanearsonate) #############    
png("Agrochemical_Review/Response_Fxs/Plots/Naqvi/Naqvi1987_msma_adult_p_clarkii.png")

plot(c(0, lc50.naq.msma.report), c(0,0.5), pch = 16, xlim = c(0, lc50.naq.msma.report*2), ylim = c(0,1),
     xlab = "MSMA (ppm)", ylab = "mortality", main = "MSMA toxicity to P. clarkii, Naqvi 1987")
  segments(x0 = 916.8, x1 = 1123.8, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,2000,20), sapply(seq(0,2000,20)*1000, muPq_msma_Naqvi87_uncertainty), pch = 5, col = 4, cex = 0.5)
  
  legend("topleft", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()
