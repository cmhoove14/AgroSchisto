#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/Revathi2010_tributyltin_predator_reproduction_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Revathi2010/revathi2010_tributyltin_pred_repro.png")

plot(revathi_dat$tbt, revathi_dat$hatch, pch = 16, xlim = c(0, 4), ylim = c(0,1),
     xlab = "tributyltin (ppm)", ylab = "viable hatching proportion",
     main = "Revathi 2010 M. rosenbergii")

  lines(seq(0,4,0.05), predict(revathi_mod, newdata = data.frame(tbt = seq(0,4,0.05))),
        lty = 2, col = 2)
  
  set.seed(43093)

  points(seq(0,4,0.05), sapply(seq(0,4,0.05)*1000, fPq_tbt_Revathi10_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomleft", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()