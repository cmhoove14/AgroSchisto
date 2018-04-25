#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/tchounwou91b_fertilizer_egg_viability_fit.R")

#Amm. sulphate model ########## 
tch91_amm_pred <- predict(tch91.egv.amm, newdata = data.frame(conc = c(0:10000)), 
                          type = "response", interval = "confidence", level = 0.95)

png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Eggs/amm_sulphate_egg_viability_data_sim.png")

plot(eggv.amm$conc, eggv.amm$per_hatch, pch = 16, 
     xlab = "Amm. P concentration (ppm)", ylab = "% hatched")
  lines(c(0:10000), tch91_amm_pred[,1], lty = 2, col = 2)
  lines(c(0:10000), tch91_amm_pred[,2], lty = 3, col = 2)
  lines(c(0:10000), tch91_amm_pred[,3], lty = 3, col = 2)

  set.seed(43093)

  points(seq(0,10000,20), sapply(seq(0,10000,20)*1000, tch91_amm_v_unc, simplify = TRUE), pch = 17, col = 4, cex = 0.6)

  legend("topright", legend = c("Observed", "Simulated"), pch = c(16, 17), col = c(1,4), cex = 0.8, bty = 'n')

dev.off()


#Urea model #############
tch91_ure_pred <- predict(tch91.egv.ure, newdata = data.frame(conc = seq(0,40000,10)), 
                          type = "response", interval = "confidence", level = 0.95)

png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Eggs/urea_egg_viability_data_sim.png")

plot(eggv.ure$conc, eggv.ure$per_hatch, pch = 16, 
     xlab = "Urea concentration (ppm)", ylab = "% hatched")
  lines(seq(0,40000,10), tch91_ure_pred[,1], lty = 2, col = 2)
  lines(seq(0,40000,10), tch91_ure_pred[,2], lty = 3, col = 2)
  lines(seq(0,40000,10), tch91_ure_pred[,3], lty = 3, col = 2)


  points(seq(0,40000,40), sapply(seq(0,40000,40)*1000, tch91_ure_v_unc, simplify = TRUE), pch = 17, col = 4, cex = 0.6)

  legend("topright", legend = c("Observed", "Simulated"), pch = c(16, 17), col = c(1,4), cex = 0.8, bty = 'n')

dev.off()  