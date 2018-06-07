#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/Lavarias2013_fenitrothion_predators_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Lavarias13/Lavarias13_fenitrothion_m_borellii.png")

plot(lav13_dat$conc, lav13_dat$mort/100, pch = 16, xlim = c(0, 20), ylim = c(0,1),
     xlab = "Fenitrothion (ppb)", ylab = "mortality",
     main = "Sarojini 1986 M. borelli")

  set.seed(43093)
  
  points(seq(0,20,0.2), sapply(seq(0,20,0.2), muPq_fenitrothion_Lavarias13_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()