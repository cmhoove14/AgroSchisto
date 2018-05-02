#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Monde et al 2016 data found in tabl2 2 for bulinus globusus 

source("Agrochemical_Review/Response_Fxs/Monde2016_endosulfan_snail_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Monde2016/endosulfan_snails_sim.png")

plot(c(807, 4160, 21457), c(.01,0.5,0.9), pch = 16, xlab = "Endosulfan (ppb)", ylab = "Snail daily mortality rate", ylim = c(0,1))
  segments(x0 = 1336, x1 = 12956, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,25000, 50), sapply(seq(0,25000, 50), muNq_endo_monde16_uncertainty, simplify = T),
         pch = 5, cex = 0.5, col = 4)
  
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  