#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(LW1949)
source("Agrochemical_Review/Response_Fxs/barbieri2016_carbofuran_predators_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Barbieri2016/carbofuran_mupq_sim_observed.png")

plotDE(barb_dat, xlab = "Carbofuran (ppm)", main = "Barbieri 2016 toxicity to Macrobrachium olfersii")

set.seed(43093)

points(seq(0,4,0.01), sapply(seq(0,4,0.01)*1000, barb_carbofuran_muPq_uncertainty)*100, pch = 5, cex = 0.6, col = 4)
points(seq(0,4,0.01), sapply(seq(0,4,0.01)*1000, muPq_carb_barb16_uncertainty)*100, pch = 5, cex = 0.6, col = 2)
  legend("bottomright", legend = c("Observed", "Simulated-LW1949", "Simulated-DRC"), pch = c(16,5,5), col = c(1,4,2), bty = 'n', cex = 0.7)

dev.off()