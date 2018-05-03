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

source("Agrochemical_Review/Response_Fxs/omkar1985_endosulfan_phosphamidon_carbaryl_predators_fit.R")

#Endosulfan plot
png("Agrochemical_Review/Response_Fxs/Plots/Omkar1985/omkar1985_predator_mortality_endosulfan.png")

plot(c(0.0052, 0.0062, 0.0073)*1000, c(.25,0.5,0.75), pch = 16, ylim = c(0,1), xlim = c(0,8),
     xlab = "Endosulfan (ppb)", ylab = "Predator daily mortality rate")
  segments(x0 = 0.0057*1000, x1 = 0.0068*1000, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,8, 0.05), sapply(seq(0,8, 0.05), muPq_endo_omkar85_uncertainty, simplify = T),
         pch = 5, cex = 0.5, col = 4)
  
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  

#Phosphamidon plot
png("Agrochemical_Review/Response_Fxs/Plots/Omkar1985/omkar1985_predator_mortality_phosphamidon.png")

plot(c(3.783, 4.825, 6.153)*1000, c(.25,0.5,0.75), pch = 16, ylim = c(0,1), xlim = c(0,7000),
     xlab = "Phosphamidon (ppb)", ylab = "Predator daily mortality rate")
  segments(x0 = 4.221*1000, x1 = 5.429*1000, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,8000, 10), sapply(seq(0,8000, 10), muPq_phos_omkar85_uncertainty, simplify = T),
         pch = 5, cex = 0.5, col = 4)
  
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  

#Carbaryl plot
png("Agrochemical_Review/Response_Fxs/Plots/Omkar1985/omkar1985_predator_mortality_carbaryl.png")

plot(c(0.0461, 0.0513, 0.0572)*1000, c(.25,0.5,0.75), pch = 16, ylim = c(0,1), xlim = c(0,60),
     xlab = "Carbaryl (ppb)", ylab = "Predator daily mortality rate")
  segments(x0 = 0.0485*1000, x1 = 0.0541*1000, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,60, 0.1), sapply(seq(0,60, 0.1), muPq_carb_omkar85_uncertainty, simplify = T),
         pch = 5, cex = 0.5, col = 4)
  
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  