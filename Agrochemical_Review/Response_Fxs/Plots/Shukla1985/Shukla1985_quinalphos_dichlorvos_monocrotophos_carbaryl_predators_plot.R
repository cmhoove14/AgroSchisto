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

source("Agrochemical_Review/Response_Fxs/shukla1985_quinalphos_dichlorvos_monocrotophos_carbaryl_predators_fit.R")

#quinalphos plot
png("Agrochemical_Review/Response_Fxs/Plots/shukla1985/shukla1985_predator_mortality_quinalphos.png")

plot(c(0.655, 0.796, 0.969)*1000, c(.25,0.5,0.75), pch = 16, ylim = c(0,1), xlim = c(0,1000),
     xlab = "quinalphos (ppb)", ylab = "Predator daily mortality rate")
  segments(x0 = 0.717*1000, x1 = 0.884*1000, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,1000, 5), sapply(seq(0,1000, 5), muPq_quin_shukla85_uncertainty, simplify = T),
         pch = 5, cex = 0.5, col = 4)
  
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  

#dichlorvos plot
png("Agrochemical_Review/Response_Fxs/Plots/shukla1985/shukla1985_predator_mortality_dichlorvos.png")

plot(c(1.194, 1.435, 1.726)*1000, c(.25,0.5,0.75), pch = 16, ylim = c(0,1), xlim = c(0,2000),
     xlab = "dichlorvos (ppb)", ylab = "Predator daily mortality rate")
  segments(x0 = 1.303*1000, x1 = 1.581*1000, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,2000, 10), sapply(seq(0,2000, 10), muPq_dich_shukla85_uncertainty, simplify = T),
         pch = 5, cex = 0.5, col = 4)
  
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  

#monocrotophos plot
png("Agrochemical_Review/Response_Fxs/Plots/shukla1985/shukla1985_predator_mortality_monocrotophos.png")

plot(c(1.963, 2.107, 2.261)*1000, c(.25,0.5,0.75), pch = 16, ylim = c(0,1), xlim = c(0,2500),
     xlab = "monocrotophos (ppb)", ylab = "Predator daily mortality rate")
  segments(x0 = 2.030*1000, x1 = 2.186*1000, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,2500, 10), sapply(seq(0,2500, 10), muPq_mono_shukla85_uncertainty, simplify = T),
         pch = 5, cex = 0.5, col = 4)
  
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  

#Carbaryl plot
png("Agrochemical_Review/Response_Fxs/Plots/shukla1985/shukla1985_predator_mortality_carbaryl.png")

plot(c(0.029, 0.033, 0.038)*1000, c(.25,0.5,0.75), pch = 16, ylim = c(0,1), xlim = c(0,40),
     xlab = "Carbaryl (ppb)", ylab = "Predator daily mortality rate")
  segments(x0 = 0.031*1000, x1 = 0.036*1000, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,40, 0.1), sapply(seq(0,40, 0.1), muPq_carb_shukla85_uncertainty, simplify = T),
         pch = 5, cex = 0.5, col = 4)
  
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  