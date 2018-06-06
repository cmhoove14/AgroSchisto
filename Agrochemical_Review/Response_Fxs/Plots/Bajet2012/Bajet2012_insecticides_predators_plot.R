#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Bajet et al 2012 uinsecticide toxicity to Macrobrachium lar##########
#Paper reports LC50, LC10, NOEC, and LOEC values. To obtain d-r function, slope is estimated
#by linear regression fit to the LC50 and LC10 values and incorporated with reported LC50 and its uncertainty
#All data points are then used in plots for qualitative validation of the resulting function

source("Agrochemical_Review/Response_Fxs/Bajet2012_insecticides_predators_fit.R")

#Lambda cyhalothrin ##############
png("Agrochemical_Review/Response_Fxs/Plots/Bajet2012/Bajet2012_lamcy_pred.png")

plot(c(noec.lamcy.baj, lc10.lamcy.baj, lc50.lamcy.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.lamcy.baj*2), ylim = c(0,1), pch = 16,
     xlab = "Lambda-Cyhalothrin (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.lamcy.baj.lo, x1 = lc50.lamcy.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.lamcy.baj*2,lc50.lamcy.baj/20), 
         sapply(seq(0,lc50.lamcy.baj*2,lc50.lamcy.baj/20), muPq_lamcy_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()

#Deltamethrin ##############
png("Agrochemical_Review/Response_Fxs/Plots/Bajet2012/Bajet2012_deltamethrin_pred.png")

plot(c(noec.delt.baj, lc10.delt.baj, lc50.delt.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.delt.baj*2), ylim = c(0,1), pch = 16,
     xlab = "Deltamethrin (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.delt.baj.lo, x1 = lc50.delt.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.delt.baj*2,lc50.delt.baj/20), 
         sapply(seq(0,lc50.delt.baj*2,lc50.delt.baj/20), muPq_deltamethrin_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

  
dev.off()

#Cypermethrin ##############
png("Agrochemical_Review/Response_Fxs/Plots/Bajet2012/Bajet2012_cypermethrin_pred.png")

plot(c(noec.cyper.baj, lc10.cyper.baj, lc50.cyper.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.cyper.baj*2), ylim = c(0,1), pch = 16,
     xlab = "cypermethrin (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.cyper.baj.lo, x1 = lc50.cyper.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.cyper.baj*2,lc50.cyper.baj/20), 
         sapply(seq(0,lc50.cyper.baj*2,lc50.cyper.baj/20), muPq_cypermethrin_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()

#Chlorpyrifos ##############
png("Agrochemical_Review/Response_Fxs/Plots/Bajet2012/Bajet2012_chlorpyrifos_pred.png")

plot(c(noec.chlor.baj, lc10.chlor.baj, lc50.chlor.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.chlor.baj*2), ylim = c(0,1), pch = 16,
     xlab = "chlorpyrifos (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.chlor.baj.lo, x1 = lc50.chlor.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.chlor.baj*2,lc50.chlor.baj/20), 
         sapply(seq(0,lc50.chlor.baj*2,lc50.chlor.baj/20), muPq_chlorpyrifos_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()

#Profenofos ##############
png("Agrochemical_Review/Response_Fxs/Plots/Bajet2012/Bajet2012_profenofos_pred.png")

plot(c(noec.prof.baj, lc10.prof.baj, lc50.prof.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.prof.baj*2), ylim = c(0,1), pch = 16,
     xlab = "profenofos (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.prof.baj.lo, x1 = lc50.prof.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.prof.baj*2,lc50.prof.baj/20), 
         sapply(seq(0,lc50.prof.baj*2,lc50.prof.baj/20), muPq_profenofos_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()

#Malathion ##############
png("Agrochemical_Review/Response_Fxs/Plots/Bajet2012/Bajet2012_malathion_pred.png")

plot(c(noec.mal.baj, lc10.mal.baj, lc50.mal.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.mal.baj*2), ylim = c(0,1), pch = 16,
     xlab = "malathion (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.mal.baj.lo, x1 = lc50.mal.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.mal.baj*2,lc50.mal.baj/20), 
         sapply(seq(0,lc50.mal.baj*2,lc50.mal.baj/20), muPq_malathion_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()

#Carbaryl ##############
png("Agrochemical_Review/Response_Fxs/Plots/Bajet2012/Bajet2012_carbaryl_pred.png")

plot(c(noec.carb.baj, lc10.carb.baj, lc50.carb.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.carb.baj*2), ylim = c(0,1), pch = 16,
     xlab = "carbaryl (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.carb.baj.lo, x1 = lc50.carb.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.carb.baj*2,lc50.carb.baj/20), 
         sapply(seq(0,lc50.carb.baj*2,lc50.carb.baj/20), muPq_carbaryl_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()

#2,4-D ##############

#Butachlor ##############
png("Agrochemical_Review/Response_Fxs/Plots/Bajet2012/Bajet2012_butachlor_pred.png")

plot(c(noec.but.baj, lc10.but.baj, lc50.but.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.but.baj*2), ylim = c(0,1), pch = 16,
     xlab = "butachlor (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.but.baj.lo, x1 = lc50.but.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.but.baj*2,lc50.but.baj/20), 
         sapply(seq(0,lc50.but.baj*2,lc50.but.baj/20), muPq_butachlor_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()