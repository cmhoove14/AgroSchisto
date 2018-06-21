#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/Rohr_unpublished_insecticides_macrobrachium_fit.R")
    
#Malathion ******************************************************************************########
png("Agrochemical_Review/Response_Fxs/Plots/Rohr_unpublished/rohr_macrobrachium_malathion.png")
  
  plot(mal_mac$Dose, mal_mac$nrx / mal_mac$ntot, pch = 16,
       ylab = "24-hr mortality rate", xlab = "malathion (ppb)",
       main = " Rohr unpublished macrobrachium toxicity data")
    lines(seq(0, max(mal_mac$Dose), max(mal_mac$Dose)/200), 
          predict(mal.mod.mac, newdata = data.frame(Dose = seq(0, max(mal_mac$Dose), max(mal_mac$Dose)/200))), 
          lty = 2, col = 2)
  set.seed(43093)
    points(seq(0, max(mal_mac$Dose), max(mal_mac$Dose)/200), 
           sapply(seq(0, max(mal_mac$Dose), max(mal_mac$Dose)/200), muPq_mal_mac_rohr_unpub_uncertainty),
           pch = 5, col = 4, cex = 0.5)
        
dev.off()

#Chlorpyrifos   *************************************************************************############
png("Agrochemical_Review/Response_Fxs/Plots/Rohr_unpublished/rohr_macrobrachium_chlorpyrifos.png")

  plot(chlor_mac$Dose, chlor_mac$nrx / chlor_mac$ntot, pch = 16, ylim = c(0,1),
       ylab = "24-hr mortality rate", xlab = "chlopyrifos (ppb)",
       main = " Rohr unpublished macrobrachium toxicity data")
    lines(seq(0, max(chlor_mac$Dose), max(chlor_mac$Dose)/200), 
          predict(chlor.mod.mac, newdata = data.frame(Dose = seq(0, max(chlor_mac$Dose), max(chlor_mac$Dose)/200))), 
          lty = 2, col = 2)
  set.seed(43093)
    points(seq(0, max(chlor_mac$Dose), max(chlor_mac$Dose)/200), 
           sapply(seq(0, max(chlor_mac$Dose), max(chlor_mac$Dose)/200), muPq_chlor_mac_rohr_unpub_uncertainty),
           pch = 5, col = 4, cex = 0.5)
      
dev.off()

#Terbufos   *****************************************************************************########
png("Agrochemical_Review/Response_Fxs/Plots/Rohr_unpublished/rohr_macrobrachium_terbufos.png")
    
  plot(terb_mac$Dose, terb_mac$nrx / terb_mac$ntot, pch = 16, ylim = c(0,1),
       ylab = "24-hr mortality rate", xlab = "Terbufos (ppb)",
       main = " Rohr unpublished macrobrachium toxicity data")
    lines(seq(0, max(terb_mac$Dose), max(terb_mac$Dose)/200), 
          predict(terb.mod.mac, newdata = data.frame(Dose = seq(0, max(terb_mac$Dose), max(terb_mac$Dose)/200))), 
          lty = 2, col = 2)
  set.seed(43093)
    points(seq(0, max(terb_mac$Dose), max(terb_mac$Dose)/200), 
           sapply(seq(0, max(terb_mac$Dose), max(terb_mac$Dose)/200), muPq_terb_mac_rohr_unpub_uncertainty),
           pch = 5, col = 4, cex = 0.5)
      
dev.off()  

#Lambda-cyhalothrin   *******************************************************************##########
png("Agrochemical_Review/Response_Fxs/Plots/Rohr_unpublished/rohr_macrobrachium_lambda-cyhalothrin.png")
    
  plot(lamcy_mac$Dose, lamcy_mac$nrx / lamcy_mac$ntot, pch = 16, ylim = c(0,1),
       ylab = "24-hr mortality rate", xlab = "Lambda-cyhalothrin (ppb)",
       main = " Rohr unpublished macrobrachium toxicity data")
    lines(seq(0, max(lamcy_mac$Dose), max(lamcy_mac$Dose)/200), 
          predict(lamcy.mod.mac, newdata = data.frame(Dose = seq(0, max(lamcy_mac$Dose), max(lamcy_mac$Dose)/200))), 
          lty = 2, col = 2)
  set.seed(43093)
    points(seq(0, max(lamcy_mac$Dose), max(lamcy_mac$Dose)/200), 
           sapply(seq(0, max(lamcy_mac$Dose), max(lamcy_mac$Dose)/200), muPq_lamcy_mac_rohr_unpub_uncertainty),
           pch = 5, col = 4, cex = 0.5)
      
dev.off()

#esfenvalerate  *************************************************************************#######
png("Agrochemical_Review/Response_Fxs/Plots/Rohr_unpublished/rohr_macrobrachium_esfenvalerate.png")
    
  plot(esfen_mac$Dose, esfen_mac$nrx / esfen_mac$ntot, pch = 16, ylim = c(0,1),
       ylab = "24-hr mortality rate", xlab = "Esfenvalerate (ppb)",
       main = " Rohr unpublished macrobrachium toxicity data")
    lines(seq(0, max(esfen_mac$Dose), max(esfen_mac$Dose)/200), 
          predict(esfen.mod.mac, newdata = data.frame(Dose = seq(0, max(esfen_mac$Dose), max(esfen_mac$Dose)/200))), 
          lty = 2, col = 2)
  set.seed(43093)
    points(seq(0, max(esfen_mac$Dose), max(esfen_mac$Dose)/200), 
           sapply(seq(0, max(esfen_mac$Dose), max(esfen_mac$Dose)/200), muPq_esfen_mac_rohr_unpub_uncertainty),
           pch = 5, col = 4, cex = 0.5)
      
dev.off()

#Permethrin *****************************************************************************##########
png("Agrochemical_Review/Response_Fxs/Plots/Rohr_unpublished/rohr_macrobrachium_permethrin.png")
    
  plot(perm_mac$Dose, perm_mac$nrx / perm_mac$ntot, pch = 16, ylim = c(0,1),
       ylab = "24-hr mortality rate", xlab = "Permethrin (ppb)",
       main = " Rohr unpublished macrobrachium toxicity data")
    lines(seq(0, max(perm_mac$Dose), max(perm_mac$Dose)/200), 
          predict(perm.mod.mac, newdata = data.frame(Dose = seq(0, max(perm_mac$Dose), max(perm_mac$Dose)/200))), 
          lty = 2, col = 2)
  set.seed(43093)
    points(seq(0, max(perm_mac$Dose), max(perm_mac$Dose)/200), 
           sapply(seq(0, max(perm_mac$Dose), max(perm_mac$Dose)/200), muPq_perm_mac_rohr_unpub_uncertainty),
           pch = 5, col = 4, cex = 0.5)
      
dev.off()