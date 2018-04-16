#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/Satapornvanit2009_insecticides_predators_fit.R")

#24-hr mortality functions for predator (macrobrachium) populations from Satapornvanit et al chemosphere paper ########
#Zinc ########    
png("Agrochemical_Review/Response_Fxs/Plots/Satapornvanit2009/Satapornvanit2009_data_functions_zinc_mupq.png")

  plot(sap.mort$conc[sap.mort$chem == 'zinc'], sap.mort$mort[sap.mort$chem == 'zinc']/100, ylim = c(0,1),
       pch = 16, xlab = 'Zinc concentration (ppb)', ylab = 'prop dead', 
       main = expression(paste('Zinc daily toxicity to post-larval ', italic('M. rosenbergii'), sep = '')))

  set.seed(43093)

  points(seq(0,900,10), sapply(seq(0,900,10), muPq_zinc_satapornvanit09_uncertainty, simplify = T), 
         pch=5, col=4, cex = 0.5)

  points(seq(0,900,10), sapply(seq(0,900,10), muPq_zinc_satapornvanit09_uncertainty_lw49),
         pch = 5, col=6, cex = 0.5)
    
    legend("bottomright", legend = c("Observed", "Simulated-DRC", "Simulated-LW1949"), pch = c(16,5,5), col = c(1,4,6), bty = 'n', cex = 0.7)

dev.off()  

#Chlorpyrifos ###########    
png("Agrochemical_Review/Response_Fxs/Plots/Satapornvanit2009/Satapornvanit2009_data_functions_chlorpyrifos_mupq.png")
  plot(sap.mort$conc[sap.mort$chem == 'chlorpyrifos'], sap.mort$mort[sap.mort$chem == 'chlorpyrifos']/100, ylim = c(0,1),
       pch = 16, xlab = 'chlorpyrifos concentration (ppb)', ylab = 'prop dead', 
       main = expression(paste('Chlorpyrifos daily toxicity to post-larval ', italic('M. rosenbergii'), sep = '')))
    
    set.seed(43093)
    
  points(seq(0,5,0.05), sapply(seq(0,5,0.05), muPq_chlor_satapornvanit09_uncertainty, simplify = T), 
          pch=5, col=4, cex = 0.5)    

  points(seq(0,5,0.02), sapply(seq(0,5,0.02), muPq_chlor_satapornvanit09_uncertainty_lw49),
         pch = 5, col=6, cex = 0.5)
  
    legend("bottomright", legend = c("Observed", "Simulated-DRC", "Simulated-LW1949"), pch = c(16,5,5), col = c(1,4,6), bty = 'n', cex = 0.7)

dev.off()

#Dimethoate ################        
png("Agrochemical_Review/Response_Fxs/Plots/Satapornvanit2009/Satapornvanit2009_data_functions_deimethoate_mupq.png")
 plot(sap.mort$conc[sap.mort$chem == 'dimethoate'], sap.mort$mort[sap.mort$chem == 'dimethoate']/100, ylim = c(0,1),
      pch = 16, xlab = 'dimethoate concentration (ppb)', ylab = 'prop dead', 
      main = expression(paste('Dimethoate daily toxicity to post-larval ', italic('M. rosenbergii'), sep = '')))
  
   set.seed(43093)

  points(seq(0,1300,10), sapply(seq(0,1300,10), muPq_dim_satapornvanit09_uncertainty, simplify = T), 
          pch=5, col=4, cex = 0.5)       
  
  points(seq(0,1300,10), sapply(seq(0,1300,10), muPq_dim_satapornvanit09_uncertainty_lw49),
         pch = 5, col=6, cex = 0.5)
    
    legend("bottomright", legend = c("Observed", "Simulated-DRC", "Simulated-LW1949"), pch = c(16,5,5), col = c(1,4,6), bty = 'n', cex = 0.7)

dev.off()

#Profenofos ################      
png("Agrochemical_Review/Response_Fxs/Plots/Satapornvanit2009/Satapornvanit2009_data_functions_profenofos_mupq.png")
      
  plot(sap.mort$conc[sap.mort$chem == 'profenofos'], sap.mort$mort[sap.mort$chem == 'profenofos']/100, ylim = c(0,1),
       pch = 16, xlab = 'profenofos concentration (ppb)', ylab = 'prop dead', 
       main = expression(paste('Profenofos daily toxicity to post-larval ', italic('M. rosenbergii'), sep = '')))  
  
  set.seed(43093)
  
  points(seq(0,60,0.5), sapply(seq(0,60,0.5), muPq_prof_satapornvanit09_uncertainty, simplify = T), 
          pch=5, col=4, cex = 0.5)

  points(seq(0,60,1), sapply(seq(0,60,1), muPq_prof_satapornvanit09_uncertainty_lw49),
         pch = 5, col=6, cex = 0.5)
  
  legend("bottomright", legend = c("Observed", "Simulated-DRC", "Simulated-LW1949"), pch = c(16,5,5), col = c(1,4,6), bty = 'n', cex = 0.7)

dev.off()  

#Feeding rate analyses ######
png("Agrochemical_Review/Response_Fxs/Plots/Satapornvanit2009/Satapornvanit2009_data_functions_zinc_psiq.png")

  plot(fr.z$conc, fr.z$feed_rate / fr.z$feed_rate[1], pch = 16, ylim = c(0,1),
       xlab = 'Zinc concentration (ppb)', ylab = 'Feeding rate (prey/prawn/hr)', 
       main = 'Prawn reduced feeding (Satapornvanit-09)')
    lines(seq(0,900,10), sapply(seq(0,900,10),psi_q_zinc_satapornvanit09, simplify = T)[1,],
          lty = 2, col = 2)
    lines(seq(0,900,10), sapply(seq(0,900,10),psi_q_zinc_satapornvanit09, simplify = T)[2,],
          lty = 3, col = 2)
    lines(seq(0,900,10), sapply(seq(0,900,10),psi_q_zinc_satapornvanit09, simplify = T)[3,],
          lty = 3, col = 2)  
    
    set.seed(43093)
    
    points(seq(0,900,10), sapply(seq(0,900,10), psi_q_zinc_satapornvanit09_uncertainty, simplify = T), 
           pch=5, col=4, cex = 0.5)   

dev.off()

#Chlorpyrifos ###################
png("Agrochemical_Review/Response_Fxs/Plots/Satapornvanit2009/Satapornvanit2009_data_functions_chlorpyrifos_psiq.png")

  plot(fr.ch$conc, fr.ch$feed_rate / fr.ch$feed_rate[1], pch = 16, ylim = c(0,1),
       xlab = 'Chlorpyrifos concentration (ppb)', ylab = 'Feeding rate (prey/prawn/hr)', 
       main = 'Prawn reduced feeding (Sata-09)')
    lines(seq(0,5,0.01), sapply(seq(0,5,0.01), psi_q_chlor_satapornvanit09, simplify = T)[1,],
          lty = 2, col = 2)
    lines(seq(0,5,0.01), sapply(seq(0,5,0.01), psi_q_chlor_satapornvanit09, simplify = T)[2,],
          lty = 3, col = 2)
    lines(seq(0,5,0.01), sapply(seq(0,5,0.01), psi_q_chlor_satapornvanit09, simplify = T)[3,],
          lty = 3, col = 2)
    
    set.seed(43093)
    
    points(seq(0,5,0.05), sapply(seq(0,5,0.05), psi_q_chlor_satapornvanit09_uncertainty, simplify = T), 
           pch=5, col=4, cex = 0.5)   

dev.off()