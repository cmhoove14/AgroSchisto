#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#herbicide toxicity to Bi. alexandrina from Omran and Salama 2013 ###############

source("Agrochemical_Review/Response_Fxs/Omran&Salama2013_atrazine_glyphosate_snails_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Omran&Salama2013/Omran&Salama2013_function_simulate_muN_atrazine.png")

  plotDE(atr.dat, xlab = "Dose (ppm)")
  predLines(fatr)
  points(seq(0,500,1), sapply(seq(0,500000,1000), ons.munq.atr)*100,
         pch = 5, col = 4, cex = 0.5)
  
  legend('bottomright', pch = c(16,5), col = c(1,4), legend = c("observed", "simulated"),
         cex = 0.8, bty = 'n')

dev.off()

#Glyphosate mortality #########
png("Agrochemical_Review/Response_Fxs/Plots/Omran&Salama2013/Omran&Salama2013_function_simulate_muN_glyphosate.png")
  
  plotDE(gly.dat)
  predLines(fgly)
  points(seq(0,500,1), sapply(seq(0,500000,1000), ons.munq.gly)*100,
         pch = 5, col = 4, cex = 0.5)
  
  legend('bottomright', pch = c(16,5), col = c(1,4), legend = c("observed", "simulated"),
         cex = 0.8, bty = 'n')

dev.off()    