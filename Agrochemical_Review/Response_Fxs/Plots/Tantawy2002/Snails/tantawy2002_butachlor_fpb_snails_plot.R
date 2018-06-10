#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
#Data extraction and model fitting to Tantawy 2002 data

source("Agrochemical_Review/Response_Fxs/tantawy2002_butachlor_fpb_snails_fit.R")

#Snail toxicity ##########
#Butachlor ##########
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Snails/tantawy2002_butachlor_muN_data_sim.png")
  plot(snail.but$conc, snail.but$mort, pch = 16, ylim = c(0,1), xlim = c(0,max(snail.but$conc)+100),
       xlab = 'Butachlor (ppb)', ylab = 'Snail mortality',
       main = expression(paste('Butachlor toxicity to ', italic('Bi. alexandrina'))))
    segments(x0 = 4060, x1 = 10400, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0, 50000, 100), sapply(seq(0, 50000, 100), muNq.tant.but_uncertainty),
         pch = 5, cex = 0.5, col = 4)

dev.off()    

#Fluazifop-p-butyl  ##########     
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Snails/tantawy2002_fpb_muN_data_sim.png")  

plot(snail.fpb$conc, snail.fpb$mort, pch = 16, ylim = c(0,1), xlim = c(0,max(snail.fpb$conc)+100),
     xlab = 'fluazifop-p-butyl (ppb)', ylab = 'Snail mortality',
     main = expression(paste('fluazifop-p-butyl toxicity to ', italic('Bi. alexandrina'))))
  segments(x0 = 11730, x1 = 26400, y0 = 0.5, y1 = 0.5)

  set.seed(43093)
  
  points(seq(0, 60000, 100), sapply(seq(0, 60000, 100), muNq.tant.fpb_uncertainty),
         pch = 5, cex = 0.5, col = 4)
dev.off()