#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Ragab 2006 data
source("Agrochemical_Review/Response_Fxs/ragab2006_fertilizers_snails_fit.R")

#Ammonium Nitrate Snail toxicity ##########
png("Agrochemical_Review/Response_Fxs/Plots/Ragab&Shoukry2006/Ragab&Shoukry2006_function_simulate_muN_ammonium_nitrate.png")

plot(c(0,lc50.amm.rag, lc90.amm.rag)*1000, c(0,0.5, 0.9), pch = 16,
     xlim = c(0, (lc90.amm.rag+100)*1000), ylim = c(0,1),
     xlab = 'amm. phosphate (ppb)', ylab = 'mortality',
     main = expression(paste('Ragab06 Ammonium phosphate toxicity to ', 
                             italic('Bi. alexandrina'))))
  segments(x0 = 435.19*1000, x1 = 507.6*1000, y0 = 0.5, y1 = 0.5)

  set.seed(43093)
  
  points(seq(0, (lc90.amm.rag+100)*1000, 2500),
         sapply(seq(0, (lc90.amm.rag+100)*1000, 2500), rag06_mun_amm),
         pch = 5, col = 4, cex = 0.5)

dev.off()

#Pottasium Sulfate Snail toxicity ##########    
png("Agrochemical_Review/Response_Fxs/Plots/Ragab&Shoukry2006/Ragab&Shoukry2006_function_simulate_muN_potassium_sulphate.png")

plot(c(0, lc50.pot.rag, lc90.pot.rag)*1000, c(0, 0.5, 0.9), pch = 16,
     xlim = c(0, (lc90.pot.rag+100)*1000), ylim = c(0,1),
     xlab = 'pot. sulfate (ppb)', ylab = 'mortality',
     main = expression(paste('Ragab06 Potassium Sulfate toxicity to ', 
                             italic('Bi. alexandrina'))))
segments(x0 = 1583.3*1000, x1 = 2280*1000, y0 = 0.5, y1 = 0.5)
set.seed(43093)

points(seq(0, (lc90.pot.rag+100)*1000, 5000),
       sapply(seq(0, (lc90.pot.rag+100)*1000, 5000), rag06_mun_pot),
       pch = 5, col = 4, cex = 0.5)

dev.off()

#Urea Snail toxicity ##########    
png("Agrochemical_Review/Response_Fxs/Plots/Ragab&Shoukry2006/Ragab&Shoukry2006_function_simulate_muN_urea.png")

plot(c(0, lc50.urea.rag, lc90.urea.rag)*1000, c(0, 0.5, 0.9), pch = 16,
     xlim = c(0, (lc90.urea.rag+100)*1000), ylim = c(0,1),
     xlab = 'urea (ppb)', ylab = 'mortality',
     main = expression(paste('Ragab06 urea toxicity to ', 
                             italic('Bi. alexandrina'))))
  segments(x0 = 24860*1000, x1 = 19469*1000, y0 = 0.5, y1 = 0.5)

  set.seed(43093)
  
  points(seq(0, (lc90.urea.rag+100)*1000, 50000),
         sapply(seq(0, (lc90.urea.rag+100)*1000, 50000), rag06_mun_urea),
         pch = 5, col = 4, cex = 0.5)

dev.off()