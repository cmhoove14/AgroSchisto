#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Data extraction and model fitting to Hussein 2011 [Bakry et al 2016](https://www.researchgate.net/publication/305149067_Effects_of_Three_Inorganic_Fertilizers_on_the_Biology_and_Histopathology_of_infected_Biomphalaria_alexandrina_snails?enrichId=rgreq-cf1a38509b7c5460cbc30f698ea58594-XXX&enrichSource=Y292ZXJQYWdlOzMwNTE0OTA2NztBUzozODI2MjU5NTgxMjE0NzJAMTQ2ODIzNjU0NTE4NQ%3D%3D&el=1_x_3&_esc=publicationCoverPdf) data
source("Agrochemical_Review/Response_Fxs/Hussein2016_fertilizer_snails_fit.R")

#balanced fertilizer reported LC50 and slope data #####################
png("Agrochemical_Review/Response_Fxs/Plots/Hussein2016/balanced_fertilizer_snail_mortality.png")

plot(c(0,505.7,851.2), c(0,.50,.90), pch = 16, cex = 1.2, ylim = c(0,1),
     xlab = "Balanced fertilizer (ppm)", ylab = "mortality from reported LC values",
     main = "B. alexandrina mortality, Hussein et al 2016")
  segments(x0 = 390.31, x1 = 625.33, y0 = 0.5,y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,900,2), sapply(seq(0,900000, 2000), muNq_balanced_hussein16_uncertainty, simplify = TRUE),
         pch = 5, cex = 0.5, col = 4)
  
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()

#high phosphorous fertilizer reported LC50 and slope data #####################
png("Agrochemical_Review/Response_Fxs/Plots/Hussein2016/highp_fertilizer_snail_mortality.png")

plot(c(0,1600,2309.3), c(0,.50,.90), pch = 16, cex = 1.2, ylim = c(0,1),
     xlab = "High P fertilizer (ppm)", ylab = "mortality from reported LC values",
     main = "B. alexandrina mortality, Hussein et al 2016")
  segments(x0 = 1356.56, x1 = 1843.4, y0 = 0.5,y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,2400,5), sapply(seq(0,2400000, 5000), muNq_highp_hussein16_uncertainty, simplify = TRUE),
         pch = 5, cex = 0.5, col = 4)
  
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()

#balanced fertilizer reported LC50 and slope data #####################
png("Agrochemical_Review/Response_Fxs/Plots/Hussein2016/highn_fertilizer_snail_mortality.png")

plot(c(0,9500,11273.3), c(0,.50,.90), pch = 16, cex = 1.2, ylim = c(0,1),
     xlab = "High N fertilizer (ppm)", ylab = "mortality from reported LC values",
     main = "B. alexandrina mortality, Hussein et al 2016")
  segments(x0 = 8891.4, x1 = 10108.5, y0 = 0.5,y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,12000,50), sapply(seq(0,12000,50)*1000, muNq_highn_hussein16_uncertainty, simplify = TRUE),
         pch = 5, cex = 0.5, col = 4)
  
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()
