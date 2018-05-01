#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Ghaffar 2016 SNAIL (B. alexandrina) data
source("Agrochemical_Review/Response_Fxs/Ghaffar2016_butralin_glyphosate_pendimethalin_snails_fit.R")

#Snail (30 B. alexandrina 6-8mm shell width) toxicity ##########
#Butralin ##############
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Snails/butralin_muN_function.png")

plot(but.dat$butralin/1000, but.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,13),
       xlab = 'Butralin (ppm)', ylab = 'snail mortality rate')
  lines(seq(0, 12, 0.01), sapply(seq(0, 12, 0.01)*1000, muNq_but_Gaf16_lin, simplify = TRUE), 
        lty = 2, col = 2)
  lines(seq(0, 12, 0.01), sapply(seq(0, 12, 0.01)*1000, muNq_but_Gaf16_lin, b = gaf.b.lin.up, simplify = TRUE), 
        lty = 3, col = 2)
  lines(seq(0, 12, 0.01), sapply(seq(0, 12, 0.01)*1000, muNq_but_Gaf16_lin, b = gaf.b.lin.lo, simplify = TRUE), 
        lty = 3, col = 2)

dev.off() 

png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Snails/butralin_muN_simulate.png")

  plot(but.dat$butralin, but.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,13000),
       xlab = 'Butralin (ppb)', ylab = 'snail mortality rate', 
       main = 'D-R function based on reported values')
    segments(x0 = 3700, x1 = 8340, y0 = 0.5, y1 = 0.5, lty = 1, col = 1)
    segments(x0 = 6220, x1 = 12800, y0 = 0.9, y1 = 0.9, lty = 1, col = 1)
  
    set.seed(43093)
    
    points(seq(0,13000,50), sapply(seq(0,13000,50), mu_Nq_butr_gaf16_uncertainty, simplify = T), 
           pch = 5, col = 4, cex = 0.5)
    
    legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()    

#Glyphosate #############
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Snails/glyphosate_muN_function.png")

plot(gly.dat$glyphosate/1000, gly.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,33.000),
     xlab = 'Glyphosate (ppm)', ylab = 'snail mortality rate',
     main = 'D-R function based on reported values')
  lines(seq(0, 33, 0.01), sapply(seq(0, 33, 0.01)*1000, muNq_gly_Gaf16_lin, simplify = TRUE), 
        lty = 2, col = 2)
  lines(seq(0, 33, 0.01), sapply(seq(0, 33, 0.01)*1000, muNq_gly_Gaf16_lin, b = gaf.b.lin.up.gly, simplify = TRUE), 
        lty = 3, col = 2)
  lines(seq(0, 33, 0.01), sapply(seq(0, 33, 0.01)*1000, muNq_gly_Gaf16_lin, b = gaf.b.lin.lo.gly, simplify = TRUE), 
        lty = 3, col = 2)

dev.off()

png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Snails/glyphosate_muN_simulate.png")

plot(gly.dat$glyphosate, gly.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,33000),
     xlab = 'Glyphosate (ppm)', ylab = 'snail mortality rate',
     main = 'D-R function based on reported values')
  segments(x0 = 9130, x1 = 16570, y0 = 0.5, y1 = 0.5, lty = 1, col = 1)
  segments(x0 = 23870, x1 = 28900, y0 = 0.9, y1 = 0.9, lty = 1, col = 1) 

    set.seed(43093)
    
  points(seq(0,33000,100), sapply(seq(0,33000,100), mu_Nq_gly_gaf16_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)
  
   legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()   

#Pendimethalin ##############

#Create function based on LC values
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Snails/pendimethalin_muN_function.png")

plot(pen.dat$pendimethalin/1000, pen.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,7),
     xlab = 'pendimethalin (ppm)', ylab = 'snail mortality rate',
     main = 'D-R function based on reported values')
  lines(seq(0, 7, 0.01), sapply(seq(0, 7, 0.01)*1000, muNq_pen_Gaf16_lin, simplify = TRUE), 
        lty = 2, col = 2)
  lines(seq(0, 7, 0.01), sapply(seq(0, 7, 0.01)*1000, muNq_pen_Gaf16_lin, b = gaf.b.lin.up.pen, simplify = TRUE), 
        lty = 3, col = 2)
  lines(seq(0, 7, 0.01), sapply(seq(0, 7, 0.01)*1000, muNq_pen_Gaf16_lin, b = gaf.b.lin.lo.pen, simplify = TRUE), 
        lty = 3, col = 2)

dev.off()

png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Snails/pendimethalin_muN_simulate.png")

  plot(pen.dat$pendimethalin, pen.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,7000),
       xlab = 'pendimethalin (ppb)', ylab = 'snail mortality rate',
       main = 'D-R function based on reported values')
    segments(x0 = 1430, x1 = 3220, y0 = 0.5, y1 = 0.5, lty = 1, col = 1)
    segments(x0 = 2400, x1 = 6420, y0 = 0.9, y1 = 0.9, lty = 1, col = 1)
  
    set.seed(43093)
    
    points(seq(0,7500,25), sapply(seq(0,7500,25), mu_Nq_pen_gaf16_uncertainty, simplify = T), 
           pch = 5, col = 4, cex = 0.5)
    
    legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()
