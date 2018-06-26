#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Data extraction and model fitting to Mohamed et al 2012 data
source("Agrochemical_Review/Response_Fxs/Mohamed2012_profenofos_diazinon_snails_fit.R")

morts <-c(0, 10, 25, 50, 90)/100
#Plot diazinon data and function fit
diaz_lcs <- c(1.41, 10.53, 12.25, 14.16, 17.78)

png("Agrochemical_Review/Response_Fxs/Plots/Mohamed2012/diazinon_snails_data_sim.png")
  plot(diaz_lcs, morts, pch = 16, ylim = c(0,1),
       xlab = "Diazinon (ppm)", ylab = "24-hr mortality")
    
    lines(seq(0,18,0.1), sapply(seq(0,18,0.1)*1000, muNq_diaz_mohamed, simplify = TRUE), lty = 2, col = 2)
    lines(seq(0,18,0.1), sapply(seq(0,18,0.1)*1000, muNq_diaz_mohamed, 
                                lc50 = 10^(log10(lc50.moh.diaz.report)-1.96*se.lc50.moh.diaz), simplify = TRUE), 
          lty = 3, col = 2)
    lines(seq(0,18,0.1), sapply(seq(0,18,0.1)*1000, muNq_diaz_mohamed, 
                                lc50 = 10^(log10(lc50.moh.diaz.report)+1.96*se.lc50.moh.diaz), simplify = TRUE), 
          lty = 3, col = 2)

  set.seed(43093)
  
  points(seq(0,18,0.1), sapply(seq(0,18,0.1)*1000, muNq_diaz_mohamed_uncertainty, simplify = TRUE), pch = 5, col = 4, cex = 0.5)
  
  legend("topleft", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()


#Plot profenofos data and function fit
prof_lcs <- c(0.402, 2.44, 3.11, 4.02, 5.57)

png("Agrochemical_Review/Response_Fxs/Plots/Mohamed2012/profenofos_snails_data_sim.png")
  plot(prof_lcs, morts, pch = 16, ylim = c(0,1),
       xlab = "Profenofos (ppm)", ylab = "24-hr mortality")
    
    lines(seq(0,6,0.1), sapply(seq(0,6,0.1)*1000, muNq_prof_mohamed, simplify = TRUE), lty = 2, col = 2)
    lines(seq(0,6,0.1), sapply(seq(0,6,0.1)*1000, muNq_prof_mohamed, 
                                lc50 = 10^(log10(lc50.moh.prof.report)-1.96*se.lc50.moh.prof), simplify = TRUE), 
          lty = 3, col = 2)
    lines(seq(0,6,0.1), sapply(seq(0,6,0.1)*1000, muNq_prof_mohamed, 
                                lc50 = 10^(log10(lc50.moh.prof.report)+1.96*se.lc50.moh.prof), simplify = TRUE), 
          lty = 3, col = 2)

  set.seed(43093)
  
  points(seq(0,6,0.1), sapply(seq(0,6,0.1)*1000, muNq_prof_mohamed_uncertainty, simplify = TRUE), pch = 5, col = 4, cex = 0.5)
  
  legend("topleft", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()

#Plot diazinon reproduction data, function, and sim
png("Agrochemical_Review/Response_Fxs/Plots/Mohamed2012/diazinon_snail_reproduction_data_sim.png")

plot(diaz_conc, diaz_rate, pch = 16, xlab = "Diazinon concentration", ylab = "Hatchlings/snail/week")
  lines(seq(0,13,0.1), sapply(seq(0,13,0.1)*1000, fNq_moh_diaz_moh12, simplify = TRUE)[1,], lty = 2, col = 2)
  lines(seq(0,13,0.1), sapply(seq(0,13,0.1)*1000, fNq_moh_diaz_moh12, simplify = TRUE)[2,], lty = 3, col = 2)
  lines(seq(0,13,0.1), sapply(seq(0,13,0.1)*1000, fNq_moh_diaz_moh12, simplify = TRUE)[3,], lty = 3, col = 2)
  
  set.seed(43093)
  
  points(seq(0,13,0.05), sapply(seq(0,13,0.05)*1000, fNq_moh_diaz_moh12_uncertainty, simplify = TRUE)*moh_repro_ref,
         pch = 5, cex = 0.5, col = 4)
  
  legend("topright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  

#plot profenofos reproduction data, funcion, and sim 
png("Agrochemical_Review/Response_Fxs/Plots/Mohamed2012/profenofos_snail_reproduction_data_sim.png")

plot(prof_conc, prof_rate, pch = 16, xlab = "Profenofos concentration", ylab = "Hatchlings/snail/week")
  lines(seq(0,3,0.01), sapply(seq(0,3,0.01)*1000, fNq_moh_prof_moh12, simplify = TRUE)[1,], lty = 2, col = 2)
  lines(seq(0,3,0.01), sapply(seq(0,3,0.01)*1000, fNq_moh_prof_moh12, simplify = TRUE)[2,], lty = 3, col = 2)
  lines(seq(0,3,0.01), sapply(seq(0,3,0.01)*1000, fNq_moh_prof_moh12, simplify = TRUE)[3,], lty = 3, col = 2)
  
  set.seed(43093)
  
  points(seq(0,3,0.025), sapply(seq(0,3,0.025)*1000, fNq_moh_prof_moh12_uncertainty, simplify = TRUE)*moh_repro_ref,
         pch = 5, cex = 0.5, col = 4)
  
  legend("topright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()
