source("Agrochemical_Review/Response_Fxs/Oliveira2009_endosulfan_snail_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Oliveira2009/endosulfan_snail_reproduction_data_sim.png")
  plot(oliv_dat$endo, oliv_dat$net, pch = 16, ylim = c(0,550),
       xlab = "Endosulfan (ppm)", ylab = "hatchlings per snail")
    
    lines(seq(0,0.12,0.001), sapply(seq(0,0.12,0.001)*1000, fNq_endo_oliv09, simplify = TRUE)[1,], lty = 2, col = 2)
    lines(seq(0,0.12,0.001), sapply(seq(0,0.12,0.001)*1000, fNq_endo_oliv09, simplify = TRUE)[2,], lty = 3, col = 2)
    lines(seq(0,0.12,0.001), sapply(seq(0,0.12,0.001)*1000, fNq_endo_oliv09, simplify = TRUE)[3,], lty = 3, col = 2)

  set.seed(43093)
  
    points(seq(0,0.12,0.001), sapply(seq(0,0.12,0.001)*1000, fNq_endo_oliv09_uncertainty, simplify = TRUE)*oliv_dat_ref, 
           pch = 5, col = 4, cex = 0.5)
  
  legend("topright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)
dev.off()