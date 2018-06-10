source("Agrochemical_Review/Response_Fxs/Browne&Moore2014_2-4D_pred_consumption_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Browne2014/pred_consumption_2-4D_data_sim.png")

plot(browne_dat$conc*1000, browne_dat$per_con / browne_ref, pch = 16, ylim = c(0,1), xlim = c(0, 35000),
     xlab = "2,4-D concentration (ppb)", ylab = "relative predator consumption rate")

  lines(seq(0,35000,50), sapply(seq(0,35000,50), psiq_24D_browne14), lty = 2, col = 2)

  set.seed(43093)
  
  points(seq(0, 35000, 50), sapply(seq(0,35000,50), psiq_24D_browne14_uncertainty), 
         pch = 5, cex = 0.5, col = 4)
  
  legend("topright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()