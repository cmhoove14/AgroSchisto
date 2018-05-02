source("Agrochemical_Review/Response_Fxs/Kristoff2011_azinphos-methyl_snail_reproduction_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Kristoff2011/kristoff2011_azinphos-methyl_snail_reproduction_data_sim.png")

plot(kristoff_dat$azmethyl, kristoff_dat$egg_mass / kristoff_ref, pch = 16, ylim = c(0,1),
     xlab = "Azinphos-methyl (ppb)", ylab = "Reproduction relative to control")

    lines(seq(0,5.1,0.1), sapply(seq(0,5.1,0.1)*1000, fNq_azmeth_kristoff11, simplify = T)[1,],
          lty = 2, col = 2)
    lines(seq(0,5.1,0.1), sapply(seq(0,5.1,0.1)*1000, fNq_azmeth_kristoff11, simplify = T)[2,],
          lty = 3, col = 2)
    lines(seq(0,5.1,0.1), sapply(seq(0,5.1,0.1)*1000, fNq_azmeth_kristoff11, simplify = T)[3,],
          lty = 3, col = 2)
    
    set.seed(43093)
    
    points(seq(0,5.1,0.05), sapply(seq(0,5.1,0.05)*1000, fNq_azmeth_kristoff11_uncertainty, simplify = T), 
           pch=5, col=4, cex = 0.5)   
    
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()