source("Agrochemical_Review/Response_Fxs/Bhavan2010_carbaryl_predator_feed_rate_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Bhavan2010/bhavan2010_carbaryl_feed_rate_data_sim.png")

plot(bhavan_dat$carbaryl, bhavan_dat$food_con / bhavan_dat$food_con[1], pch = 16, xlab = "Carbaryl (ppb)", ylab = "Consumption relative to control")

    lines(seq(0,17,0.1), sapply(seq(0,17,0.1), psi_q_carb_bhavan10, simplify = T)[1,],
          lty = 2, col = 2)
    lines(seq(0,17,0.1), sapply(seq(0,17,0.1), psi_q_carb_bhavan10, simplify = T)[2,],
          lty = 3, col = 2)
    lines(seq(0,17,0.1), sapply(seq(0,17,0.1), psi_q_carb_bhavan10, simplify = T)[3,],
          lty = 3, col = 2)
    
    set.seed(43093)
    
    points(seq(0,17,0.05), sapply(seq(0,17,0.05), psi_q_carb_bhavan10_uncertainty, simplify = T), 
           pch=5, col=4, cex = 0.5)   

dev.off()