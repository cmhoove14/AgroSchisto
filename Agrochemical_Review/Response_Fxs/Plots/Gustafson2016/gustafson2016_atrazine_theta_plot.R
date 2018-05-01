source("Agrochemical_Review/Response_Fxs/Gustafson2016_atrazine_theta_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Gustafson2016/atrazine_theta_data_sim.png")

plot(gust_dat$atr, gust_dat$cerc, pch = 16, xlab = "Atrazine (ppb)", ylab = "cercariae produced (theta)", ylim = c(0,1000))
  lines(c(0:30), gust_theta_pred[,1], lty = 2, col = 2)
  lines(c(0:30), gust_theta_pred[,2], lty = 3, col = 2)
  lines(c(0:30), gust_theta_pred[,3], lty = 3, col = 2)
  
  set.seed(43093)
  
  points(seq(0,30,0.1), sapply(seq(0,30,0.1), gust_theta_uncertainty, simplify = TRUE)*gust_int, pch = 5, cex = 0.5, col = 4)
  
dev.off()  