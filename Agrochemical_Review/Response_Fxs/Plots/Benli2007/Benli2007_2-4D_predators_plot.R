source("Agrochemical_Review/Response_Fxs/Benli2007_2-4D_predators_fit.R")

png("Agrochemical_Review/Response_Fxs/Plots/Benli2007/benli2007_2-4D_predator_mortality.png")

plot(benli_lcs, benli_morts, ylim = c(0,1), xlim = c(0,400), pch = 16,
     xlab = "2,4 D concentration (ppm)", ylab = "predator mortality", main = "Benli et al 2007")
  segments(x0 = 15.1, x1 = 327.16, y0=0.5, y1=0.5)
  lines(seq(0, 400, 2), sapply(seq(0,400,2)*1000, muPq_24d_benli07, simplify = TRUE), lty = 2, col = 2)
  lines(seq(0, 400, 2), sapply(seq(0,400,2)*1000, muPq_24d_benli07, lc50 = 327.16, simplify = TRUE), lty = 2, col = 2)
  lines(seq(0, 400, 2), sapply(seq(0,400,2)*1000, muPq_24d_benli07, lc50 = 15.1, simplify = TRUE), lty = 2, col = 2)
  
set.seed(43093)
  
  points(seq(0, 400, 2), sapply(seq(0,400,2)*1000, muPq_24d_benli07_uncertainty, simplify = TRUE),
         pch = 5, cex = 0.5, col = 4)

dev.off()  