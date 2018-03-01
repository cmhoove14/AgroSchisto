#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Data extraction and model fitting to Bakry 2011 [Bakry et al 2011](https://www.sciencedirect.com/science/article/pii/S0048357511001283) data
source("Agrochemical_Review/Response_Fxs/bakry2011_malathion_deltamethrin_snails_fit.R")

#Snail (H. duryi) toxicity ##########
#Malathion data extracted from table 1 #####################
mun.mal = data.frame(conc = c(0, .176, .480, .830, 1.760, 3.360)*1000,
                     mort = c(0, 0,   .10, .25, .50 , .90),
                     surv = 0)
  mun.mal$surv = 1 - mun.mal$mort
  mun.mal$ppm = mun.mal$conc/1000
  mun.mal$ppmlog10 = log10(mun.mal$ppm)
  mun.mal$probit = qnorm(mun.mal$mort, mean = 5)
  
png("Agrochemical_Review/Response_Fxs/Plots/Bakry_2011/bakry2011_malathion_snail_mortality.png")

plot(mun.mal$conc, mun.mal$mort, pch = 16, 
     xlab = 'Malathion (ppb)', ylab = expression(italic(mu[N])), ylim = c(0,1), xlim = c(0,5000),
     main = 'Bakry et al 2011 malathion snail mortality')
  segments(x0 = 3120, y0 = 0.5, x1 = 990, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,5000,25), sapply(seq(0,5000,25), muNq_mal_Bakry11_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)
  
dev.off()  

#Deltamethrin  #############    
mun.del = data.frame(conc = c(0, .482, 1.21, 2.034, 4.82, 7.26)*1000,
                     mort = c(0, 0,   .10, .25, .50 , .90),
                     surv = 0)
  mun.del$surv = 1 - mun.del$mort
  mun.del$ppm = mun.del$conc/1000
  mun.del$ppmlog10 = log10(mun.del$ppm)
  mun.del$probit = qnorm(mun.del$mort, mean = 5)

png("Agrochemical_Review/Response_Fxs/Plots/Bakry_2011/bakry2011_deltamethrin_snail_mortality.png")  

plot(mun.del$conc, mun.del$mort, pch = 16, ylim = c(0,1), xlim = c(0,10000),
     xlab = 'Deltamethrin (ppb)', ylab = expression(italic(mu[N])), 
     main = 'Bakry et al 2011 deltamethrin snail mortality')
  segments(x0 = 3100, y0 = 0.5, x1 = 7700, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,10000,50), sapply(seq(0,10000,50), muNq_del_Bakry11_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)
  
dev.off()  
#Reproduction over time ###########
png("Agrochemical_Review/Response_Fxs/Plots/Bakry_2011/bakry2011_snail_reproduction.png")  
  
  plot(bakry11_ctrl$time, bakry11_ctrl$eggs, type = 'b', pch = 16, ylim = c(0,max(bakry11_ctrl$eggs)), xlim = c(0,30),
       xlab = 'time (days)', ylab = "eggs laid", 
       main = 'Snail reproduction over time')
    lines(bakry11_mal$time, bakry11_mal$eggs, type = 'b', pch = 16, col = 2)
    lines(bakry11_del$time, bakry11_del$eggs, type = 'b', pch = 16, col = 4)

  legend("topright", bty="n", lty = rep(1,3), col = c(1,2,4), legend = c("Control", "Malathion", "Deltamethrin"), cex = 0.75)
  
dev.off()  