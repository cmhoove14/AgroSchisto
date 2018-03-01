#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Bakry 2012 data
source("Agrochemical_Review/Response_Fxs/bakry2012_atrazine_glyphosate_snails_fit.R")

#direct snail toxicity from Atrazine ############  
mun.bak.atr = data.frame(conc = c(.330, 1.250, 4.750)*1000,
                         mort = c(0.1, 0.5, 0.9)) 
  mun.bak.atr$ppm = mun.bak.atr$conc/1000
  mun.bak.atr$log10 = log10(mun.bak.atr$ppm)
  mun.bak.atr$probit = qnorm(mun.bak.atr$mort, mean = 5)

png("Agrochemical_Review/Response_Fxs/Plots/Bakry_2012/bakry2012_atrazine_snail_mortality.png")  
  
plot(mun.bak.atr$conc, mun.bak.atr$mort, ylim = c(0,1), xlim = c(0,5000),
       pch = 16, xlab = 'atrazine (ppb)', ylab = expression(italic(mu[N])),
     main = 'Bakry et al 2012 atrazine snail mortality')
    segments(y0 = 0.5, x0 = 830, y1 = 0.5, x1 = 1880)

  set.seed(43093)
    
  points(seq(0,5000,20), sapply(seq(0,5000,20), muNq_atr_Bakry12_uncertainty), 
       pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)
  
dev.off()

#direct snail toxicity from Glyphosate ############  
mun.bak.gly = data.frame(conc = c(.840, 3.150, 12.600)*1000,
                         mort = c(0.1, 0.5, 0.9))
  
  mun.bak.gly$ppm = mun.bak.gly$conc/1000
  mun.bak.gly$log10 = log10(mun.bak.gly$ppm)
  mun.bak.gly$probit = qnorm(mun.bak.gly$mort, mean = 5)

png("Agrochemical_Review/Response_Fxs/Plots/Bakry_2012/bakry2012_glyphosate_snail_mortality.png")  
  
  plot(mun.bak.gly$conc, mun.bak.gly$mort, ylim = c(0,1), xlim = c(0,13000),
       pch = 16, xlab = 'glyphosate (ppb)', ylab = expression(italic(mu[N])),
       main = 'Bakry et al 2012 glyphosate snail mortality')  
    segments(y0 = 0.5, x0 = 890, y1 = 0.5, x1 = 4820)
      
  set.seed(43093)
      
  points(seq(0,13000,50), sapply(seq(0,13000,50), muNq_gly_Bakry12_uncertainty), 
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)
  
dev.off()

#Reproduction over time ###########
bakry12<-read.csv('Agrochemical_Review/Response_Fxs/Data/bakry2012.csv')
  bakry12_ctrl<-subset(bakry12, chem == "control")
  bakry12_atr<-subset(bakry12, chem == "atrazine")
  bakry12_gly<-subset(bakry12, chem == "glyphosate")
  
png("Agrochemical_Review/Response_Fxs/Plots/Bakry_2012/bakry2012_snail_reproduction.png")  
  
plot(bakry12_ctrl$time, bakry12_ctrl$hatch, type = 'b', pch = 16, ylim = c(0,max(bakry12_ctrl$hatch)), xlim = c(0,42),
     xlab = 'time (days)', ylab = "eggs hatched", 
     main = 'Snail reproduction over time')
  lines(bakry12_atr$time, bakry12_atr$hatch, type = 'b', pch = 16, col = "gold")
  lines(bakry12_gly$time, bakry12_gly$hatch, type = 'b', pch = 16, col = "green")
  
  legend("topright", bty="n", lty = rep(1,3), col = c(1,"gold","green"), legend = c("Control", "Atrazine", "Glyphosate"), cex = 0.75)
  
dev.off()  