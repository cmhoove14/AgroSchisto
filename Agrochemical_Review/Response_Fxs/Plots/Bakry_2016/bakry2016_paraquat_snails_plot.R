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
source("Agrochemical_Review/Response_Fxs/bakry2016_paraquat_snails_fit.R")

#direct snail toxicity from Atrazine ############  
mun.bak.prq = data.frame(conc = c(0, 1.1, 2.3, 4.4)*1000,
                         mort = c(0, 0.25, 0.5, 0.9)) 
  mun.bak.prq$ppm = mun.bak.prq$conc/1000
  mun.bak.prq$log10 = log10(mun.bak.prq$ppm)
  mun.bak.prq$probit = qnorm(mun.bak.prq$mort, mean = 5)

png("Agrochemical_Review/Response_Fxs/Plots/Bakry_2016/bakry2016_paraquat_snail_mortality.png")  
  plot(mun.bak.prq$conc, mun.bak.prq$mort, ylim = c(0,1), xlim = c(0,5000),
       pch = 16, xlab = 'paraquat (ppb)', ylab = expression(italic(mu[N])),
       main = 'Bakry et al 2016 paraquat snail mortality')
  #segments(x0 = 1100, x1 = 4400, y0 = 0.25, y1 = 0.9) #Reported LC values look suspiciously linear...
  set.seed(43093)
  
  points(seq(0,5000,10), sapply(seq(0,5000,10), muNq_prq_Bakry16_uncertainty), 
         pch = 5, col = 4, cex = 0.5)

  dev.off()