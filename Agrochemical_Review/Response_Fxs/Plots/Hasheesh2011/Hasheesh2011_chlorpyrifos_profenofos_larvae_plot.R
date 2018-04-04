#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
source("Agrochemical_Review/Response_Fxs/Hasheesh2011_chlorpyrifos_profenofos_larvae_fit.R")
   
#Toxicity to miracidia from table 5 ###############
#ChlorP miracidia #########
piM.hash.ch = data.frame(conc = c(0, 0.78, 1.8)*1000,
                        mort = c(0, .50 , .90),
                        surv = 0)
  piM.hash.ch$surv = 1- piM.hash.ch$mort

  png("Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/hasheesh2011_chlorpyrifos_miracidia_mortality.png")
  
  plot(piM.hash.ch$conc, piM.hash.ch$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
       xlab = 'Chlorpyrifos (ppb)', ylab = '8-hr miracidial mortality',
       main = expression(paste(pi[M], ' estimate from Hasheesh 2011')))
    segments(x0 = 560, y0 = 0.5, x1 = 1100, y1 = 0.5)
    
      set.seed(43093)
      
   points(seq(0,4000,10), sapply(seq(0,4000,10), piM_ch_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)

      legend("bottomleft", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

  dev.off()    

#Profenofos miracidia #########
piM.hash.prof = data.frame(conc = c(0, 1.5, 2.51)*1000,
                           mort = c(0, .50 , .90),
                           surv = 0)
  piM.hash.prof$surv = 1- piM.hash.prof$mort
  
png("Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/hasheesh2011_profenofos_miracidia_mortality.png")
  
  plot(piM.hash.prof$conc, piM.hash.prof$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
       xlab = 'Profenofos (ppb)', ylab = '8-hr miracidial mortality',
       main = expression(paste(pi[M], ' estimate from Hasheesh 2011')))
    segments(x0 = 1150, y0 = 0.5, x1 = 1950, y1 = 0.5)
    
        set.seed(43093)
      
   points(seq(0,4000,10), sapply(seq(0,4000,10), piM_pr_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)

      legend("bottomleft", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  

#Toxicity to cercariae from table 5 ###############
#chlorpyrifos cercariae #######
piC.hash.ch = data.frame(conc = c(0, 0.96, 2.1)*1000,
                        mort = c(0,.50 , .90),
                        surv = 0)
  piC.hash.ch$surv = 1 - piC.hash.ch$mort

png("Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/hasheesh2011_chlorpyrifos_cercariae_mortality.png")
  
  plot(piC.hash.ch$conc, piC.hash.ch$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
       xlab = 'Chlorpyrifos (ppb)', ylab = '8-hr cercarial mortality',
       main = expression(paste(pi[C], ' estimate from Hasheesh 2011')))
    segments(x0 = 620, y0 = 0.5, x1 = 1440, y1 = 0.5)
 
        set.seed(43093)
      
  points(seq(0,4000,10), sapply(seq(0,4000,10), piC_ch_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)
  
    legend("bottomleft", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  
  

#profenofos cercariae #######
piC.hash.prof = data.frame(conc = c(0, 1.85, 2.85)*1000,
                          mort = c(0, .50, .90),
                          surv = 0)
  piC.hash.prof$surv = 1 - piC.hash.prof$mort

png("Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/hasheesh2011_profenofos_cercariae_mortality.png")
  
  plot(piC.hash.prof$conc, piC.hash.prof$surv, pch = 16, ylim = c(0,1), xlim = c(0,3500),
       xlab = 'Profenofos (ppb)', ylab = '8-hr cercarial mortality',
       main = expression(paste(pi[C], ' estimate from Hasheesh 2011')))
    segments(x0 = 2590, y0 = 0.5, x1 = 1320, y1 = 0.5)

        set.seed(43093)
        
   points(seq(0,4000,10), sapply(seq(0,4000,10), piC_pr_Hash11_uncertainty, simplify = T), 
         pch = 5, col = 4, cex = 0.5)  
    
    legend("bottomleft", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

dev.off()  
