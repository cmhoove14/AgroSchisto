#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/Halstead2015_insecticides_predators_fit.R")

  predict_LL2.2 <- function(b, lc50, lc50_se, unc = 0, In){
    if(unc == 1){
      e = exp(rnorm(1, lc50, lc50_se))
    }
    else{
      e = lc50
    }
    
    mort <- 1 / (1 + exp(b*(log(In)-e)))
    
    mort
  }
  
#Malathion ******************************************************************************########
  png("Agrochemical_Review/Response_Fxs/Plots/Halstead2015/halstead2015_malathion_predator_mortality.png")
  
     plot(mal$Conc, mal$nrx/mal$ntot, pch = 16, ylim = c(0,1), xlim = c(0, 100000),
          xlab = 'Malathion', ylab = 'Mortality', 
          main = expression(paste('Malathion toxicity to ', italic('P. clarkii')))) 
     
      lines(seq(0, 100000, 100), 
            sapply(seq(0, 100000, 100), predict_LL2.2, b = halstead_mal_b, lc50 = halstead_mal_lc50, lc50_se = halstead_mal_lc50_se, unc = 0), 
            col = 2, lty = 2)
      
      set.seed(43093)
      
      points(seq(0, 100000, 500), 
             sapply(seq(0, 100000, 500), muPq_mal_Halstead_uncertainty),
             pch = 5, col = 4, cex = 0.5)
      
      legend("topleft", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

  dev.off()    
#Chlorpyrifos   *************************************************************************############
   png("Agrochemical_Review/Response_Fxs/Plots/Halstead2015/halstead2015_chlorpyrifos_predator_mortality.png")
  
     plot(chlor$Conc, chlor$nrx/chlor$ntot, pch = 16, ylim = c(0,1), xlim = c(0, 100),
           xlab = 'Chlorpyrifos', ylab = 'Mortality', 
           main = expression(paste('Chlorpyrifos toxicity to ', italic('P. clarkii')))) 
      
      lines(seq(0, 100, 1), 
            sapply(seq(0, 100, 1), predict_LL2.2, b = halstead_chlor_b, lc50 = halstead_chlor_lc50, lc50_se = halstead_chlor_lc50_se, unc = 0), 
            col = 2, lty = 2)
      
      set.seed(43093)      
      
      points(seq(0, 100, 0.5), 
             sapply(seq(0, 100, 0.5), muPq_chlor_Halstead_uncertainty),
             pch = 5, col = 4, cex = 0.5)
      
      legend("right", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

    dev.off()  
    
#Terbufos   *****************************************************************************########
  png("Agrochemical_Review/Response_Fxs/Plots/Halstead2015/halstead2015_terbufos_predator_mortality.png")
    
      plot(terb$Conc, terb$nrx/terb$ntot, pch = 16, ylim = c(0,1), xlim = c(0, 200),
           xlab = 'terbufos', ylab = 'Mortality', 
           main = expression(paste('terbufos toxicity to ', italic('P. clarkii')))) 
      
      lines(seq(0, 200, 1), 
            sapply(seq(0, 200, 1), predict_LL2.2, b = halstead_terb_b, lc50 = halstead_terb_lc50, lc50_se = halstead_terb_lc50_se, unc = 0), 
            col = 2, lty = 2)
      
      set.seed(43093)

      points(seq(0, 200, 0.5), 
             sapply(seq(0, 200, 0.5), muPq_terb_Halstead_uncertainty),
             pch = 5, col = 4, cex = 0.5)
            
      legend("right", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)


  dev.off()      

#Lambda-cyhalothrin   *******************************************************************##########
  png("Agrochemical_Review/Response_Fxs/Plots/Halstead2015/halstead2015_lambda-cyhalothrin_predator_mortality.png")
  
      plot(lamcy$Conc, lamcy$nrx/lamcy$ntot, pch = 16, ylim = c(0,1), xlim = c(0, 3),
           xlab = expression(paste(lambda,'-cyhalothrin (ppb)')), ylab = 'Mortality', 
           main = expression(paste(lambda,'-cyhalothrin toxicity to ', italic('P. clarkii')))) 
      
      lines(seq(0, 3, 0.01), 
            sapply(seq(0, 3, 0.01)*1000, predict_LL2.2, b = halstead_lamcy_b, lc50 = halstead_lamcy_lc50, lc50_se = halstead_lamcy_lc50_se, unc = 0), 
            col = 2, lty = 2)
      
      set.seed(43093)
      
      points(seq(0, 3, 0.05), 
             sapply(seq(0, 3, 0.05), muPq_lamcy_Halstead_uncertainty),
             pch = 5, col = 4, cex = 0.5)

      legend("right", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

    dev.off()  
#esfenvalerate  *************************************************************************#######
  png("Agrochemical_Review/Response_Fxs/Plots/Halstead2015/halstead2015_esfenvalerate_predator_mortality.png")
      
      plot(esfen$Conc, esfen$nrx/esfen$ntot, pch = 16, ylim = c(0,1), xlim = c(0, 25),
           xlab = 'esfenvalerate', ylab = 'Mortality', 
           main = expression(paste('esfenvalerate toxicity to ', italic('P. clarkii')))) 
      
      lines(seq(0, 25, 0.1), 
            sapply(seq(0, 25, 0.1)*1000, predict_LL2.2, b = halstead_esfen_b, lc50 = halstead_esfen_lc50, lc50_se = halstead_esfen_lc50_se, unc = 0), 
            col = 2, lty = 2)
      
      set.seed(43093)
      
      points(seq(0, 25, 0.5), 
             sapply(seq(0, 25, 0.5), muPq_esfen_Halstead_uncertainty),
             pch = 5, col = 4, cex = 0.5)

      legend("right", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

  dev.off() 
  
#Permethrin *****************************************************************************##########
  png("Agrochemical_Review/Response_Fxs/Plots/Halstead2015/halstead2015_permethrin_predator_mortality.png")
      
      plot(perm$Conc, perm$nrx/perm$ntot, pch = 16, ylim = c(0,1), xlim = c(0, 6),
           xlab = 'permethrin', ylab = 'Mortality', 
           main = expression(paste('permethrin toxicity to ', italic('P. clarkii')))) 
      
      lines(seq(0, 6, 0.1), 
            sapply(seq(0, 6, 0.1)*1000, predict_LL2.2, b = halstead_perm_b, lc50 = halstead_perm_lc50, lc50_se = halstead_perm_lc50_se, unc = 0), 
            col = 2, lty = 2)
      
    set.seed(43093)

      points(seq(0, 6, 0.05), 
             sapply(seq(0, 6, 0.05), muPq_perm_Halstead_uncertainty),
             pch = 5, col = 4, cex = 0.5)

      legend("right", bty="n", pch = c(16, 5), col = c(1,4), legend = c("Reported", "Sampled"), cex = 0.75)

    dev.off()