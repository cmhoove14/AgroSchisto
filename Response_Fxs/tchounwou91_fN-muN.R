#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(drc)

#Load data #############
sn = read.csv('C:/Users/chris_hoover/Google Drive/Remais_Group_Shared_Drive/EEID_Schisto/Data/Response_Functions/Extracted_Data/Snail Mortality/tchounwou1991.csv')
  sn$conc_ppb = sn$conc*1000
  sn$mort = 20 - sn$surv
eg = read.csv('C:/Users/chris_hoover/Google Drive/Remais_Group_Shared_Drive/EEID_Schisto/Data/Response_Functions/Extracted_Data/Snail reproduction/tchounwou1991.csv')
  eg$conc_ppb = eg$conc*1000
  
#Fit model to adult daily mortality ###########
plot(sn$conc_ppb, sn$prop_dead, pch = 16, xlim = c(0,5e5),
     xlab = 'Malathion (ppb)', ylab = 'Prop dead', main = 'Tchounwou 91 B. havenensis mortality data')

sn.mod = drm(mort/20 ~ conc_ppb, weights = rep(20, length(mort)), data = sn, type = 'binomial',
             fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                        fixed = c(NA, 0, 1, NA)))

tch.mal.df = data.frame(conc = seq(0, max(sn$conc_ppb), max(sn$conc_ppb)/10000),
                        Prediction = 0,
                        Lower = 0,
                        Upper = 0)

  tch.mal.df[,2:4] <- predict(sn.mod, newdata = tch.mal.df, 
                              interval = 'confidence', level = 0.95)
  
    lines(tch.mal.df$conc, tch.mal.df$Prediction, col = 2, lty=2)
    lines(tch.mal.df$conc, tch.mal.df$Lower, col = 2, lty=3)
    lines(tch.mal.df$conc, tch.mal.df$Upper, col = 2, lty=3)

muNq_mal_tch91_uncertainty<-function(In){
  rdrm(1, LL.2(), coef(sn.mod), In, yerror = 'rbinom', ypar = 100)$y / 100 #estimate deaths / live 
}

  points(seq(0, max(sn$conc_ppb), max(sn$conc_ppb)/1000), 
         sapply(seq(0, max(sn$conc_ppb), max(sn$conc_ppb)/1000), muNq_mal_tch91_uncertainty),
         pch = 17, col=4, cex = 0.5)
  legend('bottomright', legend = c('Observed points', 'Simulated points'),
         pch = c(16,17), col = c(1,4), cex = 0.6, bty='n')
  legend('right', legend = c('Prediction/ 95%CI'),
         lty = 2, col = 2, cex = 0.6, bty='n')
  
#Fit model to egg daily mortality ###########
plot(eg$conc_ppb, eg$prop_surv, pch = 16, xlim = c(0,5e5),
     xlab = 'Malathion (ppb)', ylab = 'Prop eggs dead', 
     main = 'Tchounwou 91 B. havenensis egg viability data')

eg.mod = drm(prop_surv ~ conc_ppb, data = eg, type = 'binomial',
             fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                        fixed = c(NA, 0, 1, NA)))
  
    tch.mal.df[,2:4] <- predict(eg.mod, newdata = tch.mal.df, 
                                interval = 'confidence', level = 0.95)
  
    lines(tch.mal.df$conc, tch.mal.df$Prediction, col = 2, lty=2)
      lines(tch.mal.df$conc, tch.mal.df$Lower, col = 2, lty=3)
      lines(tch.mal.df$conc, tch.mal.df$Upper, col = 2, lty=3)
  
  fNq_mal_tch91_uncertainty<-function(In){
    rdrm(1, LL.2(), coef(eg.mod), In, yerror = 'rbinom', ypar = 100)$y / 100 #estimate deaths / live 
  }
  
  points(seq(0, 5e5, 5e5/1000), 
         sapply(seq(0, 5e5, 5e5/1000), fNq_mal_tch91_uncertainty),
         pch = 17, col=4, cex = 0.5)  
  legend('topright', legend = c('Observed points', 'Simulated points'),
         pch = c(16,17), col = c(1,4), cex = 0.6, bty='n')
  legend('right', legend = c('Prediction/ 95%CI'),
         lty = 2, col = 2, cex = 0.6, bty='n')
#Keep vector
  keep.tch91.snail = c('sn', 'eg', 'sn.mod', 'eg.mod', 
                       'fNq_mal_tch91_uncertainty', 'muNq_mal_tch91_uncertainty')