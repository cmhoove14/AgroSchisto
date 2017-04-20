require(ggplot2)
require(drc)

#Insecticide toxicity to crustaceans from Halstead et al; 4-day mortality endpoints ################
  data<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/Cray.LC50.2012.csv')
  
#Malathion ********************************************************************************************
  mal<-subset(data, Chem == 'Mal')
    mal.sum = data.frame(mal.c = unique(mal$Conc),
                         mal.total = 5,
                         mal.d = 0,
                         mort = 0)
    
    for(i in 1:length(mal.sum$mal.d)){
      mal.sum$mal.d[i] = sum(mal$Dead[mal$Conc == mal.sum$mal.c[i]])
    }
    
    mal.sum$mort = mal.sum$mal.d / mal.sum$mal.total
    
      mal.mod<-drm(mal.d / mal.total ~ mal.c, weights = mal.total, data = mal.sum, type = 'binomial',  
                   fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                              fixed = c(NA, 0, 1, NA)))

     muPq_mal_Halstead<-function(In){
        predict(mal.mod, data.frame(mal.c = In), interval = 'confidence', level = 0.95)
      }  
     
     plot(mal.sum$mal.c, mal.sum$mort, pch = 16, ylim = c(0,1),
          xlab = 'Malathion', ylab = 'Mortality', 
          main = expression(paste('Malathion toxicity to ', italic('P. clarkii')))) 
     
     
      lines(seq(0, max(mal$Conc), 100), sapply(seq(0, max(mal$Conc), 100), muPq_mal_Halstead, simplify = T)[1,],
            lty = 2, col = 2)
      lines(seq(0, max(mal$Conc), 100), sapply(seq(0, max(mal$Conc), 100), muPq_mal_Halstead, simplify = T)[2,],
            lty = 3, col = 2)
      lines(seq(0, max(mal$Conc), 100), sapply(seq(0, max(mal$Conc), 100), muPq_mal_Halstead, simplify = T)[3,],
            lty = 3, col = 2)
        
      muPq_mal_Halstead_uncertainty<-function(In){
        rdrm(1, LL.2(), coef(mal.mod), In, yerror = 'rbinom', ypar = 5)$y / 5 #estimate deaths / live 
      }
        points(seq(0, 40000, 40), sapply(seq(0, 40000, 40), muPq_mal_Halstead_uncertainty),
             pch = 5, col=4, cex = 0.5)
    
      
#Chlorpyrifos   ********************************************************************************************
  chlor<-subset(data, Chem == 'Chlor')
    chlor.sum = data.frame(chlor.c = unique(chlor$Conc),
                         chlor.total = 5,
                         chlor.d = 0,
                         mort = 0)
      
    for(i in 1:length(chlor.sum$chlor.d)){
        chlor.sum$chlor.d[i] = sum(chlor$Dead[chlor$Conc == chlor.sum$chlor.c[i]])
    }
      
    chlor.sum$mort = chlor.sum$chlor.d / chlor.sum$chlor.total
      
    chlor.mod<-drm(chlor.d / chlor.total ~ chlor.c, weights = chlor.total, data = chlor.sum, type = 'binomial',  
                   fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                              fixed = c(NA, 0, 1, NA)))

    muPq_chlor_Halstead<-function(In){
      predict(chlor.mod, data.frame(chlor.c = In), interval = 'confidence', level = 0.95)
    }  
      
    plot(chlor.sum$chlor.c, chlor.sum$mort, pch = 16, ylim = c(0,1),
         xlab = 'Chlorpyrifos', ylab = 'Mortality', 
         main = expression(paste('Chlorpyrifos toxicity to ', italic('P. clarkii')))) 
    
    lines(seq(0, max(chlor$Conc), 1), sapply(seq(0, max(chlor$Conc), 1), muPq_chlor_Halstead, simplify = T)[1,],
          lty = 2, col = 2)
    lines(seq(0, max(chlor$Conc), 1), sapply(seq(0, max(chlor$Conc), 1), muPq_chlor_Halstead, simplify = T)[2,],
          lty = 3, col = 2)
    lines(seq(0, max(chlor$Conc), 1), sapply(seq(0, max(chlor$Conc), 1), muPq_chlor_Halstead, simplify = T)[3,],
          lty = 3, col = 2)
      
      muPq_chlor_Halstead_uncertainty<-function(In){
        rdrm(1, LL.2(), coef(chlor.mod), In, yerror = 'rbinom', ypar = 5)$y / 5
      } 
        
      points(c(0:100), sapply(c(0:100), muPq_chlor_Halstead_uncertainty),
           pch = 5, col=4, cex = 0.5)
      
#Terbufos   ********************************************************************************************
  terb<-subset(data, Chem == 'Terb')
    terb.sum = data.frame(terb.c = unique(terb$Conc),
                             terb.total = 5,
                             terb.d = 0,
                             mort = 0)
      
      for(i in 1:length(terb.sum$terb.d)){
        terb.sum$terb.d[i] = sum(terb$Dead[terb$Conc == terb.sum$terb.c[i]])
      }
      
      terb.sum$mort = terb.sum$terb.d / terb.sum$terb.total
      
  terb.mod<-drm(terb.d / terb.total ~ terb.c, weights = terb.total, data = terb.sum, type = 'binomial', 
                fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                           fixed = c(NA, 0, 1, NA)))
#    summary(terb.mod)
#    plot(terb.mod, ylim = c(0,1))
      
#      muPq_terb_Halstead<-function(In){
#        predict(terb.mod, data.frame(terb.c = In))
#      }  
      
#      terb.hal.df = data.frame(terb.c = seq(0,200,0.1),
#                                Prediction = 0,
#                                Lower = 0,
#                                Upper = 0)
      
#      terb.hal.df[,2:4] <- predict(terb.mod, newdata = terb.hal.df, 
#                                    interval = 'confidence', level = 0.95)
      
#    plot(terb.sum$terb.c, terb.sum$mort, pch = 16, xlab = 'Terbufos Concentration (ppb)', 
#         ylab = 'prop dead', ylim = c(0,1))
#      lines(terb.hal.df$terb.c, terb.hal.df$Prediction, col = 2, lty=2)
#      lines(terb.hal.df$terb.c, terb.hal.df$Lower, col = 2, lty=3)
#      lines(terb.hal.df$terb.c, terb.hal.df$Upper, col = 2, lty=3)
      
    muPq_terb_Halstead_uncertainty<-function(In){
      rdrm(1, LL.2(), coef(terb.mod), In, yerror = 'rbinom', ypar = 5)$y / 5
    } 
    
    plot(seq(0, 50, 0.1), sapply(seq(0, 50, 0.1), muPq_terb_Halstead_uncertainty),
         pch = 17, col=2, cex = 0.5, xlab = 'Terbufos (ppb)', ylab = expression(paste(mu[P])),
         main = 'Sample output of predator mortality')
        
#Lambda-cyhalothrin   ********************************************************************************************
  lamcy<-subset(data, Chem == 'Lambda')
  lamcy.sum = data.frame(lamcy.c = unique(lamcy$Conc),
                        lamcy.total = 5,
                        lamcy.d = 0,
                        mort = 0)
  
    for(i in 1:length(lamcy.sum$lamcy.d)){
      lamcy.sum$lamcy.d[i] = sum(lamcy$Dead[lamcy$Conc == lamcy.sum$lamcy.c[i]])
    }
    
    lamcy.sum$mort = lamcy.sum$lamcy.d / lamcy.sum$lamcy.total
  
  lamcy.mod<-drm(lamcy.d / lamcy.total ~ lamcy.c, weights = lamcy.total, data = lamcy.sum, type = 'binomial', 
                fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                           fixed = c(NA, 0, 1, NA)))
#    summary(lamcy.mod)
#    plot(lamcy.mod, ylim = c(0,1))
  
#    muPq_lamcy_Halstead<-function(In){
#      predict(lamcy.mod, data.frame(lamcy.c = In))
#    }  
  
#  lamcy.hal.df = data.frame(lamcy.c = seq(0,5,0.01),
#                           Prediction = 0,
#                           Lower = 0,
#                           Upper = 0)
  
#    lamcy.hal.df[,2:4] <- predict(lamcy.mod, newdata = lamcy.hal.df, 
#                                 interval = 'confidence', level = 0.95)
  
#  plot(lamcy.sum$lamcy.c, lamcy.sum$mort, pch = 16, xlab = 'Lambda-Cyhalothrin Concentration (ppb)', 
#       ylab = 'prop dead', ylim = c(0,1))
#    lines(lamcy.hal.df$lamcy.c, lamcy.hal.df$Prediction, col = 2, lty=2)
#    lines(lamcy.hal.df$lamcy.c, lamcy.hal.df$Lower, col = 2, lty=3)
#    lines(lamcy.hal.df$lamcy.c, lamcy.hal.df$Upper, col = 2, lty=3)
  
  muPq_lamcy_Halstead_uncertainty<-function(In){
    rdrm(1, LL.2(), coef(lamcy.mod), In, yerror = 'rbinom', ypar = 5)$y / 5
  } 
  
  plot(seq(0, 1, 0.001), sapply(seq(0, 1, 0.001), muPq_lamcy_Halstead_uncertainty),
       pch = 17, col=2, cex = 0.5, 
       xlab = expression(paste(lambda, '-cyhalothrin (ppb)')), ylab = expression(paste(mu[P])),
       main = 'Sample output of predator mortality')
    
#esfenvalerate  ********************************************************************************************
  esfen<-subset(data, Chem == 'Esfen')
    esfen.sum = data.frame(esfen.c = unique(esfen$Conc),
                           esfen.total = 5,
                           esfen.d = 0,
                           mort = 0)
  
      for(i in 1:length(esfen.sum$esfen.d)){
        esfen.sum$esfen.d[i] = sum(esfen$Dead[esfen$Conc == esfen.sum$esfen.c[i]])
      }
      
      esfen.sum$mort = esfen.sum$esfen.d / esfen.sum$esfen.total
  
  esfen.mod<-drm(esfen.d / esfen.total ~ esfen.c, weights = esfen.total, data = esfen.sum, type = 'binomial', 
                 fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                            fixed = c(NA, 0, 1, NA)))
#    summary(esfen.mod)
#    plot(esfen.mod, ylim = c(0,1))
  
#    muPq_esfen_Halstead<-function(In){
#      predict(esfen.mod, data.frame(esfen.c = In))
#    }  
  
#  esfen.hal.df = data.frame(esfen.c = seq(0,50,0.1),
#                            Prediction = 0,
#                            Lower = 0,
#                            Upper = 0)
  
#    esfen.hal.df[,2:4] <- predict(esfen.mod, newdata = esfen.hal.df, 
#                                  interval = 'confidence', level = 0.95)
#  
#  plot(esfen.sum$esfen.c, esfen.sum$mort, pch = 16, xlab = 'Esfenvalerate Concentration (ppb)', 
#       ylab = 'prop dead', ylim = c(0,1))
#    lines(esfen.hal.df$esfen.c, esfen.hal.df$Prediction, col = 2, lty=2)
#    lines(esfen.hal.df$esfen.c, esfen.hal.df$Lower, col = 2, lty=3)
#    lines(esfen.hal.df$esfen.c, esfen.hal.df$Upper, col = 2, lty=3)
    
  muPq_esfen_Halstead_uncertainty<-function(In){
    rdrm(1, LL.2(), coef(esfen.mod), In, yerror = 'rbinom', ypar = 5)$y / 5 
  } 
  
  plot(seq(0, 1, 0.001), sapply(seq(0, 1, 0.001), muPq_esfen_Halstead_uncertainty),
       pch = 17, col=2, cex = 0.5, xlab = 'esfenvalerate', ylab = expression(paste(mu[P])),
       main = 'Sample output of predator mortality')
  
#Permethrin ********************************************************************************************
  perm<-subset(data, Chem == 'Perm')
    perm.sum = data.frame(perm.c = unique(perm$Conc),
                           perm.total = 5,
                           perm.d = 0,
                           mort = 0)
  
      for(i in 1:length(perm.sum$perm.d)){
        perm.sum$perm.d[i] = sum(perm$Dead[perm$Conc == perm.sum$perm.c[i]])
      }
      
      perm.sum$mort = perm.sum$perm.d / perm.sum$perm.total
  
  perm.mod<-drm(perm.d / perm.total ~ perm.c, weights = perm.total, data = perm.sum, type = 'binomial', 
                 fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                            fixed = c(NA, 0, 1, NA)))
#    summary(perm.mod)
#    plot(perm.mod, ylim = c(0,1))
  
#  muPq_perm_Halstead<-function(In){
#    predict(perm.mod, data.frame(perm.c = In))
#  }  
  
#  perm.hal.df = data.frame(perm.c = seq(0,10,0.01),
#                            Prediction = 0,
#                            Lower = 0,
#                            Upper = 0)
  
#    perm.hal.df[,2:4] <- predict(perm.mod, newdata = perm.hal.df, 
#                                  interval = 'confidence', level = 0.95)
  
#  plot(perm.sum$perm.c, perm.sum$mort, pch = 16, xlab = 'Permethrin Concentration (ppb)', 
#       ylab = 'prop dead', ylim = c(0,1))
#    lines(perm.hal.df$perm.c, perm.hal.df$Prediction, col = 2, lty=2)
#    lines(perm.hal.df$perm.c, perm.hal.df$Lower, col = 2, lty=3)
#    lines(perm.hal.df$perm.c, perm.hal.df$Upper, col = 2, lty=3)
  
  muPq_perm_Halstead_uncertainty<-function(In){
    rdrm(1, LL.2(), coef(perm.mod), In, yerror = 'rbinom', ypar = 5)$y / 5
  }
  
  plot(seq(0, 1, 0.001), sapply(seq(0, 1, 0.001), muPq_perm_Halstead_uncertainty),
       pch = 17, col=2, cex = 0.5, xlab = 'Permethrin (ppb)', ylab = expression(paste(mu[P])),
       main = 'Sample output of predator mortality')
  
#Insecticide info from Halstead et al chemosphere paper ################
#calculate decay rates (k) from hydrolysis half lives for each chemical in Halstead 2015 (from table S1 and pmep.cce.cornell.edu) 
  mal.k = -log(0.5)/6.2         #from table s1 and in agreement of "less than 1 week in raw river water" from Cornell
  chlor.k = -log(0.5)/25.5      #from table S1; within the range of reported half life from cornell
  terb.k = -log(0.5)/6.5        #from table s1; in agreement with Cornell estimate of 5.5 days at pH of 7; degrades into formaldehyde
  lamcy.k = -log(0.5)/1         #very fast; 0 in table S1; "Not expected to be prevalent in surface waters" according to cornell website
  esfen.k = -log(0.5)/10        #cornell 4-15 days half life in water
  perm.k = -log(0.5)/2          #cornell "half life of less than 2.5 days"
  
  #median EECs from Halstead 2015
  med.mal = 0.778
  med.chlor = 5.810
  med.terb = 1.435
  med.lamcy = 0.649
  med.esfen = 0.311
  med.perm = 1.420
  
  #suggested application intervals for each agrochemical from Halstead 2015
  mal.days = c(1,6)             #two applications 5 days apart
  chlor.days = c(1,11,21)       #three applications 10 days apart
  terb.days = c(1)              #single application
  lamcy.days = seq(1, 61, by=4) #16 applications 4 days apart
  esfen.days = seq(1, 21, by=5) #5 applications 5 days apart
  perm.days = seq(1, 36, by=5)  #8 applications 5 days apart
  

#Character vector of objects to keep from this script ###############
  keep.hal15.muP = c('mal.sum', 'chlor.sum', 'terb.sum', 'esfen.sum', 'lamcy.sum', 'perm.sum',
                     'mal.mod', 'chlor.mod', 'terb.mod', 'esfen.mod', 'lamcy.mod', 'perm.mod',
                     'muPq_mal_Halstead_uncertainty', 'muPq_chlor_Halstead_uncertainty', 'muPq_terb_Halstead_uncertainty',
                     'muPq_esfen_Halstead_uncertainty', 'muPq_lamcy_Halstead_uncertainty', 'muPq_perm_Halstead_uncertainty')