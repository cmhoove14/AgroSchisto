require(ggplot2)
require(drc)

#24-hr mortality functions for predator (macrobrachium) populations from Satapornvanit et al chemosphere paper ########
  sap.mort<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_mort.csv")
    sap.mort$dead = round(sap.mort$dead)
  mort.sub = subset(sap.mort, chem !='carbendazim') 
  
  sap.mupq<-drm(dead/total ~ conc, chem, weights = total,  data = mort.sub, fct = LL.2(), type = 'binomial')
    summary(sap.mupq)
    plot(sap.mupq)
    
#Zinc ########    
  muPq_zinc_satapornvanit09<-function(In){
    predict(sap.mupq, data.frame(conc=In, chem = 'zinc'))
  }  
    
    zinc.sat09.df = data.frame(conc = c(0:1000),
                               chem = 'zinc',
                               Prediction = 0,
                               Lower = 0,
                               Upper = 0)
    
    zinc.sat09.df[,3:5] <- predict(sap.mupq, newdata = zinc.sat09.df, 
                                   interval = 'confidence', level = 0.95)
    
  plot(sap.mort$conc[sap.mort$chem == 'zinc'], sap.mort$mort[sap.mort$chem == 'zinc']/100, ylim = c(0,1),
       pch = 16, xlab = 'Zinc concentration (ppb)', ylab = 'prop dead', main = 'Zinc daily toxicity to M. rosenbergii')
    lines(zinc.sat09.df$conc, zinc.sat09.df$Prediction, lty=2, col=2)
    lines(zinc.sat09.df$conc, zinc.sat09.df$Lower, lty=3, col=2)
    lines(zinc.sat09.df$conc, zinc.sat09.df$Upper, lty=3, col=2)
    
#Chlorpyrifos ###########    
  muPq_chlor_satapornvanit09<-function(In){
    predict(sap.mupq, data.frame(conc=In, chem = 'chlorpyrifos'))
  }  
  
      chlor.sat09.df = data.frame(conc = seq(0,10,0.01),
                                 chem = 'chlorpyrifos',
                                 Prediction = 0,
                                 Lower = 0,
                                 Upper = 0)
    
    chlor.sat09.df[,3:5] <- predict(sap.mupq, newdata = chlor.sat09.df, 
                                   interval = 'confidence', level = 0.95)
    
    plot(sap.mort$conc[sap.mort$chem == 'chlorpyrifos'], sap.mort$mort[sap.mort$chem == 'chlorpyrifos']/100, ylim = c(0,1),
         pch = 16, xlab = 'chlorpyrifos concentration (ppb)', ylab = 'prop dead', main = 'chlorpyrifos daily toxicity to M. rosenbergii')
      lines(chlor.sat09.df$conc, chlor.sat09.df$Prediction, lty=2, col=2)
      lines(chlor.sat09.df$conc, chlor.sat09.df$Lower, lty=3, col=2)
      lines(chlor.sat09.df$conc, chlor.sat09.df$Upper, lty=3, col=2)

#Dimethoate ################        
  muPq_dim_satapornvanit09<-function(In){
    predict(sap.mupq, data.frame(conc=In, chem = 'dimethoate'))
  }  
      
    dim.sat09.df = data.frame(conc = c(0:1300),
                                  chem = 'dimethoate',
                                  Prediction = 0,
                                  Lower = 0,
                                  Upper = 0)
      
      dim.sat09.df[,3:5] <- predict(sap.mupq, newdata = dim.sat09.df, 
                                      interval = 'confidence', level = 0.95)
      
    plot(sap.mort$conc[sap.mort$chem == 'dimethoate'], sap.mort$mort[sap.mort$chem == 'dimethoate']/100, ylim = c(0,1),
         pch = 16, xlab = 'dimethoate concentration (ppb)', ylab = 'prop dead', main = 'dimethoate daily toxicity to M. rosenbergii')
      lines(dim.sat09.df$conc, dim.sat09.df$Prediction, lty=2, col=2)
      lines(dim.sat09.df$conc, dim.sat09.df$Lower, lty=3, col=2)
      lines(dim.sat09.df$conc, dim.sat09.df$Upper, lty=3, col=2)

#Profenofos ################      
  muPq_pr_satapornvanit09<-function(In){
    predict(sap.mupq, data.frame(conc=In, chem = 'profenofos'))
  }  
      
    prof.sat09.df = data.frame(conc = seq(0,60,0.1),
                                chem = 'profenofos',
                                Prediction = 0,
                                Lower = 0,
                                Upper = 0)
      
      prof.sat09.df[,3:5] <- predict(sap.mupq, newdata = prof.sat09.df, 
                                    interval = 'confidence', level = 0.95)
      
    plot(sap.mort$conc[sap.mort$chem == 'profenofos'], sap.mort$mort[sap.mort$chem == 'profenofos']/100, ylim = c(0,1),
         pch = 16, xlab = 'profenofos concentration (ppb)', ylab = 'prop dead', main = 'profenofos daily toxicity to M. rosenbergii')
      lines(prof.sat09.df$conc, prof.sat09.df$Prediction, lty=2, col=2)
      lines(prof.sat09.df$conc, prof.sat09.df$Lower, lty=3, col=2)
      lines(prof.sat09.df$conc, prof.sat09.df$Upper, lty=3, col=2)
  
  muPq_carb_satapornvanit09<-function(In){ #Paper found no consistent effect of carbendazim on mortality, even
      0*In                                 #at levels higher than the solubility of carbendazim in water
  }  
  

#Reduced feeding rate from zinc and Chlorpyrifos from Satapornvanit et al chemosphere paper ########
sap.fr<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_feed_rate.csv")
  fr.z = subset(sap.fr, chem == 'zinc')
  fr.ch = subset(sap.fr, chem == 'chlorpyrifos')
  
#Zinc ###################
  zinc.fr= drm(feed_rate ~ conc, data = fr.z, type = 'continuous',
               fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                          fixed = c(NA, 0, max(fr.z$feed_rate), NA)))
    summary(zinc.fr)
    plot(zinc.fr)
    
  psi_q_zinc_satapornvanit09<-function(In){
    predict(zinc.fr, data.frame(conc = In)) / predict(zinc.fr, data.frame(conc = 0))
  }  
    
  plot(fr.z$conc, fr.z$feed_rate / fr.z$feed_rate[1], pch = 16, ylim = c(0,1),
       xlab = 'Zinc concentration (ppb)', ylab = 'Feeding rate (prey/prawn/hr)', main = 'Prawn reduced feeding (Sata-09)')
  
  fr.z.df = data.frame(conc = c(0:1100),
                       Prediction = 0,
                       Lower = 0,
                       Upper = 0)
  
  fr.z.df[,2:4] <- predict(zinc.fr, newdata = fr.z.df, 
                           interval = 'confidence', level = 0.95)
  
  
  lines(fr.z.df$conc, fr.z.df$Prediction / fr.z.df$Prediction[1], col = 2, lty=2)
  lines(fr.z.df$conc, fr.z.df$Lower / fr.z.df$Prediction[1], col = 2, lty=3)
  lines(fr.z.df$conc, fr.z.df$Upper / fr.z.df$Prediction[1], col = 2, lty=3)  
  
#Chlorpyrifos ###################
  chlor.fr= drm(feed_rate ~ conc, data = fr.ch, type = 'continuous',
               fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                          fixed = c(NA, 0, max(fr.ch$feed_rate), NA)))
    summary(chlor.fr)
    plot(chlor.fr)
  
  psi_q_chlor_satapornvanit09<-function(In){
    predict(chlor.fr, data.frame(conc = In)) / predict(chlor.fr, data.frame(conc = 0))
  }  
  
  plot(fr.ch$conc, fr.ch$feed_rate / fr.ch$feed_rate[1], pch = 16, ylim = c(0,1),
       xlab = 'Chlorpyrifos concentration (ppb)', ylab = 'Feeding rate (prey/prawn/hr)', main = 'Prawn reduced feeding (Sata-09)')
  
  fr.ch.df = data.frame(conc = seq(0,6,0.01),
                       Prediction = 0,
                       Lower = 0,
                       Upper = 0)
  
  fr.ch.df[,2:4] <- predict(chlor.fr, newdata = fr.ch.df, 
                           interval = 'confidence', level = 0.95)
  
  
  lines(fr.ch.df$conc, fr.ch.df$Prediction / fr.ch.df$Prediction[1], col = 2, lty=2)
  lines(fr.ch.df$conc, fr.ch.df$Lower / fr.ch.df$Prediction[1], col = 2, lty=3)
  lines(fr.ch.df$conc, fr.ch.df$Upper / fr.ch.df$Prediction[1], col = 2, lty=3)  
    
#*Note: checked that this produced the desired relationship on predator feeding rate and it does