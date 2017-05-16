require(ggplot2)
require(drc)

#24-hr mortality functions for predator (macrobrachium) populations from Satapornvanit et al chemosphere paper ########
  sap.mort<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_mort.csv")
    sap.mort$dead = round(sap.mort$dead)
  mort.sub = subset(sap.mort, chem !='carbendazim') 
  
  sap.mupq<-drm(dead/total ~ conc, chem, weights = total,  data = mort.sub, type = 'binomial', 
                fct = L.4(names = c('b', 'c', 'd', 'e'),
                          fixed = c(NA, 0, 1, NA)))
    summary(sap.mupq)
    plot(sap.mupq)
    
#Zinc ########    
  muPq_zinc_satapornvanit09<-function(In){
    predict(sap.mupq, data.frame(conc=In, chem = 'zinc'), interval = 'confidence', level = 0.95)
  }  
    
  plot(sap.mort$conc[sap.mort$chem == 'zinc'], sap.mort$mort[sap.mort$chem == 'zinc']/100, ylim = c(0,1),
       pch = 16, xlab = 'Zinc concentration (ppb)', ylab = 'prop dead', 
       main = expression(paste('Zinc daily toxicity to post-larval ', italic('M. rosenbergii'), sep = '')))
    lines(seq(0,900,10), sapply(seq(0,900,10), muPq_zinc_satapornvanit09, simplify = T)[1,], lty=2, col=2)
    lines(seq(0,900,10), sapply(seq(0,900,10), muPq_zinc_satapornvanit09, simplify = T)[2,], lty=3, col=2)
    lines(seq(0,900,10), sapply(seq(0,900,10), muPq_zinc_satapornvanit09, simplify = T)[3,], lty=3, col=2)
    
  muPq_zinc_satapornvanit09_uncertainty<-function(In){
    rdrm(1, L.4(), coef(sap.mupq)[c(1,5)], In, yerror = 'rbinom', ypar = 30)$y / 30
  }
    
  points(seq(0,900,10), sapply(seq(0,900,10), muPq_zinc_satapornvanit09_uncertainty, simplify = T), 
         pch=5, col=4, cex = 0.5)
  
keep.zinc.sat09 = c('muPq_zinc_satapornvanit09_uncertainty', 'sap.mupq')    
#Chlorpyrifos ###########    
  muPq_chlor_satapornvanit09<-function(In){
    predict(sap.mupq, data.frame(conc=In, chem = 'chlorpyrifos'), 
            interval = 'confidence', level = 0.95)
  }  
    
    plot(sap.mort$conc[sap.mort$chem == 'chlorpyrifos'], sap.mort$mort[sap.mort$chem == 'chlorpyrifos']/100, ylim = c(0,1),
         pch = 16, xlab = 'chlorpyrifos concentration (ppb)', ylab = 'prop dead', 
         main = expression(paste('Chlorpyrifos daily toxicity to post-larval ', italic('M. rosenbergii'), sep = '')))
      lines(seq(0,5,0.01), sapply(seq(0,5,0.01), muPq_chlor_satapornvanit09, simplify = T)[1,], lty=2, col=2)
      lines(seq(0,5,0.01), sapply(seq(0,5,0.01), muPq_chlor_satapornvanit09, simplify = T)[2,], lty=3, col=2)
      lines(seq(0,5,0.01), sapply(seq(0,5,0.01), muPq_chlor_satapornvanit09, simplify = T)[3,], lty=3, col=2)
    
  muPq_chlor_satapornvanit09_uncertainty<-function(In){
    if(In == 0) mup = 0 else{
      init = predict(sap.mupq, data.frame(conc=In, chem = 'chlorpyrifos'), se.fit = T)
      mup = rnorm(1, init[1], init[2]) - muPq_chlor_satapornvanit09(0)[1] #normalize to 0
    }
    while(mup < 0 || mup > 1){
      mup = rnorm(1, init[1], init[2]) - muPq_chlor_satapornvanit09(0)[1] #normalize to 0
    }
    return(mup)
  }
      
  points(seq(0,5,0.01), sapply(seq(0,5,0.01), muPq_chlor_satapornvanit09_uncertainty, simplify = T), 
          pch=5, col=4, cex = 0.5)    
      
keep.chlor.sat09 = c('muPq_chlor_satapornvanit09_uncertainty', 'muPq_chlor_satapornvanit09','sap.mupq')  
#Dimethoate ################        
  muPq_dim_satapornvanit09<-function(In){
    predict(sap.mupq, data.frame(conc=In, chem = 'dimethoate'), 
            interval = 'confidence', level = 0.95)
  }  
  
  plot(sap.mort$conc[sap.mort$chem == 'dimethoate'], sap.mort$mort[sap.mort$chem == 'dimethoate']/100, ylim = c(0,1),
       pch = 16, xlab = 'dimethoate concentration (ppb)', ylab = 'prop dead', 
       main = expression(paste('Dimethoate daily toxicity to post-larval ', italic('M. rosenbergii'), sep = '')))
    lines(seq(0,1300,10), sapply(seq(0,1300,10), muPq_dim_satapornvanit09, simplify = T)[1,], lty=2, col=2)
    lines(seq(0,1300,10), sapply(seq(0,1300,10), muPq_dim_satapornvanit09, simplify = T)[2,], lty=3, col=2)
    lines(seq(0,1300,10), sapply(seq(0,1300,10), muPq_dim_satapornvanit09, simplify = T)[3,], lty=3, col=2)
      
  muPq_dim_satapornvanit09_uncertainty<-function(In){
    if(In == 0) mup = 0 else{
      init = predict(sap.mupq, data.frame(conc=In, chem = 'dimethoate'), se.fit = T)
      mup = rnorm(1, init[1], init[2]) - muPq_dim_satapornvanit09(0)[1] #normalize to 0
    }
    while(mup < 0 || mup > 1){
      mup = rnorm(1, init[1], init[2]) - muPq_dim_satapornvanit09(0)[1] #normalize to 0
    }
    return(mup)
  }
    
  points(seq(0,1300,10), sapply(seq(0,1300,10), muPq_dim_satapornvanit09_uncertainty, simplify = T), 
          pch=5, col=4, cex = 0.5)       
keep.dim.sat09 = c('muPq_dim_satapornvanit09_uncertainty', 'muPq_dim_satapornvanit09','sap.mupq')
#Profenofos ################      
  muPq_pr_satapornvanit09<-function(In){
    predict(sap.mupq, data.frame(conc=In, chem = 'profenofos'), interval = 'confidence', level = 0.95)
  }  
      
  plot(sap.mort$conc[sap.mort$chem == 'profenofos'], sap.mort$mort[sap.mort$chem == 'profenofos']/100, ylim = c(0,1),
       pch = 16, xlab = 'profenofos concentration (ppb)', ylab = 'prop dead', 
       main = expression(paste('Profenofos daily toxicity to post-larval ', italic('M. rosenbergii'), sep = '')))      
    lines(seq(0,60,0.1), sapply(seq(0,60,0.1), muPq_pr_satapornvanit09, simplify = T)[1,], lty=2, col=2)
    lines(seq(0,60,0.1), sapply(seq(0,60,0.1), muPq_pr_satapornvanit09, simplify = T)[2,], lty=3, col=2)
    lines(seq(0,60,0.1), sapply(seq(0,60,0.1), muPq_pr_satapornvanit09, simplify = T)[3,], lty=3, col=2)
 
  muPq_prof_satapornvanit09_uncertainty<-function(In){
    if(In == 0) mup = 0 else{
      init = predict(sap.mupq, data.frame(conc=In, chem = 'profenofos'), se.fit = T)
      mup = rnorm(1, init[1], init[2]) - muPq_pr_satapornvanit09(0)[1] #normalize to 0
    }
    while(mup < 0 || mup > 1){
      mup = rnorm(1, init[1], init[2]) - muPq_pr_satapornvanit09(0)[1] #normalize to 0
    }
    return(mup)
  }
    
  points(seq(0,60,0.5), sapply(seq(0,60,0.5), muPq_prof_satapornvanit09_uncertainty, simplify = T), 
          pch=5, col=4, cex = 0.5) 
  
keep.prof.sat09 = c('muPq_prof_satapornvanit09_uncertainty', 'muPq_pr_satapornvanit09','sap.mupq')
  
#Carbendazim ###########       
  muPq_carb_satapornvanit09<-function(In){ #Paper found no consistent effect of carbendazim on mortality, even
      0*In                                 #at levels higher than the solubility of carbendazim in water
  }  
  
keep.carb.sat09 = c('muPq_carb_satapornvanit09')
#Reduced feeding rate from zinc and Chlorpyrifos from Satapornvanit et al chemosphere paper ########
sap.fr<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_feed_rate.csv")
  fr.z = subset(sap.fr, chem == 'zinc')
  fr.ch = subset(sap.fr, chem == 'chlorpyrifos')
  
#Zinc ###################
  zinc.fr= drm(feed_rate ~ conc, data = fr.z, type = 'continuous',
               fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                          fixed = c(NA, 0, max(fr.z$feed_rate), NA)))
    summary(zinc.fr)
    
  psi_q_zinc_satapornvanit09<-function(In){
    predict(zinc.fr, data.frame(conc = In), interval = 'confidence', level = 0.95) / fr.z$feed_rate[1]
  }  
    
  plot(fr.z$conc, fr.z$feed_rate / fr.z$feed_rate[1], pch = 16, ylim = c(0,1),
       xlab = 'Zinc concentration (ppb)', ylab = 'Feeding rate (prey/prawn/hr)', 
       main = 'Prawn reduced feeding (Satapornvanit-09)')
    lines(seq(0,900,10), sapply(seq(0,900,10),psi_q_zinc_satapornvanit09, simplify = T)[1,],
          lty = 2, col = 2)
    lines(seq(0,900,10), sapply(seq(0,900,10),psi_q_zinc_satapornvanit09, simplify = T)[2,],
          lty = 3, col = 2)
    lines(seq(0,900,10), sapply(seq(0,900,10),psi_q_zinc_satapornvanit09, simplify = T)[3,],
          lty = 3, col = 2)  
    
  par.tricksz = c(coef(zinc.fr), 'Upper Limit:(Intercept)' = max(fr.z$feed_rate))[c(1,3,2)]
    
  psi_q_zinc_satapornvanit09_uncertainty<-function(In){
    rdrm(nosim = 1, fct = LL.3(), mpar = par.tricksz, yerror = 'rnorm', xerror = In,
         ypar = c(0, predict(zinc.fr, data.frame(dose = In), se.fit = T)[2]))$y / fr.z$feed_rate[1]
  }
    
    points(seq(0,900,10), sapply(seq(0,900,10), psi_q_zinc_satapornvanit09_uncertainty, simplify = T), 
           pch=5, col=4, cex = 0.5)   

keep.zinc.sat09 = c(keep.zinc.sat09, 'psi_q_zinc_satapornvanit09_uncertainty', 
                    'par.tricksz', 'zinc.fr', 'fr.z')    
#Chlorpyrifos ###################
  chlor.fr= drm(feed_rate ~ conc, data = fr.ch, type = 'continuous',
               fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                          fixed = c(NA, 0, max(fr.ch$feed_rate), NA)))
    summary(chlor.fr)

  psi_q_chlor_satapornvanit09<-function(In){
    predict(chlor.fr, data.frame(conc = In), interval = 'confidence', level = 0.95) / fr.ch$feed_rate[1]
  }  
  
  plot(fr.ch$conc, fr.ch$feed_rate / fr.ch$feed_rate[1], pch = 16, ylim = c(0,1),
       xlab = 'Chlorpyrifos concentration (ppb)', ylab = 'Feeding rate (prey/prawn/hr)', 
       main = 'Prawn reduced feeding (Sata-09)')
    lines(seq(0,5,0.01), sapply(seq(0,5,0.01), psi_q_chlor_satapornvanit09, simplify = T)[1,],
          lty = 2, col = 2)
    lines(seq(0,5,0.01), sapply(seq(0,5,0.01), psi_q_chlor_satapornvanit09, simplify = T)[2,],
          lty = 3, col = 2)
    lines(seq(0,5,0.01), sapply(seq(0,5,0.01), psi_q_chlor_satapornvanit09, simplify = T)[3,],
          lty = 3, col = 2)
    
  par.tricksc = c(coef(chlor.fr), 'Upper Limit:(Intercept)' = max(fr.ch$feed_rate))[c(1,3,2)]
    
    psi_q_chlor_satapornvanit09_uncertainty<-function(In){
      rdrm(nosim = 1, fct = LL.3(), mpar = par.tricksc, yerror = 'rnorm', xerror = In,
           ypar = c(0, predict(chlor.fr, data.frame(dose = In), se.fit = T)[2]))$y / fr.ch$feed_rate[1]
    }
    
    points(seq(0,5,0.05), sapply(seq(0,5,0.05), psi_q_chlor_satapornvanit09_uncertainty, simplify = T), 
           pch=5, col=4, cex = 0.5)   

keep.chlor.sat09 = c(keep.chlor.sat09, 'psi_q_chlor_satapornvanit09_uncertainty', 
                    'par.tricksc', 'chlor.fr', 'fr.ch')    
#*Note: checked that this produced the desired relationship on predator feeding rate and it does