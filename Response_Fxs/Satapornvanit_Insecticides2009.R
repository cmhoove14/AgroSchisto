require(ggplot2)
require(drc)

#24-hr mortality functions for predator (macrobrachium) populations from Satapornvanit et al chemosphere paper ########
  sap.mort<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_mort.csv")
  
#Zinc ########
  plot(sap.mort$conc[sap.mort$chem == 'zinc'], sap.mort$mort[sap.mort$chem == 'zinc']/100,
       pch = 16, xlab = 'Zinc concentration (ppb)', ylab = 'prop dead')
  
  z.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'zinc']/100 ~ sap.mort$conc[sap.mort$chem == 'zinc'],
                  data = sap.mort, type = 'binomial', fct = LL.2())
    
    muPq_zinc_satapornvanit09<-function(In){
      1/(1+exp(z.mod.mupq$coefficients[1]*(log(In)-log(z.mod.mupq$coefficients[2]))))
    }  
    
  lines(c(0:max(sap.mort$conc[sap.mort$chem == 'zinc'])), 
        muPq_zinc_satapornvanit09(c(0:max(sap.mort$conc[sap.mort$chem == 'zinc']))), 
        lty=2, col='red')  
  
#Chlorpyrifos ###########
  ch.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'chlorpyrifos']/100 ~ sap.mort$conc[sap.mort$chem == 'chlorpyrifos'],
                    data = sap.mort, type = 'binomial', fct = LL.2())
  
  plot(sap.mort$conc[sap.mort$chem == 'chlorpyrifos'], sap.mort$mort[sap.mort$chem == 'chlorpyrifos']/100,
       pch = 16, xlab = 'Chlorpyrifos concentration (ppb)', ylab = 'prop dead', ylim = c(0,1))
    
    muPq_chlor_satapornvanit09<-function(In){
      1/(1+exp(ch.mod.mupq$coefficients[1]*(log(In)-log(ch.mod.mupq$coefficients[2]))))
    }  
    
  lines(seq(0, max(sap.mort$conc[sap.mort$chem == 'chlorpyrifos']), 0.1), 
        muPq_chlor_satapornvanit09(seq(0, max(sap.mort$conc[sap.mort$chem == 'chlorpyrifos']), 0.1)), 
        lty=2, col='red')  
    
#Dimethoate #########
  dim.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'dimethoate']/100 ~ sap.mort$conc[sap.mort$chem == 'dimethoate'],
                    data = sap.mort, type = 'binomial', fct = LL.2())
  
  plot(sap.mort$conc[sap.mort$chem == 'dimethoate'], sap.mort$mort[sap.mort$chem == 'dimethoate']/100,
       pch = 16, xlab = 'dimethoate concentration (ppb)', ylab = 'prop dead', ylim = c(0,1))
    
    muPq_dim_satapornvanit09<-function(In){
      1/(1+exp(dim.mod.mupq$coefficients[1]*(log(In)-log(dim.mod.mupq$coefficients[2]))))
    }  
    
  lines(c(0:max(sap.mort$conc[sap.mort$chem == 'dimethoate'])), 
        muPq_dim_satapornvanit09(c(0:max(sap.mort$conc[sap.mort$chem == 'dimethoate']))), 
        lty=2, col='red') 
    
#profenofos #############
  pr.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'profenofos']/100 ~ sap.mort$conc[sap.mort$chem == 'profenofos'],
                   data = sap.mort, type = 'binomial', fct = LL.2())
  
  plot(sap.mort$conc[sap.mort$chem == 'profenofos'], sap.mort$mort[sap.mort$chem == 'profenofos']/100,
       pch = 16, xlab = 'profenofos concentration (ppb)', ylab = 'prop dead', ylim = c(0,1))
  
    muPq_pr_satapornvanit09<-function(In){
      1/(1+exp(pr.mod.mupq$coefficients[1]*(log(In)-log(pr.mod.mupq$coefficients[2]))))
    }  
    
  lines(c(0:max(sap.mort$conc[sap.mort$chem == 'profenofos'])), 
        muPq_pr_satapornvanit09(c(0:max(sap.mort$conc[sap.mort$chem == 'profenofos']))), 
        lty=2, col='red') 
    
#NOTE: carbendazim not modeled because it has no effect on toxicity, but will be included as such in simulations #####
    muPq_carb_satapornvanit09<-function(In){
      In*0
    }

#Reduced feeding rate from zinc and Chlorpyrifos from Satapornvanit et al chemosphere paper ########
sap.fr<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_feed_rate.csv")
  fr.z = subset(sap.fr, chem == 'zinc' & feed_rate > 0 )
  fr.ch = subset(sap.fr, chem == 'chlorpyrifos' & feed_rate > 0)
  
#Zinc
  plot(fr.z$conc, fr.z$per_control, pch = 16, ylim = c(0,1),
       xlab = 'Zinc concentration (ppb)', ylab = 'Feeding rate relative to control')
  
  z.mod.fr = nls(per_control ~ exp(-b*conc), data = fr.z, start = list(b=0.01))
    summary(z.mod.fr)
    
  psi_q_zinc_satapornvanit09<-function(In){
    exp(-summary(z.mod.fr)$parameters[1]*In)
  }  
    
  lines(c(0:501), psi_q_zinc_satapornvanit09(c(0:501)), lty=2, col='red')  
  
#Chlorpyrifos
  plot(fr.ch$conc, fr.ch$per_control, pch = 16, ylim = c(0,1),
       xlab = 'Chlorpyrifos concentration (ppb)', ylab = 'Feeding rate relative to control')
  
  ch.mod.fr = nls(per_control ~ exp(-b*conc), data = fr.ch, start = list(b=0.01))
    summary(ch.mod.fr)
  
  psi_q_chlor_satapornvanit09<-function(In){
    exp(-summary(ch.mod.fr)$parameters[1]*In)
  }  
  
  lines(seq(0, max(sap.mort$conc[sap.mort$chem == 'chlorpyrifos']), 0.1), 
        psi_q_chlor_satapornvanit09(seq(0, max(sap.mort$conc[sap.mort$chem == 'chlorpyrifos']), 0.1)), 
        lty=2, col='red') 

    
#*Note: checked that this produced the desired relationship on predator feeding rate and it does