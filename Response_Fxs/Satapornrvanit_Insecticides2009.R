require(ggplot2)
require(drc)

#24-hr mortality functions for predator (macrobrachium) populations from Satapornvanit et al chemosphere paper ########
  sap.mort<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_mort.csv")
    
  #Zinc
  z.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'zinc']/100 ~ sap.mort$conc[sap.mort$chem == 'zinc'],
                  data = sap.mort, type = 'binomial', fct = LL.2())
    
    muPq_zinc_satapornvanit09<-function(In){
      
      1/(1+exp(z.mod.mupq$coefficients[1]*(log(In)-log(z.mod.mupq$coefficients[2]))))
      
    }  
  #Chlorpyrifos  
  ch.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'chlorpyrifos']/100 ~ sap.mort$conc[sap.mort$chem == 'chlorpyrifos'],
                    data = sap.mort, type = 'binomial', fct = LL.2())
    
    muPq_chlor_satapornvanit09<-function(In){
      
      1/(1+exp(ch.mod.mupq$coefficients[1]*(log(In)-log(ch.mod.mupq$coefficients[2]))))
      
    }  
  #Dimethoate
  dim.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'dimethoate']/100 ~ sap.mort$conc[sap.mort$chem == 'dimethoate'],
                    data = sap.mort, type = 'binomial', fct = LL.2())
    
    muPq_dim_satapornvanit09<-function(In){
      
      1/(1+exp(dim.mod.mupq$coefficients[1]*(log(In)-log(dim.mod.mupq$coefficients[2]))))
      
    }  
  #profenofos
  pr.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'profenofos']/100 ~ sap.mort$conc[sap.mort$chem == 'profenofos'],
                   data = sap.mort, type = 'binomial', fct = LL.2())
  
    muPq_pr_satapornvanit09<-function(In){
      
      1/(1+exp(pr.mod.mupq$coefficients[1]*(log(In)-log(pr.mod.mupq$coefficients[2]))))
      
    }  
    
    
  #NOTE: carbendazim not modeled because it has no effect on toxicity, but will be included as such in simulations
    muPq_carb_satapornvanit09<-function(In){
      In*0
    }
#Reduced feeding rate from zinc and Chlorpyrifos from Satapornvanit et al chemosphere paper ########
sap.fr<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_feed_rate.csv")
  fr.z = subset(sap.fr, chem == 'zinc' & feed_rate > 0 )
  fr.ch = subset(sap.fr, chem == 'chlorpyrifos' & feed_rate > 0)

  #Get slope only and fix y intercept at 1 (i.e. no reduction in feeding rate)
    z.mod.fr<-lm(log(per_control) ~ conc + 0, data = fr.z) 
    
      pred_q_z<-function(In){
        exp(z.mod.fr$coefficients[1]*In)
      }
      #*Note: checked that this produced the desired relationship on predator feeding rate and it does
    
    ch.mod.fr<-lm(log(per_control) ~ conc + 0, data = fr.ch) 

    pred_q_ch<-function(In){
      exp(ch.mod.fr$coefficients[1]*In)
    }
    #*Note: checked that this produced the desired relationship on predator feeding rate and it does
#Combine predator mortality and feeding rate reduction from Satapornvanit into single function ##############
  Satapornvanit2009.z<-function(In){
    return(c(f_muPq_z(In), pred_q_z(In)))
  }
    
  Satapornvanit2009.ch<-function(In){
    return(c(f_muPq_ch(In), pred_q_ch(In)))
  }
  
  Satapornvanit2009.dim<-function(In){
    return(c(f_muPq_dim(In), 1)) #No effect on pred feeding rate
  }
  
  Satapornvanit2009.pr<-function(In){
    return(c(f_muPq_pr(In), 1)) #No effect on pred feeding rate
  }
  
  Satapornvanit2009.carb<-function(In){
    return(c(0, 1)) #No effect on predator mortality or feeding rate
  }