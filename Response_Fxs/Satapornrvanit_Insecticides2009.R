require(ggplot2)
require(drc)

#24-hr mortality functions for predator (macrobrachium) populations from Satapornvanit et al chemosphere paper ########
  sap.mort<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_mort.csv")
    
    z.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'zinc']/100 ~ sap.mort$conc[sap.mort$chem == 'zinc'],
                    data = sap.mort, type = 'binomial', fct = LL.4(fixed = c(NA, NA, 1, NA)))
    
    f_muPq_z_Satapornvanit<-function(In){
      muplus = predict(z.mod.mupq, newdata = data.frame(conc = In), type = 'response')
      muuse = parameters['mu_P'] + muplus
      if(muuse > 1){
        muuse = 1
      }
      return(muuse)
    }
    
    ch.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'chlorpyrifos']/100 ~ sap.mort$conc[sap.mort$chem == 'chlorpyrifos'],
                     data = sap.mort, type = 'binomial', fct = LL.2())
    
    f_muPq_ch_Satapornvanit<-function(In){
      muplus = predict(ch.mod.mupq, newdata = data.frame(conc = In), type = 'response')
      muuse = parameters['mu_P'] + muplus
      if(muuse > 1){
        muuse = 1
      }
      return(muuse)
    }
    
    dim.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'dimethoate']/100 ~ sap.mort$conc[sap.mort$chem == 'dimethoate'],
                      data = sap.mort, type = 'binomial', fct = LL.2())
    
    f_muPq_dim_Satapornvanit<-function(In){
      muplus = predict(dim.mod.mupq, newdata = data.frame(conc = In), type = 'response')
      muuse = parameters['mu_P'] + muplus
      if(muuse > 1){
        muuse = 1
      }
      return(muuse)
    }
    
    pr.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'profenofos']/100 ~ sap.mort$conc[sap.mort$chem == 'profenofos'],
                     data = sap.mort, type = 'binomial', fct = LL.2())
    
    f_muPq_pr_Satapornvanit<-function(In){
      muplus = predict(pr.mod.mupq, newdata = data.frame(conc = In), type = 'response')
      muuse = parameters['mu_P'] + muplus
      if(muuse > 1){
        muuse = 1
      }
      return(muuse)
    }
    
    #NOTE: carbendazim not modeled because it has no effect on toxicity, but will be included as such in simulations
    f_muPq_carb_Satapornvanit<-function(In){
      muplus = 0
      muuse = parameters['mu_P'] + muplus
      return(muuse)
    }
#Reduced feeding rate from zinc and Chlorpyrifos from Satapornvanit et al chemosphere paper ########
sap.fr<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_feed_rate.csv")
  fr.z = subset(sap.fr, chem == 'zinc' & feed_rate > 0 )
  fr.ch = subset(sap.fr, chem == 'chlorpyrifos' & feed_rate > 0)

  #Get slope only and fix y intercept at 1 (i.e. no reduction in feeding rate)
    z.mod.fr<-lm(log(per_control) ~ conc + 0, data = fr.z) 
    
      pred_q_z<-function(In){
        pred_red = exp(predict(z.mod.fr, newdata = data.frame(conc = In), type = 'response'))
        return(pred_red)
      }
      #*Note: checked that this produced the desired relationship on predator feeding rate and it does
    
    ch.mod.fr<-lm(log(per_control) ~ conc + 0, data = fr.ch) 

    pred_q_ch<-function(In){
      pred_red = exp(predict(ch.mod.fr, newdata = data.frame(conc = In), type = 'response'))
      return(pred_red)
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