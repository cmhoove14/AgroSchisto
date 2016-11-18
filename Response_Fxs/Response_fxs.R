require(ggplot2)
require(drc)

#Atrazine effect on phi_Nq (snail carrying capacity) ################
atra.df<-data.frame('atra' = c(0,1,10,30,100),              #Raw atrazine concentration (ppb)
                    'logatra' = log(c(0,1,10,30,100)+1),    #Log atrazine concentration (ppb)
                    'phiNq' = c(0,0.2888,0.6535,0,1.3215))  #Snail population response measured as peak snail growth rate and 
                                                            #interpreted as changes in snail carrying capacity

plot(atra.df$logatra, atra.df$phiNq, pch = 16)

atra_mod<-glm(phiNq ~ logatra+0, data=atra.df)

atra.slope = atra_mod$coefficients[1]

f_phi_Nq_at = function(He){ #Function to use in model
  phi_Nq = parameters['phi_N'] + parameters['phi_N'] * (atra.slope*log(He+1))
  return(phi_Nq)
}

#24-hr mortality functions for predator (macrobrachium) populations from Satapornvanit et al chemosphere paper ########
  sap.mort<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_mort.csv")
    
    z.mod.mupq<-drm(sap.mort$mort[sap.mort$chem == 'zinc']/100 ~ sap.mort$conc[sap.mort$chem == 'zinc'],
                    data = sap.mort, type = 'binomial', fct = LL.2())
    
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
#Insecticide info from Halstead et al chemosphere paper ################
#calculate decay rates (k) from hydrolysis half lives for each chemical in Halstead 2015 (from table S1 and pmep.cce.cornell.edu) #############
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
  
#Insecticide toxicity to crustaceans from Halstead et al ################
  data<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/Cray.LC50.2012.csv')
  #Malathion
  mal<-subset(data, Chem == 'Mal')
    mal.c<-unique(mal$Conc)
    mal.d<-as.numeric()
    for(i in 1:length(mal.c)){
      mal.d[i] = sum(mal$Dead[mal$Conc == mal.c[i]])/5
    }
      mal.mod<-drm(mal.d ~ mal.c, type = 'binomial', fct = LL2.2())
      
      f_muPq_mal_Halstead<-function(In){
        muplus = predict(mal.mod, newdata = data.frame(conc = In), type = 'response')/4
        muuse = parameters['mu_P'] + muplus
        if(muuse > 1){
          muuse = 1
        }
        return(muuse)
      }  

  #Chlorpyrifos   
  chlor<-subset(data, Chem == 'Chlor')
    chlor.c<-unique(chlor$Conc)
    chlor.d<-as.numeric()
    for(i in 1:length(chlor.c)){
      chlor.d[i] = sum(chlor$Dead[chlor$Conc == chlor.c[i]])/5
    }
      chlor.mod<-drm(chlor.d ~ chlor.c, type = 'binomial', fct = LL2.2())
  
      f_muPq_ch_Halstead<-function(In){
        muplus = predict(chlor.mod, newdata = data.frame(conc = In), type = 'response')/4
        muuse = parameters['mu_P'] + muplus
        if(muuse > 1){
          muuse = 1
        }
        return(muuse)
      }
      
  #Terbufos   
  terb<-subset(data, Chem == 'Terb')
    terb.c<-unique(terb$Conc)
    terb.d<-as.numeric()
    for(i in 1:length(terb.c)){
      terb.d[i] = sum(terb$Dead[terb$Conc == terb.c[i]])/5
    }
      terb.mod<-drm(terb.d ~ terb.c, type = 'binomial', fct = LL2.2())
  
      f_muPq_terb_Halstead<-function(In){
        muplus = predict(terb.mod, newdata = data.frame(conc = In), type = 'response')/4
        muuse = parameters['mu_P'] + muplus
        if(muuse > 1){
          muuse = 1
        }
        return(muuse)
      }
      
  
  #Lambda-cyhalothrin   
  lamcy<-subset(data, Chem == 'Lambda')
    lamcy.c<-unique(lamcy$Conc)
    lamcy.d<-as.numeric()
    for(i in 1:length(lamcy.c)){
      lamcy.d[i] = sum(lamcy$Dead[lamcy$Conc == lamcy.c[i]])/5
    }
      lamcy.mod<-drm(lamcy.d ~ lamcy.c, type = 'binomial', fct = LL2.2())
  
      f_muPq_lamcy_Halstead<-function(In){
        muplus = predict(lamcy.mod, newdata = data.frame(conc = In), type = 'response')/4
        muuse = parameters['mu_P'] + muplus
        if(muuse > 1){
          muuse = 1
        }
        return(muuse)
      }
      
  #esfenvalerate  
  esfen<-subset(data, Chem == 'Esfen')
    esfen.c<-unique(esfen$Conc)
    esfen.d<-as.numeric()
    for(i in 1:length(esfen.c)){
      esfen.d[i] = sum(esfen$Dead[esfen$Conc == esfen.c[i]])/5
    }
      esfen.mod<-drm(esfen.d ~ esfen.c, type = 'binomial', fct = LL2.2())
  
      f_muPq_esfen_Halstead<-function(In){
        muplus = predict(esfen.mod, newdata = data.frame(conc = In), type = 'response')/4
        muuse = parameters['mu_P'] + muplus
        if(muuse > 1){
          muuse = 1
        }
        return(muuse)
      }
      
  #Permethrin
  perm<-subset(data, Chem == 'Perm')
    perm.c<-unique(perm$Conc)
    perm.d<-as.numeric()
    for(i in 1:length(perm.c)){
      perm.d[i] = sum(perm$Dead[perm$Conc == perm.c[i]])/5
    }
      perm.mod<-drm(perm.d ~ perm.c, type = 'binomial', fct = LL2.2())
      
      f_muPq_perm_Halstead<-function(In){
        muplus = predict(perm.mod, newdata = data.frame(conc = In), type = 'response')/4
        muuse = parameters['mu_P'] + muplus
        if(muuse > 1){
          muuse = 1
        }
        return(muuse)
      }
      