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
require(dplyr)

#Insecticide toxicity to crustaceans from Halstead et al; 4-day mortality endpoints ################
  data<-read.csv('Agrochemical_Review/Response_Fxs/Data/Halstead_Cray_LC50_2012.csv')
  
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
mal <- data %>%
    filter(Chem == "Mal") %>%
    group_by(Conc) %>%
    summarise(nrx = sum(Dead)) %>%
    mutate(ntot = 5)
  
#DRC model  
  mal.mod<-drm(nrx / ntot ~ Conc, weights = ntot, data = mal, type = 'binomial',  
               fct = LL2.2())
  
    halstead_mal_lc50 <- summary(mal.mod)$coef[2,1]   # Lc50 estimate from model
    halstead_mal_lc50_se <- summary(mal.mod)$coef[2,2]  #SE of lc50 estimate from model
    halstead_mal_b <- summary(mal.mod)$coef[1,1]  #slope parameter of log-logistic fitted here
    
    muPq_mal_Halstead_uncertainty<-function(In){
        
        b = halstead_mal_b
        lc50 = halstead_mal_lc50
        lc50_se = halstead_mal_lc50_se

          e = rnorm(1, lc50, lc50_se)
      
        mort <- 1 / (1 + exp(b*(log(In)-e)))
        
        mort
      }
      
  keep.hal15.mal = c('mal', 'halstead_mal_lc50', 'halstead_mal_lc50_se', 'halstead_mal_b', 'muPq_mal_Halstead_uncertainty')
  
#Chlorpyrifos   *************************************************************************############
chlor <- data %>%
    filter(Chem == "Chlor") %>%
    group_by(Conc) %>%
    summarise(nrx = sum(Dead)) %>%
    mutate(ntot = 5)
      
#DRC model  
  chlor.mod<-drm(nrx / ntot ~ Conc, weights = ntot, data = chlor, type = 'binomial',  
               fct = LL2.2())
      
      halstead_chlor_lc50 <- summary(chlor.mod)$coef[2,1]   # Lc50 estimate from model
      halstead_chlor_lc50_se <- summary(chlor.mod)$coef[2,2]   # SE of LC50 estimate from model
      halstead_chlor_b <- summary(chlor.mod)$coef[1,1]  #slope parameter of log-logistic fitted here
      
      muPq_chlor_Halstead_uncertainty<-function(In){ #Same structure as general predict_LL2.2 function, but specific to chlorpyrifos parameters
        
        b = halstead_chlor_b
        lc50 = halstead_chlor_lc50
        lc50_se = halstead_chlor_lc50_se

          e = rnorm(1, lc50, lc50_se)
        
        mort <- 1 / (1 + exp(b*(log(In)-e)))
        
        mort
      }
  
  keep.hal15.chlor = c('chlor', 'halstead_chlor_lc50', 'halstead_chlor_lc50_se', 'halstead_chlor_b', 'muPq_chlor_Halstead_uncertainty')

#Terbufos   *****************************************************************************########
terb <- data %>%
    filter(Chem == "Terb") %>%
    group_by(Conc) %>%
    summarise(nrx = sum(Dead)) %>%
    mutate(ntot = 5)
      
#DRC model  
  terb.mod<-drm(nrx / ntot ~ Conc, weights = ntot, data = terb, type = 'binomial',  
                 fct = LL2.2())
      
      halstead_terb_lc50 <- summary(terb.mod)$coef[2,1]   #Lc50 estimate from model
      halstead_terb_lc50_se <- summary(terb.mod)$coef[2,2]   #SE of lc50 estimate from model
      halstead_terb_b <- summary(terb.mod)$coef[1,1]  #slope parameter of log-logistic fitted here
      
      muPq_terb_Halstead_uncertainty<-function(In){ #Same structure as general predict_LL2.2 function, but specific to terbufos parameters
        
        b = halstead_terb_b
        lc50 = halstead_terb_lc50
        lc50_se = halstead_terb_lc50_se
        unc = 1
        
        if(unc == 1){
          e = rnorm(1, lc50, lc50_se)
        }
        else{
          e = lc50
        }
        
        mort <- 1 / (1 + exp(b*(log(In)-e)))
        
        mort
      }

  keep.hal15.terb = c('terb', 'halstead_terb_lc50', 'halstead_terb_lc50_se', 'halstead_terb_b', 'muPq_terb_Halstead_uncertainty')
      

#Lambda-cyhalothrin   *******************************************************************##########
lamcy <- data %>%
    filter(Chem == "Lambda") %>%
    group_by(Conc) %>%
    summarise(nrx = sum(Dead)) %>%
    mutate(ntot = 5,
           Conc_ppt = Conc*1000)
      
#DRC model  
  lamcy.mod<-drm(nrx / ntot ~ Conc_ppt, weights = ntot, data = lamcy, type = 'binomial',  
                 fct = LL2.2())
      
      halstead_lamcy_lc50 <- summary(lamcy.mod)$coef[2,1]   # lc50 estimate from model
      halstead_lamcy_lc50_se <- summary(lamcy.mod)$coef[2,2]   #Reported 95%Ci of lc50 in halstead et al paper converted to SE
      halstead_lamcy_b <- summary(lamcy.mod)$coef[1,1]  #slope parameter of log-logistic fitted here
      
      muPq_lamcy_Halstead_uncertainty<-function(In){ #Same structure as general predict_LL2.2 function, but specific to lamcyufos parameters
        Ins = In*1000   #Transform from ppb to ppt
        
        b = halstead_lamcy_b
        lc50 = log(halstead_lamcy_lc50)
        lc50_se = halstead_lamcy_lc50_se

          e = rnorm(1, lc50, lc50_se)

        mort <- 1 / (1 + exp(b*(log(Ins)-e)))
        
        mort
      }

  keep.hal15.lamcy = c('lamcy', 'halstead_lamcy_lc50', 'halstead_lamcy_lc50_se', 'halstead_lamcy_b', 'muPq_lamcy_Halstead_uncertainty')
  
#esfenvalerate  *************************************************************************#######
esfen <- data %>%
    filter(Chem == "Esfen") %>%
    group_by(Conc) %>%
    summarise(nrx = sum(Dead)) %>%
    mutate(ntot = 5,
           Conc_ppt = Conc*1000)
      
#DRC model  
  esfen.mod<-drm(nrx / ntot ~ Conc_ppt, weights = ntot, data = esfen, type = 'binomial',  
                 fct = LL2.2())
      
      halstead_esfen_lc50 <- summary(esfen.mod)$coef[2,1]   #Lc50 estimate from model
      halstead_esfen_lc50_se <- summary(esfen.mod)$coef[2,2]   #SE of lc50 estimate from model
      halstead_esfen_b <- summary(esfen.mod)$coef[1,1]  #slope parameter of log-logistic fitted here
      
      muPq_esfen_Halstead_uncertainty<-function(In){ #Same structure as general predict_LL2.2 function, but specific to esfenvalerate parameters
        Ins = In*1000 #Transform from ppb to ppt
        
        b = halstead_esfen_b
        lc50 = halstead_esfen_lc50
        lc50_se = halstead_esfen_lc50_se

          e = rnorm(1, lc50, lc50_se)

        mort <- 1 / (1 + exp(b*(log(Ins)-e)))
        
        mort
      }

  keep.hal15.esfen = c('esfen', 'halstead_esfen_lc50', 'halstead_esfen_lc50_se', 'halstead_esfen_b', 'muPq_esfen_Halstead_uncertainty')
      
#Permethrin *****************************************************************************##########
perm <- data %>%
    filter(Chem == "Perm") %>%
    group_by(Conc) %>%
    summarise(nrx = sum(Dead)) %>%
    mutate(ntot = 5,
           Conc_ppt = Conc*1000)
      
#DRC model  
  perm.mod<-drm(nrx / ntot ~ Conc_ppt, weights = ntot, data = perm, type = 'binomial',  
                 fct = LL2.2())
      
      halstead_perm_lc50 <- summary(perm.mod)$coef[2,1]   #Lc50 estimate from model
      halstead_perm_lc50_se <- summary(perm.mod)$coef[2,2]   #SE of lc50 estimate from model
      halstead_perm_b <- summary(perm.mod)$coef[1,1]  #slope parameter of log-logistic fitted here
      
      muPq_perm_Halstead_uncertainty<-function(In){ #Same structure as general predict_LL2.2 function, but specific to permethrin parameters
        Ins = In*1000 #Transform from ppb to ppt
        
        b = halstead_perm_b
        lc50 = halstead_perm_lc50
        lc50_se = halstead_perm_lc50_se
        unc = 1
        
        if(unc == 1){
          e = rnorm(1, lc50, lc50_se)
        }
        else{
          e = lc50
        }
        
        mort <- 1 / (1 + exp(b*(log(Ins)-e)))
        
        mort
      }

  keep.hal15.perm = c('perm', 'halstead_perm_lc50', 'halstead_perm_lc50_se', 'halstead_perm_b', 'muPq_perm_Halstead_uncertainty')

#Character vector of (all) objects to keep from this script ###############
  keep.hal15.all = c(keep.hal15.mal, keep.hal15.chlor, keep.hal15.terb, keep.hal15.esfen, keep.hal15.lamcy, keep.hal15.perm)