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
  data<-read_csv('Agrochemical_Review/Response_Fxs/Data/rohr_macrobrachium_tox_unpublished.csv') %>% 
    mutate(dead_24hr = ifelse(Time <= 24 & Time > 0, 1, 0),
#Change dose to dose of active ingredient (treats acetone concentration as 0)
    Dose = ifelse(Treatment == "Cont_Ace", 0, Dose)) %>%  
    filter(Replicate == 2)
  
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
mal_mac <- data %>%
    dplyr::filter(Treatment %in% c("Cont_Ace", "Cont_H20", "Mal")) %>%
    group_by(Dose) %>%
    summarise(nrx = sum(dead_24hr),
              ntot = n())
  
#DRC model  
  mal.mod.mac<-drm(nrx / ntot ~ Dose, weights = ntot, data = mal_mac, type = 'binomial',  
                   fct = LL2.2())
  
    rohr_unpub_mal_mac_lc50 <- summary(mal.mod.mac)$coef[2,1]   # Lc50 estimate from model
    rohr_unpub_mal_mac_lc50_se <- summary(mal.mod.mac)$coef[2,2]  #SE of lc50 estimate from model
    rohr_unpub_mal_mac_b <- summary(mal.mod.mac)$coef[1,1]  #slope parameter of log-logistic fitted here
    
    muPq_mal_mac_rohr_unpub_uncertainty<-function(In){
        
        e = rnorm(1, rohr_unpub_mal_mac_lc50, rohr_unpub_mal_mac_lc50_se)
      
        mort <- 1 / (1 + exp(rohr_unpub_mal_mac_b*(log(In)-e)))
        
        mort
      }
  
  keep.rohr.unpublished.mal = c('mal_mac', 'rohr_unpub_mal_mac_lc50', 'rohr_unpub_mal_mac_lc50_se', 
                                'rohr_unpub_mal_mac_b', 'muPq_mal_mac_rohr_unpub_uncertainty')
  
#Chlorpyrifos   *************************************************************************############
chlor_mac <- data %>%
    dplyr::filter(Treatment %in% c("Cont_Ace", "Cont_H20", "Chlor")) %>%
    group_by(Dose) %>%
    summarise(nrx = sum(dead_24hr),
              ntot = n())
  
#DRC model  
  chlor.mod.mac<-drm(nrx / ntot ~ Dose, weights = ntot, data = chlor_mac, type = 'binomial',  
                   fct = LL2.2())
  
    rohr_unpub_chlor_mac_lc50 <- summary(chlor.mod.mac)$coef[2,1]   # Lc50 estimate from model
    rohr_unpub_chlor_mac_lc50_se <- summary(chlor.mod.mac)$coef[2,2]  #SE of lc50 estimate from model
    rohr_unpub_chlor_mac_b <- summary(chlor.mod.mac)$coef[1,1]  #slope parameter of log-logistic fitted here
    
    muPq_chlor_mac_rohr_unpub_uncertainty<-function(In){
        
        e = rnorm(1, rohr_unpub_chlor_mac_lc50, rohr_unpub_chlor_mac_lc50_se)
      
        mort <- 1 / (1 + exp(rohr_unpub_chlor_mac_b*(log(In)-e)))
        
        mort
    }
    
  keep.rohr.unpublished.chlor = c('chlor_mac', 'rohr_unpub_chlor_mac_lc50', 'rohr_unpub_chlor_mac_lc50_se', 
                                'rohr_unpub_chlor_mac_b', 'muPq_chlor_mac_rohr_unpub_uncertainty')

#Terbufos   *****************************************************************************########
terb_mac <- data %>%
    dplyr::filter(Treatment %in% c("Cont_Ace", "Cont_H20", "Terb")) %>%
    group_by(Dose) %>%
    summarise(nrx = sum(dead_24hr),
              ntot = n())
  
#DRC model  
  terb.mod.mac<-drm(nrx / ntot ~ Dose, weights = ntot, data = terb_mac, type = 'binomial',  
                   fct = LL2.2())
  
    rohr_unpub_terb_mac_lc50 <- summary(terb.mod.mac)$coef[2,1]   # Lc50 estimate from model
    rohr_unpub_terb_mac_lc50_se <- summary(terb.mod.mac)$coef[2,2]  #SE of lc50 estimate from model
    rohr_unpub_terb_mac_b <- summary(terb.mod.mac)$coef[1,1]  #slope parameter of log-logistic fitted here
    
    muPq_terb_mac_rohr_unpub_uncertainty<-function(In){
        
        e = rnorm(1, rohr_unpub_terb_mac_lc50, rohr_unpub_terb_mac_lc50_se)
      
        mort <- 1 / (1 + exp(rohr_unpub_terb_mac_b*(log(In)-e)))
        
        mort
    }
    
  keep.rohr.unpublished.terb = c('terb_mac', 'rohr_unpub_terb_mac_lc50', 'rohr_unpub_terb_mac_lc50_se', 
                                'rohr_unpub_terb_mac_b', 'muPq_terb_mac_rohr_unpub_uncertainty')

#Lambda-cyhalothrin   *******************************************************************##########
lamcy_mac <- data %>%
    dplyr::filter(Treatment %in% c("Cont_Ace", "Cont_H20", "Lamda")) %>%
    group_by(Dose) %>%
    summarise(nrx = sum(dead_24hr),
              ntot = n())
  
#DRC model  
  lamcy.mod.mac<-drm(nrx / ntot ~ Dose, weights = ntot, data = lamcy_mac, type = 'binomial',  
                   fct = LL2.2())
  
    rohr_unpub_lamcy_mac_lc50 <- summary(lamcy.mod.mac)$coef[2,1]   # Lc50 estimate from model
    rohr_unpub_lamcy_mac_lc50_se <- summary(lamcy.mod.mac)$coef[2,2]  #SE of lc50 estimate from model
    rohr_unpub_lamcy_mac_b <- summary(lamcy.mod.mac)$coef[1,1]  #slope parameter of log-logistic fitted here
    
    muPq_lamcy_mac_rohr_unpub_uncertainty<-function(In){
        
        e = rnorm(1, rohr_unpub_lamcy_mac_lc50, rohr_unpub_lamcy_mac_lc50_se)
      
        mort <- 1 / (1 + exp(rohr_unpub_lamcy_mac_b*(log(In)-e)))
        
        mort
    }
    
  keep.rohr.unpublished.lamcy = c('lamcy_mac', 'rohr_unpub_lamcy_mac_lc50', 'rohr_unpub_lamcy_mac_lc50_se', 
                                'rohr_unpub_lamcy_mac_b', 'muPq_lamcy_mac_rohr_unpub_uncertainty')
  
#esfenvalerate  *************************************************************************#######
esfen_mac <- data %>%
    dplyr::filter(Treatment %in% c("Cont_Ace", "Cont_H20", "Esfen")) %>%
    group_by(Dose) %>%
    summarise(nrx = sum(dead_24hr),
              ntot = n())
  
#DRC model  
  esfen.mod.mac<-drm(nrx / ntot ~ Dose, weights = ntot, data = esfen_mac, type = 'binomial',  
                   fct = LL2.2())
  
    rohr_unpub_esfen_mac_lc50 <- summary(esfen.mod.mac)$coef[2,1]   # Lc50 estimate from model
    rohr_unpub_esfen_mac_lc50_se <- summary(esfen.mod.mac)$coef[2,2]  #SE of lc50 estimate from model
    rohr_unpub_esfen_mac_b <- summary(esfen.mod.mac)$coef[1,1]  #slope parameter of log-logistic fitted here
    
    muPq_esfen_mac_rohr_unpub_uncertainty<-function(In){
        
        e = rnorm(1, rohr_unpub_esfen_mac_lc50, rohr_unpub_esfen_mac_lc50_se)
      
        mort <- 1 / (1 + exp(rohr_unpub_esfen_mac_b*(log(In)-e)))
        
        mort
    }
    
  keep.rohr.unpublished.esfen = c('esfen_mac', 'rohr_unpub_esfen_mac_lc50', 'rohr_unpub_esfen_mac_lc50_se', 
                                'rohr_unpub_esfen_mac_b', 'muPq_esfen_mac_rohr_unpub_uncertainty')
      
#Permethrin *****************************************************************************##########
perm_mac <- data %>%
    dplyr::filter(Treatment %in% c("Cont_Ace", "Cont_H20", "Perm")) %>%
    group_by(Dose) %>%
    summarise(nrx = sum(dead_24hr),
              ntot = n())
  
#DRC model  
  perm.mod.mac<-drm(nrx / ntot ~ Dose, weights = ntot, data = perm_mac, type = 'binomial',  
                   fct = LL2.2())
  
    rohr_unpub_perm_mac_lc50 <- summary(perm.mod.mac)$coef[2,1]   # Lc50 estimate from model
    rohr_unpub_perm_mac_lc50_se <- summary(perm.mod.mac)$coef[2,2]  #SE of lc50 estimate from model
    rohr_unpub_perm_mac_b <- summary(perm.mod.mac)$coef[1,1]  #slope parameter of log-logistic fitted here
    
    muPq_perm_mac_rohr_unpub_uncertainty<-function(In){
        
        e = rnorm(1, rohr_unpub_perm_mac_lc50, rohr_unpub_perm_mac_lc50_se)
      
        mort <- 1 / (1 + exp(rohr_unpub_perm_mac_b*(log(In)-e)))
        
        mort
    }
    
  keep.rohr.unpublished.perm = c('perm_mac', 'rohr_unpub_perm_mac_lc50', 'rohr_unpub_perm_mac_lc50_se', 
                                'rohr_unpub_perm_mac_b', 'muPq_perm_mac_rohr_unpub_uncertainty')
