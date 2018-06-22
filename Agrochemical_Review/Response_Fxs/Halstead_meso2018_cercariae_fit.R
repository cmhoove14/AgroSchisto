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

L.3.fx = function(t, lc50 = lc50, slp = slp){
  1 / (1+exp(slp*log(t / lc50)))
}

halstead_cerc <- read_csv("Agrochemical_Review/Response_Fxs/Data/Halstead_meso_cercariae.csv") %>% 
  group_by(Treatment, Time) %>% 
  summarise(ntot = sum(Total),
            nfx = sum(Dead)) %>% 
  mutate(per_mort = nfx/ntot)
  
halstead_cerc_ts <- drm(nfx/ntot ~ Time, weights = ntot, type = 'binomial', fct = LL2.2(), curveid = Treatment,
                        data = halstead_cerc)


#Control parameters
  lc50_halstead18_ctrl <- summary(halstead_cerc_ts)$coef[which(grepl( "S", rownames(summary(halstead_cerc_ts)$coef))),][2,1]
  lc50_halstead18_ctrl_se <- summary(halstead_cerc_ts)$coef[which(grepl( "S", rownames(summary(halstead_cerc_ts)$coef))),][2,2]
  b_halstead18_ctrl <- summary(halstead_cerc_ts)$coef[which(grepl( "S", rownames(summary(halstead_cerc_ts)$coef))),][1,1]
  b_halstead18_ctrl_se <- summary(halstead_cerc_ts)$coef[which(grepl( "S", rownames(summary(halstead_cerc_ts)$coef))),][1,2]
  
piC_ctrl_halstead18_ts_uncertainty <- function(...){
  e0 = rnorm(1, lc50_halstead18_ctrl, lc50_halstead18_ctrl_se)
  b0 = rnorm(1, b_halstead18_ctrl, b_halstead18_ctrl_se)
    
    auc0 = integrate(L.3.fx, lc50 = e0, slp = b0, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value
    
  return(auc0)  
}

#Atrazine parameters
  lc50_halstead18_atr <- summary(halstead_cerc_ts)$coef[which(grepl( "A", rownames(summary(halstead_cerc_ts)$coef))),][2,1]
  lc50_halstead18_atr_se <- summary(halstead_cerc_ts)$coef[which(grepl( "A", rownames(summary(halstead_cerc_ts)$coef))),][2,2]
  b_halstead18_atr <- summary(halstead_cerc_ts)$coef[which(grepl( "A", rownames(summary(halstead_cerc_ts)$coef))),][1,1]
  b_halstead18_atr_se <- summary(halstead_cerc_ts)$coef[which(grepl( "A", rownames(summary(halstead_cerc_ts)$coef))),][1,2]
  
piC_atr102_halstead18_ts_uncertainty <- function(...){
  e = rnorm(1, lc50_halstead18_atr, lc50_halstead18_atr_se)
  b = rnorm(1, b_halstead18_atr, b_halstead18_atr_se)
    
    auc = integrate(L.3.fx, lc50 = e, slp = b, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value
    
  return(auc / piC_ctrl_halstead18_ts_uncertainty())  
}

#Chlorpyrifos parameters
  lc50_halstead18_chlor <- summary(halstead_cerc_ts)$coef[which(grepl( "C", rownames(summary(halstead_cerc_ts)$coef))),][2,1]
  lc50_halstead18_chlor_se <- summary(halstead_cerc_ts)$coef[which(grepl( "C", rownames(summary(halstead_cerc_ts)$coef))),][2,2]
  b_halstead18_chlor <- summary(halstead_cerc_ts)$coef[which(grepl( "C", rownames(summary(halstead_cerc_ts)$coef))),][1,1]
  b_halstead18_chlor_se <- summary(halstead_cerc_ts)$coef[which(grepl( "C", rownames(summary(halstead_cerc_ts)$coef))),][1,2]
  
piC_chlor64_halstead18_ts_uncertainty <- function(...){
  e = rnorm(1, lc50_halstead18_chlor, lc50_halstead18_chlor_se)
  b = rnorm(1, b_halstead18_chlor, b_halstead18_chlor_se)
    
    auc = integrate(L.3.fx, lc50 = e, slp = b, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value
    
  return(auc / piC_ctrl_halstead18_ts_uncertainty())  
}

#Fertilizer parameters
  lc50_halstead18_fer <- summary(halstead_cerc_ts)$coef[which(grepl( "F", rownames(summary(halstead_cerc_ts)$coef))),][2,1]
  lc50_halstead18_fer_se <- summary(halstead_cerc_ts)$coef[which(grepl( "F", rownames(summary(halstead_cerc_ts)$coef))),][2,2]
  b_halstead18_fer <- summary(halstead_cerc_ts)$coef[which(grepl( "F", rownames(summary(halstead_cerc_ts)$coef))),][1,1]
  b_halstead18_fer_se <- summary(halstead_cerc_ts)$coef[which(grepl( "F", rownames(summary(halstead_cerc_ts)$coef))),][1,2]
  
piC_fer4400_halstead18_ts_uncertainty <- function(...){
  e = rnorm(1, lc50_halstead18_fer, lc50_halstead18_fer_se)
  b = rnorm(1, b_halstead18_fer, b_halstead18_fer_se)
    
    auc = integrate(L.3.fx, lc50 = e, slp = b, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value
    
  return(auc / piC_ctrl_halstead18_ts_uncertainty())  
}
