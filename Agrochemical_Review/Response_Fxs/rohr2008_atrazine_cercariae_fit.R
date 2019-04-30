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


atr = read.csv('Agrochemical_Review/Response_Fxs/Data/rohr08_atr.csv')

atr.lin = lm(surv ~ log_conc, weights = st_err^-1, data = atr)

rohr.atr.fx = function(He){
  heu = log(He+1)
  predict(atr.lin, newdata = data.frame(log_conc = heu), interval = 'confidence', level = 0.95)
}

#Final function estimates relative mortality to He=0 and is interpreted as pi_C parameter      
  piC.atr.rohr08.lin = function(He){
      ts = predict(atr.lin, newdata = data.frame(log_conc = log(He+1)), se.fit = TRUE)[1:2]
      piC = rnorm(1, ts$fit, ts$se.fit) / atr.lin$coefficients[1]
    
    return(piC)
  } 
  
#Time series data at single concnetrations  
rohr_cerc <- read_csv("Agrochemical_Review/Response_Fxs/Data/rohr2008.csv") %>% 
  group_by(chem, time_hrs) %>% 
  summarise(ntot = sum(total),
            nfx = sum(round(dead)),
            conc = mean(conc)) %>% 
  mutate(per_mort = nfx/ntot)

rohr_cerc_ts <- drm(nfx/ntot ~ time_hrs, weights = ntot, type = 'binomial', fct = LL2.2(), curveid = chem,
                    data = rohr_cerc)
#Control parameters
  lc50_rohr08_ctrl <- summary(rohr_cerc_ts)$coef[which(grepl( "control", rownames(summary(rohr_cerc_ts)$coef))),][2,1]
  lc50_rohr08_ctrl_se <- summary(rohr_cerc_ts)$coef[which(grepl( "control", rownames(summary(rohr_cerc_ts)$coef))),][2,2]
  b_rohr08_ctrl <- summary(rohr_cerc_ts)$coef[which(grepl( "control", rownames(summary(rohr_cerc_ts)$coef))),][1,1]
  b_rohr08_ctrl_se <- summary(rohr_cerc_ts)$coef[which(grepl( "control", rownames(summary(rohr_cerc_ts)$coef))),][1,2]
  
piC_ctrl_rohr08_ts_uncertainty <- function(...){
  e0 = rnorm(1, lc50_rohr08_ctrl, lc50_rohr08_ctrl_se)
  b0 = rnorm(1, b_rohr08_ctrl, b_rohr08_ctrl_se)
    
    auc0 = integrate(L.3.fx, lc50 = e0, slp = b0, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value
    
  return(auc0)  
}
  
#Atrazine parameters and function 
  lc50_rohr08_atr201 <- summary(rohr_cerc_ts)$coef[which(grepl( "atrazine", rownames(summary(rohr_cerc_ts)$coef))),][2,1]
  lc50_rohr08_atr201_se <- summary(rohr_cerc_ts)$coef[which(grepl( "atrazine", rownames(summary(rohr_cerc_ts)$coef))),][2,2]
  b_rohr08_atr201 <- summary(rohr_cerc_ts)$coef[which(grepl( "atrazine", rownames(summary(rohr_cerc_ts)$coef))),][1,1]
  b_rohr08_atr201_se <- summary(rohr_cerc_ts)$coef[which(grepl( "atrazine", rownames(summary(rohr_cerc_ts)$coef))),][1,2]
  
piC_atr201_rohr08_ts_uncertainty <- function(...){
  e = rnorm(1, lc50_rohr08_atr201, lc50_rohr08_atr201_se)
  b = rnorm(1, b_rohr08_atr201, b_rohr08_atr201_se)
    
    auc = integrate(L.3.fx, lc50 = e, slp = b, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value
    
  return(auc / piC_ctrl_rohr08_ts_uncertainty())  
}

#carbaryl parameters and function 
  lc50_rohr08_carb33.5 <- summary(rohr_cerc_ts)$coef[which(grepl( "carbaryl", rownames(summary(rohr_cerc_ts)$coef))),][2,1]
  lc50_rohr08_carb33.5_se <- summary(rohr_cerc_ts)$coef[which(grepl( "carbaryl", rownames(summary(rohr_cerc_ts)$coef))),][2,2]
  b_rohr08_carb33.5 <- summary(rohr_cerc_ts)$coef[which(grepl( "carbaryl", rownames(summary(rohr_cerc_ts)$coef))),][1,1]
  b_rohr08_carb33.5_se <- summary(rohr_cerc_ts)$coef[which(grepl( "carbaryl", rownames(summary(rohr_cerc_ts)$coef))),][1,2]
  
piC_carb33.5_rohr08_ts_uncertainty <- function(...){
        
  e = rnorm(1, lc50_rohr08_carb33.5, lc50_rohr08_carb33.5_se)
  b = rnorm(1, b_rohr08_carb33.5, b_rohr08_carb33.5_se)
    
    auc = integrate(L.3.fx, lc50 = e, slp = b, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value
    
  return(auc / piC_ctrl_rohr08_ts_uncertainty())  
}

#glyphosate parameters and function 
  lc50_rohr08_gly3700 <- summary(rohr_cerc_ts)$coef[which(grepl( "glyphosate", rownames(summary(rohr_cerc_ts)$coef))),][2,1]
  lc50_rohr08_gly3700_se <- summary(rohr_cerc_ts)$coef[which(grepl( "glyphosate", rownames(summary(rohr_cerc_ts)$coef))),][2,2]
  b_rohr08_gly3700 <- summary(rohr_cerc_ts)$coef[which(grepl( "glyphosate", rownames(summary(rohr_cerc_ts)$coef))),][1,1]
  b_rohr08_gly3700_se <- summary(rohr_cerc_ts)$coef[which(grepl( "glyphosate", rownames(summary(rohr_cerc_ts)$coef))),][1,2]
  
piC_gly3700_rohr08_ts_uncertainty <- function(...){
        
  e = rnorm(1, lc50_rohr08_gly3700, lc50_rohr08_gly3700_se)
  b = rnorm(1, b_rohr08_gly3700, b_rohr08_gly3700_se)
    
    auc = integrate(L.3.fx, lc50 = e, slp = b, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value
    
  return(auc / piC_ctrl_rohr08_ts_uncertainty())  
}

#malathion parameters and function 
  lc50_rohr08_mal9.6 <- summary(rohr_cerc_ts)$coef[which(grepl( "malathion", rownames(summary(rohr_cerc_ts)$coef))),][2,1]
  lc50_rohr08_mal9.6_se <- summary(rohr_cerc_ts)$coef[which(grepl( "malathion", rownames(summary(rohr_cerc_ts)$coef))),][2,2]
  b_rohr08_mal9.6 <- summary(rohr_cerc_ts)$coef[which(grepl( "malathion", rownames(summary(rohr_cerc_ts)$coef))),][1,1]
  b_rohr08_mal9.6_se <- summary(rohr_cerc_ts)$coef[which(grepl( "malathion", rownames(summary(rohr_cerc_ts)$coef))),][1,2]

piC_mal9.6_rohr08_ts_uncertainty <- function(...){
        
  e = rnorm(1, lc50_rohr08_mal9.6, lc50_rohr08_mal9.6_se)
  b = rnorm(1, b_rohr08_mal9.6, b_rohr08_mal9.6_se)
    
    auc = integrate(L.3.fx, lc50 = e, slp = b, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value
    
  return(auc / piC_ctrl_rohr08_ts_uncertainty())  
}

#Effects on snail reproduction at tested concentrations from same study
rohr08_ctrl_eggs_per_day <- 27.83
rohr08_ctrl_eggs_per_day_se <- 2.60

rohr08_atr_eggs_per_day <- 23.79
rohr08_atr_eggs_per_day_se <- 3.90

fNq_atr201_rohr08_uncertainty <- function(...){
  rnorm(1, rohr08_atr_eggs_per_day, rohr08_atr_eggs_per_day_se) / 
    rnorm(1, rohr08_ctrl_eggs_per_day, rohr08_ctrl_eggs_per_day_se)
}

rohr08_gly_eggs_per_day <- 29.93
rohr08_gly_eggs_per_day_se <- 3.90

fNq_gly3700_rohr08_uncertainty <- function(...){
  rnorm(1, rohr08_gly_eggs_per_day, rohr08_gly_eggs_per_day_se) / 
    rnorm(1, rohr08_ctrl_eggs_per_day, rohr08_ctrl_eggs_per_day_se)
}

rohr08_carb_eggs_per_day <- 28.53
rohr08_carb_eggs_per_day_se <- 3.79

fNq_carb33.5_rohr08_uncertainty <- function(...){
  rnorm(1, rohr08_carb_eggs_per_day, rohr08_carb_eggs_per_day_se) / 
    rnorm(1, rohr08_ctrl_eggs_per_day, rohr08_ctrl_eggs_per_day_se)
}

rohr08_mal_eggs_per_day <- 28.40
rohr08_mal_eggs_per_day_se <- 3.49

fNq_mal9.6_rohr08_uncertainty <- function(...){
  rnorm(1, rohr08_mal_eggs_per_day, rohr08_mal_eggs_per_day_se) / 
    rnorm(1, rohr08_ctrl_eggs_per_day, rohr08_ctrl_eggs_per_day_se)
}