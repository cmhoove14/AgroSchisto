#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Data extraction and model fitting to Hussein 2011 [Bakry et al 2016](https://www.researchgate.net/publication/305149067_Effects_of_Three_Inorganic_Fertilizers_on_the_Biology_and_Histopathology_of_infected_Biomphalaria_alexandrina_snails?enrichId=rgreq-cf1a38509b7c5460cbc30f698ea58594-XXX&enrichSource=Y292ZXJQYWdlOzMwNTE0OTA2NztBUzozODI2MjU5NTgxMjE0NzJAMTQ2ODIzNjU0NTE4NQ%3D%3D&el=1_x_3&_esc=publicationCoverPdf) data
source("Agrochemical_Review/Models/litchfield_wilcoxon_get_b1_from_slope.R")

#balanced fertilizer reported LC50 and slope data #####################
  lc50.hus.balanced = 505.7
  slp.hus.balanced = 1.89
  b1.hus.balanced = get_b1(slp.hus.balanced)
  #get standard error from reported 95% CIs of lc50
      se.lc50.hus.balanced = mean(c(log10(625.33/lc50.hus.balanced), log10(lc50.hus.balanced/390.31))) / 1.96
  
#Create function based on reverse of litchfield and wilcoxon      
muNq_balanced_hussein16_uncertainty = function(Fe){
  fer = (Fe/1000) #Parameters based on ppm, data input as ppb
  lc50 = 10^(rnorm(1, log10(lc50.hus.balanced), se.lc50.hus.balanced)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
    mun = pnorm(b1.hus.balanced * log10(fer/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }

#high phosphorous fertilizer reported LC50 and slope data #####################
  lc50.hus.highp = 1600
  slp.hus.highp = 1.45
  b1.hus.highp = get_b1(slp.hus.highp)
  #get standard error from reported 95% CIs of lc50
      se.lc50.hus.highp = mean(c(log10(1843.4/lc50.hus.highp), log10(lc50.hus.highp/1356.56))) / 1.96
  
#Create function based on reverse of litchfield and wilcoxon      
muNq_highp_hussein16_uncertainty = function(Fe){
  fer = (Fe/1000) #Parameters based on ppm, data input as ppb
  lc50 = 10^(rnorm(1, log10(lc50.hus.highp), se.lc50.hus.highp)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
    mun = pnorm(b1.hus.highp * log10(fer/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }

#balanced fertilizer reported LC50 and slope data #####################
  lc50.hus.highn = 9500
  slp.hus.highn = 1.16
  b1.hus.highn = get_b1(slp.hus.highn)
  #get standard error from reported 95% CIs of lc50
      se.lc50.hus.highn = mean(c(log10(10108.5/lc50.hus.highn), log10(lc50.hus.highn/8891.4))) / 1.96
  
#Create function based on reverse of litchfield and wilcoxon      
muNq_highn_hussein16_uncertainty = function(Fe){
  fer = (Fe/1000) #Parameters based on ppm, data input as ppb
  lc50 = 10^(rnorm(1, log10(lc50.hus.highn), se.lc50.hus.highn)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
    mun = pnorm(b1.hus.highn * log10(fer/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }


#Keep vector
keep.hussein16 <- c("lc50.hus.balanced", "b1.hus.balanced", "se.lc50.hus.balanced", "muNq_balanced_hussein16_uncertainty",
                    "lc50.hus.highp", "b1.hus.highp", "se.lc50.hus.highp", "muNq_highp_hussein16_uncertainty",
                    "lc50.hus.highn", "b1.hus.highn", "se.lc50.hus.highn", "muNq_highn_hussein16_uncertainty")