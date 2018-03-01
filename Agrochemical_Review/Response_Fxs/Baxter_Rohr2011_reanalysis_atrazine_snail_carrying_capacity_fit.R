#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Atrazine effect on phi_Nq (snail carrying capacity) from Baxter et al 2011 data with Rohr et al analysis####
atra.df<-data.frame('atra' = c(0,1,10,30,100),              #Raw atrazine concentration (ppb)
                    'logatra' = log(c(0,1,10,30,100)+1),    #Log+1 atrazine concentration (ppb)
                    'growthrate' = c(0.119406, 0.153891, 0.19744, 0.118918, 0.27719), #Peak growth rate in Rohr reanalysis of Baxter data interpreted as changes in snail carrying capacity 
                    'st.err' = c(0.015949, 0.026329, 0.035467, 0.016286, 0.047113))    

bax.mod = lm(growthrate ~ logatra, weights = st.err^-1, data = atra.df)
  base.growth = predict(bax.mod, newdata = data.frame(logatra = 0), 
                  type = 'response')
  
#Function to estimate change in carrying capacity with uncertainty  
phi_Nq_atr_baxrohr = function(He){
    u = predict(bax.mod, newdata = data.frame(logatra = log(He+1)), type = 'response',
                se.fit = TRUE)[1:2]
    phiNq = rnorm(1, u$fit, u$se.fit) / base.growth   #Get proportional increase in peak growth as increase in carrying capacity

  return(phiNq)
  
}

#Exclude the 30ppb data point as it had high block variability in original study
atra.df.no30 = subset(atra.df, atra != 30)

bax.mod.no30 = lm(growthrate ~ logatra, weights = st.err^-1, data = atra.df.no30)
  base.growth.no30 = predict(bax.mod, newdata = data.frame(logatra = 0), 
                        type = 'response')
  
#Function to estimate change in carrying capacity with uncertainty  
phi_Nq_atr_baxrohr.no30 = function(He){
    u = predict(bax.mod.no30, newdata = data.frame(logatra = log(He+1)), type = 'response',
                se.fit = TRUE)[1:2]
    phiNq = rnorm(1, u$fit, u$se.fit) / base.growth.no30
  
  return(phiNq)
  
}


#Function to use this info for any herbicide given its eec, assuming a proportional response of carrying capacity to atrazine conc
  eec.atr = 102
  phi_Nq_rel_baxrohr = function(He, eec){
      prop = He / eec
      rel = prop*eec.atr
      
      u = predict(bax.mod, newdata = data.frame(logatra = log(rel+1)), type = 'response',
                  se.fit = TRUE)[1:2]
      phiNq = rnorm(1, u$fit, u$se.fit) / base.growth
    
    return(phiNq)
    
  }
  
 
  phi_Nq_rel_baxrohr.no30 = function(He, eec){
      prop = He / eec    #get proportion of input concentration to cehm's EEC
      rel = prop*eec.atr #get atrazine concentration equivalent  
      
      u = predict(bax.mod.no30, newdata = data.frame(logatra = log(rel+1)), type = 'response',
                  se.fit = TRUE)[1:2]
      phiNq = rnorm(1, u$fit, u$se.fit) / base.growth.no30

    return(phiNq)
    
  }
  
  keep.baxrohr = c('atra.df', 'phi_Nq_atr_baxrohr', 'bax.mod', 'base.growth', 'phi_Nq_atr_baxrohr.no30', 'bax.mod.no30', 'base.growth.no30', 'phi_Nq_rel_baxrohr', 'eec.atr')