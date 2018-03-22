#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

##################################################################################################

require(deSolve)
require(graphics)
require(ggplot2)

source("Elimination_Feasibility/Organize/Models/Reff_BBR_fns.R")
source("Elimination_Feasibility/Organize/Models/schisto_mods_pdd_nopdd.R")

#Get transmission parameter sets
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/model_fit_profile_likelihood_parameters.Rdata")
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/best_fit_params.Rdata")

hist(fin_pars95ci$lamda_twa, breaks = 20)
fin_pars95ci %>% 
  ggplot(aes(x = lamda_twa, y = negLL)) + geom_point(pch = 16) + theme_bw()

#Reduce parameter sets to explore based on stricter fit
shortlist_pars <- fin_pars95ci %>%
  arrange(lamda_twa) %>% 
  head(100) #keep 100 lowest parameter sets
  
hist(shortlist_pars$lamda_twa, breaks = 20)

#eqb, dfs to fill  
eqbm = data.frame(S = 0, E = 0, I = 0,
                       Wt = 0, Wu = 0)

eqbm.pdd = data.frame(S = 0, E = 0, I = 0,
                      Wt = 0, Wu = 0)

eqbm.addNDD = data.frame(S = 0, E = 0, I = 0,
                         Wt = 0, Wu = 0)

eqbm.pdd.addNDD = data.frame(S = 0, E = 0, I = 0,
                             Wt = 0, Wu = 0)

#Fill shortlist parameter vector with equilibrium estimates
  k.fit = 0.08 #clumping parameter derived from fit to epi data
  covrg = 0.8 #80% coverage assumed
  
  #time and starting conditions for each model
    time<-seq(from=0, to=200*365, by=10)
    nstart=c(S=3892,E=3750,I=1200, Wt=50, Wu=50)
    params<-pars_Chris1
    params['cov']<-covrg
    params['k'] = k.fit  

for(i in 1:nrow(shortlist_pars)){
  #Change parameter values to reflect parameter set
    params["lamda"]<-shortlist_pars$lamda_twa[i] 

  #Run each model to equilibrium  
      output<-as.data.frame(ode(nstart,time,schisto_mod_nopdd,params))
      output.addNDD<-as.data.frame(ode(nstart,time,schisto_mod_nopdd_add_ndds,params))
      
      output.PDD<-as.data.frame(ode(nstart,time,schisto_mod_pdd,params))
      output.PDD.addNDD<-as.data.frame(ode(nstart,time,schisto_mod_pdd_add_ndds,params))

  #Store equilibrium values    
    eqbm[i,] = output[dim(output)[1],c(2:6)]
    eqbm.addNDD[i,] = output.addNDD[dim(output.addNDD)[1],c(2:6)]
    eqbm.pdd[i,] = output.PDD[dim(output.PDD)[1],c(2:6)]
    eqbm.pdd.addNDD[i,] = output.PDD.addNDD[dim(output.PDD.addNDD)[1],c(2:6)]
    
    print(i)

  }

  save(eqbm, file = paste0("Elimination_Feasibility/Organize/Models/Outputs_Refs/", 
                           "trans_pars_eqbms_k_", 
                           as.character(k.fit),"noPdd_profileFit.Rdata"))
  
  save(eqbm.addNDD, file = paste0("Elimination_Feasibility/Organize/Models/Outputs_Refs/", 
                           "trans_pars_eqbms_k_", 
                           as.character(k.fit),"noPdd.addNDD_profileFit.Rdata"))
   
  save(eqbm.pdd, file = paste0("Elimination_Feasibility/Organize/Models/Outputs_Refs/", 
                               "trans_pars_eqbms_k_", 
                               as.character(k.fit),"withPdd_profileFit.Rdata"))
  
  save(eqbm.pdd.addNDD, file = paste0("Elimination_Feasibility/Organize/Models/Outputs_Refs/", 
                                      "trans_pars_eqbms_k_", 
                                      as.character(k.fit),"withPdd.addNDD_profileFit.Rdata"))