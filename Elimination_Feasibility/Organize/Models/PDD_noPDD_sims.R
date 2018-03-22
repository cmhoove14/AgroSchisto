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

#Load parameter sets and eqbm values from Get_eqbm_for_trans_parameter_sets.R script
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/model_fit_profile_likelihood_parameters.Rdata")
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/best_fit_params.Rdata")
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/trans_pars_eqbms_k_0.08noPdd_profileFit.Rdata")
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/trans_pars_eqbms_k_0.08withPdd_profileFit.Rdata")
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/trans_pars_eqbms_k_0.08noPdd.addNDD_profileFit.Rdata")
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/trans_pars_eqbms_k_0.08withPdd.addNDD_profileFit.Rdata")

shortlist_pars <- fin_pars95ci %>%
  arrange(lamda_twa) %>% 
  head(100) #keep 100 lowest parameter sets

Prevalence <- function(W, k) {
  p=1 - (1/(1+W/k)^(k))*(1+W/(1+W/k)) # fraction of humans with at least 2 parasites
  return(p)
}

# next create 100 models with these beta, lamda combinations and run with DD, get 100 BBrate curves
covrg = 0.8  # 80% coverage
eff = 0.94 # 94% efficacy from the Egg reduction rate reported for S. haematobium in <http://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0003286>
k.fit = 0.08

mda.years = c(1:20) # annual MDA for 20 years
time = c(1:(365*61)) #1 year run up, 20 years annual mda followed by 40 years rebound with no MDA

#1 annual round of MDA
mda.events = data.frame(var = rep('Wt', length(mda.years)),
                         time = c(mda.years*365),
                         value = rep((1 - eff), length(mda.years)),
                         method = rep('mult', length(mda.years)))

#Objects to fill in loop over 100 parameter sets
#array with all 100 60-year runs with and w/o PDD
  det.runs.k008.mda1 = array(data = NA, dim = c(length(time), ncol(eqbm)+1, nrow(shortlist_pars), 4))
  
#Vector of W-pre and W-pos values for non-NDD added models
  w.pre.k008.mda1 = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  w.pos.k008.mda1 = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  w.pre.k008.mda1.pdd = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  w.pos.k008.mda1.pdd = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))

#Vector of W-pre and W-pos values for NDD added models
  w.pre.k008.mda1.addNDD = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  w.pos.k008.mda1.addNDD = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  w.pre.k008.mda1.pdd.addNDD = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  w.pos.k008.mda1.pdd.addNDD = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  
#Vector of Prev-pre and Prev-pos values for non-NDD added models
  prev.pre.k008.mda1 = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  prev.pos.k008.mda1 = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  prev.pre.k008.mda1.pdd = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  prev.pos.k008.mda1.pdd = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))

#Vector of Prev-pre and Prev-pos values for NDD added models
  prev.pre.k008.mda1.addNDD = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  prev.pos.k008.mda1.addNDD = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  prev.pre.k008.mda1.pdd.addNDD = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  prev.pos.k008.mda1.pdd.addNDD = array(data = NA, dim = c(nrow(shortlist_pars), length(mda.years)))
  
params<-pars_Chris1
  params["cov"]<-covrg
  params["k"]<-k.fit

for(y in 1:nrow(shortlist_pars)){
    params["lamda"]<-shortlist_pars$lamda_twa[y] 

#No PDD, no added NDDs    
  nstart = c(S=eqbm$S[y], 
             E=eqbm$E[y], 
             I=eqbm$I[y], 
             Wt=eqbm$Wt[y], 
             Wu=eqbm$Wu[y])
  
  det.runs.k008.mda1[, , y, 1] = ode(nstart, time, schisto_mod_nopdd, params,
                           events = list(data = mda.events))
    for(i in mda.years){
      w.pre.k008.mda1[y,i] = covrg*det.runs.k008.mda1[(i*365), 5, y, 1] + (1-covrg) * det.runs.k008.mda1[(i*365), 6, y, 1]
      w.pos.k008.mda1[y,i] = covrg*det.runs.k008.mda1[(i*365+1), 5, y, 1] + (1-covrg) * det.runs.k008.mda1[(i*365+1), 6, y, 1]
      prev.pre.k008.mda1[y,i] = Prevalence(w.pre.k008.mda1[y,i], k.fit)
      prev.pos.k008.mda1[y,i] = Prevalence(w.pos.k008.mda1[y,i], k.fit)
    } 
  
#No PDD add NDDs  
  nstart.addNDD = c(S=eqbm.addNDD$S[y], 
                    E=eqbm.addNDD$E[y], 
                    I=eqbm.addNDD$I[y], 
                    Wt=eqbm.addNDD$Wt[y], 
                    Wu=eqbm.addNDD$Wu[y])
  
  det.runs.k008.mda1[, , y, 2] = ode(nstart.addNDD, time, schisto_mod_nopdd_add_ndds, params,
                           events = list(data = mda.events))
    for(i in mda.years){
      w.pre.k008.mda1.addNDD[y,i] = covrg*det.runs.k008.mda1[(i*365), 5, y, 2] + (1-covrg) * det.runs.k008.mda1[(i*365), 6, y, 2]
      w.pos.k008.mda1.addNDD[y,i] = covrg*det.runs.k008.mda1[(i*365+1), 5, y, 2] + (1-covrg) * det.runs.k008.mda1[(i*365+1), 6, y, 2]
      prev.pre.k008.mda1.addNDD[y,i] = Prevalence(w.pre.k008.mda1.addNDD[y,i], k.fit)
      prev.pos.k008.mda1.addNDD[y,i] = Prevalence(w.pos.k008.mda1.addNDD[y,i], k.fit)
    } 

#With PDD, no added NDDs    
  nstart.pdd = c(S=eqbm.pdd$S[y], 
                 E=eqbm.pdd$E[y], 
                 I=eqbm.pdd$I[y], 
                 Wt=eqbm.pdd$Wt[y], 
                 Wu=eqbm.pdd$Wu[y])
  
  det.runs.k008.mda1[, , y, 3] = ode(nstart.pdd, time, schisto_mod_pdd, params,
                           events = list(data = mda.events))
  
  for(i in mda.years){
    w.pre.k008.mda1.pdd[y,i] = covrg*det.runs.k008.mda1[(i*365), 5, y, 3] + (1-covrg) * det.runs.k008.mda1[(i*365), 6, y, 3]
    w.pos.k008.mda1.pdd[y,i] = covrg*det.runs.k008.mda1[(i*365+1), 5, y, 3] + (1-covrg) * det.runs.k008.mda1[(i*365+1), 6, y, 3]
    prev.pre.k008.mda1.pdd[y,i] = Prevalence(w.pre.k008.mda1.pdd[y,i], k.fit)
    prev.pos.k008.mda1.pdd[y,i] = Prevalence(w.pre.k008.mda1.pdd[y,i], k.fit)
  } 
  
#With PDD And added NDDs    
  nstart.pdd.addNDD = c(S=eqbm.pdd.addNDD$S[y], 
                        E=eqbm.pdd.addNDD$E[y], 
                        I=eqbm.pdd.addNDD$I[y], 
                        Wt=eqbm.pdd.addNDD$Wt[y], 
                        Wu=eqbm.pdd.addNDD$Wu[y])
  
  det.runs.k008.mda1[, , y, 4] = ode(nstart.pdd.addNDD, time, schisto_mod_pdd_add_ndds, params,
                           events = list(data = mda.events))
  
  for(i in mda.years){
    w.pre.k008.mda1.pdd.addNDD[y,i] = covrg*det.runs.k008.mda1[(i*365), 5, y, 4] + (1-covrg) * det.runs.k008.mda1[(i*365), 6, y, 4]
    w.pos.k008.mda1.pdd.addNDD[y,i] = covrg*det.runs.k008.mda1[(i*365+1), 5, y, 4] + (1-covrg) * det.runs.k008.mda1[(i*365+1), 6, y, 4]
    prev.pre.k008.mda1.pdd.addNDD[y,i] = Prevalence(w.pos.k008.mda1.pdd.addNDD[y,i], k.fit)
    prev.pos.k008.mda1.pdd.addNDD[y,i] = Prevalence(w.pos.k008.mda1.pdd.addNDD[y,i], k.fit)
  } 
  
  if(y %in% seq(1, 101, 9)){
    plot(det.runs.k008.mda1[, 1, y, 4], det.runs.k008.mda1[, 5, y, 4], type = "l", lwd = 2)
      lines(det.runs.k008.mda1[, 1, y, 2], det.runs.k008.mda1[, 5, y, 2], type = "l", lwd = 2, col = 2)
  }

  print(y)
}
  
#Some post processing #########
  mean.w.k008.mda1 = matrix(ncol = 3, nrow = length(time))
    mean.w.k008.mda1[,1] = rowMeans(det.runs.k008.mda1[ , 5, , 1])
    mean.w.k008.mda1[,2] = rowMeans(det.runs.k008.mda1[ , 6, , 1])
    mean.w.k008.mda1[,3] = covrg * mean.w.k008.mda1[,1] + (1-covrg) * mean.w.k008.mda1[,2]
   
  mean.w.k008.mda1.addNDD = matrix(ncol = 3, nrow = length(time))
    mean.w.k008.mda1.addNDD[,1] = rowMeans(det.runs.k008.mda1[ , 5, , 2])
    mean.w.k008.mda1.addNDD[,2] = rowMeans(det.runs.k008.mda1[ , 6, , 2])
    mean.w.k008.mda1.addNDD[,3] = covrg * mean.w.k008.mda1.addNDD[,1] + (1-covrg) * mean.w.k008.mda1.addNDD[,2]
     
  mean.w.pdd.k008.mda1 = matrix(ncol = 3, nrow = length(time))
    mean.w.pdd.k008.mda1[,1] = rowMeans(det.runs.k008.mda1[ , 5, , 3])
    mean.w.pdd.k008.mda1[,2] = rowMeans(det.runs.k008.mda1[ , 6, , 3])
    mean.w.pdd.k008.mda1[,3] = covrg * mean.w.pdd.k008.mda1[,1] + (1-covrg) * mean.w.pdd.k008.mda1[,2]
   
  mean.w.pdd.k008.mda1.addNDD = matrix(ncol = 3, nrow = length(time))
    mean.w.pdd.k008.mda1.addNDD[,1] = rowMeans(det.runs.k008.mda1[ , 5, , 4])
    mean.w.pdd.k008.mda1.addNDD[,2] = rowMeans(det.runs.k008.mda1[ , 6, , 4])
    mean.w.pdd.k008.mda1.addNDD[,3] = covrg * mean.w.pdd.k008.mda1.addNDD[,1] + (1-covrg) * mean.w.pdd.k008.mda1.addNDD[,2]
  
#Save outputs ###########
#Mean w matrices    
  save(mean.w.k008.mda1, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/mean_w_mda_noPdd_noAddNDD_k008_profileFit.Rdata")  
  save(mean.w.k008.mda1.addNDD, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/mean_w_mda_noPdd_AddNDD_k008_profileFit.Rdata")  
  save(mean.w.pdd.k008.mda1, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/mean_w_mda_Pdd_noAddNDD_k008_profileFit.Rdata")  
  save(mean.w.pdd.k008.mda1.addNDD, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/mean_w_mda_Pdd_AddNDD_k008_profileFit.Rdata")  

# Pre and most MDA W estimates
  save(w.pre.k008.mda1, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/w_pre_mda_noPdd_noAddNDD_k008_profileFit.Rdata")  
  save(w.pos.k008.mda1, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/w_post_mda_noPdd_noAddNDD_k008_profileFit.Rdata")  
  
  save(w.pre.k008.mda1.pdd, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/w_pre_mda_Pdd_noAddNDD_k008_profileFit.Rdata") 
  save(w.pos.k008.mda1.pdd, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/w_pos_mda_Pdd_noAddNDD_k008_profileFit.Rdata") 

  save(w.pre.k008.mda1.addNDD, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/w_pre_mda_noPdd_AddNDD_k008_profileFit.Rdata")  
  save(w.pos.k008.mda1.addNDD, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/w_pos_mda_noPdd_AddNDD_k008_profileFit.Rdata")  
       
  save(w.pre.k008.mda1.pdd.addNDD, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/w_pre_mda_Pdd_AddNDD_k008_profileFit.Rdata")  
  save(w.pos.k008.mda1.pdd.addNDD, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/w_pos_mda_Pdd_AddNDD_k008_profileFit.Rdata")  
  
# Pre and most MDA Prevalence estimates
  save(prev.pre.k008.mda1, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/prev_pre_mda_noPdd_noAddNDD_k008_profileFit.Rdata")  
  save(prev.pos.k008.mda1, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/prev_post_mda_noPdd_noAddNDD_k008_profileFit.Rdata")  
  
  save(prev.pre.k008.mda1.pdd, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/prev_pre_mda_Pdd_noAddNDD_k008_profileFit.Rdata") 
  save(prev.pos.k008.mda1.pdd, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/prev_pos_mda_Pdd_noAddNDD_k008_profileFit.Rdata") 

  save(prev.pre.k008.mda1.addNDD, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/prev_pre_mda_noPdd_AddNDD_k008_profileFit.Rdata")  
  save(prev.pos.k008.mda1.addNDD, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/prev_pos_mda_noPdd_AddNDD_k008_profileFit.Rdata")  
       
  save(prev.pre.k008.mda1.pdd.addNDD, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/prev_pre_mda_Pdd_AddNDD_k008_profileFit.Rdata")  
  save(prev.pos.k008.mda1.pdd.addNDD, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/prev_pos_mda_Pdd_AddNDD_k008_profileFit.Rdata")  

#Full simulation array
  save(det.runs.k008.mda1, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/k008_mda_fullarray_profileFit.Rdata")