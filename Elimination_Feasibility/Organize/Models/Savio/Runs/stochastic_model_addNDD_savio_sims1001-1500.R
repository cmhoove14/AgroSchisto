#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Load packages and other files ##########
source("~/ElimFeas_StochMod/PostReview/Refs/stochastic_model_addNDD_savio.R")
source("~/ElimFeas_StochMod/PostReview/Refs/schisto_mods_pdd_nopdd_savio.R")
library(parallel)

load("~/ElimFeas_StochMod/PostReview/Outputs/model_fit_profile_likelihood_parameters.Rdata")
load("~/ElimFeas_StochMod/PostReview/Refs/best_fit_params.Rdata")

par.sims = 50
  lam.range = seq(min(fin_pars95ci$lamda_twa), max(fin_pars95ci$lamda_twa), length.out = par.sims)  #transmission intensity range
  kap.range = seq(0, 2, length.out = par.sims)        #Pos. density dependence range

n.cores = detectCores() - 1

min.sim = 1001
max.sim = 1500

#Set parameters #######
covrg = 0.8
eff = 0.94
  mda.years = c(2:21)
  years = 61

#days where MDA is applied
  year.days = as.numeric()
  for(i in 1:20){
    year.days[i] = 365*i + (i-1)
  }
  
params = as.list(params)
  params$cov = covrg
  
#objects to fill
  fill = list()           #list to fill with simulations
  pe1 = as.numeric()      #elimination binary vector to fill for each sim
  w.pre = as.numeric()    #vector to fill with w_pre values
  w.pos = as.numeric()    #vector to fill with w_pos values
  bbr.w = as.numeric()    #vector to fill with bbr values from w_pre and w_pos  
  eps.w = as.numeric()    #vector to fill with epsilon (elim. feas. estimator) vals from slope of bbr
  prev.pre = as.numeric() #vector to fill with prevalence_pre values
  prev.pos = as.numeric() #vector to fill with prevalence_pos values
  bbr.prev = as.numeric() #vector to fill with prevalence based bbr values from w_pre and w_pos  
  eps.prev = as.numeric() #vector to fill with prevalence based epsilon (elim. feas. estimator) vals from slope of bbr


#Run simulations #######

#get equilibrium state variables    
load("~/ElimFeas_StochMod/PostReview/Refs/eqbm_stoch_sims2.Rdata")  
stoch_sims_eqbm_vals_fin = stoch_sims_eqbm_vals_fin[c(min.sim:max.sim),]
  
stoch.sims = 1000  #number of simulations for each parameter set

#Final values array to fill with p(e), eps, eps.sd
fill.arr = array(data = NA, dim = c(max.sim - min.sim + 1, 5, stoch.sims))

#Make cluster ######
clust = makeCluster(n.cores)
clusterExport(cl = clust, 
              varlist = c('params','lam.range', 'kap.range', 'years', 'stoch.sim.noise', 'stoch_sims_eqbm_vals_fin',
                          'mda.years', 'year.days', 'par.sims', 'stoch.sims','fill', 'transitions', 'sfx', 'eff',
                          'fill.arr', 'ssa.adaptivetau', 'covrg', 'w.pre', 'w.pos', 
                          'Prevalence', 'prev.pre', 'prev.pos', 'bbr.prev', 'eps.prev',
                          'phi_Wk', 'f_Wgk', 'R_Wv', 'pe1', 'bbr.w', 'eps.w', 'min.sim', 'max.sim'))  


#run all simulations######
for(s in c(min.sim:max.sim)){
  s2 = s - min.sim + 1

#Run sims for parameter set  
  fill.arr[s2, , ] = 
  parSapply(clust, c(1:stoch.sims), stoch.sim.noise, init = stoch_sims_eqbm_vals_fin[s2,c(3:7)], 
                                                lam = stoch_sims_eqbm_vals_fin[s2,2], 
                                                k = stoch_sims_eqbm_vals_fin[s2,1], simplify = T)[c(1,2,4,5,6),]
  
  print(s2)
}  

stopCluster(clust)

arr.name = paste('~/ElimFeas_StochMod/PostReview/Outputs/stoch_sims', min.sim, '_', max.sim, '.Rdata', sep = '')
save(fill.arr, file = arr.name)
