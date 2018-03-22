require(deSolve)
require(graphics)
require(ggplot2)
require(parallel)

source("~/ElimFeas_StochMod/PostReview/Refs/stochastic_model_addNDD_savio.R")
source("~/ElimFeas_StochMod/PostReview/Refs/schisto_mods_pdd_nopdd_savio.R")

load("~/ElimFeas_StochMod/PostReview/Outputs/model_fit_profile_likelihood_parameters.Rdata")
load("~/ElimFeas_StochMod/PostReview/Refs/best_fit_params.Rdata")

n.cores = detectCores() - 1

par.sims = 50
  lam.range = seq(min(fin_pars95ci$lamda_twa), max(fin_pars95ci$lamda_twa), length.out = par.sims)  #transmission intensity range
  kap.range = seq(0, 2, length.out = par.sims)        #Pos. density dependence range

par_grid <- expand.grid(k = kap.range, lambda = lam.range)

t = seq(0, 365*200, 30)

start=c(S=5000,E=1000,I=500, Wt=50, Wu=50)

#Customize model function to take inputs of kappa and lambda and return equilibrium values
schisto_mod_pdd_add_ndds_k_lam_eqbm <- function(k, lam, time = t, nstart = start){
  
  params["k"] <- k
  params["lamda"] <- lam
  
  run <- ode(nstart, time, schisto_mod_pdd_add_ndds, params)
  
  if(k==0) print(run[dim(run)[1],c(2:6)])

  return(run[dim(run)[1],c(2:6)])
}

clust = makeCluster(n.cores)
clusterExport(cl = clust, 
              varlist = c('params','lam.range', 'kap.range', 'schisto_mod_pdd_add_ndds_k_lam_eqbm', 
                          'par.sims', 'start','t', 'par_grid'))  

stoch_sims_eqbm_vals <- mcmapply(schisto_mod_pdd_add_ndds_k_lam_eqbm, k = par_grid[,1], lam = par_grid[,2])


stoch_sims_eqbm_vals_t <- t(stoch_sims_eqbm_vals)

stoch_sims_eqbm_vals_fin <- cbind(par_grid, stoch_sims_eqbm_vals_t)

stopCluster(clust)

save(stoch_sims_eqbm_vals_fin, file = "~/ElimFeas_StochMod/PostReview/Outputs/eqbm_stoch_sims2.Rdata")