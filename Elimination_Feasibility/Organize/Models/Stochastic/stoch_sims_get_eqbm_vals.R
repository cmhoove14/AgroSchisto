require(deSolve)
require(graphics)
require(ggplot2)

source("Elimination_Feasibility/Organize/Models/Reff_BBR_fns.R")
source("Elimination_Feasibility/Organize/Models/schisto_mods_pdd_nopdd.R")

load("Elimination_Feasibility/Organize/Models/Outputs_Refs/model_fit_profile_likelihood_parameters.Rdata")
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/best_fit_params.Rdata")

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
  
  return(run[dim(run)[1],c(2:6)])
}

stoch_sims_eqbm_vals <- mapply(schisto_mod_pdd_add_ndds_k_lam_eqbm, k = par_grid[,1], lam = par_grid[,2])

stoch_sims_eqbm_vals_t <- t(stoch_sims_eqbm_vals)

stoch_sims_eqbm_vals_fin <- cbind(par_grid, stoch_sims_eqbm_vals_t)

save(stoch_sims_eqbm_vals_fin, file = "Elimination_Feasibility/Organize/Models/Outputs_Refs/eqbm_stoch_sims.Rdata")