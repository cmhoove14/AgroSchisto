# Estimating the elimination feasibility in the ‘end game’ of control efforts for parasites subjected to regular mass drug administration: a schistosomiasis case study  

## `Outputs_Refs` folder  
Contains outputs of simulations and parameter sets used in calls in plotting scripts used to produce manuscript figures  

## `Savio` folder  
Contains scripts with minor modifications that allow them to be run on savio, the Berkeley supercomputing cluster

## `Stochastic` folder  
Stochastic versions of the models 

## Scripts contained in this folder (in order that they would be run to reproduce results presented in the manuscript)  
+ `schisto_mods_pdd_nopdd.R` -- Deterministic models used to simulate schisto dynamics in the presence of different density dependencies. Also contain best-fitting parameter set   
+ `model_fit_profile.R` -- Fits the model with density dependencies to observed data from a study site in Senegal and stores the resulting transmission parameter sets in `Outputs_Refs` folder as `model_fit_profile_likelihood_parameters.Rdata` (run on savio with savio version)  
+ `Reff_BBR_fns` -- Functions to estimate Reff and BBR, derived from models in `schisto_mods_pdd_nopdd.R`  
+ `Reff_addNDDs_BBR_plots.R` -- Estimates Reff, BBR, and dW/dt across a range of worm burdens for the best fitting parameter set. Used to produce Figures 1 and 2 that demonstrate difference between model predictions with and without influence of positive density dependence  
+ `Get_eqbm_for_trans_parameter_sets.R` -- For each transmission parameter set, simulates pdd and pdd-free models over long time period to reach equilibrium and saves the resulting equilibrium values in `Outputs_Refs` folder  
+ `PDD_noPDD_sims.R` -- For each transmission parameter set and each model: Gets equilibrium estimates produced with `Get_eqbm_for_trans_parameter_sets.R` and simulates 20 years of annual MDA, then simulates 40 years of intervention-free transmission. During the 20 years of MDA, estimates BBR between each MDA, epsilon over entire MDA sequence for that transmission parameter set. Stores all outputs in `Outputs_Refs` folder  
+ `PDD_noPDD_plots_analyze` -- Processes simulations produced in `PDD_noPDD_sims.R` to produce Figure 3 showing (a) worm burden trajectory for the pdd and pdd-free models over the course of simulated 20 years MDA and 40 years intervention-free transmission (b) estimates of BBR in each MDA interval for each model and (c) Estimate of epsilon for each year >3 in the 20 years of MDA. Also produces Figure X which shows the same three panels but based on prevalence data rather than worm burden data  
+ `stochModSims_postProcess` -- takes stochastic simulations run on Savio and process them, plots results shown in Figure 4 of relationship between probabilioty of elimination (P(e)) from stochastic simulations to elimination feasibility coefficient/estimator (epsilon)  
+ `eps_surface_kap_lam.R` -- For a range of positive density dependence (k) and transmission intensity (lambda) estimates epsilon over 20 years of simulated MDA  
+ `eps_surface_kap_lam_plots.R` -- Takes output from `eps_surface_kap_lam.R` and produces a surface plot (Figure 5) showing how epsilon varies across the two key parameters  
