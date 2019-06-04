#agrochemical forced model with NO PREDATORS for dynamic agrochemical concentration through time
agrochem_forced_mda_mod_no_preds = function(t, n, parameters, force_fx, 
                                            f.Knq = nil1, f.fnq = nil1, f.munq = nil0,  
                                            f.vq = nil1, f.pimq = nil1, f.thetaq = nil1, 
                                            f.picq = nil1) {
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
  
  #Dynamic variables
    N = S+E+I  #Total number of snails
    
    W = parameters["cvrg"]*Wt + (1-parameters["cvrg"])*Wu
    
    phi = phi_Wk(W, parameters["kappa"])           

  #Agrochemical concentration from forcing function    
    q = force_fx(t)
#agrochemical parameters adjusted
    Knadj = parameters["K_N"] * f.Knq(q)
    fnadj = parameters["f_N"] * f.fnq(q)
    munadj = parameters["mu_N"] + f.munq(q)
    vadj = parameters["v"] * f.vq(q)
    pimadj = parameters["pi_M"] * f.pimq(q)
    thetaadj = parameters["theta"] * f.thetaq(q)
    picadj = parameters["pi_C"] * f.picq(q)

#Schistosome larval concentration equations
    M = 0.5*W*parameters["H"]*phi*parameters["m"]*parameters["v"]*parameters["pi_M"] #- (mu_M + mu_Mq)*M - beta*M*N 
    C = parameters["theta"]*I*parameters["pi_C"] #- (mu_C + mu_Cq)*C - lambda*H*C
    
    #differential equations
    
    dSdt= fnadj*(1-(N/Knadj))*(S+E) - munadj*S - parameters["beta"]*M*S  #Susceptible snails
    
    dEdt= parameters["beta"]*M*S - munadj*E - parameters["sigma"]*E                                     #Exposed snails
    
    dIdt= parameters["sigma"]*E - (munadj + parameters["mu_I"])*I                                       #Infected snails
    
    dWtdt= parameters["lambda"]*C - (parameters["mu_W"]+parameters["mu_H"])*Wt  #mean worm burden in treated human population
    dWudt= parameters["lambda"]*C - (parameters["mu_W"]+parameters["mu_H"])*Wu  #mean worm burden in untreated human population

    
    return(list(c(dSdt, dEdt, dIdt, dWtdt, dWudt)))
}

#agrochemical forced model for dynamic agrochemical concentration through time
agrochem_forced_mda_mod = function(t, n, parameters, force_fx, 
                                   f.Knq = nil1, f.fnq = nil1, f.munq = nil0,  
                                   f.vq = nil1, f.pimq = nil1, f.thetaq = nil1, 
                                   f.picq = nil1, f.mupq = nil0, f.psiq = nil1) {
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    P=n[6]
    
  #Dynamic variables
    N = S+E+I  #Total number of snails
    
    W = parameters["cvrg"]*Wt + (1-parameters["cvrg"])*Wu
    
    phi = phi_Wk(W, parameters["kappa"])           

    pred_S = ((parameters["alpha"]*parameters["eps"])*(S/parameters["A"])^parameters["nn"]) / 
      (1 + sum((parameters["alpha"]*parameters["eps"])*parameters["Th"]*(S/parameters["A"])^parameters["nn"],
               (parameters["alpha"]*parameters["eps"])*parameters["Th"]*(E/parameters["A"])^parameters["nn"],
               (parameters["alpha"]*parameters["eps"])*parameters["Th"]*(I/parameters["A"])^parameters["nn"]))
    
    pred_E = ((parameters["alpha"]*parameters["eps"])*(E/parameters["A"])^parameters["nn"]) / 
      (1 + sum((parameters["alpha"]*parameters["eps"])*parameters["Th"]*(S/parameters["A"])^parameters["nn"],
               (parameters["alpha"]*parameters["eps"])*parameters["Th"]*(E/parameters["A"])^parameters["nn"],
               (parameters["alpha"]*parameters["eps"])*parameters["Th"]*(I/parameters["A"])^parameters["nn"]))
    
    pred_I = ((parameters["alpha"]*parameters["eps"])*(I/parameters["A"])^parameters["nn"]) / 
      (1 + sum((parameters["alpha"]*parameters["eps"])*parameters["Th"]*(S/parameters["A"])^parameters["nn"],
               (parameters["alpha"]*parameters["eps"])*parameters["Th"]*(E/parameters["A"])^parameters["nn"],
               (parameters["alpha"]*parameters["eps"])*parameters["Th"]*(I/parameters["A"])^parameters["nn"]))    
  
  #Agrochemical concentration from forcing function    
    q = force_fx(t)
#agrochemical parameters adjusted
    Knadj = parameters["K_N"] * f.Knq(q)
    fnadj = parameters["f_N"] * f.fnq(q)
    munadj = parameters["mu_N"] + f.munq(q)
    vadj = parameters["v"] * f.vq(q)
    pimadj = parameters["pi_M"] * f.pimq(q)
    thetaadj = parameters["theta"] * f.thetaq(q)
    picadj = parameters["pi_C"] * f.picq(q)
    mupadj = parameters["mu_P"] + f.mupq(q)
    psiadj = f.psiq(q)
        
    pred_Sadj = pred_S * psiadj
    pred_Eadj = pred_E * psiadj
    pred_Iadj = pred_I * psiadj
    
#Schistosome larval concentration equations
    M = 0.5*W*parameters["H"]*phi*parameters["m"]*parameters["v"]*parameters["pi_M"] #- (mu_M + mu_Mq)*M - beta*M*N 
    C = parameters["theta"]*I*parameters["pi_C"] #- (mu_C + mu_Cq)*C - lambda*H*C
    
    #differential equations
    
    dSdt= fnadj*(1-(N/Knadj))*(S+E) - munadj*S - pred_Sadj*P - parameters["beta"]*M*S  #Susceptible snails
    
    dEdt= parameters["beta"]*M*S - munadj*E - pred_Eadj*P - parameters["sigma"]*E                                     #Exposed snails
    
    dIdt= parameters["sigma"]*E - (munadj + parameters["mu_I"])*I - pred_Iadj*P                                       #Infected snails
    
    dWtdt= parameters["lambda"]*C - (parameters["mu_W"]+parameters["mu_H"])*Wt  #mean worm burden in treated human population
    dWudt= parameters["lambda"]*C - (parameters["mu_W"]+parameters["mu_H"])*Wu  #mean worm burden in untreated human population

    dPdt= parameters["f_P"]*(1-(P/parameters["K_P"]))*P - mupadj*P                                      #prawn population (number individuals)

    
    return(list(c(dSdt, dEdt, dIdt, dWtdt, dWudt, dPdt)))
}

#Function to simulate 4 scenarios: MDA+agrochem pollution, MDA only, Agrochem pollution only, and neither MDA or agrochemical pollution and compare DALYs accumulated in each
sim_scenarios <- function(pars, years, force_fx, 
                          f.Knq = nil1, f.fnq = nil1, f.munq = nil0,  
                          f.vq = nil1, f.pimq = nil1, f.thetaq = nil1, 
                          f.picq = nil1, f.mupq = nil0, f.psiq = nil1){
  
#Run to equilibirum with input parameter set and initial conditions ############
  eqbm_vals <- runsteady(y = init_vals, func = mda_mod, 
                         parms = pars)$y
    names(eqbm_vals) <- c("S", "E", "I", "Wt", "Wu", "P")

  eqbm_vals_no_preds <- runsteady(y = init_vals_no_preds, func = mda_pred_free_mod, 
                         parms = pars)$y
    names(eqbm_vals_no_preds) <- c("S", "E", "I", "Wt", "Wu")
    
  #Simulation with NOTHING #############
    sim_nothing <- as.data.frame(ode(y = eqbm_vals_no_preds,
                                     times = t_sim,
                                     func = mda_pred_free_mod,
                                     parms = pars,
                                     method = "euler")) %>% 
      mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
             prev = get_prev(pars["kappa"], W),
             DALYs_Wt = map_dbl(Wt, est_dalys, 
                                pars["kappa"], weight_lo, weight_hi, epmL, 
                                round(pars["H"] * pars["cvrg"])),
             DALYs_Wu = map_dbl(Wu, est_dalys, 
                                pars["kappa"], weight_lo, weight_hi, epmL, 
                                round(pars["H"] * (1-pars["cvrg"]))),
             DALYs_W = DALYs_Wt + DALYs_Wu)

  #Simulation with PREDS ONLY ###############
    sim_preds <- as.data.frame(ode(y = eqbm_vals,
                                   times = t_sim,
                                   func = mda_mod,
                                   parms = pars,
                                   method = "euler")) %>% 
      mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
             prev = get_prev(pars["kappa"], W),
             DALYs_Wt = map_dbl(Wt, est_dalys, 
                                pars["kappa"], weight_lo, weight_hi, epmL, 
                                round(pars["H"] * pars["cvrg"])),
             DALYs_Wu = map_dbl(Wu, est_dalys, 
                                pars["kappa"], weight_lo, weight_hi, epmL, 
                                round(pars["H"] * (1-pars["cvrg"]))),
             DALYs_W = DALYs_Wt + DALYs_Wu)
    
  #Simulation with ANNUAL MDA ONLY ##############
    sim_mda <- as.data.frame(ode(y = eqbm_vals_no_preds,
                                 times = t_sim,
                                 func = mda_pred_free_mod,
                                 parms = pars,
                                 events = list(data = mdas),
                                 method = "euler")) %>% 
      mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
             prev = get_prev(pars["kappa"], W),
             DALYs_Wt = map_dbl(Wt, est_dalys, 
                                pars["kappa"], weight_lo, weight_hi, epmL, 
                                round(pars["H"] * pars["cvrg"])),
             DALYs_Wu = map_dbl(Wu, est_dalys, 
                                pars["kappa"], weight_lo, weight_hi, epmL, 
                                round(pars["H"] * (1-pars["cvrg"]))),
             DALYs_W = DALYs_Wt + DALYs_Wu)
    
  #Simulation with ANNUAL MDA AND PREDATORS  ###############
    sim_preds_mda <- as.data.frame(ode(y = eqbm_vals,
                                       times = t_sim,
                                       func = mda_mod,
                                       parms = pars,
                                       events = list(data = mdas),
                                       method = "euler")) %>% 
      mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
             prev = get_prev(pars["kappa"], W),
             DALYs_Wt = map_dbl(Wt, est_dalys, 
                                pars["kappa"], weight_lo, weight_hi, epmL, 
                                round(pars["H"] * pars["cvrg"])),
             DALYs_Wu = map_dbl(Wu, est_dalys, 
                                pars["kappa"], weight_lo, weight_hi, epmL, 
                                round(pars["H"] * (1-pars["cvrg"]))),
             DALYs_W = DALYs_Wt + DALYs_Wu)
    
  #simulation with AGRO ONLY ##################
    sim_agro <- as.data.frame(ode(y = eqbm_vals_no_preds,
                              times = t_sim,
                              func = agrochem_forced_mda_mod_no_preds,
                              force_fx = force_fx,
                              f.Knq = f.Knq, f.fnq = f.fnq, f.munq = f.munq,  
                              f.vq = f.vq, f.pimq = f.pimq, f.thetaq = f.thetaq, 
                              f.picq = f.picq,
                              parms = pars,
                              events = list(data = data.frame(var = "S", #Events dataframe that doesn't do anything
                                                              time = 1, 
                                                              value = 1, 
                                                              method = "multiply")),
                              method = "euler")) %>% 
    mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
           prev = get_prev(pars["kappa"], W),
           DALYs_Wt = map_dbl(Wt, est_dalys, 
                              pars["kappa"], weight_lo, weight_hi, epmL, 
                              round(pars["H"] * pars["cvrg"])),
           DALYs_Wu = map_dbl(Wu, est_dalys, 
                              pars["kappa"], weight_lo, weight_hi, epmL, 
                              round(pars["H"] * (1-pars["cvrg"]))),
           DALYs_W = DALYs_Wt + DALYs_Wu)
    
  #Simulation with AGRO AND MDA #####################  
    sim_agro_mda <- as.data.frame(ode(y = eqbm_vals_no_preds,
                              times = t_sim,
                              func = agrochem_forced_mda_mod_no_preds,
                              force_fx = force_fx,
                              f.Knq = f.Knq, f.fnq = f.fnq, f.munq = f.munq,  
                              f.vq = f.vq, f.pimq = f.pimq, f.thetaq = f.thetaq, 
                              f.picq = f.picq,
                              parms = pars,
                              events = list(data = mdas),
                              method = "euler")) %>% 
    mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
           prev = get_prev(pars["kappa"], W),
           DALYs_Wt = map_dbl(Wt, est_dalys, 
                              pars["kappa"], weight_lo, weight_hi, epmL, 
                              round(pars["H"] * pars["cvrg"])),
           DALYs_Wu = map_dbl(Wu, est_dalys, 
                              pars["kappa"], weight_lo, weight_hi, epmL, 
                              round(pars["H"] * (1-pars["cvrg"]))),
           DALYs_W = DALYs_Wt + DALYs_Wu)
    
  #Simulation with AGRO and PREDS   #################
    sim_agro_preds <- as.data.frame(ode(y = eqbm_vals,
                              times = t_sim,
                              func = agrochem_forced_mda_mod,
                              force_fx = force_fx,
                              f.Knq = f.Knq, f.fnq = f.fnq, f.munq = f.munq,  
                              f.vq = f.vq, f.pimq = f.pimq, f.thetaq = f.thetaq, 
                              f.picq = f.picq, f.mupq = f.mupq, f.psiq = f.psiq,
                              parms = pars,
                              events = list(data = data.frame(var = "S", #Events dataframe that doesn't do anything
                                                           time = 1, 
                                                           value = 1, 
                                                           method = "multiply")),
                              method = "euler")) %>% 
    mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
           prev = get_prev(pars["kappa"], W),
           DALYs_Wt = map_dbl(Wt, est_dalys, 
                              pars["kappa"], weight_lo, weight_hi, epmL, 
                              round(pars["H"] * pars["cvrg"])),
           DALYs_Wu = map_dbl(Wu, est_dalys, 
                              pars["kappa"], weight_lo, weight_hi, epmL, 
                              round(pars["H"] * (1-pars["cvrg"]))),
           DALYs_W = DALYs_Wt + DALYs_Wu)
    
  #Simulation with AGRO PREDS AND MDA #################
    sim_agro_preds_mda <- as.data.frame(ode(y = eqbm_vals,
                              times = t_sim,
                              func = agrochem_forced_mda_mod,
                              force_fx = force_fx,
                              f.Knq = f.Knq, f.fnq = f.fnq, f.munq = f.munq,  
                              f.vq = f.vq, f.pimq = f.pimq, f.thetaq = f.thetaq, 
                              f.picq = f.picq, f.mupq = f.mupq, f.psiq = f.psiq,
                              parms = pars,
                              events = list(data = mdas),
                              method = "euler")) %>% 
    mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
           prev = get_prev(pars["kappa"], W),
           DALYs_Wt = map_dbl(Wt, est_dalys, 
                              pars["kappa"], weight_lo, weight_hi, epmL, 
                              round(pars["H"] * pars["cvrg"])),
           DALYs_Wu = map_dbl(Wu, est_dalys, 
                              pars["kappa"], weight_lo, weight_hi, epmL, 
                              round(pars["H"] * (1-pars["cvrg"]))),
           DALYs_W = DALYs_Wt + DALYs_Wu)
  
    
    return(list("Nothing" = as.list(sim_nothing),
                "Preds_only" = as.list(sim_preds),
                "MDA_only" = as.list(sim_mda),
                "MDA_Preds" = as.list(sim_preds_mda),
                "Agro_only" = as.list(sim_agro),
                "Agro_MDA" = as.list(sim_agro_mda),
                "Agro_Preds" = as.list(sim_agro_preds),
                "Agro_Preds_MDA" = as.list(sim_agro_preds_mda)))
}

#Function to get total dalys accumulated in a sim from full simulation list
get_dalys_from_sim <- function(sim_list, sim_type, dalys_col){
  sapply(1:length(sim_list), function(x){
    sim_df <- as.data.frame(sim_list[[x]][[sim_type]])
    sum(sim_df[dalys_col][which(sim_df["time"] >365),])
  })
}

# Function that uses get_dalys_from_sim above to go through all simulations to get DALYs estimates and return dataframe
comp_dalys_across_sims <- function(simulation_list){
  #Estimate cumulative DALYs during ten year simulation periods    ################
    dalys_nothing <- get_dalys_from_sim(simulation_list, "Nothing", "DALYs_W")
    dalys_preds <- get_dalys_from_sim(simulation_list, "Preds_only", "DALYs_W")
    dalys_mda <- get_dalys_from_sim(simulation_list, "MDA_only", "DALYs_W")
    dalys_preds_mda <-  get_dalys_from_sim(simulation_list, "MDA_Preds", "DALYs_W")
    dalys_agro <- get_dalys_from_sim(simulation_list, "Agro_only", "DALYs_W")
    dalys_agro_mda <- get_dalys_from_sim(simulation_list, "Agro_MDA", "DALYs_W")
    dalys_agro_preds <- get_dalys_from_sim(simulation_list, "Agro_Preds", "DALYs_W")
    dalys_agro_preds_mda <- get_dalys_from_sim(simulation_list, "Agro_Preds_MDA", "DALYs_W")

    return(data.frame("Nothing" = dalys_nothing,
                      "Preds_only" = dalys_preds,
                      "MDA_only" = dalys_mda,
                      "MDA_Preds" = dalys_preds_mda,
                      "Agro_only" = dalys_agro,
                      "Agro_MDA" = dalys_agro_mda,
                      "Agro_Preds" = dalys_agro_preds,
                      "Agro_Preds_MDA" = dalys_agro_preds_mda))  
}

# Function to get time series of state variables from full simulation list
get_ts_from_sim <- function(simulation_list, simulation_run){
  #Get names of all scenarios
    scenario_names <- names(simulation_list[[simulation_run]])
  
  # Apply function that puts data into long format to every simulation type then bind together  
  scenario_ts <- lapply(scenario_names, function(name){
    as.data.frame(simulation_list[[simulation_run]][[name]]) %>% 
      gather("Variable", "Value", S:DALYs_W) %>% 
      mutate(sim_type = name,
             sim_run = simulation_run)
  })
  
  sim_full_df <- do.call(rbind, scenario_ts)
  
  return(sim_full_df)
}

# Function that uses above function that returns time series from single simulation to take simulation list and return dataframe of time series from all simulations
comp_ts_across_sims <- function(simulation_list){
  do.call(rbind, lapply(1:length(simulation_list),
                        get_ts_from_sim, simulation_list = simulation_list))
}
