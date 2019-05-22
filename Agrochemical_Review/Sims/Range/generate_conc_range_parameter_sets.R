require(tidyverse)

load("Agrochemical_Review/Models/trans_pars_fit.Rdata")
load("Agrochemical_Review/Models/fit_pars.Rdata")

source("Agrochemical_Review/Models/r0_functions.R")
source("Agrochemical_Review/Models/model_helper_functions.R")
source("Agrochemical_Review/Models/initial_parameters.R")

nil1 <- function(...){
  return(1)
}

nil0 <-function(...){
  return(0)
}

pars_df <- as.data.frame(as.list(fit_pars))

#Get csv file containing summary of all response functions
RFxSum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv") %>% 
  arrange(Study) %>% 
  group_by(Study) %>% 
  mutate(ID = row_number(),
         Study_ID_Unique = paste0(Study, "_", ID)) %>% ungroup()

#Source all .R scripts that contain derivation of response functions
RFx_files <- paste0("Agrochemical_Review/Response_Fxs/", 
                    list.files("Agrochemical_Review/Response_Fxs/")[grepl(".R", list.files("Agrochemical_Review/Response_Fxs/"))])

sapply(RFx_files, source)

#Function to get concentration range to test from eec
get_conc_range <- function(eec, n_range){
  seq(0, sqrt(eec*2), by = sqrt(eec*2)/n_range)^2
}

#Function to generate parameter sets to simulate net agrochemical effects resulting from multiple responses for particular agrochemical
gen_conc_range_par_set <- function(conc_range, nsims, 
                                   f.Knq = nil1, f.fnq = nil1, f.munq = nil0,  
                                   f.vq = nil1, f.pimq = nil1, f.thetaq = nil1, 
                                   f.picq = nil1, f.mupq = nil0, f.psiq = nil1){
  
  #create nsims duplicate rows of base parameters
  par_set <- pars_df %>% slice(rep(1:n(), each = nsims))
  
  par_set <- par_set %>% slice(rep(row_number(), length(conc_range))) %>% 
    mutate(conc = rep(conc_range, each = nsims),
           logp1_conc = log(conc+1),
        #Replace lambda transmission parameter with sample from vals within 95%CI
           lambda = sample(fit_1par_mat$lambda_1par, nrow(.), 
                           replace = T, prob = fit_1par_mat$weight),
        #Replace affected agrochemical parameters with estimates from input rfx
           Knq = map_dbl(conc, f.Knq),
           fnq = map_dbl(conc, f.fnq),
           munq = map_dbl(conc, f.munq),
           vq = map_dbl(conc, f.vq),
           pimq = map_dbl(conc, f.pimq),
           thetaq = map_dbl(conc, f.thetaq),
           picq = map_dbl(conc, f.picq),
           mupq = map_dbl(conc, f.mupq),
           psiq = map_dbl(conc, f.psiq))
  
  return(par_set)
  
}

#Global parameters for all chemicals
num_sims = 1000
n_conc = 200

set.seed(43093)

# Generate parameter set for Atrazine #############
atr_conc_range_pars <- gen_conc_range_par_set(conc_range = get_conc_range(RFxSum %>% 
                                                                            filter(Chemical == "Atrazine") %>% 
                                                                            slice(1) %>% 
                                                                            pull(eec), n_range = n_conc),
                                              nsims = num_sims,
                                              f.Knq = phiNq_atr_baxrohr.no30,
                                              f.munq = ons.muNq.atr,
                                              f.vq = halstead_meso18_atr_mans_v_uncertainty,
                                              f.picq = piC.atr.rohr08.lin)

#Save data frame of simulation parameter set
saveRDS(atr_conc_range_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                           "atrazine_net_parameters.rds"))

# Generate parameter set for Chlorpyrifos #############
chlor_conc_range_pars <- gen_conc_range_par_set(conc_range = get_conc_range(RFxSum %>% 
                                                                              filter(Chemical == "Chlorpyrifos") %>% 
                                                                              slice(1) %>% 
                                                                              pull(eec), n_range = n_conc),
                                                nsims = num_sims,
                                                f.fnq = fNq_chlor_ibr92_uncertainty, 
                                                f.munq = muNq_ch_hash11_uncertainty,  
                                                f.vq = halstead_meso18_chlor_mans_v_uncertainty, 
                                                f.pimq = piM_ch_Hash11_uncertainty,
                                                f.picq = piC_ch_Hash11_uncertainty, 
                                                f.mupq = muPq_chlor_mac_rohr_unpub_uncertainty)

#Save data frame of simulation parameter set
saveRDS(chlor_conc_range_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                             "chlorpyrifos_net_parameters.rds"))

# Generate parameter set for Glyphosate #############
gly_conc_range_pars <- gen_conc_range_par_set(conc_range = get_conc_range(RFxSum %>% 
                                                                            filter(Chemical == "Glyphosate") %>% 
                                                                            slice(1) %>% 
                                                                            pull(eec), n_range = n_conc),
                                              nsims = num_sims,
                                              f.Knq = phiNq_gly_baxrohr.no30,
                                              f.fnq = fNq.gly.fx.uncertainty, 
                                              f.munq = ons.muNq.gly,  
                                              f.pimq = piM.ghaf_gly.exp_unc,
                                              f.picq = piC.ghaf_gly.exp_unc)

#Save data frame of simulation parameter set
saveRDS(gly_conc_range_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                             "glyphosate_net_parameters.rds"))

# Generate parameter set for Malathion #############
mal_conc_range_pars <- gen_conc_range_par_set(conc_range = get_conc_range(RFxSum %>% 
                                                                            filter(Chemical == "Malathion") %>% 
                                                                            slice(1) %>% 
                                                                            pull(eec), n_range = n_conc),
                                              nsims = num_sims,
                                              f.fnq = fNq_mal_tch91_uncertainty, 
                                              f.munq = muNq_mal_tch91_uncertainty,  
                                              f.pimq = piM.tch91_mal_unc,
                                              f.picq = piC.tch92_mal_unc,
                                              f.mupq = muPq_mal_mac_rohr_unpub_uncertainty)
 
#Save data frame of simulation parameter set
saveRDS(mal_conc_range_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                             "malathion_net_parameters.rds"))

# Generate parameter set for Profenofos #############
prof_conc_range_pars <- gen_conc_range_par_set(conc_range = get_conc_range(RFxSum %>% 
                                                                            filter(Chemical == "Profenofos") %>% 
                                                                            slice(1) %>% 
                                                                            pull(eec), n_range = n_conc),
                                              nsims = num_sims,
                                              f.fnq = fNq_moh_prof_moh12_uncertainty, 
                                              f.munq = muNq_prof_mohamed_uncertainty,  
                                              f.pimq = piM_pr_Hash11_uncertainty,
                                              f.picq = piC_pr_Hash11_uncertainty,
                                              f.mupq = muPq_profenofos_Bajet12_uncertainty)

#Save data frame of simulation parameter set
saveRDS(prof_conc_range_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                             "profenofos_net_parameters.rds"))
