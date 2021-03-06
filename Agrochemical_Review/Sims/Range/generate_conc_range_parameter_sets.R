require(tidyverse)

load("Agrochemical_Review/Models/trans_pars_fit.Rdata")
load("Agrochemical_Review/Models/fit_pars.Rdata")

source("Agrochemical_Review/Models/r0_functions.R")
source("Agrochemical_Review/Models/model_helper_functions.R")
source("Agrochemical_Review/Models/initial_parameters.R")
source("Agrochemical_Review/Sims/Data/Obs_data_functions.R")

msqa2013_fin <- readRDS("Agrochemical_Review/Sims/Data/midwest_pesticides_2013.rds")
pestdat <- readRDS("Agrochemical_Review/Sims/Data/NAWQA_pesticides.rds")
ca_surf_clean <- readRDS("Agrochemical_Review/Sims/Data/CA_SURF_clean.rds")

#Get EEC values from EPA models
ins_eecs <- read_csv("Agrochemical_Review/Sims/Data/Insecticides_EECs.csv") %>% 
  dplyr::filter(X1 == "Peak EEC in a reservoir (MS) (ug/L)") %>% 
  gather("Chemical", "EPA_EEC_ppb", Carbaryl:Terbufos) %>% 
  dplyr::select(-X1)
  
hrb_eecs <- read_csv("Agrochemical_Review/Sims/Data/Herbicides_EECs.csv") %>% 
  dplyr::filter(X1 == "Peak EEC in a reservoir (MS) (ug/L)") %>% 
  gather("Chemical", "EPA_EEC_ppb", `2,4-D`:Trifluralin) %>% 
  dplyr::select(-X1)


nil1 <- function(...){
  return(1)
}

nil0 <-function(...){
  return(0)
}

pars_df <- as.data.frame(as.list(fit_pars))

#Get csv file containing summary of all response functions
RFxSum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv") %>% 
  left_join(bind_rows(ins_eecs, hrb_eecs), by = "Chemical") %>% 
  mutate(nawqa_peak = map_dbl(Chemical, get_nawqa_max),
         msqa2013_peak = map_dbl(Chemical, get_msqa2013_max)/1000,
         ca_surf_peak = map_dbl(Chemical, get_CA_SURF_max),
         peak_obs = pmax(nawqa_peak, msqa2013_peak, ca_surf_peak, na.rm = T),
         peak_eec_poc = pmax(peak_obs, as.numeric(EPA_EEC_ppb))) %>% 
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

# Generate parameter sets for Atrazine #############
# With EEC as peak concentration
atr_conc_range_pars <- gen_conc_range_par_set(conc_range = get_conc_range(RFxSum %>% 
                                                                            filter(Chemical == "Atrazine") %>% 
                                                                            slice(1) %>% 
                                                                            pull(peak_eec_poc), n_range = n_conc),
                                              nsims = num_sims,
                                              f.Knq = phiNq_atr_baxrohr.no30,
                                              f.munq = ons.muNq.atr,
                                              f.vq = halstead_meso18_atr_mans_v_uncertainty,
                                              f.picq = piC.atr.rohr08.lin)

#Save data frame of simulation parameter set
saveRDS(atr_conc_range_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                           "atrazine_net_parameters_eec.rds"))

# With EEC as ONLY concentration
atr_eec_pars <- gen_conc_range_par_set(conc_range = as.numeric(RFxSum %>% 
                                                filter(Chemical == "Atrazine") %>% 
                                                slice(1) %>% 
                                                pull(EPA_EEC_ppb)),
                                              nsims = num_sims,
                                              f.Knq = phiNq_atr_baxrohr.no30,
                                              f.munq = ons.muNq.atr,
                                              f.vq = halstead_meso18_atr_mans_v_uncertainty,
                                              f.picq = piC.atr.rohr08.lin)

#Save data frame of simulation parameter set
saveRDS(atr_eec_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                           "atrazine_eec_pars.rds"))

# With POC as ONLY concentration
atr_poc_pars <- gen_conc_range_par_set(conc_range = RFxSum %>% 
                                                filter(Chemical == "Atrazine") %>% 
                                                slice(1) %>% 
                                                pull(peak_obs),
                                              nsims = num_sims,
                                              f.Knq = phiNq_atr_baxrohr.no30,
                                              f.munq = ons.muNq.atr,
                                              f.vq = halstead_meso18_atr_mans_v_uncertainty,
                                              f.picq = piC.atr.rohr08.lin)

#Save data frame of simulation parameter set
saveRDS(atr_poc_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                           "atrazine_poc_pars.rds"))


# Generate parameter set for Chlorpyrifos #############
# With EEC as peak concentration
chlor_conc_range_pars <- gen_conc_range_par_set(conc_range = get_conc_range(RFxSum %>% 
                                                                              filter(Chemical == "Chlorpyrifos") %>% 
                                                                              slice(1) %>% 
                                                                              pull(peak_eec_poc), n_range = n_conc),
                                                nsims = num_sims,
                                                f.fnq = fNq_chlor_ibr92_uncertainty, 
                                                f.munq = muNq_ch_hash11_uncertainty,  
                                                f.vq = halstead_meso18_chlor_mans_v_uncertainty, 
                                                f.pimq = piM_ch_Hash11_uncertainty,
                                                f.picq = piC_ch_Hash11_uncertainty, 
                                                f.mupq = muPq_chlor_mac_rohr_unpub_uncertainty)

#Save data frame of simulation parameter set
saveRDS(chlor_conc_range_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                             "chlorpyrifos_net_parameters_eec.rds"))

# With EEC as ONLY concentration
chlor_eec_pars <- gen_conc_range_par_set(conc_range = as.numeric(RFxSum %>% 
                                           filter(Chemical == "Chlorpyrifos") %>% 
                                           slice(1) %>% 
                                           pull(EPA_EEC_ppb)),
                                         nsims = num_sims,
                                         f.fnq = fNq_chlor_ibr92_uncertainty, 
                                         f.munq = muNq_ch_hash11_uncertainty,  
                                         f.vq = halstead_meso18_chlor_mans_v_uncertainty, 
                                         f.pimq = piM_ch_Hash11_uncertainty,
                                         f.picq = piC_ch_Hash11_uncertainty, 
                                         f.mupq = muPq_chlor_mac_rohr_unpub_uncertainty)

#Save data frame of simulation parameter set
saveRDS(chlor_eec_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                      "chlorpyrifos_eec_pars.rds"))

# With POC as ONLY concentration
chlor_poc_pars <- gen_conc_range_par_set(conc_range = RFxSum %>% 
                                           filter(Chemical == "Chlorpyrifos") %>% 
                                           slice(1) %>% 
                                           pull(peak_obs),
                                         nsims = num_sims,
                                         f.fnq = fNq_chlor_ibr92_uncertainty, 
                                         f.munq = muNq_ch_hash11_uncertainty,  
                                         f.vq = halstead_meso18_chlor_mans_v_uncertainty, 
                                         f.pimq = piM_ch_Hash11_uncertainty,
                                         f.picq = piC_ch_Hash11_uncertainty, 
                                         f.mupq = muPq_chlor_mac_rohr_unpub_uncertainty)

#Save data frame of simulation parameter set
saveRDS(chlor_poc_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                      "chlorpyrifos_poc_pars.rds"))

# Generate parameter set for Glyphosate #############
#With EEC as peak concentration 
gly_conc_range_pars <- gen_conc_range_par_set(conc_range = get_conc_range(RFxSum %>% 
                                                                            filter(Chemical == "Glyphosate") %>% 
                                                                            slice(1) %>% 
                                                                            pull(peak_eec_poc), n_range = n_conc),
                                              nsims = num_sims,
                                              f.Knq = phiNq_gly_baxrohr.no30,
                                              f.fnq = fNq.gly.fx.uncertainty, 
                                              f.munq = ons.muNq.gly,  
                                              f.pimq = piM.ghaf_gly.exp_unc,
                                              f.picq = piC.ghaf_gly.exp_unc)

#Save data frame of simulation parameter set
saveRDS(gly_conc_range_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                             "glyphosate_net_parameters_eec.rds"))

# With EEC as ONLY concentration
gly_eec_pars <- gen_conc_range_par_set(conc_range = as.numeric(RFxSum %>% 
                                         filter(Chemical == "Glyphosate") %>% 
                                         slice(1) %>% 
                                         pull(EPA_EEC_ppb)),
                                       nsims = num_sims,
                                       f.Knq = phiNq_gly_baxrohr.no30,
                                       f.fnq = fNq.gly.fx.uncertainty, 
                                       f.munq = ons.muNq.gly,  
                                       f.pimq = piM.ghaf_gly.exp_unc,
                                       f.picq = piC.ghaf_gly.exp_unc)

#Save data frame of simulation parameter set
saveRDS(gly_eec_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                    "glyphosate_eec_pars.rds"))

# With POC as ONLY concentration
gly_poc_pars <- gen_conc_range_par_set(conc_range = RFxSum %>% 
                                         filter(Chemical == "Glyphosate") %>% 
                                         slice(1) %>% 
                                         pull(peak_obs),
                                       nsims = num_sims,
                                       f.Knq = phiNq_gly_baxrohr.no30,
                                       f.fnq = fNq.gly.fx.uncertainty, 
                                       f.munq = ons.muNq.gly,  
                                       f.pimq = piM.ghaf_gly.exp_unc,
                                       f.picq = piC.ghaf_gly.exp_unc)

#Save data frame of simulation parameter set
saveRDS(gly_poc_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                    "glyphosate_poc_pars.rds"))

# Generate parameter set for Malathion #############
#With EEC as peak concentration 
mal_conc_range_pars <- gen_conc_range_par_set(conc_range = get_conc_range(RFxSum %>% 
                                                                            filter(Chemical == "Malathion") %>% 
                                                                            slice(1) %>% 
                                                                            pull(peak_eec_poc), n_range = n_conc),
                                              nsims = num_sims,
                                              f.fnq = fNq_mal_tch91_uncertainty, 
                                              f.munq = muNq_mal_tch91_uncertainty,  
                                              f.pimq = piM.tch91_mal_unc,
                                              f.picq = piC.tch92_mal_unc,
                                              f.mupq = muPq_mal_mac_rohr_unpub_uncertainty)
 
#Save data frame of simulation parameter set
saveRDS(mal_conc_range_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                             "malathion_net_parameters_eec.rds"))

# With EEC as ONLY concentration
mal_eec_pars <- gen_conc_range_par_set(conc_range = as.numeric(RFxSum %>% 
                                         filter(Chemical == "Malathion") %>% 
                                         slice(1) %>% 
                                         pull(EPA_EEC_ppb)),
                                       nsims = num_sims,
                                       f.fnq = fNq_mal_tch91_uncertainty, 
                                       f.munq = muNq_mal_tch91_uncertainty,  
                                       f.pimq = piM.tch91_mal_unc,
                                       f.picq = piC.tch92_mal_unc,
                                       f.mupq = muPq_mal_mac_rohr_unpub_uncertainty)

#Save data frame of simulation parameter set
saveRDS(mal_eec_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                    "malathion_eec_pars.rds"))

# With POC as ONLY concentration
mal_poc_pars <- gen_conc_range_par_set(conc_range = RFxSum %>% 
                                         filter(Chemical == "Malathion") %>% 
                                         slice(1) %>% 
                                         pull(peak_obs),
                                       nsims = num_sims,
                                       f.fnq = fNq_mal_tch91_uncertainty, 
                                       f.munq = muNq_mal_tch91_uncertainty,  
                                       f.pimq = piM.tch91_mal_unc,
                                       f.picq = piC.tch92_mal_unc,
                                       f.mupq = muPq_mal_mac_rohr_unpub_uncertainty)

#Save data frame of simulation parameter set
saveRDS(mal_poc_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                    "malathion_poc_pars.rds"))

# Generate parameter set for Profenofos #############
# With EEC as peak concentration
prof_conc_range_pars <- gen_conc_range_par_set(conc_range = get_conc_range(RFxSum %>% 
                                                                            filter(Chemical == "Profenofos") %>% 
                                                                            slice(1) %>% 
                                                                            pull(peak_eec_poc), n_range = n_conc),
                                              nsims = num_sims,
                                              f.fnq = fNq_moh_prof_moh12_uncertainty, 
                                              f.munq = muNq_prof_mohamed_uncertainty,  
                                              f.pimq = piM_pr_Hash11_uncertainty,
                                              f.picq = piC_pr_Hash11_uncertainty,
                                              f.mupq = muPq_profenofos_Bajet12_uncertainty)

#Save data frame of simulation parameter set
saveRDS(prof_conc_range_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                             "profenofos_net_parameters_eec.rds"))

# With EEC as ONLY concentration
prof_eec_pars <- gen_conc_range_par_set(conc_range = as.numeric(RFxSum %>% 
                                          filter(Chemical == "Profenofos") %>% 
                                          slice(1) %>% 
                                          pull(EPA_EEC_ppb)),
                                        nsims = num_sims,
                                        f.fnq = fNq_moh_prof_moh12_uncertainty, 
                                        f.munq = muNq_prof_mohamed_uncertainty,  
                                        f.pimq = piM_pr_Hash11_uncertainty,
                                        f.picq = piC_pr_Hash11_uncertainty,
                                        f.mupq = muPq_profenofos_Bajet12_uncertainty)

#Save data frame of simulation parameter set
saveRDS(prof_eec_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                     "profenofos_eec_pars.rds"))

# With POC as ONLY concentration
prof_poc_pars <- gen_conc_range_par_set(conc_range = RFxSum %>% 
                                          filter(Chemical == "Profenofos") %>% 
                                          slice(1) %>% 
                                          pull(peak_obs),
                                        nsims = num_sims,
                                        f.fnq = fNq_moh_prof_moh12_uncertainty, 
                                        f.munq = muNq_prof_mohamed_uncertainty,  
                                        f.pimq = piM_pr_Hash11_uncertainty,
                                        f.picq = piC_pr_Hash11_uncertainty,
                                        f.mupq = muPq_profenofos_Bajet12_uncertainty)

#Save data frame of simulation parameter set
saveRDS(prof_poc_pars, file = paste0("Agrochemical_Review/Sims/Range/Parameter_Sets/", 
                                     "profenofos_poc_pars.rds"))
