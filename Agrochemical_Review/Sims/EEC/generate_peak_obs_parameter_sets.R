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

pars_df <- as.data.frame(as.list(fit_pars))

#Get csv file containing summary of all response functions
RFxSum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

#Add estimates of peak concentrations for each chem found in each observed dataset and add unique study_IDs
RFxSum <- RFxSum %>% 
  mutate(nawqa_peak = map_dbl(Chemical, get_nawqa_max),
         msqa2013_peak = map_dbl(Chemical, get_msqa2013_max)/1000,
         ca_surf_peak = map_dbl(Chemical, get_CA_SURF_max),
         peak_observed = pmax(nawqa_peak, msqa2013_peak, ca_surf_peak, na.rm = T)) %>% 
  arrange(Study) %>% 
  group_by(Study) %>% 
  mutate(ID = row_number(),
         Study_ID_Unique = paste0(Study, "_", ID)) %>% ungroup()

# Some tweaks to remove NAs and allow ammonium fertilizer studies that were performed at mesocosm concentrations
RFxSum$peak_observed[RFxSum$Chemical == "Ammonium Fertilizer"] <- RFxSum$eec[RFxSum$Chemical == "Ammonium Fertilizer"][1]

#Source all .R scripts that contain derivation of response functions
RFx_files <- paste0("Agrochemical_Review/Response_Fxs/", 
                    list.files("Agrochemical_Review/Response_Fxs/")[grepl(".R", list.files("Agrochemical_Review/Response_Fxs/"))])

sapply(RFx_files, source)

#Function to generate parameter sets to simulate agrochemical effects at EEC or other concentration for each response function
gen_peak_conc_par_set <- function(peak_conc, nsims, rfx, par, study){
  
  #create nsims duplicate rows of base parameters
  par_set <- pars_df %>% slice(rep(1:n(), each = nsims))
  
  #Replace lambda transmission parameter with sample from vals within 95%CI
  par_set["lambda"] <- sample(fit_1par_mat$lambda_1par, nsims, 
                              replace = T, prob = fit_1par_mat$weight)
  
  #Replace affected agrochemical parameter with estimates from rfx
  par_set[par] <- sapply(rep(peak_conc, nsims), rfx)
  
  par_set["Study_ID_Unique"] <- study

  #Save data frame of simulation parameter set
  saveRDS(par_set, file = paste0("Agrochemical_Review/Sims/EEC/Parameter_Sets/", 
                                 study, 
                                 "_peak_conc_parameters.rds"))
  
}

#Do this for all response functions with relevant parameters affected
RFxSum_Sub <- RFxSum %>% filter(parameter %in% colnames(pars_df) & !is.na(rfx) & !is.na(peak_observed))
 
num_sims = 1000

set.seed(43093)

for(i in 1:nrow(RFxSum_Sub)){
  gen_peak_conc_par_set(peak_conc = RFxSum_Sub$peak_observed[i], 
                        nsims = num_sims, 
                        rfx = RFxSum_Sub$rfx[i], 
                        par = RFxSum_Sub$parameter[i],
                        study = RFxSum_Sub$Study_ID_Unique[i])
  
  if(i %% 20 == 0) print(i)
} 