require(tidyverse)

load("Agrochemical_Review/Models/trans_pars_fit.Rdata")
load("Agrochemical_Review/Models/fit_pars.Rdata")

source("Agrochemical_Review/Models/r0_functions.R")
source("Agrochemical_Review/Models/model_helper_functions.R")
source("Agrochemical_Review/Models/initial_parameters.R")

pars_df <- as.data.frame(as.list(fit_pars))

#Get csv file containing summary of all response functions
RFxSum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

RFxSum <- RFxSum %>% 
  arrange(Study) %>% 
  group_by(Study) %>% 
  mutate(ID = row_number(),
         Study_ID_Unique = paste0(Study, "_", ID)) %>% ungroup()

#Source all .R scripts that contain derivation of response functions
RFx_files <- paste0("Agrochemical_Review/Response_Fxs/", 
                    list.files("Agrochemical_Review/Response_Fxs/")[grepl(".R", list.files("Agrochemical_Review/Response_Fxs/"))])

sapply(RFx_files, source)

#Function to generate parameter sets to simulate agrochemical effects at EEC for each response function
gen_eec_par_set <- function(eec, nsims, rfx, par, study){
  
  #create nsims duplicate rows of base parameters
  par_set <- pars_df %>% slice(rep(1:n(), each = nsims))
  
  #Replace lambda transmission parameter with sample from vals within 95%CI
  par_set["lambda"] <- sample(fit_1par_mat$lambda_1par, nsims, 
                              replace = T, prob = fit_1par_mat$weight)
  
  #Replace affected agrochemical parameter with estimates from rfx
  par_set[par] <- sapply(rep(eec, nsims), rfx)
  
  par_set["Study_ID_Unique"] <- study

  #Save data frame of simulation parameter set
  saveRDS(par_set, file = paste0("Agrochemical_Review/Sims/EEC/Parameter_Sets/", 
                                 study, 
                                 "_eec_parameters.rds"))
  
}

#Do this for all response functions with relevant parameters affected
RFxSum_Sub <- RFxSum %>% filter(parameter %in% colnames(pars_df) & !is.na(rfx))
 
num_sims = 1000

set.seed(43093)

for(i in 1:nrow(RFxSum_Sub)){
  gen_eec_par_set(eec = RFxSum_Sub$eec[i], 
                  nsims = num_sims, 
                  rfx = RFxSum_Sub$rfx[i], 
                  par = RFxSum_Sub$parameter[i],
                  study = RFxSum_Sub$Study_ID_Unique[i])
  
  if(i %% 20 == 0) print(i)
} 

#Also add dataset with no response function for baseline compariosn ############
  #create nsims duplicate rows of base parameters
  par_set_base <- pars_df %>% slice(rep(1:n(), each = num_sims))
  
  #Replace lambda transmission parameter with sample from vals within 95%CI
  par_set_base["lambda"] <- sample(fit_1par_mat$lambda_1par, num_sims, 
                                   replace = T, prob = fit_1par_mat$weight)

  par_set_base["Study_ID_Unique"] <- "baseline_eec_parameters"

  #Save data frame of simulation parameter set
  saveRDS(par_set_base, file = paste0("Agrochemical_Review/Sims/EEC/Parameter_Sets/", "baseline_eec_parameters.rds"))

#Check that each unique response function has a unique parameter set ###########
#Right now, this is not true. Need to update naming convention
nrow(RFxSum_Sub) == length(list.files("Agrochemical_Review/Sims/EEC/Parameter_Sets/"))-1
