require(tidyverse)

load("Agrochemical_Review/Models/trans_pars_fit.Rdata")
load("Agrochemical_Review/Models/fit_pars.Rdata")

source("Agrochemical_Review/Models/r0_functions.R")
source("Agrochemical_Review/Models/model_helper_functions.R")
source("Agrochemical_Review/Models/initial_parameters.R")

pars_df <- as.data.frame(as.list(fit_pars))

#Get csv file containing summary of all response functions
RFxSum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

#Source all .R scripts that contain derivation of response functions
RFx_files <- paste0("Agrochemical_Review/Response_Fxs/", 
                    list.files("Agrochemical_Review/Response_Fxs/")[grepl(".R", list.files("Agrochemical_Review/Response_Fxs/"))])

sapply(RFx_files, source)

#Function to generate parameter sets to simulate agrochemical effects at EEC for each response function, then apply the r0 function to every realization of the parameter set and save the resulting dataframe
gen_eec_par_set <- function(eec, nsims, rfx, par, out_name){
  
  #create nsims duplicate rows of base parameters
  par_set <- pars_df %>% slice(rep(1:n(), each = nsims))
  
  #Replace lambda transmission parameter with sample from vals within 95%CI
  par_set["lambda"] <- sample(fit_1par_mat$lambda_1par, nsims, 
                              replace = T, prob = fit_1par_mat$weight)
  
  #Replace affected agrochemical parameter with estimates from rfx
  par_set[par] <- sapply(rep(eec, nsims), rfx)
  
  #Apply r0 function to parameter set
  r0s <- pmap_df(par_set, r0.Ag.pars)
  
  #bind parameters and rsults together
  r0_out <- cbind(par_set, r0s)
  
  r0_out["Study"] <- out_name

  #Save data frame of simulation parameter set
  save(r0_out, file = paste0("Agrochemical_Review/Sims/EEC/Parameter_Sets/", out_name, "_eec_parameters_r0s.Rdata"))
  
}

#Do this for all response functions with relevant parameters affected
RFxSum_Sub <- RFxSum %>% filter(parameter %in% colnames(pars_df) & !is.na(eec) & Class != "Fertilizer")
 
set.seed(43093)

for(i in 1:nrow(RFxSum_Sub)){
  gen_eec_par_set(RFxSum_Sub$eec[i], 1000, 
                  RFxSum_Sub$rfx[i], RFxSum_Sub$parameter[i], 
                  paste(RFxSum_Sub$Study[i], RFxSum_Sub$Chemical[i], sep = "_"))
  
  if(i %% 20 == 0) print(i)
} 

#Check that each unique response function has a unique parameter set
#Right now, this is not true. Need to update naming convention
nrow(RFxSum_Sub) == length(list.files("Agrochemical_Review/Sims/EEC/Parameter_Sets/"))
