#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
require(tidyverse)

#load R0 function
source("Agrochemical_Review/Models/r0_of_q.R")

#Load profenofos concentration values from NAWQA and sample range for r0 simulation
load("Agrochemical_Review/Sims/Range/Profenofos/prof_range.RData")

#Response functions summary
rfx_sum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

rfx_dir <- "Agrochemical_Review/Response_Fxs/"

rfx_files <- list.files(path = "Agrochemical_Review/Response_Fxs",
                        pattern = "_fit.R", recursive = TRUE)

#Profenofos response functions for S. mansoni
rfx_profenofos <- rfx_sum %>% filter(Chemical == "Profenofos" & 
                                     best == 1) %>% 
  rename(Parameter = parameter,
         Study = study_long,
         study_abrev = Study)

#Load relevant Profenofos response functions ######
  profenofos_studies <- unique(rfx_profenofos$study_abrev)
  profenofos_functions <- rfx_profenofos$rfx
  profenofos_files <- rfx_files[unlist(sapply(profenofos_studies, grep, rfx_files, ignore.case = TRUE))]

sapply(paste0(rfx_dir, profenofos_files), source)

#Function to simulate n times parameter draws from each relevant response function, estimate net r0 incorporating all response fxs, and retunr median and IQR
nsims = 5000

require(parallel)

n_cores <- detectCores() - 1

clust <- makeCluster(n_cores)
clusterExport(clust, varlist = c(ls(), "uniroot.all"), envir = .GlobalEnv)

r0_prof <- function(prof){
  r0_profs <- parSapply(cl = clust, 
                       rep(prof, nsims), 
                       r0.In,
                       f.f_Nq = fNq_prof_tch91_uncertainty, 
                       f.mu_Pq = muPq_prof_mac_rohr_unpub_uncertainty,
                       f.phi_Nq = nil1,
                       f.mu_Nq = muNq_prof_tch91_uncertainty,
                       f.alpha_q = nil1,
                       f.theta_q = nil1, 
                       f.pi_Mq = piM.tch91_prof_unc, 
                       f.pi_Cq = piC.tch92_prof_unc,
                       f.v_q = nil1)[3,]

  return(c(median(r0_profs),
           quantile(r0_profs, 0.25),
           quantile(r0_profs, 0.75)))
}

clusterSetRNGStream(clust, 43093)

prof_sims <- data.frame(t(sapply(prof_range, r0_prof, simplify = T)))

stopCluster(clust)

prof_sims <- cbind(prof_sims, prof = prof_range)
colnames(prof_sims)[c(1:3)] <- c("r0_med", "r0_25", "r0_75")

save(prof_sims, file = "Agrochemical_Review/Sims/Range/Profenofos/prof_r0_sims.RData")
