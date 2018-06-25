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

#Load glyphosate concentration values from NAWQA and sample range for r0 simulation
load("Agrochemical_Review/Sims/Range/Glyphosate/chlor_range.RData")

#Response functions summary
rfx_sum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

rfx_dir <- "Agrochemical_Review/Response_Fxs/"

rfx_files <- list.files(path = "Agrochemical_Review/Response_Fxs",
                        pattern = "_fit.R", recursive = TRUE)

#Glyphosate response functions for S. mansoni
rfx_glyphosate <- rfx_sum %>% filter(Chemical == "Glyphosate" & 
                                     best == 1) %>% 
  rename(Parameter = parameter,
         Study = study_long,
         study_abrev = Study)

#Load relevant Glyphosate response functions ######
  glyphosate_studies <- unique(rfx_glyphosate$study_abrev)
  glyphosate_functions <- rfx_glyphosate$rfx
  glyphosate_files <- rfx_files[unlist(sapply(glyphosate_studies, grep, rfx_files, ignore.case = TRUE))]

sapply(paste0(rfx_dir, glyphosate_files), source)

#Function to simulate n times parameter draws from each relevant response function, estimate net r0 incorporating all response fxs, and retunr median and IQR
nsims = 5000

require(parallel)

n_cores <- detectCores() - 1

clust <- makeCluster(n_cores)
clusterExport(clust, varlist = c(ls(), "uniroot.all"), envir = .GlobalEnv)

r0_chlor <- function(chlor){
  r0_chlors <- parSapply(cl = clust, 
                       rep(chlor, nsims), 
                       r0.In,
                       f.f_Nq = nil1, 
                       f.mu_Pq = muPq_chlor_mac_rohr_unpub_uncertainty,
                       f.phi_Nq = nil1,
                       f.mu_Nq = muNq_chlor_ibr92_uncertainty,
                       f.alpha_q = nil1,
                       f.theta_q = nil1, 
                       f.pi_Mq = piM_ch_Hash11_uncertainty, 
                       f.pi_Cq = piC_ch_Hash11_uncertainty,
                       f.v_q = halstead_meso18_chlor_mans_v_uncertainty)[3,]

  return(c(median(r0_chlors),
           quantile(r0_chlors, 0.25),
           quantile(r0_chlors, 0.75)))
}

clusterSetRNGStream(clust, 43093)

chlor_sims <- data.frame(t(sapply(chlor_range, r0_chlor, simplify = T)))

stopCluster(clust)

chlor_sims <- cbind(chlor_sims, chlor = chlor_range)
colnames(chlor_sims)[c(1:3)] <- c("r0_med", "r0_25", "r0_75")

save(chlor_sims, file = "Agrochemical_Review/Sims/Range/glyphosate/chlor_r0_sims.RData")
