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

#Load butralin concentration values from NAWQA and sample range for r0 simulation
load("Agrochemical_Review/Sims/Range/Butralin/btr_range.RData")

#Response functions summary
rfx_sum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

rfx_dir <- "Agrochemical_Review/Response_Fxs/"

rfx_files <- list.files(path = "Agrochemical_Review/Response_Fxs",
                        pattern = "_fit.R", recursive = TRUE)

#Butralin response functions for S. mansoni
rfx_butralin <- rfx_sum %>% filter(Chemical == "Butralin") %>% 
  rename(Parameter = parameter,
         Study = study_long,
         study_abrev = Study)

#Load relevant butralin response functions ######
  butralin_studies <- unique(rfx_butralin$study_abrev)
  butralin_functions <- rfx_butralin$rfx
  butralin_files <- rfx_files[unlist(sapply(butralin_studies, grep, rfx_files, ignore.case = TRUE))]

sapply(paste0(rfx_dir, butralin_files), source)

#Get relative bottom up function
source("Agrochemical_Review/Response_Fxs/Baxter_Rohr2011_reanalysis_atrazine_snail_carrying_capacity_fit.R")

#Function to simulate n times parameter draws from each relevant response function, estimate net r0 incorporating all response fxs, and retunr median and IQR
nsims = 5000

require(parallel)

n_cores <- detectCores() - 1

clust <- makeCluster(n_cores)
clusterExport(clust, varlist = c(ls(), "uniroot.all"), envir = .GlobalEnv)

r0_btr <- function(btr){
  r0_btrs <- parSapply(cl = clust, 
                       rep(btr, nsims), 
                       r0.He,
                       f.f_Nq = fNq.butr.fx.uncertainty, 
                       f.mu_Pq = nil0,
                       f.phi_Nq = phiNq_btr_baxrohr.no30,
                       f.mu_Nq = muNq_butr_gaf16_uncertainty,
                       f.alpha_q = nil1,
                       f.theta_q = nil1, 
                       f.pi_Mq = piM.ghaf_butr.exp_unc, 
                       f.pi_Cq = piC.ghaf_butr.exp_unc,
                       f.v_q = nil1)[3,]

  return(c(median(r0_btrs),
           quantile(r0_btrs, 0.25),
           quantile(r0_btrs, 0.75)))
}

clusterSetRNGStream(clust, 43093)

btr_sims <- data.frame(t(sapply(btr_range, r0_btr, simplify = T)))

stopCluster(clust)

btr_sims <- cbind(btr_sims, btr = btr_range)
colnames(btr_sims)[c(1:3)] <- c("r0_med", "r0_25", "r0_75")

save(btr_sims, file = "Agrochemical_Review/Sims/Range/Butralin/btr_r0_sims.RData")
