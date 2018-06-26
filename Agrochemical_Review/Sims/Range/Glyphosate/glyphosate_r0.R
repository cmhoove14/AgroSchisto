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
load("Agrochemical_Review/Sims/Range/Glyphosate/gly_range.RData")

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

#Get phiNq function from baxter et al reanalysis
source("Agrochemical_Review/Response_Fxs/Baxter_Rohr2011_reanalysis_atrazine_snail_carrying_capacity_fit.R")

#Function to simulate n times parameter draws from each relevant response function, estimate net r0 incorporating all response fxs, and retunr median and IQR
nsims = 5000

require(parallel)

n_cores <- detectCores() - 1

clust <- makeCluster(n_cores)
clusterExport(clust, varlist = c(ls(), "uniroot.all"), envir = .GlobalEnv)

r0_gly <- function(gly){
  r0_glys <- parSapply(cl = clust, 
                       rep(gly, nsims), 
                       r0.He,
                       f.f_Nq = fNq.gly.fx.uncertainty, 
                       f.mu_Pq = nil0,
                       f.phi_Nq = phiNq_gly_baxrohr.no30,
                       f.mu_Nq = ons.muNq.gly,
                       f.alpha_q = nil1,
                       f.theta_q = nil1, 
                       f.pi_Mq = piM.ghaf_gly.exp_unc, 
                       f.pi_Cq = piC.ghaf_gly.exp_unc,
                       f.v_q = nil1)[3,]

  return(c(median(r0_glys),
           quantile(r0_glys, 0.25),
           quantile(r0_glys, 0.75)))
}

clusterSetRNGStream(clust, 43093)

gly_sims <- data.frame(t(sapply(gly_range, r0_gly, simplify = T)))

stopCluster(clust)

gly_sims <- cbind(gly_sims, gly = gly_range)
colnames(gly_sims)[c(1:3)] <- c("r0_med", "r0_25", "r0_75")

save(gly_sims, file = "Agrochemical_Review/Sims/Range/Glyphosate/gly_r0_sims.RData")
