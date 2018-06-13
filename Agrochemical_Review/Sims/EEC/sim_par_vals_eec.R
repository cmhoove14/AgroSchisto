#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

library(tidyverse)

#load R0 function
source("Agrochemical_Review/Models/r0_of_q.R")

#Response functions summary
rfx_sum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

rfx_dir <- "Agrochemical_Review/Response_Fxs/"

rfx_files <- list.files(path = "Agrochemical_Review/Response_Fxs",
                        pattern = "_fit.R", recursive = TRUE)

#Initially was trying to work on full data frame, but introduces issues since there are so many objects required, therefore separating jobs by pathway. Also makse sense since we're constructing separate forestplots for each pathway
rfx_topdown <- rfx_sum %>% filter(Pathway == "top-down")

#Load all top-down response functions ######
topdown_studies <- unique(rfx_topdown$Study)

for(i in 1:length(topdown_studies)){
  source(paste0(rfx_dir, rfx_files[grep(topdown_studies[i], rfx_files, ignore.case = TRUE)]))
  print(i)
}

#Function to simulate 5000 parameter values, estimate r0 for each, return median and IQR for muPq response functions
muPq_r0 <- function(rfx, eec){
  print(rfx)
  muPq_r0s <- sapply(replicate(5000, do.call(rfx, list(eec))), r0.fix, 
                     fNqx = 1, 
                     #muPqx = 0,
                     fPqx = 1,
                     phiNqx = 1, 
                     muNqx = 0, 
                     psiqx = 1,
                     thetaqx = 1, 
                     piMqx = 1, 
                     piCqx = 1, 
                     vqx = 1)[3,]
  
  return(data.frame(r0_med = median(muPq_r0s),
                    r0_025 = quantile(muPq_r0s, 0.25),
                    r0_975 = quantile(muPq_r0s, 0.75)))
}

rfx_topdown_muPq <- rfx_topdown %>% filter(parameter == "muPq" & !is.na(eec)) %>% 
  cbind(map2_df(.x = .$rfx, .y = .$eec, muPq_r0))

#Function to simulate 5000 parameter values, estimate r0 for each, return median and IQR for psiq response functions
psiq_r0 <- function(rfx, eec){
  print(rfx)
  psiq_r0s <- sapply(replicate(5000, do.call(rfx, list(eec))), r0.fix, 
                     fNqx = 1, 
                     muPqx = 0,
                     fPqx = 1,
                     phiNqx = 1, 
                     muNqx = 0, 
                     #psiqx = 1,
                     thetaqx = 1, 
                     piMqx = 1, 
                     piCqx = 1, 
                     vqx = 1)[3,]
  
  return(data.frame(r0_med = median(psiq_r0s),
                    r0_025 = quantile(psiq_r0s, 0.25),
                    r0_975 = quantile(psiq_r0s, 0.75)))
}

rfx_topdown_psiq <- rfx_topdown %>% filter(parameter == "psiq" & !is.na(eec)) %>% 
  cbind(map2_df(.x = .$rfx, .y = .$eec, psiq_r0))

#Function to simulate 5000 parameter values, estimate r0 for each, return median and IQR for fPq response functions
fPq_r0 <- function(rfx, eec){
  print(rfx)
  fPq_r0s <- sapply(replicate(5000, do.call(rfx, list(eec))), r0.fix, 
                     fNqx = 1, 
                     muPqx = 0,
                     #fPqx = 1,
                     phiNqx = 1, 
                     muNqx = 0, 
                     psiqx = 1,
                     thetaqx = 1, 
                     piMqx = 1, 
                     piCqx = 1, 
                     vqx = 1)[3,]
  
  return(data.frame(r0_med = median(fPq_r0s),
                    r0_025 = quantile(fPq_r0s, 0.25),
                    r0_975 = quantile(fPq_r0s, 0.75)))
}

rfx_topdown_fPq <- rfx_topdown %>% filter(parameter == "fPq" & !is.na(eec)) %>% 
  cbind(map2_df(.x = .$rfx, .y = .$eec, psiq_r0))

#Put dfs back together
rfx_topdown_all <- rbind(rfx_topdown_muPq, rfx_topdown_psiq, rfx_topdown_fPq) %>% 
  mutate(r0_med_rel = (r0_med / r0.fix()[3]) * 100,
         r0_025_rel = (r0_025 / r0.fix()[3]) * 100,
         r0_975_rel = (r0_975 / r0.fix()[3]) * 100)

rfx_topdown_all %>% 
  ggplot(aes(x = Chemical, y = r0_med, shape = parameter, color = study_long)) + 
  geom_hline(yintercept = r0.fix()[3], lty = 2) +
  geom_point(size = 3, position = position_dodge(1)) + 
  geom_errorbar(aes(ymin = r0_025, ymax = r0_975, x = Chemical, width = 0.01), position = position_dodge(1)) +
  theme_bw() + ylim(0,4) + coord_flip()

