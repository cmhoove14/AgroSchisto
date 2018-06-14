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
library(pals)

#load R0 function
source("Agrochemical_Review/Models/r0_of_q.R")

#Response functions summary
rfx_sum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

rfx_dir <- "Agrochemical_Review/Response_Fxs/"

rfx_files <- list.files(path = "Agrochemical_Review/Response_Fxs",
                        pattern = "_fit.R", recursive = TRUE)

#Initially was trying to work on full data frame, but introduces issues since there are so many objects required, therefore separating jobs by pathway. Also makse sense since we're constructing separate forestplots for each pathway
rfx_dirsnail <- rfx_sum %>% filter(Pathway == "direct snail") %>% 
  rename(Parameter = parameter,
         Study = study_long,
         study_abrev = Study)

#Load all top-down response functions ######
dirsnail_studies <- unique(rfx_dirsnail$study_abrev)

sapply(paste0(rfx_dir, rfx_files[unlist(sapply(dirsnail_studies, grep, rfx_files, ignore.case = TRUE))]), source)

#Function to simulate 5000 parameter values, estimate r0 for each, return median and IQR for muNq response functions
muNq_r0 <- function(rfx, eec){
  print(rfx)
  muNq_r0s <- sapply(replicate(5000, do.call(rfx, list(eec))), r0.fix, 
                     fNqx = 1, 
                     muPqx = 0,
                     fPqx = 1,
                     phiNqx = 1, 
                     #muNqx = 0, 
                     psiqx = 1,
                     thetaqx = 1, 
                     piMqx = 1, 
                     piCqx = 1, 
                     vqx = 1)[3,]
  
  return(data.frame(r0_med = median(muNq_r0s),
                    r0_025 = quantile(muNq_r0s, 0.25),
                    r0_975 = quantile(muNq_r0s, 0.75)))
}

rfx_dirsnail_muNq <- rfx_dirsnail %>% filter(Parameter == "muNq" & !is.na(eec)) %>% 
  cbind(map2_df(.x = .$rfx, .y = .$eec, muNq_r0))

#Function to simulate 5000 parameter values, estimate r0 for each, return median and IQR for fNq response functions
fNq_r0 <- function(rfx, eec){
  print(rfx)
  fNq_r0s <- sapply(replicate(5000, do.call(rfx, list(eec))), r0.fix, 
                     #fNqx = 1, 
                     muPqx = 0,
                     fPqx = 1,
                     phiNqx = 1, 
                     muNqx = 0, 
                     psiqx = 1,
                     thetaqx = 1, 
                     piMqx = 1, 
                     piCqx = 1, 
                     vqx = 1)[3,]
  
  return(data.frame(r0_med = median(fNq_r0s),
                    r0_025 = quantile(fNq_r0s, 0.25),
                    r0_975 = quantile(fNq_r0s, 0.75)))
}

rfx_dirsnail_fNq <- rfx_dirsnail %>% filter(Parameter == "fNq" & !is.na(eec)) %>% 
  cbind(map2_df(.x = .$rfx, .y = .$eec, fNq_r0))

#Put dfs back together
rfx_dirsnail_all <- rbind(rfx_dirsnail_muNq, rfx_dirsnail_fNq) %>% 
  mutate(r0_med_rel = (r0_med / r0.fix()[3]) * 100,
         r0_025_rel = (r0_025 / r0.fix()[3]) * 100,
         r0_975_rel = (r0_975 / r0.fix()[3]) * 100)

#Forestplotish thing with GGPlot
my_labs <- list(bquote(mu[N]), bquote(f[N]))

tiff(paste('~/RemaisWork/Schisto/Agro_Review/Figures/EEC_forest/ggplot_forest_dirsnail', Sys.Date(), '.tiff', sep = ''),
     width = 2480, height = 3508*0.5, res = 300)
rfx_dirsnail_all %>% #mutate(axis_lab = paste(study_long, Species, sep = "  ")) %>% 
  ggplot(aes(x = Chemical, y = r0_med, shape = Parameter, col = Study)) + 
    geom_hline(yintercept = r0.fix()[3], lty = 2) +
    geom_point(size = 2, position = position_dodge(0.5)) + 
    geom_errorbar(aes(ymin = r0_025, ymax = r0_975, x = Chemical, width = 0.01), 
                  position = position_dodge(0.5)) +
    theme_bw() + ylim(0,4) + coord_flip() + ylab(expression(paste(R['0']))) +
    ggtitle("Direct snail effects") +
    scale_color_manual(values = glasbey()) + 
    scale_shape_manual(values = c(16,15), breaks = c("muNq", "fNq"),
                       labels = my_labs)
dev.off()
