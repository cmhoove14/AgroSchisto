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
library(forestplot)

#load R0 function
source("Agrochemical_Review/Models/r0_of_q.R")

#GGplot theme for manuscripts
source("Agrochemical_Review/Sims/ggplot_theme.R")

#Response functions summary
rfx_sum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

rfx_dir <- "Agrochemical_Review/Response_Fxs/"

rfx_files <- list.files(path = "Agrochemical_Review/Response_Fxs",
                        pattern = "_fit.R", recursive = TRUE)

#All direct larvae studies
rfx_dirlarv <- rfx_sum %>% filter(Pathway == "direct larvae") %>% 
  rename(Parameter = parameter,
         Study = study_long,
         study_abrev = Study)

#Load all top-down response functions ######
dirlarv_studies <- unique(rfx_dirlarv$study_abrev)
dirlarv_files <- rfx_files[unlist(sapply(dirlarv_studies, grep, rfx_files, ignore.case = TRUE))]
dirlarv_files <- dirlarv_files[-unlist(sapply("snail", grep, dirlarv_files, ignore.case = TRUE))]

sapply(paste0(rfx_dir, dirlarv_files), source)

#Function to simulate 5000 parameter values, estimate r0 for each, return median and IQR for piC response functions
piC_r0 <- function(rfx, eec){
  print(rfx)
  piC_r0s <- sapply(replicate(5000, do.call(rfx, list(eec))), r0.fix, 
                     fNqx = 1, 
                     muPqx = 0,
                     fPqx = 1,
                     phiNqx = 1, 
                     muNqx = 0, 
                     psiqx = 1,
                     thetaqx = 1, 
                     piMqx = 1, 
                     #piCqx = 1, 
                     vqx = 1)[3,]
  
  return(data.frame(r0_med = median(piC_r0s),
                    r0_025 = quantile(piC_r0s, 0.25),
                    r0_975 = quantile(piC_r0s, 0.75)))
}

rfx_dirlarv_piC <- rfx_dirlarv %>% filter(Parameter == "piC" & !is.na(eec)) %>% 
  cbind(map2_df(.x = .$rfx, .y = .$eec, piC_r0))

#Function to simulate 5000 parameter values, estimate r0 for each, return median and IQR for piM response functions
piM_r0 <- function(rfx, eec){
  print(rfx)
  piM_r0s <- sapply(replicate(5000, do.call(rfx, list(eec))), r0.fix, 
                     fNqx = 1, 
                     muPqx = 0,
                     fPqx = 1,
                     phiNqx = 1, 
                     muNqx = 0, 
                     psiqx = 1,
                     thetaqx = 1, 
                     #piMqx = 1, 
                     piCqx = 1, 
                     vqx = 1)[3,]
  
  return(data.frame(r0_med = median(piM_r0s),
                    r0_025 = quantile(piM_r0s, 0.25),
                    r0_975 = quantile(piM_r0s, 0.75)))
}

rfx_dirlarv_piM <- rfx_dirlarv %>% filter(Parameter == "piM" & !is.na(eec)) %>% 
  cbind(map2_df(.x = .$rfx, .y = .$eec, piM_r0))

#Function to simulate 5000 parameter values, estimate r0 for each, return median and IQR for piM response functions
v_r0 <- function(rfx, eec){
  print(rfx)
  v_r0s <- sapply(replicate(5000, do.call(rfx, list(eec))), r0.fix, 
                     fNqx = 1, 
                     muPqx = 0,
                     fPqx = 1,
                     phiNqx = 1, 
                     muNqx = 0, 
                     psiqx = 1,
                     thetaqx = 1, 
                     piMqx = 1, 
                     #vqx = 1,
                     piCqx = 1)[3,]
  
  return(data.frame(r0_med = median(v_r0s),
                    r0_025 = quantile(v_r0s, 0.25),
                    r0_975 = quantile(v_r0s, 0.75)))
}

rfx_dirlarv_v <- rfx_dirlarv %>% filter(Parameter == "vq" & !is.na(eec)) %>% 
  cbind(map2_df(.x = .$rfx, .y = .$eec, v_r0))

#Put dfs back together
rfx_dirlarv_all <- rbind(rfx_dirlarv_piC, rfx_dirlarv_piM, rfx_dirlarv_v) %>% 
  mutate(r0_med_rel = (r0_med / r0.fix()[3]) * 100,
         r0_025_rel = (r0_025 / r0.fix()[3]) * 100,
         r0_975_rel = (r0_975 / r0.fix()[3]) * 100)

save(rfx_dirlarv_all, file = "Agrochemical_Review/Sims/EEC/all_dirlarv_summary.RData")

#Forestplotish thing with GGPlot
my_labs <- list(bquote(pi[C]), bquote(pi[M]), bquote(v))

rfx_dirlarv_all %>% 
  ggplot(aes(x = reorder(Chemical, desc(Chemical)), y = r0_med, shape = Parameter, col = Study)) + 
    geom_hline(yintercept = r0.fix()[3], lty = 2) +
    geom_point(size = 3, position = position_dodge(0.9)) + 
    geom_errorbar(aes(ymin = r0_025, ymax = r0_975, x = Chemical, width = 0.01), 
                  position = position_dodge(0.9)) +
    theme_ms() + ylim(0,4) + coord_flip() + ylab(expression(paste(R['0']))) +
    ggtitle("Direct larval effects") + xlab("Chemical") +
    scale_color_manual(values = glasbey()) + 
    scale_shape_manual(values = c(15,16,17), breaks = c("piC", "piM", "vq"),
                       labels = my_labs)

ggsave(paste('~/RemaisWork/Schisto/Agro_Review/Figures/EEC_forest/ggplot_forest_dirlarv', Sys.Date(), '.tiff', sep = ''),
       width = 7.3, height = 7.3, dpi = 600)
