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

#Get NAWQA data and function
load("~/RemaisWork/Schisto/R Codes/ag_schist/Agrochemical_Review/Sims/Data/NAWQA_dat_functions.RData")

#Get range of atrazine concentration to focus on based on NAWQA values
atr_vals <- as.numeric(get_nawqa_dat("Atrazine") %>% pull(CONCENTRATION))

atr_range <- seq(0, max(atr_vals), length.out = 30)

#GGplot theme for manuscripts
theme_ms <- function(base_size=12, base_family="Helvetica") {
  library(grid)
  (theme_bw(base_size = base_size, base_family = base_family)+
      theme(text=element_text(color="black"),
            axis.title=element_text(face="bold", size = rel(1.3)),
            axis.text=element_text(size = rel(1), color = "black"),
            legend.title=element_text(face="bold"),
            legend.text=element_text(face="bold"),
            legend.background=element_rect(fill="transparent"),
            legend.key.size = unit(0.8, 'lines'),
            panel.border=element_rect(color="black",size=1),
            panel.grid=element_blank()
    ))
}

#Response functions summary
rfx_sum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")

rfx_dir <- "Agrochemical_Review/Response_Fxs/"

rfx_files <- list.files(path = "Agrochemical_Review/Response_Fxs",
                        pattern = "_fit.R", recursive = TRUE)

#Atrazine response functions for S. mansoni
rfx_atrazine <- rfx_sum %>% filter(Chemical == "Atrazine" & 
                                     best == 1 & 
                                     System %in% c("Mansoni", "Any")) %>% 
  rename(Parameter = parameter,
         Study = study_long,
         study_abrev = Study)

#Load relevant atrazine response functions ######
  atrazine_studies <- unique(rfx_atrazine$study_abrev)
  atrazine_functions <- rfx_atrazine$rfx
  atrazine_files <- rfx_files[unlist(sapply(atrazine_studies, grep, rfx_files, ignore.case = TRUE))]

sapply(paste0(rfx_dir, atrazine_files), source)
#manually get right cercarial die off function
source("Agrochemical_Review/Response_Fxs/rohr2008_atrazine_cercariae_fit.R")

#Function to simulate 1000 parameter draws from each relevant response function, estimate net r0 incorporating all response fxs, and retunr median and IQR
atr_list <- list()

r0_atr <- function(atr){
  r0_atrs <- sapply(rep(atr, 1000), r0.He,
                    f.f_Nq = nil1, 
                    f.mu_Pq = muPq_halstead18_atr102_uncertainty,
                    f.phi_Nq = phiNq_atr_baxrohr.no30,
                    f.mu_Nq = ons.muNq.atr,
                    f.alpha_q = nil1,
                    f.theta_q = nil1, 
                    f.pi_Mq = nil1, 
                    f.pi_Cq = piC.atr.rohr08.lin,
                    f.v_q = halstead_meso18_atr_mans_v_uncertainty)[3,]

  return(data.frame(r0_med = median(r0_atrs),
                    r0_25 = quantile(r0_atrs, 0.25),
                    r0_75 = quantile(r0_atrs, 0.75)))
}

set.seed(43093)
atr_sims <- data.frame(t(sapply(atr_range, r0_atr, simplify = T)))

atr_sims <- cbind(atr_sims, atr = atr_range)

atr_sims %>% ggplot(aes(x = atr, y = unlist(r0_med))) + stat_smooth()
