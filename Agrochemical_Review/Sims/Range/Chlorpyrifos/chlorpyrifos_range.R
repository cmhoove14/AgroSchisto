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

#Get NAWQA data and function
load("~/RemaisWork/Schisto/R Codes/ag_schist/Agrochemical_Review/Sims/Data/NAWQA_dat_functions.RData")

#Get chlorpyrifos EEC
rfx_sum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")
  chlor_eec <- unique(rfx_sum %>% filter(Chemical == "Chlorpyrifos") %>% pull(eec))

#Get range of Chlorpyrifos concentration to focus on based on NAWQA values
chlor_vals <- as.numeric(get_nawqa_dat("Chlorpyrifos") %>% 
                           filter(SITE_TYPE == "Agriculture" & REMARK != "<") %>% 
                           pull(CONCENTRATION))

chlor_range <- get_range(chlor_vals, chlor_eec)

save(chlor_vals, chlor_range, chlor_eec, file = "Agrochemical_Review/Sims/Range/Chlorpyrifos/chlor_range.RData")
