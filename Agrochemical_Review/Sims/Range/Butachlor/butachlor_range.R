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

#Get butachlor EEC
rfx_sum <- read_csv("Agrochemical_Review/Response_Fxs/Summary/Response_Fx_Summary.csv")
  but_eec <- unique(rfx_sum %>% filter(Chemical == "Butachlor") %>% pull(eec))

#Get range of butralin concentration to focus on based on NAWQA values
but_vals <- as.numeric(get_nawqa_dat("Butachlor") %>% 
                         filter(REMARK != "<") %>% 
                         pull(CONCENTRATION))

but_range <- get_range(but_vals, but_eec)

save(but_vals, but_range, but_eec, file = "Agrochemical_Review/Sims/Range/Butachlor/but_range.RData")
