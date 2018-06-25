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

#Get range of Glyphosate concentration to focus on based on NAWQA values
gly_vals <- as.numeric(get_nawqa_dat("glyphosate") %>% 
                           filter(GLYPHOSATE_RMK_ELISA != "<") %>% 
                           pull(GLYPHOSATE_ELISA_ug_L))

gly_range <- c(0, exp(seq(log(min(gly_vals)), log(max(gly_vals)), length.out = 100)), seq(28, 40, by = 2))

save(gly_vals, gly_range, file = "Agrochemical_Review/Sims/Range/Glyphosate/gly_range.RData")
