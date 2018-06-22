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

#Estimate 1-hr Schistosoma egg viability
meso_dat_eggs <- read_csv('Agrochemical_Review/Response_Fxs/Data/Halstead_meso_egg_viability.csv') %>% 
  mutate(atr_conc = 102*Atrazine,
         chlor_conc = 64*Chlorpyrifos,
         fert_conc = 4400*Fertilizer)

#Check for block and replicate effects 
eggs_check <- glm(miracidia ~ as.factor(Plate) + as.factor(Block) + atr_conc + chlor_conc + fert_conc + offset(log(eggs)), 
                  family = "poisson", data = meso_dat_eggs)

# Looks like only changes in egg viability are due to block/replicate effects, which agrees with analysis in supp info for Nature comms paper

#No evidence of excess mortality due to chlorpyrifos at EEC
vq_halstead18_chlor64_uncertainty <- function(...){
  return(0)
}

#No evidence of excess mortality due to atrazine at EEC
vq_halstead18_atr102_uncertainty <- function(...){
  return(0)
}

#No evidence of excess mortality due to fertilizer at EEC
vq_halstead18_fert4400_uncertainty <- function(...){
  return(0)
}