#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Estimate 24-hr P. alleni mortality from mesocosm data
meso_dat_preds<-read_csv('Agrochemical_Review/Response_Fxs/Data/Halstead_meso_data.csv') %>%
  mutate(atr_conc = 102 * A2,
         chlor_conc = 64 * C2,
         fert_conc = 4400 * F2,
         p.all_tot = 3)

#Chlorpyrifos
chlor_muPq_rate <- meso_dat_preds$P.all.24[which(meso_dat_preds$Treat1 == "C 1x")] / 
  meso_dat_preds$p.all_tot[which(meso_dat_preds$Treat1 == "C 1x")]
  
muPq_halstead18_chlor64_uncertainty <- function(...){
  sample(chlor_muPq_rate, 1)
}

#Atrazine
atr_muPq_rate <- meso_dat_preds$P.all.24[which(meso_dat_preds$Treat1 == "A 1x")] / 
  meso_dat_preds$p.all_tot[which(meso_dat_preds$Treat1 == "A 1x")]

muPq_halstead18_atr102_uncertainty <- function(...){
  sample(atr_muPq_rate, 1)
}

#Fertilizer
fert_muPq_rate <- meso_dat_preds$P.all.24[which(meso_dat_preds$Treat1 == "F 1x")] / 
  meso_dat_preds$p.all_tot[which(meso_dat_preds$Treat1 == "F 1x")]

muPq_halstead18_fert4400_uncertainty <- function(...){
  sample(fert_muPq_rate, 1)
}