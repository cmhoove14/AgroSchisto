#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
source("Agrochemical_Review/Response_Fxs/Halstead_meso2018_eggs_fit.R")

#Atrazine mansoni
meso_dat_mans_eggs %>% filter(Treatment %in% c("S", "A1x", "A2x", "H2O", "W")) %>% 
  ggplot(aes(x = atr_conc_logp1, y = probit_per_hatch)) + geom_point() + theme_bw() +
  geom_smooth(method = "lm") + ylab("probit hatching percent") + xlab("log+1 Atrazine (ppb)") +
  geom_point(data = data.frame(atr_conc_logp1 = log(c(0:102)+1),
                        prediction = qnorm(sapply(c(0:102), halstead_meso18_atr_mans_v_uncertainty)*meso_mans_egg_atr_int)), 
             aes(x = atr_conc_logp1, y = prediction), pch = 5, col = 4, cex = 0.7)

ggsave("Agrochemical_Review/Response_Fxs/Plots/Halstead_meso2018/mansoni_egg_v_atrazine.png", device = "png")

#Chlorpyrifos mansoni
meso_dat_mans_eggs %>% filter(Treatment %in% c("S", "C1x", "C1X", "C2x", "H2O", "W")) %>% 
  ggplot(aes(x = chlor_conc_logp1, y = probit_per_hatch)) + geom_point() + theme_bw() +
  geom_smooth(method = "lm") + ylab("probit hatching percent") + xlab("log+1 Chlorpyrifos (ppb)") +
  geom_point(data = data.frame(chlor_conc_logp1 = log(c(0:64)+1),
                        prediction = qnorm(sapply(c(0:64), halstead_meso18_chlor_mans_v_uncertainty)*meso_mans_egg_chlor_int)), 
             aes(x = chlor_conc_logp1, y = prediction), pch = 5, col = 4, cex = 0.7)

ggsave("Agrochemical_Review/Response_Fxs/Plots/Halstead_meso2018/mansoni_egg_v_chlorpyrifos.png", device = "png")

#Fertilizer mansoni
meso_dat_mans_eggs %>% filter(Treatment %in% c("S", "F1x", "F1X", "F2x", "H2O", "W")) %>% 
  ggplot(aes(x = fert_conc_logp1, y = probit_per_hatch)) + geom_point() + theme_bw() +
  geom_smooth(method = "lm") + ylab("probit hatching percent") + xlab("log+1 Fertilizer (N-ppb)") +
  geom_point(data = data.frame(fert_conc_logp1 = log(seq(0,4400,40)+1),
                        prediction = qnorm(sapply(seq(0,4400,40), halstead_meso18_fert_mans_v_uncertainty)*meso_mans_egg_fert_int)), 
             aes(x = fert_conc_logp1, y = prediction), pch = 5, col = 4, cex = 0.7)

ggsave("Agrochemical_Review/Response_Fxs/Plots/Halstead_meso2018/mansoni_egg_v_mansoni.png", device = "png")

#haematobium Atrazine
meso_dat_haem_eggs %>% filter(Treatment %in% c("S", "A1x", "A2x", "H2O", "W")) %>% 
  ggplot(aes(x = atr_conc_logp1, y = probit_per_hatch)) + geom_point() + theme_bw() +
  geom_smooth(method = "lm") + ylab("probit hatching percent") + xlab("log+1 Atrazine (ppb)") +
  geom_point(data = data.frame(atr_conc_logp1 = log(c(0:102)+1),
                        prediction = qnorm(sapply(c(0:102), halstead_meso18_atr_haem_v_uncertainty)*meso_haem_egg_atr_int)), 
             aes(x = atr_conc_logp1, y = prediction), pch = 5, col = 4, cex = 0.7)

ggsave("Agrochemical_Review/Response_Fxs/Plots/Halstead_meso2018/haematobium_egg_v_atrazine.png", device = "png")

# haematobium Chlorpyrifos
meso_dat_haem_eggs %>% filter(Treatment %in% c("S", "C1x", "C1X", "C2x", "H2O", "W")) %>% 
  ggplot(aes(x = chlor_conc_logp1, y = probit_per_hatch)) + geom_point() + theme_bw() +
  geom_smooth(method = "lm") + ylab("probit hatching percent") + xlab("log+1 Chlorpyrifos (ppb)") +
  geom_point(data = data.frame(chlor_conc_logp1 = log(c(0:64)+1),
                        prediction = qnorm(sapply(c(0:64), halstead_meso18_chlor_haem_v_uncertainty)*meso_haem_egg_chlor_int)), 
             aes(x = chlor_conc_logp1, y = prediction), pch = 5, col = 4, cex = 0.7)

ggsave("Agrochemical_Review/Response_Fxs/Plots/Halstead_meso2018/haematobium_egg_v_chlorpyrifos.png", device = "png")

# Hameatobium Fertilizer
meso_dat_haem_eggs %>% filter(Treatment %in% c("S", "F1x", "F1X", "F2x", "H2O", "W")) %>% 
  ggplot(aes(x = fert_conc_logp1, y = probit_per_hatch)) + geom_point() + theme_bw() +
  geom_smooth(method = "lm") + ylab("probit hatching percent") + xlab("log+1 Fertilizer (N-ppb)") +
  geom_point(data = data.frame(fert_conc_logp1 = log(seq(0,4400,40)+1),
                        prediction = qnorm(sapply(seq(0,4400,40), halstead_meso18_fert_haem_v_uncertainty)*meso_haem_egg_fert_int)), 
             aes(x = fert_conc_logp1, y = prediction), pch = 5, col = 4, cex = 0.7)

ggsave("Agrochemical_Review/Response_Fxs/Plots/Halstead_meso2018/haematobium_egg_v_fertilizer.png", device = "png")
