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
  mutate(probit_per_hatch = qnorm(per_hatch),
         atr_conc = 102*Atrazine,
         atr_conc_logp1 = log(atr_conc+1),
         chlor_conc = 64*Chlorpyrifos,
         chlor_conc_logp1 = log(chlor_conc+1),
         fert_conc = 4400*Fertilizer,
         fert_conc_logp1 = log(fert_conc+1))

#Mansoni data***********################################
meso_dat_mans_eggs <- meso_dat_eggs %>% filter(species == "mansoni")

#Atrazine ##########
#model
meso_mans_egg_atr <- lm(probit_per_hatch ~ atr_conc_logp1, weights = eggs,
                        data = meso_dat_mans_eggs %>% 
                          filter(Treatment %in% c("S", "A1x", "A2x", "H2O", "W") & 
                                   is.finite(probit_per_hatch)))

meso_mans_egg_atr_int <- pnorm(coef(meso_mans_egg_atr)[1])
  
#Atrazine function returns relative proportion change in viability from baseline
halstead_meso18_atr_mans_v_uncertainty<-function(He){
  init = predict(meso_mans_egg_atr, newdata = data.frame(atr_conc_logp1 = log(He+1)), se.fit = T)
  
  pred = pnorm(rnorm(1, init$`fit`, init$`se.fit`))
  
  pred / meso_mans_egg_atr_int
}

#plot viability over atrazine concentration with predictions
meso_dat_mans_eggs %>% filter(Treatment %in% c("S", "A1x", "A2x", "H2O", "W")) %>% 
  ggplot(aes(x = atr_conc_logp1, y = probit_per_hatch)) + geom_point() + theme_bw() +
  geom_smooth(method = "lm") + ylab("probit hatching percent") + xlab("log+1 Atrazine (ppb)") +
  geom_point(data = data.frame(atr_conc_logp1 = log(c(0:102)+1),
                        prediction = qnorm(sapply(c(0:102), halstead_meso18_atr_mans_v_uncertainty)*meso_mans_egg_atr_int)), 
             aes(x = atr_conc_logp1, y = prediction), pch = 5, col = 4, cex = 0.7)

#Chlorpyrifos ##########
#Chlorpyrifos model
meso_mans_egg_chlor<-lm(probit_per_hatch ~ chlor_conc_logp1, weights = eggs,
                        data = meso_dat_mans_eggs %>% 
                          filter(Treatment %in% c("S", "C1x", "C1X", "C2x", "H2O", "W") & 
                                   is.finite(probit_per_hatch)))

meso_mans_egg_chlor_int <- pnorm(coef(meso_mans_egg_chlor)[1])

#Chlorpyrifos function
halstead_meso18_chlor_mans_v_uncertainty<-function(In){
  init = predict(meso_mans_egg_chlor, newdata = data.frame(chlor_conc_logp1 = log(In+1)), se.fit = T)
  
  pred = pnorm(rnorm(1, init$`fit`, init$`se.fit`))
  
  pred / meso_mans_egg_chlor_int
}

#plot viability over chlorpyrifos concentration
meso_dat_mans_eggs %>% filter(Treatment %in% c("S", "C1x", "C1X", "C2x", "H2O", "W")) %>% 
  ggplot(aes(x = chlor_conc_logp1, y = probit_per_hatch)) + geom_point() + theme_bw() +
  geom_smooth(method = "lm") + ylab("probit hatching percent") + xlab("log+1 Chlorpyrifos (ppb)") +
  geom_point(data = data.frame(chlor_conc_logp1 = log(c(0:64)+1),
                        prediction = qnorm(sapply(c(0:64), halstead_meso18_chlor_mans_v_uncertainty)*meso_mans_egg_chlor_int)), 
             aes(x = chlor_conc_logp1, y = prediction), pch = 5, col = 4, cex = 0.7)

#Fertilizer ##########
#Fertilizer model
meso_mans_egg_fert<-lm(probit_per_hatch ~ fert_conc_logp1, weights = eggs,
                        data = meso_dat_mans_eggs %>% 
                          filter(Treatment %in% c("S", "F1x", "F1X", "F2x", "H2O", "W") & 
                                   is.finite(probit_per_hatch)))

meso_mans_egg_fert_int <- pnorm(coef(meso_mans_egg_fert)[1])

#Fertilizer function
halstead_meso18_fert_mans_v_uncertainty<-function(Fe){
  init = predict(meso_mans_egg_fert, newdata = data.frame(fert_conc_logp1 = log(Fe+1)), se.fit = T)
  
  pred = pnorm(rnorm(1, init$`fit`, init$`se.fit`))
  
  pred / meso_mans_egg_fert_int
}

#plot viability over fertilizer concentration
meso_dat_mans_eggs %>% filter(Treatment %in% c("S", "F1x", "F1X", "F2x", "H2O", "W")) %>% 
  ggplot(aes(x = fert_conc_logp1, y = probit_per_hatch)) + geom_point() + theme_bw() +
  geom_smooth(method = "lm") + ylab("probit hatching percent") + xlab("log+1 Fertilizer (N-ppb)") +
  geom_point(data = data.frame(fert_conc_logp1 = log(seq(0,4400,40)+1),
                        prediction = qnorm(sapply(seq(0,4400,40), halstead_meso18_fert_mans_v_uncertainty)*meso_mans_egg_fert_int)), 
             aes(x = fert_conc_logp1, y = prediction), pch = 5, col = 4, cex = 0.7)


#haematobium data ************################################
meso_dat_haem_eggs <- meso_dat_eggs %>% filter(species == "haematobium")

#Atrazine ##########
#model
meso_haem_egg_atr <- lm(probit_per_hatch ~ atr_conc_logp1, weights = eggs,
                        data = meso_dat_haem_eggs %>% 
                          filter(Treatment %in% c("S", "A1x", "A2x", "H2O", "W") & 
                                   is.finite(probit_per_hatch)))

meso_haem_egg_atr_int <- pnorm(coef(meso_haem_egg_atr)[1])
  
#Atrazine function returns relative proportion change in viability from baseline
halstead_meso18_atr_haem_v_uncertainty<-function(He){
  init = predict(meso_haem_egg_atr, newdata = data.frame(atr_conc_logp1 = log(He+1)), se.fit = T)
  
  pred = pnorm(rnorm(1, init$`fit`, init$`se.fit`))
  
  pred / meso_haem_egg_atr_int
}

#plot viability over atrazine concentration with predictions
meso_dat_haem_eggs %>% filter(Treatment %in% c("S", "A1x", "A2x", "H2O", "W")) %>% 
  ggplot(aes(x = atr_conc_logp1, y = probit_per_hatch)) + geom_point() + theme_bw() +
  geom_smooth(method = "lm") + ylab("probit hatching percent") + xlab("log+1 Atrazine (ppb)") +
  geom_point(data = data.frame(atr_conc_logp1 = log(c(0:102)+1),
                        prediction = qnorm(sapply(c(0:102), halstead_meso18_atr_haem_v_uncertainty)*meso_haem_egg_atr_int)), 
             aes(x = atr_conc_logp1, y = prediction), pch = 5, col = 4, cex = 0.7)

#Chlorpyrifos ##########
#Chlorpyrifos model
meso_haem_egg_chlor<-lm(probit_per_hatch ~ chlor_conc_logp1, weights = eggs,
                        data = meso_dat_haem_eggs %>% 
                          filter(Treatment %in% c("S", "C1x", "C1X", "C2x", "H2O", "W") & 
                                   is.finite(probit_per_hatch)))

meso_haem_egg_chlor_int <- pnorm(coef(meso_haem_egg_chlor)[1])

#Chlorpyrifos function
halstead_meso18_chlor_haem_v_uncertainty<-function(In){
  init = predict(meso_haem_egg_chlor, newdata = data.frame(chlor_conc_logp1 = log(In+1)), se.fit = T)
  
  pred = pnorm(rnorm(1, init$`fit`, init$`se.fit`))
  
  pred / meso_haem_egg_chlor_int
}

#plot viability over chlorpyrifos concentration
meso_dat_haem_eggs %>% filter(Treatment %in% c("S", "C1x", "C1X", "C2x", "H2O", "W")) %>% 
  ggplot(aes(x = chlor_conc_logp1, y = probit_per_hatch)) + geom_point() + theme_bw() +
  geom_smooth(method = "lm") + ylab("probit hatching percent") + xlab("log+1 Chlorpyrifos (ppb)") +
  geom_point(data = data.frame(chlor_conc_logp1 = log(c(0:64)+1),
                        prediction = qnorm(sapply(c(0:64), halstead_meso18_chlor_haem_v_uncertainty)*meso_haem_egg_chlor_int)), 
             aes(x = chlor_conc_logp1, y = prediction), pch = 5, col = 4, cex = 0.7)

#Fertilizer ##########
#Fertilizer model
meso_haem_egg_fert<-lm(probit_per_hatch ~ fert_conc_logp1, weights = eggs,
                        data = meso_dat_haem_eggs %>% 
                          filter(Treatment %in% c("S", "F1x", "F1X", "F2x", "H2O", "W") & 
                                   is.finite(probit_per_hatch)))

meso_haem_egg_fert_int <- pnorm(coef(meso_haem_egg_fert)[1])

#Fertilizer function
halstead_meso18_fert_haem_v_uncertainty<-function(Fe){
  init = predict(meso_haem_egg_fert, newdata = data.frame(fert_conc_logp1 = log(Fe+1)), se.fit = T)
  
  pred = pnorm(rnorm(1, init$`fit`, init$`se.fit`))
  
  pred / meso_haem_egg_fert_int
}

#plot viability over fertilizer concentration
meso_dat_haem_eggs %>% filter(Treatment %in% c("S", "F1x", "F1X", "F2x", "H2O", "W")) %>% 
  ggplot(aes(x = fert_conc_logp1, y = probit_per_hatch)) + geom_point() + theme_bw() +
  geom_smooth(method = "lm") + ylab("probit hatching percent") + xlab("log+1 Fertilizer (N-ppb)") +
  geom_point(data = data.frame(fert_conc_logp1 = log(seq(0,4400,40)+1),
                        prediction = qnorm(sapply(seq(0,4400,40), halstead_meso18_fert_haem_v_uncertainty)*meso_haem_egg_fert_int)), 
             aes(x = fert_conc_logp1, y = prediction), pch = 5, col = 4, cex = 0.7)
