#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(LW1949)
require(drc)
require(tidyverse)

hua_dat <- read.csv("Agrochemical_Review/Response_Fxs/Data/hua_etal_insecticide_cercariae_data.csv")
  hua_dat$ntot <- 40
  hua_dat$nfx <- round(hua_dat$mort/100 * hua_dat$ntot)
  hua_dat$probit <- qnorm(hua_dat$mort/100)
  hua_dat$log10_conc <- log10(hua_dat$conc)

ggplot(data = hua_dat, aes(x = log10_conc, y = probit, col = chem)) + 
  geom_point() + 
  theme_bw() + 
  geom_smooth(method = "lm") +
  facet_grid(~pop)

#Add mean estimate of mortality between the two populations
hua_dat <- hua_dat %>% 
  group_by(chem, conc) %>% 
  mutate(mean_mort = mean(mort),
         mean_nfx = round(mean_mort/100 * ntot),
         mean_prob = qnorm(mean_mort/100)) %>% 
  ungroup()

ggplot(data = hua_dat, aes(x = log10_conc, y = mean_prob, col = chem)) + 
  geom_point() + 
  theme_bw() + 
  geom_smooth(method = "lm")

#they tested insecticide toxicity to two populations: PWA is from a pond far from any agricultural influence and ICP is from a pond next to AG. We'll use PWA since we don't want to consider the role of adaptation in any quantitative analysis here. Also have a feeling that within the small concentration they investigated, they may be picking up more natural variability and/or error in their actual concentrations than any sort of direct effect of the insecticides themselves. In addition, they should really test a range of concnetrations that spans from no effect to complete effect (e.g. 100% mortality within 6 hours) in order to accurately estimate a dose-response relationship. Will proceed noneltheless

hua_pwa <- subset(hua_dat, pop == "PWA")

#carbaryl function NOT CONVERGING
hua_carb <- subset(hua_pwa, chem == "carbaryl")

#lw1949 method  
  hua_carb_lw49 <- dataprep(hua_carb$conc, hua_carb$ntot, hua_carb$nfx)
  hua_carb_fit <- fitLWauto(hua_carb_lw49)
  hua_carb_est <- LWestimate(hua_carb_fit, hua_carb_lw49)

    hua_carb_lc50 <- as.numeric(hua_carb_est$LWest[1])
      hua_carb_lc50_se <- as.numeric(log10(hua_carb_est$LWest[3]/hua_carb_est$LWest[1])) / qnorm(0.975)
    hua_carb_slp <- as.numeric(hua_carb_est$params[2])
    
  piCq_hua_carb_uncertainty <- function(In){
    lc50 = 10^rnorm(1, log10(hua_carb_lc50), hua_carb_lc50_se)
    
    1 - pnorm(hua_carb_slp * log10(In/lc50))
  }  

plot(hua_carb$conc, 1-hua_carb$mort/100, ylim = c(0,1), xlim = c(0,100), pch = 16,
     xlab = "carbaryl (ppb)", ylab = "relative 6-hr survival")    
  points(seq(0,100,0.5), sapply(seq(0,100,0.5), piCq_hua_carb_uncertainty), pch = 5, cex = 0.5, col = 4)

#malathion function
hua_mal <- subset(hua_pwa, chem == "malathion")

#lw1949 method  
  hua_mal_lw49 <- dataprep(hua_mal$conc, hua_mal$ntot, hua_mal$nfx)
  hua_mal_fit <- fitLWauto(hua_mal_lw49)
  hua_mal_est <- LWestimate(hua_mal_fit, hua_mal_lw49)
  
    hua_mal_lc50 <- as.numeric(hua_mal_est$LWest[1])
      hua_mal_lc50_se <- as.numeric(log10(hua_mal_est$LWest[3]/hua_mal_est$LWest[1])) / qnorm(0.975)
    hua_mal_slp <- as.numeric(hua_mal_est$params[2])
    
  piCq_hua_mal_uncertainty <- function(In){
    lc50 = 10^rnorm(1, log10(hua_mal_lc50), hua_mal_lc50_se)
    
    1 - pnorm(hua_mal_slp * log10(In/lc50))
  }  

plot(hua_mal$conc, 1-hua_mal$mort/100, ylim = c(0,1), xlim = c(0,100), pch = 16,
     xlab = "malathion (ppb)", ylab = "relative 6-hr survival")    
  points(seq(0,100,0.5), sapply(seq(0,100,0.5), piCq_hua_mal_uncertainty), pch = 5, cex = 0.5, col = 4)

#permethrin function
hua_perm <- subset(hua_pwa, chem == "permethrin")

#lw1949 method  
  hua_perm_lw49 <- dataprep(hua_perm$conc, hua_perm$ntot, hua_perm$nfx)
  hua_perm_fit <- fitLWauto(hua_perm_lw49)
  hua_perm_est <- LWestimate(hua_perm_fit, hua_perm_lw49)
  
    hua_perm_lc50 <- as.numeric(hua_perm_est$LWest[1])
      hua_perm_lc50_se <- as.numeric(log10(hua_perm_est$LWest[3]/hua_perm_est$LWest[1])) / qnorm(0.975)
    hua_perm_slp <- as.numeric(hua_perm_est$params[2])
    
  piCq_hua_perm_uncertainty <- function(In){
    lc50 = 10^rnorm(1, log10(hua_perm_lc50), hua_perm_lc50_se)
    
    1 - pnorm(hua_perm_slp * log10(In/lc50))
  }  

plot(hua_perm$conc, 1-hua_perm$mort/100, ylim = c(0,1), xlim = c(0,100), pch = 16,
     xlab = "permethrin (ppb)", ylab = "relative 6-hr survival")    
  points(seq(0,100,0.5), sapply(seq(0,100,0.5), piCq_hua_perm_uncertainty), pch = 5, cex = 0.5, col = 4)
  
#cypermethrin function
hua_cyp <- subset(hua_pwa, chem == "cypermethrin")

#lw1949 method  
  hua_cyp_lw49 <- dataprep(hua_cyp$conc, hua_cyp$ntot, hua_cyp$nfx)
  hua_cyp_fit <- fitLWauto(hua_cyp_lw49)
  hua_cyp_est <- LWestimate(hua_cyp_fit, hua_cyp_lw49)
  
    hua_cyp_lc50 <- as.numeric(hua_cyp_est$LWest[1])
      hua_cyp_lc50_se <- as.numeric(log10(hua_cyp_est$LWest[3]/hua_cyp_est$LWest[1])) / qnorm(0.975)
    hua_cyp_slp <- as.numeric(hua_cyp_est$params[2])
    
  piCq_hua_cyp_uncertainty <- function(In){
    lc50 = 10^rnorm(1, log10(hua_cyp_lc50), hua_cyp_lc50_se)
    
    1 - pnorm(hua_cyp_slp * log10(In/lc50))
  }  

plot(hua_cyp$conc, 1-hua_cyp$mort/100, ylim = c(0,1), xlim = c(0,100), pch = 16,
     xlab = "cypermethrin (ppb)", ylab = "relative 6-hr survival")    
  points(seq(0,100,0.5), sapply(seq(0,100,0.5), piCq_hua_cyp_uncertainty), pch = 5, cex = 0.5, col = 4)
  
#imidacloprid function
hua_imd <- subset(hua_pwa, chem == "imidacloprid")

#lw1949 method  
  hua_imd_lw49 <- dataprep(hua_imd$conc, hua_imd$ntot, hua_imd$nfx)
  hua_imd_fit <- fitLWauto(hua_imd_lw49)
  hua_imd_est <- LWestimate(hua_imd_fit, hua_imd_lw49)
  
    hua_imd_lc50 <- as.numeric(hua_imd_est$LWest[1])
      hua_imd_lc50_se <- as.numeric(log10(hua_imd_est$LWest[3]/hua_imd_est$LWest[1])) / qnorm(0.975)
    hua_imd_slp <- as.numeric(hua_imd_est$params[2])
    
  piCq_hua_imd_uncertainty <- function(In){
    lc50 = 10^rnorm(1, log10(hua_imd_lc50), hua_imd_lc50_se)
    
    1 - pnorm(hua_imd_slp * log10(In/lc50))
  }  

plot(hua_imd$conc, 1-hua_imd$mort/100, ylim = c(0,1), xlim = c(0,100), pch = 16,
     xlab = "imidacloprid (ppb)", ylab = "relative 6-hr survival")    
  points(seq(0,100,0.5), sapply(seq(0,100,0.5), piCq_hua_imd_uncertainty), pch = 5, cex = 0.5, col = 4)
   
#thiamethoxam function
hua_thi <- subset(hua_pwa, chem == "thiamethoxam")

#lw1949 method  
  hua_thi_lw49 <- dataprep(hua_thi$conc, hua_thi$ntot, hua_thi$nfx)
  hua_thi_fit <- fitLWauto(hua_thi_lw49)
  hua_thi_est <- LWestimate(hua_thi_fit, hua_thi_lw49)
  
    hua_thi_lc50 <- as.numeric(hua_thi_est$LWest[1])
      hua_thi_lc50_se <- as.numeric(log10(hua_thi_est$LWest[3]/hua_thi_est$LWest[1])) / qnorm(0.975)
    hua_thi_slp <- as.numeric(hua_thi_est$params[2])
    
  piCq_hua_thi_uncertainty <- function(In){
    lc50 = 10^rnorm(1, log10(hua_thi_lc50), hua_thi_lc50_se)
    
    1 - pnorm(hua_thi_slp * log10(In/lc50))
  }  

plot(hua_thi$conc, 1-hua_thi$mort/100, ylim = c(0,1), xlim = c(0,100), pch = 16,
     xlab = "thiamethoxam (ppb)", ylab = "relative 6-hr survival")    
  points(seq(0,100,0.5), sapply(seq(0,100,0.5), piCq_hua_thi_uncertainty), pch = 5, cex = 0.5, col = 4)
   