#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Seeland et al 2013 pyrimethanil toxicity to snails#########
#Paper reports LC50, LC10, NOEC, and LOEC values. To obtain d-r function, slope is estimated
#by linear regression fit to the LC50 and LC10 values and incorporated with reported LC50 and its uncertainty
#All data points are then used in plots for qualitative validation of the resulting function

require(drc)

#Reproduction ############
Seel13_repro_dat <- read.csv("Agrochemical_Review/Response_Fxs/Data/Seeland2013_snail_repro.csv")
  Seel13_repro_dat$probit <- qnorm(Seel13_repro_dat$mort)
  Seel13_repro_dat$log10ppm <- log10(Seel13_repro_dat$conc)
  
#y=mx+b  
Seel13_repro_mod <- lm(probit ~ log10ppm, data = Seel13_repro_dat[-c(1,2),])

  Seel13_repro_b <- summary(Seel13_repro_mod)$coef[1,1]   # Intercept
  Seel13_repro_b_se <- summary(Seel13_repro_mod)$coef[1,2]  #SE of intercept estimate from model
  Seel13_repro_m <- summary(Seel13_repro_mod)$coef[2,1]  #slope parameter
    
fNq_pyrimethanil_Seeland13_uncertainty<-function(Fg){

  b = rnorm(1, Seel13_repro_b, Seel13_repro_b_se)
  
  mort <- 1 - pnorm(Seel13_repro_m*log10(Fg/1000) + b)
  
  mort
}
  
#Mortality ############
Seel13_mort_dat <- data.frame(conc = c(0,0.06,0.12,0.25,0.5,1.0), #ppm
                              mort = c(6.6, 13.2, 6.6, 6.6, 19.8, 46.4)/100,
                              tot = rep(15, 6))

Seel13_mort_dat$dead <- round(Seel13_mort_dat$mort*Seel13_mort_dat$tot)

Seel13_pyr_mod <- drm(dead / tot ~ conc, 
                      weights = tot, type = 'binomial',  
                      data = Seel13_mort_dat, fct = LL2.2())
  
  Seel13_pyr_lc50 <- summary(Seel13_pyr_mod)$coef[2,1]   # Lc50 estimate from model
  Seel13_pyr_lc50_se <- summary(Seel13_pyr_mod)$coef[2,2]  #SE of lc50 estimate from model
  Seel13_pyr_b <- summary(Seel13_pyr_mod)$coef[1,1]  #slope parameter of log-logistic fitted here
    
muNq_pyrimethanil_Seeland13_uncertainty<-function(Fg){

  e = rnorm(1, Seel13_pyr_lc50, Seel13_pyr_lc50_se)
  
  mort <- 1 / (1 + exp(Seel13_pyr_b*(log(Fg/1000)-e)))
  
  mort
}

#keep vector
keep.Seeland.pyri = c('muNq_pyrimethanil_Seeland13_uncertainty', 'Seel13_pyr_lc50', 'Seel13_pyr_lc50_se', 'Seel13_pyr_b',
                      'fNq_pyrimethanil_Seeland13_uncertainty', 'Seel13_repro_m', 'Seel13_repro_b', 'Seel13_repro_b_se')    
