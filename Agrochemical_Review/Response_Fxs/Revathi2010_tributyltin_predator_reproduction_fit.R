#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(drc)

revathi_dat <- data.frame(hatch = c(96.17, 64.83, 38.30, 17.36)/100,
                          se = c(1.9, 1.85, 1.38, 1.4)/100,
                          tbt = c(0, 0.78, 1.56, 3.12))

revathi_ref <- revathi_dat$hatch[1]
  
revathi_mod= drm(hatch ~ tbt, data = revathi_dat, type = 'binomial',
                 fct = LL2.2())

  revathi_tbt_lc50 <- summary(revathi_mod)$coef[2,1]   # Lc50 estimate from model
  revathi_tbt_lc50_se <- summary(revathi_mod)$coef[2,2]  #SE of lc50 estimate from model
  revathi_tbt_b <- summary(revathi_mod)$coef[1,1]  #slope parameter of log-logistic fitted here
    
fPq_tbt_Revathi10_uncertainty<-function(He){
  e = rnorm(1, revathi_tbt_lc50, revathi_tbt_lc50_se)
  
  mort <- 1 / (1 + exp(revathi_tbt_b*(log(He/1000)-e)))
  
  mort
}
    
keep.tbt.Revathi10 = c('fPq_tbt_Revathi10_uncertainty', 'revathi_tbt_lc50', 'revathi_tbt_lc50_se', 'revathi_tbt_b')    
