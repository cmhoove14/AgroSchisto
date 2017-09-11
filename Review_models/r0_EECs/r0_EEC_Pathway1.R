#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Load response and R0 functions #######
source('Response_Fxs/Baxter_Rohr_Atrazine2011.R')
source('Response_Fxs/Johnson07_theta_fN.R')
source('Response_Fxs/rohr08_nature_fN.R')

source('Review_models/r0_of_q.R')

library(parallel)
library(fBasics)

keep.fin.p1 = c(keep.baxrohr, keep.johnson07, 'rohr08_fN_uncertainty',
                'r0.He', 'r0.Fe', 'r0.fix', 'parameters', 'nil0', 'nil1', 'keep.fin.p1')

rm(list = setdiff(ls(), keep.fin.p1))
dev.off()

no.cores = detectCores() - 1

#Run simulations of malathion concentrations, start with individual functions then combine ################
nsims = 1000         #Number of simulations to run
#Expected environmental concentrations for agrochemicals (NEED TO BE FINALIZED)
  eec.atr = rep(102, nsims)     #atrazine
  blnk.vec = c(1:nsims)
  
#All pathway 2 response functions  
parfx = c(phi_Nq_atr_baxrohr.no30, johnson07_fN_uncertainty, johnson07_theta_uncertainty, rohr08_fN_uncertainty)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.p1, 'uniroot.all', 'rdrm', 'LL.2', 'L.4'))

r0.eec.p1 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.eec.p1 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.0.5eec.p1 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.0.5eec.p1 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.0.1eec.p1 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.0.1eec.p1 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

set.seed(0)

#Fill r0 estimates for EEC #####################
  #individual parameters
    r0.eec.p1[, 1] = parSapply(clus1, eec.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.eec.p1[, 2] = parSapply(clus1, blnk.vec, r0.He, f.f_Nq = johnson07_fN_uncertainty)[3,]
    r0.eec.p1[, 3] = parSapply(clus1, blnk.vec, r0.He, f.theta_q = johnson07_theta_uncertainty)[3,]
    r0.eec.p1[, 4] = parSapply(clus1, blnk.vec, r0.He, f.f_Nq = rohr08_fN_uncertainty)[3,]
    
    
  #fill parameter values for EEC values
    par.eec.p1[, 1] = parSapply(clus1, eec.atr, phi_Nq_atr_baxrohr.no30) 
    par.eec.p1[, 2] = parSapply(clus1, blnk.vec, johnson07_fN_uncertainty) 
    par.eec.p1[, 3] = parSapply(clus1, blnk.vec, johnson07_theta_uncertainty) 
    par.eec.p1[, 4] = parSapply(clus1, blnk.vec, rohr08_fN_uncertainty) 
    
    
#Fill r0 estimates for 50% EEC ########################
  #individual parameters
    r0.0.5eec.p1[, 1] = parSapply(clus1, 0.5*eec.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.0.5eec.p1[, 2] = parSapply(clus1, blnk.vec, r0.He, f.f_Nq = johnson07_fN_uncertainty)[3,]
    r0.0.5eec.p1[, 3] = parSapply(clus1, blnk.vec, r0.He, f.theta_q = johnson07_theta_uncertainty)[3,]
    r0.0.5eec.p1[, 4] = parSapply(clus1, blnk.vec, r0.He, f.f_Nq = rohr08_fN_uncertainty)[3,]
    
    
  #fill parameter values for EEC values
    par.0.5eec.p1[, 1] = parSapply(clus1, 0.5*eec.atr, phi_Nq_atr_baxrohr.no30) 
    par.0.5eec.p1[, 2] = parSapply(clus1, blnk.vec, johnson07_fN_uncertainty) 
    par.0.5eec.p1[, 3] = parSapply(clus1, blnk.vec, johnson07_theta_uncertainty) 
    par.0.5eec.p1[, 4] = parSapply(clus1, blnk.vec, rohr08_fN_uncertainty) 
    
    
    
#Fill r0 estimates for 10% EEC #####################
  #individual parameters
    r0.0.1eec.p1[, 1] = parSapply(clus1, 0.1*eec.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.0.1eec.p1[, 2] = parSapply(clus1, blnk.vec, r0.He, f.f_Nq = johnson07_fN_uncertainty)[3,]
    r0.0.1eec.p1[, 3] = parSapply(clus1, blnk.vec, r0.He, f.theta_q = johnson07_theta_uncertainty)[3,]
    r0.0.1eec.p1[, 4] = parSapply(clus1, blnk.vec, r0.He, f.f_Nq = rohr08_fN_uncertainty)[3,]
    
    
  #fill parameter values for EEC values
    par.0.1eec.p1[, 1] = parSapply(clus1, 0.1*eec.atr, phi_Nq_atr_baxrohr.no30) 
    par.0.1eec.p1[, 2] = parSapply(clus1, blnk.vec, johnson07_fN_uncertainty) 
    par.0.1eec.p1[, 3] = parSapply(clus1, blnk.vec, johnson07_theta_uncertainty) 
    par.0.1eec.p1[, 4] = parSapply(clus1, blnk.vec, rohr08_fN_uncertainty) 
    
    
stopCluster(clus1) 
#Post process ############ 
#EEC runs ################
eec.p1.df = data.frame(chem = c('Atrazine', 'Ammonium Nitrate & Phosphoric Acid', 'Ammonium Nitrate & Phosphoric Acid', 'Atrazine'),
                       study = c('Rohr et al, 2012', 'Johnson et al 2007', 'Johnson et al 2007', 'Rohr et al 2008'),
                       Species = c('Physella spp', 'Planorbella trivolvis', 'Planorbella trivolvis', 'Planorbella trivolvis'),
                       Parameter = c('phi_N', 'fN', 'theta', 'fN'),
                       r0 = colMeans(r0.eec.p1),
                       r0.sd = apply(r0.eec.p1, 2, sd),
                       par.mean = colMeans(par.eec.p1))

eec.p1.df$r0.up = eec.p1.df$r0 + eec.p1.df$r0.sd
eec.p1.df$r0.lo = eec.p1.df$r0 - eec.p1.df$r0.sd

save(eec.p1.df, file = 'Review_models/r0_EECs/eec.p1.df.RData')

#50% EEC values #################
eec0.5.p1.df = data.frame(chem = c('Atrazine', 'Ammonium Nitrate & Phosphoric Acid', 'Ammonium Nitrate & Phosphoric Acid', 'Atrazine'),
                          study = c('Rohr et al, 2012', 'Johnson et al 2007', 'Johnson et al 2007', 'Rohr et al 2008'),
                          Species = c('Physella spp', 'Planorbella trivolvis', 'Planorbella trivolvis', 'Planorbella trivolvis'),
                          Parameter = c('phi_N', 'fN', 'theta', 'fN'),
                          r0 = colMeans(r0.0.5eec.p1),
                          r0.sd = apply(r0.0.5eec.p1, 2, sd),
                          par.mean = colMeans(par.0.5eec.p1))

eec0.5.p1.df$r0.up = eec0.5.p1.df$r0 + eec0.5.p1.df$r0.sd
eec0.5.p1.df$r0.lo = eec0.5.p1.df$r0 - eec0.5.p1.df$r0.sd

save(eec0.5.p1.df, file = 'Review_models/r0_EECs/eec0.5.p1.df.RData')


#10% EEC values ######################
eec0.1.p1.df = data.frame(chem = c('Atrazine', 'Ammonium Nitrate & Phosphoric Acid', 'Ammonium Nitrate & Phosphoric Acid', 'Atrazine'),
                          study = c('Rohr et al, 2012', 'Johnson et al 2007', 'Johnson et al 2007', 'Rohr et al 2008'),
                          Species = c('Physella spp', 'Planorbella trivolvis', 'Planorbella trivolvis', 'Planorbella trivolvis'),
                          Parameter = c('phi_N', 'fN', 'theta', 'fN'),
                          r0 = colMeans(r0.0.1eec.p1),
                          r0.sd = apply(r0.0.1eec.p1, 2, sd),
                          par.mean = colMeans(par.0.1eec.p1))

eec0.1.p1.df$r0.up = eec0.1.p1.df$r0 + eec0.1.p1.df$r0.sd
eec0.1.p1.df$r0.lo = eec0.1.p1.df$r0 - eec0.1.p1.df$r0.sd

save(eec0.1.p1.df, file = 'Review_models/r0_EECs/eec0.1.p1.df.RData')
