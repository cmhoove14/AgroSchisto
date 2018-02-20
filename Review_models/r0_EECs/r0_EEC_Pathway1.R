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
source('Response_Fxs/fin/Baxter_Rohr_Atrazine2011_fin.R')
source('Response_Fxs/fin/Johnson07_theta_fN_fin.R')
source('Response_Fxs/fin/rohr08_nature_fN_fin.R')
source('Response_Fxs/fin/Halstead_meso2017_fin.R')

source('Review_models/fin/r0_of_q_fin.R')

today = Sys.Date()

library(parallel)
library(fBasics)
library(drc)

keep.fin.p1 = c(keep.baxrohr, keep.johnson07, 'rohr08_fN_uncertainty', 'rohr08_fN_uncertainty2', keep.halstead17,
                'r0.He', 'r0.Fe', 'r0.fix', 'parameters', 'nil0', 'nil1', 'keep.fin.p1', 'today', 'kmat')

rm(list = setdiff(ls(), keep.fin.p1))
dev.off()

no.cores = detectCores() - 1

#Run simulations for pathway 1 response functions ################
nsims = 5000         #Number of simulations to run
#Expected environmental concentrations for agrochemicals (NEED TO BE FINALIZED)
  eec.atr = rep(102, nsims)     #atrazine, EEC from Halstead et al 2017
  blnk.vec = c(1:nsims)
  
#All pathway 1 response functions  
parfx = c(phi_Nq_atr_baxrohr.no30, phi_Nq_atr_baxrohr, johnson07_theta_uncertainty, johnson07_phin_uncertainty, halstead17_phiN_at_uncertainty, halstead17_phiN_fe_uncertainty, rohr08_fN_uncertainty2)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.p1, 'uniroot.all', 'rdrm', 'LL.2', 'L.4'))

r0.eec.p1 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.eec.p1 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.0.5eec.p1 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.0.5eec.p1 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.0.1eec.p1 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.0.1eec.p1 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

clusterSetRNGStream(cl = clus1, iseed = 43093) #set cluster seed for reproducibility

#Fill r0 estimates for EEC #####################
  #individual parameters
    r0.eec.p1[, 1] = parSapply(clus1, eec.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.eec.p1[, 2] = parSapply(clus1, eec.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr)[3,]
    r0.eec.p1[, 3] = parSapply(clus1, blnk.vec, r0.He, f.phi_Nq = johnson07_phin_uncertainty)[3,]
    r0.eec.p1[, 4] = parSapply(clus1, blnk.vec, r0.He, f.theta_q = johnson07_theta_uncertainty)[3,]
    r0.eec.p1[, 5] = parSapply(clus1, blnk.vec, r0.He, f.phi_Nq = rohr08_fN_uncertainty2)[3,]
    r0.eec.p1[, 6] = parSapply(clus1, blnk.vec, r0.He, f.phi_Nq = halstead17_phiN_at_uncertainty)[3,]
    r0.eec.p1[, 7] = parSapply(clus1, blnk.vec, r0.He, f.phi_Nq = halstead17_phiN_fe_uncertainty)[3,]
    
    
  #fill parameter values for EEC values
    par.eec.p1[, 1] = parSapply(clus1, eec.atr, phi_Nq_atr_baxrohr.no30) 
    par.eec.p1[, 2] = parSapply(clus1, eec.atr, phi_Nq_atr_baxrohr) 
    par.eec.p1[, 3] = parSapply(clus1, blnk.vec, johnson07_phin_uncertainty) 
    par.eec.p1[, 4] = parSapply(clus1, blnk.vec, johnson07_theta_uncertainty) 
    par.eec.p1[, 5] = parSapply(clus1, blnk.vec, rohr08_fN_uncertainty2) 
    par.eec.p1[, 6] = parSapply(clus1, blnk.vec, halstead17_phiN_at_uncertainty) 
    par.eec.p1[, 7] = parSapply(clus1, blnk.vec, halstead17_phiN_fe_uncertainty) 
    
    
#Fill r0 estimates for 50% EEC ########################
  #individual parameters
    r0.0.5eec.p1[, 1] = parSapply(clus1, 0.5*eec.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.0.5eec.p1[, 2] = parSapply(clus1, 0.5*eec.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr)[3,]
    r0.0.5eec.p1[, 3] = parSapply(clus1, blnk.vec, r0.He, f.phi_Nq = johnson07_phin_uncertainty)[3,]
    r0.0.5eec.p1[, 4] = parSapply(clus1, blnk.vec, r0.He, f.theta_q = johnson07_theta_uncertainty)[3,]
    r0.0.5eec.p1[, 5] = parSapply(clus1, blnk.vec, r0.He, f.phi_Nq = rohr08_fN_uncertainty2)[3,]
    r0.0.5eec.p1[, 6] = parSapply(clus1, blnk.vec, r0.He, f.phi_Nq = halstead17_phiN_at_uncertainty)[3,]
    r0.0.5eec.p1[, 7] = parSapply(clus1, blnk.vec, r0.He, f.phi_Nq = halstead17_phiN_fe_uncertainty)[3,]
    
  #fill parameter values for EEC values
    par.0.5eec.p1[, 1] = parSapply(clus1, 0.5*eec.atr, phi_Nq_atr_baxrohr.no30) 
    par.0.5eec.p1[, 2] = parSapply(clus1, 0.5*eec.atr, phi_Nq_atr_baxrohr) 
    par.0.5eec.p1[, 3] = parSapply(clus1, blnk.vec, johnson07_phin_uncertainty) 
    par.0.5eec.p1[, 4] = parSapply(clus1, blnk.vec, johnson07_theta_uncertainty) 
    par.0.5eec.p1[, 5] = parSapply(clus1, blnk.vec, rohr08_fN_uncertainty2) 
    par.0.5eec.p1[, 6] = parSapply(clus1, blnk.vec, halstead17_phiN_at_uncertainty) 
    par.0.5eec.p1[, 7] = parSapply(clus1, blnk.vec, halstead17_phiN_fe_uncertainty) 
    
    
    
#Fill r0 estimates for 10% EEC #####################
  #individual parameters
    r0.0.1eec.p1[, 1] = parSapply(clus1, 0.1*eec.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.0.1eec.p1[, 2] = parSapply(clus1, 0.1*eec.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr)[3,]
    r0.0.1eec.p1[, 3] = parSapply(clus1, blnk.vec, r0.He, f.phi_Nq = johnson07_phin_uncertainty)[3,]
    r0.0.1eec.p1[, 4] = parSapply(clus1, blnk.vec, r0.He, f.theta_q = johnson07_theta_uncertainty)[3,]
    r0.0.1eec.p1[, 5] = parSapply(clus1, blnk.vec, r0.He, f.phi_Nq = rohr08_fN_uncertainty2)[3,]
    r0.0.1eec.p1[, 6] = parSapply(clus1, blnk.vec, r0.He, f.phi_Nq = halstead17_phiN_at_uncertainty)[3,]
    r0.0.1eec.p1[, 7] = parSapply(clus1, blnk.vec, r0.He, f.phi_Nq = halstead17_phiN_fe_uncertainty)[3,]
    
    
  #fill parameter values for EEC values
    par.0.1eec.p1[, 1] = parSapply(clus1, 0.1*eec.atr, phi_Nq_atr_baxrohr.no30) 
    par.0.1eec.p1[, 2] = parSapply(clus1, 0.1*eec.atr, phi_Nq_atr_baxrohr) 
    par.0.1eec.p1[, 3] = parSapply(clus1, blnk.vec, johnson07_phin_uncertainty) 
    par.0.1eec.p1[, 4] = parSapply(clus1, blnk.vec, johnson07_theta_uncertainty) 
    par.0.1eec.p1[, 5] = parSapply(clus1, blnk.vec, rohr08_fN_uncertainty2) 
    par.0.1eec.p1[, 6] = parSapply(clus1, blnk.vec, halstead17_phiN_at_uncertainty) 
    par.0.1eec.p1[, 7] = parSapply(clus1, blnk.vec, halstead17_phiN_fe_uncertainty) 
    
    
stopCluster(clus1) 
#Post process ############ 
#EEC runs ################
eec.p1.df = data.frame(chem = c('Atrazine', 'Atrazine','Ammonium Nitrate & Phosphoric Acid', 
                                'Ammonium Nitrate & Phosphoric Acid', 'Atrazine', 'Atrazine', 'Sodium Nitrate & Sodium phosphate'),
                       study = c('Rohr et al, 2012', 'Rohr et al, 2012', 'Johnson et al 2007', 
                                 'Johnson et al 2007', 'Rohr et al 2008', 'Halstead et al 2017', 'Halstead et al 2017'),
                       Species = c('Physella spp', 'Physella spp', 'Planorbella trivolvis', 
                                   'Planorbella trivolvis', 'Planorbella trivolvis', 'Bulinus truncatus', 'Bulinus truncatus'),
                       Parameter = c('phi_N', 'phi_N', 'phi_N', 'theta', 'phi_N', 'phi_N', 'phi_N'),
                       r0 = colMeans(r0.eec.p1),
                       r0.sd = apply(r0.eec.p1, 2, sd),
                       r0.med = apply(r0.eec.p1, 2, median),
                       r0.25 = apply(r0.eec.p1, 2, quantile, prob = 0.25),
                       r0.75 = apply(r0.eec.p1, 2, quantile, prob = 0.75),
                       par.mean = colMeans(par.eec.p1))

eec.p1.df$r0.up = eec.p1.df$r0 + eec.p1.df$r0.sd
eec.p1.df$r0.lo = eec.p1.df$r0 - eec.p1.df$r0.sd

eec.p1.df$deltar0 = eec.p1.df$r0 - r0.He()[3]
eec.p1.df$deltar0.up = (eec.p1.df$r0 + eec.p1.df$r0.sd) - r0.He()[3]
eec.p1.df$deltar0.lo = (eec.p1.df$r0 - eec.p1.df$r0.sd) - r0.He()[3]

eec.p1.df$relr0 = eec.p1.df$r0 / r0.He()[3] * 100 - 100
eec.p1.df$relr0.up = (eec.p1.df$r0 + eec.p1.df$r0.sd) / r0.He()[3] * 100 - 100
eec.p1.df$relr0.lo = (eec.p1.df$r0 - eec.p1.df$r0.sd) / r0.He()[3] * 100 - 100

eec.p1.df$relr0.med = eec.p1.df$r0.med / r0.He()[3] * 100 - 100
eec.p1.df$relr0.25 = eec.p1.df$r0.25 / r0.He()[3] * 100 - 100
eec.p1.df$relr0.75 = eec.p1.df$r0.75 / r0.He()[3] * 100 - 100

save(eec.p1.df, file = paste('Review_models/r0_EECs/eec.p1.df', today, '.RData', sep = ''))

#50% EEC values #################
eec0.5.p1.df = data.frame(chem = c('Atrazine', 'Atrazine','Ammonium Nitrate & Phosphoric Acid', 
                                   'Ammonium Nitrate & Phosphoric Acid', 'Atrazine', 'Atrazine', 'Sodium Nitrate & Sodium phosphate'),
                          study = c('Rohr et al, 2012', 'Rohr et al, 2012', 'Johnson et al 2007', 
                                    'Johnson et al 2007', 'Rohr et al 2008', 'Halstead et al 2017', 'Halstead et al 2017'),
                          Species = c('Physella spp', 'Physella spp', 'Planorbella trivolvis', 
                                      'Planorbella trivolvis', 'Planorbella trivolvis', 'Bulinus truncatus', 'Bulinus truncatus'),
                          Parameter = c('phi_N', 'phi_N', 'phi_N', 'theta', 'phi_N', 'phi_N', 'phi_N'),
                          r0 = colMeans(r0.0.5eec.p1),
                          r0.sd = apply(r0.0.5eec.p1, 2, sd),
                          r0.med = apply(r0.0.5eec.p1, 2, median),
                          r0.25 = apply(r0.0.5eec.p1, 2, quantile, prob = 0.25),
                          r0.75 = apply(r0.0.5eec.p1, 2, quantile, prob = 0.75),
                          par.mean = colMeans(par.0.5eec.p1))

eec0.5.p1.df$r0.up = eec0.5.p1.df$r0 + eec0.5.p1.df$r0.sd
eec0.5.p1.df$r0.lo = eec0.5.p1.df$r0 - eec0.5.p1.df$r0.sd

eec0.5.p1.df$deltar0 = eec0.5.p1.df$r0 - r0.He()[3]
eec0.5.p1.df$deltar0.up = (eec0.5.p1.df$r0 + eec0.5.p1.df$r0.sd) - r0.He()[3]
eec0.5.p1.df$deltar0.lo = (eec0.5.p1.df$r0 - eec0.5.p1.df$r0.sd) - r0.He()[3]

eec0.5.p1.df$relr0 = eec0.5.p1.df$r0 / r0.He()[3] * 100 - 100
eec0.5.p1.df$relr0.up = (eec0.5.p1.df$r0 + eec0.5.p1.df$r0.sd) / r0.He()[3] * 100 - 100
eec0.5.p1.df$relr0.lo = (eec0.5.p1.df$r0 - eec0.5.p1.df$r0.sd) / r0.He()[3] * 100 - 100

eec0.5.p1.df$relr0.med = eec0.5.p1.df$r0.med / r0.He()[3] * 100 - 100
eec0.5.p1.df$relr0.25 = eec0.5.p1.df$r0.25 / r0.He()[3] * 100 - 100
eec0.5.p1.df$relr0.75 = eec0.5.p1.df$r0.75 / r0.He()[3] * 100 - 100

save(eec0.5.p1.df, file = paste('Review_models/r0_EECs/eec0.5.p1.df', today, '.RData', sep = ''))


#10% EEC values ######################
eec0.1.p1.df = data.frame(chem = c('Atrazine', 'Atrazine','Ammonium Nitrate & Phosphoric Acid', 
                                   'Ammonium Nitrate & Phosphoric Acid', 'Atrazine', 'Atrazine', 'Sodium Nitrate & Sodium phosphate'),
                          study = c('Rohr et al, 2012', 'Rohr et al, 2012', 'Johnson et al 2007', 
                                    'Johnson et al 2007', 'Rohr et al 2008', 'Halstead et al 2017', 'Halstead et al 2017'),
                          Species = c('Physella spp', 'Physella spp', 'Planorbella trivolvis', 
                                      'Planorbella trivolvis', 'Planorbella trivolvis', 'Bulinus truncatus', 'Bulinus truncatus'),
                          Parameter = c('phi_N', 'phi_N', 'phi_N', 'theta', 'phi_N', 'phi_N', 'phi_N'),
                          r0 = colMeans(r0.0.1eec.p1),
                          r0.sd = apply(r0.0.1eec.p1, 2, sd),
                          r0.med = apply(r0.0.1eec.p1, 2, median),
                          r0.25 = apply(r0.0.1eec.p1, 2, quantile, prob = 0.25),
                          r0.75 = apply(r0.0.1eec.p1, 2, quantile, prob = 0.75),
                          par.mean = colMeans(par.0.1eec.p1))

eec0.1.p1.df$r0.up = eec0.1.p1.df$r0 + eec0.1.p1.df$r0.sd
eec0.1.p1.df$r0.lo = eec0.1.p1.df$r0 - eec0.1.p1.df$r0.sd

eec0.1.p1.df$deltar0 = eec0.1.p1.df$r0 - r0.He()[3]
eec0.1.p1.df$deltar0.up = (eec0.1.p1.df$r0 + eec0.1.p1.df$r0.sd) - r0.He()[3]
eec0.1.p1.df$deltar0.lo = (eec0.1.p1.df$r0 - eec0.1.p1.df$r0.sd) - r0.He()[3]

eec0.1.p1.df$relr0 = eec0.1.p1.df$r0 / r0.He()[3] * 100 - 100
eec0.1.p1.df$relr0.up = (eec0.1.p1.df$r0 + eec0.1.p1.df$r0.sd) / r0.He()[3] * 100 - 100
eec0.1.p1.df$relr0.lo = (eec0.1.p1.df$r0 - eec0.1.p1.df$r0.sd) / r0.He()[3] * 100 - 100

eec0.1.p1.df$relr0.med = eec0.1.p1.df$r0.med / r0.He()[3] * 100 - 100
eec0.1.p1.df$relr0.25 = eec0.1.p1.df$r0.25 / r0.He()[3] * 100 - 100
eec0.1.p1.df$relr0.75 = eec0.1.p1.df$r0.75 / r0.He()[3] * 100 - 100

save(eec0.1.p1.df, file = paste('Review_models/r0_EECs/eec0.1.p1.df', today, '.RData', sep = ''))
