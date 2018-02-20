#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Load parameter sets and r0 functions #######
load("Review_models/r0_EECs/eec_par_values_2-19-18.RData")
load("Review_models/r0_EECs/eec10percent_par_values_2-19-18.RData")  

source('Review_models/fin/r0_of_q_fin.R')
  #r0.In()    #Check to make sure baseline r0 = 2.95

today = Sys.Date()

library(parallel)
library(fBasics)
library(drc)

dev.off()

no.cores = detectCores() - 1

#Run simulations for pathway 1 response functions ################
clus1 = makeCluster(no.cores)
clusterExport(clus1, c('uniroot.all', 'rdrm', 'LL.2', 'L.4','r0.He', 'r0.Fe', 'r0.fix', 'parameters', 'nil0', 'nil1', 'par.eec.all', 'par.0.1eec.all', 'today', 'kmat'))

r0.eec.p1 = matrix(data = NA, nrow = nrow(par.eec.all), ncol = 7)

r0.0.1eec.p1 = matrix(data = NA, nrow = nrow(par.eec.all), ncol = 7)

clusterSetRNGStream(cl = clus1, iseed = 43093) #set cluster seed for reproducibility

#Fill r0 estimates for EEC #####################
  #individual parameters
    r0.eec.p1[, 1] = parSapply(clus1, par.eec.all[,1], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               #phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p1[, 2] = parSapply(clus1, par.eec.all[,2], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               #phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p1[, 3] = parSapply(clus1, par.eec.all[,3], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               #phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p1[, 4] = parSapply(clus1, par.eec.all[,4], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               #theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p1[, 5] = parSapply(clus1, par.eec.all[,5], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               #phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p1[, 6] = parSapply(clus1, par.eec.all[,6], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               #phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p1[, 7] = parSapply(clus1, par.eec.all[,7], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               #phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    
    
#Fill r0 estimates for 10% EEC #####################
  #individual parameters
   r0.0.1eec.p1[, 1] = parSapply(clus1, par.0.1eec.all[,1], FUN = r0.fix, 
                                 f_Nqx = 1, 
                                 mu_Pqx = 0,
                                 #phi_Nqx = 1, 
                                 mu_Nqx = 0, 
                                 alpha_qx = 1,
                                 theta_qx = 1, 
                                 pi_Mqx = 1, 
                                 pi_Cqx = 1, 
                                 v_qx = 1)[3,]
   r0.0.1eec.p1[, 2] = parSapply(clus1, par.0.1eec.all[,2], FUN = r0.fix, 
                                 f_Nqx = 1, 
                                 mu_Pqx = 0,
                                 #phi_Nqx = 1, 
                                 mu_Nqx = 0, 
                                 alpha_qx = 1,
                                 theta_qx = 1, 
                                 pi_Mqx = 1, 
                                 pi_Cqx = 1, 
                                 v_qx = 1)[3,]
   r0.0.1eec.p1[, 3] = parSapply(clus1, par.0.1eec.all[,3], FUN = r0.fix, 
                                 f_Nqx = 1, 
                                 mu_Pqx = 0,
                                 #phi_Nqx = 1, 
                                 mu_Nqx = 0, 
                                 alpha_qx = 1,
                                 theta_qx = 1, 
                                 pi_Mqx = 1, 
                                 pi_Cqx = 1, 
                                 v_qx = 1)[3,]
   r0.0.1eec.p1[, 4] = parSapply(clus1, par.0.1eec.all[,4], FUN = r0.fix, 
                                 f_Nqx = 1, 
                                 mu_Pqx = 0,
                                 phi_Nqx = 1, 
                                 mu_Nqx = 0, 
                                 alpha_qx = 1,
                                 #theta_qx = 1, 
                                 pi_Mqx = 1, 
                                 pi_Cqx = 1, 
                                 v_qx = 1)[3,]
   r0.0.1eec.p1[, 5] = parSapply(clus1, par.0.1eec.all[,5], FUN = r0.fix, 
                                 f_Nqx = 1, 
                                 mu_Pqx = 0,
                                 #phi_Nqx = 1, 
                                 mu_Nqx = 0, 
                                 alpha_qx = 1,
                                 theta_qx = 1, 
                                 pi_Mqx = 1, 
                                 pi_Cqx = 1, 
                                 v_qx = 1)[3,]
   r0.0.1eec.p1[, 6] = parSapply(clus1, par.0.1eec.all[,6], FUN = r0.fix, 
                                 f_Nqx = 1, 
                                 mu_Pqx = 0,
                                 #phi_Nqx = 1, 
                                 mu_Nqx = 0, 
                                 alpha_qx = 1,
                                 theta_qx = 1, 
                                 pi_Mqx = 1, 
                                 pi_Cqx = 1, 
                                 v_qx = 1)[3,]
   r0.0.1eec.p1[, 7] = parSapply(clus1, par.0.1eec.all[,7], FUN = r0.fix, 
                                 f_Nqx = 1, 
                                 mu_Pqx = 0,
                                 #phi_Nqx = 1, 
                                 mu_Nqx = 0, 
                                 alpha_qx = 1,
                                 theta_qx = 1, 
                                 pi_Mqx = 1, 
                                 pi_Cqx = 1, 
                                 v_qx = 1)[3,]    
    
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
                       par.mean = colMeans(par.eec.all[,1:7]))

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
                          par.mean = colMeans(par.0.1eec.all[,1:7]))

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
