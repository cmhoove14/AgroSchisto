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
#r0.In()    #Check to make sure baseline r0 = 1.86

today = Sys.Date()

library(parallel)
library(fBasics)
library(drc)

dev.off()

no.cores = detectCores() - 1

#Run simulations for pathway 3 response functions ################
clus1 = makeCluster(no.cores)
clusterExport(clus1, c('uniroot.all', 'rdrm', 'LL.2', 'LL.3','L.4','r0.He', 'r0.Fe', 'r0.fix', 'parameters', 'nil0', 'nil1', 'par.eec.all', 'par.0.1eec.all', 'today', 'kmat'))

r0.eec.p3 = matrix(data = NA, nrow = nrow(par.eec.all), ncol = 12)

r0.0.1eec.p3 = matrix(data = NA, nrow = nrow(par.eec.all), ncol = 12)

clusterSetRNGStream(cl = clus1, iseed = 43093) #set cluster seed for reproducibility

#Fill r0 estimates for EEC #####################
  #individual parameters
    r0.eec.p3[, 1] = parSapply(clus1, par.eec.all[,37], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p3[, 2] = parSapply(clus1, par.eec.all[,38], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p3[, 3] = parSapply(clus1, par.eec.all[,39], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p3[, 4] = parSapply(clus1, par.eec.all[,40], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p3[, 5] = parSapply(clus1, par.eec.all[,41], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p3[, 6] = parSapply(clus1, par.eec.all[,42], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p3[, 7] = parSapply(clus1, par.eec.all[,43], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p3[, 8] = parSapply(clus1, par.eec.all[,44], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p3[, 9] = parSapply(clus1, par.eec.all[,45], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p3[, 10] = parSapply(clus1, par.eec.all[,46], FUN = r0.fix, 
                                f_Nqx = 1, 
                                #mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p3[, 11] = parSapply(clus1, par.eec.all[,47], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                #alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p3[, 12] = parSapply(clus1, par.eec.all[,48], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                #alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    
#Fill r0 estimates for 10% EEC #####################
  #individual parameters
    r0.0.1eec.p3[, 1] = parSapply(clus1, par.0.1eec.all[,37], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p3[, 2] = parSapply(clus1, par.0.1eec.all[,38], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p3[, 3] = parSapply(clus1, par.0.1eec.all[,39], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p3[, 4] = parSapply(clus1, par.0.1eec.all[,40], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p3[, 5] = parSapply(clus1, par.0.1eec.all[,41], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p3[, 6] = parSapply(clus1, par.0.1eec.all[,42], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p3[, 7] = parSapply(clus1, par.0.1eec.all[,43], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p3[, 8] = parSapply(clus1, par.0.1eec.all[,44], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p3[, 9] = parSapply(clus1, par.0.1eec.all[,45], FUN = r0.fix, 
                               f_Nqx = 1, 
                               #mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p3[, 10] = parSapply(clus1, par.0.1eec.all[,46], FUN = r0.fix, 
                                f_Nqx = 1, 
                                #mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p3[, 11] = parSapply(clus1, par.0.1eec.all[,47], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                #alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p3[, 12] = parSapply(clus1, par.0.1eec.all[,48], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                #alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    
  stopCluster(clus1) 
#Post process ############ 
#EEC runs ################
eec.p3.df = data.frame(chem = c('Chlorpyrifos', 'Chlorpyrifos', 'Dimethoate', 'Esfenvalerate', 'Lamda-Cyhalothrin',
                                'Malathion', 'Permethrin', 'Profenofos', 'Terbufos', 'Zinc', 'Chlorpyrifos', 'Zinc'),
                       study = c('Satapornvanit et al 2009', 'Halstead et al 2015', 'Satapornvanit et al 2009', 
                                 rep('Halstead et al 2015', 4), 'Satapornvanit et al 2009', 'Halstead et al 2015',
                                 rep('Satapornvanit et al 2009', 3)),
                       Species = c('Macrobrachium rosenbergii', 'Procambarus clarkii', 'Macrobrachium rosenbergii', 
                                   rep('Procambarus clarkii', 4), 'Macrobrachium rosenbergii', 'Procambarus clarkii',
                                   rep('Macrobrachium rosenbergii', 3)),
                       Parameter = c(rep('muP', 10), rep('psi', 2)),
                       r0 = colMeans(r0.eec.p3),
                       r0.sd = apply(r0.eec.p3, 2, sd),
                       r0.med = apply(r0.eec.p3, 2, median),
                       r0.25 = apply(r0.eec.p3, 2, quantile, prob = 0.25),
                       r0.75 = apply(r0.eec.p3, 2, quantile, prob = 0.75),
                       par.mean = colMeans(par.eec.all[,37:48]))

eec.p3.df$r0.up = eec.p3.df$r0 + eec.p3.df$r0.sd
eec.p3.df$r0.lo = eec.p3.df$r0 - eec.p3.df$r0.sd

eec.p3.df$deltar0 = eec.p3.df$r0 - r0.In()[3]
eec.p3.df$deltar0.up = (eec.p3.df$r0 + eec.p3.df$r0.sd) - r0.In()[3]
eec.p3.df$deltar0.lo = (eec.p3.df$r0 - eec.p3.df$r0.sd) - r0.In()[3]

eec.p3.df$relr0 = eec.p3.df$r0 / r0.In()[3] * 100 - 100
eec.p3.df$relr0.up = (eec.p3.df$r0 + eec.p3.df$r0.sd) / r0.In()[3] * 100 - 100
eec.p3.df$relr0.lo = (eec.p3.df$r0 - eec.p3.df$r0.sd) / r0.In()[3] * 100 - 100

eec.p3.df$relr0.med = eec.p3.df$r0.med / r0.In()[3] * 100 - 100
eec.p3.df$relr0.25 = eec.p3.df$r0.25 / r0.In()[3] * 100 - 100
eec.p3.df$relr0.75 = eec.p3.df$r0.75 / r0.In()[3] * 100 - 100

save(eec.p3.df, file = paste('Review_models/r0_EECs/eec.p3.df', today, '.RData', sep = ''))

#10% EEC values ######################
eec0.1.p3.df = data.frame(chem = c('Chlorpyrifos', 'Chlorpyrifos', 'Dimethoate', 'Esfenvalerate', 'Lamda-Cyhalothrin',
                                   'Malathion', 'Permethrin', 'Profenofos', 'Terbufos', 'Zinc', 'Chlorpyrifos', 'Zinc'),
                          study = c('Satapornvanit et al 2009', 'Halstead et al 2015', 'Satapornvanit et al 2009', 
                                    rep('Halstead et al 2015', 4), 'Satapornvanit et al 2009', 'Halstead et al 2015',
                                    rep('Satapornvanit et al 2009', 3)),
                          Species = c('Macrobrachium rosenbergii', 'Procambarus clarkii', 'Macrobrachium rosenbergii', 
                                      rep('Procambarus clarkii', 4), 'Macrobrachium rosenbergii', 'Procambarus clarkii',
                                      rep('Macrobrachium rosenbergii', 3)),
                          Parameter = c(rep('muP', 10), rep('psi', 2)),
                          r0 = colMeans(r0.0.1eec.p3),
                          r0.sd = apply(r0.0.1eec.p3, 2, sd),
                          r0.med = apply(r0.0.1eec.p3, 2, median),
                          r0.25 = apply(r0.0.1eec.p3, 2, quantile, prob = 0.25),
                          r0.75 = apply(r0.0.1eec.p3, 2, quantile, prob = 0.75),
                          par.mean = colMeans(par.0.1eec.all[,37:48]))

eec0.1.p3.df$r0.up = eec0.1.p3.df$r0 + eec0.1.p3.df$r0.sd
eec0.1.p3.df$r0.lo = eec0.1.p3.df$r0 - eec0.1.p3.df$r0.sd

eec0.1.p3.df$deltar0 = eec0.1.p3.df$r0 - r0.In()[3]
eec0.1.p3.df$deltar0.up = (eec0.1.p3.df$r0 + eec0.1.p3.df$r0.sd) - r0.In()[3]
eec0.1.p3.df$deltar0.lo = (eec0.1.p3.df$r0 - eec0.1.p3.df$r0.sd) - r0.In()[3]

eec0.1.p3.df$relr0 = eec0.1.p3.df$r0 / r0.In()[3] * 100 - 100
eec0.1.p3.df$relr0.up = (eec0.1.p3.df$r0 + eec0.1.p3.df$r0.sd) / r0.In()[3] * 100 - 100
eec0.1.p3.df$relr0.lo = (eec0.1.p3.df$r0 - eec0.1.p3.df$r0.sd) / r0.In()[3] * 100 - 100

eec0.1.p3.df$relr0.med = eec0.1.p3.df$r0.med / r0.In()[3] * 100 - 100
eec0.1.p3.df$relr0.25 = eec0.1.p3.df$r0.25 / r0.In()[3] * 100 - 100
eec0.1.p3.df$relr0.75 = eec0.1.p3.df$r0.75 / r0.In()[3] * 100 - 100

save(eec0.1.p3.df, file = paste('Review_models/r0_EECs/eec0.1.p3.df', today, '.RData', sep = ''))
