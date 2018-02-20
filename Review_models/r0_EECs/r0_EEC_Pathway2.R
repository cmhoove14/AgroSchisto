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

#Run simulations of pf pathway 2 functional responses ################
clus1 = makeCluster(no.cores)
clusterExport(clus1, c('uniroot.all', 'rdrm', 'LL.2', 'L.4','r0.He', 'r0.Fe', 'r0.fix', 'parameters', 'nil0', 'nil1', 'par.eec.all', 'par.0.1eec.all', 'today', 'kmat'))

r0.eec.p2 = matrix(data = NA, nrow = nrow(par.eec.all), ncol = 29)

r0.0.1eec.p2 = matrix(data = NA, nrow = nrow(par.eec.all), ncol = 29)

clusterSetRNGStream(cl = clus1, iseed = 43093) #set cluster seed for reproducibility

#Fill r0 estimates for EEC #####################
  #individual parameters
    r0.eec.p2[, 1] = parSapply(clus1, par.eec.all[,8], FUN = r0.fix, 
                               #f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p2[, 2] = parSapply(clus1, par.eec.all[,9], FUN = r0.fix, 
                               #f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p2[, 3] = parSapply(clus1, par.eec.all[,10], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               #mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p2[, 4] = parSapply(clus1, par.eec.all[,11], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               #mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p2[, 5] = parSapply(clus1, par.eec.all[,12], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               #mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p2[, 6] = parSapply(clus1, par.eec.all[,13], FUN = r0.fix, 
                               #f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p2[, 7] = parSapply(clus1, par.eec.all[,14], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               #mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p2[, 8] = parSapply(clus1, par.eec.all[,15], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               #mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p2[, 9] = parSapply(clus1, par.eec.all[,16], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               #mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p2[, 10] = parSapply(clus1, par.eec.all[,17], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 11] = parSapply(clus1, par.eec.all[,18], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 12] = parSapply(clus1, par.eec.all[,19], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 13] = parSapply(clus1, par.eec.all[,20], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 14] = parSapply(clus1, par.eec.all[,21], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 15] = parSapply(clus1, par.eec.all[,22], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 16] = parSapply(clus1, par.eec.all[,23], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 17] = parSapply(clus1, par.eec.all[,24], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 18] = parSapply(clus1, par.eec.all[,25], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 19] = parSapply(clus1, par.eec.all[,26], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 20] = parSapply(clus1, par.eec.all[,27], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 21] = parSapply(clus1, par.eec.all[,28], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 22] = parSapply(clus1, par.eec.all[,29], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 23] = parSapply(clus1, par.eec.all[,30], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 24] = parSapply(clus1, par.eec.all[,31], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 25] = parSapply(clus1, par.eec.all[,32], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 26] = parSapply(clus1, par.eec.all[,33], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 27] = parSapply(clus1, par.eec.all[,34], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 28] = parSapply(clus1, par.eec.all[,35], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p2[, 29] = parSapply(clus1, par.eec.all[,36], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    
#Fill r0 estimates for 10% EEC #####################
  #individual parameters
    r0.0.1eec.p2[, 1] = parSapply(clus1, par.0.1eec.all[,8], FUN = r0.fix, 
                               #f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p2[, 2] = parSapply(clus1, par.0.1eec.all[,9], FUN = r0.fix, 
                               #f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p2[, 3] = parSapply(clus1, par.0.1eec.all[,10], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               #mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p2[, 4] = parSapply(clus1, par.0.1eec.all[,11], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               #mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p2[, 5] = parSapply(clus1, par.0.1eec.all[,12], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               #mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p2[, 6] = parSapply(clus1, par.0.1eec.all[,13], FUN = r0.fix, 
                               #f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p2[, 7] = parSapply(clus1, par.0.1eec.all[,14], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               #mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p2[, 8] = parSapply(clus1, par.0.1eec.all[,15], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               #mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p2[, 9] = parSapply(clus1, par.0.1eec.all[,16], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               #mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p2[, 10] = parSapply(clus1, par.0.1eec.all[,17], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 11] = parSapply(clus1, par.0.1eec.all[,18], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 12] = parSapply(clus1, par.0.1eec.all[,19], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 13] = parSapply(clus1, par.0.1eec.all[,20], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 14] = parSapply(clus1, par.0.1eec.all[,21], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 15] = parSapply(clus1, par.0.1eec.all[,22], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 16] = parSapply(clus1, par.0.1eec.all[,23], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 17] = parSapply(clus1, par.0.1eec.all[,24], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 18] = parSapply(clus1, par.0.1eec.all[,25], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 19] = parSapply(clus1, par.0.1eec.all[,26], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 20] = parSapply(clus1, par.0.1eec.all[,27], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 21] = parSapply(clus1, par.0.1eec.all[,28], FUN = r0.fix, 
                                #f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 22] = parSapply(clus1, par.0.1eec.all[,29], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 23] = parSapply(clus1, par.0.1eec.all[,30], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 24] = parSapply(clus1, par.0.1eec.all[,31], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 25] = parSapply(clus1, par.0.1eec.all[,32], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 26] = parSapply(clus1, par.0.1eec.all[,33], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 27] = parSapply(clus1, par.0.1eec.all[,34], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 28] = parSapply(clus1, par.0.1eec.all[,35], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p2[, 29] = parSapply(clus1, par.0.1eec.all[,36], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                #mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]    
    
  
stopCluster(clus1) 
#Post process ############ 
#EEC runs ################
eec.p2.df = data.frame(chem = c('Atrazine', 'Glyphosate', 'Glyphosate', 'Atrazine', 'Malathion', 'Malathion',
                                'Butralin', 'Glyphosate', 'Pendimethalin', 'Butralin', 'Glyphosate', 'Pendimethalin', 
                                'Malathion', 'Deltamethrin', 'Malathion', 'Deltamethrin', 'Chlorpyrifos', 'Profenofos',
                                'Chlorpyrifos', 'Profenofos', 'Chlorpyrifos', 'Chlorpyrifos', 'Atrazine', 'Glyphosate',
                                'Butachlor', 'Fluazifop-p-butyl', 'Urea', 'Potassium Sulphate', 'Ammonium Nitrate'),
                       study = c(rep('Bakry et al 2012', 4), rep('Tchounwou et al 1991', 2), rep('Abdel-Ghaffar et al 2016', 6),
                                 rep('Bakry et al 2011', 4), rep('Hasheesh & Mohamed 2011', 4), rep('Ibrahim et al 1992', 2),
                                 rep('Omran & Salama 2013', 2), rep('Tantawy 2002', 2), rep('Ragab & Shoukry 2006', 3)),
                       Species = c(rep('Biomphalaria alexandrina', 4), rep('Bulinus havenensis',2), rep('Biomphalaria alexandrina', 6),
                                   rep('Helisoma duryi', 4), rep('Bulinus truncatus', 4), rep('Biomphalaria alexandrina', 9)),
                       Parameter = c(rep('fN', 2), rep('muN', 3), 'fN', rep('muN', 3), rep('fN', 3), rep('muN', 2),
                                     rep('fN', 4), rep('muN', 2), 'fN', rep('muN', 8)),
                       r0 = colMeans(r0.eec.p2),
                       r0.sd = apply(r0.eec.p2, 2, sd),
                       r0.med = apply(r0.eec.p2, 2, median),
                       r0.25 = apply(r0.eec.p2, 2, quantile, prob = 0.25),
                       r0.75 = apply(r0.eec.p2, 2, quantile, prob = 0.75),
                       par.mean = colMeans(par.eec.all[,8:36]))

eec.p2.df$r0.up = eec.p2.df$r0 + eec.p2.df$r0.sd
eec.p2.df$r0.lo = eec.p2.df$r0 - eec.p2.df$r0.sd

eec.p2.df$deltar0 = eec.p2.df$r0 - r0.He()[3]
eec.p2.df$deltar0.up = (eec.p2.df$r0 + eec.p2.df$r0.sd) - r0.He()[3]
eec.p2.df$deltar0.lo = (eec.p2.df$r0 - eec.p2.df$r0.sd) - r0.He()[3]

eec.p2.df$relr0 = eec.p2.df$r0 / r0.He()[3] * 100 - 100
eec.p2.df$relr0.up = (eec.p2.df$r0 + eec.p2.df$r0.sd) / r0.He()[3] * 100 - 100
eec.p2.df$relr0.lo = (eec.p2.df$r0 - eec.p2.df$r0.sd) / r0.He()[3] * 100 - 100

eec.p2.df$relr0.med = eec.p2.df$r0.med / r0.He()[3] * 100 - 100
eec.p2.df$relr0.25 = eec.p2.df$r0.25 / r0.He()[3] * 100 - 100
eec.p2.df$relr0.75 = eec.p2.df$r0.75 / r0.He()[3] * 100 - 100

save(eec.p2.df, file = paste('Review_models/r0_EECs/eec.p2.df', today, '.RData', sep = ''))

#10% EEC values ######################
eec0.1.p2.df = data.frame(chem = c('Atrazine', 'Glyphosate', 'Glyphosate', 'Atrazine', 'Malathion', 'Malathion',
                                   'Butralin', 'Glyphosate', 'Pendimethalin', 'Butralin', 'Glyphosate', 'Pendimethalin', 
                                   'Malathion', 'Deltamethrin', 'Malathion', 'Deltamethrin', 'Chlorpyrifos', 'Profenofos',
                                   'Chlorpyrifos', 'Profenofos', 'Chlorpyrifos', 'Chlorpyrifos', 'Atrazine', 'Glyphosate',
                                   'Butachlor', 'Fluazifop-p-butyl', 'Urea', 'Potassium Sulphate', 'Ammonium Nitrate'),
                          study = c(rep('Bakry et al 2012', 4), rep('Tchounwou et al 1991', 2), rep('Abdel-Ghaffar et al 2016', 6),
                                    rep('Bakry et al 2011', 4), rep('Hasheesh & Mohamed 2011', 4), rep('Ibrahim et al 1992', 2),
                                    rep('Omran & Salama 2013', 2), rep('Tantawy 2002', 2), rep('Ragab & Shoukry 2006', 3)),
                          Species = c(rep('Biomphalaria alexandrina', 4), rep('Bulinus havenensis',2), rep('Biomphalaria alexandrina', 6),
                                      rep('Helisoma duryi', 4), rep('Bulinus truncatus', 4), rep('Biomphalaria alexandrina', 9)),
                          Parameter = c(rep('fN', 2), rep('muN', 3), 'fN', rep('muN', 3), rep('fN', 3), rep('muN', 2),
                                        rep('fN', 4), rep('muN', 2), 'fN', rep('muN', 8)),
                          r0 = colMeans(r0.0.1eec.p2),
                          r0.sd = apply(r0.0.1eec.p2, 2, sd),
                          r0.med = apply(r0.0.1eec.p2, 2, median),
                          r0.25 = apply(r0.0.1eec.p2, 2, quantile, prob = 0.25),
                          r0.75 = apply(r0.0.1eec.p2, 2, quantile, prob = 0.75),
                          par.mean = colMeans(par.0.1eec.all[,8:36]))

eec0.1.p2.df$r0.up = eec0.1.p2.df$r0 + eec0.1.p2.df$r0.sd
eec0.1.p2.df$r0.lo = eec0.1.p2.df$r0 - eec0.1.p2.df$r0.sd

eec0.1.p2.df$deltar0 = eec0.1.p2.df$r0 - r0.He()[3]
eec0.1.p2.df$deltar0.up = (eec0.1.p2.df$r0 + eec0.1.p2.df$r0.sd) - r0.He()[3]
eec0.1.p2.df$deltar0.lo = (eec0.1.p2.df$r0 - eec0.1.p2.df$r0.sd) - r0.He()[3]

eec0.1.p2.df$relr0 = eec0.1.p2.df$r0 / r0.He()[3] * 100 - 100
eec0.1.p2.df$relr0.up = (eec0.1.p2.df$r0 + eec0.1.p2.df$r0.sd) / r0.He()[3] * 100 - 100
eec0.1.p2.df$relr0.lo = (eec0.1.p2.df$r0 - eec0.1.p2.df$r0.sd) / r0.He()[3] * 100 - 100

eec0.1.p2.df$relr0.med = eec0.1.p2.df$r0.med / r0.He()[3] * 100 - 100
eec0.1.p2.df$relr0.25 = eec0.1.p2.df$r0.25 / r0.He()[3] * 100 - 100
eec0.1.p2.df$relr0.75 = eec0.1.p2.df$r0.75 / r0.He()[3] * 100 - 100

save(eec0.1.p2.df, file = paste('Review_models/r0_EECs/eec0.1.p2.df', today, '.RData', sep = ''))
