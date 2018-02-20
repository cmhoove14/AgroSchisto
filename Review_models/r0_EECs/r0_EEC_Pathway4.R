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

#Run simulations of malathion concentrations, start with individual functions then combine ################
clus1 = makeCluster(no.cores)
clusterExport(clus1, c('uniroot.all', 'rdrm', 'LL.2', 'L.4','r0.He', 'r0.Fe', 'r0.fix', 'parameters', 'nil0', 'nil1', 'par.eec.all', 'par.0.1eec.all', 'today', 'kmat'))

r0.eec.p4 = matrix(data = NA, nrow = nrow(par.eec.all), ncol = 24)

r0.0.1eec.p4 = matrix(data = NA, nrow = nrow(par.eec.all), ncol = 24)

clusterSetRNGStream(cl = clus1, iseed = 43093) #set cluster seed for reproducibility

#Fill r0 estimates for EEC #####################
  #individual parameters
    r0.eec.p4[, 1] = parSapply(clus1, par.eec.all[,49], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p4[, 2] = parSapply(clus1, par.eec.all[,50], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p4[, 3] = parSapply(clus1, par.eec.all[,51], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p4[, 4] = parSapply(clus1, par.eec.all[,52], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p4[, 5] = parSapply(clus1, par.eec.all[,53], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p4[, 6] = parSapply(clus1, par.eec.all[,54], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p4[, 7] = parSapply(clus1, par.eec.all[,55], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p4[, 8] = parSapply(clus1, par.eec.all[,56], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p4[, 9] = parSapply(clus1, par.eec.all[,57], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.eec.p4[, 10] = parSapply(clus1, par.eec.all[,58], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                #pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p4[, 11] = parSapply(clus1, par.eec.all[,59], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                #pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p4[, 12] = parSapply(clus1, par.eec.all[,60], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                #pi_Cqx = 1, 
                                v_qx = 1)[3,]
    
    r0.eec.p4[, 13] = parSapply(clus1, par.eec.all[,61], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p4[, 14] = parSapply(clus1, par.eec.all[,62], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p4[, 15] = parSapply(clus1, par.eec.all[,63], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p4[, 16] = parSapply(clus1, par.eec.all[,64], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p4[, 17] = parSapply(clus1, par.eec.all[,65], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p4[, 18] = parSapply(clus1, par.eec.all[,66], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p4[, 19] = parSapply(clus1, par.eec.all[,67], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p4[, 20] = parSapply(clus1, par.eec.all[,68], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p4[, 21] = parSapply(clus1, par.eec.all[,69], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.eec.p4[, 22] = parSapply(clus1, par.eec.all[,70], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    
    r0.eec.p4[, 23] = parSapply(clus1, par.eec.all[,71], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1)[3,] 
                                #v_qx = 1
                                #
    r0.eec.p4[, 24] = parSapply(clus1, par.eec.all[,72], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1)[3,] 
                                #v_qx = 1
    
#Fill r0 estimates for 10% EEC #####################
  #individual parameters
    r0.0.1eec.p4[, 1] = parSapply(clus1, par.eec.all[,49], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p4[, 2] = parSapply(clus1, par.eec.all[,50], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p4[, 3] = parSapply(clus1, par.eec.all[,51], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p4[, 4] = parSapply(clus1, par.eec.all[,52], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p4[, 5] = parSapply(clus1, par.eec.all[,53], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p4[, 6] = parSapply(clus1, par.eec.all[,54], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p4[, 7] = parSapply(clus1, par.eec.all[,55], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p4[, 8] = parSapply(clus1, par.eec.all[,56], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p4[, 9] = parSapply(clus1, par.eec.all[,57], FUN = r0.fix, 
                               f_Nqx = 1, 
                               mu_Pqx = 0,
                               phi_Nqx = 1, 
                               mu_Nqx = 0, 
                               alpha_qx = 1,
                               theta_qx = 1, 
                               pi_Mqx = 1, 
                               #pi_Cqx = 1, 
                               v_qx = 1)[3,]
    r0.0.1eec.p4[, 10] = parSapply(clus1, par.eec.all[,58], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                #pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p4[, 11] = parSapply(clus1, par.eec.all[,59], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                #pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p4[, 12] = parSapply(clus1, par.eec.all[,60], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                #pi_Cqx = 1, 
                                v_qx = 1)[3,]
    
    r0.0.1eec.p4[, 13] = parSapply(clus1, par.eec.all[,61], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p4[, 14] = parSapply(clus1, par.eec.all[,62], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p4[, 15] = parSapply(clus1, par.eec.all[,63], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p4[, 16] = parSapply(clus1, par.eec.all[,64], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p4[, 17] = parSapply(clus1, par.eec.all[,65], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p4[, 18] = parSapply(clus1, par.eec.all[,66], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p4[, 19] = parSapply(clus1, par.eec.all[,67], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p4[, 20] = parSapply(clus1, par.eec.all[,68], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p4[, 21] = parSapply(clus1, par.eec.all[,69], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    r0.0.1eec.p4[, 22] = parSapply(clus1, par.eec.all[,70], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                #pi_Mqx = 1, 
                                pi_Cqx = 1, 
                                v_qx = 1)[3,]
    
    r0.0.1eec.p4[, 23] = parSapply(clus1, par.eec.all[,71], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1)[3,] 
                                #v_qx = 1
    
    r0.0.1eec.p4[, 24] = parSapply(clus1, par.eec.all[,72], FUN = r0.fix, 
                                f_Nqx = 1, 
                                mu_Pqx = 0,
                                phi_Nqx = 1, 
                                mu_Nqx = 0, 
                                alpha_qx = 1,
                                theta_qx = 1, 
                                pi_Mqx = 1, 
                                pi_Cqx = 1)[3,] 
                                #v_qx = 1
    
stopCluster(clus1) 
#Post process ############ 
#EEC runs ################
eec.p4.df = data.frame(chem = c(rep('Atrazine',4), 'Butachlor', 'Fluazifop-p-butyl', 'Butralin', 'Glyphosate', 'Pendimethalin',
                                'Malathion', 'Chlorpyrifos', 'Profenofos', 'Chlorpyrifos', 'Profenofos', 'Butralin', 'Glyphosate', 
                                'Pendimethalin', 'Malathion', 'Butachlor', 'Fluazifop-p-butyl', 'Ammonium Sulphate', 'Urea',
                                'Ammonium Sulphate', 'Urea'),
                       study = c('Griggs et al 2008', 'Koprivnikar et al 2006', 'Rohr et al 2008', 'Meta',
                                 rep('Tantawy 2002', 2), rep('Abdel-Ghaffar et al 2016', 3), 'Tchounwou et al 1992',
                                 rep('Hasheesh & Mohamed 2011', 4), rep('Abdel-Ghaffar et al 2016', 3), 'Tchounwou et al 1991a',
                                 rep('Tantawy 2002', 2), rep('Tchounwou et al 1991b', 4)),
                       Species = c(rep('Echinistoma trivolvis', 4), rep('Schistosoma mansoni', 6), rep('Schistosoma haemotobium', 4),
                                   rep('Schistosoma mansoni', 10)),
                       Parameter = c(rep('pi_C', 12), rep('pi_M', 10), rep('v', 2)),
                       r0 = colMeans(r0.eec.p4),
                       r0.sd = apply(r0.eec.p4, 2, sd),
                       r0.med = apply(r0.eec.p4, 2, median),
                       r0.25 = apply(r0.eec.p4, 2, quantile, prob = 0.25),
                       r0.75 = apply(r0.eec.p4, 2, quantile, prob = 0.75),
                       par.mean = colMeans(par.eec.all[,49:72]))

eec.p4.df$r0.up = eec.p4.df$r0 + eec.p4.df$r0.sd
eec.p4.df$r0.lo = eec.p4.df$r0 - eec.p4.df$r0.sd

eec.p4.df$deltar0 = eec.p4.df$r0 - r0.He()[3]
eec.p4.df$deltar0.up = (eec.p4.df$r0 + eec.p4.df$r0.sd) - r0.He()[3]
eec.p4.df$deltar0.lo = (eec.p4.df$r0 - eec.p4.df$r0.sd) - r0.He()[3]

eec.p4.df$relr0 = eec.p4.df$r0 / r0.He()[3] * 100 - 100
eec.p4.df$relr0.up = (eec.p4.df$r0 + eec.p4.df$r0.sd) / r0.He()[3] * 100 - 100
eec.p4.df$relr0.lo = (eec.p4.df$r0 - eec.p4.df$r0.sd) / r0.He()[3] * 100 - 100

eec.p4.df$relr0.med = eec.p4.df$r0.med / r0.He()[3] * 100 - 100
eec.p4.df$relr0.25 = eec.p4.df$r0.25 / r0.He()[3] * 100 - 100
eec.p4.df$relr0.75 = eec.p4.df$r0.75 / r0.He()[3] * 100 - 100

save(eec.p4.df, file = paste('Review_models/r0_EECs/eec.p4.df', today, '.RData', sep = ''))

#10% EEC values ######################
eec0.1.p4.df = data.frame(chem = c(rep('Atrazine',4), 'Butachlor', 'Fluazifop-p-butyl', 'Butralin', 'Glyphosate', 'Pendimethalin',
                                   'Malathion', 'Chlorpyrifos', 'Profenofos', 'Chlorpyrifos', 'Profenofos', 'Butralin', 'Glyphosate', 
                                   'Pendimethalin', 'Malathion', 'Butachlor', 'Fluazifop-p-butyl', 'Ammonium Sulphate', 'Urea',
                                   'Ammonium Sulphate', 'Urea'),
                          study = c('Griggs et al 2008', 'Koprivnikar et al 2006', 'Rohr et al 2008', 'Meta',
                                    rep('Tantawy 2002', 2), rep('Abdel-Ghaffar et al 2016', 3), 'Tchounwou et al 1992',
                                    rep('Hasheesh & Mohamed 2011', 4), rep('Abdel-Ghaffar et al 2016', 3), 'Tchounwou et al 1991a',
                                    rep('Tantawy 2002', 2), rep('Tchounwou et al 1991b', 4)),
                          Species = c(rep('Echinistoma trivolvis', 4), rep('Schistosoma mansoni', 6), rep('Schistosoma haemotobium', 4),
                                      rep('Schistosoma mansoni', 10)),
                          Parameter = c(rep('pi_C', 12), rep('pi_M', 10), rep('v', 2)),
                          r0 = colMeans(r0.0.1eec.p4),
                          r0.sd = apply(r0.0.1eec.p4, 2, sd),
                          r0.med = apply(r0.0.1eec.p4, 2, median),
                          r0.25 = apply(r0.0.1eec.p4, 2, quantile, prob = 0.25),
                          r0.75 = apply(r0.0.1eec.p4, 2, quantile, prob = 0.75),
                          par.mean = colMeans(par.0.1eec.all[,49:72]))

eec0.1.p4.df$r0.up = eec0.1.p4.df$r0 + eec0.1.p4.df$r0.sd
eec0.1.p4.df$r0.lo = eec0.1.p4.df$r0 - eec0.1.p4.df$r0.sd

eec0.1.p4.df$deltar0 = eec0.1.p4.df$r0 - r0.He()[3]
eec0.1.p4.df$deltar0.up = (eec0.1.p4.df$r0 + eec0.1.p4.df$r0.sd) - r0.He()[3]
eec0.1.p4.df$deltar0.lo = (eec0.1.p4.df$r0 - eec0.1.p4.df$r0.sd) - r0.He()[3]

eec0.1.p4.df$relr0 = eec0.1.p4.df$r0 / r0.He()[3] * 100 - 100
eec0.1.p4.df$relr0.up = (eec0.1.p4.df$r0 + eec0.1.p4.df$r0.sd) / r0.He()[3] * 100 - 100
eec0.1.p4.df$relr0.lo = (eec0.1.p4.df$r0 - eec0.1.p4.df$r0.sd) / r0.He()[3] * 100 - 100

eec0.1.p4.df$relr0.med = eec0.1.p4.df$r0.med / r0.He()[3] * 100 - 100
eec0.1.p4.df$relr0.25 = eec0.1.p4.df$r0.25 / r0.He()[3] * 100 - 100
eec0.1.p4.df$relr0.75 = eec0.1.p4.df$r0.75 / r0.He()[3] * 100 - 100

save(eec0.1.p4.df, file = paste('Review_models/r0_EECs/eec0.1.p4.df', today, '.RData', sep = ''))
