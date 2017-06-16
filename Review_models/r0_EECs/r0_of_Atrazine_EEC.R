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
source('Response_Fxs/bakry2012.R')
source('Response_Fxs/piC_atr_meta_beq.R')
source('Response_Fxs/koprivnikar06_piC_beq.R')
source('Response_Fxs/griggs08_piC_beq.R')
source('Response_Fxs/Baxter_Rohr_Atrazine2011.R')
source('Response_Fxs/rohr08_piC_atrOnly.R')
source('Response_Fxs/Omran&Salama_snails.R')

source('Review_models/r0_of_q.R')

library(parallel)
library(reshape2)

keep.fin.atr = c(keep.bak.atr, keep.grg08, keep.kop06.beq, keep.meta.piC, keep.baxrohr, keep.atr.rohr08, keep.ons.atr,
                 'r0.He', 'r0.fix', 'parameters', 'nil0', 'nil1', 'keep.fin.atr')

rm(list = setdiff(ls(), keep.fin.atr))
dev.off()

no.cores = detectCores() - 1

#Run simulations of atrazine concentrations from 0 - 2000, start with individual functions then combine ################
eec.atr = 102
nsims = 1000       #Number of simulations to run
  conc.atr = rep(eec.atr, nsims)
parfx = c(piC.meta_atr_unc, muNq_atr_Bakry12_uncertainty, phi_Nq_atr_baxrohr, #Functions corresponding to affected parameters
          piC.grg08_atr_unc, piC.kop_atr_unc2, piC.atr.rohr08.lin, ons.munq.atr)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.atr, 'uniroot.all', 'rdrm', 'LL.2'))

r0.atr.eec = matrix(data = NA, nrow = nsims, ncol = length(parfx)+1)
par.atr.eec = matrix(data = NA, nrow = nsims, ncol = length(parfx))

set.seed(0)

   #Fill r0 estimates ######
    #individual parameters
    r0.atr.eec[, 1] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.meta_atr_unc)[3,]
    r0.atr.eec[, 2] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.grg08_atr_unc)[3,]
    r0.atr.eec[, 3] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.kop_atr_unc2)[3,]
    r0.atr.eec[, 4] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.atr.rohr08.lin)[3,]
    r0.atr.eec[, 5] = parSapply(clus1, conc.atr, r0.He, f.mu_Nq = muNq_atr_Bakry12_uncertainty)[3,]
    r0.atr.eec[, 6] = parSapply(clus1, conc.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.atr.eec[, 7] = parSapply(clus1, conc.atr, r0.He, f.mu_Nq = ons.munq.atr)[3,]
  
    r0.atr.eec[, 8] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = ons.munq.atr,
                                     f.pi_Cq = piC.atr.rohr08.lin,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
  #Store parameter values   #######
    par.atr.eec[, 1] = parSapply(clus1, conc.atr, piC.meta_atr_unc)
    par.atr.eec[, 2] = parSapply(clus1, conc.atr, piC.grg08_atr_unc)
    par.atr.eec[, 3] = parSapply(clus1, conc.atr, piC.kop_atr_unc2)
    par.atr.eec[, 4] = parSapply(clus1, conc.atr, piC.atr.rohr08.lin)
    par.atr.eec[, 5] = parSapply(clus1, conc.atr, muNq_atr_Bakry12_uncertainty)
    par.atr.eec[, 6] = parSapply(clus1, conc.atr, phi_Nq_atr_baxrohr.no30)
    par.atr.eec[, 7] = parSapply(clus1, conc.atr, ons.munq.atr)
    

stopCluster(clus1)
 
########## Post process ############
atr.eec.df = data.frame(chem = rep('Atrazine', length(parfx)+1),
                        study = c('piC meta', 'Griggs & Belden, 2008', 'Koprivnikar et al, 2007',
                                  'Rohr et al, 2008', 'Bakry et al, 2012', 'Rohr et al, 2012',
                                  'Omran & Salama, 2013', 'Combined'),
                        Species = c(rep('Echinistoma trivolvis', 4), 'Biomphalaria alexandrina',
                                    'Physella spp.', 'Biomphalaria alexandrina', NA), 
                        Parameter = c('piC', 'piC', 'piC', 'piC',
                                      'muN', 'phiN', 'muN', 'Combined'),
                        r0 = c(mean(r0.atr.eec[, 1]), mean(r0.atr.eec[, 2]), mean(r0.atr.eec[, 3]),
                               mean(r0.atr.eec[, 4]), mean(r0.atr.eec[, 5]), mean(r0.atr.eec[, 6]),
                               mean(r0.atr.eec[, 7]), mean(r0.atr.eec[, 8])),
                        r0.sd = c(sd(r0.atr.eec[, 1]), sd(r0.atr.eec[, 2]), sd(r0.atr.eec[, 3]),
                                  sd(r0.atr.eec[, 4]), sd(r0.atr.eec[, 5]), sd(r0.atr.eec[, 6]),
                                  sd(r0.atr.eec[, 7]), sd(r0.atr.eec[, 8])),
                        par.mean = c(mean(par.atr.eec[, 1]), mean(par.atr.eec[, 2]), mean(par.atr.eec[, 3]),
                                     mean(par.atr.eec[, 4]), mean(par.atr.eec[, 5]), mean(par.atr.eec[, 6]),
                                     mean(par.atr.eec[, 7]), 0))

atr.eec.df$r0.up = atr.eec.df$r0 + atr.eec.df$r0.sd
atr.eec.df$r0.lo = atr.eec.df$r0 - atr.eec.df$r0.sd

save(atr.eec.df, file = 'Review_models/r0_EECs/atr.eec.df.RData')