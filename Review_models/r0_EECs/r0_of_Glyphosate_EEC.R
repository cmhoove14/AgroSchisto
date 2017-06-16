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
source('Response_Fxs/Ghaffar2016_cercariae.R')
source('Response_Fxs/Ghaffar2016_miracidia.R')
source('Response_Fxs/Ghaffar2016_snails.R')
source('Response_Fxs/bakry2012.R')
source('Response_Fxs/Omran&Salama_snails.R')

source('Review_models/r0_of_q.R')

library(parallel)
library(reshape2)

keep.fin.gly = c(keep.gaf.gly.cer, keep.gaf.gly.mir, keep.gaf.gly, keep.gaf.gly.fn,
                 keep.bak.gly, keep.ons.gly,
                 'r0.He', 'r0.fix', 'parameters', 'nil0', 'nil1', 'keep.fin.gly')

rm(list = setdiff(ls(), keep.fin.gly))
dev.off()

no.cores = detectCores() - 1

#Run simulations of atrazine concentrations from 0 - 2000, start with individual functions then combine ################
eec.gly = 3700  #Concentration range to test
nsims = 1000       #Number of simulations to run
  conc.gly = rep(eec.gly, nsims)
parfx = c(piM.ghaf_gly.exp_unc, piC.ghaf_gly.exp_unc, mu_Nq_gly_gaf16_uncertainty, #Functions corresponding to affected parameters
          fN.gly.fx.uncertainty, muNq_gly_Bakry12_uncertainty, fN.gly.bak.uncertainty, ons.munq.gly)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.gly, 'uniroot.all', 'rdrm', 'LL.2', 'LL.4'))

r0.gly.eec = matrix(data = NA, nrow = nsims, ncol = length(parfx)+2)
par.gly.eec = matrix(data = NA, nrow = nsims, ncol = length(parfx))

set.seed(0)

  #Fill r0 estimates ######
    #individual parameters
    r0.gly.eec[, 1] = parSapply(clus1, conc.gly, r0.He, f.pi_Mq = piM.ghaf_gly.exp_unc)[3,]
    r0.gly.eec[, 2] = parSapply(clus1, conc.gly, r0.He, f.pi_Cq = piC.ghaf_gly.exp_unc)[3,]
    r0.gly.eec[, 3] = parSapply(clus1, conc.gly, r0.He, f.mu_Nq = mu_Nq_gly_gaf16_uncertainty)[3,]
    r0.gly.eec[, 4] = parSapply(clus1, conc.gly, r0.He, f.f_Nq = fN.gly.fx.uncertainty)[3,]
    r0.gly.eec[, 5] = parSapply(clus1, conc.gly, r0.He, f.mu_Nq = muNq_gly_Bakry12_uncertainty)[3,]
    r0.gly.eec[, 6] = parSapply(clus1, conc.gly, r0.He, f.f_Nq = fN.gly.bak.uncertainty)[3,]
    r0.gly.eec[, 7] = parSapply(clus1, conc.gly, r0.He, f.mu_Nq = ons.munq.gly)[3,]
  
    r0.gly.eec[, 8] = parSapply(clus1, conc.gly, r0.He,
                                f.pi_Mq = piM.ghaf_gly.exp_unc,
                                f.pi_Cq = piC.ghaf_gly.exp_unc,
                                f.mu_Nq = mu_Nq_gly_gaf16_uncertainty,
                                f.f_Nq = fN.gly.fx.uncertainty)[3,]
    r0.gly.eec[, 9] = parSapply(clus1, conc.gly, r0.He,
                                f.pi_Mq = piM.ghaf_gly.exp_unc,
                                f.pi_Cq = piC.ghaf_gly.exp_unc,
                                f.mu_Nq = mu_Nq_gly_gaf16_uncertainty)[3,]

  #Store parameter values   #######
    par.gly.eec[, 1] = parSapply(clus1, conc.gly, piM.ghaf_gly.exp_unc)
    par.gly.eec[, 2] = parSapply(clus1, conc.gly, piC.ghaf_gly.exp_unc)
    par.gly.eec[, 3] = parSapply(clus1, conc.gly, mu_Nq_gly_gaf16_uncertainty)
    par.gly.eec[, 4] = parSapply(clus1, conc.gly, fN.gly.fx.uncertainty)
    par.gly.eec[, 5] = parSapply(clus1, conc.gly, muNq_gly_Bakry12_uncertainty)
    par.gly.eec[, 6] = parSapply(clus1, conc.gly, fN.gly.bak.uncertainty)
    par.gly.eec[, 7] = parSapply(clus1, conc.gly, ons.munq.gly)

stopCluster(clus1)
 
########## Post process ############
gly.eec.df = data.frame(chem = rep('Glyphosate', length(parfx)+2),
                        study = c(rep('Abdel-Ghaffar et al, 2016',4), rep('Bakry et al, 2012',2),
                                  'Omran & Salama 2013', 'Combined', 'Combined2'),
                        Species = c(rep('Schistosoma mansoni',2), rep('Biomphalaria alexandrina',5), NA, NA),
                        Parameter = c('piM', 'piC', 'muN', 'fN', 'muN', 'fN', 'muN', 'Combined', 'Combined2'),
                        r0 = c(mean(r0.gly.eec[, 1]), mean(r0.gly.eec[, 2]), mean(r0.gly.eec[, 3]),
                               mean(r0.gly.eec[, 4]), mean(r0.gly.eec[, 5]), mean(r0.gly.eec[, 6]),
                               mean(r0.gly.eec[, 7]), mean(r0.gly.eec[, 8]), mean(r0.gly.eec[, 9])),
                        r0.sd = c(sd(r0.gly.eec[, 1]), sd(r0.gly.eec[, 2]), sd(r0.gly.eec[, 3]),
                                  sd(r0.gly.eec[, 4]), sd(r0.gly.eec[, 5]), sd(r0.gly.eec[, 6]),
                                  sd(r0.gly.eec[, 7]), sd(r0.gly.eec[, 8]), sd(r0.gly.eec[, 9])),
                        par.mean = c(mean(par.gly.eec[, 1]), mean(par.gly.eec[, 2]), mean(par.gly.eec[, 3]),
                                     mean(par.gly.eec[, 4]), mean(par.gly.eec[, 5]), mean(par.gly.eec[, 6]),
                                     mean(par.gly.eec[, 7]), 0, 0))

gly.eec.df$r0.up = gly.eec.df$r0 + gly.eec.df$r0.sd
gly.eec.df$r0.lo = gly.eec.df$r0 - gly.eec.df$r0.sd

save(gly.eec.df, file = 'Review_models/r0_EECs/gly.eec.df.RData')