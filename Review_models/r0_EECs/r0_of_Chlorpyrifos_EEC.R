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
source('Response_Fxs/Halstead_Insecticides2015.R')
source('Response_Fxs/Hasheesh2011_larval.R')
source('Response_Fxs/Hasheesh2011_snails.R')
source('Response_Fxs/Ibrahim92_fN_muN.R')
source('Response_Fxs/Satapornvanit_Insecticides2009.R')

source('Review_models/r0_of_q.R')

library(parallel)
library(fBasics)
library(drc)
 
keep.fin.ch = c(keep.hal15.muP[c(1,7,12,18,24)+1], keep.hsh.ch.pim, keep.hsh.ch.pic,
                keep.hsh.ch, keep.ibr.ch, keep.chlor.sat09,
                'r0.In', 'r0.fix', 'parameters', 'nil0', 'nil1', 'keep.fin.ch')

rm(list = setdiff(ls(), keep.fin.ch))
dev.off()

no.cores = detectCores() - 1

#Run simulations of chlorpyrifos concentrations, start with individual functions then combine ################
ch.eec = 64  #Concentration range to test
nsims = 1000         #Number of simulations to run
  conc.ch = rep(ch.eec, nsims)
parfx = c(muPq_chlor_Halstead_uncertainty, muPq_chlor_satapornvanit09_uncertainty,
          psi_q_chlor_satapornvanit09_uncertainty,
          piC_ch_Hash11_uncertainty, piM_ch_Hash11_uncertainty,
          muNq_ch_hash11_uncertainty, mu_N_chlor_ibr92_uncertainty,
          fN.hash.chlor.uncertainty, f_N_chlor_ibr92_uncertainty)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.ch, 'uniroot.all', 'rdrm', 'LL.2', 'L.4', 'LL.3'))

r0.ch.eec = matrix(data = NA, nrow = nsims, ncol = length(parfx)+1)
par.ch.eec = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.ch.0.5eec = matrix(data = NA, nrow = nsims, ncol = length(parfx)+1)
par.ch.0.5eec = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.ch.0.1eec = matrix(data = NA, nrow = nsims, ncol = length(parfx)+1)
par.ch.0.1eec = matrix(data = NA, nrow = nsims, ncol = length(parfx))

set.seed(0)

  #Fill r0 estimates for EEC runs   ########
    #individual parameters #########
      r0.ch.eec[, 1] = parSapply(clus1, conc.ch, r0.In, f.mu_Pq = muPq_chlor_Halstead_uncertainty)[3,]
      r0.ch.eec[, 2] = parSapply(clus1, conc.ch, r0.In, f.mu_Pq = muPq_chlor_satapornvanit09_uncertainty)[3,]
      r0.ch.eec[, 3] = parSapply(clus1, conc.ch, r0.In, f.alpha_q = psi_q_chlor_satapornvanit09_uncertainty)[3,]
      r0.ch.eec[, 4] = parSapply(clus1, conc.ch, r0.In, f.pi_Cq = piC_ch_Hash11_uncertainty)[3,]
      r0.ch.eec[, 5] = parSapply(clus1, conc.ch, r0.In, f.pi_Mq = piM_ch_Hash11_uncertainty)[3,]
      r0.ch.eec[, 6] = parSapply(clus1, conc.ch, r0.In, f.mu_Nq = muNq_ch_hash11_uncertainty)[3,]
      r0.ch.eec[, 7] = parSapply(clus1, conc.ch, r0.In, f.mu_Nq = mu_N_chlor_ibr92_uncertainty)[3,]
      r0.ch.eec[, 8] = parSapply(clus1, conc.ch, r0.In, f.f_Nq = fN.hash.chlor.uncertainty)[3,]
      r0.ch.eec[, 9] = parSapply(clus1, conc.ch, r0.In, f.f_Nq = f_N_chlor_ibr92_uncertainty)[3,]
      
    #final run with chosen parameters ######
      r0.ch.eec[, 10] = parSapply(clus1, conc.ch, r0.In,
                                      f.mu_Nq = muNq_ch_hash11_uncertainty,
                                      f.f_Nq = f_N_chlor_ibr92_uncertainty,
                                      f.pi_Mq = piM_ch_Hash11_uncertainty,
                                      f.pi_Cq = piC_ch_Hash11_uncertainty,
                                      f.mu_Pq = muPq_chlor_Halstead_uncertainty)[3,]
  
  #fill parameter values at EEC#######
    par.ch.eec[, 1] = parSapply(clus1, conc.ch, muPq_chlor_Halstead_uncertainty)
    par.ch.eec[, 2] = parSapply(clus1, conc.ch, muPq_chlor_satapornvanit09_uncertainty)
    par.ch.eec[, 3] = parSapply(clus1, conc.ch, psi_q_chlor_satapornvanit09_uncertainty)
    par.ch.eec[, 4] = parSapply(clus1, conc.ch, piC_ch_Hash11_uncertainty)
    par.ch.eec[, 5] = parSapply(clus1, conc.ch, piM_ch_Hash11_uncertainty)
    par.ch.eec[, 6] = parSapply(clus1, conc.ch, muNq_ch_hash11_uncertainty)
    par.ch.eec[, 7] = parSapply(clus1, conc.ch, mu_N_chlor_ibr92_uncertainty)
    par.ch.eec[, 8] = parSapply(clus1, conc.ch, fN.hash.chlor.uncertainty)
    par.ch.eec[, 9] = parSapply(clus1, conc.ch, f_N_chlor_ibr92_uncertainty)
    
#Fill r0 estimates for 50% EEC runs   ########
  #individual parameters #########
    r0.ch.0.5eec[, 1] = parSapply(clus1, conc.ch*0.5, r0.In, f.mu_Pq = muPq_chlor_Halstead_uncertainty)[3,]
    r0.ch.0.5eec[, 2] = parSapply(clus1, conc.ch*0.5, r0.In, f.mu_Pq = muPq_chlor_satapornvanit09_uncertainty)[3,]
    r0.ch.0.5eec[, 3] = parSapply(clus1, conc.ch*0.5, r0.In, f.alpha_q = psi_q_chlor_satapornvanit09_uncertainty)[3,]
    r0.ch.0.5eec[, 4] = parSapply(clus1, conc.ch*0.5, r0.In, f.pi_Cq = piC_ch_Hash11_uncertainty)[3,]
    r0.ch.0.5eec[, 5] = parSapply(clus1, conc.ch*0.5, r0.In, f.pi_Mq = piM_ch_Hash11_uncertainty)[3,]
    r0.ch.0.5eec[, 6] = parSapply(clus1, conc.ch*0.5, r0.In, f.mu_Nq = muNq_ch_hash11_uncertainty)[3,]
    r0.ch.0.5eec[, 7] = parSapply(clus1, conc.ch*0.5, r0.In, f.mu_Nq = mu_N_chlor_ibr92_uncertainty)[3,]
    r0.ch.0.5eec[, 8] = parSapply(clus1, conc.ch*0.5, r0.In, f.f_Nq = fN.hash.chlor.uncertainty)[3,]
    r0.ch.0.5eec[, 9] = parSapply(clus1, conc.ch*0.5, r0.In, f.f_Nq = f_N_chlor_ibr92_uncertainty)[3,]

  #final run with chosen parameters ######
    r0.ch.0.5eec[, 10] = parSapply(clus1, conc.ch*0.5, r0.In,
                                   f.mu_Nq = muNq_ch_hash11_uncertainty,
                                   f.f_Nq = f_N_chlor_ibr92_uncertainty,
                                   f.pi_Mq = piM_ch_Hash11_uncertainty,
                                   f.pi_Cq = piC_ch_Hash11_uncertainty,
                                   f.mu_Pq = muPq_chlor_Halstead_uncertainty)[3,]

  #fill parameter values for 50% EEC runs #######
    par.ch.0.5eec[, 1] = parSapply(clus1, conc.ch*0.5, muPq_chlor_Halstead_uncertainty)
    par.ch.0.5eec[, 2] = parSapply(clus1, conc.ch*0.5, muPq_chlor_satapornvanit09_uncertainty)
    par.ch.0.5eec[, 3] = parSapply(clus1, conc.ch*0.5, psi_q_chlor_satapornvanit09_uncertainty)
    par.ch.0.5eec[, 4] = parSapply(clus1, conc.ch*0.5, piC_ch_Hash11_uncertainty)
    par.ch.0.5eec[, 5] = parSapply(clus1, conc.ch*0.5, piM_ch_Hash11_uncertainty)
    par.ch.0.5eec[, 6] = parSapply(clus1, conc.ch*0.5, muNq_ch_hash11_uncertainty)
    par.ch.0.5eec[, 7] = parSapply(clus1, conc.ch*0.5, mu_N_chlor_ibr92_uncertainty)
    par.ch.0.5eec[, 8] = parSapply(clus1, conc.ch*0.5, fN.hash.chlor.uncertainty)
    par.ch.0.5eec[, 9] = parSapply(clus1, conc.ch*0.5, f_N_chlor_ibr92_uncertainty)


      
#Fill r0 estimates for 10% EEC runs   ########
  #individual parameters #########
    r0.ch.0.1eec[, 1] = parSapply(clus1, conc.ch*0.1, r0.In, f.mu_Pq = muPq_chlor_Halstead_uncertainty)[3,]
    r0.ch.0.1eec[, 2] = parSapply(clus1, conc.ch*0.1, r0.In, f.mu_Pq = muPq_chlor_satapornvanit09_uncertainty)[3,]
    r0.ch.0.1eec[, 3] = parSapply(clus1, conc.ch*0.1, r0.In, f.alpha_q = psi_q_chlor_satapornvanit09_uncertainty)[3,]
    r0.ch.0.1eec[, 4] = parSapply(clus1, conc.ch*0.1, r0.In, f.pi_Cq = piC_ch_Hash11_uncertainty)[3,]
    r0.ch.0.1eec[, 5] = parSapply(clus1, conc.ch*0.1, r0.In, f.pi_Mq = piM_ch_Hash11_uncertainty)[3,]
    r0.ch.0.1eec[, 6] = parSapply(clus1, conc.ch*0.1, r0.In, f.mu_Nq = muNq_ch_hash11_uncertainty)[3,]
    r0.ch.0.1eec[, 7] = parSapply(clus1, conc.ch*0.1, r0.In, f.mu_Nq = mu_N_chlor_ibr92_uncertainty)[3,]
    r0.ch.0.1eec[, 8] = parSapply(clus1, conc.ch*0.1, r0.In, f.f_Nq = fN.hash.chlor.uncertainty)[3,]
    r0.ch.0.1eec[, 9] = parSapply(clus1, conc.ch*0.1, r0.In, f.f_Nq = f_N_chlor_ibr92_uncertainty)[3,]
  
  #final run with chosen parameters ######
    r0.ch.0.1eec[, 10] = parSapply(clus1, conc.ch*0.1, r0.In,
                                f.mu_Nq = muNq_ch_hash11_uncertainty,
                                f.f_Nq = f_N_chlor_ibr92_uncertainty,
                                f.pi_Mq = piM_ch_Hash11_uncertainty,
                                f.pi_Cq = piC_ch_Hash11_uncertainty,
                                f.mu_Pq = muPq_chlor_Halstead_uncertainty)[3,]
  
  #fill parameter values for 10% EEC runs #######
    par.ch.0.1eec[, 1] = parSapply(clus1, conc.ch*0.1, muPq_chlor_Halstead_uncertainty)
    par.ch.0.1eec[, 2] = parSapply(clus1, conc.ch*0.1, muPq_chlor_satapornvanit09_uncertainty)
    par.ch.0.1eec[, 3] = parSapply(clus1, conc.ch*0.1, psi_q_chlor_satapornvanit09_uncertainty)
    par.ch.0.1eec[, 4] = parSapply(clus1, conc.ch*0.1, piC_ch_Hash11_uncertainty)
    par.ch.0.1eec[, 5] = parSapply(clus1, conc.ch*0.1, piM_ch_Hash11_uncertainty)
    par.ch.0.1eec[, 6] = parSapply(clus1, conc.ch*0.1, muNq_ch_hash11_uncertainty)
    par.ch.0.1eec[, 7] = parSapply(clus1, conc.ch*0.1, mu_N_chlor_ibr92_uncertainty)
    par.ch.0.1eec[, 8] = parSapply(clus1, conc.ch*0.1, fN.hash.chlor.uncertainty)
    par.ch.0.1eec[, 9] = parSapply(clus1, conc.ch*0.1, f_N_chlor_ibr92_uncertainty)
  
stopCluster(clus1)
########## Post process ############ 
#EEC runs
ch.eec.df = data.frame(chem = rep('Chlorpyrifos', length(parfx)+1),
                       study = c('Halstead et al, 2015', rep('Satapornvanit et al, 2009', 2),
                                 rep('Hasheesh & Mohamed, 2011', 3), 'Ibrahim et al, 1992', 
                                 'Hasheesh & Mohamed, 2011', 'Ibrahim et al, 1992', 'Combined'),
                       Species = c('Procambarus clarkii', rep('Macrobrachium rosenbergii', 2), 
                                   rep('schistosoma haematobium', 2), 'Bulinus truncatus', 'Biomphalaria alexandrina',
                                   'Bulinus truncatus', 'Biomphalaria alexandrina', NA),
                       Parameter = c('muP',  'muP',  'psiP', 'piC', 'piM',
                                      'muN', 'muN', 'fN', 'fN', 'Combined'),
                       r0 = c(mean(r0.ch.eec[, 1]), mean(r0.ch.eec[, 2]), mean(r0.ch.eec[, 3]),
                              mean(r0.ch.eec[, 4]), mean(r0.ch.eec[, 5]), mean(r0.ch.eec[, 6]),
                              mean(r0.ch.eec[, 7]), mean(r0.ch.eec[, 8]), mean(r0.ch.eec[, 9]),
                              mean(r0.ch.eec[, 10])),
                       r0.sd = c(sd(r0.ch.eec[, 1]), sd(r0.ch.eec[, 2]), sd(r0.ch.eec[, 3]),
                                 sd(r0.ch.eec[, 4]), sd(r0.ch.eec[, 5]), sd(r0.ch.eec[, 6]),
                                 sd(r0.ch.eec[, 7]), sd(r0.ch.eec[, 8]), sd(r0.ch.eec[, 9]),
                                 sd(r0.ch.eec[, 10])),
                       par.mean = c(mean(par.ch.eec[, 1]), mean(par.ch.eec[, 2]), mean(par.ch.eec[, 3]),
                                    mean(par.ch.eec[, 4]), mean(par.ch.eec[, 5]), mean(par.ch.eec[, 6]),
                                    mean(par.ch.eec[, 7]), mean(par.ch.eec[, 8]), mean(par.ch.eec[, 9]),0))

ch.eec.df$r0.up = ch.eec.df$r0 + ch.eec.df$r0.sd
ch.eec.df$r0.lo = ch.eec.df$r0 - ch.eec.df$r0.sd

save(ch.eec.df, file = 'Review_models/r0_EECs/ch.eec.df.RData')

#10% EEC runs
ch.0.5eec.df = data.frame(chem = rep('Chlorpyrifos', length(parfx)+1),
                          study = c('Halstead et al, 2015', rep('Satapornvanit et al, 2009', 2),
                                    rep('Hasheesh & Mohamed, 2011', 3), 'Ibrahim et al, 1992', 
                                    'Hasheesh & Mohamed, 2011', 'Ibrahim et al, 1992', 'Combined'),
                          Species = c('Procambarus clarkii', rep('Macrobrachium rosenbergii', 2), 
                                      rep('schistosoma haematobium', 2), 'Bulinus truncatus', 'Biomphalaria alexandrina',
                                      'Bulinus truncatus', 'Biomphalaria alexandrina', NA),
                          Parameter = c('muP',  'muP',  'psiP', 'piC', 'piM',
                                        'muN', 'muN', 'fN', 'fN', 'Combined'),
                          r0 = c(mean(r0.ch.0.5eec[, 1]), mean(r0.ch.0.5eec[, 2]), mean(r0.ch.0.5eec[, 3]),
                                 mean(r0.ch.0.5eec[, 4]), mean(r0.ch.0.5eec[, 5]), mean(r0.ch.0.5eec[, 6]),
                                 mean(r0.ch.0.5eec[, 7]), mean(r0.ch.0.5eec[, 8]), mean(r0.ch.0.5eec[, 9]),
                                 mean(r0.ch.0.5eec[, 10])),
                          r0.sd = c(sd(r0.ch.0.5eec[, 1]), sd(r0.ch.0.5eec[, 2]), sd(r0.ch.0.5eec[, 3]),
                                    sd(r0.ch.0.5eec[, 4]), sd(r0.ch.0.5eec[, 5]), sd(r0.ch.0.5eec[, 6]),
                                    sd(r0.ch.0.5eec[, 7]), sd(r0.ch.0.5eec[, 8]), sd(r0.ch.0.5eec[, 9]),
                                    sd(r0.ch.0.5eec[, 10])),
                          par.mean = c(mean(par.ch.0.5eec[, 1]), mean(par.ch.0.5eec[, 2]), mean(par.ch.0.5eec[, 3]),
                                       mean(par.ch.0.5eec[, 4]), mean(par.ch.0.5eec[, 5]), mean(par.ch.0.5eec[, 6]),
                                       mean(par.ch.0.5eec[, 7]), mean(par.ch.0.5eec[, 8]), mean(par.ch.0.5eec[, 9]),0))

ch.0.5eec.df$r0.up = ch.0.5eec.df$r0 + ch.0.5eec.df$r0.sd
ch.0.5eec.df$r0.lo = ch.0.5eec.df$r0 - ch.0.5eec.df$r0.sd

save(ch.0.5eec.df, file = 'Review_models/r0_EECs/ch.0.5eec.df.RData')

#10% EEC runs
ch.0.1eec.df = data.frame(chem = rep('Chlorpyrifos', length(parfx)+1),
                       study = c('Halstead et al, 2015', rep('Satapornvanit et al, 2009', 2),
                                 rep('Hasheesh & Mohamed, 2011', 3), 'Ibrahim et al, 1992', 
                                 'Hasheesh & Mohamed, 2011', 'Ibrahim et al, 1992', 'Combined'),
                       Species = c('Procambarus clarkii', rep('Macrobrachium rosenbergii', 2), 
                                   rep('schistosoma haematobium', 2), 'Bulinus truncatus', 'Biomphalaria alexandrina',
                                   'Bulinus truncatus', 'Biomphalaria alexandrina', NA),
                       Parameter = c('muP',  'muP',  'psiP', 'piC', 'piM',
                                     'muN', 'muN', 'fN', 'fN', 'Combined'),
                       r0 = c(mean(r0.ch.0.1eec[, 1]), mean(r0.ch.0.1eec[, 2]), mean(r0.ch.0.1eec[, 3]),
                              mean(r0.ch.0.1eec[, 4]), mean(r0.ch.0.1eec[, 5]), mean(r0.ch.0.1eec[, 6]),
                              mean(r0.ch.0.1eec[, 7]), mean(r0.ch.0.1eec[, 8]), mean(r0.ch.0.1eec[, 9]),
                              mean(r0.ch.0.1eec[, 10])),
                       r0.sd = c(sd(r0.ch.0.1eec[, 1]), sd(r0.ch.0.1eec[, 2]), sd(r0.ch.0.1eec[, 3]),
                                 sd(r0.ch.0.1eec[, 4]), sd(r0.ch.0.1eec[, 5]), sd(r0.ch.0.1eec[, 6]),
                                 sd(r0.ch.0.1eec[, 7]), sd(r0.ch.0.1eec[, 8]), sd(r0.ch.0.1eec[, 9]),
                                 sd(r0.ch.0.1eec[, 10])),
                       par.mean = c(mean(par.ch.0.1eec[, 1]), mean(par.ch.0.1eec[, 2]), mean(par.ch.0.1eec[, 3]),
                                    mean(par.ch.0.1eec[, 4]), mean(par.ch.0.1eec[, 5]), mean(par.ch.0.1eec[, 6]),
                                    mean(par.ch.0.1eec[, 7]), mean(par.ch.0.1eec[, 8]), mean(par.ch.0.1eec[, 9]),0))

ch.0.1eec.df$r0.up = ch.0.1eec.df$r0 + ch.0.1eec.df$r0.sd
ch.0.1eec.df$r0.lo = ch.0.1eec.df$r0 - ch.0.1eec.df$r0.sd

save(ch.0.1eec.df, file = 'Review_models/r0_EECs/ch.0.1eec.df.RData')