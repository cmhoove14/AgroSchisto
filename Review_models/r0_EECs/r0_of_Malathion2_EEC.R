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
source('Response_Fxs/tchounwou92_piC2_beq.R')
source('Response_Fxs/tchounwou91_piM_beq.R')
source('Response_Fxs/Halstead_Insecticides2015.R')
source('Response_Fxs/bakry2011.R')
source('Response_Fxs/tchounwou91_fN-muN.R')

source('Review_models/r0_of_q.R')

library(parallel)
library(fBasics)

keep.fin.mal = c(keep.tch92.beq, keep.tch91.beq, keep.tch91.snail, 
                 keep.hal15.muP[c(1,7,18,24)], keep.bak.mal, 
                 'r0.In', 'r0.fix', 'parameters', 'nil0', 'nil1', 'keep.fin.mal')

rm(list = setdiff(ls(), keep.fin.mal))
dev.off()

no.cores = detectCores() - 1

#Run simulations of malathion concentrations, start with individual functions then combine ################
eec.mal = 583
nsims = 1000         #Number of simulations to run
  conc.mal = rep(eec.mal, nsims)
  
parfx = c(piC.tch92_mal_unc, piM.tch91_mal_unc, #Functions corresponding to affected parameters
          fN.mal.fx.uncertainty, muNq_mal_Bakry11_uncertainty, muPq_mal_Halstead_uncertainty,
          muNq_mal_tch91_uncertainty, fNq_mal_tch91_uncertainty)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.mal, 'uniroot.all', 'rdrm', 'LL.2', 'L.4'))

r0.mal.eec = matrix(data = NA, nrow = nsims, ncol = length(parfx)+1)
par.mal.eec = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.mal.0.5eec = matrix(data = NA, nrow = nsims, ncol = length(parfx)+1)
par.mal.0.5eec = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.mal.0.1eec = matrix(data = NA, nrow = nsims, ncol = length(parfx)+1)
par.mal.0.1eec = matrix(data = NA, nrow = nsims, ncol = length(parfx))

set.seed(0)

  #Fill r0 estimates for EEC
  #individual parameters
    r0.mal.eec[, 1] = parSapply(clus1, conc.mal, r0.In, f.pi_Cq = piC.tch92_mal_unc)[3,]
    r0.mal.eec[, 2] = parSapply(clus1, conc.mal, r0.In, f.pi_Mq = piM.tch91_mal_unc)[3,]
    r0.mal.eec[, 3] = parSapply(clus1, conc.mal, r0.In, f.f_Nq = fN.mal.fx.uncertainty)[3,]
    r0.mal.eec[, 4] = parSapply(clus1, conc.mal, r0.In, f.mu_Nq = muNq_mal_Bakry11_uncertainty)[3,]
    r0.mal.eec[, 5] = parSapply(clus1, conc.mal, r0.In, f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    r0.mal.eec[, 6] = parSapply(clus1, conc.mal, r0.In, f.mu_Nq = muNq_mal_tch91_uncertainty)[3,]
    r0.mal.eec[, 7] = parSapply(clus1, conc.mal, r0.In, f.f_Nq = fNq_mal_tch91_uncertainty)[3,]
    
  #Full simulation w/all parameters
    r0.mal.eec[, 8] = parSapply(clus1, conc.mal, r0.In,
                                    f.mu_Nq = muNq_mal_tch91_uncertainty,
                                    f.f_Nq = fNq_mal_tch91_uncertainty,
                                    f.pi_Mq = piM.tch91_mal_unc,
                                    f.pi_Cq = piC.tch92_mal_unc,
                                    f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]

  #fill parameter values for EEC values
    par.mal.eec[, 1] = parSapply(clus1, conc.mal, piC.tch92_mal_unc)
    par.mal.eec[, 2] = parSapply(clus1, conc.mal, piM.tch91_mal_unc)
    par.mal.eec[, 3] = parSapply(clus1, conc.mal, fN.mal.fx.uncertainty)
    par.mal.eec[, 4] = parSapply(clus1, conc.mal, muNq_mal_Bakry11_uncertainty)
    par.mal.eec[, 5] = parSapply(clus1, conc.mal, muPq_mal_Halstead_uncertainty)
    par.mal.eec[, 6] = parSapply(clus1, conc.mal, muNq_mal_tch91_uncertainty)
    par.mal.eec[, 7] = parSapply(clus1, conc.mal, fNq_mal_tch91_uncertainty)
    
  #Fill r0 estimates for 50% EEC
    #individual parameters
    r0.mal.0.5eec[, 1] = parSapply(clus1, conc.mal*0.5, r0.In, f.pi_Cq = piC.tch92_mal_unc)[3,]
    r0.mal.0.5eec[, 2] = parSapply(clus1, conc.mal*0.5, r0.In, f.pi_Mq = piM.tch91_mal_unc)[3,]
    r0.mal.0.5eec[, 3] = parSapply(clus1, conc.mal*0.5, r0.In, f.f_Nq = fN.mal.fx.uncertainty)[3,]
    r0.mal.0.5eec[, 4] = parSapply(clus1, conc.mal*0.5, r0.In, f.mu_Nq = muNq_mal_Bakry11_uncertainty)[3,]
    r0.mal.0.5eec[, 5] = parSapply(clus1, conc.mal*0.5, r0.In, f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    r0.mal.0.5eec[, 6] = parSapply(clus1, conc.mal*0.5, r0.In, f.mu_Nq = muNq_mal_tch91_uncertainty)[3,]
    r0.mal.0.5eec[, 7] = parSapply(clus1, conc.mal*0.5, r0.In, f.f_Nq = fNq_mal_tch91_uncertainty)[3,]
    
  #Full simulation w/all parameters
    r0.mal.0.5eec[, 8] = parSapply(clus1, conc.mal*0.5, r0.In,
                                   f.mu_Nq = muNq_mal_tch91_uncertainty,
                                   f.f_Nq = fNq_mal_tch91_uncertainty,
                                   f.pi_Mq = piM.tch91_mal_unc,
                                   f.pi_Cq = piC.tch92_mal_unc,
                                   f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    
  #fill parameter values for 50% EEC values
    par.mal.0.5eec[, 1] = parSapply(clus1, conc.mal*0.5, piC.tch92_mal_unc)
    par.mal.0.5eec[, 2] = parSapply(clus1, conc.mal*0.5, piM.tch91_mal_unc)
    par.mal.0.5eec[, 3] = parSapply(clus1, conc.mal*0.5, fN.mal.fx.uncertainty)
    par.mal.0.5eec[, 4] = parSapply(clus1, conc.mal*0.5, muNq_mal_Bakry11_uncertainty)
    par.mal.0.5eec[, 5] = parSapply(clus1, conc.mal*0.5, muPq_mal_Halstead_uncertainty)
    par.mal.0.5eec[, 6] = parSapply(clus1, conc.mal*0.5, muNq_mal_tch91_uncertainty)
    par.mal.0.5eec[, 7] = parSapply(clus1, conc.mal*0.5, fNq_mal_tch91_uncertainty)
    
    
  #Fill r0 estimates for 10% EEC
  #individual parameters
    r0.mal.0.1eec[, 1] = parSapply(clus1, conc.mal*0.1, r0.In, f.pi_Cq = piC.tch92_mal_unc)[3,]
    r0.mal.0.1eec[, 2] = parSapply(clus1, conc.mal*0.1, r0.In, f.pi_Mq = piM.tch91_mal_unc)[3,]
    r0.mal.0.1eec[, 3] = parSapply(clus1, conc.mal*0.1, r0.In, f.f_Nq = fN.mal.fx.uncertainty)[3,]
    r0.mal.0.1eec[, 4] = parSapply(clus1, conc.mal*0.1, r0.In, f.mu_Nq = muNq_mal_Bakry11_uncertainty)[3,]
    r0.mal.0.1eec[, 5] = parSapply(clus1, conc.mal*0.1, r0.In, f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    r0.mal.0.1eec[, 6] = parSapply(clus1, conc.mal*0.1, r0.In, f.mu_Nq = muNq_mal_tch91_uncertainty)[3,]
    r0.mal.0.1eec[, 7] = parSapply(clus1, conc.mal*0.1, r0.In, f.f_Nq = fNq_mal_tch91_uncertainty)[3,]
  
  #Full simulation w/all parameters
    r0.mal.0.1eec[, 8] = parSapply(clus1, conc.mal*0.1, r0.In,
                                f.mu_Nq = muNq_mal_tch91_uncertainty,
                                f.f_Nq = fNq_mal_tch91_uncertainty,
                                f.pi_Mq = piM.tch91_mal_unc,
                                f.pi_Cq = piC.tch92_mal_unc,
                                f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
  
  #fill parameter values for 10% EEC values
  par.mal.0.1eec[, 1] = parSapply(clus1, conc.mal*0.1, piC.tch92_mal_unc)
  par.mal.0.1eec[, 2] = parSapply(clus1, conc.mal*0.1, piM.tch91_mal_unc)
  par.mal.0.1eec[, 3] = parSapply(clus1, conc.mal*0.1, fN.mal.fx.uncertainty)
  par.mal.0.1eec[, 4] = parSapply(clus1, conc.mal*0.1, muNq_mal_Bakry11_uncertainty)
  par.mal.0.1eec[, 5] = parSapply(clus1, conc.mal*0.1, muPq_mal_Halstead_uncertainty)
  par.mal.0.1eec[, 6] = parSapply(clus1, conc.mal*0.1, muNq_mal_tch91_uncertainty)
  par.mal.0.1eec[, 7] = parSapply(clus1, conc.mal*0.1, fNq_mal_tch91_uncertainty)

stopCluster(clus1)
########## Post process ############ 
#EEC runs
mal.eec.df = data.frame(chem = rep('Malathion', length(parfx)+1),
                       study = c('Tchounwou et al, 1992', 'Tchounwou et al, 1991a', 'Bakry et al, 2011',
                                 'Bakry et al, 2011', 'Halstead et al, 2015', 
                                 'Tchounwou et al, 1991b', 'Tchounwou et al, 1991b', 'Combined'),
                       Species = c('Schistosoma mansoni', 'Schistosoma mansoni', 'Helisoma duryi','Helisoma duryi',
                                   'Procambarus clarkii', 'Bulinus havenensis', 'Bulinus havenensis', NA),
                       Parameter = c('piC', 'piM', 'fN','muN', 'muP',
                                     'muN', 'fN', 'Combined'),
                       r0 = c(mean(r0.mal.eec[, 1]), mean(r0.mal.eec[, 2]), mean(r0.mal.eec[, 3]),
                              mean(r0.mal.eec[, 4]), mean(r0.mal.eec[, 5]), mean(r0.mal.eec[, 6]),
                              mean(r0.mal.eec[, 7]), mean(r0.mal.eec[, 8])),
                       r0.sd = c(sd(r0.mal.eec[, 1]), sd(r0.mal.eec[, 2]), sd(r0.mal.eec[, 3]),
                                 sd(r0.mal.eec[, 4]), sd(r0.mal.eec[, 5]), sd(r0.mal.eec[, 6]),
                                 sd(r0.mal.eec[, 7]), sd(r0.mal.eec[, 8])),
                       par.mean = c(mean(par.mal.eec[, 1]), mean(par.mal.eec[, 2]), mean(par.mal.eec[, 3]),
                                    mean(par.mal.eec[, 4]), mean(par.mal.eec[, 5]), mean(par.mal.eec[, 6]),
                                    mean(par.mal.eec[, 7]), 0))

mal.eec.df$r0.up = mal.eec.df$r0 + mal.eec.df$r0.sd
mal.eec.df$r0.lo = mal.eec.df$r0 - mal.eec.df$r0.sd

save(mal.eec.df, file = 'Review_models/r0_EECs/mal.eec.df.RData')

#50% EEC values
mal.0.5eec.df = data.frame(chem = rep('Malathion', length(parfx)+1),
                           study = c('Tchounwou et al, 1992', 'Tchounwou et al, 1991a', 'Bakry et al, 2011',
                                     'Bakry et al, 2011', 'Halstead et al, 2015', 
                                     'Tchounwou et al, 1991b', 'Tchounwou et al, 1991b', 'Combined'),
                           Species = c('Schistosoma mansoni', 'Schistosoma mansoni', 'Helisoma duryi','Helisoma duryi',
                                       'Procambarus clarkii', 'Bulinus havenensis', 'Bulinus havenensis', NA),
                           Parameter = c('piC', 'piM', 'fN','muN', 'muP',
                                         'muN', 'fN', 'Combined'),
                           r0 = c(mean(r0.mal.0.5eec[, 1]), mean(r0.mal.0.5eec[, 2]), mean(r0.mal.0.5eec[, 3]),
                                  mean(r0.mal.0.5eec[, 4]), mean(r0.mal.0.5eec[, 5]), mean(r0.mal.0.5eec[, 6]),
                                  mean(r0.mal.0.5eec[, 7]), mean(r0.mal.0.5eec[, 8])),
                           r0.sd = c(sd(r0.mal.0.5eec[, 1]), sd(r0.mal.0.5eec[, 2]), sd(r0.mal.0.5eec[, 3]),
                                     sd(r0.mal.0.5eec[, 4]), sd(r0.mal.0.5eec[, 5]), sd(r0.mal.0.5eec[, 6]),
                                     sd(r0.mal.0.5eec[, 7]), sd(r0.mal.0.5eec[, 8])),
                           par.mean = c(mean(par.mal.0.5eec[, 1]), mean(par.mal.0.5eec[, 2]), mean(par.mal.0.5eec[, 3]),
                                        mean(par.mal.0.5eec[, 4]), mean(par.mal.0.5eec[, 5]), mean(par.mal.0.5eec[, 6]),
                                        mean(par.mal.0.5eec[, 7]), 0))

mal.0.5eec.df$r0.up = mal.0.5eec.df$r0 + mal.0.5eec.df$r0.sd
mal.0.5eec.df$r0.lo = mal.0.5eec.df$r0 - mal.0.5eec.df$r0.sd

save(mal.0.5eec.df, file = 'Review_models/r0_EECs/mal.0.5eec.df.RData')


#10% EEC values
mal.0.1eec.df = data.frame(chem = rep('Malathion', length(parfx)+1),
                        study = c('Tchounwou et al, 1992', 'Tchounwou et al, 1991a', 'Bakry et al, 2011',
                                  'Bakry et al, 2011', 'Halstead et al, 2015', 
                                  'Tchounwou et al, 1991b', 'Tchounwou et al, 1991b', 'Combined'),
                        Species = c('Schistosoma mansoni', 'Schistosoma mansoni', 'Helisoma duryi','Helisoma duryi',
                                    'Procambarus clarkii', 'Bulinus havenensis', 'Bulinus havenensis', NA),
                        Parameter = c('piC', 'piM', 'fN','muN', 'muP',
                                      'muN', 'fN', 'Combined'),
                        r0 = c(mean(r0.mal.0.1eec[, 1]), mean(r0.mal.0.1eec[, 2]), mean(r0.mal.0.1eec[, 3]),
                               mean(r0.mal.0.1eec[, 4]), mean(r0.mal.0.1eec[, 5]), mean(r0.mal.0.1eec[, 6]),
                               mean(r0.mal.0.1eec[, 7]), mean(r0.mal.0.1eec[, 8])),
                        r0.sd = c(sd(r0.mal.0.1eec[, 1]), sd(r0.mal.0.1eec[, 2]), sd(r0.mal.0.1eec[, 3]),
                                  sd(r0.mal.0.1eec[, 4]), sd(r0.mal.0.1eec[, 5]), sd(r0.mal.0.1eec[, 6]),
                                  sd(r0.mal.0.1eec[, 7]), sd(r0.mal.0.1eec[, 8])),
                        par.mean = c(mean(par.mal.0.1eec[, 1]), mean(par.mal.0.1eec[, 2]), mean(par.mal.0.1eec[, 3]),
                                     mean(par.mal.0.1eec[, 4]), mean(par.mal.0.1eec[, 5]), mean(par.mal.0.1eec[, 6]),
                                     mean(par.mal.0.1eec[, 7]), 0))

mal.0.1eec.df$r0.up = mal.0.1eec.df$r0 + mal.0.1eec.df$r0.sd
mal.0.1eec.df$r0.lo = mal.0.1eec.df$r0 - mal.0.1eec.df$r0.sd

save(mal.0.1eec.df, file = 'Review_models/r0_EECs/mal.0.1eec.df.RData')
