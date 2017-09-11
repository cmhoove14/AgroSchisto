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
source('Response_Fxs/Satapornvanit_Insecticides2009.R')

source('Review_models/r0_of_q.R')

library(parallel)
library(fBasics)

keep.fin.p3 = c(keep.hal15.muP, keep.all.sat09,
                 'r0.In', 'r0.fix', 'parameters', 'nil0', 'nil1', 'keep.fin.p3')

rm(list = setdiff(ls(), keep.fin.p3))
dev.off()

no.cores = detectCores() - 1

#Run simulations of malathion concentrations, start with individual functions then combine ################
nsims = 1000         #Number of simulations to run
#Expected environmental concentrations for agrochemicals (NEED TO BE FINALIZED)
  eec.ch = rep(64, nsims)       #chlorpyrifos
  eec.dim = rep(100, nsims)     #dimethoate
  eec.esf = rep(1.03, nsims)    #esfenvalerate
  eec.lcy = rep(1.77, nsims)    #lamda-cyhalothrin
  eec.mal = rep(583, nsims)     #malathion
  eec.per = rep(5.98, nsims)    #permethrin 
  eec.pro = rep(100, nsims)     #profenofos
  eec.trb = rep(36.6, nsims)    #terbufos
  eec.znc = rep(100, nsims)     #zinc sulfate
  
#All pathway 3 response functions  
parfx = c(muPq_chlor_satapornvanit09_uncertainty, muPq_chlor_Halstead_uncertainty, muPq_dim_satapornvanit09_uncertainty,
          muPq_esfen_Halstead_uncertainty, muPq_lamcy_Halstead_uncertainty, muPq_mal_Halstead_uncertainty,
          muPq_perm_Halstead_uncertainty, muPq_prof_satapornvanit09_uncertainty, muPq_terb_Halstead_uncertainty,
          muPq_zinc_satapornvanit09_uncertainty, psi_q_chlor_satapornvanit09_uncertainty, psi_q_zinc_satapornvanit09_uncertainty)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.p3, 'uniroot.all', 'rdrm', 'LL.2', 'LL.3','L.4'))

r0.eec.p3 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.eec.p3 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.0.5eec.p3 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.0.5eec.p3 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.0.1eec.p3 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.0.1eec.p3 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

set.seed(0)

#Fill r0 estimates for EEC #####################
  #individual parameters
    r0.eec.p3[, 1] = parSapply(clus1, eec.ch, r0.In, f.mu_Pq = muPq_chlor_satapornvanit09_uncertainty)[3,]
    r0.eec.p3[, 2] = parSapply(clus1, eec.ch, r0.In, f.mu_Pq = muPq_chlor_Halstead_uncertainty)[3,]
    r0.eec.p3[, 3] = parSapply(clus1, eec.dim, r0.In, f.mu_Pq = muPq_dim_satapornvanit09_uncertainty)[3,]
    r0.eec.p3[, 4] = parSapply(clus1, eec.esf, r0.In, f.mu_Pq = muPq_esfen_Halstead_uncertainty)[3,]
    r0.eec.p3[, 5] = parSapply(clus1, eec.lcy, r0.In, f.mu_Pq = muPq_lamcy_Halstead_uncertainty)[3,]
    r0.eec.p3[, 6] = parSapply(clus1, eec.mal, r0.In, f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    r0.eec.p3[, 7] = parSapply(clus1, eec.per, r0.In, f.mu_Pq = muPq_perm_Halstead_uncertainty)[3,]
    r0.eec.p3[, 8] = parSapply(clus1, eec.pro, r0.In, f.mu_Pq = muPq_prof_satapornvanit09_uncertainty)[3,]
    r0.eec.p3[, 9] = parSapply(clus1, eec.trb, r0.In, f.mu_Pq = muPq_terb_Halstead_uncertainty)[3,]
    r0.eec.p3[, 10] = parSapply(clus1, eec.znc, r0.In, f.mu_Pq = muPq_zinc_satapornvanit09_uncertainty)[3,]
    r0.eec.p3[, 11] = parSapply(clus1, eec.ch, r0.In, f.alpha_q = psi_q_chlor_satapornvanit09_uncertainty)[3,]
    r0.eec.p3[, 12] = parSapply(clus1, eec.znc, r0.In, f.alpha_q = psi_q_zinc_satapornvanit09_uncertainty)[3,]
    
    
    
  #fill parameter values for EEC values
    par.eec.p3[, 1] = parSapply(clus1, eec.ch, muPq_chlor_satapornvanit09_uncertainty) 
    par.eec.p3[, 2] = parSapply(clus1, eec.ch, muPq_chlor_Halstead_uncertainty) 
    par.eec.p3[, 3] = parSapply(clus1, eec.dim, muPq_dim_satapornvanit09_uncertainty) 
    par.eec.p3[, 4] = parSapply(clus1, eec.esf, muPq_esfen_Halstead_uncertainty) 
    par.eec.p3[, 5] = parSapply(clus1, eec.lcy, muPq_lamcy_Halstead_uncertainty) 
    par.eec.p3[, 6] = parSapply(clus1, eec.mal, muPq_mal_Halstead_uncertainty) 
    par.eec.p3[, 7] = parSapply(clus1, eec.per, muPq_perm_Halstead_uncertainty) 
    par.eec.p3[, 8] = parSapply(clus1, eec.pro, muPq_prof_satapornvanit09_uncertainty) 
    par.eec.p3[, 9] = parSapply(clus1, eec.trb, muPq_terb_Halstead_uncertainty) 
    par.eec.p3[, 10] = parSapply(clus1, eec.znc, muPq_zinc_satapornvanit09_uncertainty) 
    par.eec.p3[, 11] = parSapply(clus1, eec.ch, psi_q_chlor_satapornvanit09_uncertainty) 
    par.eec.p3[, 12] = parSapply(clus1, eec.znc, psi_q_zinc_satapornvanit09_uncertainty) 
    
    
#Fill r0 estimates for 50% EEC ########################
  #individual parameters
     r0.0.5eec.p3[, 1] = parSapply(clus1, 0.5*eec.ch, r0.In, f.mu_Pq = muPq_chlor_satapornvanit09_uncertainty)[3,]
     r0.0.5eec.p3[, 2] = parSapply(clus1, 0.5*eec.ch, r0.In, f.mu_Pq = muPq_chlor_Halstead_uncertainty)[3,]
     r0.0.5eec.p3[, 3] = parSapply(clus1, 0.5*eec.dim, r0.In, f.mu_Pq = muPq_dim_satapornvanit09_uncertainty)[3,]
     r0.0.5eec.p3[, 4] = parSapply(clus1, 0.5*eec.esf, r0.In, f.mu_Pq = muPq_esfen_Halstead_uncertainty)[3,]
     r0.0.5eec.p3[, 5] = parSapply(clus1, 0.5*eec.lcy, r0.In, f.mu_Pq = muPq_lamcy_Halstead_uncertainty)[3,]
     r0.0.5eec.p3[, 6] = parSapply(clus1, 0.5*eec.mal, r0.In, f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
     r0.0.5eec.p3[, 7] = parSapply(clus1, 0.5*eec.per, r0.In, f.mu_Pq = muPq_perm_Halstead_uncertainty)[3,]
     r0.0.5eec.p3[, 8] = parSapply(clus1, 0.5*eec.pro, r0.In, f.mu_Pq = muPq_prof_satapornvanit09_uncertainty)[3,]
     r0.0.5eec.p3[, 9] = parSapply(clus1, 0.5*eec.trb, r0.In, f.mu_Pq = muPq_terb_Halstead_uncertainty)[3,]
     r0.0.5eec.p3[, 10] = parSapply(clus1, 0.5*eec.znc, r0.In, f.mu_Pq = muPq_zinc_satapornvanit09_uncertainty)[3,]
     r0.0.5eec.p3[, 11] = parSapply(clus1, 0.5*eec.ch, r0.In, f.alpha_q = psi_q_chlor_satapornvanit09_uncertainty)[3,]
     r0.0.5eec.p3[, 12] = parSapply(clus1, 0.5*eec.znc, r0.In, f.alpha_q = psi_q_zinc_satapornvanit09_uncertainty)[3,]
    
    
    
  #fill parameter values for EEC values
    par.0.5eec.p3[, 1] = parSapply(clus1, 0.5*eec.ch, muPq_chlor_satapornvanit09_uncertainty) 
    par.0.5eec.p3[, 2] = parSapply(clus1, 0.5*eec.ch, muPq_chlor_Halstead_uncertainty) 
    par.0.5eec.p3[, 3] = parSapply(clus1, 0.5*eec.dim, muPq_dim_satapornvanit09_uncertainty) 
    par.0.5eec.p3[, 4] = parSapply(clus1, 0.5*eec.esf, muPq_esfen_Halstead_uncertainty) 
    par.0.5eec.p3[, 5] = parSapply(clus1, 0.5*eec.lcy, muPq_lamcy_Halstead_uncertainty) 
    par.0.5eec.p3[, 6] = parSapply(clus1, 0.5*eec.mal, muPq_mal_Halstead_uncertainty) 
    par.0.5eec.p3[, 7] = parSapply(clus1, 0.5*eec.per, muPq_perm_Halstead_uncertainty) 
    par.0.5eec.p3[, 8] = parSapply(clus1, 0.5*eec.pro, muPq_prof_satapornvanit09_uncertainty) 
    par.0.5eec.p3[, 9] = parSapply(clus1, 0.5*eec.trb, muPq_terb_Halstead_uncertainty) 
    par.0.5eec.p3[, 10] = parSapply(clus1, 0.5*eec.znc, muPq_zinc_satapornvanit09_uncertainty) 
    par.0.5eec.p3[, 11] = parSapply(clus1, 0.5*eec.ch, psi_q_chlor_satapornvanit09_uncertainty) 
    par.0.5eec.p3[, 12] = parSapply(clus1, 0.5*eec.znc, psi_q_zinc_satapornvanit09_uncertainty) 
    
    
#Fill r0 estimates for 10% EEC #####################
  #individual parameters
    r0.0.1eec.p3[, 1] = parSapply(clus1, 0.1*eec.ch, r0.In, f.mu_Pq = muPq_chlor_satapornvanit09_uncertainty)[3,]
    r0.0.1eec.p3[, 2] = parSapply(clus1, 0.1*eec.ch, r0.In, f.mu_Pq = muPq_chlor_Halstead_uncertainty)[3,]
    r0.0.1eec.p3[, 3] = parSapply(clus1, 0.1*eec.dim, r0.In, f.mu_Pq = muPq_dim_satapornvanit09_uncertainty)[3,]
    r0.0.1eec.p3[, 4] = parSapply(clus1, 0.1*eec.esf, r0.In, f.mu_Pq = muPq_esfen_Halstead_uncertainty)[3,]
    r0.0.1eec.p3[, 5] = parSapply(clus1, 0.1*eec.lcy, r0.In, f.mu_Pq = muPq_lamcy_Halstead_uncertainty)[3,]
    r0.0.1eec.p3[, 6] = parSapply(clus1, 0.1*eec.mal, r0.In, f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    r0.0.1eec.p3[, 7] = parSapply(clus1, 0.1*eec.per, r0.In, f.mu_Pq = muPq_perm_Halstead_uncertainty)[3,]
    r0.0.1eec.p3[, 8] = parSapply(clus1, 0.1*eec.pro, r0.In, f.mu_Pq = muPq_prof_satapornvanit09_uncertainty)[3,]
    r0.0.1eec.p3[, 9] = parSapply(clus1, 0.1*eec.trb, r0.In, f.mu_Pq = muPq_terb_Halstead_uncertainty)[3,]
    r0.0.1eec.p3[, 10] = parSapply(clus1, 0.1*eec.znc, r0.In, f.mu_Pq = muPq_zinc_satapornvanit09_uncertainty)[3,]
    r0.0.1eec.p3[, 11] = parSapply(clus1, 0.1*eec.ch, r0.In, f.alpha_q = psi_q_chlor_satapornvanit09_uncertainty)[3,]
    r0.0.1eec.p3[, 12] = parSapply(clus1, 0.1*eec.znc, r0.In, f.alpha_q = psi_q_zinc_satapornvanit09_uncertainty)[3,]
    
    
    
  #fill parameter values for EEC values
    par.0.1eec.p3[, 1] = parSapply(clus1, 0.1*eec.ch, muPq_chlor_satapornvanit09_uncertainty) 
    par.0.1eec.p3[, 2] = parSapply(clus1, 0.1*eec.ch, muPq_chlor_Halstead_uncertainty) 
    par.0.1eec.p3[, 3] = parSapply(clus1, 0.1*eec.dim, muPq_dim_satapornvanit09_uncertainty) 
    par.0.1eec.p3[, 4] = parSapply(clus1, 0.1*eec.esf, muPq_esfen_Halstead_uncertainty) 
    par.0.1eec.p3[, 5] = parSapply(clus1, 0.1*eec.lcy, muPq_lamcy_Halstead_uncertainty) 
    par.0.1eec.p3[, 6] = parSapply(clus1, 0.1*eec.mal, muPq_mal_Halstead_uncertainty) 
    par.0.1eec.p3[, 7] = parSapply(clus1, 0.1*eec.per, muPq_perm_Halstead_uncertainty) 
    par.0.1eec.p3[, 8] = parSapply(clus1, 0.1*eec.pro, muPq_prof_satapornvanit09_uncertainty) 
    par.0.1eec.p3[, 9] = parSapply(clus1, 0.1*eec.trb, muPq_terb_Halstead_uncertainty) 
    par.0.1eec.p3[, 10] = parSapply(clus1, 0.1*eec.znc, muPq_zinc_satapornvanit09_uncertainty) 
    par.0.1eec.p3[, 11] = parSapply(clus1, 0.1*eec.ch, psi_q_chlor_satapornvanit09_uncertainty) 
    par.0.1eec.p3[, 12] = parSapply(clus1, 0.1*eec.znc, psi_q_zinc_satapornvanit09_uncertainty) 

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
                       par.mean = colMeans(par.eec.p3))

eec.p3.df$r0.up = eec.p3.df$r0 + eec.p3.df$r0.sd
eec.p3.df$r0.lo = eec.p3.df$r0 - eec.p3.df$r0.sd

save(eec.p3.df, file = 'Review_models/r0_EECs/eec.p3.df.RData')

#50% EEC values #################
eec0.5.p3.df = data.frame(chem = c('Chlorpyrifos', 'Chlorpyrifos', 'Dimethoate', 'Esfenvalerate', 'Lamda-Cyhalothrin',
                                   'Malathion', 'Permethrin', 'Profenofos', 'Terbufos', 'Zinc', 'Chlorpyrifos', 'Zinc'),
                          study = c('Satapornvanit et al 2009', 'Halstead et al 2015', 'Satapornvanit et al 2009', 
                                    rep('Halstead et al 2015', 4), 'Satapornvanit et al 2009', 'Halstead et al 2015',
                                    rep('Satapornvanit et al 2009', 3)),
                          Species = c('Macrobrachium rosenbergii', 'Procambarus clarkii', 'Macrobrachium rosenbergii', 
                                      rep('Procambarus clarkii', 4), 'Macrobrachium rosenbergii', 'Procambarus clarkii',
                                      rep('Macrobrachium rosenbergii', 3)),
                          Parameter = c(rep('muP', 10), rep('psi', 2)),
                          r0 = colMeans(r0.0.5eec.p3),
                          r0.sd = apply(r0.0.5eec.p3, 2, sd),
                          par.mean = colMeans(par.0.5eec.p3))

eec0.5.p3.df$r0.up = eec0.5.p3.df$r0 + eec0.5.p3.df$r0.sd
eec0.5.p3.df$r0.lo = eec0.5.p3.df$r0 - eec0.5.p3.df$r0.sd

save(eec0.5.p3.df, file = 'Review_models/r0_EECs/eec0.5.p3.df.RData')


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
                          par.mean = colMeans(par.0.1eec.p3))

eec0.1.p3.df$r0.up = eec0.1.p3.df$r0 + eec0.1.p3.df$r0.sd
eec0.1.p3.df$r0.lo = eec0.1.p3.df$r0 - eec0.1.p3.df$r0.sd

save(eec0.1.p3.df, file = 'Review_models/r0_EECs/eec0.1.p3.df.RData')
