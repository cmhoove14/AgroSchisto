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
source('Response_Fxs/tchounwou92_piM_beq.R')
source('Response_Fxs/Halstead_Insecticides2015.R')
source('Response_Fxs/bakry2011.R')

source('Review_models/r0_of_q.R')

library(parallel)

keep.fin.mal = c(keep.tch92.beq, keep.tch91.beq, keep.hal15.muP[c(1,7,13)], keep.bak11.N[c(1,2,3,6,7)], 
                 'r0.In', 'r0.fix', 'parameters', 'nil0', 'nil1', 'keep.fin.mal')

rm(list = setdiff(ls(), keep.fin.mal))
dev.off()

no.cores = detectCores() - 1

#Fill data frame with point estimates of R0 driven by data in each study ################
r0.mal.fix = data.frame(mal = c(cerc.auc.mal$mal,
                                mir.auc.mal$mal,
                                mal.sum$mal.c,
                                mun.mal$conc*1000),
                        par = c(rep('piC', length(cerc.auc.mal$mal)), 
                                rep('piM', length(mir.auc.mal$mal)),
                                rep('muP', length(mal.sum$mal.c)),
                                rep('muN', length(mun.mal$conc))),
                        study = c(rep('Tch92', length(cerc.auc.mal$mal)), 
                                  rep('Tch91', length(mir.auc.mal$mal)),
                                  rep('Halstead15', length(mal.sum$mal.c)),
                                  rep('Bakry11', length(mun.mal$conc))),
                        r0 = c(r0.fix(pi_Cqx = cerc.auc.mal$piC[1])[3],
                               r0.fix(pi_Cqx = cerc.auc.mal$piC[2])[3],
                               r0.fix(pi_Cqx = cerc.auc.mal$piC[3])[3],
                               r0.fix(pi_Cqx = cerc.auc.mal$piC[4])[3],
                               r0.fix(pi_Cqx = cerc.auc.mal$piC[5])[3],
                               r0.fix(pi_Cqx = cerc.auc.mal$piC[6])[3],
                               
                               r0.fix(pi_Mqx = mir.auc.mal$piM[1])[3],
                               r0.fix(pi_Mqx = mir.auc.mal$piM[2])[3],
                               r0.fix(pi_Mqx = mir.auc.mal$piM[3])[3],
                               r0.fix(pi_Mqx = mir.auc.mal$piM[4])[3],
                               r0.fix(pi_Mqx = mir.auc.mal$piM[5])[3],
                               r0.fix(pi_Mqx = mir.auc.mal$piM[6])[3],
                               
                               r0.fix(mu_Pqx = mal.sum$mort[1])[3],
                               r0.fix(mu_Pqx = mal.sum$mort[2])[3],
                               r0.fix(mu_Pqx = mal.sum$mort[3])[3],
                               r0.fix(mu_Pqx = mal.sum$mort[4])[3],
                               r0.fix(mu_Pqx = mal.sum$mort[5])[3],
                               r0.fix(mu_Pqx = mal.sum$mort[6])[3],
                               
                               r0.fix(mu_Nqx = mun.mal$mort[1])[3],
                               r0.fix(mu_Nqx = mun.mal$mort[2])[3],
                               r0.fix(mu_Nqx = mun.mal$mort[3])[3],
                               r0.fix(mu_Nqx = mun.mal$mort[4])[3],
                               r0.fix(mu_Nqx = mun.mal$mort[5])[3]))

r0.mal.fix.n0 = subset(r0.mal.fix, mal != 0)  


#Run simulations of malathion concentrations, start with individual functions then combine ################
conc.mal = c(seq(0,2000,2), seq(2000, 200000, 200))  #Concentration range to test
nsims = 30         #Number of simulations to run
parfx = c(piC.tch92_mal_unc, piM.tch91_mal_unc, #Functions corresponding to affected parameters
          fNq_mal_Bakry11_uncertainty, muNq_mal_Bakry11_uncertainty, muPq_mal_Halstead_uncertainty)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.mal, 'uniroot.all', 'rdrm', 'LL.2'))

r0.fill.mal = array(data = NA, dim = c(length(conc.mal), nsims, length(parfx)+6))
par.fill.mal = array(data = NA, dim = c(length(conc.mal), nsims, length(parfx)))
Neq.fill.mal = array(data = NA, dim = c(length(conc.mal), nsims, length(parfx)+6))

set.seed(0)

  for(i in 1:nsims){
  #Fill r0 estimates  
  #individual parameters
    r0.fill.mal[, i, 1] = parSapply(clus1, conc.mal, r0.In, f.pi_Cq = piC.tch92_mal_unc)[3,]
    r0.fill.mal[, i, 2] = parSapply(clus1, conc.mal, r0.In, f.pi_Mq = piM.tch91_mal_unc)[3,]
    r0.fill.mal[, i, 3] = parSapply(clus1, conc.mal, r0.In, f.f_Nq = fNq_mal_Bakry11_uncertainty)[3,]
    r0.fill.mal[, i, 4] = parSapply(clus1, conc.mal, r0.In, f.mu_Nq = muNq_mal_Bakry11_uncertainty)[3,]
    r0.fill.mal[, i, 5] = parSapply(clus1, conc.mal, r0.In, f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    
  #parameter pairs and triplets
    #larval toxicity to both miracidia and cercariae
    r0.fill.mal[, i, 6] = parSapply(clus1, conc.mal, r0.In,
                                    f.pi_Mq = piM.tch91_mal_unc,
                                    f.pi_Cq = piC.tch92_mal_unc)[3,]
    #both snail effects (mortality and reproduction)
    r0.fill.mal[, i, 7] = parSapply(clus1, conc.mal, r0.In,
                                    f.mu_Nq = muNq_mal_Bakry11_uncertainty,
                                    f.f_Nq = fNq_mal_Bakry11_uncertainty)[3,]
    #larval mortality with predator mortality
    r0.fill.mal[, i, 8] = parSapply(clus1, conc.mal, r0.In,
                                    f.pi_Mq = piM.tch91_mal_unc,
                                    f.pi_Cq = piC.tch92_mal_unc,
                                    f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    #snail effects with predator mortality
    r0.fill.mal[, i, 9] = parSapply(clus1, conc.mal, r0.In,
                                    f.mu_Nq = muNq_mal_Bakry11_uncertainty,
                                    f.f_Nq = fNq_mal_Bakry11_uncertainty,
                                    f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    #Full simulation w/all parameters
    r0.fill.mal[, i, 10] = parSapply(clus1, conc.mal, r0.In,
                                              f.mu_Nq = muNq_mal_Bakry11_uncertainty,
                                              f.f_Nq = fNq_mal_Bakry11_uncertainty,
                                              f.pi_Mq = piM.tch91_mal_unc,
                                              f.pi_Cq = piC.tch92_mal_unc,
                                              f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    #Full simulation w/all parameters EXCEPT snail reproduction since function is questionable
    r0.fill.mal[, i, 11] = parSapply(clus1, conc.mal, r0.In,
                                     f.mu_Nq = muNq_mal_Bakry11_uncertainty,
                                     f.pi_Mq = piM.tch91_mal_unc,
                                     f.pi_Cq = piC.tch92_mal_unc,
                                     f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
  #fill parameter values
    par.fill.mal[, i, 1] = parSapply(clus1, conc.mal, piC.tch92_mal_unc)
    par.fill.mal[, i, 2] = parSapply(clus1, conc.mal, piM.tch91_mal_unc)
    par.fill.mal[, i, 3] = parSapply(clus1, conc.mal, fNq_mal_Bakry11_uncertainty)
    par.fill.mal[, i, 4] = parSapply(clus1, conc.mal, muNq_mal_Bakry11_uncertainty)
    par.fill.mal[, i, 5] = parSapply(clus1, conc.mal, muPq_mal_Halstead_uncertainty)
    
  #fill equilibrium snail pop estimates  
    #individual parameters
    Neq.fill.mal[, i, 1] = parSapply(clus1, conc.mal, r0.In, f.pi_Cq = piC.tch92_mal_unc)[1,]
    Neq.fill.mal[, i, 2] = parSapply(clus1, conc.mal, r0.In, f.pi_Mq = piM.tch91_mal_unc)[1,]
    Neq.fill.mal[, i, 3] = parSapply(clus1, conc.mal, r0.In, f.f_Nq = fNq_mal_Bakry11_uncertainty)[1,]
    Neq.fill.mal[, i, 4] = parSapply(clus1, conc.mal, r0.In, f.mu_Nq = muNq_mal_Bakry11_uncertainty)[1,]
    Neq.fill.mal[, i, 5] = parSapply(clus1, conc.mal, r0.In, f.mu_Pq = muPq_mal_Halstead_uncertainty)[1,]
    
    #parameter pairs and triplets
    #larval toxicity to both miracidia and cercariae
    Neq.fill.mal[, i, 6] = parSapply(clus1, conc.mal, r0.In,
                                    f.pi_Mq = piM.tch91_mal_unc,
                                    f.pi_Cq = piC.tch92_mal_unc)[1,]
    #both snail effects (mortality and reproduction)
    Neq.fill.mal[, i, 7] = parSapply(clus1, conc.mal, r0.In,
                                    f.mu_Nq = muNq_mal_Bakry11_uncertainty,
                                    f.f_Nq = fNq_mal_Bakry11_uncertainty)[1,]
    #larval mortality with predator mortality
    Neq.fill.mal[, i, 8] = parSapply(clus1, conc.mal, r0.In,
                                    f.pi_Mq = piM.tch91_mal_unc,
                                    f.pi_Cq = piC.tch92_mal_unc,
                                    f.mu_Pq = muPq_mal_Halstead_uncertainty)[1,]
    #snail effects with predator mortality
    Neq.fill.mal[, i, 9] = parSapply(clus1, conc.mal, r0.In,
                                    f.mu_Nq = muNq_mal_Bakry11_uncertainty,
                                    f.f_Nq = fNq_mal_Bakry11_uncertainty,
                                    f.mu_Pq = muPq_mal_Halstead_uncertainty)[1,]
    #Full simulation w/all parameters
    Neq.fill.mal[, i, 10] = parSapply(clus1, conc.mal, r0.In,
                                     f.mu_Nq = muNq_mal_Bakry11_uncertainty,
                                     f.f_Nq = fNq_mal_Bakry11_uncertainty,
                                     f.pi_Mq = piM.tch91_mal_unc,
                                     f.pi_Cq = piC.tch92_mal_unc,
                                     f.mu_Pq = muPq_mal_Halstead_uncertainty)[1,]
    #Full simulation w/all parameters EXCEPT snail reproduction since function is questionable
    Neq.fill.mal[, i, 11] = parSapply(clus1, conc.mal, r0.In,
                                     f.mu_Nq = muNq_mal_Bakry11_uncertainty,
                                     f.pi_Mq = piM.tch91_mal_unc,
                                     f.pi_Cq = piC.tch92_mal_unc,
                                     f.mu_Pq = muPq_mal_Halstead_uncertainty)[1,]
  } #End function

stopCluster(clus1)
########## Post process ############ 
#Get means and sds of all simulations  #########
#mean r0
conc.mal.means.r0 = matrix(nrow = length(conc.mal), ncol = length(parfx)+6)
  for(i in 1:(length(parfx)+6)){
    conc.mal.means.r0[,i] = rowMeans(r0.fill.mal[ , , i])
  }

#sd r0
conc.mal.sds.r0 = matrix(nrow = length(conc.mal), ncol = length(parfx)+6)
  for(i in 1:(length(parfx)+6)){
    conc.mal.sds.r0[,i] = fBasics::rowSds(r0.fill.mal[ , , i])
  }

#mean Neq
conc.mal.means.Neq = matrix(nrow = length(conc.mal), ncol = length(parfx)+6)
  for(i in 1:(length(parfx)+6)){
    conc.mal.means.Neq[,i] = rowMeans(Neq.fill.mal[ , , i])
  }

#sd Neq
conc.mal.sds.Neq = matrix(nrow = length(conc.mal), ncol = length(parfx)+6)
  for(i in 1:(length(parfx)+6)){
    conc.mal.sds.Neq[,i] = fBasics::rowSds(Neq.fill.mal[ , , i])
  }

#mean parameter values
conc.mal.means.par = matrix(nrow = length(conc.mal), ncol = length(parfx))
  for(i in 1:length(parfx)){
    conc.mal.means.par[,i] = rowMeans(par.fill.mal[ , , i])
  }
  colnames(conc.mal.means.par) = c('piC_tch92', 'piM_tch92', 'f_Nq_bak11', 'mu_Nq_bak11', 'mu_Pq_hal15')
