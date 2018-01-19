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

#Run simulations of glyphosate concentrations from 0 - 2000, start with individual functions then combine ################
conc.gly = seq(0, 3700*2, length.out = 200)  #Concentration range to test
nsims = 100       #Number of simulations to run
parfx = c(piC.meta_atr_unc, muNq_atr_Bakry12_uncertainty, phi_Nq_atr_baxrohr, #Functions corresponding to affected parameters
          piC.grg08_atr_unc, piC.kop_atr_unc2, piC.atr.rohr08.lin, ons.munq.atr)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.atr, 'uniroot.all', 'rdrm', 'LL.2'))

r0.atr.fill = array(data = NA, dim = c(length(conc.atr), nsims, length(parfx)+10))

clusterSetRNGStream(cl = clus1, iseed = 43093) #set cluster seed for reproducibility

set.seed(0)

  for(i in 1:nsims){
  #Fill r0 estimates ######
    #individual parameters
    r0.atr.fill[, i, 1] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.meta_atr_unc)[3,]
    r0.atr.fill[, i, 2] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.grg08_atr_unc)[3,]
    r0.atr.fill[, i, 3] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.kop_atr_unc2)[3,]
    r0.atr.fill[, i, 4] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.atr.rohr08.lin)[3,]
    r0.atr.fill[, i, 5] = parSapply(clus1, conc.atr, r0.He, f.mu_Nq = muNq_atr_Bakry12_uncertainty)[3,]
    r0.atr.fill[, i, 6] = parSapply(clus1, conc.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.atr.fill[, i, 7] = parSapply(clus1, conc.atr, r0.He, f.mu_Nq = ons.munq.atr)[3,]
  #pairwise parameter combinations
    #meta cercarial mortality and rohr analysis of Baxter carrying capcity  
    r0.atr.fill[, i, 8] = parSapply(clus1, conc.atr, r0.He,
                                    f.pi_Cq = piC.meta_atr_unc,
                                    f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    #meta cercarial mortality and bakry 2012 snail mortality 
    r0.atr.fill[, i, 9] = parSapply(clus1, conc.atr, r0.He,
                                    f.mu_Nq = muNq_atr_Bakry12_uncertainty,
                                    f.pi_Cq = piC.meta_atr_unc)[3,]
    #Rohr08 cercarial mortality and rohr analysis of Baxter carrying capcity 
    r0.atr.fill[, i, 10] = parSapply(clus1, conc.atr, r0.He,
                                    f.pi_Cq = piC.atr.rohr08.lin,
                                    f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    #Rohr cercarial mortality and bakry 2012 snail mortality 
    r0.atr.fill[, i, 11] = parSapply(clus1, conc.atr, r0.He,
                                    f.mu_Nq = muNq_atr_Bakry12_uncertainty,
                                    f.pi_Cq = piC.atr.rohr08.lin)[3,]
    #Bakry 2012 snail mortality and rohr analysis of Baxter carrying capcity 
    r0.atr.fill[, i, 12] = parSapply(clus1, conc.atr, r0.He,
                                    f.mu_Nq = muNq_atr_Bakry12_uncertainty,
                                    f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    #Rohr cercarial mortality and Omran&Salama snail mortality
    r0.atr.fill[, i, 13] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = ons.munq.atr,
                                     f.pi_Cq = piC.atr.rohr08.lin)[3,]
    #Omran&Salama snail mortality and rohr analysis of Baxter carrying capcity 
    r0.atr.fill[, i, 14] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = ons.munq.atr,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    
    #triplicate parameter combinations
    r0.atr.fill[, i, 15] = parSapply(clus1, conc.atr, r0.He,
                                              f.mu_Nq = muNq_atr_Bakry12_uncertainty,
                                              f.pi_Cq = piC.meta_atr_unc,
                                              f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.atr.fill[, i, 16] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = muNq_atr_Bakry12_uncertainty,
                                     f.pi_Cq = piC.atr.rohr08.lin,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.atr.fill[, i, 17] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = ons.munq.atr,
                                     f.pi_Cq = piC.atr.rohr08.lin,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
  #Fill equilibrium snail pop estimates #######
    #individual parameters
    Neq.atr.fill[, i, 1] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.meta_atr_unc)[1,]
    Neq.atr.fill[, i, 2] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.grg08_atr_unc)[1,]
    Neq.atr.fill[, i, 3] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.kop_atr_unc2)[1,]
    Neq.atr.fill[, i, 4] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.atr.rohr08.lin)[1,]
    Neq.atr.fill[, i, 5] = parSapply(clus1, conc.atr, r0.He, f.mu_Nq = muNq_atr_Bakry12_uncertainty)[1,]
    Neq.atr.fill[, i, 6] = parSapply(clus1, conc.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
    Neq.atr.fill[, i, 7] = parSapply(clus1, conc.atr, r0.He, f.mu_Nq = ons.munq.atr)[1,]
    #pairwise parameter combinations
    #meta cercarial mortality and rohr analysis of Baxter carrying capcity  
    Neq.atr.fill[, i, 8] = parSapply(clus1, conc.atr, r0.He,
                                    f.pi_Cq = piC.meta_atr_unc,
                                    f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
    #meta cercarial mortality and bakry 2012 snail mortality 
    Neq.atr.fill[, i, 9] = parSapply(clus1, conc.atr, r0.He,
                                    f.mu_Nq = muNq_atr_Bakry12_uncertainty,
                                    f.pi_Cq = piC.meta_atr_unc)[1,]
    #Rohr08 cercarial mortality and rohr analysis of Baxter carrying capcity 
    Neq.atr.fill[, i, 10] = parSapply(clus1, conc.atr, r0.He,
                                     f.pi_Cq = piC.atr.rohr08.lin,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
    #Rohr cercarial mortality and bakry 2012 snail mortality 
    Neq.atr.fill[, i, 11] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = muNq_atr_Bakry12_uncertainty,
                                     f.pi_Cq = piC.atr.rohr08.lin)[1,]
    #Bakry 2012 snail mortality and rohr analysis of Baxter carrying capcity 
    Neq.atr.fill[, i, 12] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = muNq_atr_Bakry12_uncertainty,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
    #Rohr cercarial mortality and Omran&Salama snail mortality
    Neq.atr.fill[, i, 13] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = ons.munq.atr,
                                     f.pi_Cq = piC.atr.rohr08.lin)[1,]
    #Omran&Salama snail mortality and rohr analysis of Baxter carrying capcity 
    Neq.atr.fill[, i, 14] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = ons.munq.atr,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
    
    #triplicate parameter combinations
    Neq.atr.fill[, i, 15] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = muNq_atr_Bakry12_uncertainty,
                                     f.pi_Cq = piC.meta_atr_unc,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
    Neq.atr.fill[, i, 16] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = muNq_atr_Bakry12_uncertainty,
                                     f.pi_Cq = piC.atr.rohr08.lin,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
    Neq.atr.fill[, i, 17] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = ons.munq.atr,
                                     f.pi_Cq = piC.atr.rohr08.lin,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
  #Store parameter values   #######
    par.atr.fill[, i, 1] = parSapply(clus1, conc.atr, piC.meta_atr_unc)
    par.atr.fill[, i, 2] = parSapply(clus1, conc.atr, piC.grg08_atr_unc)
    par.atr.fill[, i, 3] = parSapply(clus1, conc.atr, piC.kop_atr_unc2)
    par.atr.fill[, i, 4] = parSapply(clus1, conc.atr, piC.atr.rohr08.lin)
    par.atr.fill[, i, 5] = parSapply(clus1, conc.atr, muNq_atr_Bakry12_uncertainty)
    par.atr.fill[, i, 6] = parSapply(clus1, conc.atr, phi_Nq_atr_baxrohr.no30)
    par.atr.fill[, i, 7] = parSapply(clus1, conc.atr, ons.munq.atr)
    
    print(i)
  }

stopCluster(clus1)
 
########## Post process ############
#mean and sd of r0 and N eq simulations #########
#mean r0
conc.atr.means.r0 = matrix(nrow = length(conc.atr), ncol = length(parfx)+10)
  for(i in 1:(length(parfx)+10)){
    conc.atr.means.r0[,i] = rowMeans(r0.atr.fill[ , , i])
  }
#sd r0
conc.atr.sds.r0 = matrix(nrow = length(conc.atr), ncol = length(parfx)+10)
  for(i in 1:(length(parfx)+10)){
    conc.atr.sds.r0[,i] = fBasics::rowSds(r0.atr.fill[ , , i])
  }
#mean Neq
conc.atr.means.Neq = matrix(nrow = length(conc.atr), ncol = length(parfx)+10)
  for(i in 1:(length(parfx)+10)){
    conc.atr.means.Neq[,i] = rowMeans(Neq.atr.fill[ , , i])
  }
#sd eq
conc.atr.sds.Neq = matrix(nrow = length(conc.atr), ncol = length(parfx)+10)
  for(i in 1:(length(parfx)+10)){
    conc.atr.sds.Neq[,i] = fBasics::rowSds(Neq.atr.fill[ , , i])
  }
#mean parameter values
conc.atr.means.par = matrix(nrow = length(conc.atr), ncol = length(parfx))
for(i in 1:length(parfx)){
  conc.atr.means.par[,i] = rowMeans(par.atr.fill[ , , i])
}
  colnames(conc.atr.means.par) = c('piC_meta', 'piC_grg', 'piC_kop', 'piC_rohr',
                                   'muN_bak', 'phiN_bax', 'muN_OnS')