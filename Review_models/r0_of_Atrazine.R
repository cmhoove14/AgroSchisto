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

source('Review_models/r0_of_q.R')

library(parallel)
library(reshape2)

keep.fin.atr = c(keep.bak12, keep.grg08, keep.kop06.beq, keep.meta.piC, keep.baxrohr, keep.atr.rohr08,
                 'r0.He', 'r0.fix', 'parameters', 'nil0', 'nil1', 'keep.fin.atr')

rm(list = setdiff(ls(), keep.fin.atr))
dev.off()

no.cores = detectCores() - 1

#Fill data frame with point estimates of R0 driven by data in each study ################
r0.atr.fix = data.frame(atr = c(atra.df$atra, 
                                kopatr.auc$atr,
                                grgatr.auc$atr,
                                atr$conc,
                                muN.bak$atr*1000),
                        par = c(rep('phi_N', length(atra.df$atra)), 
                                rep('piC', length(kopatr.auc$atr)),
                                rep('piC', length(grgatr.auc$atr)),
                                rep('piC', length(atr$conc)),
                                rep('muN', length(muN.bak$atr))),
                        study = c(rep('Baxter11', length(atra.df$atra)), 
                                  rep('Koprivnikar06', length(kopatr.auc$atr)),
                                  rep('Griggs08', length(grgatr.auc$atr)),
                                  rep('Rohr08', length(atr$conc)),
                                  rep('Bakry12', length(muN.bak$atr))),
                        r0 = c(r0.fix(phi_Nqx = atra.df$growthrate[1] / atra.df$growthrate[1])[3],
                               r0.fix(phi_Nqx = atra.df$growthrate[2] / atra.df$growthrate[1])[3],
                               r0.fix(phi_Nqx = atra.df$growthrate[3] / atra.df$growthrate[1])[3],
                               r0.fix(phi_Nqx = atra.df$growthrate[4] / atra.df$growthrate[1])[3],
                               r0.fix(phi_Nqx = atra.df$growthrate[5] / atra.df$growthrate[1])[3],
                               r0.fix(pi_Cqx = kopatr.auc$piC[1])[3],
                               r0.fix(pi_Cqx = kopatr.auc$piC[2])[3],
                               r0.fix(pi_Cqx = kopatr.auc$piC[3])[3],
                               r0.fix(pi_Cqx = grgatr.auc$piC[1])[3],
                               r0.fix(pi_Cqx = grgatr.auc$piC[2])[3],
                               r0.fix(pi_Cqx = grgatr.auc$piC[3])[3],
                               r0.fix(pi_Cqx = (atr$surv[1] / atr$surv[1]))[3],
                               r0.fix(pi_Cqx = (atr$surv[2] / atr$surv[1]))[3],
                               r0.fix(pi_Cqx = (atr$surv[3] / atr$surv[1]))[3],
                               r0.fix(pi_Cqx = (atr$surv[4] / atr$surv[1]))[3],
                               r0.fix(pi_Cqx = (atr$surv[5] / atr$surv[1]))[3],
                               r0.fix(pi_Cqx = (atr$surv[6] / atr$surv[1]))[3],
                               r0.fix(mu_Nqx = muN.bak$mort[1])[3],
                               r0.fix(mu_Nqx = muN.bak$mort[2])[3],
                               r0.fix(mu_Nqx = muN.bak$mort[3])[3]))

r0.atr.fix.n0 = subset(r0.atr.fix, atr != 0)  

#Run simulations of atrazine concentrations from 0 - 200, start with individual functions then combine ################
conc.atr = c(0:500)  #Concentration range to test
nsims = 30       #Number of simulations to run
parfx = c(piC.meta_atr_unc, muNq_atr_bak12, phi_Nq_atr_baxrohr, #Functions corresponding to affected parameters
          piC.grg08_atr_unc, piC.kop_atr_unc2, piC.atr.rohr08.lin)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.atr, 'uniroot.all', 'rdrm', 'LL.2'))

r0.atr.fill = array(data = NA, dim = c(length(conc.atr), nsims, length(parfx)+7))
par.atr.fill = array(data = NA, dim = c(length(conc.atr), nsims, length(parfx)))
Neq.atr.fill = array(data = NA, dim = c(length(conc.atr), nsims, length(parfx)+7))

set.seed(0)

  for(i in 1:nsims){
  #Fill r0 estimates
    #individual parameters
    r0.atr.fill[, i, 1] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.meta_atr_unc)[3,]
    r0.atr.fill[, i, 2] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.grg08_atr_unc)[3,]
    r0.atr.fill[, i, 3] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.kop_atr_unc2)[3,]
    r0.atr.fill[, i, 4] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.atr.rohr08.lin)[3,]
    r0.atr.fill[, i, 5] = parSapply(clus1, conc.atr, r0.He, f.mu_Nq = muNq_atr_bak12)[3,]
    r0.atr.fill[, i, 6] = parSapply(clus1, conc.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    #pairwise parameter combinations
    r0.atr.fill[, i, 7] = parSapply(clus1, conc.atr, r0.He,
                                    f.pi_Cq = piC.meta_atr_unc,
                                    f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.atr.fill[, i, 8] = parSapply(clus1, conc.atr, r0.He,
                                    f.mu_Nq = muNq_atr_bak12,
                                    f.pi_Cq = piC.meta_atr_unc)[3,]
    r0.atr.fill[, i, 9] = parSapply(clus1, conc.atr, r0.He,
                                    f.pi_Cq = piC.atr.rohr08.lin,
                                    f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.atr.fill[, i, 10] = parSapply(clus1, conc.atr, r0.He,
                                    f.mu_Nq = muNq_atr_bak12,
                                    f.pi_Cq = piC.atr.rohr08.lin)[3,]
    r0.atr.fill[, i, 11] = parSapply(clus1, conc.atr, r0.He,
                                    f.mu_Nq = muNq_atr_bak12,
                                    f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    #triplicate parameter combinations
    r0.atr.fill[, i, 12] = parSapply(clus1, conc.atr, r0.He,
                                              f.mu_Nq = muNq_atr_bak12,
                                              f.pi_Cq = piC.meta_atr_unc,
                                              f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.atr.fill[, i, 13] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = muNq_atr_bak12,
                                     f.pi_Cq = piC.atr.rohr08.lin,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
  #Fill equilibrium snail pop estimates  
    #individual parameters
    Neq.atr.fill[, i, 1] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.meta_atr_unc)[1,]
    Neq.atr.fill[, i, 2] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.grg08_atr_unc)[1,]
    Neq.atr.fill[, i, 3] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.kop_atr_unc2)[1,]
    Neq.atr.fill[, i, 4] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.atr.rohr08.lin)[1,]
    Neq.atr.fill[, i, 5] = parSapply(clus1, conc.atr, r0.He, f.mu_Nq = muNq_atr_bak12)[1,]
    Neq.atr.fill[, i, 6] = parSapply(clus1, conc.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
    #pairwise parameter combinations
    Neq.atr.fill[, i, 7] = parSapply(clus1, conc.atr, r0.He,
                                    f.pi_Cq = piC.meta_atr_unc,
                                    f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
    Neq.atr.fill[, i, 8] = parSapply(clus1, conc.atr, r0.He,
                                    f.mu_Nq = muNq_atr_bak12,
                                    f.pi_Cq = piC.meta_atr_unc)[1,]
    Neq.atr.fill[, i, 9] = parSapply(clus1, conc.atr, r0.He,
                                    f.pi_Cq = piC.atr.rohr08.lin,
                                    f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
    Neq.atr.fill[, i, 10] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = muNq_atr_bak12,
                                     f.pi_Cq = piC.atr.rohr08.lin)[1,]
    Neq.atr.fill[, i, 11] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = muNq_atr_bak12,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
    #triplicate parameter combinations
    Neq.atr.fill[, i, 12] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = muNq_atr_bak12,
                                     f.pi_Cq = piC.meta_atr_unc,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
    Neq.atr.fill[, i, 13] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = muNq_atr_bak12,
                                     f.pi_Cq = piC.atr.rohr08.lin,
                                     f.phi_Nq = phi_Nq_atr_baxrohr.no30)[1,]
  #Store parameter values  
    par.atr.fill[, i, 1] = parSapply(clus1, conc.atr, piC.meta_atr_unc)
    par.atr.fill[, i, 2] = parSapply(clus1, conc.atr, piC.grg08_atr_unc)
    par.atr.fill[, i, 3] = parSapply(clus1, conc.atr, piC.kop_atr_unc2)
    par.atr.fill[, i, 4] = parSapply(clus1, conc.atr, piC.atr.rohr08.lin)
    par.atr.fill[, i, 5] = parSapply(clus1, conc.atr, muNq_atr_bak12)
    par.atr.fill[, i, 6] = parSapply(clus1, conc.atr, phi_Nq_atr_baxrohr.no30)
  }

stopCluster(clus1)
 
########## Post process ############
#mean and sd of r0 and N eq simulations #########
conc.atr.means.r0 = matrix(nrow = length(conc.atr), ncol = length(parfx)+7)
  for(i in 1:(length(parfx)+7)){
    conc.atr.means.r0[,i] = rowMeans(r0.atr.fill[ , , i])
  }

conc.atr.sds.r0 = matrix(nrow = length(conc.atr), ncol = length(parfx)+7)
  for(i in 1:(length(parfx)+7)){
    conc.atr.sds.r0[,i] = fBasics::rowSds(r0.atr.fill[ , , i])
  }

conc.atr.means.Neq = matrix(nrow = length(conc.atr), ncol = length(parfx)+7)
  for(i in 1:(length(parfx)+7)){
    conc.atr.means.Neq[,i] = rowMeans(Neq.atr.fill[ , , i])
  }

conc.atr.sds.Neq = matrix(nrow = length(conc.atr), ncol = length(parfx)+7)
  for(i in 1:(length(parfx)+7)){
    conc.atr.sds.Neq[,i] = fBasics::rowSds(Neq.atr.fill[ , , i])
  }

conc.atr.means.par = matrix(nrow = length(conc.atr), ncol = length(parfx))
for(i in 1:length(parfx)){
  conc.atr.means.par[,i] = rowMeans(par.atr.fill[ , , i])
}
  colnames(conc.atr.means.par) = c('piC_meta', 'piC_grg', 'piC_kop', 'piC_rohr',
                                   'muN_bak', 'phiN_bax')

#Loess regressions for each parameter #####################
d1.atr = as.data.frame(r0.atr.fill[ , , 1])
  d1.atr2 = reshape(d1.atr, varying = c(colnames(d1.atr)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d1.atr)), direction = 'long')
#Fit loess to output of all simulations
  l1.atr = loess(r0 ~ id, data = d1.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l1.atr.2 = loess(r0.atr.fill[ , 1, 1] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l1.atr.df = data.frame(atr = conc.atr,
                     pred = predict(l1.atr, newdata=data.frame(id = conc.atr)),
                     se = predict(l1.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                     pred.2 = predict(l1.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                     se.2 = predict(l1.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)
  
#*********************************************************************  
  
d2.atr = as.data.frame(r0.atr.fill[ , , 2])
  d2.atr2 = reshape(d2.atr, varying = c(colnames(d2.atr)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d2.atr)), direction = 'long')
#Fit loess to output of all simulations
  l2.atr = loess(r0 ~ id, data = d2.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l2.atr.2 = loess(r0.atr.fill[ , 1, 2] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l2.atr.df = data.frame(atr = conc.atr,
                     pred = predict(l2.atr, newdata=data.frame(id = conc.atr)),
                     se = predict(l2.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                     pred.2 = predict(l2.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                     se.2 = predict(l2.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)

#*********************************************************************

d3.atr = as.data.frame(r0.atr.fill[ , , 3])
  d3.atr2 = reshape(d3.atr, varying = c(colnames(d3.atr)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d3.atr)), direction = 'long')
#Fit loess to output of all simulations
  l3.atr = loess(r0 ~ id, data = d3.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l3.atr.2 = loess(r0.atr.fill[ , 1, 3] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l3.atr.df = data.frame(atr = conc.atr,
                     pred = predict(l3.atr, newdata=data.frame(id = conc.atr)),
                     se = predict(l3.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                     pred.2 = predict(l3.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                     se.2 = predict(l3.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)

#*********************************************************************
  
d4.atr = as.data.frame(r0.atr.fill[ , , 4])
  d4.atr2 = reshape(d4.atr, varying = c(colnames(d4.atr)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d4.atr)), direction = 'long')
  
#Fit loess to output of all simulations
  l4.atr = loess(r0 ~ id, data = d4.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l4.atr.2 = loess(r0.atr.fill[ , 1, 4] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l4.atr.df = data.frame(atr = conc.atr,
                     pred = predict(l4.atr, newdata=data.frame(id = conc.atr)),
                     se = predict(l4.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                     pred.2 = predict(l4.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                     se.2 = predict(l4.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)
  
#*********************************************************************  
  
d5.atr = as.data.frame(r0.atr.fill[ , , 5])
  d5.atr2 = reshape(d5.atr, varying = c(colnames(d5.atr)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d5.atr)), direction = 'long')
  
#Fit loess to output of all simulations
  l5.atr = loess(r0 ~ id, data = d5.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l5.atr.2 = loess(r0.atr.fill[ , 1, 5] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l5.atr.df = data.frame(atr = conc.atr,
                     pred = predict(l5.atr, newdata=data.frame(id = conc.atr)),
                     se = predict(l5.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                     pred.2 = predict(l5.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                     se.2 = predict(l5.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)
  
#*********************************************************************
  
d6.atr = as.data.frame(r0.atr.fill[ , , 6])
  d6.atr2 = reshape(d6.atr, varying = c(colnames(d6.atr)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d6.atr)), direction = 'long')

#Fit loess to output of all simulations
  l6.atr = loess(r0 ~ id, data = d6.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l6.atr.2 = loess(r0.atr.fill[ , 1, 6] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l6.atr.df = data.frame(atr = conc.atr,
                     pred = predict(l6.atr, newdata=data.frame(id = conc.atr)),
                     se = predict(l6.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                     pred.2 = predict(l6.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                     se.2 = predict(l6.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)
  
#*********************************************************************  
  
d7.atr = as.data.frame(r0.atr.fill[ , , 7])
  d7.atr2 = reshape(d7.atr, varying = c(colnames(d7.atr)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d7.atr)), direction = 'long')
#Fit loess to output of all simulations
  l7.atr = loess(r0 ~ id, data = d7.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l7.atr.2 = loess(r0.atr.fill[ , 1, 7] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l7.atr.df = data.frame(atr = conc.atr,
                         pred = predict(l7.atr, newdata=data.frame(id = conc.atr)),
                         se = predict(l7.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                         pred.2 = predict(l7.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                         se.2 = predict(l7.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)
  
#*********************************************************************  
  
d8.atr = as.data.frame(r0.atr.fill[ , , 8])
  d8.atr2 = reshape(d8.atr, varying = c(colnames(d8.atr)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d8.atr)), direction = 'long')
#Fit loess to output of all simulations
  l8.atr = loess(r0 ~ id, data = d8.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l8.atr.2 = loess(r0.atr.fill[ , 1, 8] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l8.atr.df = data.frame(atr = conc.atr,
                         pred = predict(l8.atr, newdata=data.frame(id = conc.atr)),
                         se = predict(l8.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                         pred.2 = predict(l8.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                         se.2 = predict(l8.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)
  
#*********************************************************************  
  
d9.atr = as.data.frame(r0.atr.fill[ , , 9])
  d9.atr2 = reshape(d9.atr, varying = c(colnames(d9.atr)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d9.atr)), direction = 'long')
#Fit loess to output of all simulations
  l9.atr = loess(r0 ~ id, data = d9.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l9.atr.2 = loess(r0.atr.fill[ , 1, 9] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l9.atr.df = data.frame(atr = conc.atr,
                         pred = predict(l9.atr, newdata=data.frame(id = conc.atr)),
                         se = predict(l9.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                         pred.2 = predict(l9.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                         se.2 = predict(l9.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)

#*********************************************************************  
  
d10.atr = as.data.frame(r0.atr.fill[ , , 10])
  d10.atr2 = reshape(d10.atr, varying = c(colnames(d10.atr)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d10.atr)), direction = 'long')
#Fit loess to output of all simulations
  l10.atr = loess(r0 ~ id, data = d10.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l10.atr.2 = loess(r0.atr.fill[ , 1, 10] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l10.atr.df = data.frame(atr = conc.atr,
                         pred = predict(l10.atr, newdata=data.frame(id = conc.atr)),
                         se = predict(l10.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                         pred.2 = predict(l10.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                         se.2 = predict(l10.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)

#*********************************************************************  
  
d11.atr = as.data.frame(r0.atr.fill[ , , 11])
  d11.atr2 = reshape(d11.atr, varying = c(colnames(d11.atr)),
                     v.names = 'r0', timevar = 'sim', 
                     times = c(colnames(d11.atr)), direction = 'long')
#Fit loess to output of all simulations
  l11.atr = loess(r0 ~ id, data = d11.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l11.atr.2 = loess(r0.atr.fill[ , 1, 11] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l11.atr.df = data.frame(atr = conc.atr,
                          pred = predict(l11.atr, newdata=data.frame(id = conc.atr)),
                          se = predict(l11.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                          pred.2 = predict(l11.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                          se.2 = predict(l11.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)
  
#*********************************************************************  
  
d12.atr = as.data.frame(r0.atr.fill[ , , 12])
  d12.atr2 = reshape(d12.atr, varying = c(colnames(d12.atr)),
                     v.names = 'r0', timevar = 'sim', 
                     times = c(colnames(d12.atr)), direction = 'long')
#Fit loess to output of all simulations
  l12.atr = loess(r0 ~ id, data = d12.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l12.atr.2 = loess(r0.atr.fill[ , 1, 12] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l12.atr.df = data.frame(atr = conc.atr,
                          pred = predict(l12.atr, newdata=data.frame(id = conc.atr)),
                          se = predict(l12.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                          pred.2 = predict(l12.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                          se.2 = predict(l12.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)

#*********************************************************************  
  
d13.atr = as.data.frame(r0.atr.fill[ , , 13])
  d13.atr2 = reshape(d13.atr, varying = c(colnames(d13.atr)),
                     v.names = 'r0', timevar = 'sim', 
                     times = c(colnames(d13.atr)), direction = 'long')
#Fit loess to output of all simulations
  l13.atr = loess(r0 ~ id, data = d13.atr2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l13.atr.2 = loess(r0.atr.fill[ , 1, 13] ~ conc.atr, degree = 2, span = 0.75, family = 'gaussian')
  
  l13.atr.df = data.frame(atr = conc.atr,
                          pred = predict(l13.atr, newdata=data.frame(id = conc.atr)),
                          se = predict(l13.atr, newdata=data.frame(id = conc.atr), se = TRUE)$se.fit,
                          pred.2 = predict(l13.atr.2, newdata=data.frame(conc.atr = conc.atr)),
                          se.2 = predict(l13.atr.2, newdata=data.frame(conc.atr = conc.atr), se = TRUE)$se.fit)
  
#Plot results ##############################
#mean of all sims  
plot(conc.atr, conc.atr.means.r0[,1], pch = 17, cex=0.5, col = 2, xlab = 'Atrazine (ppb)', ylab = expression(paste(R[0], '(atrazine)')),
     ylim = c(2.5,6.5), xlim = c(0, max(conc.atr)))
  for(i in 2:(length(parfx)+4)){
    points(conc.atr, conc.atr.means.r0[,i], pch = 17, cex=0.5, col = i+1)
  }
  
plot(conc.atr, conc.atr.means.Neq[,1], pch = 17, cex=0.5, col = 2, xlab = 'Atrazine (ppb)', ylab = expression(N[eq]),
     xlim = c(0, max(conc.atr)), ylim = c(8000, 28000))
  for(i in 2:(length(parfx)+4)){
    points(conc.atr, conc.atr.means.Neq[,i], pch = 17, cex=0.5, col = i+1)
  }  