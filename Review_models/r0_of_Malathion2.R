#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Run this and skip to plotting functions unless manipulating parameters, simulations, etc. #############
  load('Review_models/r0_of_malathion2_ws.RData')

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
conc.mal = c(c(0:999), seq(1000, 200000, length.out = 1001))  #Concentration range to test
nsims = 10         #Number of simulations to run
parfx = c(piC.tch92_mal_unc, piM.tch91_mal_unc, #Functions corresponding to affected parameters
          fNq_mal_Bakry11_uncertainty, muNq_mal_Bakry11_uncertainty, muPq_mal_Halstead_uncertainty)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.mal, 'uniroot.all', 'rdrm', 'LL.2'))

r0.fill.mal = array(data = NA, dim = c(length(conc.mal), nsims, length(parfx)+2))
par.fill.mal = array(data = NA, dim = c(length(conc.mal), nsims, length(parfx)))
Neq.fill.mal = array(data = NA, dim = c(length(conc.mal), nsims, length(parfx)+2))

set.seed(0)

  for(i in 1:nsims){
  #Fill r0 estimates  
    r0.fill.mal[, i, 1] = parSapply(clus1, conc.mal, r0.In, f.pi_Cq = piC.tch92_mal_unc)[3,]
    r0.fill.mal[, i, 2] = parSapply(clus1, conc.mal, r0.In, f.pi_Mq = piM.tch91_mal_unc)[3,]
    r0.fill.mal[, i, 3] = parSapply(clus1, conc.mal, r0.In, f.f_Nq = fNq_mal_Bakry11_uncertainty)[3,]
    r0.fill.mal[, i, 4] = parSapply(clus1, conc.mal, r0.In, f.mu_Nq = muNq_mal_Bakry11_uncertainty)[3,]
    r0.fill.mal[, i, 5] = parSapply(clus1, conc.mal, r0.In, f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    r0.fill.mal[, i, 6] = parSapply(clus1, conc.mal, r0.In,
                                              f.mu_Nq = muNq_mal_Bakry11_uncertainty,
                                              f.f_Nq = fNq_mal_Bakry11_uncertainty,
                                              f.pi_Mq = piM.tch91_mal_unc,
                                              f.pi_Cq = piC.tch92_mal_unc,
                                              f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    r0.fill.mal[, i, 7] = parSapply(clus1, conc.mal, r0.In,
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
    Neq.fill.mal[, i, 1] = parSapply(clus1, conc.mal, r0.In, f.pi_Cq = piC.tch92_mal_unc)[1,]
    Neq.fill.mal[, i, 2] = parSapply(clus1, conc.mal, r0.In, f.pi_Mq = piM.tch91_mal_unc)[1,]
    Neq.fill.mal[, i, 3] = parSapply(clus1, conc.mal, r0.In, f.f_Nq = fNq_mal_Bakry11_uncertainty)[1,]
    Neq.fill.mal[, i, 4] = parSapply(clus1, conc.mal, r0.In, f.mu_Nq = muNq_mal_Bakry11_uncertainty)[1,]
    Neq.fill.mal[, i, 5] = parSapply(clus1, conc.mal, r0.In, f.mu_Pq = muPq_mal_Halstead_uncertainty)[1,]
    Neq.fill.mal[, i, 6] = parSapply(clus1, conc.mal, r0.In,
                                    f.mu_Nq = muNq_mal_Bakry11_uncertainty,
                                    f.f_Nq = fNq_mal_Bakry11_uncertainty,
                                    f.pi_Mq = piM.tch91_mal_unc,
                                    f.pi_Cq = piC.tch92_mal_unc,
                                    f.mu_Pq = muPq_mal_Halstead_uncertainty)[1,]
    Neq.fill.mal[, i, 7] = parSapply(clus1, conc.mal, r0.In,
                                    f.pi_Mq = piM.tch91_mal_unc,
                                    f.pi_Cq = piC.tch92_mal_unc,
                                    f.mu_Pq = muPq_mal_Halstead_uncertainty)[1,]
  }

stopCluster(clus1)
########## Post process ############ 
#Get means of all simulations  #########
conc.mal.means.r0 = matrix(nrow = length(conc.mal), ncol = length(parfx)+2)
  for(i in 1:(length(parfx)+2)){
    conc.mal.means.r0[,i] = rowMeans(r0.fill.mal[ , , i])
  }

conc.mal.means.Neq = matrix(nrow = length(conc.mal), ncol = length(parfx)+2)
  for(i in 1:(length(parfx)+2)){
    conc.mal.means.Neq[,i] = rowMeans(Neq.fill.mal[ , , i])
  }


#Format simulations into long format for scatterplots; fit loess regressions for smoothed lines #######
d1.mal = as.data.frame(r0.fill.mal[ , , 1])
  d1.mal2 = reshape(d1.mal, varying = c(colnames(d1.mal)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d1.mal)), direction = 'long')
#Fit loess to output of all simulations
  l1.mal = loess(r0 ~ id, data = d1.mal2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l1.mal.2 = loess(r0.fill.mal[ , 1, 1] ~ conc.mal, degree = 2, span = 0.75, family = 'gaussian')

  l1.mal.df = data.frame(mal = conc.mal,
                     pred = predict(l1.mal, newdata=data.frame(id = conc.mal)),
                     se = predict(l1.mal, newdata=data.frame(id = conc.mal), se = TRUE)$se.fit,
                     pred.2 = predict(l1.mal.2, newdata=data.frame(conc.mal = conc.mal)),
                     se.2 = predict(l1.mal.2, newdata=data.frame(conc.mal = conc.mal), se = TRUE)$se.fit)

#*********************************************************************  

d2.mal = as.data.frame(r0.fill.mal[ , , 2])
  d2.mal2 = reshape(d2.mal, varying = c(colnames(d2.mal)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d2.mal)), direction = 'long')
#Fit loess to output of all simulations
  l2.mal = loess(r0 ~ id, data = d2.mal2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l2.mal.2 = loess(r0.fill.mal[ , 1, 2] ~ conc.mal, degree = 2, span = 0.75, family = 'gaussian')

  l2.mal.df = data.frame(mal = conc.mal,
                     pred = predict(l2.mal, newdata=data.frame(id = conc.mal)),
                     se = predict(l2.mal, newdata=data.frame(id = conc.mal), se = TRUE)$se.fit,
                     pred.2 = predict(l2.mal.2, newdata=data.frame(conc.mal = conc.mal)),
                     se.2 = predict(l2.mal.2, newdata=data.frame(conc.mal = conc.mal), se = TRUE)$se.fit)

#*********************************************************************

d3.mal = as.data.frame(r0.fill.mal[ , , 3])
  d3.mal2 = reshape(d3.mal, varying = c(colnames(d3.mal)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d3.mal)), direction = 'long')
#Fit loess to output of all simulations
  l3.mal = loess(r0 ~ id, data = d3.mal2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l3.mal.2 = loess(r0.fill.mal[ , 1, 3] ~ conc.mal, degree = 2, span = 0.75, family = 'gaussian')

  l3.mal.df = data.frame(mal = conc.mal,
                     pred = predict(l3.mal, newdata=data.frame(id = conc.mal)),
                     se = predict(l3.mal, newdata=data.frame(id = conc.mal), se = TRUE)$se.fit,
                     pred.2 = predict(l3.mal.2, newdata=data.frame(conc.mal = conc.mal)),
                     se.2 = predict(l3.mal.2, newdata=data.frame(conc.mal = conc.mal), se = TRUE)$se.fit)

#*********************************************************************

d4.mal = as.data.frame(r0.fill.mal[ , , 4])
  d4.mal2 = reshape(d4.mal, varying = c(colnames(d4.mal)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d4.mal)), direction = 'long')

#Fit loess to output of all simulations
  l4.mal = loess(r0 ~ id, data = d4.mal2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l4.mal.2 = loess(r0.fill.mal[ , 1, 4] ~ conc.mal, degree = 2, span = 0.75, family = 'gaussian')

  l4.mal.df = data.frame(mal = conc.mal,
                     pred = predict(l4.mal, newdata=data.frame(id = conc.mal)),
                     se = predict(l4.mal, newdata=data.frame(id = conc.mal), se = TRUE)$se.fit,
                     pred.2 = predict(l4.mal.2, newdata=data.frame(conc.mal = conc.mal)),
                     se.2 = predict(l4.mal.2, newdata=data.frame(conc.mal = conc.mal), se = TRUE)$se.fit)

#*********************************************************************  

d5.mal = as.data.frame(r0.fill.mal[ , , 5])
  d5.mal2 = reshape(d5.mal, varying = c(colnames(d5.mal)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d5.mal)), direction = 'long')

#Fit loess to output of all simulations
  l5.mal = loess(r0 ~ id, data = d5.mal2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l5.mal.2 = loess(r0.fill.mal[ , 1, 5] ~ conc.mal, degree = 2, span = 0.75, family = 'gaussian')

  l5.mal.df = data.frame(mal = conc.mal,
                     pred = predict(l5.mal, newdata=data.frame(id = conc.mal)),
                     se = predict(l5.mal, newdata=data.frame(id = conc.mal), se = TRUE)$se.fit,
                     pred.2 = predict(l5.mal.2, newdata=data.frame(conc.mal = conc.mal)),
                     se.2 = predict(l5.mal.2, newdata=data.frame(conc.mal = conc.mal), se = TRUE)$se.fit)

#*********************************************************************

d6.mal = as.data.frame(r0.fill.mal[ , , 6])
  d6.mal2 = reshape(d6.mal, varying = c(colnames(d6.mal)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d6.mal)), direction = 'long')

#Fit loess to output of all simulations
  l6.mal = loess(r0 ~ id, data = d6.mal2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l6.mal.2 = loess(r0.fill.mal[ , 1, 6] ~ conc.mal, degree = 2, span = 0.75, family = 'gaussian')

  l6.mal.df = data.frame(mal = conc.mal,
                     pred = predict(l6.mal, newdata=data.frame(id = conc.mal)),
                     se = predict(l6.mal, newdata=data.frame(id = conc.mal), se = TRUE)$se.fit,
                     pred.2 = predict(l6.mal.2, newdata=data.frame(conc.mal = conc.mal)),
                     se.2 = predict(l6.mal.2, newdata=data.frame(conc.mal = conc.mal), se = TRUE)$se.fit)
  
#*********************************************************************
d7.mal = as.data.frame(r0.fill.mal[ , , 7])
  d7.mal2 = reshape(d7.mal, varying = c(colnames(d7.mal)),
                    v.names = 'r0', timevar = 'sim', 
                    times = c(colnames(d7.mal)), direction = 'long')
  
#Fit loess to output of all simulations
  l7.mal = loess(r0 ~ id, data = d7.mal2, degree = 2, span = 0.75, family = 'gaussian') 
#Fit loess to output of 1 simulation
  l7.mal.2 = loess(r0.fill.mal[ , 1, 7] ~ conc.mal, degree = 2, span = 0.75, family = 'gaussian')
  
  l7.mal.df = data.frame(mal = conc.mal,
                     pred = predict(l7.mal, newdata=data.frame(id = conc.mal)),
                     se = predict(l7.mal, newdata=data.frame(id = conc.mal), se = TRUE)$se.fit,
                     pred.2 = predict(l7.mal.2, newdata=data.frame(conc.mal = conc.mal)),
                     se.2 = predict(l7.mal.2, newdata=data.frame(conc.mal = conc.mal), se = TRUE)$se.fit)
  

#Plot results ###################
#mean of all sims as points  
plot(conc.mal, conc.mal.means.r0[,1], pch = 17, cex=0.5, col = 2, xlab = 'Malathion (ppb)', ylab = expression(paste(R[0], '(q)')),
     ylim = c(0,7), xlim = c(0, max(conc.mal)))
  for(i in 2:(length(parfx)+2)){
    points(conc.mal, conc.mal.means.r0[,i], pch = 17, cex=0.5, col = i+1)
  }
  legend('topright', pch = 17, col = c(2,3,6,7,8), cex = 0.7,
         legend = c('piC', 'piM', 'muP', 'muN', 'piC & piM & muP'))
  
plot(conc.mal, conc.mal.means.r0[,1], pch = 17, cex=0.5, col = 2, xlab = 'Malathion (ppb)', ylab = expression(paste(R[0], '(q)')),
     ylim = c(0,4), xlim = c(0, 1000))
  for(i in 2:(length(parfx)+2)){
    points(conc.mal, conc.mal.means.r0[,i], pch = 17, cex=0.5, col = i+1)
  }  
  legend('bottomleft', pch = 17, col = c(4,5,7), cex = 0.7,
         legend = c('fN', 'muN', 'fN & muN'))

#loess regression lines of single sims
plot(c(0,max(conc.mal)), c(rep(r0.In(0)[3],2) - r0.In(0)[3]), type='l', lwd = 2, ylim = c(-4,3), xlim = c(0, max(conc.mal)),
     xlab = 'Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  lines(l1.mal.df$mal, l1.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=2)
    lines(l1.mal.df$mal, l1.mal.df$pred.2 + l1.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=2)
    lines(l1.mal.df$mal, l1.mal.df$pred.2 - l1.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=2)
  lines(l2.mal.df$mal, l2.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=5)
    lines(l2.mal.df$mal, l2.mal.df$pred.2 + l2.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=5)
    lines(l2.mal.df$mal, l2.mal.df$pred.2 - l2.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=5)
  lines(l3.mal.df$mal, l3.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=7)
    lines(l3.mal.df$mal, l3.mal.df$pred.2 + l3.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=7)
    lines(l3.mal.df$mal, l3.mal.df$pred.2 - l3.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=7)
  lines(l4.mal.df$mal, l4.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=4)
    lines(l4.mal.df$mal, l4.mal.df$pred.2 + l4.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=4)
    lines(l4.mal.df$mal, l4.mal.df$pred.2 - l4.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=4)
  lines(l5.mal.df$mal, l5.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=6)
    lines(l5.mal.df$mal, l5.mal.df$pred.2 + l5.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=6)
    lines(l5.mal.df$mal, l5.mal.df$pred.2 - l5.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=6)
  lines(l6.mal.df$mal, l6.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=8)
    lines(l6.mal.df$mal, l6.mal.df$pred.2 + l6.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=8)
    lines(l6.mal.df$mal, l6.mal.df$pred.2 - l6.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=8)
  lines(l7.mal.df$mal, l7.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=3)
    lines(l7.mal.df$mal, l7.mal.df$pred.2 + l7.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=3)
    lines(l7.mal.df$mal, l7.mal.df$pred.2 - l7.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=3)
    
      #Point estimates from studies
  points(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'phi_N'], r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'phi_N'] - r0.In(0)[3],
         pch = 17, cex = 0.8, col = 6)  
  points(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'piC'], r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'piC'] - r0.In(0)[3],
         pch = 17, cex = 0.8, col = 2)  
  points(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'muN'], r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'muN'] - r0.In(0)[3],
         pch = 17, cex = 0.8, col = 4)  
  legend('bottomleft', title = 'Parameters', legend = c('Baseline',
                                                        expression(paste(pi[C], ' meta')),
                                                        expression(paste(mu[N], ' Bakry 2012')),
                                                        expression(paste(Phi[N], ' Baxter 2011')),
                                                        'Combined functions'),
         lty = c(1, rep(2,4)), col = c(1,2,4,6,3), cex = 0.8)  
  title(expression(paste('Malathion influence on ', R[0], sep ='')))
  

  legend('topright', legend = c(expression(paste(pi[C], '(q)')),
                                expression(paste(pi[M], '(q)')),
                                expression(paste(f[N], '(q)')),
                                expression(paste(mu[N], '(q)')),
                                expression(paste(mu[P], '(q)'))), pch = 17, col = c(2:6), cex = 0.7)