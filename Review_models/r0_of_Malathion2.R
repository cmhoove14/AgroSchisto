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
source('Response_Fxs/tchounwou92_piC2_b_e_q_function.R')
source('Response_Fxs/tchounwou92_piM_b_e_q.R')
source('Response_Fxs/Halstead_Insecticides2015.R')
source('Response_Fxs/bakry2011.R')

source('Review_models/r0_of_q.R')

library(parallel)

keep.fin = c(keep.tch92.beq, keep.tch91.beq, keep.hal15.muP[c(1,7)], keep.bak11.N[c(1,2,5,6)], 
             'r0.In', 'parameters', 'nil0', 'nil1', 'keep.fin')

rm(list = setdiff(ls(), keep.fin))
dev.off()

no.cores = detectCores() - 1

#Run 100 simulations of malathion concentrations from 0 - 10000, start with individual functions then combine ################
conc = c(c(0:999), seq(1000, 200000, length.out = 1001))  #Concentration range to test
nsims = 100         #Number of simulations to run
parfx = c(piC.tch92_mal_unc, piM.tch91_mal_unc, #Functions corresponding to affected parameters
          fNq_mal_Bakry11_uncertainty, muNq_mal_Bakry11_uncertainty, muPq_mal_Halstead_uncertainty)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin, 'uniroot.all', 'rdrm', 'LL.2'))

ar.fill = array(data = NA, dim = c(length(conc), nsims, length(parfx)+1))

set.seed(0)

  for(i in 1:nsims){
    ar.fill[, i, 1] = parSapply(clus1, conc, r0.In, f.pi_Cq = piC.tch92_mal_unc)[3,]
    ar.fill[, i, 2] = parSapply(clus1, conc, r0.In, f.pi_Mq = piM.tch91_mal_unc)[3,]
    ar.fill[, i, 3] = parSapply(clus1, conc, r0.In, f.f_Nq = fNq_mal_Bakry11_uncertainty)[3,]
    ar.fill[, i, 4] = parSapply(clus1, conc, r0.In, f.mu_Nq = muNq_mal_Bakry11_uncertainty)[3,]
    ar.fill[, i, 5] = parSapply(clus1, conc, r0.In, f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
    ar.fill[, i, 6] = parSapply(clus1, conc, r0.In,
                                              f.mu_Nq = muNq_mal_Bakry11_uncertainty,
                                              f.f_Nq = fNq_mal_Bakry11_uncertainty,
                                              f.pi_Mq = piM.tch91_mal_unc,
                                              f.pi_Cq = piC.tch92_mal_unc,
                                              f.mu_Pq = muPq_mal_Halstead_uncertainty)[3,]
  }

stopCluster(clus1)
 
conc.means = matrix(nrow = length(conc), ncol = length(parfx)+1)

  for(i in 1:length(parfx)+1){
    conc.means[,i] = rowMeans(ar.fill[ , , i])
  }

#Plot results of individual parameter responses as points
plot(conc, conc.means[,1], pch = 17, cex=0.5, col = 2, xlab = 'Malathion (ppb)', ylab = expression(paste(R[0], '(q)')),
     ylim = c(0,10), xlim = c(0, max(conc)))
  for(i in 2:length(parfx)){
    points(conc, conc.means[,i], pch = 17, cex=0.5, col = i+1)
  }

#Plot results as lines
plot(conc, conc.means[,1], pch = 17, cex=0.5, col = 2, ylim = c(0,10), xlim = c(0, max(conc)), type = 'l',
     xlab = 'Malathion (ppb)', ylab = expression(paste(R[0], '(q)')))
  for(i in 2:length(parfx)){
    lines(conc, conc.means[,i], cex=0.5, col = i+1)
  }

#Plot results of combined function as points
plot(conc, conc.means[,6], pch = 17, cex=0.5, col = 2, ylim = c(0,5), xlim = c(0, 1000),
     xlab = 'Malathion (ppb)', ylab = expression(paste(R[0], '(q) combined')))

#Plot results with smooth scatter; only good to view a single parameter, can't do multiple parameters at once
pal1 = colorRampPalette(c('white', 'red'))
smoothScatter(conc, conc.means[,1], xlab = 'Malathion (ppb)', ylab = expression(paste(R[0], '(q)')),
              colramp = c('red'), ylim = c(0,10), xlim = c(0, max(conc)), pch = NA)


  legend('topright', legend = c(expression(paste(pi[C], '(q)')),
                                expression(paste(pi[M], '(q)')),
                                expression(paste(f[N], '(q)')),
                                expression(paste(mu[N], '(q)')),
                                expression(paste(mu[P], '(q)'))), pch = 17, col = c(2:6), cex = 0.7)