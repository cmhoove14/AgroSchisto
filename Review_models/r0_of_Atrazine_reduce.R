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
source('Response_Fxs/fin/bakry2012_fin.R')
source('Response_Fxs/fin/piC_atr_meta_beq_fin.R')
source('Response_Fxs/fin/koprivnikar06_piC_beq_fin.R')
source('Response_Fxs/fin/griggs08_piC_beq_fin.R')
source('Response_Fxs/fin/Baxter_Rohr_Atrazine2011_fin.R')
source('Response_Fxs/fin/rohr08_piC_atrOnly_fin.R')
source('Response_Fxs/fin/Omran&Salama_snails_finlw49.R')

source('Review_models/fin/r0_of_q_fin.R')

library(parallel)
library(reshape2)

keep.fin.atr = c(keep.bak.atr, keep.grg08, keep.kop06.beq, keep.meta.piC, keep.baxrohr, keep.atr.rohr08, keep.ons.atr,
                 'r0.He', 'r0.fix', 'kmat', 'parameters', 'nil0', 'nil1', 'keep.fin.atr')

rm(list = setdiff(ls(), keep.fin.atr))
dev.off()

no.cores = detectCores() - 1

#Fill data frame with point estimates of R0 driven by data in each study ################
r0.atr.fix = data.frame(atr = c(atra.df$atra, 
                                kopatr.auc$atr,
                                grgc.df$atr,
                                atr$conc,
                                mun.bak.atr$conc,
                                ons.atr$conc),
                        par = c(rep('phi_N', length(atra.df$atra)), 
                                rep('piC', length(kopatr.auc$atr)),
                                rep('piC', length(grgc.df$atr)),
                                rep('piC', length(atr$conc)),
                                rep('muN', length(mun.bak.atr$conc)),
                                rep('muN', length(ons.atr$conc))),
                        study = c(rep('Baxter11', length(atra.df$atra)), 
                                  rep('Koprivnikar06', length(kopatr.auc$atr)),
                                  rep('Griggs08', length(grgc.df$atr)),
                                  rep('Rohr08', length(atr$conc)),
                                  rep('Bakry12', length(mun.bak.atr$conc)),
                                  rep('Omran&Salama', length(ons.atr$conc))),
                        r0 = c(r0.fix(phi_Nqx = atra.df$growthrate[1] / atra.df$growthrate[1])[3],
                               r0.fix(phi_Nqx = atra.df$growthrate[2] / atra.df$growthrate[1])[3],
                               r0.fix(phi_Nqx = atra.df$growthrate[3] / atra.df$growthrate[1])[3],
                               r0.fix(phi_Nqx = atra.df$growthrate[4] / atra.df$growthrate[1])[3],
                               r0.fix(phi_Nqx = atra.df$growthrate[5] / atra.df$growthrate[1])[3],
                               r0.fix(pi_Cqx = kopatr.auc$piC[1])[3],
                               r0.fix(pi_Cqx = kopatr.auc$piC[2])[3],
                               r0.fix(pi_Cqx = kopatr.auc$piC[3])[3],
                               r0.fix(pi_Cqx = grg08_atr_aucs[1]/grg08_atr_aucs[1])[3],
                               r0.fix(pi_Cqx = grg08_atr_aucs[2]/grg08_atr_aucs[1])[3],
                               r0.fix(pi_Cqx = grg08_atr_aucs[3]/grg08_atr_aucs[1])[3],
                               r0.fix(pi_Cqx = (atr$surv[1] / atr$surv[1]))[3],
                               r0.fix(pi_Cqx = (atr$surv[2] / atr$surv[1]))[3],
                               r0.fix(pi_Cqx = (atr$surv[3] / atr$surv[1]))[3],
                               r0.fix(pi_Cqx = (atr$surv[4] / atr$surv[1]))[3],
                               r0.fix(pi_Cqx = (atr$surv[5] / atr$surv[1]))[3],
                               r0.fix(pi_Cqx = (atr$surv[6] / atr$surv[1]))[3],
                               r0.fix(mu_Nqx = mun.bak.atr$mort[1])[3],
                               r0.fix(mu_Nqx = mun.bak.atr$mort[2])[3],
                               r0.fix(mu_Nqx = mun.bak.atr$mort[3])[3],
                               r0.fix(mu_Nqx = ons.atr$mort[1])[3],
                               r0.fix(mu_Nqx = ons.atr$mort[2])[3],
                               r0.fix(mu_Nqx = ons.atr$mort[3])[3],
                               r0.fix(mu_Nqx = ons.atr$mort[4])[3],
                               r0.fix(mu_Nqx = ons.atr$mort[5])[3],
                               r0.fix(mu_Nqx = ons.atr$mort[6])[3],
                               r0.fix(mu_Nqx = ons.atr$mort[7])[3]))

r0.atr.fix.n0 = subset(r0.atr.fix, atr != 0)  

#Run simulations of atrazine concentrations from 0 - 2000, start with individual functions then combine ################
conc.atr = seq(0, 204, length.out = 100)  #Concentration range to test (atrazine EEC = 102)
nsims = 100       #Number of simulations to run
parfx = c(piC.meta_atr_unc, fN.atr.bak.uncertainty, muNq_atr_Bakry12_uncertainty, phi_Nq_atr_baxrohr, #Functions corresponding to affected parameters
          piC.grg08_atr_unc, auc.kop.lin.atr0, piC.atr.rohr08.lin, ons.munq.atr)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.atr, 'uniroot.all', 'rdrm', 'LL.2', 'kmat'))

r0.atr.fill = array(data = NA, dim = c(length(conc.atr), nsims, length(parfx)+3))

clusterSetRNGStream(cl = clus1, iseed = 43093) #set cluster seed for reproducibility

  for(i in 1:nsims){
  #Fill r0 estimates ######
    #individual parameters
    r0.atr.fill[, i, 1] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.meta_atr_unc)[3,]
    r0.atr.fill[, i, 2] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.grg08_atr_unc)[3,]
    r0.atr.fill[, i, 3] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = auc.kop.lin.atr0)[3,]
    r0.atr.fill[, i, 4] = parSapply(clus1, conc.atr, r0.He, f.pi_Cq = piC.atr.rohr08.lin)[3,]
    r0.atr.fill[, i, 5] = parSapply(clus1, conc.atr, r0.He, f.f_Nq = fN.atr.bak.uncertainty)[3,]
    r0.atr.fill[, i, 6] = parSapply(clus1, conc.atr, r0.He, f.mu_Nq = muNq_atr_Bakry12_uncertainty)[3,]
    r0.atr.fill[, i, 7] = parSapply(clus1, conc.atr, r0.He, f.phi_Nq = phi_Nq_atr_baxrohr.no30)[3,]
    r0.atr.fill[, i, 8] = parSapply(clus1, conc.atr, r0.He, f.mu_Nq = ons.munq.atr)[3,]
  
    #combined estimate 1: Bakry et al snail data
    r0.atr.fill[, i, 9] = parSapply(clus1, conc.atr, r0.He,
                                    f.mu_Nq = muNq_atr_Bakry12_uncertainty,
                                    f.f_Nq = fN.atr.bak.uncertainty,
                                    f.pi_Cq = piC.atr.rohr08.lin,
                                    f.phi_Nq = phi_Nq_atr_baxrohr)[3,]
    #combined estimate 2: Omran and Salama snail data
    r0.atr.fill[, i, 10] = parSapply(clus1, conc.atr, r0.He,
                                    f.mu_Nq = ons.munq.atr,
                                    f.pi_Cq = piC.atr.rohr08.lin,
                                    f.phi_Nq = phi_Nq_atr_baxrohr)[3,]
    
    #combined estimate 3: Omran and Salama snail data, no cercariae effect
    r0.atr.fill[, i, 11] = parSapply(clus1, conc.atr, r0.He,
                                     f.mu_Nq = ons.munq.atr,
                                     f.phi_Nq = phi_Nq_atr_baxrohr)[3,]
  
    print(i)
  }

stopCluster(clus1)
 
########## Post process ############
#median and IQR of r0 simulations #########
#median r0
  #raw r0 value
    atr.med.r0 = matrix(nrow = length(conc.atr), ncol = length(parfx)+3)
      for(i in 1:(length(parfx)+2)){
        atr.med.r0[,i] = apply(r0.atr.fill[ , , i], 1, median)
      }
    
  #change in r0 as a percent from baseline
    atr.med.r0_percent = matrix(nrow = length(conc.atr), ncol = length(parfx)+3)
    for(i in 1:(length(parfx)+2)){
      atr.med.r0_percent[,i] = apply(r0.atr.fill[ , , i], 1, median) / r0.He(0)[3] *100-100
    }

  #25 %ile r0
    #raw r0 value
      atr.med.r0.25 = matrix(nrow = length(conc.atr), ncol = length(parfx)+3)
      for(i in 1:(length(parfx)+2)){
        atr.med.r0.25[,i] = apply(r0.atr.fill[ , , i], 1, quantile, probs = 0.25)
      }

    #change in r0 as a percent from baseline
      atr.med.r0.25_percent = matrix(nrow = length(conc.atr), ncol = length(parfx)+3)
      for(i in 1:(length(parfx)+2)){
        atr.med.r0.25_percent[,i] = apply(r0.atr.fill[ , , i], 1, quantile, probs = 0.25) / r0.He(0)[3] *100-100
      }

  #75 %ile r0
    #raw r0 value
      atr.med.r0.75 = matrix(nrow = length(conc.atr), ncol = length(parfx)+3)
      for(i in 1:(length(parfx)+2)){
        atr.med.r0.75[,i] = apply(r0.atr.fill[ , , i], 1, quantile, probs = 0.75)
      }
    #change in r0 as a percent from baseline
      atr.med.r0.75_percent = matrix(nrow = length(conc.atr), ncol = length(parfx)+3)
      for(i in 1:(length(parfx)+2)){
        atr.med.r0.75_percent[,i] = apply(r0.atr.fill[ , , i], 1, quantile, probs = 0.75) / r0.He(0)[3] *100-100
      }

plot(conc.atr, atr.med.r0[,10], type= 'l', lwd = 2, xlab = 'atrazine (ppb)', ylab = 'delta r0 (%)', ylim = c(1.5, 2.5))
      lines(conc.atr, atr.med.r0.25[,10], lty = 2)
      lines(conc.atr, atr.med.r0.75[,10], lty = 2)
      
plot(conc.atr, atr.med.r0_percent[,10], type= 'l', lwd = 2, xlab = 'atrazine (ppb)', ylab = 'delta r0 (%)',
     ylim = c(-50,50))
  lines(conc.atr, atr.med.r0.25_percent[,10], lty = 2)
  lines(conc.atr, atr.med.r0.75_percent[,10], lty = 2)