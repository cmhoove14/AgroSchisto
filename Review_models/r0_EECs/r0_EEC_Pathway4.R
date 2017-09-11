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
source('Response_Fxs/griggs08_piC_beq.R')
source('Response_Fxs/koprivnikar06_piC_beq.R')
source('Response_Fxs/rohr08_piC_atrOnly.R')
source('Response_Fxs/piC_atr_meta_beq.R')
source('Response_Fxs/tantawy2002_piC.R')
source('Response_Fxs/Ghaffar2016_cercariae.R')
source('Response_Fxs/tchounwou92_piC2_beq.R')
source('Response_Fxs/Hasheesh2011_larval.R')
source('Response_Fxs/Ghaffar2016_miracidia.R')
source('Response_Fxs/tchounwou91_piM_beq.R')
source('Response_Fxs/tantawy2002_piM.R')
source('Response_Fxs/tchounwou91b_Fe_piM_beq.R')
source('Response_Fxs/tchounwou91b_Fe_egg_viability.R')

source('Review_models/r0_of_q.R')

library(parallel)
library(fBasics)

keep.fin.p4 = c(keep.grg08, keep.kop06.beq, keep.atr.rohr08, keep.meta.piC, keep.tantawy.piC, keep.gaf.piC,
                keep.tch92.beq, keep.hsh.pic, keep.hsh.pim, keep.gaf.piM, keep.tch91.beq, keep.tantawy.piM,
                keep.tch91.Fe, keep.tch91.egv, 'r0.In', 'r0.He', 'r0.Fe', 'r0.fix', 'parameters', 
                'nil0', 'nil1', 'keep.fin.p4')

rm(list = setdiff(ls(), keep.fin.p4))
dev.off()

no.cores = detectCores() - 1

#Run simulations of malathion concentrations, start with individual functions then combine ################
nsims = 1000         #Number of simulations to run
#Expected environmental concentrations for agrochemicals (NEED TO BE FINALIZED)
  eec.atr = rep(102, nsims)     #atrazine
  eec.but = rep(100, nsims)     #butachlor
  eec.btr = rep(100, nsims)     #butralin
  eec.fpb = rep(100, nsims)     #fluazifop-p-butyl
  eec.gly = rep(3700, nsims)    #Glyphosate
  eec.pen = rep(100, nsims)     #pendimethalin
  eec.mal = rep(583, nsims)     #malathion
  eec.ch = rep(64, nsims)       #chlorpyrifos
  eec.prof = rep(100, nsims)    #profenofos
  eec.urea = rep(100, nsims)    #urea
  eec.amm = rep(100, nsims)     #Ammonium sulphate
  
#All pathway 2 response functions  
parfx = c(piC.grg08_atr_unc2, piC.kop_atr_unc, piC.atr.rohr08.lin, piC.meta_atr_unc, piC.tant02_but.lin_unc,
          piC.tant02_fpb.exp_unc, piC.ghaf_butr.exp_unc, piC.ghaf_gly.exp_unc, piC.ghaf_pen.exp_unc, 
          piC.tch92_mal_unc, piC_ch_Hash11_uncertainty, piC_pr_Hash11_uncertainty, piM_ch_Hash11_uncertainty,
          piM_pr_Hash11_uncertainty, piM.ghaf_butr.exp_unc, piM.ghaf_gly.exp_unc, piM.ghaf_pen.exp_unc, piM.tch91_mal_unc, 
          piM.tant02_but.exp_unc, piM.tant02_fpb.lin_unc, piM.tch91_amm_unc, piM.tch91_ure_unc, tch91.egv.amm_unc, 
          tch91.egv.ure_unc)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.p4, 'uniroot.all', 'rdrm', 'LL.2', 'L.4'))

r0.eec.p4 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.eec.p4 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.0.5eec.p4 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.0.5eec.p4 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.0.1eec.p4 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.0.1eec.p4 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

set.seed(0)

#Fill r0 estimates for EEC #####################
  #individual parameters
    r0.eec.p4[, 1] = parSapply(clus1, eec.atr, r0.He, f.pi_Cq = piC.grg08_atr_unc2)[3,]
    r0.eec.p4[, 2] = parSapply(clus1, eec.atr, r0.He, f.pi_Cq = piC.kop_atr_unc)[3,]
    r0.eec.p4[, 3] = parSapply(clus1, eec.atr, r0.He, f.pi_Cq = piC.atr.rohr08.lin)[3,]
    r0.eec.p4[, 4] = parSapply(clus1, eec.atr, r0.He, f.pi_Cq = piC.meta_atr_unc)[3,]
    r0.eec.p4[, 5] = parSapply(clus1, eec.but, r0.He, f.pi_Cq = piC.tant02_but.lin_unc)[3,]
    r0.eec.p4[, 6] = parSapply(clus1, eec.fpb, r0.He, f.pi_Cq = piC.tant02_fpb.exp_unc)[3,]
    r0.eec.p4[, 7] = parSapply(clus1, eec.btr, r0.He, f.pi_Cq = piC.ghaf_butr.exp_unc)[3,]
    r0.eec.p4[, 8] = parSapply(clus1, eec.gly, r0.He, f.pi_Cq = piC.ghaf_gly.exp_unc)[3,]
    r0.eec.p4[, 9] = parSapply(clus1, eec.pen, r0.He, f.pi_Cq = piC.ghaf_pen.exp_unc)[3,]
    r0.eec.p4[, 10] = parSapply(clus1, eec.mal, r0.In, f.pi_Cq = piC.tch92_mal_unc)[3,]
    r0.eec.p4[, 11] = parSapply(clus1, eec.ch, r0.In, f.pi_Cq = piC_ch_Hash11_uncertainty)[3,]
    r0.eec.p4[, 12] = parSapply(clus1, eec.prof, r0.In, f.pi_Cq = piC_pr_Hash11_uncertainty)[3,]
    
    r0.eec.p4[, 13] = parSapply(clus1, eec.ch, r0.In, f.pi_Mq = piM_ch_Hash11_uncertainty)[3,]
    r0.eec.p4[, 14] = parSapply(clus1, eec.prof, r0.In, f.pi_Mq = piM_pr_Hash11_uncertainty)[3,]
    r0.eec.p4[, 15] = parSapply(clus1, eec.btr, r0.He, f.pi_Mq = piM.ghaf_butr.exp_unc)[3,]
    r0.eec.p4[, 16] = parSapply(clus1, eec.gly, r0.He, f.pi_Mq = piM.ghaf_gly.exp_unc)[3,]
    r0.eec.p4[, 17] = parSapply(clus1, eec.pen, r0.He, f.pi_Mq = piM.ghaf_pen.exp_unc)[3,]
    r0.eec.p4[, 18] = parSapply(clus1, eec.mal, r0.In, f.pi_Mq = piM.tch91_mal_unc)[3,]
    r0.eec.p4[, 19] = parSapply(clus1, eec.but, r0.He, f.pi_Mq = piM.tant02_but.exp_unc)[3,]
    r0.eec.p4[, 20] = parSapply(clus1, eec.fpb, r0.He, f.pi_Mq = piM.tant02_fpb.lin_unc)[3,]
    r0.eec.p4[, 21] = parSapply(clus1, eec.amm, r0.Fe, f.pi_Mq = piM.tch91_amm_unc)[3,]
    r0.eec.p4[, 22] = parSapply(clus1, eec.urea, r0.Fe, f.pi_Mq = piM.tch91_ure_unc)[3,]
    
    r0.eec.p4[, 23] = parSapply(clus1, eec.amm, r0.Fe, f.v_q = tch91.egv.amm_unc)[3,]
    r0.eec.p4[, 24] = parSapply(clus1, eec.urea, r0.Fe, f.v_q = tch91.egv.ure_unc)[3,]
    
    
  #fill parameter values for EEC values
    par.eec.p4[, 1] = parSapply(clus1, eec.atr, piC.grg08_atr_unc2) 
    par.eec.p4[, 2] = parSapply(clus1, eec.atr, piC.kop_atr_unc) 
    par.eec.p4[, 3] = parSapply(clus1, eec.atr, piC.atr.rohr08.lin) 
    par.eec.p4[, 4] = parSapply(clus1, eec.atr, piC.meta_atr_unc) 
    par.eec.p4[, 5] = parSapply(clus1, eec.but, piC.tant02_but.lin_unc) 
    par.eec.p4[, 6] = parSapply(clus1, eec.fpb, piC.tant02_fpb.exp_unc) 
    par.eec.p4[, 7] = parSapply(clus1, eec.btr, piC.ghaf_butr.exp_unc) 
    par.eec.p4[, 8] = parSapply(clus1, eec.gly, piC.ghaf_gly.exp_unc) 
    par.eec.p4[, 9] = parSapply(clus1, eec.pen, piC.ghaf_pen.exp_unc) 
    par.eec.p4[, 10] = parSapply(clus1, eec.mal, piC.tch92_mal_unc) 
    par.eec.p4[, 11] = parSapply(clus1, eec.ch, piC_ch_Hash11_uncertainty) 
    par.eec.p4[, 12] = parSapply(clus1, eec.prof, piC_pr_Hash11_uncertainty) 
    
    par.eec.p4[, 13] = parSapply(clus1, eec.ch, piM_ch_Hash11_uncertainty) 
    par.eec.p4[, 14] = parSapply(clus1, eec.prof, piM_pr_Hash11_uncertainty) 
    par.eec.p4[, 15] = parSapply(clus1, eec.btr, piM.ghaf_butr.exp_unc) 
    par.eec.p4[, 16] = parSapply(clus1, eec.gly, piM.ghaf_gly.exp_unc) 
    par.eec.p4[, 17] = parSapply(clus1, eec.pen, piM.ghaf_pen.exp_unc) 
    par.eec.p4[, 18] = parSapply(clus1, eec.mal, piM.tch91_mal_unc) 
    par.eec.p4[, 19] = parSapply(clus1, eec.but, piM.tant02_but.exp_unc) 
    par.eec.p4[, 20] = parSapply(clus1, eec.fpb, piM.tant02_fpb.lin_unc) 
    par.eec.p4[, 21] = parSapply(clus1, eec.amm, piM.tch91_amm_unc) 
    par.eec.p4[, 22] = parSapply(clus1, eec.urea, piM.tch91_ure_unc) 
    
    par.eec.p4[, 23] = parSapply(clus1, eec.amm, tch91.egv.amm_unc) 
    par.eec.p4[, 24] = parSapply(clus1, eec.urea, tch91.egv.ure_unc) 
    
#Fill r0 estimates for 50% EEC ########################
  #individual parameters
    r0.0.5eec.p4[, 1] = parSapply(clus1, 0.5*eec.atr, r0.He, f.pi_Cq = piC.grg08_atr_unc2)[3,]
    r0.0.5eec.p4[, 2] = parSapply(clus1, 0.5*eec.atr, r0.He, f.pi_Cq = piC.kop_atr_unc)[3,]
    r0.0.5eec.p4[, 3] = parSapply(clus1, 0.5*eec.atr, r0.He, f.pi_Cq = piC.atr.rohr08.lin)[3,]
    r0.0.5eec.p4[, 4] = parSapply(clus1, 0.5*eec.atr, r0.He, f.pi_Cq = piC.meta_atr_unc)[3,]
    r0.0.5eec.p4[, 5] = parSapply(clus1, 0.5*eec.but, r0.He, f.pi_Cq = piC.tant02_but.lin_unc)[3,]
    r0.0.5eec.p4[, 6] = parSapply(clus1, 0.5*eec.fpb, r0.He, f.pi_Cq = piC.tant02_fpb.exp_unc)[3,]
    r0.0.5eec.p4[, 7] = parSapply(clus1, 0.5*eec.btr, r0.He, f.pi_Cq = piC.ghaf_butr.exp_unc)[3,]
    r0.0.5eec.p4[, 8] = parSapply(clus1, 0.5*eec.gly, r0.He, f.pi_Cq = piC.ghaf_gly.exp_unc)[3,]
    r0.0.5eec.p4[, 9] = parSapply(clus1, 0.5*eec.pen, r0.He, f.pi_Cq = piC.ghaf_pen.exp_unc)[3,]
    r0.0.5eec.p4[, 10] = parSapply(clus1, 0.5*eec.mal, r0.In, f.pi_Cq = piC.tch92_mal_unc)[3,]
    r0.0.5eec.p4[, 11] = parSapply(clus1, 0.5*eec.ch, r0.In, f.pi_Cq = piC_ch_Hash11_uncertainty)[3,]
    r0.0.5eec.p4[, 12] = parSapply(clus1, 0.5*eec.prof, r0.In, f.pi_Cq = piC_pr_Hash11_uncertainty)[3,]
    
    r0.0.5eec.p4[, 13] = parSapply(clus1, 0.5*eec.ch, r0.In, f.pi_Mq = piM_ch_Hash11_uncertainty)[3,]
    r0.0.5eec.p4[, 14] = parSapply(clus1, 0.5*eec.prof, r0.In, f.pi_Mq = piM_pr_Hash11_uncertainty)[3,]
    r0.0.5eec.p4[, 15] = parSapply(clus1, 0.5*eec.btr, r0.He, f.pi_Mq = piM.ghaf_butr.exp_unc)[3,]
    r0.0.5eec.p4[, 16] = parSapply(clus1, 0.5*eec.gly, r0.He, f.pi_Mq = piM.ghaf_gly.exp_unc)[3,]
    r0.0.5eec.p4[, 17] = parSapply(clus1, 0.5*eec.pen, r0.He, f.pi_Mq = piM.ghaf_pen.exp_unc)[3,]
    r0.0.5eec.p4[, 18] = parSapply(clus1, 0.5*eec.mal, r0.In, f.pi_Mq = piM.tch91_mal_unc)[3,]
    r0.0.5eec.p4[, 19] = parSapply(clus1, 0.5*eec.but, r0.He, f.pi_Mq = piM.tant02_but.exp_unc)[3,]
    r0.0.5eec.p4[, 20] = parSapply(clus1, 0.5*eec.fpb, r0.He, f.pi_Mq = piM.tant02_fpb.lin_unc)[3,]
    r0.0.5eec.p4[, 21] = parSapply(clus1, 0.5*eec.amm, r0.Fe, f.pi_Mq = piM.tch91_amm_unc)[3,]
    r0.0.5eec.p4[, 22] = parSapply(clus1, 0.5*eec.urea, r0.Fe, f.pi_Mq = piM.tch91_ure_unc)[3,]
    
    r0.0.5eec.p4[, 23] = parSapply(clus1, 0.5*eec.amm, r0.Fe, f.v_q = tch91.egv.amm_unc)[3,]
    r0.0.5eec.p4[, 24] = parSapply(clus1, 0.5*eec.urea, r0.Fe, f.v_q = tch91.egv.ure_unc)[3,]
    
    
  #fill parameter values for EEC values
    par.0.5eec.p4[, 1] = parSapply(clus1, 0.5*eec.atr, piC.grg08_atr_unc2) 
    par.0.5eec.p4[, 2] = parSapply(clus1, 0.5*eec.atr, piC.kop_atr_unc) 
    par.0.5eec.p4[, 3] = parSapply(clus1, 0.5*eec.atr, piC.atr.rohr08.lin) 
    par.0.5eec.p4[, 4] = parSapply(clus1, 0.5*eec.atr, piC.meta_atr_unc) 
    par.0.5eec.p4[, 5] = parSapply(clus1, 0.5*eec.but, piC.tant02_but.lin_unc) 
    par.0.5eec.p4[, 6] = parSapply(clus1, 0.5*eec.fpb, piC.tant02_fpb.exp_unc) 
    par.0.5eec.p4[, 7] = parSapply(clus1, 0.5*eec.btr, piC.ghaf_butr.exp_unc) 
    par.0.5eec.p4[, 8] = parSapply(clus1, 0.5*eec.gly, piC.ghaf_gly.exp_unc) 
    par.0.5eec.p4[, 9] = parSapply(clus1, 0.5*eec.pen, piC.ghaf_pen.exp_unc) 
    par.0.5eec.p4[, 10] = parSapply(clus1, 0.5*eec.mal, piC.tch92_mal_unc) 
    par.0.5eec.p4[, 11] = parSapply(clus1, 0.5*eec.ch, piC_ch_Hash11_uncertainty) 
    par.0.5eec.p4[, 12] = parSapply(clus1, 0.5*eec.prof, piC_pr_Hash11_uncertainty) 
    
    par.0.5eec.p4[, 13] = parSapply(clus1, 0.5*eec.ch, piM_ch_Hash11_uncertainty) 
    par.0.5eec.p4[, 14] = parSapply(clus1, 0.5*eec.prof, piM_pr_Hash11_uncertainty) 
    par.0.5eec.p4[, 15] = parSapply(clus1, 0.5*eec.btr, piM.ghaf_butr.exp_unc) 
    par.0.5eec.p4[, 16] = parSapply(clus1, 0.5*eec.gly, piM.ghaf_gly.exp_unc) 
    par.0.5eec.p4[, 17] = parSapply(clus1, 0.5*eec.pen, piM.ghaf_pen.exp_unc) 
    par.0.5eec.p4[, 18] = parSapply(clus1, 0.5*eec.mal, piM.tch91_mal_unc) 
    par.0.5eec.p4[, 19] = parSapply(clus1, 0.5*eec.but, piM.tant02_but.exp_unc) 
    par.0.5eec.p4[, 20] = parSapply(clus1, 0.5*eec.fpb, piM.tant02_fpb.lin_unc) 
    par.0.5eec.p4[, 21] = parSapply(clus1, 0.5*eec.amm, piM.tch91_amm_unc) 
    par.0.5eec.p4[, 22] = parSapply(clus1, 0.5*eec.urea, piM.tch91_ure_unc) 
    
    par.0.5eec.p4[, 23] = parSapply(clus1, 0.5*eec.amm, tch91.egv.amm_unc) 
    par.0.5eec.p4[, 24] = parSapply(clus1, 0.5*eec.urea, tch91.egv.ure_unc) 
    
#Fill r0 estimates for 10% EEC #####################
  #individual parameters
    r0.0.1eec.p4[, 1] = parSapply(clus1, 0.1*eec.atr, r0.He, f.pi_Cq = piC.grg08_atr_unc2)[3,]
    r0.0.1eec.p4[, 2] = parSapply(clus1, 0.1*eec.atr, r0.He, f.pi_Cq = piC.kop_atr_unc)[3,]
    r0.0.1eec.p4[, 3] = parSapply(clus1, 0.1*eec.atr, r0.He, f.pi_Cq = piC.atr.rohr08.lin)[3,]
    r0.0.1eec.p4[, 4] = parSapply(clus1, 0.1*eec.atr, r0.He, f.pi_Cq = piC.meta_atr_unc)[3,]
    r0.0.1eec.p4[, 5] = parSapply(clus1, 0.1*eec.but, r0.He, f.pi_Cq = piC.tant02_but.lin_unc)[3,]
    r0.0.1eec.p4[, 6] = parSapply(clus1, 0.1*eec.fpb, r0.He, f.pi_Cq = piC.tant02_fpb.exp_unc)[3,]
    r0.0.1eec.p4[, 7] = parSapply(clus1, 0.1*eec.btr, r0.He, f.pi_Cq = piC.ghaf_butr.exp_unc)[3,]
    r0.0.1eec.p4[, 8] = parSapply(clus1, 0.1*eec.gly, r0.He, f.pi_Cq = piC.ghaf_gly.exp_unc)[3,]
    r0.0.1eec.p4[, 9] = parSapply(clus1, 0.1*eec.pen, r0.He, f.pi_Cq = piC.ghaf_pen.exp_unc)[3,]
    r0.0.1eec.p4[, 10] = parSapply(clus1, 0.1*eec.mal, r0.In, f.pi_Cq = piC.tch92_mal_unc)[3,]
    r0.0.1eec.p4[, 11] = parSapply(clus1, 0.1*eec.ch, r0.In, f.pi_Cq = piC_ch_Hash11_uncertainty)[3,]
    r0.0.1eec.p4[, 12] = parSapply(clus1, 0.1*eec.prof, r0.In, f.pi_Cq = piC_pr_Hash11_uncertainty)[3,]
    
    r0.0.1eec.p4[, 13] = parSapply(clus1, 0.1*eec.ch, r0.In, f.pi_Mq = piM_ch_Hash11_uncertainty)[3,]
    r0.0.1eec.p4[, 14] = parSapply(clus1, 0.1*eec.prof, r0.In, f.pi_Mq = piM_pr_Hash11_uncertainty)[3,]
    r0.0.1eec.p4[, 15] = parSapply(clus1, 0.1*eec.btr, r0.He, f.pi_Mq = piM.ghaf_butr.exp_unc)[3,]
    r0.0.1eec.p4[, 16] = parSapply(clus1, 0.1*eec.gly, r0.He, f.pi_Mq = piM.ghaf_gly.exp_unc)[3,]
    r0.0.1eec.p4[, 17] = parSapply(clus1, 0.1*eec.pen, r0.He, f.pi_Mq = piM.ghaf_pen.exp_unc)[3,]
    r0.0.1eec.p4[, 18] = parSapply(clus1, 0.1*eec.mal, r0.In, f.pi_Mq = piM.tch91_mal_unc)[3,]
    r0.0.1eec.p4[, 19] = parSapply(clus1, 0.1*eec.but, r0.He, f.pi_Mq = piM.tant02_but.exp_unc)[3,]
    r0.0.1eec.p4[, 20] = parSapply(clus1, 0.1*eec.fpb, r0.He, f.pi_Mq = piM.tant02_fpb.lin_unc)[3,]
    r0.0.1eec.p4[, 21] = parSapply(clus1, 0.1*eec.amm, r0.Fe, f.pi_Mq = piM.tch91_amm_unc)[3,]
    r0.0.1eec.p4[, 22] = parSapply(clus1, 0.1*eec.urea, r0.Fe, f.pi_Mq = piM.tch91_ure_unc)[3,]
    
    r0.0.1eec.p4[, 23] = parSapply(clus1, 0.1*eec.amm, r0.Fe, f.v_q = tch91.egv.amm_unc)[3,]
    r0.0.1eec.p4[, 24] = parSapply(clus1, 0.1*eec.urea, r0.Fe, f.v_q = tch91.egv.ure_unc)[3,]
    
    
  #fill parameter values for EEC values
    par.0.1eec.p4[, 1] = parSapply(clus1, 0.1*eec.atr, piC.grg08_atr_unc2) 
    par.0.1eec.p4[, 2] = parSapply(clus1, 0.1*eec.atr, piC.kop_atr_unc) 
    par.0.1eec.p4[, 3] = parSapply(clus1, 0.1*eec.atr, piC.atr.rohr08.lin) 
    par.0.1eec.p4[, 4] = parSapply(clus1, 0.1*eec.atr, piC.meta_atr_unc) 
    par.0.1eec.p4[, 5] = parSapply(clus1, 0.1*eec.but, piC.tant02_but.lin_unc) 
    par.0.1eec.p4[, 6] = parSapply(clus1, 0.1*eec.fpb, piC.tant02_fpb.exp_unc) 
    par.0.1eec.p4[, 7] = parSapply(clus1, 0.1*eec.btr, piC.ghaf_butr.exp_unc) 
    par.0.1eec.p4[, 8] = parSapply(clus1, 0.1*eec.gly, piC.ghaf_gly.exp_unc) 
    par.0.1eec.p4[, 9] = parSapply(clus1, 0.1*eec.pen, piC.ghaf_pen.exp_unc) 
    par.0.1eec.p4[, 10] = parSapply(clus1, 0.1*eec.mal, piC.tch92_mal_unc) 
    par.0.1eec.p4[, 11] = parSapply(clus1, 0.1*eec.ch, piC_ch_Hash11_uncertainty) 
    par.0.1eec.p4[, 12] = parSapply(clus1, 0.1*eec.prof, piC_pr_Hash11_uncertainty) 
    
    par.0.1eec.p4[, 13] = parSapply(clus1, 0.1*eec.ch, piM_ch_Hash11_uncertainty) 
    par.0.1eec.p4[, 14] = parSapply(clus1, 0.1*eec.prof, piM_pr_Hash11_uncertainty) 
    par.0.1eec.p4[, 15] = parSapply(clus1, 0.1*eec.btr, piM.ghaf_butr.exp_unc) 
    par.0.1eec.p4[, 16] = parSapply(clus1, 0.1*eec.gly, piM.ghaf_gly.exp_unc) 
    par.0.1eec.p4[, 17] = parSapply(clus1, 0.1*eec.pen, piM.ghaf_pen.exp_unc) 
    par.0.1eec.p4[, 18] = parSapply(clus1, 0.1*eec.mal, piM.tch91_mal_unc) 
    par.0.1eec.p4[, 19] = parSapply(clus1, 0.1*eec.but, piM.tant02_but.exp_unc) 
    par.0.1eec.p4[, 20] = parSapply(clus1, 0.1*eec.fpb, piM.tant02_fpb.lin_unc) 
    par.0.1eec.p4[, 21] = parSapply(clus1, 0.1*eec.amm, piM.tch91_amm_unc) 
    par.0.1eec.p4[, 22] = parSapply(clus1, 0.1*eec.urea, piM.tch91_ure_unc) 
    
    par.0.1eec.p4[, 23] = parSapply(clus1, 0.1*eec.amm, tch91.egv.amm_unc) 
    par.0.1eec.p4[, 24] = parSapply(clus1, 0.1*eec.urea, tch91.egv.ure_unc) 

stopCluster(clus1) 
#Post process ############ 
#EEC runs ################
eec.p4.df = data.frame(chem = c(rep('Atrazine',4), 'Butachlor', 'Fluazifop-p-butyl', 'Butralin', 'Glyphosate', 'Pendimethalin',
                                'Malathion', 'Chlorpyrifos', 'Profenofos', 'Chlorpyrifos', 'Profenofos', 'Butralin', 'Glyphosate', 
                                'Pendimethalin', 'Malathion', 'Butachlor', 'Fluazifop-p-butyl', 'Ammonium Sulphate', 'Urea',
                                'Ammonium Sulphate', 'Urea'),
                       study = c('Griggs et al 2008', 'Koprivnikar et al 2006', 'Rohr et al 2008', 'Meta',
                                 rep('Tantawy 2002', 2), rep('Abdel-Ghaffar et al 2016', 3), 'Tchounwou et al 1992',
                                 rep('Hasheesh & Mohamed 2011', 4), rep('Abdel-Ghaffar et al 2016', 3), 'Tchounwou et al 1991a',
                                 rep('Tantawy 2002', 2), rep('Tchounwou et al 1991b', 4)),
                       Species = c(rep('Echinistoma trivolvis', 4), rep('Schistosoma mansoni', 6), rep('Schistosoma haemotobium', 4),
                                   rep('Schistosoma mansoni', 10)),
                       Parameter = c(rep('pi_C', 12), rep('pi_M', 10), rep('v', 2)),
                       r0 = colMeans(r0.eec.p4),
                       r0.sd = apply(r0.eec.p4, 2, sd),
                       par.mean = colMeans(par.eec.p4))

eec.p4.df$r0.up = eec.p4.df$r0 + eec.p4.df$r0.sd
eec.p4.df$r0.lo = eec.p4.df$r0 - eec.p4.df$r0.sd

save(eec.p4.df, file = 'Review_models/r0_EECs/eec.p4.df.RData')

#50% EEC values #################
eec0.5.p4.df = data.frame(chem = c(rep('Atrazine',4), 'Butachlor', 'Fluazifop-p-butyl', 'Butralin', 'Glyphosate', 'Pendimethalin',
                                   'Malathion', 'Chlorpyrifos', 'Profenofos', 'Chlorpyrifos', 'Profenofos', 'Butralin', 'Glyphosate', 
                                   'Pendimethalin', 'Malathion', 'Butachlor', 'Fluazifop-p-butyl', 'Ammonium Sulphate', 'Urea',
                                   'Ammonium Sulphate', 'Urea'),
                          study = c('Griggs et al 2008', 'Koprivnikar et al 2006', 'Rohr et al 2008', 'Meta',
                                    rep('Tantawy 2002', 2), rep('Abdel-Ghaffar et al 2016', 3), 'Tchounwou et al 1992',
                                    rep('Hasheesh & Mohamed 2011', 4), rep('Abdel-Ghaffar et al 2016', 3), 'Tchounwou et al 1991a',
                                    rep('Tantawy 2002', 2), rep('Tchounwou et al 1991b', 4)),
                          Species = c(rep('Echinistoma trivolvis', 4), rep('Schistosoma mansoni', 6), rep('Schistosoma haemotobium', 4),
                                      rep('Schistosoma mansoni', 10)),
                          Parameter = c(rep('pi_C', 12), rep('pi_M', 10), rep('v', 2)),
                          r0 = colMeans(r0.0.5eec.p4),
                          r0.sd = apply(r0.0.5eec.p4, 2, sd),
                          par.mean = colMeans(par.0.5eec.p4))

eec0.5.p4.df$r0.up = eec0.5.p4.df$r0 + eec0.5.p4.df$r0.sd
eec0.5.p4.df$r0.lo = eec0.5.p4.df$r0 - eec0.5.p4.df$r0.sd

save(eec0.5.p4.df, file = 'Review_models/r0_EECs/eec0.5.p4.df.RData')


#10% EEC values ######################
eec0.1.p4.df = data.frame(chem = c(rep('Atrazine',4), 'Butachlor', 'Fluazifop-p-butyl', 'Butralin', 'Glyphosate', 'Pendimethalin',
                                   'Malathion', 'Chlorpyrifos', 'Profenofos', 'Chlorpyrifos', 'Profenofos', 'Butralin', 'Glyphosate', 
                                   'Pendimethalin', 'Malathion', 'Butachlor', 'Fluazifop-p-butyl', 'Ammonium Sulphate', 'Urea',
                                   'Ammonium Sulphate', 'Urea'),
                          study = c('Griggs et al 2008', 'Koprivnikar et al 2006', 'Rohr et al 2008', 'Meta',
                                    rep('Tantawy 2002', 2), rep('Abdel-Ghaffar et al 2016', 3), 'Tchounwou et al 1992',
                                    rep('Hasheesh & Mohamed 2011', 4), rep('Abdel-Ghaffar et al 2016', 3), 'Tchounwou et al 1991a',
                                    rep('Tantawy 2002', 2), rep('Tchounwou et al 1991b', 4)),
                          Species = c(rep('Echinistoma trivolvis', 4), rep('Schistosoma mansoni', 6), rep('Schistosoma haemotobium', 4),
                                      rep('Schistosoma mansoni', 10)),
                          Parameter = c(rep('pi_C', 12), rep('pi_M', 10), rep('v', 2)),
                          r0 = colMeans(r0.0.1eec.p4),
                          r0.sd = apply(r0.0.1eec.p4, 2, sd),
                          par.mean = colMeans(par.0.1eec.p4))

eec0.1.p4.df$r0.up = eec0.1.p4.df$r0 + eec0.1.p4.df$r0.sd
eec0.1.p4.df$r0.lo = eec0.1.p4.df$r0 - eec0.1.p4.df$r0.sd

save(eec0.1.p4.df, file = 'Review_models/r0_EECs/eec0.1.p4.df.RData')
