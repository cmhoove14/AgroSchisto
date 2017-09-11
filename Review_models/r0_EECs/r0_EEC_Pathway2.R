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
source('Response_Fxs/tchounwou91_fN-muN.R')
source('Response_Fxs/Ghaffar2016_snails.R')
source('Response_Fxs/bakry2011.R')
source('Response_Fxs/Hasheesh2011_snails.R')
source('Response_Fxs/Ibrahim92_fN_muN.R')
source('Response_Fxs/Omran&Salama_snails.R')
source('Response_Fxs/tantawy2002_muN.R')
source('Response_Fxs/ragab2006.R')


source('Review_models/r0_of_q.R')

library(parallel)
library(fBasics)

keep.fin.p2 = c(keep.bak.all, keep.tch91.snail, keep.gaf.all, keep.bak11.N,
                keep.hsh.all, keep.ibr.ch, keep.ons.all, keep.muN.tantawy, keep.ragab.mun,
                 'r0.In', 'r0.He', 'r0.Fe', 'r0.fix', 'parameters', 'nil0', 'nil1', 'keep.fin.p2')

rm(list = setdiff(ls(), keep.fin.p2))
dev.off()

no.cores = detectCores() - 1

#Run simulations of malathion concentrations, start with individual functions then combine ################
nsims = 1000         #Number of simulations to run
#Expected environmental concentrations for agrochemicals (NEED TO BE FINALIZED)
  eec.atr = rep(102, nsims)     #atrazine
  eec.gly = rep(3700, nsims)    #Glyphosate
  eec.ch = rep(64, nsims)       #chlorpyrifos
  eec.mal = rep(583, nsims)     #malathion
  eec.but = rep(100, nsims)     #butralin
  eec.pen = rep(100, nsims)     #pendimethalin
  eec.del = rep(100, nsims)     #deltamethrin
  eec.prof = rep(100, nsims)    #profenofos
  eec.fpb = rep(100, nsims)     #fluazifop-p-butyl
  eec.urea = rep(100, nsims)    #urea
  eec.pot = rep(100, nsims)     #potassium sulphate
  eec.amm = rep(100, nsims)     #Ammonium phosphate
  
#All pathway 2 response functions  
parfx = c(fN.atr.bak.uncertainty, fN.gly.bak.uncertainty, muNq_gly_Bakry12_uncertainty,
          muNq_atr_Bakry12_uncertainty, muNq_mal_tch91_uncertainty, fNq_mal_tch91_uncertainty,
          mu_Nq_butr_gaf16_uncertainty, mu_Nq_gly_gaf16_uncertainty, mu_Nq_pen_gaf16_uncertainty,
          fN.butr.fx.uncertainty, fN.gly.fx.uncertainty, fN.pen.fx.uncertainty,
          muNq_mal_Bakry11_uncertainty, muNq_del_Bakry11_uncertainty, fNq_mal_Bakry11_uncertainty, 
          fNq_del_Bakry11_uncertainty, fN.hash.chlor.uncertainty, fN.hash.prof.uncertainty,
          muNq_ch_hash11_uncertainty, muNq_prof_hash11_uncertainty, f_N_chlor_ibr92_uncertainty,
          mu_N_chlor_ibr92_uncertainty, ons.munq.atr, ons.munq.gly, muN.tant.but_uncertainty, 
          muN.tant.fpb_uncertainty, rag06_mun_urea, rag06_mun_pot, rag06_mun_amm)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.fin.p2, 'uniroot.all', 'rdrm', 'LL.2', 'L.4'))

r0.eec.p2 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.eec.p2 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.0.5eec.p2 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.0.5eec.p2 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

r0.0.1eec.p2 = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.0.1eec.p2 = matrix(data = NA, nrow = nsims, ncol = length(parfx))

set.seed(0)

#Fill r0 estimates for EEC #####################
  #individual parameters
    r0.eec.p2[, 1] = parSapply(clus1, eec.atr, r0.He, f.f_Nq = fN.atr.bak.uncertainty)[3,]
    r0.eec.p2[, 2] = parSapply(clus1, eec.gly, r0.He, f.f_Nq = fN.gly.bak.uncertainty)[3,]
    r0.eec.p2[, 3] = parSapply(clus1, eec.gly, r0.He, f.mu_Nq = muNq_gly_Bakry12_uncertainty)[3,]
    r0.eec.p2[, 4] = parSapply(clus1, eec.atr, r0.He, f.mu_Nq = muNq_atr_Bakry12_uncertainty)[3,]
    r0.eec.p2[, 5] = parSapply(clus1, eec.mal, r0.In, f.mu_Nq = muNq_mal_tch91_uncertainty)[3,]
    r0.eec.p2[, 6] = parSapply(clus1, eec.mal, r0.In, f.f_Nq = fNq_mal_tch91_uncertainty)[3,]
    r0.eec.p2[, 7] = parSapply(clus1, eec.but, r0.He, f.mu_Nq = mu_Nq_butr_gaf16_uncertainty)[3,]
    r0.eec.p2[, 8] = parSapply(clus1, eec.gly, r0.He, f.mu_Nq = mu_Nq_gly_gaf16_uncertainty)[3,]
    r0.eec.p2[, 9] = parSapply(clus1, eec.pen, r0.He, f.mu_Nq = mu_Nq_pen_gaf16_uncertainty)[3,]
    r0.eec.p2[, 10] = parSapply(clus1, eec.but, r0.He, f.f_Nq = fN.butr.fx.uncertainty)[3,]
    r0.eec.p2[, 11] = parSapply(clus1, eec.gly, r0.He, f.f_Nq = fN.gly.fx.uncertainty)[3,]
    r0.eec.p2[, 12] = parSapply(clus1, eec.pen, r0.He, f.f_Nq = fN.pen.fx.uncertainty)[3,]
    r0.eec.p2[, 13] = parSapply(clus1, eec.mal, r0.In, f.mu_Nq = muNq_mal_Bakry11_uncertainty)[3,]
    r0.eec.p2[, 14] = parSapply(clus1, eec.del, r0.In, f.mu_Nq = muNq_del_Bakry11_uncertainty)[3,]
    r0.eec.p2[, 15] = parSapply(clus1, eec.mal, r0.In, f.f_Nq = fNq_mal_Bakry11_uncertainty)[3,]
    r0.eec.p2[, 16] = parSapply(clus1, eec.del, r0.In, f.f_Nq = fNq_del_Bakry11_uncertainty)[3,]
    r0.eec.p2[, 17] = parSapply(clus1, eec.ch, r0.In, f.f_Nq = fN.hash.chlor.uncertainty)[3,]
    r0.eec.p2[, 18] = parSapply(clus1, eec.prof, r0.In, f.f_Nq = fN.hash.prof.uncertainty)[3,]
    r0.eec.p2[, 19] = parSapply(clus1, eec.ch, r0.In, f.mu_Nq = muNq_ch_hash11_uncertainty)[3,]
    r0.eec.p2[, 20] = parSapply(clus1, eec.prof, r0.In, f.mu_Nq = muNq_prof_hash11_uncertainty)[3,]
    r0.eec.p2[, 21] = parSapply(clus1, eec.ch, r0.In, f.f_Nq = f_N_chlor_ibr92_uncertainty)[3,]
    r0.eec.p2[, 22] = parSapply(clus1, eec.ch, r0.In, f.mu_Nq = mu_N_chlor_ibr92_uncertainty)[3,]
    r0.eec.p2[, 23] = parSapply(clus1, eec.atr, r0.He, f.mu_Nq = ons.munq.atr)[3,]
    r0.eec.p2[, 24] = parSapply(clus1, eec.gly, r0.He, f.mu_Nq = ons.munq.gly)[3,]
    r0.eec.p2[, 25] = parSapply(clus1, eec.but, r0.He, f.mu_Nq = muN.tant.but_uncertainty)[3,]
    r0.eec.p2[, 26] = parSapply(clus1, eec.fpb, r0.He, f.mu_Nq = muN.tant.fpb_uncertainty)[3,]
    r0.eec.p2[, 27] = parSapply(clus1, eec.urea, r0.Fe, f.mu_Nq = rag06_mun_urea)[3,]
    r0.eec.p2[, 28] = parSapply(clus1, eec.pot, r0.Fe, f.mu_Nq = rag06_mun_pot)[3,]
    r0.eec.p2[, 29] = parSapply(clus1, eec.amm, r0.Fe, f.mu_Nq = rag06_mun_amm)[3,]
    
    
  #fill parameter values for EEC values
    par.eec.p2[, 1] = parSapply(clus1, eec.atr, fN.atr.bak.uncertainty) 
    par.eec.p2[, 2] = parSapply(clus1, eec.gly, fN.gly.bak.uncertainty) 
    par.eec.p2[, 3] = parSapply(clus1, eec.gly, muNq_gly_Bakry12_uncertainty) 
    par.eec.p2[, 4] = parSapply(clus1, eec.atr, muNq_atr_Bakry12_uncertainty) 
    par.eec.p2[, 5] = parSapply(clus1, eec.mal, muNq_mal_tch91_uncertainty) 
    par.eec.p2[, 6] = parSapply(clus1, eec.mal, fNq_mal_tch91_uncertainty) 
    par.eec.p2[, 7] = parSapply(clus1, eec.but, mu_Nq_butr_gaf16_uncertainty) 
    par.eec.p2[, 8] = parSapply(clus1, eec.gly, mu_Nq_gly_gaf16_uncertainty) 
    par.eec.p2[, 9] = parSapply(clus1, eec.pen, mu_Nq_pen_gaf16_uncertainty) 
    par.eec.p2[, 10] = parSapply(clus1, eec.but, fN.butr.fx.uncertainty) 
    par.eec.p2[, 11] = parSapply(clus1, eec.gly, fN.gly.fx.uncertainty) 
    par.eec.p2[, 12] = parSapply(clus1, eec.pen, fN.pen.fx.uncertainty) 
    par.eec.p2[, 13] = parSapply(clus1, eec.mal, muNq_mal_Bakry11_uncertainty) 
    par.eec.p2[, 14] = parSapply(clus1, eec.del, muNq_del_Bakry11_uncertainty) 
    par.eec.p2[, 15] = parSapply(clus1, eec.mal, fNq_mal_Bakry11_uncertainty) 
    par.eec.p2[, 16] = parSapply(clus1, eec.del, fNq_del_Bakry11_uncertainty) 
    par.eec.p2[, 17] = parSapply(clus1, eec.ch, fN.hash.chlor.uncertainty) 
    par.eec.p2[, 18] = parSapply(clus1, eec.prof, fN.hash.prof.uncertainty) 
    par.eec.p2[, 19] = parSapply(clus1, eec.ch, muNq_ch_hash11_uncertainty) 
    par.eec.p2[, 20] = parSapply(clus1, eec.prof, muNq_prof_hash11_uncertainty) 
    par.eec.p2[, 21] = parSapply(clus1, eec.ch, f_N_chlor_ibr92_uncertainty) 
    par.eec.p2[, 22] = parSapply(clus1, eec.ch, mu_N_chlor_ibr92_uncertainty) 
    par.eec.p2[, 23] = parSapply(clus1, eec.atr, ons.munq.atr) 
    par.eec.p2[, 24] = parSapply(clus1, eec.gly, ons.munq.gly) 
    par.eec.p2[, 25] = parSapply(clus1, eec.but, muN.tant.but_uncertainty) 
    par.eec.p2[, 26] = parSapply(clus1, eec.fpb, muN.tant.fpb_uncertainty) 
    par.eec.p2[, 27] = parSapply(clus1, eec.urea, rag06_mun_urea) 
    par.eec.p2[, 28] = parSapply(clus1, eec.pot, rag06_mun_pot) 
    par.eec.p2[, 29] = parSapply(clus1, eec.amm, rag06_mun_amm) 
    
#Fill r0 estimates for 50% EEC ########################
  #individual parameters
     r0.0.5eec.p2[, 1] = parSapply(clus1, 0.5*eec.atr, r0.He, f.f_Nq = fN.atr.bak.uncertainty)[3,]
     r0.0.5eec.p2[, 2] = parSapply(clus1, 0.5*eec.gly, r0.He, f.f_Nq = fN.gly.bak.uncertainty)[3,]
     r0.0.5eec.p2[, 3] = parSapply(clus1, 0.5*eec.gly, r0.He, f.mu_Nq = muNq_gly_Bakry12_uncertainty)[3,]
     r0.0.5eec.p2[, 4] = parSapply(clus1, 0.5*eec.atr, r0.He, f.mu_Nq = muNq_atr_Bakry12_uncertainty)[3,]
     r0.0.5eec.p2[, 5] = parSapply(clus1, 0.5*eec.mal, r0.In, f.mu_Nq = muNq_mal_tch91_uncertainty)[3,]
     r0.0.5eec.p2[, 6] = parSapply(clus1, 0.5*eec.mal, r0.In, f.f_Nq = fNq_mal_tch91_uncertainty)[3,]
     r0.0.5eec.p2[, 7] = parSapply(clus1, 0.5*eec.but, r0.He, f.mu_Nq = mu_Nq_butr_gaf16_uncertainty)[3,]
     r0.0.5eec.p2[, 8] = parSapply(clus1, 0.5*eec.gly, r0.He, f.mu_Nq = mu_Nq_gly_gaf16_uncertainty)[3,]
     r0.0.5eec.p2[, 9] = parSapply(clus1, 0.5*eec.pen, r0.He, f.mu_Nq = mu_Nq_pen_gaf16_uncertainty)[3,]
     r0.0.5eec.p2[, 10] = parSapply(clus1, 0.5*eec.but, r0.He, f.f_Nq = fN.butr.fx.uncertainty)[3,]
     r0.0.5eec.p2[, 11] = parSapply(clus1, 0.5*eec.gly, r0.He, f.f_Nq = fN.gly.fx.uncertainty)[3,]
     r0.0.5eec.p2[, 12] = parSapply(clus1, 0.5*eec.pen, r0.He, f.f_Nq = fN.pen.fx.uncertainty)[3,]
     r0.0.5eec.p2[, 13] = parSapply(clus1, 0.5*eec.mal, r0.In, f.mu_Nq = muNq_mal_Bakry11_uncertainty)[3,]
     r0.0.5eec.p2[, 14] = parSapply(clus1, 0.5*eec.del, r0.In, f.mu_Nq = muNq_del_Bakry11_uncertainty)[3,]
     r0.0.5eec.p2[, 15] = parSapply(clus1, 0.5*eec.mal, r0.In, f.f_Nq = fNq_mal_Bakry11_uncertainty)[3,]
     r0.0.5eec.p2[, 16] = parSapply(clus1, 0.5*eec.del, r0.In, f.f_Nq = fNq_del_Bakry11_uncertainty)[3,]
     r0.0.5eec.p2[, 17] = parSapply(clus1, 0.5*eec.ch, r0.In, f.f_Nq = fN.hash.chlor.uncertainty)[3,]
     r0.0.5eec.p2[, 18] = parSapply(clus1, 0.5*eec.prof, r0.In, f.f_Nq = fN.hash.prof.uncertainty)[3,]
     r0.0.5eec.p2[, 19] = parSapply(clus1, 0.5*eec.ch, r0.In, f.mu_Nq = muNq_ch_hash11_uncertainty)[3,]
     r0.0.5eec.p2[, 20] = parSapply(clus1, 0.5*eec.prof, r0.In, f.mu_Nq = muNq_prof_hash11_uncertainty)[3,]
     r0.0.5eec.p2[, 21] = parSapply(clus1, 0.5*eec.ch, r0.In, f.f_Nq = f_N_chlor_ibr92_uncertainty)[3,]
     r0.0.5eec.p2[, 22] = parSapply(clus1, 0.5*eec.ch, r0.In, f.mu_Nq = mu_N_chlor_ibr92_uncertainty)[3,]
     r0.0.5eec.p2[, 23] = parSapply(clus1, 0.5*eec.atr, r0.He, f.mu_Nq = ons.munq.atr)[3,]
     r0.0.5eec.p2[, 24] = parSapply(clus1, 0.5*eec.gly, r0.He, f.mu_Nq = ons.munq.gly)[3,]
     r0.0.5eec.p2[, 25] = parSapply(clus1, 0.5*eec.but, r0.He, f.mu_Nq = muN.tant.but_uncertainty)[3,]
     r0.0.5eec.p2[, 26] = parSapply(clus1, 0.5*eec.fpb, r0.He, f.mu_Nq = muN.tant.fpb_uncertainty)[3,]
     r0.0.5eec.p2[, 27] = parSapply(clus1, 0.5*eec.urea, r0.Fe, f.mu_Nq = rag06_mun_urea)[3,]
     r0.0.5eec.p2[, 28] = parSapply(clus1, 0.5*eec.pot, r0.Fe, f.mu_Nq = rag06_mun_pot)[3,]
     r0.0.5eec.p2[, 29] = parSapply(clus1, 0.5*eec.amm, r0.Fe, f.mu_Nq = rag06_mun_amm)[3,]
    
    
  #fill parameter values for 50% EEC values
    par.0.5eec.p2[, 1] = parSapply(clus1, 0.5*eec.atr, fN.atr.bak.uncertainty) 
    par.0.5eec.p2[, 2] = parSapply(clus1, 0.5*eec.gly, fN.gly.bak.uncertainty) 
    par.0.5eec.p2[, 3] = parSapply(clus1, 0.5*eec.gly, muNq_gly_Bakry12_uncertainty) 
    par.0.5eec.p2[, 4] = parSapply(clus1, 0.5*eec.atr, muNq_atr_Bakry12_uncertainty) 
    par.0.5eec.p2[, 5] = parSapply(clus1, 0.5*eec.mal, muNq_mal_tch91_uncertainty) 
    par.0.5eec.p2[, 6] = parSapply(clus1, 0.5*eec.mal, fNq_mal_tch91_uncertainty) 
    par.0.5eec.p2[, 7] = parSapply(clus1, 0.5*eec.but, mu_Nq_butr_gaf16_uncertainty) 
    par.0.5eec.p2[, 8] = parSapply(clus1, 0.5*eec.gly, mu_Nq_gly_gaf16_uncertainty) 
    par.0.5eec.p2[, 9] = parSapply(clus1, 0.5*eec.pen, mu_Nq_pen_gaf16_uncertainty) 
    par.0.5eec.p2[, 10] = parSapply(clus1, 0.5*eec.but, fN.butr.fx.uncertainty) 
    par.0.5eec.p2[, 11] = parSapply(clus1, 0.5*eec.gly, fN.gly.fx.uncertainty) 
    par.0.5eec.p2[, 12] = parSapply(clus1, 0.5*eec.pen, fN.pen.fx.uncertainty) 
    par.0.5eec.p2[, 13] = parSapply(clus1, 0.5*eec.mal, muNq_mal_Bakry11_uncertainty) 
    par.0.5eec.p2[, 14] = parSapply(clus1, 0.5*eec.del, muNq_del_Bakry11_uncertainty) 
    par.0.5eec.p2[, 15] = parSapply(clus1, 0.5*eec.mal, fNq_mal_Bakry11_uncertainty) 
    par.0.5eec.p2[, 16] = parSapply(clus1, 0.5*eec.del, fNq_del_Bakry11_uncertainty) 
    par.0.5eec.p2[, 17] = parSapply(clus1, 0.5*eec.ch, fN.hash.chlor.uncertainty) 
    par.0.5eec.p2[, 18] = parSapply(clus1, 0.5*eec.prof, fN.hash.prof.uncertainty) 
    par.0.5eec.p2[, 19] = parSapply(clus1, 0.5*eec.ch, muNq_ch_hash11_uncertainty) 
    par.0.5eec.p2[, 20] = parSapply(clus1, 0.5*eec.prof, muNq_prof_hash11_uncertainty) 
    par.0.5eec.p2[, 21] = parSapply(clus1, 0.5*eec.ch, f_N_chlor_ibr92_uncertainty) 
    par.0.5eec.p2[, 22] = parSapply(clus1, 0.5*eec.ch, mu_N_chlor_ibr92_uncertainty) 
    par.0.5eec.p2[, 23] = parSapply(clus1, 0.5*eec.atr, ons.munq.atr) 
    par.0.5eec.p2[, 24] = parSapply(clus1, 0.5*eec.gly, ons.munq.gly) 
    par.0.5eec.p2[, 25] = parSapply(clus1, 0.5*eec.but, muN.tant.but_uncertainty) 
    par.0.5eec.p2[, 26] = parSapply(clus1, 0.5*eec.fpb, muN.tant.fpb_uncertainty) 
    par.0.5eec.p2[, 27] = parSapply(clus1, 0.5*eec.urea, rag06_mun_urea) 
    par.0.5eec.p2[, 28] = parSapply(clus1, 0.5*eec.pot, rag06_mun_pot) 
    par.0.5eec.p2[, 29] = parSapply(clus1, 0.5*eec.amm, rag06_mun_amm) 
    
#Fill r0 estimates for 10% EEC #####################
  #individual parameters
    r0.0.1eec.p2[, 1] = parSapply(clus1, 0.1*eec.atr, r0.He, f.f_Nq = fN.atr.bak.uncertainty)[3,]
    r0.0.1eec.p2[, 2] = parSapply(clus1, 0.1*eec.gly, r0.He, f.f_Nq = fN.gly.bak.uncertainty)[3,]
    r0.0.1eec.p2[, 3] = parSapply(clus1, 0.1*eec.gly, r0.He, f.mu_Nq = muNq_gly_Bakry12_uncertainty)[3,]
    r0.0.1eec.p2[, 4] = parSapply(clus1, 0.1*eec.atr, r0.He, f.mu_Nq = muNq_atr_Bakry12_uncertainty)[3,]
    r0.0.1eec.p2[, 5] = parSapply(clus1, 0.1*eec.mal, r0.In, f.mu_Nq = muNq_mal_tch91_uncertainty)[3,]
    r0.0.1eec.p2[, 6] = parSapply(clus1, 0.1*eec.mal, r0.In, f.f_Nq = fNq_mal_tch91_uncertainty)[3,]
    r0.0.1eec.p2[, 7] = parSapply(clus1, 0.1*eec.but, r0.He, f.mu_Nq = mu_Nq_butr_gaf16_uncertainty)[3,]
    r0.0.1eec.p2[, 8] = parSapply(clus1, 0.1*eec.gly, r0.He, f.mu_Nq = mu_Nq_gly_gaf16_uncertainty)[3,]
    r0.0.1eec.p2[, 9] = parSapply(clus1, 0.1*eec.pen, r0.He, f.mu_Nq = mu_Nq_pen_gaf16_uncertainty)[3,]
    r0.0.1eec.p2[, 10] = parSapply(clus1, 0.1*eec.but, r0.He, f.f_Nq = fN.butr.fx.uncertainty)[3,]
    r0.0.1eec.p2[, 11] = parSapply(clus1, 0.1*eec.gly, r0.He, f.f_Nq = fN.gly.fx.uncertainty)[3,]
    r0.0.1eec.p2[, 12] = parSapply(clus1, 0.1*eec.pen, r0.He, f.f_Nq = fN.pen.fx.uncertainty)[3,]
    r0.0.1eec.p2[, 13] = parSapply(clus1, 0.1*eec.mal, r0.In, f.mu_Nq = muNq_mal_Bakry11_uncertainty)[3,]
    r0.0.1eec.p2[, 14] = parSapply(clus1, 0.1*eec.del, r0.In, f.mu_Nq = muNq_del_Bakry11_uncertainty)[3,]
    r0.0.1eec.p2[, 15] = parSapply(clus1, 0.1*eec.mal, r0.In, f.f_Nq = fNq_mal_Bakry11_uncertainty)[3,]
    r0.0.1eec.p2[, 16] = parSapply(clus1, 0.1*eec.del, r0.In, f.f_Nq = fNq_del_Bakry11_uncertainty)[3,]
    r0.0.1eec.p2[, 17] = parSapply(clus1, 0.1*eec.ch, r0.In, f.f_Nq = fN.hash.chlor.uncertainty)[3,]
    r0.0.1eec.p2[, 18] = parSapply(clus1, 0.1*eec.prof, r0.In, f.f_Nq = fN.hash.prof.uncertainty)[3,]
    r0.0.1eec.p2[, 19] = parSapply(clus1, 0.1*eec.ch, r0.In, f.mu_Nq = muNq_ch_hash11_uncertainty)[3,]
    r0.0.1eec.p2[, 20] = parSapply(clus1, 0.1*eec.prof, r0.In, f.mu_Nq = muNq_prof_hash11_uncertainty)[3,]
    r0.0.1eec.p2[, 21] = parSapply(clus1, 0.1*eec.ch, r0.In, f.f_Nq = f_N_chlor_ibr92_uncertainty)[3,]
    r0.0.1eec.p2[, 22] = parSapply(clus1, 0.1*eec.ch, r0.In, f.mu_Nq = mu_N_chlor_ibr92_uncertainty)[3,]
    r0.0.1eec.p2[, 23] = parSapply(clus1, 0.1*eec.atr, r0.He, f.mu_Nq = ons.munq.atr)[3,]
    r0.0.1eec.p2[, 24] = parSapply(clus1, 0.1*eec.gly, r0.He, f.mu_Nq = ons.munq.gly)[3,]
    r0.0.1eec.p2[, 25] = parSapply(clus1, 0.1*eec.but, r0.He, f.mu_Nq = muN.tant.but_uncertainty)[3,]
    r0.0.1eec.p2[, 26] = parSapply(clus1, 0.1*eec.fpb, r0.He, f.mu_Nq = muN.tant.fpb_uncertainty)[3,]
    r0.0.1eec.p2[, 27] = parSapply(clus1, 0.1*eec.urea, r0.Fe, f.mu_Nq = rag06_mun_urea)[3,]
    r0.0.1eec.p2[, 28] = parSapply(clus1, 0.1*eec.pot, r0.Fe, f.mu_Nq = rag06_mun_pot)[3,]
    r0.0.1eec.p2[, 29] = parSapply(clus1, 0.1*eec.amm, r0.Fe, f.mu_Nq = rag06_mun_amm)[3,]
    
    
  #fill parameter values for 50% EEC values
    par.0.1eec.p2[, 1] = parSapply(clus1, 0.1*eec.atr, fN.atr.bak.uncertainty) 
    par.0.1eec.p2[, 2] = parSapply(clus1, 0.1*eec.gly, fN.gly.bak.uncertainty) 
    par.0.1eec.p2[, 3] = parSapply(clus1, 0.1*eec.gly, muNq_gly_Bakry12_uncertainty) 
    par.0.1eec.p2[, 4] = parSapply(clus1, 0.1*eec.atr, muNq_atr_Bakry12_uncertainty) 
    par.0.1eec.p2[, 5] = parSapply(clus1, 0.1*eec.mal, muNq_mal_tch91_uncertainty) 
    par.0.1eec.p2[, 6] = parSapply(clus1, 0.1*eec.mal, fNq_mal_tch91_uncertainty) 
    par.0.1eec.p2[, 7] = parSapply(clus1, 0.1*eec.but, mu_Nq_butr_gaf16_uncertainty) 
    par.0.1eec.p2[, 8] = parSapply(clus1, 0.1*eec.gly, mu_Nq_gly_gaf16_uncertainty) 
    par.0.1eec.p2[, 9] = parSapply(clus1, 0.1*eec.pen, mu_Nq_pen_gaf16_uncertainty) 
    par.0.1eec.p2[, 10] = parSapply(clus1, 0.1*eec.but, fN.butr.fx.uncertainty) 
    par.0.1eec.p2[, 11] = parSapply(clus1, 0.1*eec.gly, fN.gly.fx.uncertainty) 
    par.0.1eec.p2[, 12] = parSapply(clus1, 0.1*eec.pen, fN.pen.fx.uncertainty) 
    par.0.1eec.p2[, 13] = parSapply(clus1, 0.1*eec.mal, muNq_mal_Bakry11_uncertainty) 
    par.0.1eec.p2[, 14] = parSapply(clus1, 0.1*eec.del, muNq_del_Bakry11_uncertainty) 
    par.0.1eec.p2[, 15] = parSapply(clus1, 0.1*eec.mal, fNq_mal_Bakry11_uncertainty) 
    par.0.1eec.p2[, 16] = parSapply(clus1, 0.1*eec.del, fNq_del_Bakry11_uncertainty) 
    par.0.1eec.p2[, 17] = parSapply(clus1, 0.1*eec.ch, fN.hash.chlor.uncertainty) 
    par.0.1eec.p2[, 18] = parSapply(clus1, 0.1*eec.prof, fN.hash.prof.uncertainty) 
    par.0.1eec.p2[, 19] = parSapply(clus1, 0.1*eec.ch, muNq_ch_hash11_uncertainty) 
    par.0.1eec.p2[, 20] = parSapply(clus1, 0.1*eec.prof, muNq_prof_hash11_uncertainty) 
    par.0.1eec.p2[, 21] = parSapply(clus1, 0.1*eec.ch, f_N_chlor_ibr92_uncertainty) 
    par.0.1eec.p2[, 22] = parSapply(clus1, 0.1*eec.ch, mu_N_chlor_ibr92_uncertainty) 
    par.0.1eec.p2[, 23] = parSapply(clus1, 0.1*eec.atr, ons.munq.atr) 
    par.0.1eec.p2[, 24] = parSapply(clus1, 0.1*eec.gly, ons.munq.gly) 
    par.0.1eec.p2[, 25] = parSapply(clus1, 0.1*eec.but, muN.tant.but_uncertainty) 
    par.0.1eec.p2[, 26] = parSapply(clus1, 0.1*eec.fpb, muN.tant.fpb_uncertainty) 
    par.0.1eec.p2[, 27] = parSapply(clus1, 0.1*eec.urea, rag06_mun_urea) 
    par.0.1eec.p2[, 28] = parSapply(clus1, 0.1*eec.pot, rag06_mun_pot) 
    par.0.1eec.p2[, 29] = parSapply(clus1, 0.1*eec.amm, rag06_mun_amm) 

stopCluster(clus1) 
#Post process ############ 
#EEC runs ################
eec.p2.df = data.frame(chem = c('Atrazine', 'Glyphosate', 'Glyphosate', 'Atrazine', 'Malathion', 'Malathion',
                                'Butraline', 'Glyphosate', 'Pendimethalin', 'Butraline', 'Glyphosate', 'Pendimethalin', 
                                'Malathion', 'Deltamethrin', 'Malathion', 'Deltamethrin', 'Chlorpyrifos', 'Profenofos',
                                'Chlorpyrifos', 'Profenofos', 'Chlorpyrifos', 'Chlorpyrifos', 'Atrazine', 'Glyphosate',
                                'Butralin', 'Fluazifop-p-butyl', 'Urea', 'Potassium Sulphate', 'Ammonium Nitrate'),
                       study = c(rep('Bakry et al 2012', 4), rep('Tchounwou et al 1991', 2), rep('Abdel-Ghaffar et al 2016', 6),
                                 rep('Bakry et al 2011', 4), rep('Hasheesh & Mohamed 2011', 4), rep('Ibrahim et al 1992', 2),
                                 rep('Omran & Salama 2013', 2), rep('Tantawy 2002', 2), rep('Ragab & Shoukry 2006', 3)),
                       Species = c(rep('Biomphalaria alexandrina', 4), rep('Bulinus havenensis',2), rep('Biomphalaria alexandrina', 6),
                                   rep('Helisoma duryi', 4), rep('Bulinus truncatus', 4), rep('Biomphalaria alexandrina', 9)),
                       Parameter = c(rep('fN', 2), rep('muN', 3), 'fN', rep('muN', 3), rep('fN', 3), rep('muN', 2),
                                     rep('fN', 4), rep('muN', 2), 'fN', rep('muN', 8)),
                       r0 = colMeans(r0.eec.p2),
                       r0.sd = apply(r0.eec.p2, 2, sd),
                       par.mean = colMeans(par.eec.p2))

eec.p2.df$r0.up = eec.p2.df$r0 + eec.p2.df$r0.sd
eec.p2.df$r0.lo = eec.p2.df$r0 - eec.p2.df$r0.sd

save(eec.p2.df, file = 'Review_models/r0_EECs/eec.p2.df.RData')

#50% EEC values #################
eec0.5.p2.df = data.frame(chem = c('Atrazine', 'Glyphosate', 'Glyphosate', 'Atrazine', 'Malathion', 'Malathion',
                                   'Butraline', 'Glyphosate', 'Pendimethalin', 'Butraline', 'Glyphosate', 'Pendimethalin', 
                                   'Malathion', 'Deltamethrin', 'Malathion', 'Deltamethrin', 'Chlorpyrifos', 'Profenofos',
                                   'Chlorpyrifos', 'Profenofos', 'Chlorpyrifos', 'Chlorpyrifos', 'Atrazine', 'Glyphosate',
                                   'Butralin', 'Fluazifop-p-butyl', 'Urea', 'Potassium Sulphate', 'Ammonium Nitrate'),
                          study = c(rep('Bakry et al 2012', 4), rep('Tchounwou et al 1991', 2), rep('Abdel-Ghaffar et al 2016', 6),
                                    rep('Bakry et al 2011', 4), rep('Hasheesh & Mohamed 2011', 4), rep('Ibrahim et al 1992', 2),
                                    rep('Omran & Salama 2013', 2), rep('Tantawy 2002', 2), rep('Ragab & Shoukry 2006', 3)),
                          Species = c(rep('Biomphalaria alexandrina', 4), rep('Bulinus havenensis',2), rep('Biomphalaria alexandrina', 6),
                                      rep('Helisoma duryi', 4), rep('Bulinus truncatus', 4), rep('Biomphalaria alexandrina', 9)),
                          Parameter = c(rep('fN', 2), rep('muN', 3), 'fN', rep('muN', 3), rep('fN', 3), rep('muN', 2),
                                        rep('fN', 4), rep('muN', 2), 'fN', rep('muN', 8)),
                          r0 = colMeans(r0.0.5eec.p2),
                          r0.sd = apply(r0.0.5eec.p2, 2, sd),
                          par.mean = colMeans(par.0.5eec.p2))

eec0.5.p2.df$r0.up = eec0.5.p2.df$r0 + eec0.5.p2.df$r0.sd
eec0.5.p2.df$r0.lo = eec0.5.p2.df$r0 - eec0.5.p2.df$r0.sd

save(eec0.5.p2.df, file = 'Review_models/r0_EECs/eec0.5.p2.df.RData')


#10% EEC values ######################
eec0.1.p2.df = data.frame(chem = c('Atrazine', 'Glyphosate', 'Glyphosate', 'Atrazine', 'Malathion', 'Malathion',
                                   'Butraline', 'Glyphosate', 'Pendimethalin', 'Butraline', 'Glyphosate', 'Pendimethalin', 
                                   'Malathion', 'Deltamethrin', 'Malathion', 'Deltamethrin', 'Chlorpyrifos', 'Profenofos',
                                   'Chlorpyrifos', 'Profenofos', 'Chlorpyrifos', 'Chlorpyrifos', 'Atrazine', 'Glyphosate',
                                   'Butralin', 'Fluazifop-p-butyl', 'Urea', 'Potassium Sulphate', 'Ammonium Nitrate'),
                          study = c(rep('Bakry et al 2012', 4), rep('Tchounwou et al 1991', 2), rep('Abdel-Ghaffar et al 2016', 6),
                                    rep('Bakry et al 2011', 4), rep('Hasheesh & Mohamed 2011', 4), rep('Ibrahim et al 1992', 2),
                                    rep('Omran & Salama 2013', 2), rep('Tantawy 2002', 2), rep('Ragab & Shoukry 2006', 3)),
                          Species = c(rep('Biomphalaria alexandrina', 4), rep('Bulinus havenensis',2), rep('Biomphalaria alexandrina', 6),
                                      rep('Helisoma duryi', 4), rep('Bulinus truncatus', 4), rep('Biomphalaria alexandrina', 9)),
                          Parameter = c(rep('fN', 2), rep('muN', 3), 'fN', rep('muN', 3), rep('fN', 3), rep('muN', 2),
                                        rep('fN', 4), rep('muN', 2), 'fN', rep('muN', 8)),
                          r0 = colMeans(r0.0.1eec.p2),
                          r0.sd = apply(r0.0.1eec.p2, 2, sd),
                          par.mean = colMeans(par.0.1eec.p2))

eec0.1.p2.df$r0.up = eec0.1.p2.df$r0 + eec0.1.p2.df$r0.sd
eec0.1.p2.df$r0.lo = eec0.1.p2.df$r0 - eec0.1.p2.df$r0.sd

save(eec0.1.p2.df, file = 'Review_models/r0_EECs/eec0.1.p2.df.RData')
