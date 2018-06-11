#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

library(parallel)
library(fBasics)
library(tidyverse)

#Load all response functions ######
rfx_files <- list.files(path = "Agrochemical_Review/Response_Fxs",
                        pattern = "_fit.R", recursive = TRUE)
rfx_dir <- "Agrochemical_Review/Response_Fxs/"

sapply(paste0(rfx_dir, rfx_files), source)

#Prep simulations ############
nsims = 5000
no.cores = detectCores() - 1

#Expected environmental concentrations for agrochemicals
eec.atr = rep(102, nsims)     #atrazine, EEC from Halstead et al 2017
blnk.vec = c(1:nsims)
eec.gly = rep(3700, nsims)    #Glyphosate eec from Raffel et al
eec.ch = rep(64, nsims)       #chlorpyrifos eec from halstead et al 2015
eec.mal = rep(18.4, nsims)     #malathion eec from hasltead et al 2015
eec.but = rep(2.12, nsims)     #butachlor eec form lucy
eec.btr = rep(16.89, nsims)     #butralin eec from epa report
eec.pen = rep(11.4, nsims)     #pendimethalin eec from lucy
eec.del = rep(0.0052, nsims)     #deltamethrin eec from lucy
eec.prof = rep(3.01, nsims)    #profenofos eec from lucy
eec.fpb = rep(4.67, nsims)     #fluazifop-p-butyl eec form lucy
eec.urea = rep(100, nsims)    #urea eec TBD
eec.pot = rep(100, nsims)     #potassium sulphate eec TBD
eec.amm = rep(100, nsims)     #Ammonium phosphate eec TBD
eec.dim = rep(4.07, nsims)     #dimethoate eec from lucy
eec.esf = rep(1.03, nsims)    #esfenvalerate eec from halstead et al 2015
eec.lcy = rep(1.77, nsims)    #lamda-cyhalothrin eec from halstead et al 2015
eec.per = rep(5.98, nsims)    #permethrin eec from halstead et al 2015
eec.pro = rep(3.01, nsims)     #profenofos eec form lucy
eec.trb = rep(36.6, nsims)    #terbufos eec from halstead et al 2015
eec.znc = rep(100, nsims)     #zinc sulfate eec TBD

parfx = c(phi_Nq_atr_baxrohr.no30, phi_Nq_atr_baxrohr, johnson07_theta_uncertainty, johnson07_phin_uncertainty, halstead17_phiN_at_uncertainty, halstead17_phiN_fe_uncertainty, rohr08_fN_uncertainty2, fN.atr.bak.uncertainty, fN.gly.bak.uncertainty, muNq_gly_Bakry12_uncertainty, muNq_atr_Bakry12_uncertainty, muNq_mal_tch91_uncertainty, fNq_mal_tch91_uncertainty, mu_Nq_butr_gaf16_uncertainty, mu_Nq_gly_gaf16_uncertainty, mu_Nq_pen_gaf16_uncertainty, fN.butr.fx.uncertainty, fN.gly.fx.uncertainty, fN.pen.fx.uncertainty, muNq_mal_Bakry11_uncertainty, muNq_del_Bakry11_uncertainty, fNq_mal_Bakry11_uncertainty, fNq_del_Bakry11_uncertainty, fN.hash.chlor.uncertainty, fN.hash.prof.uncertainty, muNq_ch_hash11_uncertainty, muNq_prof_hash11_uncertainty, f_N_chlor_ibr92_uncertainty, mu_N_chlor_ibr92_uncertainty, ons.munq.atr, ons.munq.gly, muN.tant.but_uncertainty, muN.tant.fpb_uncertainty, rag06_mun_urea, rag06_mun_pot, rag06_mun_amm, muPq_chlor_satapornvanit09_uncertainty, muPq_chlor_Halstead_uncertainty, muPq_dim_satapornvanit09_uncertainty, muPq_esfen_Halstead_uncertainty, muPq_lamcy_Halstead_uncertainty, muPq_mal_Halstead_uncertainty, muPq_perm_Halstead_uncertainty, muPq_prof_satapornvanit09_uncertainty, muPq_terb_Halstead_uncertainty, muPq_zinc_satapornvanit09_uncertainty, psi_q_chlor_satapornvanit09_uncertainty, psi_q_zinc_satapornvanit09_uncertainty, piC.grg08_atr_unc2, piC.kop_atr_unc, piC.atr.rohr08.lin, piC.meta_atr_unc, piC.tant02_but.exp_unc, piC.tant02_fpb.exp_unc, piC.ghaf_butr.exp_unc, piC.ghaf_gly.exp_unc, piC.ghaf_pen.exp_unc, piC.tch92_mal_unc, piC_ch_Hash11_uncertainty, piC_pr_Hash11_uncertainty, piM_ch_Hash11_uncertainty, piM_pr_Hash11_uncertainty, piM.ghaf_butr.exp_unc, piM.ghaf_gly.exp_unc, piM.ghaf_pen.exp_unc, piM.tch91_mal_unc,piM.tant02_but.exp_unc, piM.tant02_fpb.exp_unc, piM.tch91_amm_unc, piM.tch91_ure_unc, tch91.egv.amm_unc, tch91.egv.ure_unc)  

clus1 = makeCluster(no.cores)
clusterExport(clus1, c(keep.all, 'uniroot.all', 'rdrm', 'LL.2', 'LL.3', 'L.4'))

#Simulation settings
nsims = 5000
seed = 43093

par.eec.all = matrix(data = NA, nrow = nsims, ncol = length(parfx))
par.0.1eec.all = matrix(data = NA, nrow = nsims, ncol = length(parfx))

clusterSetRNGStream(cl = clus1, iseed = seed) #set cluster seed for reproducibility

#fill parameter values for each response function #########
#Pathway 1
# EEC
  par.eec.all[, 1] = parSapply(clus1, eec.atr, phi_Nq_atr_baxrohr.no30) 
  par.eec.all[, 2] = parSapply(clus1, eec.atr, phi_Nq_atr_baxrohr) 
  par.eec.all[, 3] = parSapply(clus1, blnk.vec, johnson07_phin_uncertainty) 
  par.eec.all[, 4] = parSapply(clus1, blnk.vec, johnson07_theta_uncertainty) 
  par.eec.all[, 5] = parSapply(clus1, blnk.vec, rohr08_fN_uncertainty2) 
  par.eec.all[, 6] = parSapply(clus1, blnk.vec, halstead17_phiN_at_uncertainty) 
  par.eec.all[, 7] = parSapply(clus1, blnk.vec, halstead17_phiN_fe_uncertainty) 
#0.1EEC  
  par.0.1eec.all[, 1] = parSapply(clus1, 0.1*eec.atr, phi_Nq_atr_baxrohr.no30) 
  par.0.1eec.all[, 2] = parSapply(clus1, 0.1*eec.atr, phi_Nq_atr_baxrohr) 
  par.0.1eec.all[, 3] = parSapply(clus1, blnk.vec, johnson07_phin_uncertainty) 
  par.0.1eec.all[, 4] = parSapply(clus1, blnk.vec, johnson07_theta_uncertainty) 
  par.0.1eec.all[, 5] = parSapply(clus1, blnk.vec, rohr08_fN_uncertainty2) 
  par.0.1eec.all[, 6] = parSapply(clus1, blnk.vec, halstead17_phiN_at_uncertainty) 
  par.0.1eec.all[, 7] = parSapply(clus1, blnk.vec, halstead17_phiN_fe_uncertainty) 
  
#Pathway 2
#EEC
  par.eec.all[, 8] = parSapply(clus1, eec.atr, fN.atr.bak.uncertainty) 
  par.eec.all[, 9] = parSapply(clus1, eec.gly, fN.gly.bak.uncertainty) 
  par.eec.all[, 10] = parSapply(clus1, eec.gly, muNq_gly_Bakry12_uncertainty) 
  par.eec.all[, 11] = parSapply(clus1, eec.atr, muNq_atr_Bakry12_uncertainty) 
  par.eec.all[, 12] = parSapply(clus1, eec.mal, muNq_mal_tch91_uncertainty) 
  par.eec.all[, 13] = parSapply(clus1, eec.mal, fNq_mal_tch91_uncertainty) 
  par.eec.all[, 14] = parSapply(clus1, eec.btr, mu_Nq_butr_gaf16_uncertainty) 
  par.eec.all[, 15] = parSapply(clus1, eec.gly, mu_Nq_gly_gaf16_uncertainty) 
  par.eec.all[, 16] = parSapply(clus1, eec.pen, mu_Nq_pen_gaf16_uncertainty) 
  par.eec.all[, 17] = parSapply(clus1, eec.btr, fN.butr.fx.uncertainty) 
  par.eec.all[, 18] = parSapply(clus1, eec.gly, fN.gly.fx.uncertainty) 
  par.eec.all[, 19] = parSapply(clus1, eec.pen, fN.pen.fx.uncertainty) 
  par.eec.all[, 20] = parSapply(clus1, eec.mal, muNq_mal_Bakry11_uncertainty) 
  par.eec.all[, 21] = parSapply(clus1, eec.del, muNq_del_Bakry11_uncertainty) 
  par.eec.all[, 22] = parSapply(clus1, eec.mal, fNq_mal_Bakry11_uncertainty) 
  par.eec.all[, 23] = parSapply(clus1, eec.del, fNq_del_Bakry11_uncertainty) 
  par.eec.all[, 24] = parSapply(clus1, eec.ch, fN.hash.chlor.uncertainty) 
  par.eec.all[, 25] = parSapply(clus1, eec.prof, fN.hash.prof.uncertainty) 
  par.eec.all[, 26] = parSapply(clus1, eec.ch, muNq_ch_hash11_uncertainty) 
  par.eec.all[, 27] = parSapply(clus1, eec.prof, muNq_prof_hash11_uncertainty) 
  par.eec.all[, 28] = parSapply(clus1, eec.ch, f_N_chlor_ibr92_uncertainty) 
  par.eec.all[, 29] = parSapply(clus1, eec.ch, mu_N_chlor_ibr92_uncertainty) 
  par.eec.all[, 30] = parSapply(clus1, eec.atr, ons.munq.atr) 
  par.eec.all[, 31] = parSapply(clus1, eec.gly, ons.munq.gly) 
  par.eec.all[, 32] = parSapply(clus1, eec.but, muN.tant.but_uncertainty) 
  par.eec.all[, 33] = parSapply(clus1, eec.fpb, muN.tant.fpb_uncertainty) 
  par.eec.all[, 34] = parSapply(clus1, eec.urea, rag06_mun_urea) 
  par.eec.all[, 35] = parSapply(clus1, eec.pot, rag06_mun_pot) 
  par.eec.all[, 36] = parSapply(clus1, eec.amm, rag06_mun_amm) 
#0.1 EEC
  par.0.1eec.all[, 8] = parSapply(clus1, 0.1*eec.atr, fN.atr.bak.uncertainty) 
  par.0.1eec.all[, 9] = parSapply(clus1, 0.1*eec.gly, fN.gly.bak.uncertainty) 
  par.0.1eec.all[, 10] = parSapply(clus1, 0.1*eec.gly, muNq_gly_Bakry12_uncertainty) 
  par.0.1eec.all[, 11] = parSapply(clus1, 0.1*eec.atr, muNq_atr_Bakry12_uncertainty) 
  par.0.1eec.all[, 12] = parSapply(clus1, 0.1*eec.mal, muNq_mal_tch91_uncertainty) 
  par.0.1eec.all[, 13] = parSapply(clus1, 0.1*eec.mal, fNq_mal_tch91_uncertainty) 
  par.0.1eec.all[, 14] = parSapply(clus1, 0.1*eec.btr, mu_Nq_butr_gaf16_uncertainty) 
  par.0.1eec.all[, 15] = parSapply(clus1, 0.1*eec.gly, mu_Nq_gly_gaf16_uncertainty) 
  par.0.1eec.all[, 16] = parSapply(clus1, 0.1*eec.pen, mu_Nq_pen_gaf16_uncertainty) 
  par.0.1eec.all[, 17] = parSapply(clus1, 0.1*eec.btr, fN.butr.fx.uncertainty) 
  par.0.1eec.all[, 18] = parSapply(clus1, 0.1*eec.gly, fN.gly.fx.uncertainty) 
  par.0.1eec.all[, 19] = parSapply(clus1, 0.1*eec.pen, fN.pen.fx.uncertainty) 
  par.0.1eec.all[, 20] = parSapply(clus1, 0.1*eec.mal, muNq_mal_Bakry11_uncertainty) 
  par.0.1eec.all[, 21] = parSapply(clus1, 0.1*eec.del, muNq_del_Bakry11_uncertainty) 
  par.0.1eec.all[, 22] = parSapply(clus1, 0.1*eec.mal, fNq_mal_Bakry11_uncertainty) 
  par.0.1eec.all[, 23] = parSapply(clus1, 0.1*eec.del, fNq_del_Bakry11_uncertainty) 
  par.0.1eec.all[, 24] = parSapply(clus1, 0.1*eec.ch, fN.hash.chlor.uncertainty) 
  par.0.1eec.all[, 25] = parSapply(clus1, 0.1*eec.prof, fN.hash.prof.uncertainty) 
  par.0.1eec.all[, 26] = parSapply(clus1, 0.1*eec.ch, muNq_ch_hash11_uncertainty) 
  par.0.1eec.all[, 27] = parSapply(clus1, 0.1*eec.prof, muNq_prof_hash11_uncertainty) 
  par.0.1eec.all[, 28] = parSapply(clus1, 0.1*eec.ch, f_N_chlor_ibr92_uncertainty) 
  par.0.1eec.all[, 29] = parSapply(clus1, 0.1*eec.ch, mu_N_chlor_ibr92_uncertainty) 
  par.0.1eec.all[, 30] = parSapply(clus1, 0.1*eec.atr, ons.munq.atr) 
  par.0.1eec.all[, 31] = parSapply(clus1, 0.1*eec.gly, ons.munq.gly) 
  par.0.1eec.all[, 32] = parSapply(clus1, 0.1*eec.but, muN.tant.but_uncertainty) 
  par.0.1eec.all[, 33] = parSapply(clus1, 0.1*eec.fpb, muN.tant.fpb_uncertainty) 
  par.0.1eec.all[, 34] = parSapply(clus1, 0.1*eec.urea, rag06_mun_urea) 
  par.0.1eec.all[, 35] = parSapply(clus1, 0.1*eec.pot, rag06_mun_pot) 
  par.0.1eec.all[, 36] = parSapply(clus1, 0.1*eec.amm, rag06_mun_amm) 
  
#Pathway 3
#EEC
  par.eec.all[, 37] = parSapply(clus1, eec.ch, muPq_chlor_satapornvanit09_uncertainty) 
  par.eec.all[, 38] = parSapply(clus1, eec.ch, muPq_chlor_Halstead_uncertainty) 
  par.eec.all[, 39] = parSapply(clus1, eec.dim, muPq_dim_satapornvanit09_uncertainty) 
  par.eec.all[, 40] = parSapply(clus1, eec.esf, muPq_esfen_Halstead_uncertainty) 
  par.eec.all[, 41] = parSapply(clus1, eec.lcy, muPq_lamcy_Halstead_uncertainty) 
  par.eec.all[, 42] = parSapply(clus1, eec.mal, muPq_mal_Halstead_uncertainty) 
  par.eec.all[, 43] = parSapply(clus1, eec.per, muPq_perm_Halstead_uncertainty) 
  par.eec.all[, 44] = parSapply(clus1, eec.pro, muPq_prof_satapornvanit09_uncertainty) 
  par.eec.all[, 45] = parSapply(clus1, eec.trb, muPq_terb_Halstead_uncertainty) 
  par.eec.all[, 46] = parSapply(clus1, eec.znc, muPq_zinc_satapornvanit09_uncertainty) 
  par.eec.all[, 47] = parSapply(clus1, eec.ch, psi_q_chlor_satapornvanit09_uncertainty) 
  par.eec.all[, 48] = parSapply(clus1, eec.znc, psi_q_zinc_satapornvanit09_uncertainty) 
#0.1 EEC
  par.0.1eec.all[, 37] = parSapply(clus1, 0.1*eec.ch, muPq_chlor_satapornvanit09_uncertainty) 
  par.0.1eec.all[, 38] = parSapply(clus1, 0.1*eec.ch, muPq_chlor_Halstead_uncertainty) 
  par.0.1eec.all[, 39] = parSapply(clus1, 0.1*eec.dim, muPq_dim_satapornvanit09_uncertainty) 
  par.0.1eec.all[, 40] = parSapply(clus1, 0.1*eec.esf, muPq_esfen_Halstead_uncertainty) 
  par.0.1eec.all[, 41] = parSapply(clus1, 0.1*eec.lcy, muPq_lamcy_Halstead_uncertainty) 
  par.0.1eec.all[, 42] = parSapply(clus1, 0.1*eec.mal, muPq_mal_Halstead_uncertainty) 
  par.0.1eec.all[, 43] = parSapply(clus1, 0.1*eec.per, muPq_perm_Halstead_uncertainty) 
  par.0.1eec.all[, 44] = parSapply(clus1, 0.1*eec.pro, muPq_prof_satapornvanit09_uncertainty) 
  par.0.1eec.all[, 45] = parSapply(clus1, 0.1*eec.trb, muPq_terb_Halstead_uncertainty) 
  par.0.1eec.all[, 46] = parSapply(clus1, 0.1*eec.znc, muPq_zinc_satapornvanit09_uncertainty) 
  par.0.1eec.all[, 47] = parSapply(clus1, 0.1*eec.ch, psi_q_chlor_satapornvanit09_uncertainty) 
  par.0.1eec.all[, 48] = parSapply(clus1, 0.1*eec.znc, psi_q_zinc_satapornvanit09_uncertainty) 
  
#Pathway 4
#EEC
  par.eec.all[, 49] = parSapply(clus1, eec.atr, piC.grg08_atr_unc2) 
  par.eec.all[, 50] = parSapply(clus1, eec.atr, piC.kop_atr_unc) 
  par.eec.all[, 51] = parSapply(clus1, eec.atr, piC.atr.rohr08.lin) 
  par.eec.all[, 52] = parSapply(clus1, eec.atr, piC.meta_atr_unc) 
  par.eec.all[, 53] = parSapply(clus1, eec.but, piC.tant02_but.exp_unc) 
  par.eec.all[, 54] = parSapply(clus1, eec.fpb, piC.tant02_fpb.exp_unc) 
  par.eec.all[, 55] = parSapply(clus1, eec.btr, piC.ghaf_butr.exp_unc) 
  par.eec.all[, 56] = parSapply(clus1, eec.gly, piC.ghaf_gly.exp_unc) 
  par.eec.all[, 57] = parSapply(clus1, eec.pen, piC.ghaf_pen.exp_unc) 
  par.eec.all[, 58] = parSapply(clus1, eec.mal, piC.tch92_mal_unc) 
  par.eec.all[, 59] = parSapply(clus1, eec.ch, piC_ch_Hash11_uncertainty) 
  par.eec.all[, 60] = parSapply(clus1, eec.prof, piC_pr_Hash11_uncertainty) 
  par.eec.all[, 61] = parSapply(clus1, eec.ch, piM_ch_Hash11_uncertainty) 
  par.eec.all[, 62] = parSapply(clus1, eec.prof, piM_pr_Hash11_uncertainty) 
  par.eec.all[, 63] = parSapply(clus1, eec.btr, piM.ghaf_butr.exp_unc) 
  par.eec.all[, 64] = parSapply(clus1, eec.gly, piM.ghaf_gly.exp_unc) 
  par.eec.all[, 65] = parSapply(clus1, eec.pen, piM.ghaf_pen.exp_unc) 
  par.eec.all[, 66] = parSapply(clus1, eec.mal, piM.tch91_mal_unc) 
  par.eec.all[, 67] = parSapply(clus1, eec.but, piM.tant02_but.exp_unc) 
  par.eec.all[, 68] = parSapply(clus1, eec.fpb, piM.tant02_fpb.exp_unc) 
  par.eec.all[, 69] = parSapply(clus1, eec.amm, piM.tch91_amm_unc) 
  par.eec.all[, 70] = parSapply(clus1, eec.urea, piM.tch91_ure_unc) 
  par.eec.all[, 71] = parSapply(clus1, eec.amm, tch91.egv.amm_unc) 
  par.eec.all[, 72] = parSapply(clus1, eec.urea, tch91.egv.ure_unc) 
#0.1 EEC
  par.0.1eec.all[, 49] = parSapply(clus1, 0.1*eec.atr, piC.grg08_atr_unc2) 
  par.0.1eec.all[, 50] = parSapply(clus1, 0.1*eec.atr, piC.kop_atr_unc) 
  par.0.1eec.all[, 51] = parSapply(clus1, 0.1*eec.atr, piC.atr.rohr08.lin) 
  par.0.1eec.all[, 52] = parSapply(clus1, 0.1*eec.atr, piC.meta_atr_unc) 
  par.0.1eec.all[, 53] = parSapply(clus1, 0.1*eec.but, piC.tant02_but.exp_unc) 
  par.0.1eec.all[, 54] = parSapply(clus1, 0.1*eec.fpb, piC.tant02_fpb.exp_unc) 
  par.0.1eec.all[, 55] = parSapply(clus1, 0.1*eec.btr, piC.ghaf_butr.exp_unc) 
  par.0.1eec.all[, 56] = parSapply(clus1, 0.1*eec.gly, piC.ghaf_gly.exp_unc) 
  par.0.1eec.all[, 57] = parSapply(clus1, 0.1*eec.pen, piC.ghaf_pen.exp_unc) 
  par.0.1eec.all[, 58] = parSapply(clus1, 0.1*eec.mal, piC.tch92_mal_unc) 
  par.0.1eec.all[, 59] = parSapply(clus1, 0.1*eec.ch, piC_ch_Hash11_uncertainty) 
  par.0.1eec.all[, 60] = parSapply(clus1, 0.1*eec.prof, piC_pr_Hash11_uncertainty) 
  par.0.1eec.all[, 61] = parSapply(clus1, 0.1*eec.ch, piM_ch_Hash11_uncertainty) 
  par.0.1eec.all[, 62] = parSapply(clus1, 0.1*eec.prof, piM_pr_Hash11_uncertainty) 
  par.0.1eec.all[, 63] = parSapply(clus1, 0.1*eec.btr, piM.ghaf_butr.exp_unc) 
  par.0.1eec.all[, 64] = parSapply(clus1, 0.1*eec.gly, piM.ghaf_gly.exp_unc) 
  par.0.1eec.all[, 65] = parSapply(clus1, 0.1*eec.pen, piM.ghaf_pen.exp_unc) 
  par.0.1eec.all[, 66] = parSapply(clus1, 0.1*eec.mal, piM.tch91_mal_unc) 
  par.0.1eec.all[, 67] = parSapply(clus1, 0.1*eec.but, piM.tant02_but.exp_unc) 
  par.0.1eec.all[, 68] = parSapply(clus1, 0.1*eec.fpb, piM.tant02_fpb.exp_unc) 
  par.0.1eec.all[, 69] = parSapply(clus1, 0.1*eec.amm, piM.tch91_amm_unc) 
  par.0.1eec.all[, 70] = parSapply(clus1, 0.1*eec.urea, piM.tch91_ure_unc) 
  par.0.1eec.all[, 71] = parSapply(clus1, 0.1*eec.amm, tch91.egv.amm_unc) 
  par.0.1eec.all[, 72] = parSapply(clus1, 0.1*eec.urea, tch91.egv.ure_unc) 
  
stopCluster(clus1) 

save(par.eec.all, file = "Review_models/r0_EECs/eec_par_values_2-19-18.RData")  
save(par.0.1eec.all, file = "Review_models/r0_EECs/eec10percent_par_values_2-19-18.RData")  