#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Increased cercarial shedding rate and snail reproduction from Johnson et al PNAS 2007
#Snail species is Planorbella trivolvis & trematode is Ribeiroa ondatrae

cerc = read.csv('Agrochemical_Review/Response_Fxs/Data/Johnson2007_Cerc_shed.csv')

cerc_ctrl <- cerc$c_hr[1]
cerc_ctrl_se <- cerc$c_hr_se[1]

cerc_fert <- cerc$c_hr[2]
cerc_fert_se <- cerc$c_hr_se[2]

johnson07_theta_uncertainty = function(){
  theta_q = rnorm(1, cerc_fert, cerc_fert_se) / rnorm(1, cerc_ctrl, cerc_ctrl_se)
  
  theta_q
} #cercariae per infected snail per hour in high nutrient group divided by control group gives estimate of rel. increase

eggs = read.csv('Agrochemical_Review/Response_Fxs/Data/Johnson2007_snail_eggs.csv')

eggs_ctrl <- eggs$eggs_lo
eggs_ctrl_se <- eggs$eggs_lo_se

eggs_fert <- eggs$eggs_hi
eggs_fert_se <- eggs$eggs_hi_se

johnson07_fNq_uncertainty = function(){

  fn = sum(rnorm(1, eggs_fert[1], eggs_fert_se[1]),
           rnorm(1, eggs_fert[2], eggs_fert_se[2]),
           rnorm(1, eggs_fert[3], eggs_fert_se[3]),
           rnorm(1, eggs_fert[4], eggs_fert_se[4])) /
        sum(rnorm(1, eggs_ctrl[1], eggs_ctrl_se[1]),
            rnorm(1, eggs_ctrl[2], eggs_ctrl_se[2]),
            rnorm(1, eggs_ctrl[3], eggs_ctrl_se[3]),
            rnorm(1, eggs_ctrl[4], eggs_ctrl_se[4]))
  
  while(fn <0) {
    fn = sum(rnorm(1, eggs_fert[1], eggs_fert_se[1]),
           rnorm(1, eggs_fert[2], eggs_fert_se[2]),
           rnorm(1, eggs_fert[3], eggs_fert_se[3]),
           rnorm(1, eggs_fert[4], eggs_fert_se[4])) /
        sum(rnorm(1, eggs_ctrl[1], eggs_ctrl_se[1]),
            rnorm(1, eggs_ctrl[2], eggs_ctrl_se[2]),
            rnorm(1, eggs_ctrl[3], eggs_ctrl_se[3]),
            rnorm(1, eggs_ctrl[4], eggs_ctrl_se[4]))
  }
  return(fn)
}

bm = read.csv('Agrochemical_Review/Response_Fxs/Data/Johnson2007_snail_biomass.csv')

  bm_ctrl = bm$bm[bm$treat == 'ctrl' & bm$Day == max(bm$Day[bm$treat == 'fert'])]
  bm_ctrl_se = bm$bmse[bm$treat == 'ctrl' & bm$Day == max(bm$Day[bm$treat == 'fert'])]
  
  bm_fert = bm$bm[bm$treat == 'fert' & bm$Day == max(bm$Day[bm$treat == 'fert'])]
  bm_fert_se = bm$bmse[bm$treat == 'fert' & bm$Day == max(bm$Day[bm$treat == 'fert'])]
  
johnson07_phiNq_uncertainty = function(){

  phi_N = rnorm(1, bm_fert, bm_fert_se) /
          rnorm(1, bm_ctrl, bm_ctrl_se)
  
  while(phi_N < 0){
    phi_N = rnorm(1, bm_fert, bm_fert_se) /
            rnorm(1, bm_ctrl, bm_ctrl_se)
  }
  
  return(phi_N)
}

keep.johnson07 = c('johnson07_theta_uncertainty', 'cerc_ctrl', 'cerc_ctrl_se', 'cerc_fert', 'cerc_fert_se', 
                   'johnson07_fNq_uncertainty', 'eggs_ctrl', 'eggs_ctrl_se', 'eggs_fert', 'eggs_fert_se',
                   'johnson07_phiNq_uncertainty', 'bm_ctrl', 'bm_ctrl_se', 'bm_fert', 'bm_fert_se')