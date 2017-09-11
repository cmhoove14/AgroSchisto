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

cerc = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Johnson2007_Cerc_shed.csv')

johnson07_theta_uncertainty = function(seed){
  set.seed(seed)
  theta_q = rnorm(1, 102.95, 15.77) / rnorm(1, 45.96, 9.13)
  theta_q
} #cercariae per infected snail per hour in high nutrient group divided by control group gives estimate of rel. increase

eggs = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Johnson2007_snail_eggs.csv')

johnson07_fN_uncertainty = function(seed){
  set.seed(seed)
  
  fN = sum(rnorm(1, eggs$eggs_hi[1], eggs$eggs_hi_se[1]),
           rnorm(1, eggs$eggs_hi[2], eggs$eggs_hi_se[2]),
           rnorm(1, eggs$eggs_hi[3], eggs$eggs_hi_se[3]),
           rnorm(1, eggs$eggs_hi[4], eggs$eggs_hi_se[4])) /
        sum(rnorm(1, eggs$eggs_lo[1], eggs$eggs_lo_se[1]),
            rnorm(1, eggs$eggs_lo[2], eggs$eggs_lo_se[2]),
            rnorm(1, eggs$eggs_lo[3], eggs$eggs_lo_se[3]),
            rnorm(1, eggs$eggs_lo[4], eggs$eggs_lo_se[4]))
  fN
}

keep.johnson07 = c('johnson07_theta_uncertainty', 'eggs', 'johnson07_fN_uncertainty')