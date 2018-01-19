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

johnson07_theta_uncertainty = function(waste){
  wst = waste
  theta_q = rnorm(1, 102.95, 15.77) / rnorm(1, 45.96, 9.13)
  
  while(theta_q < 1) theta_q = rnorm(1, 102.95, 15.77) / rnorm(1, 45.96, 9.13)
  
  theta_q
} #cercariae per infected snail per hour in high nutrient group divided by control group gives estimate of rel. increase

eggs = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Johnson2007_snail_eggs.csv')

johnson07_fN_uncertainty = function(waste){
  wst = waste

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

bm = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Johnson2007_snail_biomass.csv')

  plot(bm$Day[bm$treat == 'ctrl'], bm$bm[bm$treat == 'ctrl'], cex = 1.2, xlab = 'day', ylab = 'biomass',
       xlim = c(0, max(bm$Day) + 2), ylim = c(8,16))
  
    for(i in 1:length(bm$Day[bm$treat == 'ctrl'])){
      segments(x0 = bm$Day[i], x1 = bm$Day[i], y0 = bm$bmup[i], y1 = bm$bmlo[i])
    }
  
  points(bm$Day[bm$treat == 'fert'], bm$bm[bm$treat == 'fert'], pch = 16, cex = 1.2)
  
    for(i in 1:length(bm$Day[bm$treat == 'fert'])){
      segments(x0 = bm$Day[bm$treat == 'fert'][i], x1 = bm$Day[bm$treat == 'fert'][i], 
               y0 = bm$bmup[bm$treat == 'fert'][i], y1 = bm$bmlo[bm$treat == 'fert'][i])
    }
  
johnson07_phin_uncertainty = function(waste){
  wst = waste
  
  phi_N = rnorm(1, bm$bm[bm$treat == 'fert' & bm$Day == max(bm$Day[bm$treat == 'fert'])],
                bm$bmse[bm$treat == 'fert' & bm$Day == max(bm$Day[bm$treat == 'fert'])]) /
          rnorm(1, bm$bm[bm$treat == 'ctrl' & bm$Day == max(bm$Day[bm$treat == 'ctrl'])],
                bm$bmse[bm$treat == 'ctrl' & bm$Day == max(bm$Day[bm$treat == 'ctrl'])])
  return(phi_N)
}

keep.johnson07 = c('johnson07_theta_uncertainty', 'eggs', 'johnson07_fN_uncertainty', 'bm', 'johnson07_phin_uncertainty')