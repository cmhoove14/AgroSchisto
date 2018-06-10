#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
#Data extraction and model fitting to Tantawy 2002 data

source("Agrochemical_Review/Models/litchfield_wilcoxon_get_b1_from_slope.R")
#Snail toxicity ##########
#Butachlor ##########
snail.but = data.frame(conc = c(0,0.65, 1.5, 4.5,6.5,44)*1000, # suspected typo in paper @4.5 (reads 45)
                       mort = c(0,0, .10, .25, .50, .90),
                       dead = c(0,0, 1, 2.5, 5, 9),
                       surv = 0)

  lc50.tant.but = 6.5
  slp.tant.but = get_b1(2.1)
  #get standard error from reported 95% CIs of lc50
  se.lc50.tant.but = mean(c(log10(10.4 / lc50.tant.but), 
                            log10(lc50.tant.but/4.06))) / 1.96 #st. err of lc50 in ppm
  
  muNq.tant.but_uncertainty<-function(He){
    Heu = (He/1000)
    lc50 = 10^(rnorm(1, log10(lc50.tant.but), se.lc50.tant.but))
    mun = pnorm(slp.tant.but * log10(Heu/lc50)) 

    return(mun)
  }
  
keep.tant.but = c('muNq.tant.but_uncertainty', 'lc50.tant.but', 'se.lc50.tant.but', 'slp.tant.but')  

#Fluazifop-p-butyl  ##########     
snail.fpb = data.frame(conc = c(0,1.76, 4.5, 9,17.6,58)*1000, 
                       mort = c(0,0, .10, .25, .50, .90),
                       dead = c(0,0, 1, 2.5, 5, 9),
                       surv = 0)

  lc50.tant.fpb = 17.6
  slp.tant.fpb = get_b1(1.8)
  #get standard error from reported 95% CIs of lc50
  se.lc50.tant.fpb = mean(c(log10(26.4 / lc50.tant.fpb), log10(lc50.tant.fpb/11.73))) / 1.96 #st. err of lc50 in ppm
  
  muNq.tant.fpb_uncertainty<-function(He){
    Heu = (He/1000)
    lc50 = 10^(rnorm(1, log10(lc50.tant.fpb), se.lc50.tant.fpb))
    mun = pnorm(slp.tant.fpb * log10(Heu/lc50)) 

    return(mun)
  }
  
keep.tant.fpb = c('muNq.tant.fpb_uncertainty', 'fx.tant.fpb', 'lc50.tant.fpb', 
                  'se.lc50.tant.fpb', 'slp.tant.fpb') 

#keep vector #########  
keep.muN.tantawy = c('muNq.tant.but_uncertainty', 'lc50.tant.but', 
                     'se.lc50.tant.but', 'slp.tant.but',
                     'muNq.tant.fpb_uncertainty', 'lc50.tant.fpb', 
                     'se.lc50.tant.fpb', 'slp.tant.fpb')  