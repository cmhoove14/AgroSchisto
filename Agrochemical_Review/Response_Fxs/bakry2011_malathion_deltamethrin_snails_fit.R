#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Data extraction and model fitting to Bakry 2011 [Bakry et al 2011](https://www.sciencedirect.com/science/article/pii/S0048357511001283) data

#Snail (H. duryi) toxicity ##########
#Malathion reported LC50 and slope data #####################
  lc50.bak.mal.report = 1.760
  slp.bak.mal.report = 2.74
  
  #get standard error from reported 95% CIs of lc50
      se.lc50.bak.mal = mean(c(log10(3.12/lc50.bak.mal.report), log10(lc50.bak.mal.report/0.99))) / 1.96
  
#Create function based on reverse of litchfield and wilcoxon      
muNq_mal_Bakry11_uncertainty = function(In){
  Ins = (In/1000) #Parameters based on ppm, data input as ppb
  lc50 = 10^(rnorm(1, log10(lc50.bak.mal.report), se.lc50.bak.mal)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
    mun = pnorm((slp.bak.mal.report) * log10(Ins/lc50)) #Estimate daily mortality (percent)

    return(mun)
  }
  
#keep vector
  keep.bak.mal = c('muNq_mal_Bakry11_uncertainty', 'lc50.bak.mal.report', 'se.lc50.bak.mal', 'slp.bak.mal.report')    
  
#Deltamethrin reported LC50 and slope data #############    
  lc50.bak.del.report = 4.82
  slp.bak.del.report = 2.74
  
  #get standard error from reported 95% CIs of lc50
    se.lc50.bak.del = mean(c(log10(7.7/lc50.bak.del.report), log10(lc50.bak.del.report/3.1))) / 1.96
  
  muNq_del_Bakry11_uncertainty = function(In){
    Ins = (In/1000) #Parameters based on ppm, data input as ppb
    lc50 = 10^(rnorm(1, log10(lc50.bak.del.report), se.lc50.bak.mal)) #Estimate lc50 with uncertainty and backtransform from log10 scale
    
    mun = pnorm((slp.bak.del.report) * log10(Ins/lc50))
    
    return(mun)
  }
  
  #keep vector
  keep.bak.del = c('muNq_del_Bakry11_uncertainty', 'lc50.bak.del.report', 'se.lc50.bak.del', 'slp.bak.del.report')    

#Snail reproduction over time, not analyzed at this point because only control and single dose group for each chemical
bakry11<-read.csv('Agrochemical_Review/Response_Fxs/Data/bakry2011.csv')
  bakry11_ctrl <-subset(bakry11, chem == "control")
  bakry11_mal <-subset(bakry11, chem == "malathion")
  bakry11_del <-subset(bakry11, chem == "deltamethrin")
  