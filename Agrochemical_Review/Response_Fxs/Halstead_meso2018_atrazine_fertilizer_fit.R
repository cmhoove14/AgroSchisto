#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Step 1.1 Estimate change in snail carrying capacity by relative increases in final snail numbers
meso_dat_botup<-read.csv('Agrochemical_Review/Response_Fxs/Data/Halstead_meso_data.csv')

    Fes<-meso_dat_botup$Snails.Live[meso_dat_botup$Ch==1 & 
                                 meso_dat_botup$At==0 &
                                 meso_dat_botup$Fe==1]
      fe.mean = mean(Fes)
      fe.sd = sd(Fes)
      
    Ats<-meso_dat_botup$Snails.Live[meso_dat_botup$Ch==1 & 
                                 meso_dat_botup$At==1 &
                                 meso_dat_botup$Fe==0]
      at.mean = mean(Ats)
      at.sd = sd(Ats)

    refs<-meso_dat_botup$Snails.Live[meso_dat_botup$Ch==1 & 
                                 meso_dat_botup$At==0 &
                                 meso_dat_botup$Fe==0]
      ref.mean<-mean(refs)
      ref.sd<-sd(refs)
    
#Fertilizer bottom up effects  
  halstead18_phiNq_fe_uncertainty = function(...){ #function based on snail hatchlings
    sample(Fes, 1) / sample(refs,1)
    
    #phin = rnorm(1, fe.mean, fe.sd) / ref.mean
    #while(phin <= 0) phin = rnorm(1, fe.mean, fe.sd) / ref.mean
    #phin
  }
      
      
#Atrazine bottom up effects  
  halstead18_phiNq_at_uncertainty = function(...){ #function based on snail hatchlings
        sample(Ats, 1) / sample(refs,1)

    #phin = rnorm(1, at.mean, at.sd) / ref.mean
    #while(phin <= 0) phin = rnorm(1, at.mean, at.sd) / ref.mean
    #phin
  }

keep.halstead18 = c('halstead18_phiNq_at_uncertainty', 'at.mean', 'at.sd', 'ref.mean',
                    'halstead18_phiNq_fe_uncertainty', 'fe.mean', 'fe.sd')    
  
#mean(sapply(runif(5000, 0, 1000), halstead17_phiN_fe_uncertainty))
#mean(sapply(runif(5000, 0, 1000), halstead17_phiN_at_uncertainty))
  