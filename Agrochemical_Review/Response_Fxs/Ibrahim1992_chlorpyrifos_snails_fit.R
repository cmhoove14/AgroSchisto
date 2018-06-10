#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Toxicity to Biomphalaria snails from Ibrahim 1992 ###################
require(drc)

# reduction in snail fecundity: from table 1, relative change in number of juveniles produced over the 8 week experiment ####
snail.repro = data.frame(dose = c(0,125,250,500), #chlorpyrifos in ppb
                         juvs = c(4225, 3005, 1749, 1313),
                        #convert 12 week juvenile cohort to juveniles per snail per day
                         juvs.sn.day = c(4225, 3005, 1749, 1313) / 30 / (7*12)) 

#model with juveniles/snail/day as outcome  
chlor.fN.predict = drm(juvs.sn.day ~ dose, data = snail.repro, type = 'continuous',
                       fct = LL.3(names = c("b", "d", "e"),
                                  fixed = c(NA, max(snail.repro$juvs.sn.day), NA)))

  fNq_chlor_ibr92_uncertainty<-function(In){
      init = predict(chlor.fN.predict, newdata = data.frame(mal = In), se.fit = T)
      fn = rnorm(1, init[1], init[2]) / snail.repro$juvs.sn.day[1]

    return(fn)
  }
  
keep.ibr.fn.ch = c('snail.repro', 'chlor.fN.predict', 'fNq_chlor_ibr92_uncertainty')  

#Snail mortality #############
  #-- relative change in survival over entire experiment period (12 weeks)
  #Not advised to use this as mortality estimate as outcomes were only assessed at 5 weeks.
      #though authors note that all snails died within the first week @ 500ppb
  snail.mort = data.frame(dose = c(0,125,250,500),
                          total = rep(30,4),
                          dead = c(5,7,8,30),
                          live = c(25,23,22,0))
  
  snail.mort$mort = snail.mort$dead / snail.mort$total
  
  ibr_muNq<-drm(dead/total ~ dose, weights = total, data = snail.mort,
                type = 'binomial', fct = LL.2(names = c('b', 'e'),
                                              fixed = c(NA, NA)))
  
  muNq_chlor_ibr92_uncertainty<-function(In){
      init = predict(ibr_muNq, data.frame(dose = In), se.fit = T)
      mun = rnorm(1, init[1], init[2])
   
    return(mun)
  }
  
keep.ibr.ch = c(keep.ibr.fn.ch, 'muNq_chlor_ibr92_uncertainty', 'ibr_muNq', 'snail.mort')      