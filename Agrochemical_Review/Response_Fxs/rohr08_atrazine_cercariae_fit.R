#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

atr = read.csv('Agrochemical_Review/Response_Fxs/Data/rohr08_atr.csv')

atr.lin = lm(surv ~ log_conc, weights = st_err^-1, data = atr)

rohr.atr.fx = function(He){
  heu = log(He+1)
  predict(atr.lin, newdata = data.frame(log_conc = heu), interval = 'confidence', level = 0.95)
}

#Final function estimates relative mortality to He=0 and is interpreted as pi_C parameter      
  piC.atr.rohr08.lin = function(He){
      ts = predict(atr.lin, newdata = data.frame(log_conc = log(He+1)), se.fit = TRUE)[1:2]
      piC = rnorm(1, ts$fit, ts$se.fit) / atr.lin$coefficients[1]
    
    return(piC)
  } 

#Vector of items to keep
  keep.atr.rohr08 = c('atr', 'atr.lin', 'piC.atr.rohr08.lin')