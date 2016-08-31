#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Toxicity to snails from Ibrahim 1992 ###################
  snail.repro = data.frame(dose = c(0,125,250,500),
                           f_red = c(1, 3005/4225, 1749/4225, 1313/4225))
  
  
  plot(snail.repro$dose, snail.repro$f_red, ylim = c(0,1), pch = 16, 
       ylab = 'reduction in snail recruitment rate', xlab = 'ChlorP ppb')
  
  mod1<-nls(f_red ~ exp(-b*(dose+0.1)), data=snail.repro, start = list(b=0.1))
    summary(mod1)
    
  dose = c(0:500)
  y = exp(-0.0028479*(dose))
  
  lines(dose, y, lty=2, col='red')
  
#toxicity to prawns from Halstead 2015 ############################
  pred.mort<-data.frame('dose'=c(0,0.64,3.2,6.4,32,64),
                        'death'=c(0,0,0,0.6,1,1))
  
  plot(pred.mort$dose, pred.mort$death, pch = 16, 
       ylab = 'predator mortality', xlab = 'chlorP ppb')
  
  mod2<-nls(death ~ exp(b*exp(-dose)), data=pred.mort, start = list(b = -4))
    summary(mod2)
    
  y2 = exp(-3.074e+02 *exp(-dose)) 
  
  lines(dose, y2, lty=2, col='red')
  