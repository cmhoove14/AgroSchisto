#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(deSolve)
source('Review_models/agroReview_mod1.1.R')
source('Review_models/Response_fxs.R')

#Set baseline ######

time = seq(0,365*5,1)
nstart3['He']<-0
base<-as.data.frame(ode(nstart3, time, mod1, parameters))
      base$N = (base$S + base$E + base$I)
      deriv = numeric(length(base$W))
      for (i in 0:(length(base$W)-1)) {
        deriv[i] = base$W[i] - base$W[i+1]
      }
      base$dW = -deriv
      
  plot(base$time, base$dW, type = 'l', lwd = 2, ylim = c(-0.1,0.1))    

#Add predators  
In_fx<-function(In){
  Satapornvanit2009.ch(In)
} #Set agrochemical response function to Chlorpyrifos influence in Satapornvanit paper
parameters['k_In'] = -log(0.5)/25.5     

       Ch.event<-data.frame(var = 'In',
                            time = 365,
                            value = 1,
                            method = 'add')  
      
base_plsCh<-as.data.frame(ode(nstart3, time, mod1, parameters, events = list(data = Ch.event)))
base_plsCh$N = (base_plsCh$S + base_plsCh$E + base_plsCh$I)
              deriv = numeric(length(base_plsCh$W))
              for (i in 0:(length(base_plsCh$W)-1)) {
                deriv[i] = base_plsCh$W[i] - base_plsCh$W[i+1]
              }
              base_plsCh$dW = -deriv  
              
    lines(base_plsCh$time[-length(base_plsCh$time)], base_plsCh$dW[-length(base_plsCh$time)], lty = 2)          
        

