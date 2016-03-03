#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

library(deSolve)
library(splines)
snails=function(t, n, parameters) { #Snail population model with birth and death rates and carrying capacity
  with(as.list(parameters),{
    N=n[1]
    dNdt= f_N*((1-(phi_N*N))*N) - (mu_N*N)
    return(list(c(dNdt)))
  })
} 

parameters=c(f_N=0.2,
             phi_N=1/250,
             mu_N=1/80)
time=c(1:(12*7)) #12 week mesocosm experiment
nstart=c(N=11) #B. truncatus initial population  

output.2=as.data.frame(ode(nstart,time,snails,parameters))

parameters['f_N']<-0.16 #Value used in Sokolow et al 2015
output.16=as.data.frame(ode(nstart,time,snails,parameters))

parameters['f_N']<-0.035 #.03 is about what was estimated from mesocosm experiments assuming exponential decay throughout
parameters['phi_N']<-0
parameters['mu_N']<-0
output.03=as.data.frame(ode(nstart,time,snails,parameters))

parameters=c(f_N=0.2,
             phi_N=1/250,
             mu_N=1/80)


  
#Find and plot derivatives (i.e. dN/dt over time) 
  dY.2<-diff(output.2[,2]/diff(output.2[,1]))
  dX.2<-rowMeans(embed(output.2[,1],2))
    
  dY.16<-diff(output.16[,2]/diff(output.16[,1]))
  dX.16<-rowMeans(embed(output.16[,1],2))
  
  dY.03<-diff(output.03[,2]/diff(output.03[,1]))
  dX.03<-rowMeans(embed(output.03[,1],2))
 
#Plot what theoretical population curves would look like given above parameters     
  par(mfrow=c(2,1), mar=c(4,4,0.4,2))
  plot(x=output.2[,1], y=output.2[,2], type='l', ylim=c(0,250), xlim=c(0, max(output.2[,1])),
       xlab='', ylab='N', col='red', lwd=2)
    lines(x=output.16[,1], y=output.16[,2], col='blue', lwd=2)
    lines(x=output.03[,1], y=output.03[,2], col='green', lwd=2)
    legend("bottomright", legend=c('f_N - 0.2','f_N - 0.16', 'f_N - 0.035'), col=c('red', 'blue', 'green'), lty=1, lwd=1.5)
    legend("topleft", legend=c(paste('Death - ', parameters['mu_N']), 
                                   paste('Capacity - ', 1/parameters['phi_N'])))
  plot(x=dX.2, y=dY.2, type='l', xlab='time (days)', ylab='dN/dt', lty=2, col='red', lwd=2)
    lines(x=dX.16, y=dY.16, lty=2, col='blue', lwd=2)
    lines(x=dX.03, y=dY.03, lty=2, col='green', lwd=2)

 

  
snailt<-data.frame('time'=c(1:(12*7)),'N'=rep(0, times=(12*7)))
  for(i in 1:(12*7)){
    snailt[i,2]=(11/250)/(11+((1/250)-11)*exp(0.16*snailt[i,1]))
  }
