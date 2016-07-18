#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#function used to model halstead et al mesocosm results impact on human infection

require(deSolve)
require(ggplot2)
require(reshape)
require(lattice)
require(ggthemes)
require(gridExtra)
st.er <- function(x) {
  sd(x)/sqrt(length(x))
} #Function to calculate standard error of the mean
cbPalette <- c("#999999", "yellow", "red", "green", "orange", "darkgreen", "pink", "blue")

#Model structure and equations ####################
#Represents a simplified version only including agroC effects found significant in mesocosms
schisto_halstead=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    P=n[5]
    N=S+E+I
    
    #Miracidia production; function of adult female worms alive in the system
    M=0.5*W*H*m 
    #egg input in mesocosm = 2276 eggs/tank/week = 325/day ; V=8.4%
    
    #Rate of predation; TOOK OUT alpha_q to simplify because of high dose of ChlorP 
    pred= (alpha*P)/(1+(alpha*N*Th)) #death rate of snails due to predators (Prawns)
    
    #TOOK OUT f_Nq and mu_Nq because they were not detectable in mesocosm experiments
    dSdt= (f_N*f_Nq)*((1-((phi_N*(1/phi_Nq))*N))*(S+(z*E))) - 
      ((mu_N)*S) - (pred*S) - (beta*(M)*S) #Susceptible snails
    
    dEdt= beta*(M)*S - ((mu_N+pred+sigma)*E) #Exposed snails
    
    dIdt= (sigma*E) - ((mu_N+pred+mu_I)*I) #Infected snails
    
    dWdt= (lamda*I) - ((mu_W+mu_H)*W) #worm burden in human population
    
    dPdt= (f_P*(1-(phi_P*P))*P)-((mu_P+muPq)*P) #prawn population
    
    
    return(list(c(dSdt,dEdt,dIdt,dWdt, dPdt)))
  }) 
} 

#Parameter values #################
parameters=c( 
  ##snail parameters 
  f_N=0.16, # recruitment rate (from sokolow et al)
  phi_N=(1-1/(80*0.16))/10000, # carrying capacity from sokolow et al
  z=0.5, #Proportion of exposed snails that reproduce from sokolow et al
  mu_N=1/80, #Mortality rate from Sokolow et al
  sigma=1/50, #Transition rate from exposed to infected from sokolow et al
  mu_I=1/20, #additional snail death due to infection from sokolow et al
  f_Nq=1,
  phi_Nq=1,
  
  #prawn parameters
  alpha=0.003, #attack rate
  Th=0.1,#~Prawn predation limit
  f_P=0.234/2, #prawn birth rate from Cervantes-Santiago Aquaculture 2010 paper (/2 for 1:1 female-male ratio)
  phi_P=1/(40*3),  #prawn carrying capacity
  mu_P= 0.03809524, #observed daily 24 hr prawn mortality rate in chlorP-free tanks of mesocosm#0.00137, previous rate
  muPq=0,
  
  #Adult Worm, Miracidia and Circariae Parameters
  lamda=0.00004, #probability of snail shedding a cercariae that infects a human host and survives to reproduction
  mu_W=1/(3*365), # death rate of adult worms
  m=0.5, #miracidial shedding rate per reproductive female divided by miracidial mortality; from sokolow et al
  
  #Human parameters
  H=300, #number of humans
  mu_H=1/(60*365), #Assumes 60 year lifespan
  beta = 2.038e-06 #max-likelihood estimate from fit to Lampsar follow-up 1 - 3
)

#Set initial values #########
  yrs=30
  time=seq(0,365*yrs,1)
  p=1 
  nstart=c(S=7000,E=3750,I=1200, W=20, P=p)
  
#Agro-free model run############  
  output.free=as.data.frame(ode(nstart,time,schisto_halstead,parameters)) 
  
  eqbm.free<-output.free[365*yrs,]
  eqbm.free #Snail extirpation and R0 below 1
  
  #plot time series of state variables
  plot(output.free[,1], output.free[,2], type='l', xlab="time",ylab="System Variables", 
       ylim=c(0,max( output.free[,2]+output.free[,3]+output.free[,4] )), col='blue', lty=2, lwd=2)
    lines(output.free[,1],output.free[,3],col='orange', lty=2, lwd=2)
    lines(output.free[,1],output.free[,4],col='red', lty=2, lwd=2)
    lines(output.free[,1],output.free[,2]+output.free[,3]+output.free[,4],col=1, lty=1, lwd=2)
    lines(output.free[,1],output.free[,5],col='green', lty=2)
    lines(output.free[,1],output.free[,6],col='purple', lty=2, lwd=2)
#Top-down model run with muPq from mesocosm chlorP########### 
  parameters['muPq']=0.73
  output.ch=as.data.frame(ode(nstart,time,schisto_halstead,parameters)) 
  
  eqbm.ch<-output.ch[365*yrs,]
  eqbm.ch #Snail extirpation and R0 below 1
  
  #plot time series of state variables
  plot(output.ch[,1], output.ch[,2], type='l', xlab="time",ylab="System Variables", 
       ylim=c(0,max( output.ch[,2]+output.ch[,3]+output.ch[,4] )), col='blue', lty=2, lwd=2)
    lines(output.ch[,1],output.ch[,3],col='orange', lty=2, lwd=2)
    lines(output.ch[,1],output.ch[,4],col='red', lty=2, lwd=2)
    lines(output.ch[,1],output.ch[,2]+output.ch[,3]+output.ch[,4],col=1, lty=1, lwd=2)
    lines(output.ch[,1],output.ch[,5],col='green', lty=2)
    lines(output.ch[,1],output.ch[,6],col='purple', lty=2, lwd=2)
  
#bottom-up model run with phi_Nq from mesocosm atrazine########### 
  parameters['muPq']=0
  parameters['phi_Nq']=1.61
  
  output.at=as.data.frame(ode(nstart,time,schisto_halstead,parameters)) 
  
  eqbm.at<-output.at[365*yrs,]
  eqbm.at #Snail extirpation and R0 below 1
  
  #plot time series of state variables
  plot(output.at[,1], output.at[,2], type='l', xlab="time",ylab="System Variables", 
       ylim=c(0,max( output.ch[,2]+output.ch[,3]+output.ch[,4] )), col='blue', lty=2, lwd=2)
    lines(output.at[,1],output.at[,3],col='orange', lty=2, lwd=2)
    lines(output.at[,1],output.at[,4],col='red', lty=2, lwd=2)
    lines(output.at[,1],output.at[,2]+output.at[,3]+output.at[,4],col=1, lty=1, lwd=2)
    lines(output.at[,1],output.at[,5],col='green', lty=2)
    lines(output.at[,1],output.at[,6],col='purple', lty=2, lwd=2)
  
#Ch:At combo model run########### 
  parameters['muPq']=0.73
  parameters['phi_Nq']=1.61
  
  output.atch=as.data.frame(ode(nstart,time,schisto_halstead,parameters)) 
  
  eqbm.atch<-output.atch[365*yrs,]
  eqbm.atch #Snail extirpation and R0 below 1
  
  #plot time series of state variables
  plot(output.atch[,1], output.atch[,2], type='l', xlab="time",ylab="System Variables", 
       ylim=c(0,max( output.ch[,2]+output.ch[,3]+output.ch[,4] )), col='blue', lty=2, lwd=2)
    lines(output.atch[,1],output.atch[,3],col='orange', lty=2, lwd=2)
    lines(output.atch[,1],output.atch[,4],col='red', lty=2, lwd=2)
    lines(output.atch[,1],output.atch[,2]+output.atch[,3]+output.atch[,4],col=1, lty=1, lwd=2)
    lines(output.atch[,1],output.atch[,5],col='green', lty=2)
    lines(output.atch[,1],output.atch[,6],col='purple', lty=2, lwd=2)
  
  
