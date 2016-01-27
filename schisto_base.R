## A schisto model
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
schisto=function(t, n, parameters) { 
  with(as.list(parameters),{
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    N=S+E+I
    pred= (alpha*P)/(1+alpha*N*Th)
    #inf= 0.5*beta*W*H*m*H
    inf= 0.5*beta*W*H*m*S
    dSdt= f_N*((1-(phi_N*N))*(S+(z*E))) - (mu_N*S) - (pred*S) - inf
    dEdt= inf - ((mu_N+pred+sigma)*E)
    dIdt= (sigma*E) - ((mu_N+mu_I)*I) - (pred*I)
    dWdt= (lamda*I) - ((mu_W+mu_H)*W)
    return(list(c(dSdt,dEdt,dIdt,dWdt)))
  }) 
} 


#Parameters from Sokolow et al ########################################
yrs=30
time=seq(0,365*yrs,1)
p=0
parameters=c(f_N=0.16, #Instantaneous snail mortality rate
             phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
             z=0.5, #Fraction of exposed snails that reproduce
             mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
             alpha=0.003, #per capita prawn attack rate on snails
             P=p, #Number of prawns in the system
             Th=0.1, #Prawn predation max
             beta=0.000004, #Snail infection probability
             H=1000, # human population
             sigma=0.02, # exposed snails becoming infected
             mu_I=0.05, #Enhanced mortality of infected snails
             lamda=0.00004, #Human infection probability from infected snails
             mu_H=0.00004566, #Worm mortality from host deaths
             mu_W=0.00091324, #Adult worm natural mortality
             m=0.5 #Miracidial production per reproductive female
             )
#Parameters we might consider changing based on other sources
    #beta = 0.000001-0.0005 from Liang et al 2002
    #f = 0.5/7 = 0.07/day snail recruitment rate from Woolhouse, 1990 fig4 @ 21 C

# Various starting values ########################
#nstart=c(S=5500,E=1500,I=1000, W=50) #VALUES USED HERE ARE Randomly contrived

nstart=c(S=2435,E=4444,I=1422, W=59) #VALUES USED HERE ARE PRAWN-FREE EQUILIBRIUM

#nstart=c(S=6097,E=346,I=45, W=26) #VALUES USED HERE ARE PRAWN=33 EQBM

#Run the model and print ~eqbm values #######################
output=as.data.frame(ode(nstart,time,schisto,parameters))
output[100:500,]
eqbm<-output[365*yrs,]
  eqbm

#Plot changes over time period
par(mfrow=c(1,1))
plot(output[,1], output[,2], type='l', xlab="time (days)",ylab="",
     ylim=c(0,max(output[,2]+output[,3]+output[,4])))
lines(output[,1],output[,3],col="red")
lines(output[,1],output[,4],col="blue")
lines(output[,1],output[,2]+output[,3]+output[,4],col="green")
lines(output[,1],output[,5],col="brown")
lines(output[,1],rep(p,times=length(output[,1])),col="purple")
legend("topright",paste(c("S","E","I","N","W","P"),
                        as.integer(c(eqbm[,2],eqbm[,3],eqbm[,4],(eqbm[,2]+eqbm[,3]+eqbm[,4]),eqbm[,5],p)),
                        sep=" - "),
       col=c("black","red","blue","green","brown","purple"),lty=1,cex=0.7)

#Plot only mean worm burden (W) over 30-yr run
plot(output[,1], output[,5], type='l', col="brown",xlab="time",ylab="Mean worm burden",
     ylim=c(0,max(output[,5])))

#Run the model with P=50 and observe dynamics over time  (Sok etal Fig5B&E) #####################
yrs=20
time=seq(0,365*yrs,1)
p=50

output=as.data.frame(ode(nstart,time,schisto,parameters))
output[100:500,]
eqbm2<-output[365*yrs,]
eqbm2

#Plots
par(mfrow=c(2,1), cex.lab=1.3, cex.axis=1.4)
par(mai=c(0.2, 1.1, 0.1,0.4))
  plot(output[,1]/365, output[,5], type='l', xlab="",ylab="", lwd=3,
        ylim=c(0,80), col="brown", xaxt='n')
    mtext("Mean worm burden (W)", side=2, line=2.3, cex=1.1)
par(mai=c(1.2, 1.1, 0.1,0.4)) 
  plot(output[,1]/365, prev(W=output[,5])*100, type='l', xlab="years",ylab="", lwd=3,
         ylim=c(0,100), col="red")
    mtext("Prevalence", side=2, line=2.3, cex=1.3)

#Function to calculate prevalence ####################################
prev<-function(W,k=0.25){
  1-(1+(W/k))^(-k)
}
prev2<-function(x){
  1-(1+(x/0.25))^(-0.25)
} #Same function, just with k already set to value cited in sokolow et al 2015
prev(W=eqbm[,5])

#Function to calculate R-0 of system give parameter values ################
r0<-function(f_N=0.16, #Instantaneous snail mortality rate
             phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
             z=0.5, #Fraction of exposed snails that reproduce
             mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
             alpha=0.003, #per capita prawn attack rate on snails
             P=p, #Number of prawns in the system
             Th=0.1, #Prawn predation max
             beta=0.000004, #Snail infection probability
             H=1000, # human population
             sigma=0.02, # exposed snails becoming infected
             mu_I=0.05, #Enhanced mortality of infected snails
             lamda=0.00004, #Human infection probability from infected snails
             mu_H=0.00004566, #Worm mortality from host deaths
             mu_W=0.00091324, #Adult worm natural mortality
             m=0.5 #Miracidial production per reproductive female
                )
  {
  pred=(alpha*P)/(1+(alpha*Th))
  print(pred)
  #top=((f+f*z)*(beta*m)*sigma*lamda)
  top=beta*0.5*H*m*(1/phi_N)*sigma*lamda
  print(top)
  bot=(mu_N+mu_I+pred)*(mu_N+sigma+pred)*(mu_H+mu_W)
  print(bot)
  r0=top/bot
  r0
}
r0()

#Function to observe instantaneous rates of change at given parameter values #############
d_all<-function(f_N=0.16, #Instantaneous snail mortality rate
                phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
                z=0.5, #Fraction of exposed snails that reproduce
                mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
                alpha=0.003, #per capita prawn attack rate on snails
                P=p, #Number of prawns in the system
                Th=0.1, #Prawn predation max
                beta=0.000004, #Snail infection probability
                H=1000, # human population
                sigma=0.02, # exposed snails becoming infected
                mu_I=0.05, #Enhanced mortality of infected snails
                lamda=0.00004, #Human infection probability from infected snails
                mu_H=0.00004566, #Worm mortality from host deaths
                mu_W=0.00091324, #Adult worm natural mortality
                m=0.5, #Miracidial production per reproductive female
                S, E, I, W){
  N=S+E+I
  pred= (alpha*P)/(1+(alpha*N*Th))
  inf= 0.5*beta*W*H*m*H
    dSdt= f_N*((1-(phi_N*N))*(S+(z*E))) - (mu_N*S) - (pred*S) - inf
  print(dSdt)
    dEdt= inf - ((mu_N+pred+sigma)*E)
  print(dEdt)
    dIdt= (sigma*E) - ((mu_N+mu_I)*I) - (pred*I)
  print(dIdt)
    dWdt= (lamda*I) - ((mu_W+mu_H)*W)
  print(dWdt)
}

#Check changes at equilibrium
d_all(S=eqbm[,2], E=eqbm[,3], I=eqbm[,4], W=eqbm[,5])

#Plot change in S,E,I,W over variable P, with runs equal to 30 years ##############
p_range<-seq(from=1, to=100, by=1)
  op<-rep(0, time=length(p_range))
  op2<-rep(0, time=length(p_range))
  op3<-rep(0, time=length(p_range))
  op4<-rep(0, time=length(p_range))
  for(p in p_range){
    ind<-which(p_range==p)
    time=seq(0,365*30,1)
    parameters=c(f_N=0.16, #Instantaneous snail mortality rate
                 phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
                 z=0.5, #Fraction of exposed snails that reproduce
                 mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
                 alpha=0.003, #per capita prawn attack rate on snails
                 P=p, #Number of prawns in the system
                 Th=0.1, #Prawn predation max
                 beta=0.000004, #Snail infection probability
                 H=1000, # human population
                 sigma=0.02, # exposed snails becoming infected
                 mu_I=0.05, #Enhanced mortality of infected snails
                 lamda=0.00004, #Human infection probability from infected snails
                 mu_H=0.00004566, #Worm mortality from host deaths
                 mu_W=0.00091324, #Adult worm natural mortality
                 m=0.5 #Miracidial production per reproductive female
    )
    nstart=c(S=4382,E=3847,I=684, W=28)  
    output=as.data.frame(ode(nstart,time,schisto,parameters)) 
    #op[ind]<-output[10950,2]+output[10950,3]+output[10950,4]
    op[ind]<-output[365*30,2]
    op2[ind]<-output[365*30,3]
    op3[ind]<-output[365*30,4]
    op4[ind]<-output[365*30,5]
  }
  #Plot results over Predator-values
  par(mfrow=c(1,1))
  plot(p_range, op/(op+op2+op3), type='l', xlab="",ylab="%Susceptible of all snails", col="black", lty=3)
  par(mfrow=c(2,2))
  plot(p_range, op, type='l', xlab="#Predators",ylab="S", col="black", ylim=c(0,max(op, na.rm=TRUE)))
  plot(p_range, op2, type="l", xlab="#Predators",ylab="E", col="red", ylim=c(0,max(op2, na.rm=TRUE)))
  plot(p_range, op3, type="l", xlab="#Predators",ylab="I", col="blue", ylim=c(0,max(op3, na.rm=TRUE)))
  plot(p_range, op4, type="l", xlab="#Predators",ylab="W", col="brown", ylim=c(0,max(op4, na.rm=TRUE)))
  
#Plot N & I over P for document#####################

par(mfrow=c(2,1), cex.lab=1.5, cex.axis=1.4)
#Plot total number of snails (N) across different values of P
par(mai=c(0.2, 1.1, 0.1,0.4))
plot(p_range, (op+op2+op3), type='l', lwd=3, xlab="",ylab="",
     col="black", xlim=c(0,100), ylim=c(0,9000), xaxt='n')
  mtext(expression(italic(N)), side=2, cex=2, las=2, line=2.6)

#Plot total number of infected snails (I) across different values of P
  par(mai=c(1.2, 1.1, 0.1,0.4))
  plot(p_range, (op3), type='l', lwd=3, xlab="Predator abundance",ylab="", 
     col="blue", xlim=c(0,100), ylim=c(0,1500))
  mtext(expression(italic(I)), side=2, cex=2, las=2, line=2.6, col="blue")

#Plot W & Prev over P for document#####################
par(mfrow=c(2,1), cex.lab=1.4, cex.axis=1.3)
  #W plot
  par(mai=c(0.2, 1.1, 0.1,0.4))
  plot(p_range, op4, type="l", ylab="", col="brown", 
       lwd=2, ylim=c(0,60), xaxt='n')
  mtext("Mean worm burden (W)", col="brown", line=2.5, cex=1.4, side=2)
  #Prev Plot
  par(mai=c(1.2, 1.1, 0.1,0.4))
  plot(p_range, prev(W=op4), type="l", xlab="Predator abundance",ylab="", 
       col="red", lwd=2, ylim=c(0,0.80))
    mtext("Prevalence", col="red", line=2.5, cex=1.5, side=2)
#########################################  
  #Plot proportion of snails that are susceptible at 30 yr equilibirum across different values of P
plot(p_range, ((op+op2)/(op+op2+op3))*100, type='l', xlab="#Predators",ylab="%NON-infected of all snails", 
     col="black", xlim=c(0,50), ylim=c(0,100))
abline(h=100, col="pink",lty=2)
