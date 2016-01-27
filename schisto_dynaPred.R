## A schisto model with dynamic predator population and agrochemical effects on predators
#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.##########

library(deSolve)
schisto_pred=function(t, n, parameters) { 
  with(as.list(parameters),{
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    P=n[5]
    N=S+E+I
    pred= (alpha*alpha_q*P)/(1+(alpha*alpha_q)*N*Th)
    #inf= 0.5*beta*M*H*m*H
    inf= 0.5*beta*W*H*m*S
    dSdt= f_N*((1-(phi_N*N))*(S+(z*E))) - (mu_N*S) - (pred*S) - inf
    dEdt= inf - ((mu_N+pred+sigma)*E)
    dIdt= (sigma*E) - ((mu_N+mu_I)*I) - (pred*I)
    dWdt= (lamda*I) - ((mu_W+mu_H)*W)
    dPdt= f_p*(1-(phi_p*P))*P-(mu_p + mu_p_q)*P
    return(list(c(dSdt,dEdt,dIdt,dWdt,dPdt)))
  }) 
} 

#Reductions in attack rate, alpha, from AgroC concentration, q #############################
#reduction in predation rate (Machrobrachium nipponense) caused by pesticide exposure taken from Yuan et al,2004
  #where q_red = fraction of the control predation rate maintained by predators exposed to sublethal pesticide doses

  #if Paraquat (herbicide)
    #dose=10 ug/L; q_red=0.60
    #dose=50 ug/L; q_red=0.34
    #dose=100 ug/L; q_red=0.17
  #if Malathion (insecticide)
    #dose=10 ug/L; q_red=0.84
    #dose=50 ug/L; q_red=0.69
    #dose=100 ug/L; q_red=0.56

alpha_q=1
  
#Increased mortality rate (Procambarus alleni) from exposure to a variety of insecticides;
#EECs and 96 hr LC50s taken from Halstead, Civitello & Rohr, 2015, table1, fig1
#In DAILY time step model, TOTAL mortality rate @these LC50s = (P*0.5)/3 
#i.e. 50% of the pop dead at 3 days (96hrs) averaged to daily rate over three days 
  #~-~-~PYRETHROIDS~-~-~#
  #if malathion (organophosphate insecticide)
    #LC50=50000 ug/L; mean EEC = ~6 ug/L; 0% of env samples within 95% CI of LC50
  #if Chlorpyrifos
    #LC50=29.3 ug/L; mean EEC = ~8 ug/L; 5% of env samples within 95% CI of LC50
  #if Terbufos
    #LC50=8.89 ug/L; mean EEC = ~2 ug/L; 22% of env samples within 95% CI of LC50  
  #~-~-~ORGANOPHOSPHATES~-~-~#
  #if Esfenvalerate
    #LC50=0.22 ug/L; mean EEC = ~0.4 ug/L; 97% of env samples within 95% CI of LC50 
  #if lambda-cyhalothrin
    #LC50=0.21 ug/L; mean EEC = ~8 ug/L; 100% of env samples within 95% CI of LC50  
  #if Permethrin
    #LC50=0.58 ug/L; mean EEC = ~2 ug/L; 99% of env samples within 95% CI of LC50  

estimate=0
#Parameters and starting values ###################################
yrs=30
time=seq(0,365*yrs,1)
p=1
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
             m=0.5, #Miracidial production per reproductive female
             f_p=0.128, #Prawn per capita daily mortality rate
             alpha_q=1, #Reduction in prawn attack rate from AgroC
             phi_p=0.02, #Predator prawn carrying capacity
             mu_p=0.00137, #Predator prawn natural per capita daily mortality rate
             mu_p_q=0 #Predator prawn enhanced morality rate from agrochemicals
)

nstart=c(S=2435,E=4445,I=1422, M=59.3, P=p) #Prawn-free equilibrium values
output=as.data.frame(ode(nstart,time,schisto_pred,parameters))
output[0:200,]
eqbm<-output[365*yrs,]
eqbm

#Plot changes over time period ##############################
par(mfrow=c(1,1))
plot(output[,1], output[,2], type='l', xlab="time",ylab="Number of Snails",
     ylim=c(0,max(output[,2]+output[10950,3]+output[10950,4])))
lines(output[,1],output[,3],col="red")
lines(output[,1],output[,4],col="blue")
lines(output[,1],output[,2]+output[,3]+output[,4],col="green")
lines(output[,1],output[,5],col="brown")
lines(output[,1],output[,6],col="purple")
legend("topright",paste(c("S","E","I","N","W","P"),
                        as.integer(c(eqbm[,2],eqbm[,3],eqbm[,4],(eqbm[,2]+eqbm[,3]+eqbm[,4]),eqbm[,5],eqbm[,6])),
                        sep=" - "),
       col=c("black","red","blue","green","brown","purple"),lty=1,cex=0.8)

#Plot showing parameters with lower magnitude ###############################
plot(output[,1], output[,5], type='l', col="brown", xlab="time",ylab="",
     main="", ylim=c(0,max(output[,5])))
lines(output[,1],output[,3],col="red")
lines(output[,1],output[,6],col="purple")
lines(output[,1],output[,4],col="blue")

#Validate that Phi_P is ~ P in Sokolow model ################
phi_P_range<-seq(from=1, to=100, by=1)
vals<-matrix(ncol=5, nrow=length(phi_P_range))
for(phi_P in phi_P_range){
  nstart=c(S=2435,E=4445,I=1422, M=59.3, P=p)
  parameters["phi_p"]<-phi_P^-1
  print(parameters["phi_p"])
  output=as.data.frame(ode(nstart,time,schisto_pred,parameters)) #params_withextra
  eqbm<-output[365*yrs,]
  eqbm
  print(eqbm)
  vals[phi_P,1]<-phi_P^-1
  vals[phi_P,2]<-as.numeric(eqbm[,2]+eqbm[,3]+eqbm[,4])
  vals[phi_P,3]<-as.numeric(eqbm[,4])
  vals[phi_P,4]<-as.numeric(eqbm[,5])
  vals[phi_P,5]<-as.numeric(prev(W=eqbm[,5]))
}

#Plot N and I over phi_p #############
par(mfrow=c(2,1), cex.lab=1.5, cex.axis=1.3)
#Plot total number of snails (N) across different values of P
par(mai=c(0.2, 1.1, 0.1,0.4))
plot(vals[,1]^-1, vals[,2], type='l', lwd=3, xlab="",ylab="",
     col="black", ylim=c(0,9000), xaxt='n')
  mtext(expression(italic(N)), side=2, cex=2, las=2, line=2.6)

#Plot total number of infected snails (I) across different values of P
par(mai=c(1.2, 1.1, 0.1,0.4))
plot(vals[,1]^-1, vals[,3], type='l', lwd=3, xlab="Prawn Carrying Capacity",ylab="", 
     col="blue", ylim=c(0,1500))
  mtext(expression(italic(I)), side=2, cex=2, las=2, line=2.6, col="blue")

#Plot showing only mean worm burden over time with the dynamic P function #################
plot(output[,1], output[,5], type='l', col="brown", xlab="time",ylab="Mean Worm Burden",
     ylim=c(0,max(output[,5])))
