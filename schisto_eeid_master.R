#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Load libraries ##########
require(deSolve)
require(ggplot2)
require(reshape)

#Model structure and equations ####################

schisto_master=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    P=n[5]
    N=S+E+I
    
    #Miracidia production; function of adult female worms alive in the system (W*H*0.5) assuming 1:1 sex ratio,
    #mating probability function (gamma) from Anderson and May,
    
    fx<-function(x, mean.worm = W, clump = k){
      (1-cos(x))/((1+((mean.worm/(clump+mean.worm))*cos(x)))^1+clump)
    }
    gamma = 1 - ((1-(W/(W+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value)
    
    #per-reproductive female egg production (m),
    #egg viability or the fraction of eggs that successfully hatch into viable miracidia (v),
    
      M=((0.5*W*H)*gamma)#*m*u_H*(v*vq)
      
    #miracidial mortality and infectivity (perhaps influenced by agrochemicals) affects beta
    
    pred= (alpha*P)/(1+(alpha*((N/200)^nn)*Th)) #per capita predation of predators on snails; a Holling's type III functional response
    
    dSdt= f_N*(1-(N/(phi_N*phi_Nq)))*(S+E) - 
      mu_N*S - pred*((S/200)^nn) - beta*M*S #Susceptible snails
    
    dEdt= beta*M*S - mu_N*E - pred*((E/200)^nn) - sigma*E #Exposed snails
    
    dIdt= sigma*E - (mu_N+mu_I)*I - pred*((I/200)^nn) #Infected snails
    
    dWdt= lamda*I - (mu_W+mu_H)*W #worm burden in human population
    
    dPdt= f_P*(1-(P/phi_P))*P-(mu_P+muPq)*P #prawn population
    
    
    return(list(c(dSdt,dEdt,dIdt,dWdt, dPdt)))
  }) 
} 

#List parameters and values #####################
parameters=c( #Updated as of 4/19 to get values directly from cited sources instead of drawing directly from PNAS 
  ##standard snail parameters 
    f_N=0.1, # recruitment rate: snails/snail/day
    phi_N=10000, # carrying capacity: max snail population (corresponds to ~50/m^2)
    z=0.5, #Proportion of exposed snails that reproduce: density dependent, but assumed constant here
    mu_N=1/60, #Mortality rate from Anderson and May (from Chu 1966): deaths/snail/day
    sigma=1/40, #Transition rate from exposed to infected; ~latent period
    mu_I=1/10 - 1/60, #additional snail death due to infection
    ## snail parameters impacted by agrochemicals
    #f_Nq=1, #Not affected in mesocosm
    phi_Nq=1, #Scalar of snail carrying capacity by chemical concentration INFORMED BY BOTTOM UP EFFECTS IN MESOCOSM
    #mu_Nq=0, #Chem concentrations too low to affect snails in mesocosom experiments
    
  
  #prawn parameters
    alpha=0.003, #attack rate
    Th=0.067,#~Prawn predation limit for procambarius clarkii: ~15 snails/day
    f_P=0.117,#prawn birth rate: prawns/prawn/day
    phi_P=120,  #prawn carrying capacity ~0.6/m^2
    mu_P= 0.038095238, #observed daily 24 hr prawn mortality rate in chlorP-free tanks of mesocosm#0.00137, previous rate
    nn=1, #exponent of the prey density in the holling's type III functional response  
  
  #prawn parameter impacted by agrochemicals
    muPq=0, #agrochemical dependent prawn mortality rate (added to baseline mortality rate) INFORMED BY MESOCOSM AND HALSTEAD 2015
    #alpha_q=1, #Scalar of predation rate due to sub-lethal toxicity; not considered in mesocosm
    
  #Worm parameters
    mu_W=1/(3.3*365), # death rate of adult worms
    m=.36, #eggs/mL urine/female adult worm
    v=0.084, #Egg viability controlling schistosome egg->infective miracidia
    vq=1, #agrochemical-caused reduction in egg viability
    
  #Human parameters
    H=300, #number of humans
    mu_H=1/(60*365), #Assumes ~60 year lifespan
    k=0.1, #clumping parameter of the negative binomial distribution
    u_H=1200, #mL urine per human/day (approximate, ranges from 800 - 2000)
    
  #Transmission parameters
    lamda=9.27e-5, #snail-to-man transmission: p(infected snail sheds cercariae that infects human and reaches adulthood)
    beta=1.63e-5 #man-to-snail transmission: p(mated female worm produces a miracidia that infects a snail)

)

#Set initial values and do some runs##########
nstart=c(S=1,E=0,I=0, W=5, P=1)
yrs=40
time=seq(0,365*yrs,1)

output=as.data.frame(ode(nstart,time,schisto_master,parameters))
  output$N = output$S+output$E+output$I
eqbm=output[dim(output)[1],]

snail.prev=eqbm$I/(eqbm$E+eqbm$I+eqbm$S)

par(mfrow=c(1,1))
plot(output$time, output$S, type='l', xlab="time",ylab="System Variables", 
     ylim=c(0,max( output$S+output$E+output$I )), col='blue', lwd=2)
  lines(output$time,output$E,col='orange', lwd=2)
  lines(output$time,output$I,col='red', lwd=2)
  lines(output$time,output$N,col='black', lwd=2)
  lines(output$time,output$P*10,col = 'brown')
  
plot(output$time, output$W, type='l', xlab="time",ylab="Worm burden (W)", 
     col='purple', lwd=2)
  
#R0 function ##########
  
  get_Ro_beta_lamda<-function(muPq = 0, phi_Nq = 1, beta, lamda, f_Nq = 1) #variable parameters to be manipulated
  { 
    f_N<-parameters["f_N"]
    phi_N<-parameters["phi_N"]
    z<-parameters["z"]
    mu_N<-parameters["mu_N"]
    sigma<-parameters["sigma"]
    mu_I<-parameters["mu_I"]
    alpha<-parameters["alpha"]
    Th<-parameters["Th"]
    f_P<-parameters["f_P"]
    phi_P<-parameters["phi_P"]
    mu_P<-parameters["mu_P"]
    mu_W<-parameters["mu_W"]
    m<-parameters["m"]
    H<-parameters["H"]
    mu_H<-parameters["mu_H"]
    
    P_eq<-(1-((muPq+mu_P)/f_P))*phi_P #Equilibrium estimate of P given prawn predator parameters
    if(P_eq<0){
      P_eq=0
    }
    #Equilibrium estimate of N given snail parameters
    #Shorthand values to use in N_eq expression
    a= -(alpha*Th*f_N*f_Nq)/(phi_N*phi_Nq)
    b= f_N*f_Nq*alpha*Th - (f_N*f_Nq)/(phi_N*phi_Nq) - mu_N*alpha*Th 
    c= f_N*f_Nq - mu_N - alpha*P_eq
    
    if((b^2-4*a*c)<0){ #If prawn population sufficient to eliminate snails, N_eq=0
      N_eq1=0
    } else {
      N_eq1 <- (-b + sqrt(b^2-4*a*c)) / (2*a) #Function to solve quadratic expression for N_eq
    }
    
    if((b^2-4*a*c)<0){ #If prawn population sufficient to eliminate snails, N_eq=0
      N_eq2=0
    } else {
      N_eq2 <- (-b - sqrt(b^2-4*a*c)) / (2*a) #Function to solve quadratic expression for N_eq
    }
    
    if(N_eq1 > N_eq2){ #If prawn population sufficient to eliminate snails, N_eq=0
      N_eq=N_eq1
    } else {
      N_eq = N_eq2 #Function to solve quadratic expression for N_eq
    }
    
    pred<-(alpha*P_eq)/(1+(alpha*N_eq*Th))#death rate of snails due to predators given equilibrium estimates of P and N
    
    T1<-0.5*beta*H*N_eq
    T2<-lamda*sigma
    T3<- (mu_W+mu_H)*(mu_N+pred+sigma)*(mu_N+pred+mu_I)
    
    Ro_est <- sqrt((T1*T2)/T3)
    
    return(c(N_eq,P_eq,Ro_est ))
    
  }  
  
  p.dead = parameters["f_P"] - parameters["mu_P"] #muPq value at which predator population is eliminated 
  
#Try values of lamda and beta to get in range of ~10% infection and R0 ~3 ##################
trans.params<-data.frame('lamda'=rep(seq(from=5e-6, to = 5e-5, by = (5e-5 - 5e-6)/20), 
                                     times = length(seq(from=5e-6, to = 5e-5, by = (5e-5 - 5e-6)/20))),
                         'beta' =rep(seq(from=5e-6, to = 5e-5, by = (5e-5 - 5e-6)/20),
                                     each = length(seq(from=5e-6, to = 5e-5, by = (5e-5 - 5e-6)/20))),
                         'W' = 0,
                         'snail.prev' = 0,
                         'R0' = 0)  

for (i in 1:nrow(trans.params)){
  params=parameters
  params['lamda'] = trans.params[i,1]
  params['beta'] = trans.params[i,2]
  nstart=c(S=4000,E=2000,I=500, W=10, P=0)
  time=seq(from=0, to=50*365, by = 1)
  
  trans.output=as.data.frame(ode(nstart,time,schisto_master,params))
  
  trans.params[i,3] = trans.output[dim(trans.output)[1],5] #mean worm burden
  
  trans.params[i,4] = trans.output[dim(trans.output)[1],4] / #shedding snail prevalence
                                       
                      (trans.output[dim(trans.output)[1],2] + 
                         trans.output[dim(trans.output)[1],3] + 
                          trans.output[dim(trans.output)[1],4])
  
  trans.params[i,5] = get_Ro_beta_lamda(beta = trans.params[i,2], 
                                        lamda = trans.params[i,1],
                                        muPq = p.dead)[3]
  
  plot(trans.output$time, trans.output$S, type='l', xlab="time",ylab="System Variables", 
       ylim=c(0,max( trans.output$S+trans.output$E+trans.output$I )), 
       col='blue', lwd=2)
    lines(trans.output$time,trans.output$E,col='orange', lwd=2)
    lines(trans.output$time,trans.output$I,col='red', lwd=2)
    lines(trans.output$time, (trans.output$S+trans.output$E+trans.output$I), 
          col='black', lwd=2)
    legend("topright", legend=round(c(trans.output[dim(trans.output)[1],5],
                                trans.output[dim(trans.output)[1],4],
                                trans.output[dim(trans.output)[1],2]), digits=2), 
           lty = c(1,1,1), lwd = c(1,1,1), col = c("purple", "red", "blue"))
    legend("bottomright", legend = c(paste('lamda=', trans.params[i,1], sep=''),
                                     paste('beta=', trans.params[i,2], sep='')))
    
  print(i)
}

  as.numeric(trans.params$snail.prev*100)  
  
  plot(trans.params$R0, trans.params$W, 
       ylab = 'Mean worm burden (W)', xlab = 'R0 estimate')
    points(trans.params$R0[trans.params$lamda==max(trans.params$lamda)], 
           trans.params$W[trans.params$lamda==max(trans.params$lamda)], pch = 16, col='red')
    points(trans.params$R0[trans.params$beta==max(trans.params$beta)], 
           trans.params$W[trans.params$beta==max(trans.params$beta)], pch = 16, col='green')
    points(trans.params$R0[trans.params$lamda==min(trans.params$lamda)], 
           trans.params$W[trans.params$lamda==min(trans.params$lamda)], pch = 16, col='yellow')
    points(trans.params$R0[trans.params$beta==min(trans.params$beta)], 
           trans.params$W[trans.params$beta==min(trans.params$beta)], pch = 16, col='blue')
    legend('bottomright', legend = c(paste('min-beta=', min(trans.params$beta), sep=''),
                                     paste('max-beta=', max(trans.params$beta), sep=''),
                                     paste('min-lamda=', min(trans.params$lamda), sep=''),
                                     paste('max-lamda=', max(trans.params$lamda), sep='')),
           col = c('blue','green','yellow','red'), pch=c(16,16,16,16), cex=0.75)
    
#Heat map of worm burden across beta and lamda values ##################
ggplot(trans.params, aes(x=lamda, y=beta, fill=W))+
  theme_bw()+
  geom_tile(color='white', size=0.1)+
  scale_fill_continuous(low='green', high='red')

#Heat map of R0 across beta and lamda values ##################
ggplot(trans.params, aes(x=lamda, y=beta, fill=R0))+
  theme_bw()+
  geom_tile(color='white', size=0.1)+
  scale_fill_continuous(low='green', high='red')

#Line plot of snail prev response to beta/lamda values ######################  
plot(trans.params$beta[trans.params$lamda==min(trans.params$lamda)], 
     trans.params$snail.prev[trans.params$lamda==min(trans.params$lamda)], 
     type='l', ylab = 'Snail.prev', xlab = 'beta', ylim = c(0,0.2))
  lines(trans.params$beta[trans.params$lamda==unique(trans.params$lamda)[6]], 
        trans.params$snail.prev[trans.params$lamda==unique(trans.params$lamda)[6]], 
        col=2)
  lines(trans.params$beta[trans.params$lamda==unique(trans.params$lamda)[11]], 
        trans.params$snail.prev[trans.params$lamda==unique(trans.params$lamda)[11]], 
        col=3)
  lines(trans.params$beta[trans.params$lamda==unique(trans.params$lamda)[16]], 
        trans.params$snail.prev[trans.params$lamda==unique(trans.params$lamda)[16]], 
        col=4)
  legend('bottomright', legend = c(paste('lamda=', min(trans.params$lamda), sep=''),
                                   paste('lamda=', unique(trans.params$lamda)[6], sep=''),
                                   paste('lamda=', unique(trans.params$lamda)[11], sep=''),
                                   paste('lamda=', unique(trans.params$lamda)[16], sep='')),
         col=c(1,2,3,4), lty=c(1,1,1,1), cex=0.75)
 
#Line plot of worm burden response to beta/lamda values ######################  
plot(trans.params$beta[trans.params$lamda==min(trans.params$lamda)], 
       trans.params$W[trans.params$lamda==min(trans.params$lamda)], 
       type='l', ylab = 'W', xlab = 'beta', ylim=c(0,65))
  lines(trans.params$beta[trans.params$lamda==unique(trans.params$lamda)[6]], 
        trans.params$W[trans.params$lamda==unique(trans.params$lamda)[6]], 
        col=2)
  lines(trans.params$beta[trans.params$lamda==unique(trans.params$lamda)[11]], 
        trans.params$W[trans.params$lamda==unique(trans.params$lamda)[11]], 
        col=3)
  lines(trans.params$beta[trans.params$lamda==unique(trans.params$lamda)[16]], 
        trans.params$W[trans.params$lamda==unique(trans.params$lamda)[16]], 
        col=4)
  legend('bottomright', legend = c(paste('lamda=', min(trans.params$lamda), sep=''),
                                   paste('lamda=', unique(trans.params$lamda)[6], sep=''),
                                   paste('lamda=', unique(trans.params$lamda)[11], sep=''),
                                   paste('lamda=', unique(trans.params$lamda)[16], sep='')),
         col=c(1,2,3,4), lty=c(1,1,1,1), cex=0.575)
  
#Line plot of R0 response to beta/lamda values ######################  
plot(trans.params$beta[trans.params$lamda==min(trans.params$lamda)], 
       trans.params$R0[trans.params$lamda==min(trans.params$lamda)], 
       type='l', ylab = 'R0', xlab = 'beta', ylim=c(0,4))
  lines(trans.params$beta[trans.params$lamda==unique(trans.params$lamda)[6]], 
        trans.params$R0[trans.params$lamda==unique(trans.params$lamda)[6]], 
        col=2)
  lines(trans.params$beta[trans.params$lamda==unique(trans.params$lamda)[11]], 
        trans.params$R0[trans.params$lamda==unique(trans.params$lamda)[11]], 
        col=3)
  lines(trans.params$beta[trans.params$lamda==unique(trans.params$lamda)[16]], 
        trans.params$R0[trans.params$lamda==unique(trans.params$lamda)[16]], 
        col=4)
  legend('bottomright', legend = c(paste('lamda=', min(trans.params$lamda), sep=''),
                                   paste('lamda=', unique(trans.params$lamda)[6], sep=''),
                                   paste('lamda=', unique(trans.params$lamda)[11], sep=''),
                                   paste('lamda=', unique(trans.params$lamda)[16], sep='')),
         col=c(1,2,3,4), lty=c(1,1,1,1), cex=0.75)

#Check out how mating probability varies with W assuming constant k ############
  Ws<-c(1:100)
  k=0.08
  mating<-rep(0,100)
  
  for(i in 1:100){
    W = Ws[i]
    mating[i] = 1 - ((1-(W/(W+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value)
  }
  
  mating2<-rep(0,100)
  k=0.15
  for(i in 1:100){
    W = Ws[i]
    mating2[i] = 1 - ((1-(W/(W+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value)
  }
  
  mating3<-rep(0,100)
  k=0.25
  for(i in 1:100){
    W = Ws[i]
    mating3[i] = 1 - ((1-(W/(W+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value)
  }
  
  plot(Ws, mating, type = 'l', xlab = 'mean worm burden (W)', ylab = 'mating probability', lwd=2)
    lines(Ws, mating2, col='red', lwd=2)
    lines(Ws, mating3, col='blue', lwd=2)
    legend('bottomright', legend = c('k=0.08', 'k=0.15', 'k=0.25'), col = c('black', 'red', 'blue'),
           lwd=2)
    
    
#Estimate z parameter from Mangal et al 2010 Fig3 #####################
  fxE<-function(dens){
    208.7*exp(-0.0797*dens)
  }
  
  fxS<-function(dens){
    231.03*exp(-0.0524*dens)
  }
  
  z = integrate(fxE, 3, 17)$value / integrate(fxS, 3, 17)$value
  
  
  