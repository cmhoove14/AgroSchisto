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
    dSdt= (f_N)*((1-((phi_N*(1/phi_Nq))*N))*(S+(z*E))) - 
      ((mu_N)*S) - (pred*S) - (beta*(M)*S) #Susceptible snails
    
    dEdt= beta*(M)*S - ((mu_N+pred+sigma)*E) #Exposed snails
    
    dIdt= (sigma*E) - ((mu_N+pred+mu_I)*I) #Infected snails
    
    dWdt= (lamda*I) - ((mu_W+mu_H)*W) #worm burden in human population
    
    dPdt= (f_P*(1-(phi_P*P))*P)-((mu_P+muPq)*P) #prawn population
    
   
    return(list(c(dSdt,dEdt,dIdt,dWdt, dPdt)))
  }) 
} 

#List parameters and values #####################
parameters=c(
  ##standard snail parameters 
    f_N=0.16, # recruitment rate (from sokolow et al)
    phi_N=(1-1/(80*0.16))/10000, # carrying capacity from sokolow et al
    z=0.5, #Proportion of exposed snails that reproduce from sokolow et al
    mu_N=1/80, #Mortality rate from Sokolow et al
    sigma=1/50, #Transition rate from exposed to infected from sokolow et al
    mu_I=1/20, #additional snail death due to infection from sokolow et al
    
    beta=0.000004, #0.0000011128-Best fit data to first two points in Lampsar 2

  ## snail parameters impacted by agrochemicals
    #f_Nq=1, #Not affected in mesocosm
    phi_Nq=1, #Scalar of snail carrying capacity by chemical concentration INFORMED BY BOTTOM UP EFFECTS IN MESOCOSM
    #mu_Nq=0, #Chem concentrations too low to affect snails in mesocosom experiments
  
  #prawn parameters
    alpha=0.003, #attack rate
    Th=0.1,#~Prawn predation limit
    f_P=0.128,#prawn birth rate
    phi_P=1/50,  #prawn carrying capacity
    mu_P= 0.00137,#daily prawn mortality
  
  #prawn parameter impacted by agrochemicals
    muPq=0, #agrochemical dependent prawn mortality rate (added to baseline mortality rate) INFORMED BY MESOCOSM AND HALSTEAD 2015
    #alpha_q=1, #Scalar of predation rate due to sub-lethal toxicity; not considered in mesocosm
  
  #Adult Worm, Miracidia and Circariae Parameters
  lamda=0.00004, #probability of snail shedding a cercariae that infects a human host and survives to reproduction
  mu_W=1/(3*365), # death rate of adult worms
  m=0.5, #miracidial shedding rate per reproductive female divided by miracidial mortality; from sokolow et al
  #V=.084, #Egg viability controlling schistosome egg->infective miracidia

  
  #Human parameters
  H=1000, #number of humans
  mu_H=1/(60*365) #Assumes 60 year lifespan
)
#Set initial values #########
    yrs=30
    time=seq(0,365*yrs,1)
    p=1 
    nstart=c(S=7000,E=3750,I=1200, W=20, P=p)
  
#Run model with base parameters and no agrochemical effects and save equilibrium values for 5 state variables #############
  output.b=as.data.frame(ode(nstart,time,schisto_halstead,parameters)) 
  
  eqbm.b<-output.b[365*yrs,]
  eqbm.b #Infection mostly prevented by prawn population

    #plot time series of state variables
      plot(output.b[,1], output.b[,2], type='l', xlab="time",ylab="System Variables", 
           ylim=c(0,max( output.b[,2]+output.b[,3]+output.b[,4] )), col=1, lty=1, lwd=2)
          lines(output.b[,1],output.b[,3],col=2, lty=2, lwd=2)
          lines(output.b[,1],output.b[,4],col=3, lty=2, lwd=2)
          lines(output.b[,1],output.b[,2]+output.b[,3]+output.b[,4],col=4, lty=2, lwd=2)
          lines(output.b[,1],output.b[,5],col=5, lty=2)
          lines(output.b[,1],output.b[,6],col=6, lty=2, lwd=2)
  
  #Plot worm burden and prawn numbers onlyonly  
    plot(output.b[,1],output.b[,5],col=5, ylab='Mean worm burden (W)', xlab='time', ylim=c(0,50))
      lines(output.b[,1],output.b[,6],col=6, lty=2, lwd=2)
      
#Add starting prawn numbers and halt reproduction to reflect mesocosm ################
  p=3*40
    nstart=c(S=7000,E=3750,I=1200, W=20, P=p)
  parameters["f_P"]<-0.128
  parameters["mu_P"]<-0.006862177 #Observed daily mortality over 12 weeks (see below for derivation)
  parameters["phi_P"]<-1/120 #Really high to eliminate its influence
  
  output.m=as.data.frame(ode(nstart,time,schisto_halstead,parameters)) 
  
  eqbm.m<-output.m[365*yrs,]
  eqbm.m
  
  #plot time series of state variables
    plot(output.m[,1], output.m[,2], type='l', xlab="time",ylab="System Variables", 
         ylim=c(0,max( output.m[,2]+output.m[,3]+output.m[,4] )), col=1, lty=1, lwd=2)
      lines(output.m[,1],output.m[,3],col=2, lty=2, lwd=2)
      lines(output.m[,1],output.m[,4],col=3, lty=2, lwd=2)
      lines(output.m[,1],output.m[,2]+output.m[,3]+output.m[,4],col=4, lty=2, lwd=2)
      lines(output.m[,1],output.m[,5],col=5, lty=2)
      lines(output.m[,1],output.m[,6],col=6, lty=2, lwd=2)
  
  #Plot worm burden and prawn numbers onlyonly  
    plot(output.m[,1],output.m[,5],col=5, ylab='Mean worm burden (W)', xlab='time', ylim=c(0,50))
      lines(output.m[,1],output.m[,6],col=6, lty=2, lwd=2)

#Add muPq from mesocosm #######################
  p=3*39
    nstart=c(S=7000,E=3750,I=1200, W=20, P=p)
  parameters["f_P"]<-0.128
  parameters["mu_P"]<-0.03883984
  parameters["phi_P"]<-1/(3*40)
  parameters["muPq"]<-1.321756-parameters["mu_P"]
    output.q=as.data.frame(ode(nstart,time,schisto_halstead,parameters)) 
    
    eqbm.q<-output.q[365*yrs,]
    eqbm.q

#plot time series of estimated values to make sure state variables reached equilibrium 
  plot(output.q[,1], output.q[,2], type='l', xlab="time",ylab="System Variables", 
       ylim=c(0,max( output.q[,2]+output.q[,3]+output.q[,4] )), col=1, lty=1, lwd=2)
    lines(output.q[,1],output.q[,3],col=2, lty=2, lwd=2)
    lines(output.q[,1],output.q[,4],col=3, lty=2, lwd=2)
    lines(output.q[,1],output.q[,2]+output.q[,3]+output.q[,4],col=4, lty=2, lwd=2)
    lines(output.q[,1],output.q[,5],col=5, lty=2)
    lines(output.q[,1],output.q[,6],col=6, lty=2, lwd=2)
    
  #Plot worm burden only  
    plot(output.q[,1],output.q[,5],col=5, ylab='Mean worm burden (W)', xlab='time') 
      lines(output.q[,1], output.q[,6], lty=2, col=6)

#A couple different options for muPq #################################
  #muPq data from halstead et al 2015 (ecotoxicology paper)
  
  ecotox<-data.frame(dose=c(0,0.64,3.2,6.4,32,64), 
                    mortality=c(0,0,0,0,0.8, 1.00))
  
  ecotox$mu_P<- -0.25*log(1-ecotox[,2])
  ecotox[6,3]= -0.25*log(1-0.95) #Assume 95% mortality instead of full mortality for conservative estimate
  
  #Probit analysis of contrived data (same pattern as in observed, but highest dose group has 95% mortality instead of 100% and used 100 in each dosage group instead of 5)
    ecotox_mod1<-data.frame('dose'=c(rep(0,100), rep(0.64,100), rep(3.2,100), rep(6.4,100), rep(32,100), rep(64,100)),
                          'response'=c(rep(0,100), rep(0,100), rep(0,100), rep(0,100), rep(0,20), rep(1,175), rep(0,5)))
      ecotox1<-glm(response ~ dose, family=binomial(link="probit"),data=ecotox_mod1)
    summary(ecotox1)
    
  #Extrapolate response to constant gradient of Chlorpyrifos concentration
    p.ecotox<-data.frame(dose=seq(from=0, to=150, by=1))
    p.ecotox[, c('mortality', 'st.er')]<-predict(ecotox1, p.ecotox, 
                                                    type = "response", se.fit=TRUE)
    
  #Actual observed data with 5 prawns in each dose group  
    ecotox_mod2<-data.frame('dose'=c(rep(0,5), rep(0.64,5), rep(3.2,5), rep(6.4,5), rep(32,5), rep(64,5)),
                           'response'=c(rep(0,5), rep(0,5), rep(0,5), rep(0,5), 0 , rep(1,9)))
    ecotox2<-glm(response ~ dose, family=binomial(link="probit"),data=ecotox_mod2)
      summary(ecotox2)
    
      #Extrapolate response to constant gradient of Chlorpyrifos concentration
      p.ecotox2<-data.frame(dose=seq(from=0, to=150, by=1))
      p.ecotox2[, c('mortality', 'st.er')]<-predict(ecotox2, p.ecotox2, 
                                                          type = "response", se.fit=TRUE)
  
  plot(x=p.ecotox$dose, y=p.ecotox$mortality, type='l', ylim=c(0,1.2), xlab='ChlorP dose', ylab='predicted % daily mortality')
    lines(x=p.ecotox$dose, y=p.ecotox$mortality+p.ecotox$st.er, lty=2, cex=0.8)
    lines(x=p.ecotox$dose, y=p.ecotox$mortality-p.ecotox$st.er, lty=2, cex=0.8)
    lines(x=p.ecotox2$dose, y=p.ecotox2$mortality, col='red')
    lines(x=p.ecotox2$dose, y=p.ecotox2$mortality+p.ecotox2$st.er, lty=2, col='red', cex=0.8)
    lines(x=p.ecotox2$dose, y=p.ecotox2$mortality-p.ecotox2$st.er, lty=2, col='red', cex=0.8)
      points(ecotox$dose, ecotox$mu_P, pch=16)
  
  #Convert %mortality to mortality rate
  muPq1<-data.frame('mu_agro_p1' = -0.25*log(1-p.ecotox[,2]),
                    'mu_agro_p1_up' = -0.25*log(1-(p.ecotox[,2]+p.ecotox[,3])),
                    'mu_agro_p1_lo' = -0.25*log(1-(p.ecotox[,2]-p.ecotox[,3])))
  
  muPq2<-data.frame('mu_agro_p2' = -0.25*log(1-p.ecotox2[,2]),
                    'mu_agro_p2_up' = -0.25*log(1-(p.ecotox2[,2]+p.ecotox2[,3])),
                    'mu_agro_p2_lo' = -0.25*log(1-(p.ecotox2[,2]-p.ecotox2[,3])))
  
  

#muPq directly from mesocosm ###############################
  dat<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/R_use.csv")
    
    #What is the observed daily prawn mortality rate in chlorP-free tanks over the full 12 weeks?
      p.0<-3*length(dat$tank[dat$chlor==0]) #Total starting number of prawns in chloP-free tanks (3 in each)
      p.12<-sum(dat$p.all_fin[dat$chlor==0]) #Total surviving prawns in chlorP tanks at end of 12-week experiment
    #daily mortality rate assuming constant death throughout 12 weeks =ln(Nt/N0)/-t with t=12 weeks *7 days=84 days
      p.r<-log(p.12/p.0)/-84 #=0.006862177
      
  #What about using the 24 hour endpoints?  
    mesotox<-data.frame('chlorP'=dat$chlor,
                          'chlorP2'=dat$treat,
                          'prawn.0'=rep(3,length(dat$chlor)),
                          'prawn.24'=dat$p.all_24)
    
      mesotox$chlorP[mesotox$chlorP2=="C 2x"]<-2 #Fill info for double EEC chlorP tanks
      
      mesotox<-mesotox[,-2] #Get rid of now unneccessary column
    
    #Prepare data for probit where each observation is a prawn, its binary survival outcome (1=dead) and a chlorPdose  
      mesotox2<-untable(mesotox[,c(1,3)], num=mesotox[,2])
      
      mesotox2$dose[mesotox2$chlorP==2]<-128
      mesotox2$dose[mesotox2$chlorP==1]<-64
      mesotox2$dose[mesotox2$chlorP==0]<-0
    
      mesotox2$outcome<-rep(0, length(mesotox2$chlorP))
      
      mesotox2$outcome[mesotox2$prawn.24==0]<-rep(c(1,1,1), length(mesotox2$prawn.24[mesotox2$prawn.24==0])/3)
      mesotox2$outcome[mesotox2$prawn.24==1]<-rep(c(1,1,0), length(mesotox2$prawn.24[mesotox2$prawn.24==1])/3)
      mesotox2$outcome[mesotox2$prawn.24==2]<-rep(c(1,0,0), length(mesotox2$prawn.24[mesotox2$prawn.24==2])/3)
      mesotox2$outcome[mesotox2$prawn.24==3]<-rep(c(0,0,0), length(mesotox2$prawn.24[mesotox2$prawn.24==3])/3)
    
    #~*~*~*~DAILY MORTALITY RATES FOR PRAWN SPECIES~*~*~*~*~
      sum(mesotox2$outcome[mesotox2$chlorP==2])/length(mesotox2$outcome[mesotox2$chlorP==2]) #60% mortality
        -log(1-0.6) #Observed daily mortality rate = 0.9162907 when ChlorP @128 ug/L is present
      
      sum(mesotox2$outcome[mesotox2$chlorP==1])/length(mesotox2$outcome[mesotox2$chlorP==1]) #76.6% mortality 
        -log(1-0.766667) #Observed daily mortality rate = 1.455289 when ChlorP @64 ug/L is present
      
      sum(mesotox2$outcome[mesotox2$chlorP!=0])/length(mesotox2$outcome[mesotox2$chlorP!=0]) #73.3% mortality when chlorP present
        -log(1-0.7333333) #Observed daily mortality rate = 1.321756 when ChlorP is present
      
      sum(mesotox2$outcome[mesotox2$chlorP==0])/length(mesotox2$outcome[mesotox2$chlorP==0]) #3.8% mortality when chlorP absent
        -log(1-0.03809524) #Observed daily mortality rate = 0.03883984 when ChlorP is absent
    
    #Probit model with three tested doses
      pr<-glm(outcome ~ dose, family=binomial(link="probit"),data=mesotox2)
      summary(pr)
      
      p.mesotox<-data.frame('dose'=seq(1,150,1),
                              'mort'=rep(0,150),
                              'st.er'=rep(0,150))
    
    p.mesotox$mort<-predict(pr, p.mesotox, type = "response", se.fit=TRUE)$fit
    p.mesotox$st.er<-predict(pr, p.mesotox, type = "response", se.fit=TRUE)$se.fit
    
    #Plot modeled dose-response curve
      plot(p.mesotox$dose, p.mesotox$mort, type='l', lwd=1.5, ylim=c(0,1),
           xlab='Dose', ylab = '% mortality')
        lines(p.mesotox$dose, p.mesotox$mort+p.mesotox$st.er, col='red', lwd=0.8, lty=2)
        lines(p.mesotox$dose, p.mesotox$mort-p.mesotox$st.er, col='red', lwd=0.8, lty=2)
        points(x=c(0,64,128), y=c(0.038,0.73333,0.60), pch=16)
    
#Values of agrochemical parameters (phi_Nq and muPq@chlorP=64 ug/L) to sample over within R0 expression #################
  #phi_Nq variability derived from mesocosom data in Halstead_et_al_bottom_up_predict code
    phi_base<-rnorm(n=10000, mean=1, sd=(0.1488876))
    phi_fert<-rnorm(n=10000, mean=1.159642, sd=(0.1330912))
    phi_atra<-rnorm(n=10000, mean=1.614304, sd=(0.1223979))
    phi_atfe<-rnorm(n=10000, mean=1.509579, sd=(0.1711861))
    
  #muPq options (derived above)  
    #Observed mortality in ecotoxicology paper was 100% (though it was 4 day mortality which is why -0.25 is in calcs below)
      eco.mort.obs<-rep(1, 10000) 
      eco.rate<-rep(9.010913347, 10000)
    #modeled with contrived prawns and 95% in last dosage group
      eco.mort.mod1<-rnorm(n=10000, mean=p.ecotox$mortality[p.ecotox$dose==64], sd=p.ecotox$st.er[p.ecotox$dose==64]) 
        eco.mort.mod1[eco.mort.mod1>1]<-1.000000e+00 #Max mortality is 1 (100%)
          eco.rate.mod1<- -0.25*log(1-eco.mort.mod1) #Convert %mortality into death rate
          eco.rate.mod1[eco.rate.mod1==Inf]<-9.010913347 #This seems to be the maximum mortality rate, maybe based on R's limit of smallest number
    #modeled with observed 5 prawns in each dose group and 100% mortality in 64 dose group      
    eco.mort.mod2<-rnorm(n=10000, mean=p.ecotox2$mortality[p.ecotox2$dose==64], sd=p.ecotox2$st.er[p.ecotox2$dose==64])
      eco.mort.mod2[eco.mort.mod2>1]<-1.000000e+00 #Max mortality is 1 (100%)
        eco.rate.mod2<- -0.25*log(1-eco.mort.mod2) #Convert %mortality into death rate
        eco.rate.mod2[eco.rate.mod1==Inf]<-9.010913347
    
    #Observed rate from mesocosm
      meso.rate<-rep(1.455289, 10000)
    #Modeled rate from dose response function to 0, 64, and 128 in mesocosm
      meso.mort.mod1<-rnorm(n=10000, mean=p.mesotox$mort[p.mesotox$dose==64], sd=p.mesotox$st.er[p.mesotox$dose==64])
        meso.mort.mod1[meso.mort.mod1>1]<-1.000000e+00 #Max mortality is 1 (100%)
          meso.rate.mod1<- -log(1-meso.mort.mod1) #This was daily mortality so no 0.25 in calculation
#R0 function *~*~*GO BACK TO ORIGINAL PARAMETER VALUES BEFORE RUNNING ANYTHING BELOW*~*~*~* #########
    
get_Ro<-function(muPq, phi_Nq)
{   #HAVE TO SET muPq and phi_Nq in function call
  f_N<-parameters["f_N"]
  phi_N<-parameters["phi_N"]
  z<-parameters["z"]
  mu_N<-parameters["mu_N"]
  beta<-parameters["beta"]
  sigma<-parameters["sigma"]
  mu_I<-parameters["mu_I"]
  alpha<-parameters["alpha"]
  Th<-parameters["Th"]
  f_P<-parameters["f_P"]
  phi_P<-parameters["phi_P"]
  mu_P<-parameters["mu_P"]
  lamda<-parameters["lamda"]
  mu_W<-parameters["mu_W"]
  m<-parameters["m"]
  H<-parameters["H"]
  mu_H<-parameters["mu_H"]
  
  muPq=muPq
  phi_Nq=phi_Nq
  
  N_eq<-(1-(mu_N/f_N))/(phi_N*(1/phi_Nq)) #Equilibrium estimate of N given snail parameters
  P_eq<-(1-((mu_P+muPq)/f_P))/phi_P #Equilibrium estimate of P given prawn predator parameters
    if(P_eq<0){
      P_eq=0
    }
  pred<-(alpha*P_eq)/(1+(alpha*N_eq*Th))#death rate of snails due to predators given equilibrium estimates of P and N
  
  T1<-0.5*beta*m*H*N_eq
  T2<-lamda*sigma
  T3<- (mu_W+mu_H)*(mu_N+pred+sigma)*(mu_N+pred+mu_I)
  
  Ro_est <- sqrt((T1*T2)/T3)
  
  Ro_est 
  
}

#Baseline R0 run ##############
  get_Ro(muPq=0, phi_Nq=1)  #R0 = 1.09706
          
#Prawn free R0 run ##########
  parameters["phi_P"]<-1e20
  get_Ro(muPq = 0, phi_Nq = 1) #R0 = 2.026636

#R0 run with prawn carrying capacity equivalent to that in mesocosm assuming scale up factor of 40 (5 sq meters to 200 sq meters) #########
  parameters["phi_P"]<-1/120
  get_Ro(muPq = 0, phi_Nq = 1) #R0 = 0.6730459
#Mesocosm informed R0 runs#############   
  #Mesocosm baseline R0 run
    R0_base<-rep(0,10000)
    for(i in 1:length(R0_base)){
      R0_base[i]<-get_Ro(muPq = 0, phi_Nq = phi_base[i])
    }
  #Mesocosm +Fertilizer R0 run
    R0_fert<-rep(0,10000)
    for(i in 1:length(R0_fert)){
      R0_fert[i]<-get_Ro(muPq = 0, phi_Nq = phi_fert[i])
    }
   
  #Mesocosm +Atrazine R0 run
    R0_atra<-rep(0,10000)
    for(i in 1:length(R0_atra)){
      R0_atra[i]<-get_Ro(muPq = 0, phi_Nq = phi_atra[i])
    }
    
  #Mesocosm +Chlorpyrifos R0 run
    R0_chlor<-rep(0,10000)
    for(i in 1:length(R0_chlor)){
      R0_chlor[i]<-get_Ro(muPq = meso.rate.mod1[i], phi_Nq = phi_base[i])
    }
    
  #Mesocosm +Fertilizer+Atrazine R0 run
    R0_atfe<-rep(0,10000)
    for(i in 1:length(R0_atfe)){
      R0_atfe[i]<-get_Ro(muPq = 0, phi_Nq = phi_atfe[i])
    }
    
  #Mesocosm +Fertilizer+Chlorpyrifos R0 run
    R0_chfe<-rep(0,10000)
    for(i in 1:length(R0_chfe)){
      R0_chfe[i]<-get_Ro(muPq = meso.rate.mod1[i], phi_Nq = phi_fert[i])
    }
    
  #Mesocosm +Atrazine+Chlorpyrifos R0 run
    R0_atch<-rep(0,10000)
    for(i in 1:length(R0_atch)){
      R0_atch[i]<-get_Ro(muPq = meso.rate.mod1[i], phi_Nq = phi_atra[i])
    }
    
  #Mesocosm +Fertilizer+Atrazine+Chlorpyrifos R0 run
    R0_tre<-rep(0,10000)
    for(i in 1:length(R0_tre)){
      R0_tre[i]<-get_Ro(muPq = meso.rate.mod1[i], phi_Nq = phi_atfe[i])
    }
  
          
#Plot results ###############

r0s<-data.frame("Treatment"=c('baseline', 'Fert only', 'Atra only', 'ChlorP Only', 
                              'Fert+Atra', 'Fert+ChlorP', 'Atra+ChlorP', 'All three'),
                      "meanR0"=c(mean(R0_base), mean(R0_fert), mean(R0_atra), mean(R0_chlor),
                                 mean(R0_atfe), mean(R0_chfe), mean(R0_atch), mean(R0_tre)),
                      "st.erR0"=c(st.er(R0_base), st.er(R0_fert), st.er(R0_atra), st.er(R0_chlor),
                                  st.er(R0_atfe), st.er(R0_chfe), st.er(R0_atch), st.er(R0_tre)))

r0s$Treatment<-factor(r0s$Treatment, levels = c('baseline', 'Fert only', 'Atra only', 'ChlorP Only', 
                                                'Fert+Atra', 'Fert+ChlorP', 'Atra+ChlorP', 'All three'))

ggplot(r0s, aes(x=Treatment, y=meanR0, fill=Treatment))+
  theme_bw()+
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=15))+
  scale_fill_manual(values=cbPalette) +
  ylab("Mean +/- SEM predicted R-0")+
  geom_bar(position=position_dodge(), stat="identity", width = .7) +
  geom_errorbar(aes(ymin=meanR0-st.erR0,
                    ymax=meanR0+st.erR0),
                width=.2, position=position_dodge(.7))
