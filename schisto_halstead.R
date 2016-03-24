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
    
    beta=0.0000011128, #-Best fit data to first two points in Lampsar 2

  ## snail parameters impacted by agrochemicals
    #f_Nq=1, #Not affected in mesocosm
    phi_Nq=1, #Scalar of snail carrying capacity by chemical concentration INFORMED BY BOTTOM UP EFFECTS IN MESOCOSM
    #mu_Nq=0, #Chem concentrations too low to affect snails in mesocosom experiments
  
  #prawn parameters
    alpha=0.003, #attack rate
    Th=0.1,#~Prawn predation limit
    f_P=0.128,#prawn birth rate
    phi_P=1/50,  #prawn carrying capacity
    mu_P= 0.03883984, #observed daily 24 hr prawn mortality rate in chlorP-free tanks of mesocosm#0.00137, previous rate
  
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
      

#Add muPq from mesocosm #######################
  p=3*39
    nstart=c(S=7000,E=3750,I=1200, W=6, P=p)
  parameters["mu_P"]<-0.03883984
    mu_P=parameters["mu_P"]
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
                          'response'=c(rep(0,100), rep(0,100), rep(0,100), rep(0,100), rep(0,20), rep(1,179), rep(0,1)))
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
    
  #What are the mean daily mortality rates in each group derived from inidividual tanks?    
    mesotox$surv.24<-mesotox$prawn.24/mesotox$prawn.0 #Equivalent to (Nt/N0)
      mesotox$surv.24[mesotox$surv.24==0]<-0.15
    mesotox$rate.24<- -log(mesotox$surv.24) #-log(Nt/N0)=r
    
      mesorate.chlorP0<-mean(mesotox$rate.24[mesotox$chlorP==0])
        mesorate.chlorP0.sd<-sd(mesotox$rate.24[mesotox$chlorP==0])
      
      mesorate.chlorP64<-mean(mesotox$rate.24[mesotox$chlorP==1])
        mesorate.chlorP64.sd<-sd(mesotox$rate.24[mesotox$chlorP==0])
      
      mesorate.chlorP128<-mean(mesotox$rate.24[mesotox$chlorP==2])
        mesorate.chlorP128.sd<-sd(mesotox$rate.24[mesotox$chlorP==2])
      
      mesorate.chlorPres<-mean(mesotox$rate.24[mesotox$chlorP>0])
        mesorate.chlorPres.sd<-sd(mesotox$rate.24[mesotox$chlorP>0])
        
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
    
        
        
#OBSOLETE - Values of agrochemical parameters (phi_Nq and muPq@chlorP=64 ug/L) to sample over within R0 expression #################
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
          
    #Mean, sd of daily prawn mortality from individual chlorP containing tanks
      meso.rate.tanks<-rnorm(10000, mean=mesorate.chlorPres, sd=mesorate.chlorPres.sd)
        meso.rate.tanks[meso.rate.tanks<(2*parameters["mu_P"])]<-(2*parameters["mu_P"])
        
        
#R0 function#########
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
  
  N_eq<-(1-(mu_N/f_N))/(phi_N*(1/phi_Nq)) #Equilibrium estimate of N given snail parameters
  P_eq<-(1-((muPq+mu_P)/f_P))/phi_P #Equilibrium estimate of P given prawn predator parameters
    if(P_eq<0){
      P_eq=0
    }
  pred<-(alpha*P_eq)/(1+(alpha*N_eq*Th))#death rate of snails due to predators given equilibrium estimates of P and N
  
  T1<-0.5*beta*m*H*N_eq
  T2<-lamda*sigma
  T3<- (mu_W+mu_H)*(mu_N+pred+sigma)*(mu_N+pred+mu_I)
  
  Ro_est <- sqrt((T1*T2)/T3)
  
  print(N_eq)
  print(P_eq)
  Ro_est 
  
}

#Baseline R0 run ##############
    parameters["phi_P"]<-1/50
    parameters["beta"]<-1.1128e-6      
  get_Ro(muPq=0, phi_Nq=1)  #R0 = 0.6684057 IF DIFFERENT, SOMETHING ISN'T RIGHT
          

#Prawn free R0 run gives our estimate of the R0 in Lampsar 2 village##########
    parameters["phi_P"]<-1e20 #Basically makes equilibrium prawns 0
      get_Ro(muPq = 0, phi_Nq = 1) #R0 = 1.068942
        R0_vil<-1.068942
  #reset beta to point estimate
    parameters["beta"]=0.0000011128
    

#R0 run with prawn carrying capacity equivalent to that in mesocosm assuming scale up factor of 40 (5 sq meters to 200 sq meters) #########
  parameters["phi_P"]<-1/120 #Also want to use this in the mesocosm runs below
  get_Ro(muPq = 0, phi_Nq = 1) #R0 = 0.4413189

#OBSOLETE - Mesocosm informed R0 runs with muPq estimated from mesocosm dose-response#############   
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
      R0_chlor[i]<-get_Ro(muPq = meso.rate.mod1[i]-mu_P, phi_Nq = phi_base[i])
    }
    
  #Mesocosm +Fertilizer+Atrazine R0 run
    R0_atfe<-rep(0,10000)
    for(i in 1:length(R0_atfe)){
      R0_atfe[i]<-get_Ro(muPq = 0, phi_Nq = phi_atfe[i])
    }
    
  #Mesocosm +Fertilizer+Chlorpyrifos R0 run
    R0_chfe<-rep(0,10000)
    for(i in 1:length(R0_chfe)){
      R0_chfe[i]<-get_Ro(muPq = meso.rate.mod1[i]-mu_P, phi_Nq = phi_fert[i])
    }
    
  #Mesocosm +Atrazine+Chlorpyrifos R0 run
    R0_atch<-rep(0,10000)
    for(i in 1:length(R0_atch)){
      R0_atch[i]<-get_Ro(muPq = meso.rate.mod1[i]-mu_P, phi_Nq = phi_atra[i])
    }
    
  #Mesocosm +Fertilizer+Atrazine+Chlorpyrifos R0 run
    R0_tre<-rep(0,10000)
    for(i in 1:length(R0_tre)){
      R0_tre[i]<-get_Ro(muPq = meso.rate.mod1[i]-mu_P, phi_Nq = phi_atfe[i])
    }
  
          
    
#OBSOLETE - Plot results ###############

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
  scale_fill_manual(values=c("#999999", 'green2', 'gold1', 'red', 'orange', 'darkgreen', 'pink', 'blue')) +
  ylab("Mean +/- SEM predicted R-0")+
  geom_bar(position=position_dodge(), stat="identity", width = .7) +
  geom_errorbar(aes(ymin=meanR0-st.erR0,
                    ymax=meanR0+st.erR0),
                width=.2, position=position_dodge(.7))

#OBSOLETE - Use box plot instead of bar plot ###########################
  r0s.box<-data.frame('At'=R0_atra,
                      'Ch'=R0_chlor,
                      'Fe'=R0_fert,  
                      'At:Fe'=R0_atfe, 
                      'Ch:Fe'=R0_chfe, 
                      'At:Ch'=R0_atch, 
                      'At:Ch:Fe'=R0_tre)

  boxplot(r0s.box)
    lines()

#OBSOLETE - Mesocosm informed R0 runs with muPq value from mesorate.chlorPres value derived above##########
  #Mesocosm baseline R0 run
    R0_base.2<-rep(0,10000)
    for(i in 1:length(R0_base)){
      R0_base[i]<-get_Ro(muPq = 0, phi_Nq = phi_base[i])
    }
    
  #Mesocosm +Fertilizer R0 run
    R0_fert.2<-rep(0,10000)
    for(i in 1:length(R0_fert)){
      R0_fert[i]<-get_Ro(muPq = 0, phi_Nq = phi_fert[i])
    }
    
  #Mesocosm +Atrazine R0 run
    R0_atra.2<-rep(0,10000)
    for(i in 1:length(R0_atra)){
      R0_atra[i]<-get_Ro(muPq = 0, phi_Nq = phi_atra[i])
    }
    
  #Mesocosm +Chlorpyrifos R0 run
    R0_chlor.2<-rep(0,10000)
    for(i in 1:length(R0_chlor)){
      R0_chlor[i]<-get_Ro(muPq = meso.rate.tanks[i]-mu_P, phi_Nq = phi_base[i])
    }
    
  #Mesocosm +Fertilizer+Atrazine R0 run
    R0_atfe.2<-rep(0,10000)
    for(i in 1:length(R0_atfe)){
      R0_atfe[i]<-get_Ro(muPq = 0, phi_Nq = phi_atfe[i])
    }
    
  #Mesocosm +Fertilizer+Chlorpyrifos R0 run
    R0_chfe.2<-rep(0,10000)
    for(i in 1:length(R0_chfe)){
      R0_chfe[i]<-get_Ro(muPq = meso.rate.tanks[i]-mu_P, phi_Nq = phi_fert[i])
    }
    
  #Mesocosm +Atrazine+Chlorpyrifos R0 run
    R0_atch.2<-rep(0,10000)
    for(i in 1:length(R0_atch)){
      R0_atch[i]<-get_Ro(muPq = meso.rate.tanks[i]-mu_P, phi_Nq = phi_atra[i])
    }
    
  #Mesocosm +Fertilizer+Atrazine+Chlorpyrifos R0 run
    R0_tre.2<-rep(0,10000)
    for(i in 1:length(R0_tre)){
      R0_tre[i]<-get_Ro(muPq = meso.rate.tanks[i]-mu_P, phi_Nq = phi_atfe[i])
    }
    

#OBSOLETE - Plot results ####################
  r0s.2<-data.frame("Treatment"=c('Fe', 'At', 'Ch', 
                                  'At:Fe', 'Ch:Fe', 'At:Ch', 'At:Ch:Fe'),
                    "meanR0"=c(mean(R0_fert), mean(R0_atra), mean(R0_chlor),
                               mean(R0_atfe), mean(R0_chfe), mean(R0_atch), mean(R0_tre)),
                    "sd.R0"=c(sd(R0_fert), sd(R0_atra), sd(R0_chlor),
                                sd(R0_atfe), sd(R0_chfe), sd(R0_atch), sd(R0_tre)))
    
  r0s.2$Treatment<-factor(r0s.2$Treatment, levels = c('Fe', 'At', 'Ch', 
                                                      'At:Fe', 'Ch:Fe', 'At:Ch', 'At:Ch:Fe'))
  
    ggplot(r0s.2, aes(x=Treatment, y=meanR0))+
      theme_bw()+
      theme(axis.title=element_text(size=17),
            axis.text=element_text(size=12))+
      xlab("Agrochemicals Present")+
      geom_hline(aes(yintercept=R0_vil), colour='red')+
      geom_hline(aes(yintercept=R0_up), colour='red', linetype=2)+
      geom_hline(aes(yintercept=R0_lo), colour='red', linetype=2)+
      #scale_alpha_manual(values=c('green2', 'gold1', 'red', 'orange', 'darkgreen', 'pink', 'blue')) +
      ylab(expression('R'[0]))+
      geom_point(size=3) +
      geom_errorbar(aes(ymin=meanR0-(1.96*sd.R0),
                        ymax=meanR0+(1.96*sd.R0)),
                    width=.1, position=position_dodge(.7))
      
    
    
    
#OBSOLETE - Another option for propagating uncertainty: just use results from individual tanks ##########################
Ch.mean.Bt=mean(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 &dat$fert==0]) #mean observed Bt. fin when only chlorP is present
  fert.vec=c(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 &dat$fert==1])/Ch.mean.Bt #estimated scalar of carrying capacity when Fe=1
  atra.vec=c(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 &dat$fert==0])/Ch.mean.Bt #estimated scalar of carrying capacity when At=1
  atfe.vec=c(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 &dat$fert==1])/Ch.mean.Bt #estimated scalar of carrying capacity when AtFe=1
  chlor.vec=c(-log(dat$p.all_24[dat$chlor==1 & dat$atra==0 &dat$fert==0]/3)) 
    chlor.vec[chlor.vec==Inf]<-1.0986123 #Replace infinite values with high prawn mortality rate that results in total prawn mortality anyways
      chlor.vec=chlor.vec-mu_P #MAKE SURE mu_P=0.03...
  chfe.vec=c(-log(dat$p.all_24[dat$chlor==1 & dat$atra==0 &dat$fert==1]/3)) 
    chfe.vec=chfe.vec-mu_P
  atch.vec=c(-log(dat$p.all_24[dat$chlor==1 & dat$atra==1 &dat$fert==0]/3)) 
    atch.vec[atch.vec==Inf]<-1.0986123
      atch.vec=atch.vec-mu_P
  tre.vec=c(-log(dat$p.all_24[dat$chlor==1 & dat$atra==1 &dat$fert==1]/3)) 
    tre.vec[tre.vec==Inf]<-1.0986123
      tre.vec=tre.vec-mu_P
#R0 runs with these vectors  
  #Mesocosm +Fertilizer R0 run
    R0_fert2<-get_Ro(muPq = 0, phi_Nq = fert.vec)
  
  #Mesocosm +Atrazine R0 run
    R0_atra2<-get_Ro(muPq = 0, phi_Nq = atra.vec)

  #Mesocosm +Chlorpyrifos R0 run
    R0_chlor2<-get_Ro(muPq = chfe.vec, phi_Nq = 1)
  
  #Mesocosm +Fertilizer+Atrazine R0 run
    R0_atfe2<-get_Ro(muPq = 0, phi_Nq = atfe.vec)
  
  #Mesocosm +Fertilizer+Chlorpyrifos R0 run
    R0_chfe2<-get_Ro(muPq = chfe.vec, phi_Nq = fert.vec)
  
  #Mesocosm +Atrazine+Chlorpyrifos R0 run
    R0_atch2<-get_Ro(muPq = atch.vec, phi_Nq = atra.vec)
  
  #Mesocosm +Fertilizer+Atrazine+Chlorpyrifos R0 run
    R0_tre2<-get_Ro(muPq = tre.vec, phi_Nq = atfe.vec)
  
#Boxplot    
r0s.box2<-data.frame('At'=R0_atra2,
                        'Ch'=R0_chlor2,
                        'Fe'=R0_fert2,  
                        'At:Fe'=R0_atfe2, 
                        'Ch:Fe'=R0_chfe2, 
                        'At:Ch'=R0_atch2, 
                        'At:Ch:Fe'=R0_tre2)
    
    boxplot(r0s.box2)
      abline(h=1.068942, col='red', lty=2)    
      
      
      
#OBSOLETE - Using regression equation for estimates #########################
#Read in data that's merged with SEModel  
  ind1<-read.csv("~/RemaisWork/Schisto/Data/Halstead_etal/ind1.csv")
    lm.bt.end1<-lm(bt_liv_fin~Pred.2+Alg2.2, data=ind1) #End B. truncatus numbers; prediction function
      summary(lm.bt.end1)
      
#Use regression equation for predicted values in final bt counts 
  bt.fin.P<-predict(lm.bt.end1, data=ind1, se.fit=T)
    ind1$bt_fin.fit<-bt.fin.P$fit
    ind1$bt_fin.fit.err<-bt.fin.P$se.fit
  #Trim the data frame a bit  
  ind1.sub<-subset(ind1, select = c(Tank,At,Ch,Fe,Treats,Alg2.2,Pred.2,bt_liv_fin, p.all_24,bt_fin.fit,bt_fin.fit.err))  
    ind1.sub$Treats[ind1.sub$bt_fin.fit<0] 
    
  mean(ind1.sub$bt_fin.fit[ind1.sub$Treats=="ChlorP"]) #155.7778 - predicted bt.fin without predation or enhanced algae
    ind1.sub$phi_Nq<-ind1.sub$bt_fin.fit/155.7778
      ind1.sub$phi_Nq[ind1.sub$Ch==0]<-1 #Not possible to assess effect of algae alone when chlorP is absent
      #So replace fert/ atra only trials with mean phi_Nq from combo treatments with chlorP
      ind1.sub$phi_Nq[ind1.sub$Treats=="Fert"]<-mean(ind1.sub$phi_Nq[ind1.sub$Treats=="Chlor_Fert"])
      ind1.sub$phi_Nq[ind1.sub$Treats=="Atra"]<-mean(ind1.sub$phi_Nq[ind1.sub$Treats=="Atra_Chlor"])
      ind1.sub$phi_Nq[ind1.sub$Treats=="Atra_Fert"]<-mean(ind1.sub$phi_Nq[ind1.sub$Treats=="All_Three"])
      
    ind1.sub$prawn.surv<-ind1.sub$p.all_24/3
      ind1.sub$prawn.surv[ind1.sub$prawn.surv==0]<-0.15#95% mortality of 3 prawns=2.85 =.15 survival
    ind1.sub$muPq<- -log(ind1.sub$prawn.surv) #fill prwan mortality variable
      ind1.sub$muPq<- ind1.sub$muPq - mu_P #want excess mortality above background rate
      ind1.sub$muPq[ind1.sub$muPq<0]<-0
      
  #How well does predicted B. truncatus match observed?
    plot(x=ind1.sub$bt_fin.fit, y=ind1.sub$bt_liv_fin, xlab='Predicted',ylab='Observed', pch=16,
         ylim=c(0,500), xlim=c(0,500))
      points(ind1.sub$bt_fin.fit[ind1.sub$Ch==1], ind1.sub$bt_liv_fin[ind1.sub$Ch==1], col='red',pch=15, cex=1.5)
      points(ind1.sub$bt_fin.fit[ind1.sub$At==1], ind1.sub$bt_liv_fin[ind1.sub$At==1], col='gold2',pch=17, cex=1.4)
      points(ind1.sub$bt_fin.fit[ind1.sub$Fe==1], ind1.sub$bt_liv_fin[ind1.sub$Fe==1], col='green',pch=16)
      legend('topright', legend=c('Ch', 'At', 'Fe'), pch=c(15,17,16), col=c('red','gold2','green'), title='agroCs present')
  
  meso.muPqs<-ind1.sub$muPq
  meso.phiNqs<-ind1.sub$phi_Nq
  
  meso.R0s<-rep(0, 60)
  for(i in 1:60){
    meso.R0s[i]=get_Ro(muPq = meso.muPqs[i], phi_Nq = meso.phiNqs[i])
  }
  
  ind1.sub$R0<-meso.R0s
  
    r0s.meso<-data.frame("At"=c(ind1.sub$R0[ind1.sub$Treats=="Atra"]),
                         "Ch"=c(ind1.sub$R0[ind1.sub$Treats=="ChlorP"]),
                         "Fe"=c(ind1.sub$R0[ind1.sub$Treats=="Fert"]),
                         "At:Ch"=c(ind1.sub$R0[ind1.sub$Treats=="Atra_Chlor"]),
                         "At:Fe"=c(ind1.sub$R0[ind1.sub$Treats=="Atra_Fert"]),
                         "Ch:Fe"=c(ind1.sub$R0[ind1.sub$Treats=="Chlor_Fert"]),
                         "At:Ch:Fe"=c(ind1.sub$R0[ind1.sub$Treats=="All_Three"]))
          
    boxplot(r0s.meso, ylab=expression('R'[0]), xlab='Agrochemicals Present')
      abline(h=R0_vil, lty=2, col='red')
      
      
#Use uncertainty in beta as our R0 uncertainty, fixed values for all else ################
  get_Ro<-function(muPq, phi_Nq, beta) #added beta to the call list so easier to change
      {   #HAVE TO SET muPq and phi_Nq in function call
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
        lamda<-parameters["lamda"]
        mu_W<-parameters["mu_W"]
        m<-parameters["m"]
        H<-parameters["H"]
        mu_H<-parameters["mu_H"]
        
        N_eq<-(1-(mu_N/f_N))/(phi_N*(1/phi_Nq)) #Equilibrium estimate of N given snail parameters
        P_eq<-(1-((muPq+mu_P)/f_P))/phi_P #Equilibrium estimate of P given prawn predator parameters
        if(P_eq<0){
          P_eq=0
        }
        pred<-(alpha*P_eq)/(1+(alpha*N_eq*Th))#death rate of snails due to predators given equilibrium estimates of P and N
        
        T1<-0.5*beta*m*H*N_eq
        T2<-lamda*sigma
        T3<- (mu_W+mu_H)*(mu_N+pred+sigma)*(mu_N+pred+mu_I)
        
        Ro_est <- sqrt((T1*T2)/T3)
        
        print(N_eq)
        print(P_eq)
        Ro_est 
        
  }
    
    parameters["phi_P"]<-1/(40*3)    
    beta_0=1.1128e-6
    beta_up=1.5484e-06
    beta_lo=4e-07#beta value corresponding to R0 of 1.01; previous value was 4e-7

#village R0        
  R0_vil=1.068942 #1.068942
  R0_up=1.260919 #1.260919
  R0_lo=0.6408785 #0.6408785
  
#village R0 with prawns
  parameters["phi_P"]<-1/(40*3) 
    R0_vil.p<-get_Ro(muPq = 0, phi_Nq = 1, beta=beta_0)
    R0_up.p<-get_Ro(muPq = 0, phi_Nq = 1, beta=beta_up)
    R0_lo.p<-get_Ro(muPq = 0, phi_Nq = 1, beta=beta_lo)

    
#atrazine only R0
  R0_atra0 = get_Ro(muPq = 0, phi_Nq = 1.614304, beta = beta_0)
  R0_atra_up = get_Ro(muPq = 0, phi_Nq = 1.614304, beta = beta_up)
  R0_atra_lo = get_Ro(muPq = 0, phi_Nq = 1.614304, beta = beta_lo)
  
#chlorpyrifos only R0
  R0_chlor0 = get_Ro(muPq = mesorate.chlorPres, phi_Nq = 1, beta = beta_0)
  R0_chlor_up = get_Ro(muPq = mesorate.chlorPres, phi_Nq = 1, beta = beta_up)
  R0_chlor_lo = get_Ro(muPq = mesorate.chlorPres, phi_Nq = 1, beta = beta_lo)
  
#fertilizer only R0
  R0_fert0 = get_Ro(muPq = 0, phi_Nq = 1.159642, beta = beta_0)
  R0_fert_up = get_Ro(muPq = 0, phi_Nq = 1.159642, beta = beta_up)
  R0_fert_lo = get_Ro(muPq = 0, phi_Nq = 1.159642, beta = beta_lo)  
  
#atrazine+chlorpyrifos R0
  R0_atch0 = get_Ro(muPq = mesorate.chlorPres, phi_Nq = 1.614304, beta = beta_0)
  R0_atch_up = get_Ro(muPq = mesorate.chlorPres, phi_Nq = 1.614304, beta = beta_up)
  R0_atch_lo = get_Ro(muPq = mesorate.chlorPres, phi_Nq = 1.614304, beta = beta_lo) 
  
#atrazine+fertilizer R0
  R0_atfe0 = get_Ro(muPq = 0, phi_Nq = 1.509579, beta = beta_0)
  R0_atfe_up = get_Ro(muPq = 0, phi_Nq = 1.509579, beta = beta_up)
  R0_atfe_lo = get_Ro(muPq = 0, phi_Nq = 1.509579, beta = beta_lo) 
  
#chlorpyrifos+fertilizer R0
  R0_chfe0 = get_Ro(muPq = mesorate.chlorPres, phi_Nq = 1.159642, beta = beta_0)
  R0_chfe_up = get_Ro(muPq = mesorate.chlorPres, phi_Nq = 1.159642, beta = beta_up)
  R0_chfe_lo = get_Ro(muPq = mesorate.chlorPres, phi_Nq = 1.159642, beta = beta_lo)  
  
#atrazine+chlorpyrifos+fertilizer R0
  R0_tre0 = get_Ro(muPq = mesorate.chlorPres, phi_Nq = 1.509579, beta = beta_0)
  R0_tre_up = get_Ro(muPq = mesorate.chlorPres, phi_Nq = 1.509579, beta = beta_up)
  R0_tre_lo = get_Ro(muPq = mesorate.chlorPres, phi_Nq = 1.509579, beta = beta_lo) 
  
r0s.3<-data.frame("Treatment"=c('At', 'Ch', 'Fe',  
                                'At:Ch', 'At:Fe', 'Ch:Fe', 'At:Ch:Fe'),
                  'r0_0'=c(R0_atra0, R0_chlor0, R0_fert0, R0_atch0, R0_atfe0, R0_chfe0, R0_tre0),
                  'r0_up'=c(R0_atra_up, R0_chlor_up, R0_fert_up, R0_atch_up, R0_atfe_up, R0_chfe_up, R0_tre_up),
                  'r0_lo'=c(R0_atra_lo,  R0_chlor_lo, R0_fert_lo, R0_atch_lo, R0_atfe_lo, R0_chfe_lo, R0_tre_lo))  

r0s.3$Treatment<-factor(r0s.3$Treatment, levels = c('Fe', 'At', 'At:Fe', #Bottom up effects only
                                                    'Ch', #Top-down effects only
                                                    'At:Ch',  'Ch:Fe',  'At:Ch:Fe')) #Both top-down and bottom-up effects

gg1<-ggplot(r0s.3, aes(x=Treatment, y=r0_0))+
  #Theme formatting
    theme_bw()+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=10))+
    scale_y_continuous(breaks=c(0, 0.44, 0.5, 1.0, 1.07, 1.5, 2.0), limits=c(0,2))+
    xlab("")+
    ylab(expression('R'[0]))+
  #Village R0 lines
    geom_hline(aes(yintercept=R0_vil), colour='grey30', size=1, linetype=3)+
    #geom_hline(aes(yintercept=R0_up), colour='red', linetype=2, alpha=0.5, size=1)+
    #geom_hline(aes(yintercept=R0_lo), colour='red', linetype=2, alpha=0.5, size=1)+
    geom_hline(aes(yintercept=R0_vil.p), colour='grey30', size=1, linetype=3)+
    #geom_hline(aes(yintercept=R0_up.p), colour='green2', linetype=2, alpha=0.5, size=1)+
    #geom_hline(aes(yintercept=R0_lo.p), colour='green2', linetype=2, alpha=0.5, size=1)+
  #Adding data from data frame
    geom_point(size=3) +
    geom_errorbar(aes(ymin=r0_lo,
                      ymax=r0_up),
                  width=.1, position=position_dodge(.7))+
  #Add labels to R0 lines
    geom_label(x=1.475, y=1.155, label="R ",  size=5, colour = 'grey30', fill='white', label.size = NA)+
      geom_text(x=1.525, y=1.125, label="0,f",  size=3, colour = 'grey30')+
    geom_label(x=5.475, y=0.53, label='R ',  size=5, colour = 'grey30', fill='white', label.size = NA)+
      geom_text(x=5.525, y=0.508, label="0,f",  size=3, colour = 'grey30')+
      geom_text(x=5.585, y=0.50, label="p",  size=3, colour = 'grey30')+
  #Add treatment labels
    geom_segment(x=0.7, xend=3.3, y=1.7, yend=1.7, colour='grey40', lineend='square')+
      geom_text(x=2, y=1.76, label='bottom-up effects', size=5, colour='grey40')+
    geom_segment(x=3.7, xend=4.3, y=1.7, yend=1.7, colour='grey40', lineend='square')+
      geom_text(x=4, y=1.84, label='top-down', size=5, colour='grey40')+
      geom_text(x=4, y=1.76, label='effects', size=5, colour='grey40')+
    geom_segment(x=4.7, xend=7.3, y=1.7, yend=1.7, colour='grey40', lineend='square')+
      geom_text(x=6, y=1.76, label='bottom-up & top-down effects', size=5, colour='grey40')+
  #Add plot label
    geom_text(label='A', x=0.57, y=2, size=10)
  

#Check out estimates of mean worm burden (W) given these paremeters ####################
  #Starting values
    p=3*39
    nstart=c(S=10000,E=0,I=0, W=6, P=p)
    parameters["mu_P"]<-0.03883984
    parameters["muPq"]<-mesorate.chlorPres
    parameters["phi_P"]<-1/(3*40)
    parameters["beta"]<-1.1128e-06
  #Run
  output.0=as.data.frame(ode(nstart,time,schisto_halstead,parameters)) 
  
  eqbm.0<-output.0[365*yrs,]
  eqbm.0
  
  mu_P=parameters["mu_P"]
  
  parameters["muPq"]<-1.321756-parameters["mu_P"]
  
  
  
#Plot R0 response across chlorP doses #####################
  parameters["muPq"]=0
  parameters["muP"]=0 #let's just make all mortality due to mortality observed in ecotox study
  p.ecotox$rate0<- -log(1-p.ecotox$mort)
    p.ecotox$rate.up<- -log(1-(p.ecotox$mort+(1.96*p.ecotox$st.er)))
      p.ecotox$rate.up[is.na(p.ecotox$rate.up)]<-max(p.ecotox$rate.up[p.ecotox$rate.up!=Inf], na.rm=T)
      p.ecotox$rate.up[p.ecotox$rate.up==Inf]<-max(p.ecotox$rate.up[p.ecotox$rate.up!=Inf], na.rm=T)
      
    p.ecotox$rate.lo<- -log(1-(p.ecotox$mort-(1.96*p.ecotox$st.er)))

  for (i in 1:nrow(p.ecotox)){
    p.ecotox$R0[i] <- get_Ro(muPq = (p.ecotox$rate0[i]), phi_Nq = 1, beta = beta_0)
  }
    
    for (i in 1:nrow(p.ecotox)){
      p.ecotox$R0_lo[i] <- get_Ro(muPq = (p.ecotox$rate.lo[i]), phi_Nq = 1, beta = beta_0)
    }
    
    for (i in 1:nrow(p.ecotox)){
      p.ecotox$R0_up[i] <- get_Ro(muPq = (p.ecotox$rate.up[i]), phi_Nq = 1, beta = beta_0)
    }
   
  gg2<-ggplot(p.ecotox, aes(x=dose, y=R0))+
    theme_bw()+
    theme(axis.text=element_text(size=10), axis.title=element_text(size=14))+
    scale_y_continuous(breaks=c(0, 0.44, 0.5,1.0,1.07,1.5,2.0), limits=c(0,2))+
    scale_x_continuous(breaks=c(0,20,40,60,64), limits=c(0,70))+
    ylab(expression('R'[0]))+
    xlab(expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep = '')))+
    geom_line()+
    geom_line(aes(y=R0_lo), linetype=2)+
    geom_line(aes(y=R0_up), linetype=2)+
    geom_text(label='B', x=0, y=2, size=10)#+
    #geom_label(x=64, y=1.15, label="concentration", colour='grey40', size=4, fill='white',label.size = NA)+
    #geom_label(x=64, y=1.225, label="tested", colour='grey40', size=4, fill='white', label.size = NA)+
    #geom_segment(x=64, xend=64, y=-Inf, yend=1.08, linetype=3, colour='grey40')
      
plot(p.ecotox$dose, p.ecotox$R0, type='l', bty='l', ylim=c(0,1.2), xlim=c(0,64),
     xlab='Chlorpyrifos dose', ylab=expression('R'[0]))
  lines(p.ecotox$dose, p.ecotox$R0_up, lty=2)
  lines(p.ecotox$dose, p.ecotox$R0_lo, lty=2)
  
  
  
#Create a heat map of R0 values across different chlorP doses and snail carrying capacity values ################
  chlorP.doses<-c(0, 0.32, 0.64, 3.2, 6.4, 32, 64)#doses used in ecotox Halstead paper
    rat.doses<-predict(ecotox1, data.frame('dose'=chlorP.doses), 
                        type = "response", se.fit=TRUE)$fit
    rate.doses<- -log(1-rat.doses)
  phiNq.doses<-seq(1,1.6,by=.1)
  phiNq.doses.percent<-seq(0,60, by=10)
  
  #Create matrix of R0 estimates given parameter values above
  r0.doses<-matrix(ncol=7, nrow=7)
  parameters["mu_P"]=0
    for (i in 1:7){
      for (j in 1:7){
        r0.doses[i,j]<-get_Ro(muPq = rate.doses[j],
                            phi_Nq = phiNq.doses[i],
                            beta = beta_0)
      }
     
    }
  
  r0.doses.df<-data.frame('R0'=c(r0.doses[,1],
                                 r0.doses[,2],
                                 r0.doses[,3],
                                 r0.doses[,4],
                                 r0.doses[,5],
                                 r0.doses[,6],
                                 r0.doses[,7]),
                          'PhiNq'=c(rep(phiNq.doses.percent, 7)),
                          'dose'=c(rep(chlorP.doses[1],7),
                                     rep(chlorP.doses[2],7),
                                     rep(chlorP.doses[3],7),
                                     rep(chlorP.doses[4],7),
                                     rep(chlorP.doses[5],7),
                                     rep(chlorP.doses[6],7),
                                     rep(chlorP.doses[7],7)))
  
  r0.doses.df$mort<-predict(ecotox1, r0.doses.df, type = "response", se.fit=TRUE)$fit
  r0.doses.df$rate<- round(-log(1-r0.doses.df$mort), digits=4)
    r0.doses.df$rate<-factor(r0.doses.df$rate, levels=c(r0.doses.df$rate[1],
                                                        r0.doses.df$rate[8],
                                                        r0.doses.df$rate[15],
                                                        r0.doses.df$rate[22],
                                                        r0.doses.df$rate[29],
                                                        r0.doses.df$rate[36],
                                                        r0.doses.df$rate[43]))
    
    r0.doses.df$dose<-factor(r0.doses.df$dose, levels=c(r0.doses.df$dose[1],
                                                        r0.doses.df$dose[8],
                                                        r0.doses.df$dose[15],
                                                        r0.doses.df$dose[22],
                                                        r0.doses.df$dose[29],
                                                        r0.doses.df$dose[36],
                                                        r0.doses.df$dose[43]))
  r0.doses.df$PhiNq<-factor(r0.doses.df$PhiNq, levels=c(phiNq.doses.percent))
    
  gg3<-ggplot(r0.doses.df, aes(x=dose, y=PhiNq, fill=R0))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='grey90', high='black')+
    coord_equal()+
    labs(x=expression(paste('Chloryrifos concentration (', mu, 'g/L)', sep = '')), 
    #(x=expression(paste('Predator mortality rate (', mu[P][,][q], ')', sep = '')), axis label for predator mortality rate
         y=expression(paste('% increase in snail carrying capacity (', phi[N][,][q], ')', sep='')))+
    theme(axis.ticks=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=14),
          legend.title=element_text(size=15), legend.text=element_text(size=12))+
    geom_text(label='C', x=0.75, y=7.25, size=10, alpha=.50)
    
  
  
#Put three plots together for final figure ###########################
  grid.arrange(gg1,gg2,gg3, ncol=2, nrow=2, layout_matrix=rbind(c(1,1),
                                                                c(2,3)))