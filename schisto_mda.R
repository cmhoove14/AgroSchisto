## A schisto model with modeled MDA interventions
library(deSolve)
#Base schisto model 
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

#W divided into treated and untreated groups
schisto_WtWu=function(t, n, parameters) { 
  with(as.list(parameters),{
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    N=S+E+I
    W=(0.8*Wt)+(0.2*Wu)
    pred= (alpha*P)/(1+alpha*N*Th)
    inf= 0.5*beta*W*H*m*S
    dSdt= f_N*((1-(phi_N*N))*(S+(z*E))) - (mu_N*S) - (pred*S) - inf
    dEdt= inf - ((mu_N+pred+sigma)*E)
    dIdt= (sigma*E) - ((mu_N+mu_I)*I) - (pred*I)
    dWtdt= (lamda*I) - ((mu_W+mu_H)*Wt)
    dWudt= (lamda*I) - ((mu_W+mu_H)*Wu)
    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt)))
  }) 
} 

#Paramaters and initial conditions ########################
time=seq(0,365*30,1)
p=0
parameters=c(f_N=0.16, #Instantaneous snail mortality rate
             phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
             z=0.5, #Fraction of exposed snails that reproduce
             mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
             alpha=0.003, #per capita prawn attack rate on humans
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
nstart=c(S=2435,E=4444,I=1422, W=59) #VALUES USED HERE ARE PRAWN-FREE EQUILIBRIUM
output=as.data.frame(ode(nstart,time,schisto,parameters)) 
output[c((max(time)-100):max(time)),2:5]
eqbm<-output[max(time),2:5]
eqbm

#MDA Model ################
MDA_5years<-function(years,coverage,efficacy, eqbm, total_years){
  eqbm<-c(S=2435,E=4445,I=1422,W=59.3)
  entryValue<-eqbm
  exitValue<-numeric()
  mda_output<-numeric() # stores all the variable values
  W_treated<-numeric()
  W_untreated<-numeric()
  start<-numeric()
  end<-numeric()
  for(y in years){
    if(y==1){
      start<-365*(y-1)
      end<-(365*y)-1
      
      time<-seq(start-start,end-start,1)
      p=0
      parameters=c(f_N=0.16, #Instantaneous snail mortality rate
                   phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
                   z=0.5, #Fraction of exposed snails that reproduce
                   mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
                   alpha=0.003, #per capita prawn attack rate on humans
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
      nstart=as.numeric(entryValue)
      output=as.data.frame(ode(nstart,time,schisto,parameters)) 
      output_correctTime<-output
      output_correctTime[,1]<-output_correctTime[,1]+start
      mda_output<-rbind(mda_output, output_correctTime)
      mda_output<-cbind(mda_output, mda_output[,5])
      colnames(mda_output)[6]<-5
      exitValue<-as.numeric(output[max(time),2:5])
      entryValue<-exitValue
      entryValue<-c(exitValue, exitValue[4])
      entryValue[4]<-(1-efficacy)*entryValue[4]
      entryValue[5]<-entryValue[5]
    }
    if(y>1){
      start<-365*(y-1)
      end<-(365*y)-1
      
      time<-seq(start-start,end-start,1)
      p=0
      parameters=c(f_N=0.16, #Instantaneous snail mortality rate
                   phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
                   z=0.5, #Fraction of exposed snails that reproduce
                   mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
                   alpha=0.003, #per capita prawn attack rate on humans
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
      nstart=as.numeric(entryValue)
      output=as.data.frame(ode(nstart,time,schisto_WtWu,parameters)) 
      output_correctTime<-output
      output_correctTime[,1]<-output_correctTime[,1]+start
      mda_output<-rbind(mda_output, output_correctTime)
      exitValue<-as.numeric(output[max(time),2:6])
      entryValue<-exitValue
      entryValue[4]<-(1-efficacy)*(entryValue[4])
      entryValue[5]<-entryValue[5]
    }
  
    print(c(y, entryValue[4], entryValue[5]))
   
   
  }
  #run till end 
  start<-end+1
  end<-365*total_years
  time<-seq(start-start,end-start,1)
  p=0
  parameters=c(f_N=0.16, #Instantaneous snail mortality rate
               phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
               z=0.5, #Fraction of exposed snails that reproduce
               mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
               alpha=0.003, #per capita prawn attack rate on humans
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
  nstart=as.numeric(entryValue)
  output=as.data.frame(ode(nstart,time,schisto_WtWu,parameters)) 
  output_correctTime<-output
  output_correctTime[,1]<-output_correctTime[,1]+start
  mda_output<-rbind(mda_output, output_correctTime)
  
 mda_output 
 
  
  
}# end of function

#Values for MDA model  ####################
years<-c(1,2,3,4,5) # years of MDA intervention
coverage<-0.8 # coverage within human population at intervention
efficacy<-0.99 # efficacy of intervetnion (% adult worms cleared)
total_years<-30 # Years of model run

mda_output<-MDA_5years(years,coverage,efficacy, eqbm, total_years)
#Plot results of MDA model ######################
par(mfrow=c(2,1), cex.lab=1.4, cex.axis=1.3)
#Mean worm burden
par(mai=c(0.2, 1.1, 0.7,0.4))
  plot(mda_output[,1], mda_output[,5], type='l', xlab="", lwd=3,
      ylab="",ylim=c(0,max(mda_output[,5])), col="blue", xaxt='n',
      main="MDA @ years 1-5")
    mtext("Mean worm burden (W)", side=2, line=2.5, cex=1.3)
    lines(mda_output[,1], mda_output[,6], type='l', col="blue", lty=2, lwd=2)
    lines(mda_output[,1], (0.2*mda_output[,6]+ 0.8*mda_output[,5]), type='l', 
          col="blue", lty=3, lwd=2)
    legend("bottomright", c("Treated pop", "Untreated pop", "Mean pop"), 
           lty=c(1,2,3), col="blue", cex=1.1)
#prevalence
  par(mai=c(1.2, 1.1, 0.1,0.4))  
  plot(mda_output[,1], prev(W=mda_output[,5]), type='l', xlab="time (days)", lwd=3,
      ylab="",ylim=c(0,1), col="blue")
    mtext("Prevalence", side=2, line=2.5, cex=1.3)
    lines(mda_output[,1], prev(W=mda_output[,6]), type='l', col="blue", lty=2, lwd=2)
    lines(mda_output[,1], prev(W=(0.2*mda_output[,6]+ 0.8*mda_output[,5])), 
          type='l', col="blue", lty=3, lwd=2)
    
  
  
  
#MDA intervention coupled with prawn stocking #########################
Prawns_MDA_2_4<-function(Prawns=50,coverage,efficacy, eqbm, total_years){
  eqbm<-c(S=2435,E=4445,I=1422,Wt=59.3,Wu=59.3)
  entryValue<-eqbm
  exitValue<-numeric()
  prawn_output<-numeric() # stores all the variable values
  W_treated<-numeric()
  W_untreated<-numeric()
  start<-numeric()
  end<-numeric()
    
    
    #first two years
    start<-0
    end<-(365*2)-1
    
    time<-seq(start-start,end-start,1)
    parameters=c(f_N=0.16, #Instantaneous snail mortality rate
                 phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
                 z=0.5, #Fraction of exposed snails that reproduce
                 mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
                 alpha=0.003, #per capita prawn attack rate on humans
                 P=Prawns, #Number of prawns in the system
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
    nstart=as.numeric(entryValue)
    output=as.data.frame(ode(nstart,time,schisto_WtWu,parameters)) 
    output_correctTime<-output
    output_correctTime[,1]<-output_correctTime[,1]+start
    prawn_output<-rbind(prawn_output, output_correctTime)
    exitValue<-as.numeric(output[max(time),2:5])
    entryValue<-exitValue
    entryValue<-c(exitValue, exitValue[4])
    entryValue[4]<-(1-efficacy)*entryValue[4]
    entryValue[5]<-entryValue[5]
    
    #first treatment at year=2
    start<-end+1
    end<-(365*4)-1
    
    time<-seq(start-start,end-start,1)
    parameters=c(f_N=0.16, #Instantaneous snail mortality rate
                 phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
                 z=0.5, #Fraction of exposed snails that reproduce
                 mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
                 alpha=0.003, #per capita prawn attack rate on humans
                 P=Prawns, #Number of prawns in the system
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
    nstart=as.numeric(entryValue)
    output=as.data.frame(ode(nstart,time,schisto_WtWu,parameters)) 
    output_correctTime<-output
    output_correctTime[,1]<-output_correctTime[,1]+start
    prawn_output<-rbind(prawn_output, output_correctTime)
    exitValue<-as.numeric(output[max(time),2:6])
    entryValue<-exitValue
    entryValue[4]<-(1-efficacy)*(entryValue[4])
    entryValue[5]<-entryValue[5]
    
    
    print(c(2, entryValue[4], entryValue[5]))
    
    # at year 4, MDA
    start<-end+1
    end<-(365*total_years)-1
    
    time<-seq(start-start,end-start,1)
    parameters=c(f_N=0.16, #Instantaneous snail mortality rate
                 phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
                 z=0.5, #Fraction of exposed snails that reproduce
                 mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
                 alpha=0.003, #per capita prawn attack rate on humans
                 P=Prawns, #Number of prawns in the system
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
    nstart=as.numeric(entryValue)
    output=as.data.frame(ode(nstart,time,schisto_WtWu,parameters)) 
    output_correctTime<-output
    output_correctTime[,1]<-output_correctTime[,1]+start
    prawn_output<-rbind(prawn_output, output_correctTime)
    exitValue<-as.numeric(output[max(time),2:6])
    entryValue<-exitValue
    entryValue[4]<-(1-efficacy)*(entryValue[4])
    entryValue[5]<-entryValue[5]
    
    print(c(4, entryValue[4], entryValue[5]))
    
    
    
    prawn_output 
    
    
    
  }# end of function
  
#Starting values and run of MDA+Prawn model ####################
Prawns=50 #Number of prawns in system (constant)
coverage<-0.8 #Coverage of MDA intervention
efficacy<-0.99 #Efficacy of MDA intervention
total_years<-30 #Time of model run
  
  prawn_output<-Prawns_MDA_2_4(Prawns,coverage,efficacy, eqbm, total_years)
  
#Plot results of MDA + Prawn model #############################  
#Mean worm burden
par(mfrow=c(2,1), cex.axis=1.3, cex.lab=1.4)
  par(mai=c(0.2, 1.1, 0.7,0.4))
    plot(prawn_output[,1], prawn_output[,5], type='l', lwd=3,
        xlab="",ylab="Mean worm burden",main="MDA at 2&4 yrs with P=50",
        ylim=c(0,max(prawn_output[,5])), col="blue", xaxt='n')
      lines(prawn_output[,1], prawn_output[,6], type='l', col="blue", lty=2, lwd=3)
      lines(prawn_output[,1], (0.2*prawn_output[,6]+ 0.8*prawn_output[,5]), type='l', 
            col="blue", lty=3, lwd=2)
      legend("topright", c("Treated pop", "Untreated pop", "Mean pop"), lty=c(1,2,3), col="blue")
#Prevalence  
  par(mai=c(1.2, 1.1, 0.1,0.4))
    plot(prawn_output[,1], prev(W=prawn_output[,5]), type='l', lwd=3,
        xlab="time (days)",ylab="Prevalence",main="",
        ylim=c(0,1), col="blue")
      lines(prawn_output[,1], prev(W=prawn_output[,6]), type='l', col="blue", lty=2, lwd=3)
      lines(prawn_output[,1], prev(W=(0.2*prawn_output[,6]+ 0.8*prawn_output[,5])), type='l', 
            col="blue", lty=3, lwd=2)
      
#Plot results of MDA alone and MDA + Prawn models together ###################
  par(mfrow=c(2,1))
  #Plot1: MDA Alone
  par(mar=c(1.5, 4.0, 2.0, 2.0))
    plot(mda_output[,1], mda_output[,5], type='l', xlab="", xaxt='n', lwd=3,
         ylab="",ylim=c(0,max(mda_output[,5])), col="blue")
    lines(mda_output[,1], mda_output[,6], type='l', col="blue", lty=2, lwd=2)

  #Plot2: MDA+Prawns
  par(mar=c(4.5, 4.0, 1.0, 2.0), cex.lab=1.2)
    plot(prawn_output[,1], prawn_output[,5], type='l', lwd=3,
         xlab="", ylab="", ylim=c(0,max(prawn_output[,5])), col="blue",
         main="MDA with Prawns")
    mtext("time (days)", side=1, line=2.5, cex=1.3)
    mtext("Mean Worm Burden (M)", side=2, at=c(-2000,75), cex=1.5, line=2.5)
    lines(prawn_output[,1], prawn_output[,6], type='l', col="blue", lty=2, lwd=3)
    legend("topright", c("Treated pop", "Untreated pop"), lty=c(1,2), col="blue")
    