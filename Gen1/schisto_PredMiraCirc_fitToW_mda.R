# This code validates the mesocosm experiment by Halstead et al.
##### Agrochem impacts the following parameters:fNq, phi_Nq, mu_Nq, alpha_q, mu_Pq, Theta_q, R_q, V_q.
require(deSolve)

######################################################################
#### This section gives the basic functions ##########################
######################################################################


#### The basic schisto function with mda where there are 2 populations, treated and untreated. #####################################
#### start of with Wt=Wu=W (eqbm), number of treted and untreated worms are the same. 
#### Assign worms as a weighted measure of the two based on coverage.
schisto_PredMiraCirc_fitToW_2pop_mda=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    P=n[6]
    N=S+E+I
    
    #Agrochemical impacted parameters
    
    #1. Rate of successful hatching of schistosome eggs
    V<-function(V_q){
      V_o * V_q
    }
    
  
    
    #3. Rate of shedding circariae by infected snails
    Theta<-function(Theta_q){
      Theta_o * Theta_q
    }
    
    #4. Attack rate of prawns
    alpha<-function(alpha_q){
      alpha_o*alpha_q
    }
    
    #Miracidia population
    W=(cov*Wt) + ((1-cov)*Wu) #weighting treated and untreated populations
    M=0.5*W*H*m*V(V_q)
    
    
    #Circariae population
    C=I*Theta(Theta_q)
    
    #Rate of predation
    pred= (alpha(alpha_q)*P)/(1+(alpha(alpha_q)*N*Th)) #death rate of snails due to predators (Prawns)
    
    
    #snail compartment model
    dSdt= (f_N*f_Nq)*((1-((phi_N*phi_Nq)*N))*(S+(z*E))) - ( ((mu_N+mu_Nq)*S) - (pred*S) - (beta*(M)*S) )
    dEdt= beta*(M)*S - ((mu_N+mu_Nq+pred+sigma)*E)
    dIdt= (sigma*E) - ((mu_N+mu_Nq+pred+mu_I)*I)
    
    #worm burden in human
    dWtdt= (lamda*C) - ((mu_W+mu_H)*Wt) - (eff*Wt*mda )
    dWudt= (lamda*C) - ((mu_W+mu_H)*Wu)
    
    #prawn population
    dPdt= (f_P*(1-(phi_P*P))*P)-((mu_P+mu_Pq)*P)
    
    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt, dPdt )))
    }) 
} 

#### The original parameters derive from Sokolow et al ###############

CC_P<-1#prawn carrying capacity

parameters_2pops_mda=c(
  ##standard snail parameters 
  f_N=0.16,
  mu_N=1/80,
  phi_N=(1-(1/80*0.16))/10000, 
  z=0.5, 
  beta=0.000004,
  sigma=0.02, 
  mu_I=0.05, #additional snail death due to infection
  
  ## snail parameters impacted by agrochemicals
  
  f_Nq=1,
  phi_Nq=1,
  mu_Nq=0,
  
  #agrochemical concentration
  #q=0,
  
  #prawn parameters
  
  alpha_o=0.003, #attack rate
  Th=0.1,
  f_P=0.128,#prawn birth rate
  phi_P=1/CC_P,  #carrying capacity
  mu_P= 0.00137, #prawn mortality
  
  #prawn parameter impacted by agrochemicals
  mu_Pq=0, #can be a function or fixed (0.02 for eg)
  alpha_q=1,
  
  #Adult Worm, Miracidia and Circariae Parameters
  lamda=1, #fraction of circariae that infect humans
  mu_W=1/(3*365), # death rate of worm
  m=0.5,
  V_o=1,
  
  Theta_o=0.00004,
  #parameters impacted by agrochemicals
  V_q=1,
 
  Theta_q=1,
  
  #Human parameters
  H=1000, #number of humans
  mu_H=1/(60*365),
  
  #treatment parameters
  cov=0.1, #coverage of treatment across the population
  eff=0.99, # efficiency of the drug
  mda=0 # flag to indicate if mda is applied or not
)

#end of parameters



Senegal_mda<-function(nstart, parameters_2pops_mda, mda_days ){
 
   #mda function indicating when PZQ is applied
  mda<-function(t){
    #ifelse( t==(28*7) | t==(31*7), 1, 0)
    ifelse( t %in% mda_days, 1, 0)
  }
  
  #### first run to equilibrium for 7 months
  output_all<-numeric()
  
  t_eqbmtoPZQ1<-28*7 # 28 weeks or 7 months
  time=seq(0,t_eqbmtoPZQ1,1)
  
  output1=as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,parameters_2pops_mda)) 
  output_7months<-output1[dim(output1)[1],]
  time_last<-output1$time[dim(output1)[1]]
  output_all<-rbind(output_all, output1)
  
  ##### apply PZQ using mda flag and then run for 1 day
  parameters_2pops_mda["mda"]<-mda(time_last+1)
  time<-seq(time_last, time_last+1,1)
  nstart<-c(S=output_7months$S, E=output_7months$E,I=output_7months$I, Wt=output_7months$Wt, Wu=output_7months$Wu, P=p)
  output2=as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,parameters_2pops_mda)) 
  output_PZQ1<-output2[dim(output2)[1],]
  output_all<-rbind(output_all, output2[-1,])
  time_last<-output2[dim(output2)[1],1]
  
  #### apply no PZQ using mda flag and run for 3 weeks
  parameters_2pops_mda["mda"]<-mda(time_last+1)
  time<-seq(time_last, time_last+(3*7),1)
  nstart<-c(S=output_PZQ1$S, E=output_PZQ1$E,I=output_PZQ1$I, Wt=output_PZQ1$Wt, Wu=output_PZQ1$Wu, P=p)
  output3=as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,parameters_2pops_mda)) 
  output_3wks<-output3[dim(output3)[1],]
  output_all<-rbind(output_all, output3[-1,])
  time_last<-output3[dim(output3)[1],1]
  
  #apply PZQ using mda flag and run for 1 day
  parameters_2pops_mda["mda"]<-mda(time_last)
  time<-seq(time_last, time_last+1,1)
  nstart<-c(S=output_3wks$S, E=output_3wks$E,I=output_3wks$I, Wt=output_3wks$Wt, Wu=output_3wks$Wu, P=p)
  output4=as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,parameters_2pops_mda)) 
  output_PZQ2<-output4[dim(output4)[1],]
  output_all<-rbind(output_all, output4[-1,])
  time_last<-output4[dim(output4)[1],1]
  
  
  #apply no PZQ using mda flag and run for 5 months (20 weeks)
  parameters_2pops_mda["mda"]<-mda(time_last)
  time<-seq(time_last, time_last+(20*7),1)
  nstart<-c(S=output_PZQ2$S, E=output_PZQ2$E,I=output_PZQ2$I, Wt=output_PZQ2$Wt, Wu=output_PZQ2$Wu, P=p)
  output5=as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,parameters_2pops_mda)) 
  output_20wks<-output5[dim(output5)[1],]
  output_all<-rbind(output_all, output5[-1,])
  time_last<-output5[dim(output5)[1],1]
  
  output_all
  
}

mda_days<-c((28*7+1),(31*7+1)) #end of 28 weeks and 31 weeks
p=1
nstart=c(S=3892,E=3750,I=1200, Wt=20, Wu=20, P=p) 
output<-Senegal_mda(nstart, parameters_2pops_mda, mda_days)
quartz()
par(mfrow=c(2,1))
plot(output$time, output$Wt, type='l', lwd=2, col=1, ylim=c(0, max(output$Wu)))
lines(output$time, output$Wu, lwd=2, col=2)
legend(x=300, y=20, legend=c("Wt", "Wu"), col=c(1,2), lwd=2)
cov<-parameters_2pops_mda["cov"]
W<-(cov*output$Wt) + ((1-cov)*output$Wu)
plot(output$time,W ,  type='l', lwd=2, col=3, ylim=c(0, max(W)))


get_Ro_mesocosm<-function(parameters, eqbm)
{
  f_N<-parameters["f_N"]
  phi_N<-parameters["phi_N"]
  z<-parameters["z"]
  mu_N<-parameters["mu_N"]
  beta<-parameters["beta"]
  sigma<-parameters["sigma"]
  mu_I<-parameters["mu_I"]
  f_Nq<-parameters["f_Nq"]
  phi_Nq<-parameters["phi_Nq"]
  mu_Nq<-parameters["mu_Nq"]
  #q<-parameters["q"]
  alpha_o<-parameters["alpha_o"]
  Th<-parameters["Th"]
  f_P<-parameters["f_P"]
  phi_P<-parameters["phi_P"]
  mu_P<-parameters["mu_P"]
  muPq<-parameters["mu_Pq"]
  alpha_q<-parameters["alpha_q"]
  lamda<-parameters["lamda"]
  mu_W<-parameters["mu_W"]
  m<-parameters["m"]
  V_o<-parameters["V_o"]
  
  Theta_o<-parameters["Theta_o"]
  V_q<-parameters["V_q"]
  
  Theta_q<-parameters["Theta_q"]
  H<-parameters["H"]
  mu_H<-parameters["mu_H"]
  
  
  
  #1. Rate of successful hatching of schistosome eggs
  V<-function(V_q){
    V_o * V_q
  }
  
  
  
  #3. Rate of shedding circariae by infected snails
  Theta<-function(Theta_q){
    Theta_o * Theta_q
  }
  
  #4. Attack rate of prawns
  alpha<-function(alpha_q){
    alpha_o*alpha_q
  }
  N<-eqbm["S"]+eqbm["E"]+eqbm["I"]
  P<-eqbm["P"]
  pred<- as.numeric((alpha(alpha_q)*P)/(1+(alpha(alpha_q)*N*Th)) )#death rate of snails due to predators (Prawns)
  
  #N_star<-(1-((mu_N + mu_Nq)/(f_N*f_Nq)))/(phi_N+phi_Nq)
  T1<-beta*0.5*H*m*V(V_q)*N
  T2<- sigma*lamda*Theta(Theta_q)
  T3<-(mu_W+mu_H)*( mu_N + mu_Nq+ sigma+pred)*(mu_N+mu_Nq+mu_I+pred)
  
  Ro_miraCirc <- sqrt( (T1*T2)/T3 ) 
  
  
  Ro_miraCirc 
  
}

#### The function for prevalence #####################################

k<-0.25 # clumping parameter of the negative binomial distribution
Prevalence <- function(W, k) {
  p=1 - (1/(1+W/k)^(k))*(1+W/(1+W/k)) # fraction of humans with at least 2 parasites
  return(p)
}


######################################################################
#### Fit the beta of the villages to their worm burden W at baseline and after PZQ #############
#### Thr baseline condition is no prawns, no agro ####################
######################################################################

# Villages are - Lampsar I, Lampsar II, Mbarigo I, Mbarigo II, Ndiol Maure, Pokhotane, Nder
W_baseline<-c(25, 6.5, 17.9, 35, 29.2, 599, 245)
W_5month_field_mean<-c(1.12,1.5,6.2,12.1,8.5,408,355)
W_5month_field_k<-c(0.02,0.02,0.18,0.23,0.04,0.17,0.78)
W_5month_field_sd<-sqrt( (W_5month_field_mean)+(W_5month_field_k*((W_5month_field_mean)^2) ) )
N<-129
W_5month_field_se<-W_5month_field_sd/sqrt(N)
timepoints<-c(0, 28*7, 51*7)
W_bestBeta<-rep(0, times=7)
#run the test only for villages 1 to 5 that have S. Mansoni
for(b in 2:5){
  
  #### No prawns #####################################################
  CC_P<-1 
  #### CC_N=10000 #####################################################
  phi_N_new<-(1-(1/80*0.16))/10000 
  #### generate a range of beta ########################################
  min_beta<-parameters_2pops_mda["beta"]/10
  max_beta<-parameters_2pops_mda["beta"]*10
  beta_range<-seq(from=min_beta, to=max_beta, by=(max_beta-min_beta)/1000  ) 
  nstart=c(S=3892,E=3750,I=1200, Wt=20, Wu=20, P=p) 
  W_baseline_range<-rep(0, times=length(beta_range))
  W_postmda_range<-rep(0, times=length(beta_range))
  sqErr_range<-rep(0, times=length(beta_range))
  
  for(i in 1:length(beta_range)){
    params<-parameters_2pops_mda
    params["phi_P"]<-1/CC_P
    params["beta"]<-beta_range[i]
    params["phi_N"]<-phi_N_new
    params["mda"]<-0 #no mda
    #### run to equilibrium
    time<-seq(from=0, to=30*365, by=1)
    output<-as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,params))
    output_eqbm<-output[dim(output)[1],]
    cov<-params["cov"]
    W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
    W_baseline_range[i]<-W_eqbm
    
    #apply mda
    nstart1<-c(S=output_eqbm$S,E=output_eqbm$E,I=output_eqbm$I, Wt=output_eqbm$Wt, Wu=output_eqbm$Wu, P=p)                   
    output<-Senegal_mda(nstart1, params, mda_days) 
    output_eqbm<-output[dim(output)[1],]
    cov<-params["cov"]
    W_postmda<-as.numeric((cov*output_eqbm$Wt ) + ((1-cov)*output_eqbm$Wu ) )
    #W_postmda_range[i]<-W_postmda
    W_postmda_range[i]<-output_eqbm$Wt
    
    sqErr<-(W_baseline_range[i]-W_baseline[b])^2 + (W_postmda_range[i]-W_5month_field_mean[b])^2
    sqErr_range[i]<-sqErr
    print(c(b,i))
  }#end of i loop
  
  
  bestBetaIndex<-which(sqErr_range==min(sqErr_range) )
  W_bestBeta[b]<-beta_range[bestBetaIndex]
  
  ### lets get all the ouput for that beta
  params<-parameters_2pops_mda
  params["phi_P"]<-1/CC_P
  params["beta"]<-beta_range[bestBetaIndex]
  params["phi_N"]<-phi_N_new
  params["mda"]<-0 #no mda
  #### run to equilibrium
  time<-seq(from=0, to=30*365, by=1)
  output<-as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,params))
  output_eqbm<-output[dim(output)[1],]
  cov<-params["cov"]
  W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )

  #apply mda
  nstart1<-c(S=output_eqbm$S,E=output_eqbm$E,I=output_eqbm$I, Wt=output_eqbm$Wt, Wu=output_eqbm$Wu, P=p)                   
  output<-Senegal_mda(nstart1, params, mda_days) 
  output_eqbm<-output[dim(output)[1],]
  cov<-params["cov"]
  W_postmda<-as.numeric((cov*output_eqbm$Wt ) + ((1-cov)*output_eqbm$Wu ) )
  W_postmda<-output_eqbm$Wt
  #lets plot and see
  datapoints<-c(W_baseline[b], W_baseline[b], W_5month_field_mean[b])
  W<-as.numeric((cov*output$Wt ) + ((1-cov)*output$Wu ) )
  W<-output$Wt
  quartz()
  plot(output$time, W , type="l", col=b, lwd=2, ylim=c(0,max(c(W_baseline[b],W) )), ylab="W", xlab="time" )
  points(timepoints, c(W_baseline[b], W_baseline[b], W_5month_field_mean[b]), pch=16 )
  segments(timepoints[3], datapoints[3]-(1.96*W_5month_field_se[b]), timepoints[3], datapoints[3]+(1.96*W_5month_field_se[b]), lwd=2)
  title(paste("Village #", b, "with best fit beta =", beta_range[bestBetaIndex], sep=" "))
 }


##### save the data
# outputfile<-paste("/Users/Arathi/Documents/2015/Monash /Data/SciencePaper_modelData/W_bestBeta", ".Rdata", sep="")
# #save(W_bestBeta, file=outputfile )
# load(outputfile)


######################################################################
#### Baseline Ro #####################################################
######################################################################
beta_villages<-W_bestBeta
CC_N_villages<-c(10000, 10000, 10000, 10000, 10000, 100000, 100000)
phi_N_villages<-(1-(1/80*0.16))/CC_N_villages
k_villages<-c(0.17, 0.08, 0.28, 0.23, 0.17, 0.29, 0.66) # from Sanna Sokolow
W_villages<-rep(0, times=7)
Ro_villages<-rep(0, times=7)
Prev_villages<-rep(0, times=7)

CC_P<-1
nstart=c(S=3892,E=3750,I=1200, Wt=20, Wu=20, P=p) 
p=1
time<-seq(from=0, to=30*365, by=1)
for(b in 2:5 ){
  params<-parameters_2pops_mda
  params["phi_P"]<-1/CC_P
  params["beta"]<-beta_villages[b]
  params["phi_N"]<-phi_N_villages[b]
  #run to equilibrium.
  output<-as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,params))
  eqbm<-output[dim(output)[1],]
  W_villages[b]<-eqbm$Wt
  Ro_villages[b]<-as.numeric(get_Ro_mesocosm(params, eqbm))
  Prev_villages[b]<-Prevalence(eqbm$Wt, k_villages[b])
}# end of b loop

#### Plot the baseline results
quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(2:5, W_villages[2:5], type='b', pch=16, lwd=2, xlab="villages", ylab="W_baseline")
plot(2:5, Ro_villages[2:5], type='b', pch=16, lwd=2, xlab="villages", ylab="Ro_baseline")
plot(2:5, Prev_villages[2:5], type='b', pch=16, lwd=2, xlab="villages", ylab="Prev_baseline")
barCenters <-barplot(mean(Ro_villages[2:5]), names.arg = "Mean Baseline Ro", width=1, ylim=c(0,5), xlim=c(0,5) )
segments(barCenters, mean(Ro_villages[2:5]) - (sqrt(var(Ro_villages[2:5])) * 2), barCenters, mean(Ro_villages[2:5]) + (sqrt(var(Ro_villages[2:5])) * 2), lwd = 1.5)

arrows(barCenters , mean(Ro_villages[2:5]) - (sqrt(var(Ro_villages[2:5])) * 2), barCenters, mean(Ro_villages[2:5]) + (sqrt(var(Ro_villages[2:5])) * 2), lwd = 1.5, angle = 90,code = 3, length = 0.05)





################################################################################
#### Fit village 2 (Lampsar II) to all 4 data points without adding seasonality to see how it does.
W_field<-c(6.5,1.5,161,17.6)
W_DP<-rep(0, times=4) #4 data points
W_model<-rep(0, times=4)
output_model<-numeric()

params<-parameters_2pops_mda
params["beta"]<-W_bestBeta[2]
params["phi_P"]<-1/CC_P
params["phi_N"]<-phi_N_villages[2]
params["mda"]<-0 #no mda
#### run to equilibrium
time<-seq(from=0, to=30*365, by=1)
nstart=c(S=3892,E=3750,I=1200, Wt=20, Wu=20, P=p) 

output<-as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,params))
output_eqbm<-output[dim(output)[1],]
cov<-params["cov"]
W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
W_model[1]<-W_eqbm
W_DP[1]<-0

#apply PZQ1 and PZQ2 mda
nstart1<-c(S=output_eqbm$S,E=output_eqbm$E,I=output_eqbm$I, Wt=output_eqbm$Wt, Wu=output_eqbm$Wu, P=p)                   
output<-Senegal_mda(nstart1, params, mda_days) 
output_model<-rbind(output_model, output)
output_DP2<-output[dim(output)[1],] #datapoint 2
cov<-params["cov"]
#W_postmda<-as.numeric((cov*output_eqbm$Wt ) + ((1-cov)*output_eqbm$Wu ) )
#W_postmda_range[i]<-W_postmda
time_last<-output$time[dim(output)[1]]
W_model[2]<-output_DP2$Wt
W_DP[2]<-time_last


#Apply PZQ 3 mda and run for 7 months
parameters_2pops_mda["mda"]<-1
time<-seq(time_last, time_last+1,1)
nstart<-c(S=output_DP2$S, E=output_DP2$E,I=output_DP2$I, Wt=output_DP2$Wt, Wu=output_DP2$Wu, P=p)
output3=as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,parameters_2pops_mda)) 
output_PZQ3<-output3[dim(output3)[1],]
output_model<-rbind(output_model, output3[-1,])
time_last<-output3[dim(output3)[1],1]

parameters_2pops_mda["mda"]<-0
time<-seq(time_last, time_last+(7*7*4),1)
nstart<-c(S=output_PZQ3$S, E=output_PZQ3$E,I=output_PZQ3$I, Wt=output_PZQ3$Wt, Wu=output_PZQ3$Wu, P=p)
output4=as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,parameters_2pops_mda)) 
output_DP3<-output4[dim(output4)[1],]
output_model<-rbind(output_model, output4[-1,])
time_last<-output4[dim(output4)[1],1]
W_model[3]<-output_DP3$Wt
W_DP[3]<-time_last

# apply PZQ 4 mda and run for 6 months

parameters_2pops_mda["mda"]<-1
time<-seq(time_last, time_last+1,1)
nstart<-c(S=output_DP3$S, E=output_DP3$E,I=output_DP3$I, Wt=output_DP3$Wt, Wu=output_DP3$Wu, P=p)
output5=as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,parameters_2pops_mda)) 
output_PZQ4<-output5[dim(output5)[1],]
output_model<-rbind(output_model, output5[-1,])
time_last<-output5[dim(output5)[1],1]

parameters_2pops_mda["mda"]<-0
time<-seq(time_last, time_last+(6*4*7),1)
nstart<-c(S=output_PZQ4$S, E=output_PZQ4$E,I=output_PZQ4$I, Wt=output_PZQ4$Wt, Wu=output_PZQ4$Wu, P=p)
output6=as.data.frame(ode(nstart,time,schisto_PredMiraCirc_fitToW_2pop_mda,parameters_2pops_mda)) 
output_DP4<-output6[dim(output6)[1],]
output_model<-rbind(output_model, output6[-1,])
time_last<-output6[dim(output6)[1],1]

W_model[4]<-output_DP4$Wt
W_DP[4]<-time_last


quartz()
par(mfrow=c(2,1))
plot(output_model$time, output_model$Wt,type="l", col=3, lwd=2, xlab="time", ylab="Wt", ylim=c(0, max(output_model$Wt)))
points(W_DP, W_model, pch=16, col=1)
title("LampsarII model with no seasonal adjustments")
plot(output_model$time, output_model$Wt,type="l", col=3, lwd=2, xlab="time", ylab="Wt", ylim=c(0, max(W_field)))
points(W_DP, W_field, pch=16, col="dark blue")
title("Lampsar II model against field data")

quartz()
par(mfrow=c(2,1))
plot(output_model$time, output_model$Wt,type="l", col=3, lwd=2, xlab="time", ylab="Wt", ylim=c(0, max(W_field)))
points(W_DP, W_field, pch=16, col="dark blue")

plot(W_DP, W_model, pch=16, col=1, type="p", xlab="time", ylab="Wt", ylim=c(0, max(W_field)) , cex=1.5)
points(W_DP, W_field, pch=16, col="dark green")
title("Comparison of model points with field points for Lamspar II")
legend(x=600, y=150, legend=c("model", "field"), col=c("black", "dark green"), pch=c(16,16))


