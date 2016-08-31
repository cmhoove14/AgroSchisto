#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Read in function
require(deSolve)
require(ggplot2)
st.er <- function(x) {
  sd(x)/sqrt(length(x))
} #Function to calculate standard error of the mean
cbPalette <- c("#999999", "yellow", "red", "green", "orange", "darkgreen", "pink", "blue")

schisto_master_mda=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    P=n[6]
    N=S+E+I
    
    W = cov*Wt + (1-cov)*Wu #weighting treated and untreated populations
    
    #Miracidia production; function of adult female worms alive in the system (W*H*0.5) assuming 1:1 sex ratio,
    #mating probability function (gamma) from Anderson and May (where it is phi),
    
    fx<-function(x, mean.worm = W, clump = k){
      (1-cos(x))/((1+((mean.worm/(clump+mean.worm))*cos(x)))^1+clump)
    }
    gamma = 1 - ((1-(W/(W+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value)
    
    #per-reproductive female egg production (m) per mL urine,
    #u_H daily mL urine production of an infected individual (u_H)
    #egg viability or the fraction of eggs that successfully hatch into viable miracidia (v) * potential ag effect(vq),
    
    M=((0.5*W*H)*gamma)#*(m*u_H)*(v*vq)
    
    pred= (alpha*P)/(1+(alpha*N*Th)) #death rate of snails due to predators (Prawns)
    
    dSdt= f_N*(1-(N/(phi_N*phi_Nq)))*(S+E) - 
      mu_N*S - pred*S - beta*M*S #Susceptible snails
    
    dEdt= beta*M*S - (mu_N+pred+sigma)*E #Exposed snails
    
    dIdt= sigma*E - (mu_N+pred+mu_I)*I #Infected snails
    
    #worm burden in human
    dWtdt= (lamda*I) - ((mu_W+mu_H)*Wt) - (eff*Wt*mda )
    dWudt= (lamda*I) - ((mu_W+mu_H)*Wu)
    
    dPdt= f_P*(1-(P/phi_P))*P-(mu_P+muPq)*P #prawn population
    
    
    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt, dPdt)))
  }) 
} 

#List parameters and values #####################
parameters=c( #Updated as of 4/19 to get values directly from cited sources instead of drawing directly from PNAS 
  ##standard snail parameters 
  f_N=0.1, # recruitment rate: snails/snail/day
  phi_N=10000, # carrying capacity: max snail population (corresponds to ~50/m^2)
  z=0.5, #Proportion of exposed snails that reproduce: density dependent, but assumed constant here
  mu_N=0.017, #Mortality rate from Anderson and May (from Chu 1966): deaths/snail/day
  sigma=0.023, #Transition rate from exposed to infected; ~latent period
  mu_I=0.083, #additional snail death due to infection
  ## snail parameters impacted by agrochemicals
  #f_Nq=1, #Not affected in mesocosm
  phi_Nq=1, #Scalar of snail carrying capacity by chemical concentration INFORMED BY BOTTOM UP EFFECTS IN MESOCOSM
  #mu_Nq=0, #Chem concentrations too low to affect snails in mesocosom experiments
  
  
  #prawn parameters
  alpha=0.003, #attack rate
  Th=0.067,#~Prawn predation limit for procambarius clarkii: ~15 snails/day
  f_P=0.117,#prawn birth rate: prawns/prawn/day
  phi_P=120,  #prawn carrying capacity ~0.6/m^2
  mu_P= 0.03883984, #observed daily 24 hr prawn mortality rate in chlorP-free tanks of mesocosm#0.00137, previous rate
  
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
  lamda=1.5e-5, #snail-to-man transmission: p(infected snail sheds cercariae that infects human and reaches adulthood)
  beta=2.5e-5, #man-to-snail transmission: p(mated female worm produces a miracidia that infects a snail)
  
  #Treatment parameters
  eff=0.99, #Efficacy of PZQ treatment (fraction of worms cleared)
  cov=0.5, #coverage of PZQ treatment in the human population
  mda=0 #binary to indicate if MDA is applied
  
)

#Set initial values ##########
nstart=c(S=8000,E=0,I=0, Wt=10, Wu=10, P=0)
yrs=30
time=seq(0,365*yrs,1)

#Run model and visualize for initial estimates #################
output=as.data.frame(ode(nstart,time,schisto_master_mda,parameters))
eqbm=output[dim(output)[1],]

snail.prev=(eqbm$I)/(eqbm$E+eqbm$I+eqbm$S)

plot(output$time, output$S, type='l', xlab="time",ylab="System Variables", 
     ylim=c(0,max( output$S+output$E+output$I )), col='blue', lwd=2)
  lines(output$time,output$E,col='orange', lwd=2)
  lines(output$time,output$I,col='red', lwd=2)
  lines(output$time,output$S+output$E+output$I,col='black', lwd=2)

# simulating mda in the Senegal site, with the simplified Halstead expt based model ####################
  Senegal_mda_halstead<-function(nstart, parameters, mda_days ){
    
    #mda function indicating when PZQ is applied
    mda<-function(t){
      #ifelse( t==(28*7) | t==(31*7), 1, 0)
      ifelse( t %in% mda_days, 1, 0)
    }
    
    # first run to equilibrium for 7 months
    output_all<-numeric()
    
    t_eqbmtoPZQ1<-28*7 # 28 weeks or 7 months
    time=seq(0,t_eqbmtoPZQ1,1)
    
    output1=as.data.frame(ode(nstart,time,schisto_master_mda,parameters)) 
    output_7months<-output1[dim(output1)[1],]
    time_last<-output1$time[dim(output1)[1]]
    output_all<-rbind(output_all, output1)
    
    ##### apply PZQ using mda flag and then run for 1 day
    parameters["mda"]<-mda(time_last+1)
    time<-seq(time_last, time_last+1,1)
    nstart<-c(S=output_7months$S, E=output_7months$E,I=output_7months$I, Wt=output_7months$Wt, Wu=output_7months$Wu, P=0)
    output2=as.data.frame(ode(nstart,time,schisto_master_mda,parameters)) 
    output_PZQ1<-output2[dim(output2)[1],]
    output_all<-rbind(output_all, output2[-1,])
    time_last<-output2[dim(output2)[1],1]
    
    #### apply no PZQ using mda flag and run for 3 weeks
    parameters["mda"]<-mda(time_last+1)
    time<-seq(time_last, time_last+(3*7),1)
    nstart<-c(S=output_PZQ1$S, E=output_PZQ1$E,I=output_PZQ1$I, Wt=output_PZQ1$Wt, Wu=output_PZQ1$Wu, P=0)
    output3=as.data.frame(ode(nstart,time,schisto_master_mda,parameters)) 
    output_3wks<-output3[dim(output3)[1],]
    output_all<-rbind(output_all, output3[-1,])
    time_last<-output3[dim(output3)[1],1]
    
    #apply PZQ using mda flag and run for 1 day
    parameters["mda"]<-mda(time_last)
    time<-seq(time_last, time_last+1,1)
    nstart<-c(S=output_3wks$S, E=output_3wks$E,I=output_3wks$I, Wt=output_3wks$Wt, Wu=output_3wks$Wu, P=0)
    output4=as.data.frame(ode(nstart,time,schisto_master_mda,parameters)) 
    output_PZQ2<-output4[dim(output4)[1],]
    output_all<-rbind(output_all, output4[-1,])
    time_last<-output4[dim(output4)[1],1]
    
    
    #apply no PZQ using mda flag and run for 5 months (20 weeks)
    parameters["mda"]<-mda(time_last)
    time<-seq(time_last, time_last+(20*7),1)
    nstart<-c(S=output_PZQ2$S, E=output_PZQ2$E,I=output_PZQ2$I, Wt=output_PZQ2$Wt, Wu=output_PZQ2$Wu, P=0)
    output5=as.data.frame(ode(nstart,time,schisto_master_mda,parameters)) 
    output_20wks<-output5[dim(output5)[1],]
    output_all<-rbind(output_all, output5[-1,])
    time_last<-output5[dim(output5)[1],1]
    
    output_all
    
  }
  
  mda_days<-c((28*7+1),(31*7+1)) #end of 28 weeks and 31 weeks
  nstart=c(S=8000,E=0,I=0, Wt=10, Wu=10, P=0) 
  output<-Senegal_mda_halstead(nstart, parameters, mda_days)
  #quartz()
  par(mfrow=c(2,1))
  plot(output$time, output$Wt, type='l', lwd=2, col=1, ylim=c(0, max(output$Wu)))
  lines(output$time, output$Wu, lwd=2, col=2)
  legend(x=300, y=20, legend=c("Wt", "Wu"), col=c(1,2), lwd=2)
  cov<-parameters["cov"]
  W<-(cov*output$Wt) + ((1-cov)*output$Wu)
  plot(output$time,W ,  type='l', lwd=2, col=3, ylim=c(0, max(W)))
  
  eqbm<-output[dim(output)[1],]
  params<-parameters

#Fit the beta of the villages to their worm burden W at baseline and after PZQ #############
# The baseline condition is no prawns, no agro; Village Epi data used is Lampsar II.
# The egg burdens -- measured as eggs/10mL urine -- are divided by 3.6 (page 3, French et al 2015) 
#  to estimate female worms/infected individual. Worm burden is then achieved by multiplying by 2 (assuming 1:1 sex ratio)
  
EggToWormConvert<-3.6
W_baseline<-(6.5/EggToWormConvert)*2
W_baseline_k<-0.08
W_baseline_sd<-sqrt( (W_baseline)+((W_baseline^2)/W_baseline_k ) )
N<-129
W_baseline_se<-W_baseline_sd/sqrt(N)
#W_limits_baseline<-c( W_baseline-(1.96*W_baseline_se), W_baseline+(1.96*W_baseline_se))

W_5month_field_mean<-(1.5/EggToWormConvert)*2
W_5month_field_k<-0.02
W_5month_field_sd<-sqrt( (W_5month_field_mean)+((W_5month_field_mean^2)/W_5month_field_k ) )
N<-129
W_5month_field_se<-W_5month_field_sd/sqrt(N)
timepoints<-c(0, 28*7, 51*7)

W_Feb13_field_mean<-(161/EggToWormConvert)*2
W_Feb13_field_k<-0.21
W_Feb13_field_sd<-sqrt( (W_Feb13_field_mean)+((W_Feb13_field_mean^2)/W_Feb13_field_k ) )
N<-129
W_Feb13_field_se<-W_Feb13_field_sd/sqrt(N)

W_Sep13_field_mean<-(17.6/EggToWormConvert)*2
W_Sep13_field_k<-0.29
W_Sep13_field_sd<-sqrt( (W_Sep13_field_mean)+((W_Sep13_field_mean^2)/W_Sep13_field_k ) )
N<-129
W_Sep13_field_se<-W_Sep13_field_sd/sqrt(N)
  
#W_limits_5months<-c(W_5month_field_mean-(1.96*W_5month_field_se), W_5month_field_mean+(1.96*W_5month_field_se))



# Obtain the negative binomial distribution for W ###########
#The ecological model for the neg prob for worm burden
# i number of parasites per person
# k clumping parameter
# m mean burden
# refer wikipedia page - https://en.wikipedia.org/wiki/Negative_binomial_distribution
# formula refer - http://influentialpoints.com/Training/negative_binomial_distribution-principles-properties-assumptions.htm
Prob_negbin<-function(i,k,m){
  p<-(gamma(k+i)/(gamma(i+1) * gamma(k))) * (1 + (m/k))^(-k-i) * ((m/k)^i)
  p
}


W_test<-0
Prob_negbin(i=W_test, k=W_baseline_k, m=W_baseline)
dnbinom(W_test, s=W_baseline_k, mu=W_baseline)

quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(1:100,Prob_negbin(i=1:100, k=W_baseline_k, m=W_baseline) , type ="l", col=1, lwd=2)
  plot(1:100,Prob_negbin(i=1:100, k=W_5month_field_k, m=W_5month_field_mean) , type ="l", col=1, lwd=2)
  plot(1:100,Prob_negbin(i=1:100, k=W_Feb13_field_k, m=W_Feb13_field_mean) , type ="l", col=1, lwd=2)
  plot(1:100,Prob_negbin(i=1:100, k=W_Sep13_field_k, m=W_Sep13_field_mean) , type ="l", col=1, lwd=2)
  title("Negative binomial disbn", outer=TRUE)

Prob_gaussian<-function(y, mu, sd){
  p<-(1/(sqrt(2*pi)*sd) ) * exp( -((y-mu)^2)/(2*sd^2) )
  p
}

range<-seq(from=0, to=100, by=100/1000)
quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(range,Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se) , type ="l", 
       col=1, lwd=2, xlab="W_range", ylab="PDF", main="Baseline")
  plot(range,Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se) , type ="l", 
       col=1, lwd=2, xlab="W_range", ylab="PDF", main="First Follow Up")
  plot(range,Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) , type ="l", 
       col=1, lwd=2, xlab="W_range", ylab="PDF", main="Second Follow Up")
  plot(range,Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se) , type ="l", 
       col=1, lwd=2, xlab="W_range", ylab="PDF", main="Third Follow Up")
  title("Distribution of the mean worm burden", outer=TRUE)

prod13<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * 
  as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) 

prod134<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * 
  as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) * 
  as.vector(Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se)) 

prod123<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * 
  as.vector(Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se)) * 
  as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) 

prod1234<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * 
  as.vector(Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se))  * 
  as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) * 
  as.vector(Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se)) 

quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(range, prod13, type='l', col=2, main="points 1+3")
  plot(range, prod134, type='l', col=2, main="points 1+3+4")
  plot(range, prod123, type='l', col=2, main="points 1+2+3")
  plot(range, prod1234, type='l', col=2, main="points 1+2+3+4")

quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(range, log(prod13), type='l', col=3, main="points 1+3")
  plot(range, log(prod134), type='l', col=3, main="points 1+3+4")
  plot(range, log(prod123), type='l', col=3, main="points 1+2+3")
  plot(range, log(prod1234), type='l', col=3, main="points 1+2+3+4")

range<-seq(from=0, to=5, by=1/1000)
  prod13<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * 
    as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) 
  
  prod134<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * 
    as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) * 
    as.vector(Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se)) 
  
  prod123<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * 
    as.vector(Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se)) * 
    as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) 
  
  prod1234<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * 
    as.vector(Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se))  * 
    as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) * 
    as.vector(Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se)) 

quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(range, prod13, type='l', col=2, main="points 1+3")
  plot(range, prod134, type='l', col=2, main="points 1+3+4")
  plot(range, prod123, type='l', col=2, main="points 1+2+3")
  plot(range, prod1234, type='l', col=2, main="points 1+2+3+4")


quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(range, log(prod13), type='l', col=3, main="points 1+3")
  plot(range, log(prod134), type='l', col=3, main="points 1+3+4")
  plot(range, log(prod123), type='l', col=3, main="points 1+2+3")
  plot(range, log(prod1234), type='l', col=3, main="points 1+2+3+4")

# Generate range of Beta values to use ##########
  min_beta<-2e-5
  max_beta<-4e-5
  beta_range<-seq(from=min_beta, to=max_beta, by=(max_beta-min_beta)/100 ) 

# generate a range of lamda ###############

  min_lamda<-1e-5
  max_lamda<-3e-5
  lamda_range<-seq(from=min_lamda, to=max_lamda, by=(max_lamda-min_lamda)/100) 
  nstart=c(S=8000,E=0,I=0, Wt=10, Wu=10, P=0) 

#Generate matrices to fill ###############    
  W_baseline_range<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  W_postmda_July12_range<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  W_postmda_Feb13_range<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  W_postmda_Sep13_range<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  
  loglikelihood_range<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range1<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range2<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range3<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range4<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range13<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range14<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range134<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range1234<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  
#Fitting function using max likelihood #####################
  output_all<-numeric()
  timepoints<-c(198, 362, 576, 760) #days when MDA was applied
  W_timepoints<-c(W_baseline, W_5month_field_mean, W_Feb13_field_mean, W_Sep13_field_mean)
  for(i in 1:length(beta_range)){
    for(j in 1:length(lamda_range)){
      
      output_all<-numeric()
      params<-parameters
      params["beta"]<-beta_range[i]
      params["lamda"]<-lamda_range[j]
      params["mda"]<-0 #no mda
      
      #### run to equilibrium
      time<-seq(from=0, to=50*365, by=1)
      nstart=c(S=8000,E=300,I=500, Wt=10, Wu=10, P=0) 
      output<-as.data.frame(ode(nstart,time,schisto_master_mda,params))
      output_eqbm<-output[dim(output)[1],]
      cov<-params["cov"]
      W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
      W_baseline_range[i,j]<-W_eqbm
      InfectedSnails<-((output_eqbm$I + output_eqbm$E)/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100
      #apply mda
      nstart1<-c(S=output_eqbm$S,E=output_eqbm$E,I=output_eqbm$I, Wt=output_eqbm$Wt, Wu=output_eqbm$Wu, P=0)                   
      output_July12<-Senegal_mda_halstead(nstart1, params, mda_days) 
      output_eqbm<-output_July12[dim(output_July12)[1],]
      cov<-params["cov"]
      W_postmda<-as.numeric((cov*output_eqbm$Wt ) + ((1-cov)*output_eqbm$Wu ) )
      #W_postmda_range[i]<-W_postmda
      W_postmda_July12_range[i,j]<-output_eqbm$Wt
      
      output_all<-rbind(output_all, output_July12 ) ##
      
      #### apply mda at first survey in July 2012.
      
      ##### apply PZQ using mda flag and then run for 1 day
      time_last<-output_July12[dim(output_July12)[1],1]
      
      params["mda"]<-1
      time<-seq(time_last, time_last+1,1)
      nstart<-c(S=output_eqbm$S, E=output_eqbm$E,I=output_eqbm$I, Wt=output_eqbm$Wt, Wu=output_eqbm$Wu, P=0)
      output_mda_July12=as.data.frame(ode(nstart,time,schisto_master_mda,params)) 
      output_PZQ2<-output_mda_July12[dim(output_mda_July12)[1],]
      
      output_all<-rbind(output_all, output_mda_July12[-1,]) ##
      time_last<-output_mda_July12[dim(output_mda_July12)[1],1]
      
      #### apply no PZQ using mda flag and run for 7 months, July to Feb = (4*31)+(3*30) = 214 days.
      params["mda"]<-0
      time<-seq(time_last, time_last+214,1)
      nstart<-c(S=output_PZQ2$S, E=output_PZQ2$E,I=output_PZQ2$I, Wt=output_PZQ2$Wt, Wu=output_PZQ2$Wu, P=0)
      output_Feb13=as.data.frame(ode(nstart,time,schisto_master_mda,params)) 
      output_Feb13_eqbm<-output_Feb13[dim(output_Feb13)[1],]
      
      W_postmda_Feb13_range[i,j]<-output_Feb13_eqbm$Wt
      output_all<-rbind(output_all, output_Feb13[-1,]) ##
      
      time_last<-output_Feb13[dim(output_Feb13)[1],1]
      
      
      
      
      #### apply mda at next survey in Feb 2013.
      ##### apply PZQ using mda flag and then run for 1 day
      
      params["mda"]<-1
      time<-seq(time_last, time_last+1,1)
      nstart<-c(S=output_Feb13_eqbm$S, E=output_Feb13_eqbm$E,I=output_Feb13_eqbm$I, Wt=output_Feb13_eqbm$Wt, Wu=output_Feb13_eqbm$Wu, P=0)
      output_mda_Feb13=as.data.frame(ode(nstart,time,schisto_master_mda,params)) 
      output_PZQ3<-output_mda_Feb13[dim(output_mda_Feb13)[1],]
      
      output_all<-rbind(output_all, output_mda_Feb13[-1,]) ##
      
      time_last<-output_mda_Feb13[dim(output_mda_Feb13)[1],1]
      
      #### apply no PZQ using mda flag and run for 6 months, Feb to Aug 2013 = (4*31)+(2*30) = 184 days.
      params["mda"]<-0
      time<-seq(time_last, time_last+184,1)
      nstart<-c(S=output_PZQ3$S, E=output_PZQ3$E,I=output_PZQ3$I, Wt=output_PZQ3$Wt, Wu=output_PZQ3$Wu, P=0)
      output_Sep13=as.data.frame(ode(nstart,time,schisto_master_mda,params)) 
      output_Sep13_eqbm<-output_Sep13[dim(output_Sep13)[1],]
      
      W_postmda_Sep13_range[i,j]<-output_Sep13_eqbm$Wt
      output_all<-rbind(output_all, output_Sep13[-1,]) ##
      
      time_last<-output_Sep13[dim(output_Sep13)[1],1]
      
      #     quartz()
      #     plot(output_all$time, output_all$Wt, col=1, lwd=2, type="l", ylim=c(0, max(output_all$Wt)), main="ind4")
      #     points(timepoints, W_timepoints, pch=13, col=2)
      #     quartz()
      #     par(mfrow=c(2,2), oma=c(0,0,2,0))
      #     plot(output_all$time, output_all$S, type="l", col=1, main="S")
      #     plot(output_all$time, output_all$E, type="l", col=2, main="E")
      #     plot(output_all$time, output_all$I, type="l", col=3, main="I")
      #     plot(output_all$time, output_all$Wt, type="l", col=4, main="Wt")
      #     title(get_Ro_mesocosm_withPrawns(params, output_eqbm), outer=TRUE)
      #     output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I)*100
      #     
      LL1<-Prob_gaussian(y=W_baseline_range[i,j], mu=W_baseline, sd=W_baseline_se)  # only baseline
      LL2<-Prob_gaussian(y=W_postmda_July12_range[i,j], mu=W_5month_field_mean, sd=W_5month_field_se)
      LL3<-Prob_gaussian(y=W_postmda_Feb13_range[i,j], mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3
      LL4<-Prob_gaussian(y=W_postmda_Sep13_range[i,j], mu=W_Sep13_field_mean, sd=W_Sep13_field_se) #only data point 4
      LL13<-LL1 * LL3
      LL14<-LL1 * LL4
      LL134<-LL1 * LL3 * LL4
      LL1234<- LL1 * LL2* LL3 * LL4
      
      likelihood_range1[i,j]<-LL1
      likelihood_range2[i,j]<-LL2
      likelihood_range3[i,j]<-LL3
      likelihood_range4[i,j]<-LL4
      likelihood_range13[i,j]<-LL13
      likelihood_range14[i,j]<-LL14
      likelihood_range134[i,j]<-LL134
      likelihood_range1234[i,j]<-LL1234
      
      print( c(i,j, W_baseline_range[i,j], InfectedSnails )) 

      
    }# end of j loop
  }# end of i loop
   
  quartz()
  par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(lamda_range, W_baseline_range, type="l", lwd=2,  xlab="lamda", ylab="W", main="Baseline")
  plot(lamda_range, W_postmda_July12_range, type="l", lwd=2, main="First Follow Up",  xlab="lamda", ylab="W")
  plot(lamda_range, W_postmda_Feb13_range, type="l", lwd=2, main="Second Follow Up",  xlab="lamda", ylab="W")
  plot(lamda_range, W_postmda_Sep13_range, type="l", lwd=2, main="Third Follow Up",  xlab="lamda", ylab="W")
  title("Variation of Worm Burden with Lamda", outer="TRUE")
#   
  quartz()
  par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(lamda_range, likelihood_range1, type="l", lwd=2, main="Baseline")
  plot(lamda_range, likelihood_range2, type="l", lwd=2, main="First follow up")
  plot(lamda_range, likelihood_range3, type="l", lwd=2, main="Second follow up")
  plot(lamda_range, likelihood_range4, type="l", lwd=2, main="Third follow up")
  title("Variation of Likelihood with Lamda", outer="TRUE")
  
  likelihood_range1[which(likelihood_range1==0)]<-1e-320
  likelihood_range2[which(likelihood_range2==0)]<-1e-320
  likelihood_range3[which(likelihood_range3==0)]<-1e-320
  likelihood_range4[which(likelihood_range4==0)]<-1e-320
  
  quartz()
  par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(lamda_range, -log(likelihood_range1), type="l", lwd=2, main="Baseline")
  plot(lamda_range, -log(likelihood_range2), type="l", lwd=2, main="First follow up")
  plot(lamda_range, -log(likelihood_range3), type="l", lwd=2, main="Second follow up")
  plot(lamda_range, -log(likelihood_range4), type="l", lwd=2, main="Third follow up")
  title("Variation of Neg Log Likelihood with Lamda", outer="TRUE")
  
  negLogLL1<- -log(likelihood_range1)
  negLogLL2<- -log(likelihood_range2)
  negLogLL3<- -log(likelihood_range3)
  negLogLL4<- -log(likelihood_range4)
  
  ind1<-which( negLogLL1==min(negLogLL1))
  ind2<-which( negLogLL2==min(negLogLL2))
  ind3<-which( negLogLL3==min(negLogLL3))
  ind4<-which( negLogLL4==min(negLogLL4))
  
  params<-parameters
  params["lamda"]<-lamda_range[ind1]
  Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
  Ro
  
  params<-parameters
  params["lamda"]<-lamda_range[ind2]
  Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
  Ro
  
  params<-parameters
  params["lamda"]<-lamda_range[ind3]
  Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
  Ro
  
  params<-parameters
  params["lamda"]<-lamda_range[ind4]
  Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
  Ro
  
  
  
  quartz()
  #par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(lamda_range[(ind1-20):(ind1+20)], -log(likelihood_range1)[(ind1-20):(ind1+20)], type="l", lwd=2, main="Baseline")
  abline(h=negLogLL1[ind1], col=3, lty=2)
  abline(h=negLogLL1[ind1]+1.92, col=3, lty=1, lwd=2)
  locator(n=2, type="p")

  quartz()
  plot(lamda_range[(ind2-20):(ind2+20)], -log(likelihood_range2)[(ind2-20):(ind2+20)], type="l", lwd=2, main="First follow up")
  abline(h=negLogLL2[ind2], col=3, lty=2)
  abline(h=negLogLL2[ind2]+1.92, col=3, lty=1, lwd=2)
  locator(n=2, type="p")
  
  
  quartz()
  plot(lamda_range[(ind3-60):(ind3+60)], -log(likelihood_range3)[(ind3-60):(ind3+60)], type="l", lwd=2, main="Second follow up")
  abline(h=negLogLL3[ind3], col=3, lty=2)
  abline(h=negLogLL3[ind3]+1.92, col=3, lty=1, lwd=2)
  locator(n=2, type="p")
  
  quartz()
  plot(lamda_range[(ind4-60):(ind4+60)], -log(likelihood_range4)[(ind4-60):(ind4+60)], type="l", lwd=2, main="Third follow up")
  abline(h=negLogLL4[ind4], col=3, lty=2)
  abline(h=negLogLL4[ind4]+1.92, col=3, lty=1, lwd=2)
  locator(n=2, type="p")
  
### combinations  
  quartz()
  par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(lamda_range, likelihood_range13, type="l", lwd=2, main="Baseline+Second follow up")
  plot(lamda_range, likelihood_range14, type="l", lwd=2, main="Baseline+Third follow up")
  plot(lamda_range, likelihood_range134, type="l", lwd=2, main="Baseline+Second+Third follow up")
  plot(lamda_range, likelihood_range1234, type="l", lwd=2, main="All 4 points")
  title("Variation of Likelihood with Lamda", outer="TRUE")
  
  likelihood_range13[which(likelihood_range13==0)]<-1e-320
  likelihood_range14[which(likelihood_range14==0)]<-1e-320
  likelihood_range134[which(likelihood_range134==0)]<-1e-320
  likelihood_range1234[which(likelihood_range1234==0)]<-1e-320
  
  negLogLL13<- -log(likelihood_range134)
  negLogLL14<- -log(likelihood_range14)
  negLogLL134<- -log(likelihood_range134)
  negLogLL1234<- -log(likelihood_range1234)
  
  ind13<-which( negLogLL13==min(negLogLL13))
  ind14<-which( negLogLL14==min(negLogLL14))
  ind134<-which( negLogLL134==min(negLogLL134))
  ind1234<-which( negLogLL1234==min(negLogLL1234))
  
  params<-parameters
  params["lamda"]<-lamda_range[ind13]
  Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
  Ro
  params<-parameters
  params["lamda"]<-lamda_range[ind14]
  Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
  Ro
  params<-parameters
  params["lamda"]<-lamda_range[ind134]
  Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
  Ro
  params<-parameters
  params["lamda"]<-lamda_range[ind1234]
  Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
  Ro
  
  
  
  
  quartz()
  plot(lamda_range[(ind13-10):(ind13+10)], negLogLL13[(ind13-10):(ind13+10)], type="l", lwd=2, main="Baseline+Second follow up")
  abline(h=min(negLogLL13), col=3, lty=2)
  abline(h=min(negLogLL13)+1.96, col=3, lty=3)
  locator(n=2, type="p")
  
  quartz()
  plot(lamda_range[(ind14-20):(ind14+20)],negLogLL14[(ind14-20):(ind14+20)], type="l", lwd=2, main="Baseline+Third follow up")
  abline(h=min(negLogLL14), col=3, lty=2)
  abline(h=min(negLogLL14)+1.96, col=3, lty=3)
  locator(n=2, type="p")
  
  quartz()
  plot(lamda_range[(ind134-20):(ind134+20)], negLogLL134[(ind134-20):(ind134+20)], type="l", lwd=2, main="Baseline+Second+Third follow up")
  abline(h=min(negLogLL134), col=3, lty=2)
  abline(h=min(negLogLL134)+1.96, col=3, lty=3)
  locator(n=2, type="p")
  
  quartz()
  plot(lamda_range[(ind1234-20):(ind1234+20)], negLogLL1234[(ind1234-20):(ind1234+20)], type="l", lwd=2, main="All 4 points")
  abline(h=min(negLogLL1234), col=3, lty=2)
  abline(h=min(negLogLL1234)+1.96, col=3, lty=3)
  locator(n=2, type="p")
  
  params<-parameters
  params["lamda"]<-lamda_range[ind3]
  Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
  Ro
 
 
#   
#   
# for(b in 1:3){
#   ### lets get the ouput curves for the 3 beta values
#   params<-parameters
#   params["phi_P"]<-1/CC_P
#   params["beta"]<-W_bestBeta[b]
#   params["phi_N"]<-phi_N_new
#   params["mda"]<-0 #no mda
#   #### run to equilibrium
#   nstart=c(S=3892,E=3750,I=1200, Wt=20, Wu=20, P=0) 
#   time<-seq(from=0, to=30*365, by=1)
#   output<-as.data.frame(ode(nstart,time,schisto_master_mda,params))
#   output_eqbm<-output[dim(output)[1],]
#   cov<-params["cov"]
#   W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
#   
#   #apply mda
#   nstart1<-c(S=output_eqbm$S,E=output_eqbm$E,I=output_eqbm$I, Wt=output_eqbm$Wt, Wu=output_eqbm$Wu, P=0)                   
#   output<-Senegal_mda_halstead(nstart1, params, mda_days) 
#   output_eqbm<-output[dim(output)[1],]
#   cov<-params["cov"]
#   W_postmda<-as.numeric((cov*output_eqbm$Wt ) + ((1-cov)*output_eqbm$Wu ) )
#   W_postmda<-output_eqbm$Wt
#   #lets plot and see
#   datapoints<-c(W_baseline_field[b], W_baseline_field[b], W_5month_field[b])
#   W<-as.numeric((cov*output$Wt ) + ((1-cov)*output$Wu ) )
#   W<-output$Wt
#   quartz()
#   plot(output$time, W , type="l", col=b, lwd=2, ylim=c(0,max(c( W_baseline_field,W) )), ylab="W", xlab="time" )
#   points(timepoints, c(W_baseline_field[1], W_baseline_field[1], W_5month_field[1]), pch=16, col=1 )
#   points(timepoints, c(W_baseline_field[2], W_baseline_field[2], W_5month_field[2]), pch=16, col=2 )
#   points(timepoints, c(W_baseline_field[3], W_baseline_field[3], W_5month_field[3]), pch=16, col=3 )
#   
#   #segments(timepoints[3], datapoints[3]-(1.96*W_5month_field_se[b]), timepoints[3], datapoints[3]+(1.96*W_5month_field_se[b]), lwd=2)
#   title(paste("Village Lampsar II ", "with best fit beta =",W_bestBeta[b] , sep=" "))
#   
#   
# }
# 
# 
# 
# 
# beta<-W_bestBeta[2]
# 
# params<-parameters
# 
# # Run the baseline model to equilibrium, with prwans and no agro
# 
# # get Ro and r(must be 1 at eqbm)
# 
# # Prawn free R0
# 
# 
# #Baseline R0 run ##############
# get_Ro(muPq=0, phi_Nq=1)  #R0 = 1.09706
# 
# #Prawn free R0 run ##########
# parameters["phi_P"]<-1e20
# get_Ro(muPq = 0, phi_Nq = 1) #R0 = 2.026636
# 
# #R0 run with prawn carrying capacity equivalent to that in mesocosm assuming scale up factor of 40 (5 sq meters to 200 sq meters) #########
# parameters["phi_P"]<-1/120
# get_Ro(muPq = 0, phi_Nq = 1) #R0 = 0.6730459
# #Mesocosm informed R0 runs#############   
# #Mesocosm baseline R0 run
# R0_base<-rep(0,10000)
# for(i in 1:length(R0_base)){
#   R0_base[i]<-get_Ro(muPq = 0, phi_Nq = phi_base[i])
# }
# #Mesocosm +Fertilizer R0 run
# R0_fert<-rep(0,10000)
# for(i in 1:length(R0_fert)){
#   R0_fert[i]<-get_Ro(muPq = 0, phi_Nq = phi_fert[i])
# }
# 
# #Mesocosm +Atrazine R0 run
# R0_atra<-rep(0,10000)
# for(i in 1:length(R0_atra)){
#   R0_atra[i]<-get_Ro(muPq = 0, phi_Nq = phi_atra[i])
# }
# 
# #Mesocosm +Chlorpyrifos R0 run
# R0_chlor<-rep(0,10000)
# for(i in 1:length(R0_chlor)){
#   R0_chlor[i]<-get_Ro(muPq = meso.rate.mod1[i], phi_Nq = phi_base[i])
# }
# 
# #Mesocosm +Fertilizer+Atrazine R0 run
# R0_atfe<-rep(0,10000)
# for(i in 1:length(R0_atfe)){
#   R0_atfe[i]<-get_Ro(muPq = 0, phi_Nq = phi_atfe[i])
# }
# 
# #Mesocosm +Fertilizer+Chlorpyrifos R0 run
# R0_chfe<-rep(0,10000)
# for(i in 1:length(R0_chfe)){
#   R0_chfe[i]<-get_Ro(muPq = meso.rate.mod1[i], phi_Nq = phi_fert[i])
# }
# 
# #Mesocosm +Atrazine+Chlorpyrifos R0 run
# R0_atch<-rep(0,10000)
# for(i in 1:length(R0_atch)){
#   R0_atch[i]<-get_Ro(muPq = meso.rate.mod1[i], phi_Nq = phi_atra[i])
# }
# 
# #Mesocosm +Fertilizer+Atrazine+Chlorpyrifos R0 run
# R0_tre<-rep(0,10000)
# for(i in 1:length(R0_tre)){
#   R0_tre[i]<-get_Ro(muPq = meso.rate.mod1[i], phi_Nq = phi_atfe[i])
# }
# 
# 
# #Plot results ###############
# 
# r0s<-data.frame("Treatment"=c('baseline', 'Fert only', 'Atra only', 'ChlorP Only', 
#                               'Fert+Atra', 'Fert+ChlorP', 'Atra+ChlorP', 'All three'),
#                 "meanR0"=c(mean(R0_base), mean(R0_fert), mean(R0_atra), mean(R0_chlor),
#                            mean(R0_atfe), mean(R0_chfe), mean(R0_atch), mean(R0_tre)),
#                 "st.erR0"=c(st.er(R0_base), st.er(R0_fert), st.er(R0_atra), st.er(R0_chlor),
#                             st.er(R0_atfe), st.er(R0_chfe), st.er(R0_atch), st.er(R0_tre)))
# 
# r0s$Treatment<-factor(r0s$Treatment, levels = c('baseline', 'Fert only', 'Atra only', 'ChlorP Only', 
#                                                 'Fert+Atra', 'Fert+ChlorP', 'Atra+ChlorP', 'All three'))
# 
# ggplot(r0s, aes(x=Treatment, y=meanR0, fill=Treatment))+
#   theme_bw()+
#   theme(axis.title=element_text(size=20),
#         axis.text=element_text(size=15))+
#   scale_fill_manual(values=cbPalette) +
#   ylab("Mean +/- SEM predicted R-0")+
#   geom_bar(position=position_dodge(), stat="identity", width = .7) +
#   geom_errorbar(aes(ymin=meanR0-st.erR0,
#                     ymax=meanR0+st.erR0),
#                 width=.2, position=position_dodge(.7))
