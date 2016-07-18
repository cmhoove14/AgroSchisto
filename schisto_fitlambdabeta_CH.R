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

#Model structure and equations ####################
schisto_halstead_2pops_mda=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    P=n[6]
    N=S+E+I
    
    W=(cov*Wt) + ((1-cov)*Wu) #weighting treated and untreated populations
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
    
    pred= (alpha*P)/(1+(alpha*N*Th)) #death rate of snails due to predators (Prawns)
    
    dSdt= f_N*(1-(N/(phi_N*phi_Nq)))*(S+E) - 
      mu_N*S - pred*S - beta*M*S #Susceptible snails
    
    dEdt= beta*M*S - (mu_N+pred+sigma)*E #Exposed snails
    
    dIdt= sigma*E - (mu_N+pred+mu_I)*I #Infected snails
    
    #worm burden in human
    dWtdt= (lamda*I) - ((mu_W+mu_H)*Wt)
    dWudt= (lamda*I) - ((mu_W+mu_H)*Wu)
    
    
    dPdt= f_P*(1-(P/phi_P))*P-(mu_P+muPq)*P #prawn population
    
    
    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt, dPdt)))
  }) 
} 

fx<-function(x, mean.worm = W, clump = k){
  (1-cos(x))/((1+((mean.worm/(clump+mean.worm))*cos(x)))^1+clump)
}

#List parameters and values #####################
parameters_2pops_mda=c(
  ##standard snail parameters 
  f_N=0.1, # recruitment rate: snails/snail/day
  phi_N=10000, # carrying capacity: max snail population (corresponds to ~50/m^2)
  z=0.5, #Proportion of exposed snails that reproduce: density dependent, but assumed constant here
  mu_N=1/60, #Mortality rate from Anderson and May (from Chu 1966): deaths/snail/day
  sigma=1/40, #Transition rate from exposed to infected; ~latent period
  mu_I=1/10 - 1/60, #additional snail death due to infection
  ## snail parameters impacted by agrochemicals
  f_Nq=1, #Not affected in mesocosm
  phi_Nq=1, #Scalar of snail carrying capacity by chemical concentration INFORMED BY BOTTOM UP EFFECTS IN MESOCOSM
  #mu_Nq=0, #Chem concentrations too low to affect snails in mesocosom experiments
  
  
  #prawn parameters
  alpha=0.003, #attack rate
  Th=0.067,#~Prawn predation limit for procambarius clarkii: ~15 snails/day
  f_P=0.117,#prawn birth rate: prawns/prawn/day
  phi_P=120,  #prawn carrying capacity ~0.6/m^2
  mu_P= 0.038095238, #observed daily 24 hr prawn mortality rate in chlorP-free tanks of mesocosm#0.00137, previous rate
  
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
  k=0.08, #clumping parameter of the negative binomial distribution
  u_H=1200, #mL urine per human/day (approximate, ranges from 800 - 2000)
  
  #Transmission parameters
  lamda=2.5e-5, #snail-to-man transmission: p(infected snail sheds cercariae that infects human and reaches adulthood)
  beta=1.75e-5, #man-to-snail transmission: p(mated female worm produces a miracidia that infects a snail)
  
  
  #treatment parameters
  cov=0.43, #coverage of treatment across the population, Lampsar I = 100/1000 =0.1 %, Lampsar II = 129/300 = 43%
  eff=0.95, # efficiency of the drug
  mda=0 # flag to indicate if mda is applied or not
)

params<-parameters_2pops_mda
parameters<-params

#R0 function ###########
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
  
  Ro_est <- sqrt((0.5*beta*H*N_eq*lamda*sigma)/((mu_W+mu_H)*(mu_N+pred+sigma)*(mu_N+pred+mu_I)))
  
  return(c(N_eq,P_eq,Ro_est ))
  
}  

p.dead = parameters["f_P"] - parameters["mu_P"] #muPq value at which predator population is eliminated 

#Run model with base parameters and no agrochemical effects and save equilibrium values for 5 state variables #############
time=seq(0,365*50,1) #50 years to ensure equilibrium is reached
nstart=c(S=4000,E=2000,I=500, Wt=10, Wu=10, P=0)

output.b=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,parameters)) 

eqbm.b<-output.b[dim(output.b)[1],]
eqbm.b
snail.prev<-eqbm.b$I/(eqbm.b$S+eqbm.b$E+eqbm.b$I)
snail.prev

plot(output.b$time, output.b$S, type='l', xlab="time",ylab="System Variables", 
     ylim=c(0,max( output.b$S+output.b$E+output.b$I )), col='blue', lwd=2)
  lines(output.b$time,output.b$E,col='orange', lwd=2)
  lines(output.b$time,output.b$I,col='red', lwd=2)
  lines(output.b$time,output.b$S+output.b$E+output.b$I,col='black', lwd=2)

# simulating mda in the Senegal site, with the simplified Halstead expt based model ###########

Senegal_mda_halstead<-function(nstart, parameters_2pops_mda, mda_days ){
  
  #mda function indicating when PZQ is applied
  mda<-function(t){
    #ifelse( t==(28*7) | t==(31*7), 1, 0)
    ifelse( t %in% mda_days, 1, 0)
  }
  
  #### first run to equilibrium for 7 months
  output_all<-numeric()
  
  t_eqbmtoPZQ1<-28*7 # 28 weeks or 7 months
  time=seq(0,t_eqbmtoPZQ1,1)
  
  output1=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,parameters_2pops_mda)) 
  
  output_7months<-output1[dim(output1)[1],]
  time_last<-output1$time[dim(output1)[1]]
  output_all<-rbind(output_all, output1)
  
  ##### apply PZQ using mda flag and then run for 1 day
  parameters_2pops_mda["mda"]<-mda(time_last+1)
  time<-seq(time_last, time_last+1,1)
  
  nstart<-c(S=output_7months$S, 
            E=output_7months$E,
            I=output_7months$I, 
            Wt=output_7months$Wt, 
            Wu=output_7months$Wu, 
            P=0)
  
  output2=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,parameters_2pops_mda)) 
  
  #Apply MDA to treated population (Wt) assuming 95% clearance of adult worms
  output2[2,5]<-output2[2,5] - output2[2,5]*.95
  
  output_PZQ1<-output2[dim(output2)[1],]
  output_all<-rbind(output_all, output2[-1,])
  time_last<-output2[dim(output2)[1],1]
  
  # apply no PZQ using mda flag and run for 3 weeks
  parameters_2pops_mda["mda"]<-mda(time_last+1)
  time<-seq(time_last, time_last+(3*7),1)
  
  nstart<-c(S=output_PZQ1$S, 
            E=output_PZQ1$E,
            I=output_PZQ1$I, 
            Wt=output_PZQ1$Wt, 
            Wu=output_PZQ1$Wu, 
            P=0)
  
  output3=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,parameters_2pops_mda)) 
  
  output_3wks<-output3[dim(output3)[1],]
  output_all<-rbind(output_all, output3[-1,])
  time_last<-output3[dim(output3)[1],1]
  
  #apply PZQ using mda flag and run for 1 day
  parameters_2pops_mda["mda"]<-mda(time_last)
  time<-seq(time_last, time_last+1,1)
  
  nstart<-c(S=output_3wks$S, 
            E=output_3wks$E,
            I=output_3wks$I, 
            Wt=output_3wks$Wt, 
            Wu=output_3wks$Wu, 
            P=0)
  
  output4=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,parameters_2pops_mda)) 
  
  #Apply MDA to treated population (Wt) assuming 95% clearance of adult worms
  output4[2,5]<-output4[2,5] - output4[2,5]*.95
  
  output_PZQ2<-output4[dim(output4)[1],]
  output_all<-rbind(output_all, output4[-1,])
  time_last<-output4[dim(output4)[1],1]
  
  
  #apply no PZQ using mda flag and run for 5 months (20 weeks)
  parameters_2pops_mda["mda"]<-mda(time_last)
  time<-seq(time_last, time_last+(20*7),1)
  
  nstart<-c(S=output_PZQ2$S, 
            E=output_PZQ2$E,
            I=output_PZQ2$I, 
            Wt=output_PZQ2$Wt, 
            Wu=output_PZQ2$Wu, 
            P=0)
  
  output5=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,parameters_2pops_mda)) 
  output_20wks<-output5[dim(output5)[1],]
  output_all<-rbind(output_all, output5[-1,])
  time_last<-output5[dim(output5)[1],1]
  
  output_all
  
}

mda_days<-c((28*7+1),(31*7+1)) #end of 28 weeks and 31 weeks
output2<-Senegal_mda_halstead(nstart, parameters_2pops_mda, mda_days)

plot(output2$time, output2$Wt, type='l', lwd=2, col=1, ylim=c(0, max(output2$Wu)))
  lines(output2$time, output2$Wu, lwd=2, col=2)
  cov<-parameters_2pops_mda["cov"]
  W<-(cov*output2$Wt) + ((1-cov)*output2$Wu)
  lines(output2$time,W , lwd=2, col=3)
  legend('bottomleft', legend=c("Wt", "Wu", "W_mean"), col=c(1,2,3), lwd=2)


#Prevalence function given W and k ##################################

k<-0.25 # clumping parameter of the negative binomial distribution
Prevalence <- function(W, k) {
  p=1 - (1/(1+W/k)^(k))*(1+W/(1+W/k)) # fraction of humans with at least 2 parasites
  return(p)
}

#Gather up epi data to fit to; Village used is Lampsar II. ##############
#egg burdens are divided by 3.6 which gives mean number of mated females per person
  #based on 3.6 eggs/10mL/mated female from cheever study
#Baseline epi data ##############
  EggToWormConvert<-3.6
  W_baseline<-6.5/EggToWormConvert
  W_baseline_k<-0.08
  W_baseline_sd<-sqrt( (W_baseline)+((W_baseline^2)/W_baseline_k ) )
  N<-129
  W_baseline_se<-W_baseline_sd/sqrt(N)

#5 month epi data point estimates ##################
  W_5month_field_mean<-1.5/EggToWormConvert
  W_5month_field_k<-0.02
  W_5month_field_sd<-sqrt( (W_5month_field_mean)+((W_5month_field_mean^2)/W_5month_field_k ) )
  N<-129
  W_5month_field_se<-W_5month_field_sd/sqrt(N)
  timepoints<-c(0, 28*7, 51*7)

#Feb 13 epi data estimates (in the middle of high transmission season) ################
  W_Feb13_field_mean<-161/EggToWormConvert
  W_Feb13_field_k<-0.21
  W_Feb13_field_sd<-sqrt( (W_Feb13_field_mean)+((W_Feb13_field_mean^2)/W_Feb13_field_k ) )
  N<-129
  W_Feb13_field_se<-W_Feb13_field_sd/sqrt(N)

#Sept 13 epi data estimates (in the middle of high transmission season) ################
  W_Sep13_field_mean<-17.6/EggToWormConvert
  W_Sep13_field_k<-0.29
  W_Sep13_field_sd<-sqrt( (W_Sep13_field_mean)+((W_Sep13_field_mean^2)/W_Sep13_field_k ) )
  N<-129
  W_Sep13_field_se<-W_Sep13_field_sd/sqrt(N)

#Get neg binomial distribution for Epi datapoints given estimated W and k #################
  Prob_negbin<-function(i,k,m){
    p<-(gamma(k+i)/(gamma(i+1) * gamma(k))) * (1 + (m/k))^(-k-i) * ((m/k)^i)
    p
  }
# i number of parasites per person
# k clumping parameter
# m mean burden
# refer wikipedia page - https://en.wikipedia.org/wiki/Negative_binomial_distribution
# formula refer - http://influentialpoints.com/Training/negative_binomial_distribution-principles-properties-assumptions.htm

  W_test<-0
  Prob_negbin(i=W_test, k=W_baseline_k, m=W_baseline)
  dnbinom(W_test, s=W_baseline_k, mu=W_baseline)

#Plot probability distributions of parasite burdens ###################
par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(1:250,Prob_negbin(i=1:250, k=W_baseline_k, m=W_baseline) , type ="l",lwd=2,
       ylab = 'prob individual has x worms', xlab = '#adult worms')
  plot(1:250,Prob_negbin(i=1:250, k=W_5month_field_k, m=W_5month_field_mean) , type ="l", lwd=2,
       ylab = 'prob individual has x worms', xlab = '#adult worms')
  plot(1:250,Prob_negbin(i=1:250, k=W_Feb13_field_k, m=W_Feb13_field_mean) , type ="l", lwd=2,
       ylab = 'prob individual has x worms', xlab = '#adult worms')
  plot(1:250,Prob_negbin(i=1:250, k=W_Sep13_field_k, m=W_Sep13_field_mean) , type ="l", lwd=2,
       ylab = 'prob individual has x worms', xlab = '#adult worms')
  title("Negative binomial disbn", outer=TRUE)
  
  
Prob_gaussian<-function(y, mu, sd){
  p<-(1/(sqrt(2*pi)*sd) ) * exp( -((y-mu)^2)/(2*sd^2) )
  p
}

#plot probability distribution functions of gaussian worm distribution? #############
#Range of 0-250 ##############
range<-seq(from=0, to=250, by=250/1000)
par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(range,Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se) , type ="l", lwd=2, 
       xlab="W_range", ylab="PDF", main="Baseline")
  plot(range,Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se) , type ="l", 
       lwd=2, xlab="W_range", ylab="PDF", main="First Follow Up")
  plot(range,Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) , type ="l", 
       lwd=2, xlab="W_range", ylab="PDF", main="Second Follow Up")
  plot(range,Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se) , type ="l", 
       lwd=2, xlab="W_range", ylab="PDF", main="Third Follow Up")
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
  as.vector(Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se)) * 
  as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) * 
  as.vector(Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se)) 

par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(range, prod13, type='l', col=2, main="points 1+3")
  plot(range, prod134, type='l', col=2, main="points 1+3+4")
  plot(range, prod123, type='l', col=2, main="points 1+2+3")
  plot(range, prod1234, type='l', col=2, main="points 1+2+3+4")

par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(range, log(prod13), type='l', col=3, main="points 1+3")
  plot(range, log(prod134), type='l', col=3, main="points 1+3+4")
  plot(range, log(prod123), type='l', col=3, main="points 1+2+3")
  plot(range, log(prod1234), type='l', col=3, main="points 1+2+3+4")

#Range of 0-10 ##############
range<-seq(from=0, to=10, by=1/1000)
  prod13<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * 
    as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) 
  prod134<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * 
    as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) * 
    as.vector(Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se)) 
  prod123<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * 
    as.vector(Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se)) * 
    as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) 
  prod1234<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * 
    as.vector(Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se)) * 
    as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) * 
    as.vector(Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se)) 

par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(range, prod13, type='l', col=2, main="points 1+3")
  plot(range, prod134, type='l', col=2, main="points 1+3+4")
  plot(range, prod123, type='l', col=2, main="points 1+2+3")
  plot(range, prod1234, type='l', col=2, main="points 1+2+3+4")

par(mfrow=c(2,2), oma=c(0,0,2,0))
  plot(range, log(prod13), type='l', col=3, main="points 1+3")
  plot(range, log(prod134), type='l', col=3, main="points 1+3+4")
  plot(range, log(prod123), type='l', col=3, main="points 1+2+3")
  plot(range, log(prod1234), type='l', col=3, main="points 1+2+3+4")

#Generate beta and lamda ranges to test #################
  beta_min<-5.0e-7
  beta_max<-5.0e-5
  
  beta_range<-seq(from=beta_min, to = beta_max, by = (beta_max-beta_min)/100)

  lamda_min<-5.0e-6
  lamda_max<-5.0e-4
  
  lamda_range<-seq(from=lamda_min, to = lamda_max, by = (lamda_max-lamda_min)/100)

#Generate matrices to fill #############  
  W_baseline_range<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  W_postmda_July12_range<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  W_postmda_Feb13_range<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  W_postmda_Sep13_range<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  
  loglikelihood_range<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range1<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range2<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range3<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range4<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range12<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range13<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range14<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range123<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range134<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  likelihood_range1234<-matrix(0, nrow = length(beta_range) , ncol = length(lamda_range))
  

  output_all<-numeric()
  #Timepoints of epi datapoint collection
    timepoints<-c(198, 362, 576, 760)
  #Estimates of measured worm burden at epi timepoints
    W_timepoints<-c(W_baseline, W_5month_field_mean, W_Feb13_field_mean, W_Sep13_field_mean)
    
    
#Generate data frame to fill ################
  fitting<-data.frame('beta' = rep(beta_range, times = length(lamda_range)),
                       'lamda' = rep(lamda_range, each = length(beta_range)),
                      'snail.prev' = 0,
                      'W1' = 0,
                      'W2' = 0,
                      'W3' = 0,
                      'W4' = 0,
                      'likelihood_range1' = 0,
                      'likelihood_range2' = 0,
                      'likelihood_range3' = 0,
                      'likelihood_range4' = 0,
                      'likelihood_range12' = 0,
                      'likelihood_range13' = 0,
                      'likelihood_range14' = 0,
                      'likelihood_range123' = 0,
                      'likelihood_range134' = 0,
                      'likelihood_range1234'= 0)  
#Run to estimate fits to epi data given tested range of beta and lamda values  #####################    
  for(i in 1:nrow(fitting)){
      
      output_all<-numeric()
      params["beta"]<-fitting[i,1]
      params["lamda"]<-fitting[i,2]
      params["mda"]<-0 #no mda
      k<-params['k']
      
      R0<-get_Ro_beta_lamda(muPq = p.dead, beta = params["beta"], lamda = params["lamda"])[3]
      
      if(R0<1){
        
        fitting[i,8] = 0
        fitting[i,9] = 0
        fitting[i,10] = 0
        fitting[i,11] = 0
        
        print(R0)
        
      } else { # run to equilibrium
              
      time<-seq(from=0, to=50*365, by=1)
      nstart=c(S=4000,E=2000,I=500, Wt=20, Wu=20, P=0)
      output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
      output_eqbm<-output[dim(output)[1],]
    
      cov<-params["cov"]
      
      W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
      
      if(W_eqbm < 1){
        
        fitting[i,8] = 0
        fitting[i,9] = 0
        fitting[i,10] = 0
        fitting[i,11] = 0
        
        print(W_eqbm)
        
      } else {
        plot(output$time, output$S, type='l', xlab="time",ylab="System Variables", 
           ylim=c(0,max( output$S+output$E+output$I )), col='blue', lwd=2)
        lines(output$time,output$E,col='orange', lwd=2)
        lines(output$time,output$I,col='red', lwd=2)
        lines(output$time,output$S+output$E+output$I,col='black', lwd=2)
        
      W1 <- output_eqbm$Wt
      
      if(W1 < 0){
        W1 = 0
      }
      
      W = W1
      
      gamma1 = (1 - ((1-(W1/(W1+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
       if(gamma1 < 0){
         gamma1 = 0
       }
            
      fitting[i,3]<-(output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100

      fitting[i,4] <- 0.5 * W1 * gamma1
      
      
      # apply mda
      nstart1<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt, 
                 Wu=output_eqbm$Wu, 
                 P=0)                   
      output_July12<-Senegal_mda_halstead(nstart1, params, mda_days) 
      output_eqbm<-output_July12[dim(output_July12)[1],]
      
      W2 = output_eqbm$Wt
      
      if(W2 < 0){
        W2 = 0
      }
      
      W = W2
      
      gamma2 = (1 - ((1-(W2/(W2+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
        if(gamma2 < 0){
          gamma2 = 0
        }
      
      fitting[i,5] <- 0.5 * W2 * gamma2
      
      output_all<-rbind(output_all, output_July12 ) #Joing dataframes for continuous time series
      
      #apply mda at first survey in July 2012 using mda flag and then run for 1 day
      time_last<-output_July12[dim(output_July12)[1],1]
      
      params["mda"]<-1
      time<-seq(time_last, time_last+1,1)
      nstart<-c(S=output_eqbm$S, 
                E=output_eqbm$E,
                I=output_eqbm$I, 
                Wt=output_eqbm$Wt, 
                Wu=output_eqbm$Wu, 
                P=0)
      
      output_mda_July12=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
      output_PZQ2<-output_mda_July12[dim(output_mda_July12)[1],]
      
      output_all<-rbind(output_all, output_mda_July12[-1,])
      time_last<-output_mda_July12[dim(output_mda_July12)[1],1]
      
      # apply no PZQ using mda flag and run for 7 months, July to Feb = (4*31)+(3*30) = 214 days.
      params["mda"]<-0
      time<-seq(time_last, time_last+214,1)
      nstart<-c(S=output_PZQ2$S, 
                E=output_PZQ2$E,
                I=output_PZQ2$I, 
                Wt=output_PZQ2$Wt, 
                Wu=output_PZQ2$Wu, 
                P=0)
      
      output_Feb13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
      output_Feb13_eqbm<-output_Feb13[dim(output_Feb13)[1],]
      
      W3 = output_Feb13_eqbm$Wt
      
      if(W3 < 0){
        W3 = 0
      }
      
      W = W3
      
      gamma3 = (1 - ((1-(W3/(W3+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
        if(gamma3 < 0){
          gamma3 = 0
        }
      
      fitting[i,6]<- 0.5 * W3 * gamma3
      output_all<-rbind(output_all, output_Feb13[-1,]) ##
      
      time_last<-output_Feb13[dim(output_Feb13)[1],1]
      
      #### apply mda at next survey in Feb 2013 using mda flag and then run for 1 day
      params["mda"]<-1
      time<-seq(time_last, time_last+1,1)
      nstart<-c(S=output_Feb13_eqbm$S, 
                E=output_Feb13_eqbm$E,
                I=output_Feb13_eqbm$I, 
                Wt=output_Feb13_eqbm$Wt, 
                Wu=output_Feb13_eqbm$Wu, 
                P=0)
      
      output_mda_Feb13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
      output_PZQ3<-output_mda_Feb13[dim(output_mda_Feb13)[1],]
      
      output_all<-rbind(output_all, output_mda_Feb13[-1,]) ##
      
      time_last<-output_mda_Feb13[dim(output_mda_Feb13)[1],1]
      
      # apply no PZQ using mda flag and run for 6 months, Feb to Aug 2013 = (4*31)+(2*30) = 184 days.
      params["mda"]<-0
      time<-seq(time_last, time_last+184,1)
      nstart<-c(S=output_PZQ3$S, 
                E=output_PZQ3$E,
                I=output_PZQ3$I, 
                Wt=output_PZQ3$Wt, 
                Wu=output_PZQ3$Wu, 
                P=0)
      output_Sep13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
      output_Sep13_eqbm<-output_Sep13[dim(output_Sep13)[1],]
      
      W4 = output_Sep13_eqbm$Wt
      
      if(W4 < 0){
        W4 = 0
      }
      
      W = W4
      
      gamma4 = (1 - ((1-(W4/(W4+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
        if(gamma4 < 0){
          gamma4 = 0
        }
      
      fitting[i,7]<- 0.5 * W4 * gamma4
      output_all<-rbind(output_all, output_Sep13[-1,]) ##
      
      time_last<-output_Sep13[dim(output_Sep13)[1],1]
      
      fitting[i,8]<-Prob_gaussian(y=fitting[i,4], mu=W_baseline, sd=W_baseline_se)  # only baseline
      fitting[i,9]<-Prob_gaussian(y=fitting[i,5], mu=W_5month_field_mean, sd=W_5month_field_se) # only data point 2
      fitting[i,10]<-Prob_gaussian(y=fitting[i,6], mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3
      fitting[i,11]<-Prob_gaussian(y=fitting[i,7], mu=W_Sep13_field_mean, sd=W_Sep13_field_se) #only data point 4
           
      print('complete')

      }
      
      
    }
      

      fitting[i,12]<-fitting[i,8] * fitting[i,9] #Fit to baseline and 2
      fitting[i,13]<-fitting[i,8] * fitting[i,10] #Fit to baseline and 3
      fitting[i,14]<-fitting[i,8] * fitting[i,1] #Fit to baseline and 4
      fitting[i,15]<-fitting[i,8] * fitting[i,9] * fitting[i,10] #Fit to baseline, 2 and 3
      fitting[i,16]<-fitting[i,8] * fitting[i,10] * fitting[i,11] #Fit to baseline, 3 and 4
      fitting[i,17]<-fitting[i,8] * fitting[i,9] * fitting[i,10] * fitting[i,11] #Fit to all epi data points
      
  }# end of i loop
  
#Save results, make subsets ################    
  write.csv(fitting, 'C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results.csv',
              row.names = F)
    
  ftng_rslt<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results.csv')  
  
  for(i in 1:nrow(ftng_rslt)){
    ftng_rslt[i,18] = get_Ro_beta_lamda(muPq = p.dead, beta = ftng_rslt[i,1], lamda = ftng_rslt[i,2])[1]
    ftng_rslt[i,19] = get_Ro_beta_lamda(muPq = p.dead, beta = ftng_rslt[i,1], lamda = ftng_rslt[i,2])[3]
    
  }
  
  colnames(ftng_rslt)[18:19]<-c('N_eq', 'R0')
  
  ftng_Ro1<-subset(ftng_rslt, R0 > 1)

  ftng_sp_15<-subset(ftng_rslt, snail.prev < 15 & snail.prev > 0)
  
  ftng_sp_10<-subset(ftng_rslt, snail.prev < 10 & snail.prev > 0)
  
  ftng_sp_10$likelihood_range1<- -log(ftng_sp_10$likelihood_range1 +1)
  ftng_sp_10$likelihood_range2<- -log(ftng_sp_10$likelihood_range2 +1)
  ftng_sp_10$likelihood_range3<- -log(ftng_sp_10$likelihood_range3 +1)
  ftng_sp_10$likelihood_range4<- -log(ftng_sp_10$likelihood_range4 +1)
  ftng_sp_10$likelihood_range12<- -log(ftng_sp_10$likelihood_range12 +1)
  ftng_sp_10$likelihood_range13<- -log(ftng_sp_10$likelihood_range13 +1)
  ftng_sp_10$likelihood_range1<- -log(ftng_sp_10$likelihood_range14 +1)
  ftng_sp_10$likelihood_range123<- -log(ftng_sp_10$likelihood_range123 +1)
  ftng_sp_10$likelihood_range134<- -log(ftng_sp_10$likelihood_range134 +1)
  ftng_sp_10$likelihood_range1234<- -log(ftng_sp_10$likelihood_range1234 +1)
  
  ftng_rslt_log<-ftng_rslt
  
  ftng_rslt_log$likelihood_range1<- -log(ftng_rslt$likelihood_range1 +1)
  ftng_rslt_log$likelihood_range2<- -log(ftng_rslt$likelihood_range2 +1)
  ftng_rslt_log$likelihood_range3<- -log(ftng_rslt$likelihood_range3 +1)
  ftng_rslt_log$likelihood_range4<- -log(ftng_rslt$likelihood_range4 +1)
  ftng_rslt_log$likelihood_range12<- -log(ftng_rslt$likelihood_range12 +1)
  ftng_rslt_log$likelihood_range13<- -log(ftng_rslt$likelihood_range13 +1)
  ftng_rslt_log$likelihood_range1<- -log(ftng_rslt$likelihood_range14 +1)
  ftng_rslt_log$likelihood_range123<- -log(ftng_rslt$likelihood_range123 +1)
  ftng_rslt_log$likelihood_range134<- -log(ftng_rslt$likelihood_range134 +1)
  ftng_rslt_log$likelihood_range1234<- -log(ftng_rslt$likelihood_range1234 +1)
  
#Get negative log likelihood results ##################
  beta1 = ftng_rslt_log$beta[ftng_rslt_log$likelihood_range1 == min(ftng_rslt_log$likelihood_range1)]
  beta2 = ftng_rslt_log$beta[ftng_rslt_log$likelihood_range2 == min(ftng_rslt_log$likelihood_range2)]
  beta3 = ftng_rslt_log$beta[ftng_rslt_log$likelihood_range3 == min(ftng_rslt_log$likelihood_range3)]
  beta4 = ftng_rslt_log$beta[ftng_rslt_log$likelihood_range4 == min(ftng_rslt_log$likelihood_range4)]
  beta12 = ftng_rslt_log$beta[ftng_rslt_log$likelihood_range12 == min(ftng_rslt_log$likelihood_range12)]
  beta13 = ftng_rslt_log$beta[ftng_rslt_log$likelihood_range13 == min(ftng_rslt_log$likelihood_range13)]
  beta14 = ftng_rslt_log$beta[ftng_rslt_log$likelihood_range14 == min(ftng_rslt_log$likelihood_range14)]
  beta123 = ftng_rslt_log$beta[ftng_rslt_log$likelihood_range123 == min(ftng_rslt_log$likelihood_range123)]
  beta134 = ftng_rslt_log$beta[ftng_rslt_log$likelihood_range134 == min(ftng_rslt_log$likelihood_range134)]
  beta1234 = ftng_rslt_log$beta[ftng_rslt_log$likelihood_range1234 == min(ftng_rslt_log$likelihood_range1234)]
  
  lamda1 = ftng_rslt_log$lamda[ftng_rslt_log$likelihood_range1 == min(ftng_rslt_log$likelihood_range1)]
  lamda2 = ftng_rslt_log$lamda[ftng_rslt_log$likelihood_range2 == min(ftng_rslt_log$likelihood_range2)]
  lamda3 = ftng_rslt_log$lamda[ftng_rslt_log$likelihood_range3 == min(ftng_rslt_log$likelihood_range3)]
  lamda4 = ftng_rslt_log$lamda[ftng_rslt_log$likelihood_range4 == min(ftng_rslt_log$likelihood_range4)]
  lamda12 = ftng_rslt_log$lamda[ftng_rslt_log$likelihood_range12 == min(ftng_rslt_log$likelihood_range12)]
  lamda13 = ftng_rslt_log$lamda[ftng_rslt_log$likelihood_range13 == min(ftng_rslt_log$likelihood_range13)]
  lamda14 = ftng_rslt_log$lamda[ftng_rslt_log$likelihood_range14 == min(ftng_rslt_log$likelihood_range14)]
  lamda123 = ftng_rslt_log$lamda[ftng_rslt_log$likelihood_range123 == min(ftng_rslt_log$likelihood_range123)]
  lamda134 = ftng_rslt_log$lamda[ftng_rslt_log$likelihood_range134 == min(ftng_rslt_log$likelihood_range134)]
  lamda1234 = ftng_rslt_log$lamda[ftng_rslt_log$likelihood_range1234 == min(ftng_rslt_log$likelihood_range1234)]
  
  
#Plotting to visualize results ############ 
  plot(ftng_rslt_log$lamda[ftng_rslt_log$beta == beta3],
       ftng_rslt_log$likelihood_range3[ftng_rslt_log$beta == beta3], type = 'l',
       xlab = 'lamda', ylab = 'neg log likelihood', lwd = 2)
  
  ggplot(ftng_sp_10, aes(x=lamda, y=beta, fill=R0))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt, aes(x=lamda, y=beta, fill=snail.prev))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt, aes(x=lamda, y=beta, fill=W1))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt, aes(x=lamda, y=beta, fill=R0))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt_log, aes(x=lamda, y=beta, fill=likelihood_range1))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt_log, aes(x=lamda, y=beta, fill=likelihood_range2))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt_log, aes(x=lamda, y=beta, fill=likelihood_range3))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt_log, aes(x=lamda, y=beta, fill=likelihood_range4))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt_log, aes(x=lamda, y=beta, fill=likelihood_range12))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt_log, aes(x=lamda, y=beta, fill=likelihood_range13))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt_log, aes(x=lamda, y=beta, fill=likelihood_range14))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt_log, aes(x=lamda, y=beta, fill=likelihood_range123))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt_log, aes(x=lamda, y=beta, fill=likelihood_range134))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
  ggplot(ftng_rslt_log, aes(x=lamda, y=beta, fill=likelihood_range1234))+
    theme_bw()+
    geom_tile(color='white', size=0.1)+
    scale_fill_continuous(low='green', high='red')
  
#Run model through mda simulations with likelihood results from high transmission#################
    test3<-data.frame('beta' = beta3,
                        'lamda' = lamda3,
                        'snail.prev' = 0,
                        'W1' = 0,
                        'W2' = 0,
                        'W3' = 0,
                        'W4' = 0,
                        'likelihood_range1' = 0,
                        'likelihood_range2' = 0,
                        'likelihood_range3' = 0,
                        'likelihood_range4' = 0,
                        'likelihood_range12' = 0,
                        'likelihood_range13' = 0,
                        'likelihood_range14' = 0,
                        'likelihood_range123' = 0,
                        'likelihood_range134' = 0,
                        'likelihood_range1234'= 0)  
    
    for(i in 1:nrow(test3)){
      
      output_all<-numeric()
      params["beta"]<-test3[i,1]
      params["lamda"]<-test3[i,2]
      params["mda"]<-0 #no mda
      k<-params['k']
      
      R0<-get_Ro_beta_lamda(muPq = p.dead, beta = params["beta"], lamda = params["lamda"])[3]
      
      if(R0<1){
        
        test3[i,8] = 0
        test3[i,9] = 0
        test3[i,10] = 0
        test3[i,11] = 0
        
        print(R0)
        
      } else { # run to equilibrium
        
        time<-seq(from=0, to=50*365, by=1)
        nstart=c(S=4000,E=2000,I=500, Wt=20, Wu=20, P=0)
        output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
        output_eqbm<-output[dim(output)[1],]
        
        cov<-params["cov"]
        
        W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
        
        if(W_eqbm < 1){
          
          test3[i,8] = 0
          test3[i,9] = 0
          test3[i,10] = 0
          test3[i,11] = 0
          
          print(W_eqbm)
          
        } else {
          plot(output$time, output$S, type='l', xlab="time",ylab="System Variables", 
               ylim=c(0,max( output$S+output$E+output$I )), col='blue', lwd=2)
          lines(output$time,output$E,col='orange', lwd=2)
          lines(output$time,output$I,col='red', lwd=2)
          lines(output$time,output$S+output$E+output$I,col='black', lwd=2)
          
          W1 <- output_eqbm$Wt
          
          if(W1 < 0){
            W1 = 0
          }
          
          W = W1
          
          gamma1 = (1 - ((1-(W1/(W1+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
          if(gamma1 < 0){
            gamma1 = 0
          }
          
          test3[i,3]<-(output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100
          
          test3[i,4] <- 0.5 * W1 * gamma1
          
          
          # apply mda
          nstart1<-c(S=output_eqbm$S, 
                     E=output_eqbm$E, 
                     I=output_eqbm$I, 
                     Wt=output_eqbm$Wt, 
                     Wu=output_eqbm$Wu, 
                     P=0)                   
          output_July12<-Senegal_mda_halstead(nstart1, params, mda_days) 
          output_eqbm<-output_July12[dim(output_July12)[1],]
          
          W2 = output_eqbm$Wt
          
          if(W2 < 0){
            W2 = 0
          }
          
          W = W2
          
          gamma2 = (1 - ((1-(W2/(W2+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
          if(gamma2 < 0){
            gamma2 = 0
          }
          
          test3[i,5] <- 0.5 * W2 * gamma2
          
          output_all<-rbind(output_all, output_July12 ) #Joing dataframes for continuous time series
          
          #apply mda at first survey in July 2012 using mda flag and then run for 1 day
          time_last<-output_July12[dim(output_July12)[1],1]
          
          params["mda"]<-1
          time<-seq(time_last, time_last+1,1)
          nstart<-c(S=output_eqbm$S, 
                    E=output_eqbm$E,
                    I=output_eqbm$I, 
                    Wt=output_eqbm$Wt, 
                    Wu=output_eqbm$Wu, 
                    P=0)
          
          output_mda_July12=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
          output_PZQ2<-output_mda_July12[dim(output_mda_July12)[1],]
          
          output_all<-rbind(output_all, output_mda_July12[-1,])
          time_last<-output_mda_July12[dim(output_mda_July12)[1],1]
          
          # apply no PZQ using mda flag and run for 7 months, July to Feb = (4*31)+(3*30) = 214 days.
          params["mda"]<-0
          time<-seq(time_last, time_last+214,1)
          nstart<-c(S=output_PZQ2$S, 
                    E=output_PZQ2$E,
                    I=output_PZQ2$I, 
                    Wt=output_PZQ2$Wt, 
                    Wu=output_PZQ2$Wu, 
                    P=0)
          
          output_Feb13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
          output_Feb13_eqbm<-output_Feb13[dim(output_Feb13)[1],]
          
          W3 = output_Feb13_eqbm$Wt
          
          if(W3 < 0){
            W3 = 0
          }
          
          W = W3
          
          gamma3 = (1 - ((1-(W3/(W3+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
          if(gamma3 < 0){
            gamma3 = 0
          }
          
          test3[i,6]<- 0.5 * W3 * gamma3
          output_all<-rbind(output_all, output_Feb13[-1,]) ##
          
          time_last<-output_Feb13[dim(output_Feb13)[1],1]
          
          #### apply mda at next survey in Feb 2013 using mda flag and then run for 1 day
          params["mda"]<-1
          time<-seq(time_last, time_last+1,1)
          nstart<-c(S=output_Feb13_eqbm$S, 
                    E=output_Feb13_eqbm$E,
                    I=output_Feb13_eqbm$I, 
                    Wt=output_Feb13_eqbm$Wt, 
                    Wu=output_Feb13_eqbm$Wu, 
                    P=0)
          
          output_mda_Feb13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
          output_PZQ3<-output_mda_Feb13[dim(output_mda_Feb13)[1],]
          
          output_all<-rbind(output_all, output_mda_Feb13[-1,]) ##
          
          time_last<-output_mda_Feb13[dim(output_mda_Feb13)[1],1]
          
          # apply no PZQ using mda flag and run for 6 months, Feb to Aug 2013 = (4*31)+(2*30) = 184 days.
          params["mda"]<-0
          time<-seq(time_last, time_last+184,1)
          nstart<-c(S=output_PZQ3$S, 
                    E=output_PZQ3$E,
                    I=output_PZQ3$I, 
                    Wt=output_PZQ3$Wt, 
                    Wu=output_PZQ3$Wu, 
                    P=0)
          output_Sep13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
          output_Sep13_eqbm<-output_Sep13[dim(output_Sep13)[1],]
          
          W4 = output_Sep13_eqbm$Wt
          
          if(W4 < 0){
            W4 = 0
          }
          
          W = W4
          
          gamma4 = (1 - ((1-(W4/(W4+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
          if(gamma4 < 0){
            gamma4 = 0
          }
          
          test3[i,7]<- 0.5 * W4 * gamma4
          output_all<-rbind(output_all, output_Sep13[-1,]) ##
          
          time_last<-output_Sep13[dim(output_Sep13)[1],1]
          
          test3[i,8]<-Prob_gaussian(y=test3[i,4], mu=W_baseline, sd=W_baseline_se)  # only baseline
          test3[i,9]<-Prob_gaussian(y=test3[i,5], mu=W_5month_field_mean, sd=W_5month_field_se) # only data point 2
          test3[i,10]<-Prob_gaussian(y=test3[i,6], mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3
          test3[i,11]<-Prob_gaussian(y=test3[i,7], mu=W_Sep13_field_mean, sd=W_Sep13_field_se) #only data point 4
          
          print('complete')
          
        }
        
        
      }
      
      
      test3[i,12]<-test3[i,8] * test3[i,9] #Fit to baseline and 2
      test3[i,13]<-test3[i,8] * test3[i,10] #Fit to baseline and 3
      test3[i,14]<-test3[i,8] * test3[i,1] #Fit to baseline and 4
      test3[i,15]<-test3[i,8] * test3[i,9] * test3[i,10] #Fit to baseline, 2 and 3
      test3[i,16]<-test3[i,8] * test3[i,10] * test3[i,11] #Fit to baseline, 3 and 4
      test3[i,17]<-test3[i,8] * test3[i,9] * test3[i,10] * test3[i,11] #Fit to all epi data points
      
    }# end of i loop
    
    plot(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), type='l', lwd = 2,
         xlab = 'time', ylab = '~Mean Worm burden')
    points(x = c(198, 362, 576, 760),
           y = c(W_baseline, W_5month_field_mean, W_Feb13_field_mean, W_Sep13_field_mean),
           pch = 1, col = 'red')
#Results from low transmission###############    
    test2<-data.frame('beta' = beta2,
                      'lamda' = lamda2,
                      'snail.prev' = 0,
                      'W1' = 0,
                      'W2' = 0,
                      'W3' = 0,
                      'W4' = 0,
                      'likelihood_range1' = 0,
                      'likelihood_range2' = 0,
                      'likelihood_range3' = 0,
                      'likelihood_range4' = 0,
                      'likelihood_range12' = 0,
                      'likelihood_range13' = 0,
                      'likelihood_range14' = 0,
                      'likelihood_range123' = 0,
                      'likelihood_range134' = 0,
                      'likelihood_range1234'= 0)  
    
    for(i in 1:nrow(test2)){
      
      output_all<-numeric()
      params["beta"]<-test2[i,1]
      params["lamda"]<-test2[i,2]
      params["mda"]<-0 #no mda
      k<-params['k']
      
      R0<-get_Ro_beta_lamda(muPq = p.dead, beta = params["beta"], lamda = params["lamda"])[3]
      
      if(R0<1){
        
        test2[i,8] = 0
        test2[i,9] = 0
        test2[i,10] = 0
        test2[i,11] = 0
        
        print(R0)
        
      } else { # run to equilibrium
        
        time<-seq(from=0, to=50*365, by=1)
        nstart=c(S=4000,E=2000,I=500, Wt=20, Wu=20, P=0)
        output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
        output_eqbm<-output[dim(output)[1],]
        
        cov<-params["cov"]
        
        W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
        
        if(W_eqbm < 1){
          
          test2[i,8] = 0
          test2[i,9] = 0
          test2[i,10] = 0
          test2[i,11] = 0
          
          print(W_eqbm)
          
        } else {
          plot(output$time, output$S, type='l', xlab="time",ylab="System Variables", 
               ylim=c(0,max( output$S+output$E+output$I )), col='blue', lwd=2)
          lines(output$time,output$E,col='orange', lwd=2)
          lines(output$time,output$I,col='red', lwd=2)
          lines(output$time,output$S+output$E+output$I,col='black', lwd=2)
          
          W1 <- output_eqbm$Wt
          
          if(W1 < 0){
            W1 = 0
          }
          
          W = W1
          
          gamma1 = (1 - ((1-(W1/(W1+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
          if(gamma1 < 0){
            gamma1 = 0
          }
          
          test2[i,3]<-(output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100
          
          test2[i,4] <- 0.5 * W1 * gamma1
          
          
          # apply mda
          nstart1<-c(S=output_eqbm$S, 
                     E=output_eqbm$E, 
                     I=output_eqbm$I, 
                     Wt=output_eqbm$Wt, 
                     Wu=output_eqbm$Wu, 
                     P=0)                   
          output_July12<-Senegal_mda_halstead(nstart1, params, mda_days) 
          output_eqbm<-output_July12[dim(output_July12)[1],]
          
          W2 = output_eqbm$Wt
          
          if(W2 < 0){
            W2 = 0
          }
          
          W = W2
          
          gamma2 = (1 - ((1-(W2/(W2+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
          if(gamma2 < 0){
            gamma2 = 0
          }
          
          test2[i,5] <- 0.5 * W2 * gamma2
          
          output_all<-rbind(output_all, output_July12 ) #Joing dataframes for continuous time series
          
          #apply mda at first survey in July 2012 using mda flag and then run for 1 day
          time_last<-output_July12[dim(output_July12)[1],1]
          
          params["mda"]<-1
          time<-seq(time_last, time_last+1,1)
          nstart<-c(S=output_eqbm$S, 
                    E=output_eqbm$E,
                    I=output_eqbm$I, 
                    Wt=output_eqbm$Wt, 
                    Wu=output_eqbm$Wu, 
                    P=0)
          
          output_mda_July12=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
          output_PZQ2<-output_mda_July12[dim(output_mda_July12)[1],]
          
          output_all<-rbind(output_all, output_mda_July12[-1,])
          time_last<-output_mda_July12[dim(output_mda_July12)[1],1]
          
          # apply no PZQ using mda flag and run for 7 months, July to Feb = (4*31)+(3*30) = 214 days.
          params["mda"]<-0
          time<-seq(time_last, time_last+214,1)
          nstart<-c(S=output_PZQ2$S, 
                    E=output_PZQ2$E,
                    I=output_PZQ2$I, 
                    Wt=output_PZQ2$Wt, 
                    Wu=output_PZQ2$Wu, 
                    P=0)
          
          output_Feb13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
          output_Feb13_eqbm<-output_Feb13[dim(output_Feb13)[1],]
          
          W3 = output_Feb13_eqbm$Wt
          
          if(W3 < 0){
            W3 = 0
          }
          
          W = W3
          
          gamma3 = (1 - ((1-(W3/(W3+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
          if(gamma3 < 0){
            gamma3 = 0
          }
          
          test2[i,6]<- 0.5 * W3 * gamma3
          output_all<-rbind(output_all, output_Feb13[-1,]) ##
          
          time_last<-output_Feb13[dim(output_Feb13)[1],1]
          
          #### apply mda at next survey in Feb 2013 using mda flag and then run for 1 day
          params["mda"]<-1
          time<-seq(time_last, time_last+1,1)
          nstart<-c(S=output_Feb13_eqbm$S, 
                    E=output_Feb13_eqbm$E,
                    I=output_Feb13_eqbm$I, 
                    Wt=output_Feb13_eqbm$Wt, 
                    Wu=output_Feb13_eqbm$Wu, 
                    P=0)
          
          output_mda_Feb13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
          output_PZQ3<-output_mda_Feb13[dim(output_mda_Feb13)[1],]
          
          output_all<-rbind(output_all, output_mda_Feb13[-1,]) ##
          
          time_last<-output_mda_Feb13[dim(output_mda_Feb13)[1],1]
          
          # apply no PZQ using mda flag and run for 6 months, Feb to Aug 2013 = (4*31)+(2*30) = 184 days.
          params["mda"]<-0
          time<-seq(time_last, time_last+184,1)
          nstart<-c(S=output_PZQ3$S, 
                    E=output_PZQ3$E,
                    I=output_PZQ3$I, 
                    Wt=output_PZQ3$Wt, 
                    Wu=output_PZQ3$Wu, 
                    P=0)
          output_Sep13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
          output_Sep13_eqbm<-output_Sep13[dim(output_Sep13)[1],]
          
          W4 = output_Sep13_eqbm$Wt
          
          if(W4 < 0){
            W4 = 0
          }
          
          W = W4
          
          gamma4 = (1 - ((1-(W4/(W4+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
          if(gamma4 < 0){
            gamma4 = 0
          }
          
          test2[i,7]<- 0.5 * W4 * gamma4
          output_all<-rbind(output_all, output_Sep13[-1,]) ##
          
          time_last<-output_Sep13[dim(output_Sep13)[1],1]
          
          test2[i,8]<-Prob_gaussian(y=test2[i,4], mu=W_baseline, sd=W_baseline_se)  # only baseline
          test2[i,9]<-Prob_gaussian(y=test2[i,5], mu=W_5month_field_mean, sd=W_5month_field_se) # only data point 2
          test2[i,10]<-Prob_gaussian(y=test2[i,6], mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3
          test2[i,11]<-Prob_gaussian(y=test2[i,7], mu=W_Sep13_field_mean, sd=W_Sep13_field_se) #only data point 4
          
          print('complete')
          
        }
        
        
      }
      
      
      test2[i,12]<-test2[i,8] * test2[i,9] #Fit to baseline and 2
      test2[i,13]<-test2[i,8] * test2[i,10] #Fit to baseline and 3
      test2[i,14]<-test2[i,8] * test2[i,1] #Fit to baseline and 4
      test2[i,15]<-test2[i,8] * test2[i,9] * test2[i,10] #Fit to baseline, 2 and 3
      test2[i,16]<-test2[i,8] * test2[i,10] * test2[i,11] #Fit to baseline, 3 and 4
      test2[i,17]<-test2[i,8] * test2[i,9] * test2[i,10] * test2[i,11] #Fit to all epi data points
      
    }# end of i loop
    
    plot(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), type='l', lwd = 2,
         xlab = 'time', ylab = '~Mean Worm burden')
    points(x = c(0, time_last2, time_last4, time_last5),
           y = c(W_baseline, W_5month_field_mean, W_Feb13_field_mean, W_Sep13_field_mean),
           pch = 1, col = 'red')
    
#Generate beta and lamda ranges to test #################
    beta_min2<-1.0e-7
    beta_max2<-1.0e-5
    
    beta_range2<-seq(from=beta_min2, to = beta_max2, by = (beta_max2-beta_min2)/50)
    
    lamda_min2<-1.0e-5
    lamda_max2<-1.0e-3
    
    lamda_range2<-seq(from=lamda_min2, to = lamda_max2, by = (lamda_max2-lamda_min2)/50)
  
    #Timepoints of epi datapoint collection
    timepoints<-c(198, 362, 576, 760)
    #Estimates of measured worm burden at epi timepoints
    W_timepoints<-c(W_baseline, W_5month_field_mean, W_Feb13_field_mean, W_Sep13_field_mean)
    
    
#Generate data frame to fill ################
    fitting2<-data.frame('beta' = rep(beta_range2, times = length(lamda_range2)),
                        'lamda' = rep(lamda_range2, each = length(beta_range2)),
                        'snail.prev' = 0,
                        'W1' = 0,
                        'W2' = 0,
                        'W3' = 0,
                        'W4' = 0,
                        'likelihood_range1' = 0,
                        'likelihood_range2' = 0,
                        'likelihood_range3' = 0,
                        'likelihood_range4' = 0,
                        'likelihood_range12' = 0,
                        'likelihood_range13' = 0,
                        'likelihood_range14' = 0,
                        'likelihood_range123' = 0,
                        'likelihood_range124' = 0,
                        'likelihood_range134' = 0,
                        'likelihood_range1234'= 0)  
    #Run to estimate fits to epi data given tested range of beta and lamda values  #####################    
    for(i in 1:nrow(fitting2)){
      
      print(i)
      
      output_all<-numeric()
      params["beta"]<-fitting2[i,1]
      params["lamda"]<-fitting2[i,2]
      params["mda"]<-0 #no mda
      k<-params['k']
      
      R0<-get_Ro_beta_lamda(muPq = p.dead, beta = params["beta"], lamda = params["lamda"])[3]
      
      if(R0<1){
        
        print(paste('R0=', R0, sep = ''))
        
      } else { # run to equilibrium
        
        time<-seq(from=0, to=50*365, by=1)
        nstart=c(S=4000,E=2000,I=500, Wt=20, Wu=20, P=0)
        output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
        output_eqbm<-output[dim(output)[1],]
        
        cov<-params["cov"]
        
        W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
        
      if (W_eqbm < .4){

          print(paste('W_eqbm=', W_eqbm, sep = ''))
          
        } else if((output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100 > 10){
          
          print('Snail prev > 10')
          
          
          
        } 
          else {
            
          plot(output$time, output$S, type='l', xlab="time",ylab="System Variables", 
               ylim=c(0,max( output$S+output$E+output$I )), col='blue', lwd=2)
          lines(output$time,output$E,col='orange', lwd=2)
          lines(output$time,output$I,col='red', lwd=2)
          lines(output$time,output$S+output$E+output$I,col='black', lwd=2)
          W1 <- output_eqbm$Wt
          
          if(W1 < 0){
            W1 = 0
          }
          
          W = W1
          
          gamma1 = (1 - ((1-(W1/(W1+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
          if(gamma1 < 0){
            gamma1 = 0
          }
          
          fitting2[i,3]<-(output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100
          
          fitting2[i,4] <- 0.5 * W1 * gamma1
          
          
          # apply mda
          nstart1<-c(S=output_eqbm$S, 
                     E=output_eqbm$E, 
                     I=output_eqbm$I, 
                     Wt=output_eqbm$Wt, 
                     Wu=output_eqbm$Wu, 
                     P=0)  
          
          output_July12<-Senegal_mda_halstead(nstart1, params, mda_days) 
          output_eqbm<-output_July12[dim(output_July12)[1],]
          
          W2 = output_eqbm$Wt
          
          if(W2 < 0){
            W2 = 0
          }
          
          W = W2
          
          gamma2 = (1 - ((1-(W2/(W2+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
          if(gamma2 < 0){
            gamma2 = 0
          }
          
          fitting2[i,5] <- 0.5 * W2 * gamma2
          
          output_all<-rbind(output_all, output_July12 ) #Joing dataframes for continuous time series
          
          #apply mda at first survey in July 2012 using mda flag and then run for 1 day
          time_last<-output_July12[dim(output_July12)[1],1]
          
          time<-seq(time_last, time_last+1,1)
          nstart<-c(S=output_eqbm$S, 
                    E=output_eqbm$E,
                    I=output_eqbm$I, 
                    Wt=output_eqbm$Wt, 
                    Wu=output_eqbm$Wu, 
                    P=0)
          
          output_mda_July12=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
          
          #Apply MDA to treated population (Wt) assuming 95% clearance of adult worms
          output_mda_July12[2,5]<-output_mda_July12[2,5] - output_mda_July12[2,5]*.95
          
          output_PZQ2<-output_mda_July12[dim(output_mda_July12)[1],]
          
          output_all<-rbind(output_all, output_mda_July12[-1,])
          time_last<-output_mda_July12[dim(output_mda_July12)[1],1]
          
          # apply no PZQ using mda flag and run for 7 months, July to Feb = (4*31)+(3*30) = 214 days.
          params["mda"]<-0
          time<-seq(time_last, time_last+214,1)
          nstart<-c(S=output_PZQ2$S, 
                    E=output_PZQ2$E,
                    I=output_PZQ2$I, 
                    Wt=output_PZQ2$Wt, 
                    Wu=output_PZQ2$Wu, 
                    P=0)
          
          output_Feb13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
          output_Feb13_eqbm<-output_Feb13[dim(output_Feb13)[1],]
          
          W3 = output_Feb13_eqbm$Wt
          
          if(W3 < 0){
            W3 = 0
          }
          
          W = W3
          
          gamma3 = (1 - ((1-(W3/(W3+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
          if(gamma3 < 0){
            gamma3 = 0
          }
          
          fitting2[i,6]<- 0.5 * W3 * gamma3
          output_all<-rbind(output_all, output_Feb13[-1,]) ##
          
          time_last<-output_Feb13[dim(output_Feb13)[1],1]
          
          #### apply mda at next survey in Feb 2013 using mda flag and then run for 1 day
          params["mda"]<-1
          time<-seq(time_last, time_last+1,1)
          nstart<-c(S=output_Feb13_eqbm$S, 
                    E=output_Feb13_eqbm$E,
                    I=output_Feb13_eqbm$I, 
                    Wt=output_Feb13_eqbm$Wt, 
                    Wu=output_Feb13_eqbm$Wu, 
                    P=0)
          
          output_mda_Feb13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
          
          #Apply MDA to treated population (Wt) assuming 95% clearance of adult worms
          output_mda_Feb13[2,5]<-output_mda_Feb13[2,5] - output_mda_Feb13[2,5]*.95
          
          output_PZQ3<-output_mda_Feb13[dim(output_mda_Feb13)[1],]
          
          output_all<-rbind(output_all, output_mda_Feb13[-1,]) ##
          
          time_last<-output_mda_Feb13[dim(output_mda_Feb13)[1],1]
          
          # apply no PZQ using mda flag and run for 6 months, Feb to Aug 2013 = (4*31)+(2*30) = 184 days.
          params["mda"]<-0
          time<-seq(time_last, time_last+184,1)
          nstart<-c(S=output_PZQ3$S, 
                    E=output_PZQ3$E,
                    I=output_PZQ3$I, 
                    Wt=output_PZQ3$Wt, 
                    Wu=output_PZQ3$Wu, 
                    P=0)
          output_Sep13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
          output_Sep13_eqbm<-output_Sep13[dim(output_Sep13)[1],]
          
          W4 = output_Sep13_eqbm$Wt
          
          if(W4 < 0){
            W4 = 0
          }
          
          W = W4
          
          gamma4 = (1 - ((1-(W4/(W4+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value))
          if(gamma4 < 0){
            gamma4 = 0
          }
          
          fitting2[i,7]<- 0.5 * W4 * gamma4
          output_all<-rbind(output_all, output_Sep13[-1,]) ##
          
          time_last<-output_Sep13[dim(output_Sep13)[1],1]
          
          fitting2[i,8]<-Prob_gaussian(y=fitting2[i,4], mu=W_baseline, sd=W_baseline_se)  # only baseline
          fitting2[i,9]<-Prob_gaussian(y=fitting2[i,5], mu=W_5month_field_mean, sd=W_5month_field_se) # only data point 2
          fitting2[i,10]<-Prob_gaussian(y=fitting2[i,6], mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3
          fitting2[i,11]<-Prob_gaussian(y=fitting2[i,7], mu=W_Sep13_field_mean, sd=W_Sep13_field_se) #only data point 4
          
          plot(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), type='l', lwd = 2,
               xlab = 'time', ylab = '~Mean Worm burden')
          points(x = c(198, 362, 576, 760),
                 y = c(W_baseline, W_5month_field_mean, W_Feb13_field_mean, W_Sep13_field_mean),
                 pch = 1, col = 'red')
          
      
      fitting2[i,12]<-fitting2[i,8] * fitting2[i,9] #Fit to baseline and 2
      fitting2[i,13]<-fitting2[i,8] * fitting2[i,10] #Fit to baseline and 3
      fitting2[i,14]<-fitting2[i,8] * fitting2[i,1] #Fit to baseline and 4
      fitting2[i,15]<-fitting2[i,8] * fitting2[i,9] * fitting2[i,10] #Fit to baseline, 2 and 3
      fitting2[i,16]<-fitting2[i,8] * fitting2[i,9] * fitting2[i,11] #Fit to baseline, 2 and 4
      fitting2[i,17]<-fitting2[i,8] * fitting2[i,10] * fitting2[i,11] #Fit to baseline, 3 and 4
      fitting2[i,18]<-fitting2[i,8] * fitting2[i,9] * fitting2[i,10] * fitting2[i,11] #Fit to all epi data points
      
          }
      
      } #Close final else statement
      
    }# end of i loop
    
#Save results, make subsets ############################    
write.csv(fitting2, 'C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results2.csv',
              row.names = F)
    
    ftng2<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results2.csv')  
    
    for(i in 1:nrow(ftng2)){
      ftng2[i,19] = get_Ro_beta_lamda(muPq = p.dead, beta = ftng2[i,1], lamda = ftng2[i,2])[1]
      ftng2[i,20] = get_Ro_beta_lamda(muPq = p.dead, beta = ftng2[i,1], lamda = ftng2[i,2])[3]
      
    }
    
    colnames(ftng2)[19:20]<-c('N_eq', 'R0')
    
    ftng2_Ro1<-subset(ftng2, R0 > 1)
    
    ftng2_log<-ftng2
    
    ftng2_log$likelihood_range1<- -log(ftng2$likelihood_range1 +1)
    ftng2_log$likelihood_range2<- -log(ftng2$likelihood_range2 +1)
    ftng2_log$likelihood_range3<- -log(ftng2$likelihood_range3 +1)
    ftng2_log$likelihood_range4<- -log(ftng2$likelihood_range4 +1)
    ftng2_log$likelihood_range12<- -log(ftng2$likelihood_range12 +1)
    ftng2_log$likelihood_range13<- -log(ftng2$likelihood_range13 +1)
    ftng2_log$likelihood_range1<- -log(ftng2$likelihood_range14 +1)
    ftng2_log$likelihood_range123<- -log(ftng2$likelihood_range123 +1)
    ftng2_log$likelihood_range134<- -log(ftng2$likelihood_range134 +1)
    ftng2_log$likelihood_range1234<- -log(ftng2$likelihood_range1234 +1)

#Get estimates for each reinfection period #################
    
    beta1 = ftng2_log$beta[ftng2_log$likelihood_range1 == min(ftng2_log$likelihood_range1)]
    beta2 = ftng2_log$beta[ftng2_log$likelihood_range2 == min(ftng2_log$likelihood_range2)]
    beta3 = ftng2_log$beta[ftng2_log$likelihood_range3 == min(ftng2_log$likelihood_range3)]
    beta4 = ftng2_log$beta[ftng2_log$likelihood_range4 == min(ftng2_log$likelihood_range4)]
    beta12 = ftng2_log$beta[ftng2_log$likelihood_range12 == min(ftng2_log$likelihood_range12)]
    beta13 = ftng2_log$beta[ftng2_log$likelihood_range13 == min(ftng2_log$likelihood_range13)]
    beta14 = ftng2_log$beta[ftng2_log$likelihood_range14 == min(ftng2_log$likelihood_range14)]
    beta123 = ftng2_log$beta[ftng2_log$likelihood_range123 == min(ftng2_log$likelihood_range123)]
    beta124 = ftng2_log$beta[ftng2_log$likelihood_range124 == min(ftng2_log$likelihood_range124)]
    beta134 = ftng2_log$beta[ftng2_log$likelihood_range134 == min(ftng2_log$likelihood_range134)]
    beta1234 = ftng2_log$beta[ftng2_log$likelihood_range1234 == min(ftng2_log$likelihood_range1234)]
    
    lamda1 = ftng2_log$lamda[ftng2_log$likelihood_range1 == min(ftng2_log$likelihood_range1)]
    lamda2 = ftng2_log$lamda[ftng2_log$likelihood_range2 == min(ftng2_log$likelihood_range2)]
    lamda3 = ftng2_log$lamda[ftng2_log$likelihood_range3 == min(ftng2_log$likelihood_range3)]
    lamda4 = ftng2_log$lamda[ftng2_log$likelihood_range4 == min(ftng2_log$likelihood_range4)]
    lamda12 = ftng2_log$lamda[ftng2_log$likelihood_range12 == min(ftng2_log$likelihood_range12)]
    lamda13 = ftng2_log$lamda[ftng2_log$likelihood_range13 == min(ftng2_log$likelihood_range13)]
    lamda14 = ftng2_log$lamda[ftng2_log$likelihood_range14 == min(ftng2_log$likelihood_range14)]
    lamda123 = ftng2_log$lamda[ftng2_log$likelihood_range123 == min(ftng2_log$likelihood_range123)]
    lamda124 = ftng2_log$lamda[ftng2_log$likelihood_range124 == min(ftng2_log$likelihood_range124)]
    lamda134 = ftng2_log$lamda[ftng2_log$likelihood_range134 == min(ftng2_log$likelihood_range134)]
    lamda1234 = ftng2_log$lamda[ftng2_log$likelihood_range1234 == min(ftng2_log$likelihood_range1234)]
    
    plot(ftng2_log$lamda[ftng2_log$beta == beta3],
         ftng2_log$likelihood_range3[ftng2_log$beta == beta3], type = 'l',
         xlab = 'lamda', ylab = 'neg log likelihood', lwd = 2)
    
#Get estimates for each reinfection period ##############
    best1<-subset(ftng2_log, beta == beta1 & lamda == lamda1)
    best2<-subset(ftng2_log, beta == beta2 & lamda == lamda2)
    best3<-subset(ftng2_log, beta == beta3 & lamda == lamda3)
    best4<-subset(ftng2_log, beta == beta4 & lamda == lamda4)
    
#Run simulation with lamda/beta parameters for each transmission season as above #############
    output_all<-numeric()
    params["beta"]<-beta1
    params["lamda"]<-lamda1
    params["mda"]<-0 #no mda
    k<-params['k']
    
    # run to equilibrium
      time<-seq(from=0, to=50*365, by=1)
      nstart=c(S=4000,E=2000,I=500, Wt=20, Wu=20, P=0)
      output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
      output_eqbm<-output[dim(output)[1],]
      
      cov<-params["cov"]
      
      W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
      
      snail.prev<-(output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100 
      
        
        plot(output$time, output$S, type='l', xlab="time",ylab="System Variables", 
             ylim=c(0,max( output$S+output$E+output$I )), col='blue', lwd=2)
        lines(output$time,output$E,col='orange', lwd=2)
        lines(output$time,output$I,col='red', lwd=2)
        lines(output$time,output$S+output$E+output$I,col='black', lwd=2)
        
        # apply mda
        params["beta"]<-beta2
        params["lamda"]<-lamda2
        
        nstart1<-c(S=output_eqbm$S, 
                   E=output_eqbm$E, 
                   I=output_eqbm$I, 
                   Wt=output_eqbm$Wt, 
                   Wu=output_eqbm$Wu, 
                   P=0)  
        
        output_July12<-Senegal_mda_halstead(nstart1, params, mda_days) 
        output_eqbm<-output_July12[dim(output_July12)[1],]
        
        
        output_all<-rbind(output_all, output_July12 ) #Join dataframes for continuous time series
        
        #apply mda at first survey in July 2012 using mda flag and then run for 1 day
        time_last<-output_July12[dim(output_July12)[1],1]
        
        time<-seq(time_last, time_last+1,1)
        nstart<-c(S=output_eqbm$S, 
                  E=output_eqbm$E,
                  I=output_eqbm$I, 
                  Wt=output_eqbm$Wt, 
                  Wu=output_eqbm$Wu, 
                  P=0)
        
        output_mda_July12=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
        
        #Apply MDA to treated population (Wt) assuming 95% clearance of adult worms
        output_mda_July12[2,5]<-output_mda_July12[2,5] - output_mda_July12[2,5]*.95
        
        output_PZQ2<-output_mda_July12[dim(output_mda_July12)[1],]
        
        output_all<-rbind(output_all, output_mda_July12[-1,])
        time_last<-output_mda_July12[dim(output_mda_July12)[1],1]
        
        # apply no PZQ using mda flag and run for 7 months, July to Feb = (4*31)+(3*30) = 214 days.
        params["mda"]<-0
        params["beta"]<-beta3
        params["lamda"]<-lamda3
        
        time<-seq(time_last, time_last+214,1)
        nstart<-c(S=output_PZQ2$S, 
                  E=output_PZQ2$E,
                  I=output_PZQ2$I, 
                  Wt=output_PZQ2$Wt, 
                  Wu=output_PZQ2$Wu, 
                  P=0)
        
        output_Feb13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
        output_Feb13_eqbm<-output_Feb13[dim(output_Feb13)[1],]
        
        output_all<-rbind(output_all, output_Feb13[-1,]) ##
        
        time_last<-output_Feb13[dim(output_Feb13)[1],1]
        
        #### apply mda at next survey in Feb 2013 using mda flag and then run for 1 day
        params["mda"]<-1
        params["beta"]<-beta4
        params["lamda"]<-lamda4
        
        time<-seq(time_last, time_last+1,1)
        nstart<-c(S=output_Feb13_eqbm$S, 
                  E=output_Feb13_eqbm$E,
                  I=output_Feb13_eqbm$I, 
                  Wt=output_Feb13_eqbm$Wt, 
                  Wu=output_Feb13_eqbm$Wu, 
                  P=0)
        
        output_mda_Feb13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
        
        #Apply MDA to treated population (Wt) assuming 95% clearance of adult worms
        output_mda_Feb13[2,5]<-output_mda_Feb13[2,5] - output_mda_Feb13[2,5]*.95
        
        output_PZQ3<-output_mda_Feb13[dim(output_mda_Feb13)[1],]
        
        output_all<-rbind(output_all, output_mda_Feb13[-1,]) ##
        
        time_last<-output_mda_Feb13[dim(output_mda_Feb13)[1],1]
        
        # apply no PZQ using mda flag and run for 6 months, Feb to Aug 2013 = (4*31)+(2*30) = 184 days.
        params["mda"]<-0
        time<-seq(time_last, time_last+184,1)
        nstart<-c(S=output_PZQ3$S, 
                  E=output_PZQ3$E,
                  I=output_PZQ3$I, 
                  Wt=output_PZQ3$Wt, 
                  Wu=output_PZQ3$Wu, 
                  P=0)
        output_Sep13=as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params)) 
        output_Sep13_eqbm<-output_Sep13[dim(output_Sep13)[1],]
        
        
        output_all<-rbind(output_all, output_Sep13[-1,]) ##
        
        time_last<-output_Sep13[dim(output_Sep13)[1],1]
        
        plot(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), type='l', lwd = 2,
             xlab = 'time', ylab = '~Mean Worm burden')
        points(x = c(198, 362, 576, 760),
               y = c(W_baseline, W_5month_field_mean, W_Feb13_field_mean, W_Sep13_field_mean),
               pch = 1, col = 'red')
        
    
    ggplot(ftng2_log, aes(x=lamda, y=beta, fill=likelihood_range1))+
      theme_bw()+
      geom_tile(color='white', size=0.1)+
      scale_fill_continuous(low='green', high='red')
    