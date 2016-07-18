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

time=seq(0,365*50,1) #50 years to ensure equilibrium is reached
nstart=c(S=4000,E=2000,I=500, Wt=10, Wu=10, P=0)

#Prevalence function given W and k ##################################

k<-0.1 # clumping parameter of the negative binomial distribution
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

#Get distributions for Epi datapoints given estimated W and k #################
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

Prob_gaussian<-function(y, mu, sd){
  p<-(1/(sqrt(2*pi)*sd) ) * exp( -((y-mu)^2)/(2*sd^2) )
  p
}

#Generate beta and lamda ranges to test #################
beta_min3<-1.0e-7
beta_max3<-1.0e-5

beta_range3<-seq(from=beta_min3, to = beta_max3, by = (beta_max3-beta_min3)/75)

lamda_min3<-1.0e-6
lamda_max3<-1.25e-3

lamda_range3<-seq(from=lamda_min3, to = lamda_max3, by = (lamda_max3-lamda_min3)/75)

lamda_min3.2<-1.0e-5
lamda_max3.2<-1.0e-3

lamda_range3.2<-seq(from=lamda_min3.2, to = lamda_max3.2, by = (lamda_max3.2-lamda_min3.2)/50)
#Prepare data frames and vectors ####################
#Timepoints of epi datapoint collection
timepoints<-c(198, 362, 576, 760)
#Estimates of measured worm burden at epi timepoints
W_timepoints<-c(W_baseline, W_5month_field_mean, W_Feb13_field_mean, W_Sep13_field_mean)

fitting1<-data.frame('beta1' = rep(beta_range3, times = length(lamda_range3)),
                     'lamda1' = rep(lamda_range3, each = length(beta_range3)),
                     'R0' = 0,
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
                     'likelihood_range1234' = 0)

fitting3<-data.frame('beta1' = rep(beta_range3, times = length(lamda_range3)),
                     'lamda1' = rep(lamda_range3, each = length(beta_range3)),
                     'lamda2' = rep(lamda_range3*10, each = length(beta_range3)),
                     'snail.prev' = 0,
                     'W1' = 0,
                     'W2' = 0,
                     'W3' = 0,
                     'likelihood_range1' = 0,
                     'likelihood_range2' = 0,
                     'likelihood_range3' = 0,
                     'likelihood_range12' = 0,
                     'likelihood_range13' = 0,
                     'likelihood_range123' = 0)
#Run to estimate fits to epi data from single lamda and beta values ####################
for(i in 1:nrow(fitting1)){
  
  print(i)
  
  output_all<-numeric()
  params["beta"]<-fitting1[i,1]
  params["lamda"]<-fitting1[i,2]
  k<-params['k']
  
  R0<-get_Ro_beta_lamda(muPq = p.dead, beta = fitting1[i,1], lamda = fitting1[i,2])[3]
  
  fitting1[i,3] = R0
  
  if(R0 < 1){
    print(paste('R0=', R0, sep = ''))
  } else {
    
    # run to equilibrium
    
    time<-seq(from=0, to=50*365, by=1)
    nstart=c(S=8000,E=0,I=0, Wt=10, Wu=10, P=0)
    output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
    output_eqbm<-output[dim(output)[1],]
    output_eqbm$time<-0
    
    cov<-params["cov"]
    
    W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
    
    if (W_eqbm < .4){
      
      print(paste('W_eqbm=', W_eqbm, sep = ''))
      
    } else if (output_eqbm$E < 0){
      
      print(paste('E_eqbm=', output_eqbm$E, sep = ''))
      
    } else if (output_eqbm$S > output$S[10000]+100){
      
      print('Eqbm not reached')
      
      plot(output$time, output$S, type='l', xlab="time",ylab="System Variables", 
           ylim=c(0,max( output$S+output$E+output$I )), col='blue', lwd=2)
      lines(output$time,output$E,col='orange', lwd=2)
      lines(output$time,output$I,col='red', lwd=2)
      lines(output$time,output$S+output$E+output$I,col='black', lwd=2)
      
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
      
      gamma1 = (1 - ((1-(W1/(W1+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma1 < 0){
        gamma1 = 0
      }
      
      fitting1[i,4]<-(output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100
      
      fitting1[i,5] <- 0.5 * W1 * gamma1
      
      output_all<-rbind(output_all, output_eqbm ) #Join dataframes for continuous time series
      
      #Apply first MDA at t=1
      
      nstart1<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=1, to=21, by=1)
      
      output_3wks<-as.data.frame(ode(nstart1,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_3wks ) #Join dataframes for continuous time series
      
      output_eqbm<-output_3wks[dim(output_3wks)[1],]
      
      #Apply second MDA at t=21
      
      nstart2<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      #Run for 5 months (150 days)
      time = seq(from=22, to=172, by=1)
      
      output2<-as.data.frame(ode(nstart2,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output2 ) #Join dataframes for continuous time series
      
      output_eqbm<-output2[dim(output2)[1],]
      
      W2 = output_eqbm$Wt
      
      if(W2 < 0){
        W2 = 0
      }
      
      W = W2
      
      gamma2 = (1 - ((1-(W2/(W2+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma2 < 0){
        gamma2 = 0
      }
      
      fitting1[i,6] <- 0.5 * W2 * gamma2
      
      # Apply MDA, run for 7 months (210 days)
      time<-seq(from=173, t=383, by=1)
      
      nstart3<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, 
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output3=as.data.frame(ode(nstart3,time,schisto_halstead_2pops_mda,params)) 
      
      output_all<-rbind(output_all, output3 ) #Join dataframes for continuous time series
      
      output_eqbm<-output3[dim(output3)[1],]
      
      W3 = output_eqbm$Wt
      
      if(W3 < 0){
        W3 = 0
      }
      
      W = W3
      
      gamma3 = (1 - ((1-(W3/(W3+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma3 < 0){
        gamma3 = 0
      }
      
      fitting1[i,7]<- 0.5 * W3 * gamma3
      
      # Apply MDA, run for 6 months reinfection (180 days)
      time<-seq(from=384, t=564, by=1)
      
      nstart4<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, 
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output4=as.data.frame(ode(nstart4,time,schisto_halstead_2pops_mda,params)) 
      
      output_all<-rbind(output_all, output4 ) #Join dataframes for continuous time series
      
      output_eqbm<-output4[dim(output4)[1],]
      
      W4 = output_eqbm$Wt
      
      if(W4 < 0){
        W4 = 0
      }
      
      W = W4
      
      gamma4 = (1 - ((1-(W4/(W4+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma4 < 0){
        gamma4 = 0
      }
      
      fitting1[i,8]<- 0.5 * W4 * gamma4
      
      fitting1[i,9]<-Prob_gaussian(y=fitting1[i,5], mu=W_baseline, sd=W_baseline_se)  # only baseline
      fitting1[i,10]<-Prob_gaussian(y=fitting1[i,6], mu=W_5month_field_mean, sd=W_5month_field_se) # only data point 2
      fitting1[i,11]<-Prob_gaussian(y=fitting1[i,7], mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3
      fitting1[i,12]<-Prob_gaussian(y=fitting1[i,8], mu=W_Sep13_field_mean, sd=W_Sep13_field_se) #only data point 3
      
      plot(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), type='l', lwd = 2,
           xlab = 'time', ylab = '~Mean Worm burden', ylim = c(0,45))
      points(x = timepoints - 198, y = W_timepoints,
             pch = 1, col = 'red')
      
      
      fitting1[i,13]<-fitting1[i,9] * fitting1[i,10] #Fit to baseline and 2
      fitting1[i,14]<-fitting1[i,9] * fitting1[i,11] #Fit to baseline and 3
      fitting1[i,15]<-fitting1[i,9] * fitting1[i,12] #Fit to baseline and 4
      fitting1[i,16]<-fitting1[i,9] * fitting1[i,10] * fitting1[i,11] #Fit to baseline, 2 and 3
      fitting1[i,17]<-fitting1[i,9] * fitting1[i,11] * fitting1[i,12] #Fir to baseline, 3 and 4 
      fitting1[i,18]<-fitting1[i,9] * fitting1[i,10] * fitting1[i,11] * fitting1[i,12] #Fit to baseline, 2, 3 and 4
      
      
    }
    
  }
  
} # end of i loop

#Save results, make subsets ############################    
write.csv(fitting1, 'C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results_beta_lamda.csv',
          row.names = F)

ftng1<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results_beta_lamda.csv')  

ftng1$likelihood_range1<- -log(ftng1$likelihood_range1 +1)
ftng1$likelihood_range2<- -log(ftng1$likelihood_range2 +1)
ftng1$likelihood_range3<- -log(ftng1$likelihood_range3 +1)
ftng1$likelihood_range4<- -log(ftng1$likelihood_range4 +1)

ftng1$likelihood_range12<- -log(ftng1$likelihood_range12 +1)
ftng1$likelihood_range13<- -log(ftng1$likelihood_range13 +1)
ftng1$likelihood_range14<- -log(ftng1$likelihood_range14 +1)

ftng1$likelihood_range123<- -log(ftng1$likelihood_range123 +1)
ftng1$likelihood_range134<- -log(ftng1$likelihood_range134 +1)

ftng1$likelihood_range1234<- -log(ftng1$likelihood_range1234 +1)

#Get estimates for each reinfection period #################

ftng1.1<-subset(ftng1, likelihood_range1 < 0 & snail.prev <10)
  unique(ftng1.1$beta1)

ftng1.2<-subset(ftng1, likelihood_range2 < 0 & snail.prev <10)
  unique(ftng1.2$beta1)
  
ftng1.3<-subset(ftng1, likelihood_range3 < 0 & snail.prev <10)
  unique(ftng1.3$beta1)

ftng1.4<-subset(ftng1, likelihood_range4 < 0 & snail.prev <10)
  unique(ftng1.4$beta1)


#Get estimates for each reinfection period #################

beta1 = ftng1$beta1[ftng1$likelihood_range1 == min(ftng1$likelihood_range1)]
beta2 = ftng1$beta1[ftng1$likelihood_range2 == min(ftng1$likelihood_range2)]
beta3 = ftng1$beta1[ftng1$likelihood_range3 == min(ftng1$likelihood_range3)]
beta4 = ftng1$beta1[ftng1$likelihood_range4 == min(ftng1$likelihood_range4)]

beta12 = ftng1$beta1[ftng1$likelihood_range12 == min(ftng1$likelihood_range12)]
beta13 = ftng1$beta1[ftng1$likelihood_range13 == min(ftng1$likelihood_range13)]
beta14 = ftng1$beta1[ftng1$likelihood_range14 == min(ftng1$likelihood_range14)]

beta123 = ftng1$beta1[ftng1$likelihood_range123 == min(ftng1$likelihood_range123)]
beta134 = ftng1$beta1[ftng1$likelihood_range134 == min(ftng1$likelihood_range134)]

lamda1 = ftng1$lamda1[ftng1$likelihood_range1 == min(ftng1$likelihood_range1)]
lamda2 = ftng1$lamda1[ftng1$likelihood_range2 == min(ftng1$likelihood_range2)]
lamda3 = ftng1$lamda1[ftng1$likelihood_range3 == min(ftng1$likelihood_range3)]
lamda4 = ftng1$lamda1[ftng1$likelihood_range4 == min(ftng1$likelihood_range4)]

lamda12 = ftng1$lamda1[ftng1$likelihood_range12 == min(ftng1$likelihood_range12)]
lamda13 = ftng1$lamda1[ftng1$likelihood_range13 == min(ftng1$likelihood_range13)]
lamda14 = ftng1$lamda1[ftng1$likelihood_range14 == min(ftng1$likelihood_range14)]

lamda123 = ftng1$lamda1[ftng1$likelihood_range123 == min(ftng1$likelihood_range123)]
lamda134 = ftng1$lamda1[ftng1$likelihood_range134 == min(ftng1$likelihood_range134)]

pt2<-list('beta' = beta2,
          'lamda' = lamda2,
          's.prev' = ftng1$snail.prev[ftng1$likelihood_range2 == min(ftng1$likelihood_range2)],
          'R0' = get_Ro_beta_lamda(muPq = p.dead, beta = beta2, lamda = lamda2))

plot(x = ftng1$lamda1[ftng1$beta1 == pt2$beta],
     y = ftng1$likelihood_range2[ftng1$beta1 == pt2$beta],
     xlab = 'lamda', ylab = '-log likelihood', type = 'l', main = 'Follow up 1')

pt3<-list('beta' = beta3,
          'lamda' = lamda3,
          's.prev' = ftng1$snail.prev[ftng1$likelihood_range3 == min(ftng1$likelihood_range3)],
          'R0' = get_Ro_beta_lamda(muPq = p.dead, beta = beta3, lamda = lamda3))

plot(x = ftng1$lamda1[ftng1$beta1 == pt3$beta],
     y = ftng1$likelihood_range3[ftng1$beta1 == pt3$beta],
     xlab = 'lamda', ylab = '-log likelihood', type = 'l', main = 'Follow up 2')

pt4<-list('beta' = beta4,
          'lamda' = lamda4,
          's.prev' = ftng1$snail.prev[ftng1$likelihood_range4 == min(ftng1$likelihood_range4)],
          'R0' = get_Ro_beta_lamda(muPq = p.dead, beta = beta4, lamda = lamda4))

plot(x = ftng1$lamda1[ftng1$beta1 == pt4$beta],
     y = ftng1$likelihood_range4[ftng1$beta1 == pt4$beta],
     xlab = 'lamda', ylab = '-log likelihood', type = 'l', main = 'Follow up 3')



#Run to estimate fits to epi data given tested range of beta and two lamda values  #####################    
for(i in 1:nrow(fitting3)){
  
  print(i)
  
  output_all<-numeric()
  params["beta"]<-fitting3[i,1]
  params["lamda"]<-fitting3[i,2]
  k<-params['k']
  
  R0<-get_Ro_beta_lamda(muPq = p.dead, beta = fitting3[i,1], lamda = fitting3[i,2])[3]
  
  if(R0 < 1){
    print(paste('R0=', R0, sep = ''))
  } else {
  
   # run to equilibrium
    
    time<-seq(from=0, to=50*365, by=1)
    nstart=c(S=8000,E=0,I=0, Wt=10, Wu=10, P=0)
    output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
    output_eqbm<-output[dim(output)[1],]
    output_eqbm$time<-0
    
    cov<-params["cov"]
    
    W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
    
    if (W_eqbm < .4){
      
      print(paste('W_eqbm=', W_eqbm, sep = ''))
      
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
      
      gamma1 = (1 - ((1-(W1/(W1+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma1 < 0){
        gamma1 = 0
      }
      
      fitting3[i,4]<-(output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100
      
      fitting3[i,5] <- 0.5 * W1 * gamma1
      
      output_all<-rbind(output_all, output_eqbm ) #Join dataframes for continuous time series
      
  #Apply first MDA at t=1
      
      nstart1<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=1, to=21, by=1)
      
      output_3wks<-as.data.frame(ode(nstart1,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_3wks ) #Join dataframes for continuous time series
      
      output_eqbm<-output_3wks[dim(output_3wks)[1],]
      
  #Apply second MDA at t=21 and then run until rainy season (june)
      
      nstart2<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=22, to=121, by=1)
      
      output_june<-as.data.frame(ode(nstart2,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_june ) #Join dataframes for continuous time series
      
      output_eqbm<-output_june[dim(output_june)[1],]
      
  #Change lamda value to correspond to high transmission season, run for 1 month
      params['lamda']<-fitting3[i,3]
      
      nstart3<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=122, to=151, by=1)
      
      output_july<-as.data.frame(ode(nstart3,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_july ) #Join dataframes for continuous time series
      
      output_eqbm<-output_july[dim(output_july)[1],]
      
      W2 = output_eqbm$Wt
      
      if(W2 < 0){
        W2 = 0
      }
      
      W = W2
      
      gamma2 = (1 - ((1-(W2/(W2+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma2 < 0){
        gamma2 = 0
      }
      
      fitting3[i,6] <- 0.5 * W2 * gamma2
      
      # Apply MDA and run through rainy season (end of october)
      time<-seq(from=152, t=273, by=1)
      
      nstart4<-c(S=output_eqbm$S, 
                E=output_eqbm$E,
                I=output_eqbm$I, 
                Wt=output_eqbm$Wt * 0.05, 
                Wu=output_eqbm$Wu, 
                P=0)
      
      output_oct=as.data.frame(ode(nstart4,time,schisto_halstead_2pops_mda,params)) 
      
      output_all<-rbind(output_all, output_oct ) #Join dataframes for continuous time series
      
      output_eqbm<-output_oct[dim(output_oct)[1],]
      
      #End of rainy season, go back to low lamda and run until Feb (end of year)
      time<-seq(274, 365, 1)
      params['lamda']<-fitting3[i,2]
      
      nstart5<-c(S=output_eqbm$S, 
                E=output_eqbm$E,
                I=output_eqbm$I, 
                Wt=output_eqbm$Wt, 
                Wu=output_eqbm$Wu, 
                P=0)
      
      output_Feb13=as.data.frame(ode(nstart5,time,schisto_halstead_2pops_mda,params)) 
      output_Feb13_eqbm<-output_Feb13[dim(output_Feb13)[1],]
      
      W3 = output_Feb13_eqbm$Wt
      
      if(W3 < 0){
        W3 = 0
      }
      
      W = W3
      
      gamma3 = (1 - ((1-(W3/(W3+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma3 < 0){
        gamma3 = 0
      }
      
      fitting3[i,7]<- 0.5 * W3 * gamma3
      output_all<-rbind(output_all, output_Feb13) ##
      
      fitting3[i,8]<-Prob_gaussian(y=fitting3[i,4], mu=W_baseline, sd=W_baseline_se)  # only baseline
      fitting3[i,9]<-Prob_gaussian(y=fitting3[i,5], mu=W_5month_field_mean, sd=W_5month_field_se) # only data point 2
      fitting3[i,10]<-Prob_gaussian(y=fitting3[i,6], mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3

      plot(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), type='l', lwd = 2,
           xlab = 'time', ylab = '~Mean Worm burden', ylim = c(0,45))
      points(x = c(0, 150, 364),
             y = c(W_baseline, W_5month_field_mean, W_Feb13_field_mean),
             pch = 1, col = 'red')
      
      
      fitting3[i,11]<-fitting3[i,8] * fitting3[i,9] #Fit to baseline and 2
      fitting3[i,12]<-fitting3[i,8] * fitting3[i,10] #Fit to baseline and 3
      fitting3[i,13]<-fitting3[i,8] * fitting3[i,9] * fitting3[i,10]#Fit to baseline 2 and 3
     
    }
    
  }
    
  } # end of i loop

#Save results, make subsets ############################    
write.csv(fitting3, 'C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results_2lamda.csv',
          row.names = F)

ftng3<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results_2lamda.csv')  

ftng3$likelihood_range1<- -log(ftng3$likelihood_range1 +1)
ftng3$likelihood_range2<- -log(ftng3$likelihood_range2 +1)
ftng3$likelihood_range3<- -log(ftng3$likelihood_range3 +1)
ftng3$likelihood_range12<- -log(ftng3$likelihood_range12 +1)
ftng3$likelihood_range13<- -log(ftng3$likelihood_range13 +1)
ftng3$likelihood_range123<- -log(ftng3$likelihood_range123 +1)

for(i in nrow(ftng3)){
  ftng3[i,14]<-get_Ro_beta_lamda(muPq = p.dead, beta = ftng3[i,1], lamda = ftng3[i,2])[3]
}

ftng3.gd<-subset(ftng3, snail.prev < 10)

  

#Get estimates for each reinfection period #################

beta1 = ftng3$beta[ftng3$likelihood_range1 == min(ftng3$likelihood_range1)]
beta2 = ftng3$beta[ftng3$likelihood_range2 == min(ftng3$likelihood_range2)]
beta3 = ftng3$beta[ftng3$likelihood_range3 == min(ftng3$likelihood_range3)]
beta12 = ftng3$beta[ftng3$likelihood_range12 == min(ftng3$likelihood_range12)]
beta13 = ftng3$beta[ftng3$likelihood_range13 == min(ftng3$likelihood_range13)]
beta123 = ftng3$beta[ftng3$likelihood_range123 == min(ftng3$likelihood_range123)]

lamda1 = ftng3$lamda1[ftng3$likelihood_range1 == min(ftng3$likelihood_range1)]
lamda2 = ftng3$lamda1[ftng3$likelihood_range2 == min(ftng3$likelihood_range2)]
lamda3 = ftng3$lamda1[ftng3$likelihood_range3 == min(ftng3$likelihood_range3)]
lamda12 = ftng3$lamda1[ftng3$likelihood_range12 == min(ftng3$likelihood_range12)]
lamda13 = ftng3$lamda1[ftng3$likelihood_range13 == min(ftng3$likelihood_range13)]
lamda123 = ftng3$lamda1[ftng3$likelihood_range123 == min(ftng3$likelihood_range123)]

lamda1.2 = ftng3$lamda2[ftng3$likelihood_range1 == min(ftng3$likelihood_range1)]
lamda2.2 = ftng3$lamda2[ftng3$likelihood_range2 == min(ftng3$likelihood_range2)]
lamda3.2 = ftng3$lamda2[ftng3$likelihood_range3 == min(ftng3$likelihood_range3)]
lamda12.2 = ftng3$lamda2[ftng3$likelihood_range12 == min(ftng3$likelihood_range12)]
lamda13.2 = ftng3$lamda2[ftng3$likelihood_range13 == min(ftng3$likelihood_range13)]
lamda123.2 = ftng3$lamda2[ftng3$likelihood_range123 == min(ftng3$likelihood_range123)]

#Generate beta and lamda ranges to test and data frame to fill#################
beta_min4<-5.0e-6
beta_max4<-5.0e-5

beta1s<-round(runif(5000, max = beta_max4, min = beta_min4), digits = 8)

lamda_min4<-1.0e-5
lamda_max4<-3.0e-5

lamda1s<-round(runif(5000, max = lamda_max4, min = lamda_min4), digits = 8)

lamda_min4.2<-3.0e-4
lamda_max4.2<-1.5e-3

lamda2s<-round(runif(5000, max = lamda_max4.2, min = lamda_min4.2), digits = 8)

fitting4<-data.frame('beta1' = beta1s,
                     'lamda1' = lamda1s,
                     'lamda2' = lamda2s,
                     'R0' = 0,
                     'snail.prev' = 0,
                     'W1' = 0,
                     'W2' = 0,
                     'W3' = 0,
                     'W4' = 0,
                     'likelihood_range1' = 0,
                     'likelihood_range2' = 0,
                     'likelihood_range3' = 0,
                     'likelihood_range4' = 0,
                     'likelihood_range234' = 0,
                     'likelihood_range1234' = 0)

hist(fitting4$beta1)
hist(fitting4$lamda1)
hist(fitting4$lamda2)



#Run to estimate fits to epi data given tested range of beta and two lamda values  ##################### 

for(i in 1:nrow(fitting4)){
  
  print(i)
  
  output_all<-numeric()
  params["beta"]<-fitting4[i,1]
  params["lamda"]<-fitting4[i,2]
  k<-params['k']
  
  R0<-get_Ro_beta_lamda(muPq = p.dead, beta = fitting4[i,1], lamda = fitting4[i,2])[3]
  
  fitting4[i,4] = R0 
  
  if(R0 < 1){
    print(paste('R0=', R0, sep = ''))
  } else {
    
    # run to equilibrium
    
    time<-seq(from=0, to=50*365, by=1)
    nstart=c(S=8000,E=0,I=0, Wt=10, Wu=10, P=0)
    output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
    output_eqbm<-output[dim(output)[1],]
    output_eqbm$time<-0
    
    cov<-params["cov"]
    
    fitting4[i,5]<-(output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100
    
    W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
    
    if (W_eqbm < .4){
      
      print(paste('W_eqbm=', W_eqbm, sep = ''))
      
    } else if (output_eqbm$E < 0){
      
      print(paste('E_eqbm=', output_eqbm$E, sep = ''))
      
    } else if (output_eqbm$S > output$S[10000]+100){
      
      print('Eqbm not reached')
      
      plot(output$time, output$S, type='l', xlab="time",ylab="System Variables", 
           ylim=c(0,max( output$S+output$E+output$I )), col='blue', lwd=2)
      lines(output$time,output$E,col='orange', lwd=2)
      lines(output$time,output$I,col='red', lwd=2)
      lines(output$time,output$S+output$E+output$I,col='black', lwd=2)
      
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
      
      gamma1 = (1 - ((1-(W1/(W1+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma1 < 0){
        gamma1 = 0
      }
      
      fitting4[i,6] <- 0.5 * W1 * gamma1
      
      output_all<-rbind(output_all, output_eqbm ) #Join dataframes for continuous time series
      
      #Apply first MDA at t=1
      
      nstart1<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=1, to=21, by=1)
      
      output_3wks<-as.data.frame(ode(nstart1,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_3wks ) #Join dataframes for continuous time series
      
      output_eqbm<-output_3wks[dim(output_3wks)[1],]
      
      #Apply second MDA at t=21 and then run until rainy season (june)
      
      nstart2<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply second MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=22, to=121, by=1)
      
      output_june<-as.data.frame(ode(nstart2,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_june ) #Join dataframes for continuous time series
      
      output_eqbm<-output_june[dim(output_june)[1],]
      
      #Change lamda value to correspond to high transmission season, run for 1 month
      params['lamda']<-fitting4[i,3]
      
      nstart3<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt,
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=122, to=151, by=1)
      
      output_july<-as.data.frame(ode(nstart3,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_july ) #Join dataframes for continuous time series
      
      output_eqbm<-output_july[dim(output_july)[1],]
      
      W2 = output_eqbm$Wt
      
      if(W2 < 0){
        W2 = 0
      }
      
      W = W2
      
      gamma2 = (1 - ((1-(W2/(W2+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma2 < 0){
        gamma2 = 0
      }
      
      fitting4[i,7] <- 0.5 * W2 * gamma2
      
      # Apply MDA and run through rainy season (end of october)
      time<-seq(from=152, t=273, by=1)
      
      nstart4<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply MDA 3
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_oct=as.data.frame(ode(nstart4,time,schisto_halstead_2pops_mda,params)) 
      
      output_all<-rbind(output_all, output_oct ) #Join dataframes for continuous time series
      
      output_eqbm<-output_oct[dim(output_oct)[1],]
      
      #End of rainy season, go back to low lamda and run until Feb (end of year)
      time<-seq(274, 365, 1)
      params['lamda']<-fitting4[i,2]
      
      nstart5<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt, 
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_Feb13=as.data.frame(ode(nstart5,time,schisto_halstead_2pops_mda,params)) 
      output_Feb13_eqbm<-output_Feb13[dim(output_Feb13)[1],]
      
      W3 = output_Feb13_eqbm$Wt
      
      if(W3 < 0){
        W3 = 0
      }
      
      W = W3
      
      gamma3 = (1 - ((1-(W3/(W3+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma3 < 0){
        gamma3 = 0
      }
      
      fitting4[i,8]<- 0.5 * W3 * gamma3
      
      output_all<-rbind(output_all, output_Feb13) 
      
      #Apply MDA and run from Feb '13 to beginning of rainy season (June)
      time<-seq(from=366, t=483, by=1)
      
      nstart6<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply MDA 4
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_june13=as.data.frame(ode(nstart6,time,schisto_halstead_2pops_mda,params)) 
      
      output_all<-rbind(output_all, output_june13 ) #Join dataframes for continuous time series
      
      output_eqbm<-output_june13[dim(output_june13)[1],]
      
      #Change beta value to correspond to high transmission season, run to Sept '13 follow up
      params['lamda']<-fitting4[i,3]
      
      nstart7<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt,
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=484, to=584, by=1)
      
      output_sept<-as.data.frame(ode(nstart7,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_sept ) #Join dataframes for continuous time series
      
      output_eqbm<-output_sept[dim(output_sept)[1],]
      
      W4 = output_eqbm$Wt
      
      if(W4 < 0){
        W4 = 0
      }
      
      W = W4
      
      gamma4 = (1 - ((1-(W4/(W4+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma4 < 0){
        gamma4 = 0
      }
      
      fitting4[i,9] <- 0.5 * W4 * gamma4
      
      fitting4[i,10]<-Prob_gaussian(y=fitting4[i,6], mu=W_baseline, sd=W_baseline_se)  # only baseline
      fitting4[i,11]<-Prob_gaussian(y=fitting4[i,7], mu=W_5month_field_mean, sd=W_5month_field_se) # only data point 2
      fitting4[i,12]<-Prob_gaussian(y=fitting4[i,8], mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3
      fitting4[i,13]<-Prob_gaussian(y=fitting4[i,9], mu=W_Sep13_field_mean, sd=W_Sep13_field_se) #only data point 3
      
      
      plot(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), type='l', lwd = 2,
           xlab = 'time', ylab = '~Mean Worm burden', ylim = c(-5,45))
      segments(x0 = 122, x1 = 274, y0 = -4, y1 = -4, col = 'blue', lwd=2)
      segments(x0 = 515, x1 = 605, y0 = -4, y1 = -4, col = 'blue', lwd=2)
      text(x = 200, y = -2.5, labels = 'Rainy season', col = 'blue', cex = 0.75)
      text(x = 555, y = -2.5, labels = 'Rainy season', col = 'blue', cex = 0.75)
      lines(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wu), col = 'red')
      points(x = c(0, 150, 364, 584),
             y = c(W_baseline, W_5month_field_mean, W_Feb13_field_mean,W_Sep13_field_mean),
             pch = 1, col = 'red')
      legend('topleft', lwd = 1, legend = c('Wt', 'Wu'), col = c(1,2))
      
      #likelihood for hitting three follow-up points
      fitting4[i,14]<-fitting4[i,11]*fitting4[i,12]*fitting4[i,13]
      
      #likelihood for hitting all 4 points
      fitting4[i,15]<-fitting4[i,10]*fitting4[i,11]*fitting4[i,12]*fitting4[i,13]
      
    }
    
  }
  
} # end of i loop

#Save results, make subsets ############################    
write.csv(fitting4, 'C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results_beta_2lamda_rand3.csv',
          row.names = F)

ftng<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results_beta_2lamda_rand3.csv')

ftng$lamda.twa<-(301/584)*ftng$lamda1+(253/584)*ftng$lamda2

for(i in 1:nrow(ftng)){
  ftng[i,17]<-get_Ro_beta_lamda(muPq = p.dead, beta = ftng[i,1], lamda = ftng[i,16])[3]
}
colnames(ftng)[17]='R0.twa'

good<-subset(ftng, likelihood_range234 >0 & snail.prev <10)

plot(x=good$beta1, y=good$lamda1, pch = 16, xlab = 'beta', ylab = 'lamda', 
     ylim = c(0, max(good$lamda2)))
  points(x = good$beta1, y = good$lamda2, pch = 16, col = 'red')
  points(x = good$beta1, y = good$lamda.twa, pch = 16, col = 'blue')
  legend('topright', legend = c('lamda1', 'lamda2', 'lamda.twa'), 
         pch = 16, col = c(1,2,4), cex=0.5)
 
hist(good$R0.twa) 

hist(good$beta1)
hist(good$lamda1)
hist(good$lamda2)

hist(-log(good$likelihood_range234))

  for(i in 1:nrow(good)){
    good[i,18]=(-log(good[i,14])/sum(-log(good[,14])))
  }

hist(good$V18)

beta.samp<-sample(good$beta1, size=1000000, replace = T, prob = good$V18)
  hist(beta.samp)

lamda1.samp<-sample(good$lamda1, size=1000000, replace = T, prob = good$V18)
  hist(lamda1.samp)

lamda2.samp<-sample(good$lamda2, size=1000000, replace = T, prob = good$V18)
  hist(lamda2.samp)
  
lamda.twa.samp<-sample(good$lamda.twa, size = 1000000, replace = T, prob = good$V18)
  hist(lamda.twa.samp)

At.Ch<-sapply(1:100000, function(i) get_Ro_beta_lamda(muPq = p.dead, phi_Nq = sample(At.mean, 1),
                    beta = sample(beta.samp, 1), lamda = sample(lamda.twa.samp, 1))[3])
  hist(At.Ch)

#First run with 1000 parameter combinations ################
ftng4<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results_beta_2lamda_rand.csv')  

ftng4$likelihood_range1<- -log(ftng4$likelihood_range1 +1)
ftng4$likelihood_range2<- -log(ftng4$likelihood_range2 +1)
ftng4$likelihood_range3<- -log(ftng4$likelihood_range3 +1)
ftng4$likelihood_range12<- -log(ftng4$likelihood_range12 +1)
ftng4$likelihood_range13<- -log(ftng4$likelihood_range13 +1)
ftng4$likelihood_range123<- -log(ftng4$likelihood_range123 +1)



ftng4.123<-subset(ftng4, likelihood_range123 != 0) #parameter combinations that fit all 3 points
ftng4.12<-subset(ftng4, likelihood_range12 != 0) #parameter combinations that fit first 2 points
ftng4.13<-subset(ftng4, likelihood_range13 != 0) #parameter combinations that fit first and third points

plot(x = c(0, 150, 365), y = c(W_baseline, W_5month_field_mean, W_Feb13_field_mean),
           pch = 1, col = 'red', ylim = c(-5,55), xlim = c(0,365), xlab = 'time', ylab = '~Mean Worm burden')
  segments(x0 = 0, x1 = 0, 
           y0 = (W_baseline - W_baseline_se), y1 = (W_baseline + W_baseline_se), col = 'red', lty = 2)
  segments(x0 = 150, x1 = 150, 
           y0 = (W_5month_field_mean - W_5month_field_se), 
           y1 = (W_5month_field_mean + W_5month_field_se), col = 'red', lty = 2)
  segments(x0 = 365, x1 = 365, 
           y0 = (W_Feb13_field_mean - W_Feb13_field_se), 
           y1 = (W_Feb13_field_mean + W_Feb13_field_se), col = 'red', lty = 2)
  segments(x0 = 122, x1 = 274, y0 = -4, y1 = -4, col = 'blue', lwd=2)
  text(x = 200, y = -2.5, labels = 'Rainy season', col = 'blue', cex = 0.75)
  legend('topleft', lwd = 1, legend = c('Wt', 'Wu'), col = c(1,2))
  
  for(i in 1:nrow(ftng4.13)){
    output_all<-numeric()
    params["beta"]<-ftng4.13[i,1]
    params["lamda"]<-ftng4.13[i,2]
    k<-params['k']
      
      # run to equilibrium
      
      time<-seq(from=0, to=50*365, by=1)
      nstart=c(S=8000,E=0,I=0, Wt=10, Wu=10, P=0)
      output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
      output_eqbm<-output[dim(output)[1],]
      output_eqbm$time<-0
      
      cov<-params["cov"]
        
        
        output_all<-rbind(output_all, output_eqbm ) #Join dataframes for continuous time series
        
        #Apply first MDA at t=1
        
        nstart1<-c(S=output_eqbm$S, 
                   E=output_eqbm$E, 
                   I=output_eqbm$I, 
                   Wt=output_eqbm$Wt * 0.05, #Apply first MDA
                   Wu=output_eqbm$Wu, 
                   P=0)  
        
        time = seq(from=1, to=21, by=1)
        
        output_3wks<-as.data.frame(ode(nstart1,time,schisto_halstead_2pops_mda,params))
        
        output_all<-rbind(output_all, output_3wks ) #Join dataframes for continuous time series
        
        output_eqbm<-output_3wks[dim(output_3wks)[1],]
        
        #Apply second MDA at t=21 and then run until rainy season (june)
        
        nstart2<-c(S=output_eqbm$S, 
                   E=output_eqbm$E, 
                   I=output_eqbm$I, 
                   Wt=output_eqbm$Wt * 0.05, #Apply first MDA
                   Wu=output_eqbm$Wu, 
                   P=0)  
        
        time = seq(from=22, to=121, by=1)
        
        output_june<-as.data.frame(ode(nstart2,time,schisto_halstead_2pops_mda,params))
        
        output_all<-rbind(output_all, output_june ) #Join dataframes for continuous time series
        
        output_eqbm<-output_june[dim(output_june)[1],]
        
        #Change beta value to correspond to high transmission season, run for 1 month
        params["lamda"]<ftng4.13[i,3]
        
        
        nstart3<-c(S=output_eqbm$S, 
                   E=output_eqbm$E, 
                   I=output_eqbm$I, 
                   Wt=output_eqbm$Wt, #Apply first MDA
                   Wu=output_eqbm$Wu, 
                   P=0)  
        
        time = seq(from=122, to=151, by=1)
        
        output_july<-as.data.frame(ode(nstart3,time,schisto_halstead_2pops_mda,params))
        
        output_all<-rbind(output_all, output_july ) #Join dataframes for continuous time series
        
        output_eqbm<-output_july[dim(output_july)[1],]
        
        
        # Apply MDA and run through rainy season (end of october)
        time<-seq(from=152, t=273, by=1)
        
        nstart4<-c(S=output_eqbm$S, 
                   E=output_eqbm$E,
                   I=output_eqbm$I, 
                   Wt=output_eqbm$Wt * 0.05, 
                   Wu=output_eqbm$Wu, 
                   P=0)
        
        output_oct=as.data.frame(ode(nstart4,time,schisto_halstead_2pops_mda,params)) 
        
        output_all<-rbind(output_all, output_oct ) #Join dataframes for continuous time series
        
        output_eqbm<-output_oct[dim(output_oct)[1],]
        
        #End of rainy season, go back to low beta and run until Feb (end of year)
        time<-seq(274, 365, 1)
        params["lamda"]<-ftng4.13[i,2]
        
        
        nstart5<-c(S=output_eqbm$S, 
                   E=output_eqbm$E,
                   I=output_eqbm$I, 
                   Wt=output_eqbm$Wt, 
                   Wu=output_eqbm$Wu, 
                   P=0)
        
        output_Feb13=as.data.frame(ode(nstart5,time,schisto_halstead_2pops_mda,params)) 
        output_Feb13_eqbm<-output_Feb13[dim(output_Feb13)[1],]
        
        
        output_all<-rbind(output_all, output_Feb13) ##
        
        
        lines(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), col = 'black')
        lines(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wu), col = 'red')
        
      
  }


beta1 = ftng4$beta[ftng4$likelihood_range1 == min(ftng4$likelihood_range1)]
beta2 = ftng4$beta[ftng4$likelihood_range2 == min(ftng4$likelihood_range2)]
beta3 = ftng4$beta[ftng4$likelihood_range3 == min(ftng4$likelihood_range3)]
beta12 = ftng4$beta[ftng4$likelihood_range12 == min(ftng4$likelihood_range12)]
beta13 = ftng4$beta[ftng4$likelihood_range13 == min(ftng4$likelihood_range13)]
beta123 = ftng4$beta[ftng4$likelihood_range123 == min(ftng4$likelihood_range123)]

lamda1 = ftng4$lamda1[ftng4$likelihood_range1 == min(ftng4$likelihood_range1)]
lamda2 = ftng4$lamda1[ftng4$likelihood_range2 == min(ftng4$likelihood_range2)]
lamda3 = ftng4$lamda1[ftng4$likelihood_range3 == min(ftng4$likelihood_range3)]
lamda12 = ftng4$lamda1[ftng4$likelihood_range12 == min(ftng4$likelihood_range12)]
lamda13 = ftng4$lamda1[ftng4$likelihood_range13 == min(ftng4$likelihood_range13)]
lamda123 = ftng4$lamda1[ftng4$likelihood_range123 == min(ftng4$likelihood_range123)]

lamda1.2 = ftng4$lamda1.1[ftng4$likelihood_range1 == min(ftng4$likelihood_range1)]
lamda2.2 = ftng4$lamda1.1[ftng4$likelihood_range2 == min(ftng4$likelihood_range2)]
lamda3.2 = ftng4$lamda1.1[ftng4$likelihood_range3 == min(ftng4$likelihood_range3)]
lamda12.2 = ftng4$lamda1.1[ftng4$likelihood_range12 == min(ftng4$likelihood_range12)]
lamda13.2 = ftng4$lamda1.1[ftng4$likelihood_range13 == min(ftng4$likelihood_range13)]
lamda123.2 = ftng4$lamda1.1[ftng4$likelihood_range123 == min(ftng4$likelihood_range123)]

lamda.weight<-function(l1, l2){
  l.w = ((213/365)*l1 + (152/365)*l2)/365
  return(l.w)
}

get_Ro_beta_lamda(muPq = p.dead, beta = ftng4.123[2,1], lamda = ftng4.123[2,2])

lamda.weight(l1 = ftng4.123[2,2], l2 = ftng4.123[2,3])

get_Ro_beta_lamda(muPq = p.dead, beta = ftng4.123[2,1], 
                  lamda = lamda.weight(l1 = ftng4.123[2,2], l2 = ftng4.123[2,3])
)

#Second run with 5000 parameter combinations ################
ftng4<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results_beta_2lamda_rand2.csv')  

ftng4$likelihood_range1<- -log(ftng4$likelihood_range1 +1)
ftng4$likelihood_range2<- -log(ftng4$likelihood_range2 +1)
ftng4$likelihood_range3<- -log(ftng4$likelihood_range3 +1)
ftng4$likelihood_range12<- -log(ftng4$likelihood_range12 +1)
ftng4$likelihood_range13<- -log(ftng4$likelihood_range13 +1)
ftng4$likelihood_range123<- -log(ftng4$likelihood_range123 +1)

for(i in 1:nrow(ftng4)){
  ftng4[i,14] = get_Ro_beta_lamda(muPq = p.dead, beta = ftng4[i,1], lamda = ftng4[i,2])[3]
  ftng4[i,15] = get_Ro_beta_lamda(muPq = p.dead, beta = ftng4[i,1], lamda = ftng4[i,3])[3]
  
}

ftng4.123<-subset(ftng4, likelihood_range1 != 0 & 
                    likelihood_range2 != 0 & 
                    likelihood_range3 != 0) #parameter combinations that fit all 3 points
beta_min = min(ftng4.123$beta1)
beta_max = max(ftng4.123$beta1)

lamda1_min = min(ftng4.123$lamda1)
lamda1_max = max(ftng4.123$lamda1)

lamda2_min = min(ftng4.123$lamda1.1)
lamda2_max = max(ftng4.123$lamda1.1)

#Refine and rerun ################
beta1s<-round(runif(500, max = beta_max, min = beta_min), digits = 8)

lamda1s<-round(runif(500, max = lamda1_max, min = lamda1_min), digits = 8)

lamda2s<-round(runif(500, max = 0.00125, min = 0.0005), digits = 8)

fitting5<-data.frame('beta1' = beta1s,
                     'lamda1' = lamda1s,
                     'lamda2' = lamda2s,
                     'snail.prev' = 0,
                     'W1' = 0,
                     'W2' = 0,
                     'W3' = 0,
                     'likelihood_range1' = 0,
                     'likelihood_range2' = 0,
                     'likelihood_range3' = 0,
                     'likelihood_range12' = 0,
                     'likelihood_range13' = 0,
                     'likelihood_range123' = 0)

for(i in 1:nrow(fitting5)){
  
  print(i)
  
  output_all<-numeric()
  params["beta"]<-fitting5[i,1]
  params["lamda"]<-fitting5[i,2]
  k<-params['k']
  
  R0<-get_Ro_beta_lamda(muPq = p.dead, beta = fitting5[i,1], lamda = fitting5[i,2])[3]
  
  if(R0 < 1){
    print(paste('R0=', R0, sep = ''))
  } else {
    
    # run to equilibrium
    
    time<-seq(from=0, to=50*365, by=1)
    nstart=c(S=8000,E=0,I=0, Wt=10, Wu=10, P=0)
    output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
    output_eqbm<-output[dim(output)[1],]
    output_eqbm$time<-0
    
    cov<-params["cov"]
    
    W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
    
    if (W_eqbm < .4){
      
      print(paste('W_eqbm=', W_eqbm, sep = ''))
      
    } else if (output_eqbm$E < 0){
      
      print(paste('E_eqbm=', output_eqbm$E, sep = ''))
      
    } else if (output_eqbm$S > output$S[10000]+100){
      
      print('Eqbm not reached')
      
      plot(output$time, output$S, type='l', xlab="time",ylab="System Variables", 
           ylim=c(0,max( output$S+output$E+output$I )), col='blue', lwd=2)
      lines(output$time,output$E,col='orange', lwd=2)
      lines(output$time,output$I,col='red', lwd=2)
      lines(output$time,output$S+output$E+output$I,col='black', lwd=2)
      
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
      
      gamma1 = (1 - ((1-(W1/(W1+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma1 < 0){
        gamma1 = 0
      }
      
      fitting5[i,4]<-(output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100
      
      fitting5[i,5] <- 0.5 * W1 * gamma1
      
      output_all<-rbind(output_all, output_eqbm ) #Join dataframes for continuous time series
      
      #Apply first MDA at t=1
      
      nstart1<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=1, to=21, by=1)
      
      output_3wks<-as.data.frame(ode(nstart1,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_3wks ) #Join dataframes for continuous time series
      
      output_eqbm<-output_3wks[dim(output_3wks)[1],]
      
      #Apply second MDA at t=21 and then run until rainy season (june)
      
      nstart2<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=22, to=121, by=1)
      
      output_june<-as.data.frame(ode(nstart2,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_june ) #Join dataframes for continuous time series
      
      output_eqbm<-output_june[dim(output_june)[1],]
      
      #Change beta value to correspond to high transmission season, run for 1 month
      params['lamda']<-fitting5[i,3]
      
      nstart3<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=122, to=151, by=1)
      
      output_july<-as.data.frame(ode(nstart3,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_july ) #Join dataframes for continuous time series
      
      output_eqbm<-output_july[dim(output_july)[1],]
      
      W2 = output_eqbm$Wt
      
      if(W2 < 0){
        W2 = 0
      }
      
      W = W2
      
      gamma2 = (1 - ((1-(W2/(W2+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma2 < 0){
        gamma2 = 0
      }
      
      fitting5[i,6] <- 0.5 * W2 * gamma2
      
      # Apply MDA and run through rainy season (end of october)
      time<-seq(from=152, t=273, by=1)
      
      nstart4<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, 
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_oct=as.data.frame(ode(nstart4,time,schisto_halstead_2pops_mda,params)) 
      
      output_all<-rbind(output_all, output_oct ) #Join dataframes for continuous time series
      
      output_eqbm<-output_oct[dim(output_oct)[1],]
      
      #End of rainy season, go back to low beta and run until Feb (end of year)
      time<-seq(274, 365, 1)
      params['lamda']<-fitting5[i,2]
      
      nstart5<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt, 
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_Feb13=as.data.frame(ode(nstart5,time,schisto_halstead_2pops_mda,params)) 
      output_Feb13_eqbm<-output_Feb13[dim(output_Feb13)[1],]
      
      W3 = output_Feb13_eqbm$Wt
      
      if(W3 < 0){
        W3 = 0
      }
      
      W = W3
      
      gamma3 = (1 - ((1-(W3/(W3+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma3 < 0){
        gamma3 = 0
      }
      
      fitting5[i,7]<- 0.5 * W3 * gamma3
      output_all<-rbind(output_all, output_Feb13) ##
      
      fitting5[i,8]<-Prob_gaussian(y=fitting5[i,5], mu=W_baseline, sd=W_baseline_se)  # only baseline
      fitting5[i,9]<-Prob_gaussian(y=fitting5[i,6], mu=W_5month_field_mean, sd=W_5month_field_se) # only data point 2
      fitting5[i,10]<-Prob_gaussian(y=fitting5[i,7], mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3
      
      plot(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), type='l', lwd = 2,
           xlab = 'time', ylab = '~Mean Worm burden', ylim = c(-5,45))
      segments(x0 = 122, x1 = 274, y0 = -4, y1 = -4, col = 'blue', lwd=2)
      text(x = 200, y = -2.5, labels = 'Rainy season', col = 'blue', cex = 0.75)
      lines(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wu), col = 'red')
      points(x = c(0, 150, 364),
             y = c(W_baseline, W_5month_field_mean, W_Feb13_field_mean),
             pch = 1, col = 'red')
      legend('topleft', lwd = 1, legend = c('Wt', 'Wu'), col = c(1,2))
      
      
      fitting5[i,11]<-fitting5[i,8] * fitting5[i,9] #Fit to baseline and 2
      fitting5[i,12]<-fitting5[i,8] * fitting5[i,10] #Fit to baseline and 3
      fitting5[i,13]<-fitting5[i,8] * fitting5[i,9] * fitting5[i,10]#Fit to baseline 2 and 3
      
    }
    
  }
  
} # end of i loop

#Save and analyze ##################
write.csv(fitting5, 'C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results_beta_2lamda_rand3.csv',
          row.names = F)

ftng5<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/R Codes/MLE_results_beta_2lamda_rand3.csv')  

ftng5$likelihood_range1<- -log(ftng5$likelihood_range1 +1)
ftng5$likelihood_range2<- -log(ftng5$likelihood_range2 +1)
ftng5$likelihood_range3<- -log(ftng5$likelihood_range3 +1)
ftng5$likelihood_range12<- -log(ftng5$likelihood_range12 +1)
ftng5$likelihood_range13<- -log(ftng5$likelihood_range13 +1)
ftng5$likelihood_range123<- -log(ftng5$likelihood_range123 +1)

for(i in 1:nrow(ftng5)){
  ftng5[i,14] = get_Ro_beta_lamda(muPq = p.dead, beta = ftng5[i,1], lamda = ftng5[i,2])[3]
  ftng5[i,15] = get_Ro_beta_lamda(muPq = p.dead, beta = ftng5[i,1], lamda = ftng5[i,3])[3]
  ftng5[i,16] = (((213/365)*ftng5[i,2]) + ((152/365) * ftng5[i,3])) #time-weighted average of lamda
  ftng5[i,17] = get_Ro_beta_lamda(muPq = p.dead, beta = ftng5[i,1], lamda = ftng5[i,16])[3]
  
}

colnames(ftng5)[c(14:17)]<-c('R0_low', 'R0_high', 'lamda_twa', 'R0_twa')

ftng5.123<-subset(ftng5, likelihood_range123 != 0 & 
                         likelihood_range2 != 0 & 
                         likelihood_range3 != 0) #parameter combinations that fit all 3 points

snail.prev<-mean(ftng5.123$snail.prev)
beta<-mean(ftng5.123$beta1)
lamda1<-mean(ftng5.123$lamda1)
lamda2<-mean(ftng5.123$lamda2)
lamda.twa<-mean(ftng5.123$lamda_twa)


#Function to observe fit to epi data ###################
fit.test<-function(b1, b2, l1, l2){

  output_all<-numeric()
  params["beta"]<-b1
  params["lamda"]<-l1
  k<-params['k']
  
  R0<-get_Ro_beta_lamda(muPq = p.dead, beta = b1, lamda = l1)[3]
  
  if(R0 < 1){
    print(paste('R0=', R0, sep = ''))
  } else {
    
    # run to equilibrium
    
    time<-seq(from=0, to=50*365, by=1)
    nstart=c(S=8000,E=0,I=0, Wt=10, Wu=10, P=0)
    output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
    output_eqbm<-output[dim(output)[1],]
    output_eqbm$time<-0
    
    cov<-params["cov"]
    
    W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
    
    if (W_eqbm < .4){
      
      print(paste('W_eqbm=', W_eqbm, sep = ''))
      
    } 
    
    else {
      
      plot(output$time, output$S, type='l', xlab="time",ylab="System Variables", 
           ylim=c(0,max( output$S+output$E+output$I )), col='blue', lwd=2)
      lines(output$time,output$E,col='orange', lwd=2)
      lines(output$time,output$I,col='red', lwd=2)
      lines(output$time,output$S+output$E+output$I,col='black', lwd=2)
      W1 <- output_eqbm$Wt
    
      
      output_all<-rbind(output_all, output_eqbm ) #Join dataframes for continuous time series
      
      #Apply first MDA at t=1
      
      nstart1<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=1, to=21, by=1)
      
      output_3wks<-as.data.frame(ode(nstart1,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_3wks ) #Join dataframes for continuous time series
      
      output_eqbm<-output_3wks[dim(output_3wks)[1],]
      
      #Apply second MDA at t=21 and then run until rainy season (june)
      
      nstart2<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=22, to=121, by=1)
      
      output_june<-as.data.frame(ode(nstart2,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_june ) #Join dataframes for continuous time series
      
      output_eqbm<-output_june[dim(output_june)[1],]
      
      #Change beta value to correspond to high transmission season, run for 1 month
      params['beta']<-b2
      params["lamda"]<-l2
      
      
      nstart3<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=122, to=151, by=1)
      
      output_july<-as.data.frame(ode(nstart3,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_july ) #Join dataframes for continuous time series
      
      output_eqbm<-output_july[dim(output_july)[1],]
      
      
      # Apply MDA and run through rainy season (end of october)
      time<-seq(from=152, t=273, by=1)
      
      nstart4<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, 
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_oct=as.data.frame(ode(nstart4,time,schisto_halstead_2pops_mda,params)) 
      
      output_all<-rbind(output_all, output_oct ) #Join dataframes for continuous time series
      
      output_eqbm<-output_oct[dim(output_oct)[1],]
      
      #End of rainy season, go back to low beta and run until Feb (end of year)
      time<-seq(274, 365, 1)
      params['beta']<-b1
      params["lamda"]<-l1
      
      
      nstart5<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt, 
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_Feb13=as.data.frame(ode(nstart5,time,schisto_halstead_2pops_mda,params)) 
      output_Feb13_eqbm<-output_Feb13[dim(output_Feb13)[1],]
      
      
      output_all<-rbind(output_all, output_Feb13) ##
      
      
      plot(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), type='l', lwd = 2,
           xlab = 'time', ylab = '~Mean Worm burden', ylim = c(-5,45))
      segments(x0 = 122, x1 = 274, y0 = -4, y1 = -4, col = 'blue', lwd=2)
      text(x = 200, y = -2.5, labels = 'Rainy season', col = 'blue', cex = 0.75)
      lines(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wu), col = 'red')
      points(x = c(0, 150, 364),
             y = c(W_baseline, W_5month_field_mean, W_Feb13_field_mean),
             pch = 1, col = 'red')
      legend('topleft', lwd = 1, legend = c('Wt', 'Wu'), col = c(1,2), cex = 0.75)
      legend('top', legend = c(paste('beta = ', b1, sep = ''),
                                       paste('lamd1 = ', l1, sep = ''),
                                       paste('lamda2 = ', l2, sep = '')),
             cex = 0.75)
      
    }
    
  }
  
}

fit.test(b1 = 6.44e-06, b2 = 6.44e-06, l1 = pt3$lamda, l2 = pt3$lamda)

Rlo<-get_Ro_beta_lamda(muPq = p.dead, beta = beta, lamda = lamda1)[3]
Rhi<-get_Ro_beta_lamda(muPq = p.dead, beta = beta, lamda = lamda2)[3]
R.twa1<-get_Ro_beta_lamda(muPq = p.dead, beta = beta, lamda = lamda.twa)[3]
R.twa2 = (213/365)*Rlo + (152/365)*Rhi


