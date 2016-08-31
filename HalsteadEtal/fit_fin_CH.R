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

#Generate beta and lamda ranges to test ####################
betas<-runif(5000, 1e-6, 5e-5)
  hist(betas)
  min(betas)
  max(betas)
  
lamda1s<-runif(5000, 1e-6, 5e-5)
  hist(lamda1s)
  min(lamda1s)
  max(lamda1s)
  
lamda2s<-runif(5000, 1e-4, 1e-3) 
  hist(lamda2s)
  min(lamda2s)
  max(lamda2s)
  
fit.df<-data.frame('beta' = betas,
                   'lamda1' = lamda1s,
                   'lamda2' = lamda2s,
                   'lamda.twa' = 0,
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

for(i in 1:nrow(fit.df)){
  
  print(i)
  
  output_all<-numeric()
  params["beta"]<-fit.2[i,1]
  params["lamda"]<-fit.2[i,2]
  k<-params['k']
  
  R0<-get_Ro_beta_lamda(muPq = p.dead, beta = fit.2[i,1], lamda = fit.2[i,2])[3]
  
  fit.2[i,4] = (375/594) * fit.2[i,3] + (219/594) * fit.2[i,3]
  fit.2[i,5] = R0 
  
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
    
    fit.2[i,6]<-(output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100
    
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
      
      fit.2[i,7] <- 0.5 * W1 * gamma1
      
      output_all<-rbind(output_all, output_eqbm ) #Join dataframes for continuous time series
      
      #Apply first MDA at t=1; assumed to be Feb 1, 2013
      
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
      
      #Apply second MDA at t=21 (Feb 22) and then run until rainy season (june 15; 112 days)
      
      nstart2<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply second MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=22, to=134, by=1)
      
      output_june<-as.data.frame(ode(nstart2,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_june ) #Join dataframes for continuous time series
      
      output_eqbm<-output_june[dim(output_june)[1],]
      
      #Change lamda value to correspond to high transmission season, run to July follow up
      params['lamda']<-fit.2[i,3]
      
      nstart3<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt,
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=135, to=156, by=1)
      
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
      
      fit.2[i,8] <- 0.5 * W2 * gamma2
      
      # Apply MDA and run through rainy season (mid-october; 15 weeks)
      time<-seq(from=157, t=262, by=1)
      
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
      time<-seq(263, 375, 1)
      params['lamda']<-fit.2[i,2]
      
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
      
      fit.2[i,9]<- 0.5 * W3 * gamma3
      
      output_all<-rbind(output_all, output_Feb13) 
      
      #Apply MDA and run from Feb '13 to beginning of rainy season (June; 18 weeks)
      time<-seq(from=376, t=502, by=1)
      
      nstart6<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply MDA 4
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_june13=as.data.frame(ode(nstart6,time,schisto_halstead_2pops_mda,params)) 
      
      output_all<-rbind(output_all, output_june13 ) #Join dataframes for continuous time series
      
      output_eqbm<-output_june13[dim(output_june13)[1],]
      
      #Change lamda value to correspond to high transmission season, run to Sept '13 follow up
      params['lamda']<-fit.2[i,3]
      
      nstart7<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt,
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=503, to=594, by=1)
      
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
      
      fit.2[i,10] <- 0.5 * W4 * gamma4
      
      fit.2[i,11]<-Prob_gaussian(y=fit.2[i,7], mu=W_baseline, sd=W_baseline_se)  # only baseline
      fit.2[i,12]<-Prob_gaussian(y=fit.2[i,8], mu=W_5month_field_mean, sd=W_5month_field_se) # only data point 2
      fit.2[i,13]<-Prob_gaussian(y=fit.2[i,9], mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3
      fit.2[i,14]<-Prob_gaussian(y=fit.2[i,10], mu=W_Sep13_field_mean, sd=W_Sep13_field_se) #only data point 3
      
      
      plot(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), type='l', lwd = 2,
           xlab = 'time', ylab = '~Mean Worm burden', ylim = c(-5,45))
      segments(x0 = 135, x1 = 262, y0 = -4, y1 = -4, col = 'blue', lwd=2)
      segments(x0 = 502, x1 = 594, y0 = -4, y1 = -4, col = 'blue', lwd=2)
      text(x = 200, y = -2.5, labels = 'Rainy season', col = 'blue', cex = 0.75)
      text(x = 525, y = -2.5, labels = 'Rainy season', col = 'blue', cex = 0.75)
      lines(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wu), col = 'red')
      points(x = c(0, 156, 375, 594),
             y = c(W_baseline, W_5month_field_mean, W_Feb13_field_mean,W_Sep13_field_mean),
             pch = 1, col = 'red')
      legend('topleft', lwd = 1, legend = c('Wt', 'Wu'), col = c(1,2))
      
      #likelihood for hitting three follow-up points
      fit.2[i,15]<-fit.2[i,12]*fit.2[i,13]*fit.2[i,14]
      
      #likelihood for hitting all 4 points
      fit.2[i,16]<-fit.2[i,11]*fit.2[i,12]*fit.2[i,13]*fit.2[i,14]
      
    }
    
  }
  
} # end of i loop


#Fit based on high transmission only between follow up 1 and 2 ####################
betas<-c(4.0e-6, 4.5e-6, 5.0e-6, 5.5e-6, 6.0e-6, 6.5e-6, 7.0e-6, 7.5e-6, 8.0e-6, 8.5e-6)
lamda1s<-c(8.0e-6, 9.0e-6, 1.0e-5, 2.0e-5, 3.0e-5, 4.0e-5, 5.0e-5, 6.0e-5, 7.0e-5, 8.0e-5)
lamda2s<-c(7.0e-5, 8.0e-5, 9.0e-5, 1.0e-4, 2.0e-4, 3.0e-4, 4.0e-4, 5.0e-4, 6.0e-4, 7.0e-4)

fit.2<-data.frame('beta' = rep(betas, each = 1000),
                   'lamda1' = rep(c(rep(lamda1s[1], times = 10),
                                  rep(lamda1s[2], times = 10),
                                  rep(lamda1s[3], times = 10),
                                  rep(lamda1s[4], times = 10),
                                  rep(lamda1s[5], times = 10),
                                  rep(lamda1s[6], times = 10),
                                  rep(lamda1s[7], times = 10),
                                  rep(lamda1s[8], times = 10),
                                  rep(lamda1s[9], times = 10),
                                  rep(lamda1s[10], times = 10)), times = 100),
                   'lamda2' = rep(lamda2s, times = 1000),
                   'lamda.twa' = 0,
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



for(i in 1:nrow(fit.2)){
  
  print(i)
  
  output_all<-numeric()
  params["beta"]<-fit.2[i,1]
  params["lamda"]<-fit.2[i,2]
  k<-params['k']
  
  R0<-get_Ro_beta_lamda(muPq = p.dead, beta = fit.2[i,1], lamda = fit.2[i,2])[3]
  
  fit.2[i,4] = (376/594) * fit.2[i,2] + (218/594) * fit.2[i,3]
  fit.2[i,5] = R0 
  
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
    
    fit.2[i,6]<-(output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I))*100
    
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
      
      fit.2[i,7] <- 0.5 * W1 * gamma1
      
      output_all<-rbind(output_all, output_eqbm ) #Join dataframes for continuous time series
      
      #Apply first MDA at t=1; assumed to be Feb 1, 2013
      
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
      
      #Apply second MDA at t=21 (Feb 22) and then run until first follow up, June 
      
      nstart2<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply second MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=22, to=156, by=1)
      
      output_fu1<-as.data.frame(ode(nstart2,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_fu1 ) #Join dataframes for continuous time series
      
      output_eqbm<-output_fu1[dim(output_fu1)[1],]
      
      W2 = output_eqbm$Wt
      
      if(W2 < 0){
        W2 = 0
      }
      
      W = W2
      
      gamma2 = (1 - ((1-(W2/(W2+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma2 < 0){
        gamma2 = 0
      }
      
      fit.2[i,8] <- 0.5 * W2 * gamma2
      
    # Apply MDA and run to second follow-up with higher transmission
      params['lamda']<-fit.2[i,3]
      
      time<-seq(from=157, t=375, by=1)
      
      nstart4<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply MDA 3
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_fu2=as.data.frame(ode(nstart4,time,schisto_halstead_2pops_mda,params)) 
      
      output_all<-rbind(output_all, output_fu2 ) #Join dataframes for continuous time series
      
      output_eqbm<-output_fu2[dim(output_fu2)[1],]
      
      W3 = output_eqbm$Wt
      
      if(W3 < 0){
        W3 = 0
      }
      
      W = W3
      
      gamma3 = (1 - ((1-(W3/(W3+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma3 < 0){
        gamma3 = 0
      }
      
      fit.2[i,9]<- 0.5 * W3 * gamma3
      
    #Apply MDA, return to low transmission, run to follow up 3
      time<-seq(376, 594, 1)
      params['lamda']<-fit.2[i,2]
      
      nstart5<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt *0.05, 
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_fu3=as.data.frame(ode(nstart5,time,schisto_halstead_2pops_mda,params)) 
      output_eqbm<-output_fu3[dim(output_fu3)[1],]
      
      output_all<-rbind(output_all, output_fu3) 
      
      W4 = output_eqbm$Wt
      
      if(W4 < 0){
        W4 = 0
      }
      
      W = W4
      
      gamma4 = (1 - ((1-(W4/(W4+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi)$value))
      if(gamma4 < 0){
        gamma4 = 0
      }
      
      fit.2[i,10] <- 0.5 * W4 * gamma4
      
      fit.2[i,11]<-Prob_gaussian(y=fit.2[i,7], mu=W_baseline, sd=W_baseline_se)  # only baseline
      fit.2[i,12]<-Prob_gaussian(y=fit.2[i,8], mu=W_5month_field_mean, sd=W_5month_field_se) # only data point 2
      fit.2[i,13]<-Prob_gaussian(y=fit.2[i,9], mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3
      fit.2[i,14]<-Prob_gaussian(y=fit.2[i,10], mu=W_Sep13_field_mean, sd=W_Sep13_field_se) #only data point 3
      
      
      plot(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), type='l', lwd = 2,
           xlab = 'time', ylab = '~Mean Worm burden', ylim = c(0,50))
      lines(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wu), col = 'red')
      lines(x = c(1:nrow(output_all)), 
            y = ((output_all$I / (output_all$I + output_all$E + output_all$S))*100), 
            col = 'purple', lty=2)
      points(x = c(0, 156, 375, 594),
             y = c(W_baseline, W_5month_field_mean, W_Feb13_field_mean,W_Sep13_field_mean),
             pch = 1, col = 'red')
      legend('topleft', lwd = 1, legend = c('Wt', 'Wu', 'S.Prev'), 
             col = c('black', 'red', 'purple'), cex=0.75)
      legend('topright', 
             legend = c(paste('beta=', fit.2[i,1], sep = ''),
                        paste('lam1=', fit.2[i,2], sep = ''),
                        paste('lam2=', fit.2[i,3], sep = '')), 
             cex = 0.75)
      
      #likelihood for hitting three follow-up points
      fit.2[i,15]<-fit.2[i,12]*fit.2[i,13]*fit.2[i,14]
      
      #likelihood for hitting all 4 points
      fit.2[i,16]<-fit.2[i,11]*fit.2[i,12]*fit.2[i,13]*fit.2[i,14]
      
    }
    
  }
  
} # end of i loop

write.csv(fit.2, 'fit.fin2.csv', row.names = FALSE)

fin<-read.csv('fit.fin2.csv')

fin<-subset(fin, snail.prev < 10 & R0 >1 & likelihood_range234 != 0)

  fin$likelihood_range234<- -log(fin$likelihood_range234)
  
  fin<-unique(fin)
  
  fin$prob = fin$likelihood_range234 / sum(fin$likelihood_range234)
  
  trips<-subset(fin, select = c(beta, lamda.twa, prob))

  tester<-data.frame('beta' = trips[sample(nrow(trips), 1000, replace = T, prob = trips$prob), 1],
                     'lamda'= trips[sample(nrow(trips), 1000, replace = T, prob = trips$prob), 2],
                     'R0' = 0)
  
  for(i in 1:nrow(tester)){
    tester[i,3] = get_Ro_beta_lamda(muPq = p.dead, beta = tester[i,1], lamda = tester[i,2])[3]
  }
    plot(density(tester$R0))
  
  
  R0.samp<-sapply(1:10000, 
                  function(i) get_Ro_beta_lamda(muPq = p.dead,
                                                phi_Nq = 1,
                                                beta = sample(fin$beta, size = 1, replace = T, prob = fin$prob),
                                                lamda = sample(fin$lamda.twa, size = 1, replace = T, prob = fin$prob))[3]) 
  plot(density(R0.samp))
  
  beta = fin$beta[fin$likelihood_range234 == min(fin$likelihood_range234)]

  lamda1 = fin$lamda1[fin$likelihood_range234 == min(fin$likelihood_range234)]

  lamda2 = fin$lamda2[fin$likelihood_range234 == min(fin$likelihood_range234)]
  
  
fit.vis<-function(beta, lamda1, lamda2){
  
  output_all<-numeric()
  params["beta"]<-beta
  params["lamda"]<-lamda1
  k<-params['k']
    
    # run to equilibrium
    
    time<-seq(from=0, to=50*365, by=1)
    nstart=c(S=8000,E=0,I=0, Wt=10, Wu=10, P=0)
    output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
    output_eqbm<-output[dim(output)[1],]
    output_eqbm$time<-0
    
    cov<-params["cov"]
    
      output_all<-rbind(output_all, output_eqbm ) #Join dataframes for continuous time series
      
      #Apply first MDA at t=1; assumed to be Feb 1, 2013
      
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
      
      #Apply second MDA at t=21 (Feb 22) and then run until first follow up, June 
      
      nstart2<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply second MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=22, to=156, by=1)
      
      output_fu1<-as.data.frame(ode(nstart2,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_fu1 ) #Join dataframes for continuous time series
      
      output_eqbm<-output_fu1[dim(output_fu1)[1],]
      
      
      # Apply MDA and run to second follow-up with higher transmission
      params['lamda']<-lamda2
      
      time<-seq(from=157, t=375, by=1)
      
      nstart4<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply MDA 3
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_fu2=as.data.frame(ode(nstart4,time,schisto_halstead_2pops_mda,params)) 
      
      output_all<-rbind(output_all, output_fu2 ) #Join dataframes for continuous time series
      
      output_eqbm<-output_fu2[dim(output_fu2)[1],]
      
      
      #Apply MDA, return to low transmission, run to follow up 3
      time<-seq(376, 594, 1)
      params['lamda']<-lamda1
      
      nstart5<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt *0.05, 
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_fu3=as.data.frame(ode(nstart5,time,schisto_halstead_2pops_mda,params)) 
      output_eqbm<-output_fu3[dim(output_fu3)[1],]
      
      output_all<-rbind(output_all, output_fu3) 
      
      
      plot(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wt), type='l', lwd = 2,
           xlab = 'time', ylab = '~Mean Worm burden', ylim = c(0,50))
      lines(x = c(1:nrow(output_all)), y = (0.5 * output_all$Wu), col = 'red')
      lines(x = c(1:nrow(output_all)), 
            y = ((output_all$I / (output_all$I + output_all$E + output_all$S))*100), 
            col = 'purple', lty=2)
      points(x = c(0, 156, 375, 594),
             y = c(W_baseline, W_5month_field_mean, W_Feb13_field_mean,W_Sep13_field_mean),
             pch = 1, col = 'red')
      legend('topleft', lwd = 1, legend = c('Wt', 'Wu', 'S.Prev'), 
             col = c('black', 'red', 'purple'), cex=0.75)
      legend('topright', 
             legend = c(paste('beta=', beta, sep = ''),
                        paste('lam1=', lamda1, sep = ''),
                        paste('lam2=', lamda2, sep = '')), 
             cex = 0.75)
      
} 

  fit.vis(beta = 7e-6, lamda1 = 3e-5, lamda2 = 7e-4)

  