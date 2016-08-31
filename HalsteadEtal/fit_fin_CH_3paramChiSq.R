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

fx<-function(x, mean.worm, clump){
  (1-cos(x))/((1+((mean.worm/(clump+mean.worm))*cos(x)))^1+clump)
}
W<-1
k<-0.08
angles<-seq(from=0, to=2*pi, by=(2*pi)/360)

plot(angles, fx(angles,W, k ), type="b")

integrate(fx, 0, 2*pi,W, k )$value

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
p.dead = parameters["f_P"] - parameters["mu_P"]

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



#Parameter ranges ####################
betas<-c(4.0e-6, 4.5e-6, 5.0e-6, 5.5e-6, 6.0e-6, 6.5e-6, 7.0e-6, 7.5e-6, 8.0e-6, 8.5e-6)
lamda1s<-c(8.0e-6, 9.0e-6, 1.0e-5, 2.0e-5, 3.0e-5, 4.0e-5, 5.0e-5, 6.0e-5, 7.0e-5, 8.0e-5)
lamda2s<-c(7.0e-5, 8.0e-5, 9.0e-5, 1.0e-4, 2.0e-4, 3.0e-4, 4.0e-4, 5.0e-4, 6.0e-4, 7.0e-4)
#Function to simulate the transmission, given a triplet.(beta, lamda1, lamda2) #############

Tx_simulate<-function(triplet){
  
  timepoints<-numeric()
  output_all<-numeric()
  params<-parameters_2pops_mda
  params["beta"]<-triplet[1]
  params["lamda"]<-( (375/594) * triplet[2]) +( (219/594) * triplet[3]) 
  k<-params['k']
  
  # run to equilibrium
    time<-seq(from=0, to=50*365, by=1)
    nstart=c(S=3892,E=3750,I=1200, Wt=20, Wu=20, P=0)
    output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
    
    output_eqbm<-output[dim(output)[1],]
    output_eqbm$time<-0
    output_eqbm[which(output_eqbm<0)]<-0
    cov<-params["cov"]
    W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
    
      output_all<-rbind(output_all, output_eqbm ) #Join dataframes for continuous time series
      timepoints<-c(timepoints, output_all$time[length(output_all$time)])
      #Apply first MDA at t=1; assumed to be Feb 1, 2013
      
      nstart1<-c(S=output_eqbm$S, 
                 E=output_eqbm$E, 
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply first MDA
                 Wu=output_eqbm$Wu, 
                 P=0)  
      
      time = seq(from=1, to=21, by=1)
      
      params["lamda"]<-triplet[2]
      
      output_3wks<-as.data.frame(ode(nstart1,time,schisto_halstead_2pops_mda,params))
      
      output_all<-rbind(output_all, output_3wks ) #Join dataframes for continuous time series
      
      output_eqbm<-output_3wks[dim(output_3wks)[1],]
      output_eqbm[which(output_eqbm<0)]<-0
      
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
      timepoints<-c(timepoints, output_all$time[length(output_all$time)])
      
      
      output_eqbm<-output_fu1[dim(output_fu1)[1],]
      output_eqbm[which(output_eqbm<0)]<-0
      
      
      W2<-output_eqbm$Wt
    
      W<-W2
      gamma2 = (1 - ((1-(W2/(W2+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi, W, k)$value))
      if(gamma2 < 0){
        gamma2 = 0
      }
      
      W2_model<- 0.5 * W2 * gamma2
      
      # Apply MDA and run to second follow-up with higher transmission
      params['lamda']<-triplet[3]
      
      time<-seq(from=157, t=375, by=1)
      
      nstart4<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt * 0.05, #Apply MDA 3
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_fu2=as.data.frame(ode(nstart4,time,schisto_halstead_2pops_mda,params)) 
      
      output_all<-rbind(output_all, output_fu2 ) #Join dataframes for continuous time series
      timepoints<-c(timepoints, output_all$time[length(output_all$time)])
      
      output_eqbm<-output_fu2[dim(output_fu2)[1],]
      output_eqbm[which(output_eqbm<0)]<-0
      
      W3<- output_eqbm$Wt
      
      if(W3 < 0){
        W3 = 0
      }
      
      W<- W3
      
      gamma3 = (1 - ((1-(W3/(W3+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi, W, k )$value))
      if(gamma3 < 0){
        gamma3 = 0
      }
      
      W3_model<- 0.5 * W3 * gamma3
      
      #Apply MDA, return to low transmission, run to follow up 3
      time<-seq(376, 594, 1)
      params['lamda']<-triplet[2]
      
      nstart5<-c(S=output_eqbm$S, 
                 E=output_eqbm$E,
                 I=output_eqbm$I, 
                 Wt=output_eqbm$Wt *0.05, 
                 Wu=output_eqbm$Wu, 
                 P=0)
      
      output_fu3=as.data.frame(ode(nstart5,time,schisto_halstead_2pops_mda,params)) 
      output_eqbm<-output_fu3[dim(output_fu3)[1],]
      output_eqbm[which(output_eqbm<0)]<-0
      
      output_all<-rbind(output_all, output_fu3) 
      timepoints<-c(timepoints, output_all$time[length(output_all$time)])
      
      W4<- output_eqbm$Wt
      
      if(W4 < 0){
        W4 = 0
      }
      
      W<- W4
      
      gamma4 = (1 - ((1-(W4/(W4+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi, W, k)$value))
      if(gamma4 < 0){
        gamma4 = 0
      }
      
      W4_model <- 0.5 * W4 * gamma4
      
      L2<-Prob_gaussian(y=W2, mu=W_5month_field_mean, sd=W_5month_field_se) # only data point 2
      L3<-Prob_gaussian(y=W3, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3
      L4<-Prob_gaussian(y=W4, mu=W_Sep13_field_mean, sd=W_Sep13_field_se) #only data point 3
      
#       Wt_old<-output_all$Wt
#       Wt_modified<-rep(0, times=length(Wt_old) )
#       for(i in 1:length(Wt_old)){
#         W<-Wt_old[i]
#         gamma = (1 - ((1-(W/(W+k)))^(1+k)/(2*pi))*(integrate(fx, 0, 2*pi, W, k)$value))
#         if(gamma < 0){
#           gamma = 0
#         }
#         
#         Wt_modified[i]<-Wt_old[i]*0.5*gamma
#       }
#       
#       
#         quartz()
#         par(mfrow=c(2,2))
#         plot(output_all$time, Wt_old, type="l", col=1, lwd=2, main="worm")
#         plot(output_all$time, Wt_modified, type="l", col=4, lty=3, main="female worms")
#         points(timepoints, W_timepoints, pch=16, col=2, cex=1)
#         
#         plot(output_all$time, output_all$I, type="l", col=4, lty=3, main="infected")
#         plot(output_all$time, output_all$I/(output_all$S+output_all$E+output_all$I)*100, type="l", col=4, lty=3, main="proportion infected")
#         
# 
      
      
     LL234<-L2*L3*L4
     if(LL234==0)
       LL234<- 1e-320
     
     negLogL234<--log(LL234)
     
     negLogL234
    }
  
  
triplet<-as.numeric(c(mean(betas), mean(lamda1s), mean(lamda2s)) )
params<-parameters_2pops_mda
Tx_simulate(triplet)

#Fit based on high transmission only between follow up 1 and 2 ####################
op_min<-optim(par=c(min(betas), min(lamda1s), min(lamda2s)), Tx_simulate, method="Nelder-Mead")
op_min_NelderMead<-op_min$par
params<-parameters_2pops_mda
beta<-op_min_NelderMead[1]
lamda<-( (375/594) * op_min_NelderMead[2]) +( (219/594) * op_min_NelderMead[3]) 
get_Ro_beta_lamda(muPq = p.dead, phi_Nq = 1, beta, lamda, f_Nq = 1)[3]

bestBeta<-op_min_NelderMead[1]
bestLamda1<-op_min_NelderMead[2]
bestLamda2<-op_min_NelderMead[3]
triplet<-c(bestBeta,bestLamda1,bestLamda2)

#Now do the grid search on a 3D grid of beta, lamda1 and lamda2 around the best fit parameter values
beta_range<-seq(from=bestBeta*0.1, to=bestBeta*1.9, by=(bestBeta*1.5-bestBeta*0.5)/50)
lamda_range_1<-seq(from=bestLamda1*0.1, to=bestLamda1*1.9, by=(bestLamda1*1.5-bestLamda1*0.5)/50)
lamda_range_2<-seq(from=bestLamda2*0.1, to=bestLamda2*1.9, by=(bestLamda2*1.5-bestLamda2*0.5)/50)

NegLLOutput_mx.test<-load("~/RemaisWork/Schisto/R Codes/ag_schist/output.Rdata")

NegLLOutput_mx<-array(0, dim=c(length(beta_range),length(lamda_range_1),length(lamda_range_2)) )
output_all<-numeric()
timepoints<-c(0,  362, 576, 760)
W_timepoints<-c(W_baseline, W_5month_field_mean, W_Feb13_field_mean, W_Sep13_field_mean)

for(i in 1:length(beta_range)){
  for(j1 in 1:length(lamda_range_1)){
    for(j2 in 1:length(lamda_range_2)){
      triplet<-c(c(beta_range[i], lamda_range_1[j1], lamda_range_2[j2]))
      #only allow triplets where Ro>1
      params<-parameters_2pops_mda
      params["phi_P"]<-1 #remember to set this, for no prawns.
      beta<-triplet[1]
      lamda<-( (375/594) * triplet[2]) +( (219/594) * triplet[3]) 
      
      Ro<-get_Ro_beta_lamda(muPq = p.dead, phi_Nq = 1, beta, lamda, f_Nq = 1)[3]
      if(Ro<1) NegLLOutput_mx[i,j1,j2]<--log(1e-320)
      if(Ro>=1){
        NegLLOutput_mx[i,j1,j2]<-Tx_simulate(triplet)
      }
      print(c(i, j1, j2, NegLLOutput_mx[i,j1,j2] ))
    }}}

NegLLOutput_mx[which(NegLLOutput_mx==0)]<--log(1e-320)
boundary<-Tx_simulate(as.numeric(op_min_NelderMead))+7.815 #chi-square with 3 degree of freedom, 95% CI
indToKeep<-which(NegLLOutput_mx2<=boundary, arr.ind=TRUE)

shortlist<-c(as.numeric(op_min_NelderMead), exp( - Tx_simulate(as.numeric(op_min_NelderMead)) ) )
for(i in 1:dim(indToKeep)[1]){
  ind<-indToKeep[i,]
  triplet<-c(beta_range[ind[1]], lamda_range_1[ind[2]], lamda_range_2[ind[3]])
  negLL<-Tx_simulate(triplet)
  tmp_row<-c(triplet,exp(-negLL))
  shortlist<-rbind(shortlist, tmp_row)
  print(i)
}
shortlist_df<-data.frame(beta=shortlist[,1], 
                         lamda1=shortlist[,2], 
                         lamda2=shortlist[,3], 
                         p=shortlist[,4])

for(i in 1:nrow(shortlist_df)){
  shortlist_df[i,5] = -log(shortlist_df[i,4])
  shortlist_df[i,6] = (375/594) * shortlist_df[i,2] +( (219/594) * shortlist_df[i,3]) 
}

for(i in 1:nrow(shortlist_df)){
  shortlist_df[i,7] = get_Ro_beta_lamda(muPq = p.dead, 
                                        beta = shortlist_df[i,1], 
                                        lamda = shortlist_df[i,6])[3]
  shortlist_df[i,8] = shortlist_df[i,4] / sum(shortlist_df[,4])
}

colnames(shortlist_df)[c(5:8)]<-c('negLL','lamda.twa', 'R0', 'prob')
plot(density(shortlist_df$R0))

write.csv(shortlist_df, 'trips10x10x10.csv', row.names = FALSE)


params<-parameters_2pops_mda
#params["phi_P"]<-1 #remember to set this, for no prawns.
#params["mu_P"]<-1
beta<-op_min_NelderMead[1]
lamda<-( (375/594) * op_min_NelderMead[2]) +( (219/594) * op_min_NelderMead[3]) 
get_Ro_beta_lamda(muPq = 0, phi_Nq = 1, beta, lamda, f_Nq = 1, params)

Ro_range<-c(as.numeric(op_min_NelderMead), as.numeric( get_Ro_beta_lamda(muPq = 0, phi_Nq = 1, beta, lamda, f_Nq = 1, params)[3]) )
for(i in 1:dim(indToKeep)[1]){
  ind<-indToKeep[i,]
  beta<-beta_range[ind[1]]
  lamda<-( (375/594) * lamda_range_1[ind[2]]) +( (219/594) * lamda_range_2[ind[3]]) 
  ro<-as.numeric( get_Ro_beta_lamda(muPq = 0, phi_Nq = 1, beta, lamda, f_Nq = 1, params)[3] )
  
  triplet<-c(beta_range[ind[1]], lamda_range_1[ind[2]], lamda_range_2[ind[3]])
  tmp_row<-c(triplet,ro)
  Ro_range<-rbind(Ro_range, tmp_row)
  print(i)
}
Ro_df<-data.frame(beta=Ro_range[,1], lamda1=Ro_range[,2], lamda2=Ro_range[,3], ro=Ro_range[,4])
opfile<-"/Users/Arathi/Documents/2015/Monash /Data/SciencePaper_modelData/"
outputfile<-paste(opfile, "shortlistTriplets_Ro.Rdata", sep="")
save(Ro_df, file=outputfile)

range( Ro_range[,4])
mean(Ro_range[,4])
Ro_list<-as.numeric(Ro_range[,4])
Ro_list<-sort(Ro_list)
quartz()
hist(Ro_list)
points(x= Ro_range[,4][1], y=0, pch=16, col=2, cex=3)

