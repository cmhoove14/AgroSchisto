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

source("lib_schistoModels.R")
source("lib_parameters.R")

# fx<-function(x, mean.worm, clump){
#   op<-(1-cos(x))/((1+((mean.worm/(clump+mean.worm))*cos(x)))^1+clump)
#   op
# }
# range<-seq(0, 2*pi, (2*pi)/360)
# quartz()
# plot(range, fx(range, 1.8, 0.08), type="l", lwd=2 )


     
#Set initial values #########
yrs=30
time=seq(0,365*yrs,1)
p=1
nstart=c(S=7000,E=3750,I=1200, Wt=20, Wu=20, P=p)
parameters<-parameters_2pops_mda_Chris1
parameters["phi_P"]<-1

#### First fix beta to get I to 10% of total sails.
yrs=30
time=seq(0,365*yrs,1)
p=1
nstart=c(S=7000,E=3750,I=1200, Wt=20, Wu=20, P=p)
parameters<-parameters_2pops_mda_Chris1
parameters["phi_P"]<-1


min_beta<-1.4e-05
max_beta<-parameters_2pops_mda_Chris1["beta"]*2
beta_range<-seq(from=min_beta, to=max_beta, by=(max_beta-min_beta)/50 ) 
W_beta<-rep(0, times=length(beta_range))
Ro_range<-rep(0, times=length(beta_range))
I_range<-rep(0, times=length(beta_range))
for(b in 1:length(beta_range) ){
  params<-parameters_2pops_mda_Chris1
  params["phi_P"]<-1
  params["beta"]<-beta_range[b]
  
  output.b=as.data.frame(lsoda(nstart,time,schisto_master_2pops_mda_seas,params)) 
  eqbm.b<-output.b[dim(output.b)[1],]
  print(eqbm.b) #Infection mostly prevented by prawn population
  
   
#   #plot time series of state variables
#   quartz()
#   par(mfrow=c(2,1), oma=c(0,0,2,0))
# #   plot(output.b$time, output.b$S, type='l', xlab="time",ylab="System Variables", 
# #        ylim=c(-500,max( output.b$S+output.b$E+output.b$I )), col=1, lty=1, lwd=2)
#   plot(output.b$time, output.b$I, type='l', xlab="time",ylab="System Variables", 
#        ylim=c(-500,max( output.b$S+output.b$E+output.b$I )), col=1, lty=1, lwd=2)
#   lines(output.b$time,output.b$E,col=2, lty=2, lwd=2)
#   lines(output.b$time,output.b$I,col=3, lty=2, lwd=2)
#   lines(output.b$time,output.b$S+output.b$E+output.b$I,col=4, lty=2, lwd=2)
#   legend(x=8000, y=10000, legend=c("S", "E", "I", "N"), col=1:4, lwd=2)
#   
#   plot(output.b$time,output.b$Wt, type='l', xlab="time",ylab="System Variables", ylim=c(0, max( max(output.b$Wu), max(output.b$P) )  ), col=5, lty=1 )
#   lines(output.b$time,output.b$Wu,col=6, lty=3)
#   lines(output.b$time,output.b$P,col=7, lty=2, lwd=2)
#   legend(x=8000, y=15, legend=c("Wt", "Wu", "P"), col=5:7, lwd=2)
  
  
  
  
  W_beta[b]<-eqbm.b$Wt * 0.5 * gamma_fn(eqbm.b$Wt, k=0.08)
  Ro_range[b]<-get_Ro_mesocosm_withPrawns_newModel(params, eqbm.b)
  I_range[b]<-(eqbm.b$I/(eqbm.b$S+eqbm.b$E+eqbm.b$I))*100

   
  print(c(b, Ro_range[b],  I_range[b]))
}#end of beta_range
bestBetaInd<-order(abs(I_range-10))[1]
bestBeta<-beta_range[bestBetaInd]

quartz()
par(mfrow=c(2,2))
plot(beta_range, I_range, type="l", col=1, lwd=2, ylab="Infected Snail Percentage", xlab="beta_range", ylim=c(0,max(I_range)), main="Infected snails vs beta" )
abline(h=10, col=1, lwd=1, lty=2)
abline(v=bestBeta, col=1, lwd=1, lty=3)
text(x=bestBeta, y=0, labels=bestBeta, pos=4, cex=1 )

plot(beta_range, W_beta, lwd=2, col=6, main="W vs beta", type="l")
abline(h=10, col=1, lwd=1, lty=2)
abline(v=bestBeta, col=1, lwd=1, lty=3)
text(x=bestBeta, y=0, labels=bestBeta, pos=4, cex=1 )

plot(beta_range, Ro_range, type="l", lwd=2, ylab="Ro",  xlab="beta_range", main="Ro vs beta")
abline(h=10, col=1, lwd=1, lty=2)
abline(v=bestBeta, col=1, lwd=1, lty=3)
text(x=bestBeta, y=0, labels=bestBeta, pos=4, cex=1 )
params<-parameters_2pops_mda_Chris1
params["phi_P"]<-1
params["beta"]<-bestBeta
Ro<-get_Ro_mesocosm_withPrawns_newModel(params, eqbm.b)
text(x=bestBeta, y=Ro, labels=round(Ro,4), pos=4, cex=1 )

#### Now look at the range of W for diff lamda values.

min_lamda<-parameters_2pops_mda_Chris1["lamda"]/2
max_lamda<-parameters_2pops_mda_Chris1["lamda"]*2
lamda_range<-seq(from=min_lamda, to=max_lamda, by=(max_lamda-min_lamda)/50) 
W_lamda<-rep(0, times=length(lamda_range))
W_lamda_1<-rep(0, times=length(lamda_range))
I_range<-rep(0, times=length(lamda_range))

Ro_range<-rep(0, times=length(lamda_range))
for(b in 1:length(lamda_range) ){
  params<-parameters_2pops_mda_Chris1
  params["beta"]<-bestBeta
  params["phi_P"]<-1
  params["lamda"]<-lamda_range[b]
  output.b=as.data.frame(ode(nstart,time,schisto_master_2pops_mda_seas,params)) 
  
  eqbm.b<-output.b[dim(output.b)[1],]
  print(eqbm.b) #Infection mostly prevented by prawn population
  W_lamda[b]<-eqbm.b$Wt
  W_lamda_1[b]<-eqbm.b$Wt*gamma_fn(eqbm.b$Wt, k=0.08)*0.5
  
  Ro_range[b]<-get_Ro_mesocosm_withPrawns_newModel(params, eqbm.b)
  I_range[b]<-(eqbm.b$I/(eqbm.b$S+eqbm.b$E+eqbm.b$I))*100
  
#   #plot time series of state variables
#   quartz()
#   par(mfrow=c(2,1), oma=c(0,0,2,0))
#   plot(output.b$time, output.b$S, type='l', xlab="time",ylab="System Variables", 
#        ylim=c(0,max( output.b$S+output.b$E+output.b$I )), col=1, lty=1, lwd=2)
#   lines(output.b$time,output.b$E,col=2, lty=2, lwd=2)
#   lines(output.b$time,output.b$I,col=3, lty=2, lwd=2)
#   lines(output.b$time,output.b$S+output.b$E+output.b$I,col=4, lty=2, lwd=2)
#   legend(x=8000, y=10000, legend=c("S", "E", "I", "N"), col=1:4, lwd=2)
#   
#   plot(output.b$time,output.b$Wt, type='l', xlab="time",ylab="System Variables", ylim=c(0, max( max(output.b$Wu), max(output.b$P) )  ), col=5, lty=1 )
#   lines(output.b$time,output.b$Wu,col=6, lty=3)
#   lines(output.b$time,output.b$P,col=7, lty=2, lwd=2)
#   legend(x=8000, y=15, legend=c("Wt", "Wu", "P"), col=5:7, lwd=2)
  
  print(c(b, Ro_range[b],  I_range[b], W_lamda[b], W_lamda_1[b]))
  

}#end of beta_range

quartz()
par(mfrow=c(2,1))
plot(lamda_range, W_lamda, type="l", lwd=2, col=6, ylab="W/I", xlab="lamda_range")
lines(lamda_range, W_lamda_1, type='l', lwd=2, col=6, lty=3)
lines(lamda_range, I_range, type='l', lwd=2, col=2)
legend(lamda_range[40], 15, legend=c("W", "I"), col=c(6,2), lwd=2)
plot(lamda_range, Ro_range, type="l", lwd=2)

# 
#   params<-parameters_2pops_mda_Chris1
#   params["phi_P"]<-1
#   params["beta"]<-1.75e-5
#   params["lamda"]<-1.5e-5
#   get_Ro_mesocosm_withPrawns_newModel(params, eqbm.b)


###################################################################



#get_Ro_mesocosm_withPrawns(parameters, eqbm.b)

p=1
nstart=c(S=7000,E=3750,I=1200, Wt=20, Wu=20, P=p)
parameters<-parameters_2pops_mda_Chris1
parameters["beta"]<-bestBeta
parameters["lamda"]<-1.47e-04
CC_P<-1
parameters["phi_P"]<-1/CC_P

#run to eqbm
time=seq(0,30*365,1)
output=as.data.frame(ode(nstart,time,schisto_master_2pops_mda_seas,parameters)) 
nstart=c(S=output$S[dim(output)[1]],E=output$E[dim(output)[1]],I=output$I[dim(output)[1]], Wt=output$Wt[dim(output)[1]], Wu=output$Wu[dim(output)[1]], P=p) 

mda_days<-c((28*7+1),(31*7+1)) #end of 28 weeks and 31 weeks
p=1
output<-Senegal_mda_halstead_1(nstart, parameters, mda_days)
quartz()
par(mfrow=c(2,1))
plot(output$time, output$Wt, type='l', lwd=2, col=1, ylim=c(0, max(output$Wu)), main="treated and untreated worms")
abline(h=0.05*output$Wt[mda_days[1]-1], lty=2)
plot(output$time, (output$I)*parameters["lamda"], lty=2, type="l", col=1, main="contribution of new worms" )

quartz()
par(mfrow=c(2,2))
plot(output$time, output$Wt, type='l', lwd=2, col=1, ylim=c(0, max(output$Wu)), main="treated and untreated worms")
lines(output$time, output$Wu, lwd=2, col=2)
legend(x=300, y=20, legend=c("Wt", "Wu"), col=c(1,2), lwd=2)
cov<-parameters["cov"]
W<-(cov*output$Wt) + ((1-cov)*output$Wu)
plot(output$time,W ,  type='l', lwd=2, col=3, ylim=c(0, max(W)), main="avg worms")

plot(output$time, output$S, type='l', lwd=2, col=3, ylim=c(0, max(output$S+output$E+output$I)), main="Snails")
lines(output$time, output$E, lwd=2, col="blue")
lines(output$time, output$I, lwd=2, col="red")
lines(output$time, output$S+output$E+output$I, lwd=2, col="black")


# eqbm<-output[dim(output)[1],]
# params<-parameters_2pops_mda_Chris
# get_Ro_mesocosm_withPrawns(params, eqbm)
# 

#### The function for prevalence #####################################


######################################################################
#### Fit the beta of the villages to their worm burden W at baseline and after PZQ #############
#### Thr baseline condition is no prawns, no agro ####################
######################################################################

# Village used is Lampsar II.
# The baseline egg burdens are divided by 3.6 (page 3, French et al 2015)
EggToWormConvert<-3.6
W_baseline<-6.5/EggToWormConvert
W_baseline_k<-0.08
W_baseline_sd<-sqrt( (W_baseline)+((W_baseline^2)/W_baseline_k ) )
N<-129
W_baseline_se<-W_baseline_sd/sqrt(N)
#W_limits_baseline<-c( W_baseline-(1.96*W_baseline_se), W_baseline+(1.96*W_baseline_se))

W_5month_field_mean<-1.5/EggToWormConvert
W_5month_field_k<-0.02
W_5month_field_sd<-sqrt( (W_5month_field_mean)+((W_5month_field_mean^2)/W_5month_field_k ) )
N<-129
W_5month_field_se<-W_5month_field_sd/sqrt(N)
timepoints<-c(0, 28*7, 51*7)

W_Feb13_field_mean<-161/EggToWormConvert
W_Feb13_field_k<-0.21
W_Feb13_field_sd<-sqrt( (W_Feb13_field_mean)+((W_Feb13_field_mean^2)/W_Feb13_field_k ) )
N<-129
W_Feb13_field_se<-W_Feb13_field_sd/sqrt(N)

W_Sep13_field_mean<-17.6/EggToWormConvert
W_Sep13_field_k<-0.29
W_Sep13_field_sd<-sqrt( (W_Sep13_field_mean)+((W_Sep13_field_mean^2)/W_Sep13_field_k ) )
N<-129
W_Sep13_field_se<-W_Sep13_field_sd/sqrt(N)
  
#W_limits_5months<-c(W_5month_field_mean-(1.96*W_5month_field_se), W_5month_field_mean+(1.96*W_5month_field_se))



### Obtain the negative binomial distribution for W #####
#The ecological model for the neg prob for worm burden
# i number of parasites per person
# k clumping parameter
# m mean burden
# refer wikipedia page - https://en.wikipedia.org/wiki/Negative_binomial_distribution
# formula refer - http://influentialpoints.com/Training/negative_binomial_distribution-principles-properties-assumptions.htm
k<-0.25 # clumping parameter of the negative binomial distribution

Prob_negbin<-function(i,k,m){
  p<-(gamma(k+i)/(gamma(i+1) * gamma(k))) * (1 + (m/k))^(-k-i) * ((m/k)^i)
  p
}


# W_test<-0
# Prob_negbin(i=W_test, k=W_baseline_k, m=W_baseline)
# dnbinom(W_test, s=W_baseline_k, mu=W_baseline)
# 
# quartz()
# par(mfrow=c(2,2), oma=c(0,0,2,0))
# plot(1:100,Prob_negbin(i=1:100, k=W_baseline_k, m=W_baseline) , type ="l", col=1, lwd=2)
# plot(1:100,Prob_negbin(i=1:100, k=W_5month_field_k, m=W_5month_field_mean) , type ="l", col=1, lwd=2)
# plot(1:100,Prob_negbin(i=1:100, k=W_Feb13_field_k, m=W_Feb13_field_mean) , type ="l", col=1, lwd=2)
# plot(1:100,Prob_negbin(i=1:100, k=W_Sep13_field_k, m=W_Sep13_field_mean) , type ="l", col=1, lwd=2)
# title("Negative binomial disbn", outer=TRUE)

Prob_gaussian<-function(y, mu, sd){
  p<-(1/(sqrt(2*pi)*sd) ) * exp( -((y-mu)^2)/(2*sd^2) )
  p
}

# range<-seq(from=0, to=70, by=70/1000)
# quartz()
# par(mfrow=c(2,2), oma=c(0,0,2,0))
# plot(range,Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se) , type ="l", col=1, lwd=2, xlab="W_range", ylab="PDF", main="Baseline")
# plot(range,Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se) , type ="l", col=1, lwd=2, xlab="W_range", ylab="PDF", main="First Follow Up")
# plot(range,Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) , type ="l", col=1, lwd=2, xlab="W_range", ylab="PDF", main="Second Follow Up")
# plot(range,Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se) , type ="l", col=1, lwd=2, xlab="W_range", ylab="PDF", main="Third Follow Up")
# title("Distribution of the mean worm burden", outer=TRUE)
# 
# prod13<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) 
# prod134<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) * as.vector(Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se)) 
# prod123<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * as.vector(Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se)) * as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) 
# prod1234<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * as.vector(Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se))  * as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) * as.vector(Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se)) 
# 
# quartz()
# par(mfrow=c(2,2), oma=c(0,0,2,0))
# plot(range, prod13, type='l', col=2, main="points 1+3")
# plot(range, prod134, type='l', col=2, main="points 1+3+4")
# plot(range, prod123, type='l', col=2, main="points 1+2+3")
# plot(range, prod1234, type='l', col=2, main="points 1+2+3+4")
# 
# quartz()
# par(mfrow=c(2,2), oma=c(0,0,2,0))
# plot(range, log(prod13), type='l', col=3, main="points 1+3")
# plot(range, log(prod134), type='l', col=3, main="points 1+3+4")
# plot(range, log(prod123), type='l', col=3, main="points 1+2+3")
# plot(range, log(prod1234), type='l', col=3, main="points 1+2+3+4")
# 
# range<-seq(from=0, to=5, by=1/1000)
# prod13<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) 
# prod134<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) * as.vector(Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se)) 
# prod123<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * as.vector(Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se)) * as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) 
# prod1234<- as.vector(Prob_gaussian(y=range, mu=W_baseline, sd=W_baseline_se)) * as.vector(Prob_gaussian(y=range, mu=W_5month_field_mean, sd=W_5month_field_se))  * as.vector(Prob_gaussian(y=range, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) ) * as.vector(Prob_gaussian(y=range, mu=W_Sep13_field_mean, sd=W_Sep13_field_se)) 
# 
# quartz()
# par(mfrow=c(2,2), oma=c(0,0,2,0))
# plot(range, prod13, type='l', col=2, main="points 1+3")
# plot(range, prod134, type='l', col=2, main="points 1+3+4")
# plot(range, prod123, type='l', col=2, main="points 1+2+3")
# plot(range, prod1234, type='l', col=2, main="points 1+2+3+4")
# 
# 
# quartz()
# par(mfrow=c(2,2), oma=c(0,0,2,0))
# plot(range, log(prod13), type='l', col=3, main="points 1+3")
# plot(range, log(prod134), type='l', col=3, main="points 1+3+4")
# plot(range, log(prod123), type='l', col=3, main="points 1+2+3")
# plot(range, log(prod1234), type='l', col=3, main="points 1+2+3+4")

##### Fit lambda to the values of W at first 2,3,4 data points

  
  #### No prawns #####################################################
  CC_P<-1 
  #### CC_N=10000 #####################################################
  phi_N_new<-(1-(1/80*0.16))/10000 
mda_days<-c((28*7+1),(31*7+1)) #end of 28 weeks and 31 weeks

  #### generate a range of beta ########################################
#   min_beta<-parameters_2pops_mda_Chris1["beta"]/2
#   max_beta<-parameters_2pops_mda_Chris1["beta"]*2
#   beta_range<-seq(from=min_beta, to=max_beta, by=(max_beta-min_beta)/50 ) 

min_beta<-bestBeta
max_beta<-bestBeta
beta_range<-seq(from=min_beta, to=max_beta, by=(max_beta-min_beta)/50 ) 
#### generate a range of lamda ########################################

  # range to fit to baseline, follow up 1 and follow up 3
 # min_lamda_1<-parameters_2pops_mda_Chris1["lamda"]/2
  min_lamda_1<-1.155e-05
  max_lamda_1<-parameters_2pops_mda_Chris1["lamda"]*2
  
  #range to fit to follow up 2
  min_lamda_2<-parameters_2pops_mda_Chris1["lamda"]*10
  max_lamda_2<-parameters_2pops_mda_Chris1["lamda"]*70
  
  lamda_range_1<-seq(from=min_lamda_1, to=max_lamda_1, by=(max_lamda_1-min_lamda_1)/50) 
  lamda_range_2<-seq(from=min_lamda_2, to=max_lamda_2, by=(max_lamda_2-min_lamda_2)/50) 
  
  lamda_range<-lamda_range_1
  
  nstart=c(S=3892,E=3750,I=1200, Wt=20, Wu=20, P=p) 
  
  W_baseline_range<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  W_postmda_July12_range<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  W_postmda_Feb13_range<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  W_postmda_Sep13_range<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  
  loglikelihood_range<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  likelihood_range1<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  likelihood_range2<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  likelihood_range3<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  likelihood_range4<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  likelihood_range13<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  likelihood_range14<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  likelihood_range134<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  likelihood_range1234<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  likelihood_range124<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))
  
  InfectedSnails_mx<-matrix(0, nrow = length(lamda_range_1) , ncol = length(lamda_range_2))

  output_all<-numeric()
  timepoints<-c(198, 362, 576, 760)
  W_timepoints<-c(W_baseline, W_5month_field_mean, W_Feb13_field_mean, W_Sep13_field_mean)
    for(j1 in 1:length(lamda_range_1)){
      for(j2 in 1:length(lamda_range_2)){
      
      output_all<-numeric()
      params<-parameters_2pops_mda_Chris1
      params["phi_P"]<-1
      params["beta"]<-bestBeta
      params["lamda"]<-lamda_range_1[j1]
      params["mda"]<-0 #no mda
      
      #### run to equilibrium
      time<-seq(from=0, to=30*365, by=1)
      nstart=c(S=3892,E=3750,I=1200, Wt=20, Wu=20, P=p) 
      output<-as.data.frame(ode(nstart,time,schisto_master_2pops_mda_seas,params))
      output_eqbm<-output[dim(output)[1],]
      cov<-params["cov"]
      W_eqbm<- mean( (cov*output$Wt[(dim(output)[1]-365): dim(output)[1] ]) + ((1-cov)*output$Wu[(dim(output)[1]-365): dim(output)[1] ] ) )
      W_baseline_range[j1,j2]<-W_eqbm
      InfectedSnails<-mean( (output$I[(dim(output)[1]-365): dim(output)[1] ]/(output$S[(dim(output)[1]-365): dim(output)[1] ]+output$E[(dim(output)[1]-365): dim(output)[1] ]+output$I[(dim(output)[1]-365): dim(output)[1] ]))*100 ) 
      InfectedSnails_mx[j1,j2]<-InfectedSnails
      output_eqbm[which(output_eqbm<0)]<-0
#       quartz()
#       plot(output$time, output$Wt, col=1, lwd=2, type="l", ylim=c(0, max(output$Wt)), main="ind4")
#       
#       quartz()
#       par(mfrow=c(2,2), oma=c(0,0,2,0))
#       plot(output$time, output$S, type="l", col=1, main="S")
#       plot(output$time, output$E, type="l", col=2, main="E")
#       plot(output$time, output$I, type="l", col=3, main="I")
#       plot(output$time, output$Wt, type="l", col=4, main="Wt")
#       title(get_Ro_mesocosm_withPrawns_newModel(params, output_eqbm), outer=TRUE)
#       output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I)*100
#       
      
      
      #apply mda
      nstart1<-c(S=output_eqbm$S,E=output_eqbm$E,I=output_eqbm$I, Wt=output_eqbm$Wt, Wu=output_eqbm$Wu, P=p)                   
      output_July12<-Senegal_mda_halstead_1(nstart1, params, mda_days) 
      output_eqbm<-output_July12[dim(output_July12)[1],]
      cov<-params["cov"]
      W_postmda<-as.numeric((cov*output_eqbm$Wt ) + ((1-cov)*output_eqbm$Wu ) )
      #W_postmda_range[i]<-W_postmda
      W_postmda_July12_range[j1,j2]<-output_eqbm$Wt
      
      output_all<-rbind(output_all, output_July12 ) ##
      
      #### apply mda at first survey in July 2012.
      
      ##### apply PZQ using mda flag and then run for 1 day
      time_last<-output_July12[dim(output_July12)[1],1]
      output_PZQ2<-output_July12[dim(output_July12)[1],]
      output_PZQ2$time<-time_last+1
      output_PZQ2$Wt<-output_eqbm$Wt *( 1- params["eff"])
      output_all<-rbind(output_all, output_PZQ2) ##
      time_last<-time_last+1
      
      #### apply no PZQ using mda flag and run for 7 months, July to Feb = (4*31)+(3*30) = 214 days.
      time<-seq(time_last, time_last+214,1)
      nstart<-c(S=output_PZQ2$S, E=output_PZQ2$E,I=output_PZQ2$I, Wt=output_PZQ2$Wt, Wu=output_PZQ2$Wu, P=p)
      params["lamda"]<-lamda_range_2[j2]
      
      output_Feb13=as.data.frame(ode(nstart,time,schisto_master_2pops_mda_seas,params)) 
      output_Feb13_eqbm<-output_Feb13[dim(output_Feb13)[1],]
      
      W_postmda_Feb13_range[j1,j2]<-output_Feb13_eqbm$Wt
      output_all<-rbind(output_all, output_Feb13[-1,]) ##
      
      time_last<-output_Feb13[dim(output_Feb13)[1],1]
      
      
      #### apply mda at next survey in Feb 2013.
      ##### apply PZQ using mda flag and then run for 1 day
      output_PZQ3<-output_Feb13_eqbm
      output_PZQ3$time<-time_last+1
      output_PZQ3$Wt<-output_Feb13_eqbm$Wt * (1-params["eff"])
      output_all<-rbind(output_all, output_PZQ3) ##
      time_last<-time_last+1
      
      #### apply no PZQ using mda flag and run for 6 months, Feb to Aug 2013 = (4*31)+(2*30) = 184 days.
      params["mda"]<-0
      time<-seq(time_last, time_last+184,1)
      params["lamda"]<-lamda_range_1[j1]
      nstart<-c(S=output_PZQ3$S, E=output_PZQ3$E,I=output_PZQ3$I, Wt=output_PZQ3$Wt, Wu=output_PZQ3$Wu, P=p)
      output_Sep13=as.data.frame(ode(nstart,time,schisto_master_2pops_mda_seas,params)) 
      output_Sep13_eqbm<-output_Sep13[dim(output_Sep13)[1],]
      
      W_postmda_Sep13_range[j1,j2]<-output_Sep13_eqbm$Wt
      output_all<-rbind(output_all, output_Sep13[-1,]) ##
      
      time_last<-output_Sep13[dim(output_Sep13)[1],1]
      
      W_modified<-output_all$Wt
      for(m in 1:length(W_modified)){
        W_modified[m]<-gamma_fn(W_baseline_range[j1,j2], k=0.08)*0.5 * output_all$Wt[m]
      }
      
          quartz()
          plot(output_all$time, W_modified, col=1, lwd=2, type="l", ylim=c(0, max(max(output_all$Wt),W_timepoints ) ), main="ind4")
           points(timepoints, W_timepoints, pch=13, col=2)
#           quartz()
#           par(mfrow=c(2,2), oma=c(0,0,2,0))
#           plot(output_all$time, output_all$S, type="l", col=1, main="S")
#           plot(output_all$time, output_all$E, type="l", col=2, main="E")
#           plot(output_all$time, output_all$I, type="l", col=3, main="I")
#           plot(output_all$time, output_all$Wt, type="l", col=4, main="Wt")
#           title(get_Ro_mesocosm_withPrawns_newModel(params, output_eqbm), outer=TRUE)
#           output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I)*100
          
      LL1<-Prob_gaussian(y=W_baseline_range[j1,j2]*gamma_fn(W_baseline_range[j1,j2], k=0.08)*0.5, mu=W_baseline, sd=W_baseline_se)  # only baseline
      LL2<-Prob_gaussian(y=W_postmda_July12_range[j1,j2]*gamma_fn(W_postmda_July12_range[j1,j2], k=0.02)*0.5, mu=W_5month_field_mean, sd=W_5month_field_se)
      LL3<-Prob_gaussian(y=W_postmda_Feb13_range[j1,j2]*gamma_fn(W_postmda_Feb13_range[j1,j2], k=0.21)*0.5, mu=W_Feb13_field_mean, sd=W_Feb13_field_se) #only data point 3
      LL4<-Prob_gaussian(y=W_postmda_Sep13_range[j1,j2]*gamma_fn(W_postmda_Sep13_range[j1,j2], k=0.29)*0.5, mu=W_Sep13_field_mean, sd=W_Sep13_field_se) #only data point 4
      LL13<-LL1 * LL3
      LL14<-LL1 * LL4
      LL134<-LL1 * LL3 * LL4
      LL1234<- LL1 * LL2* LL3 * LL4
      LL124<-LL1 * LL2 * LL4
      
      likelihood_range1[i,j]<-LL1
      likelihood_range2[i,j]<-LL2
      likelihood_range3[i,j]<-LL3
      likelihood_range4[i,j]<-LL4
      likelihood_range13[i,j]<-LL13
      likelihood_range14[i,j]<-LL14
      likelihood_range134[i,j]<-LL134
      likelihood_range1234[i,j]<-LL1234
      likelihood_range124[i,j]<-LL124
      

      #save all the matrices
      opfile<-"/Users/Arathi/Documents/2015/Monash /Data/SciencePaper_modelData/"
      outputfile<-paste(opfile, "betalamdaFit_LL1_v2.Rdata", sep="")
      save(likelihood_range1, file=outputfile)

      outputfile<-paste(opfile, "betalamdaFit_LL2_v2.Rdata", sep="")
      save(likelihood_range2, file=outputfile)
      
      outputfile<-paste(opfile, "betalamdaFit_LL3_v2.Rdata", sep="")
      save(likelihood_range3, file=outputfile)
      
      outputfile<-paste(opfile, "betalamdaFit_LL4_v2.Rdata", sep="")
      save(likelihood_range4, file=outputfile)
      
      outputfile<-paste(opfile, "betalamdaFit_LL13_v2.Rdata", sep="")
      save(likelihood_range13, file=outputfile)
      
      outputfile<-paste(opfile, "betalamdaFit_LL14_v2.Rdata", sep="")
      save(likelihood_range14, file=outputfile)
      
      outputfile<-paste(opfile, "betalamdaFit_LL134_v2.Rdata", sep="")
      save(likelihood_range134, file=outputfile)
      
      outputfile<-paste(opfile, "betalamdaFit_LL1234_v2.Rdata", sep="")
      save(likelihood_range1234, file=outputfile)
      
      outputfile<-paste(opfile, "betalamdaFit_LL124_v2.Rdata", sep="")
      save(likelihood_range124, file=outputfile)
      
      
      outputfile<-paste(opfile, "betalamdaFit_InfectedSnails_mx_v2.Rdata", sep="")
      save(InfectedSnails_mx, file=outputfile)
      
      print( c(i,j, W_baseline_range[i,j], InfectedSnails )) 
      
      }# end of j2 loop   
    }# end of j1 loop

  #load all the matrices
  opfile<-"/Users/Arathi/Documents/2015/Monash /Data/SciencePaper_modelData/"
  outputfile<-paste(opfile, "betalamdaFit_LL1_v2.Rdata", sep="")
  load(outputfile)
  
  outputfile<-paste(opfile, "betalamdaFit_LL2_v2.Rdata", sep="")
  load(outputfile)
  
  outputfile<-paste(opfile, "betalamdaFit_LL3_v2.Rdata", sep="")
  load(outputfile)
  
  outputfile<-paste(opfile, "betalamdaFit_LL4_v2.Rdata", sep="")
  load(outputfile)
  
  outputfile<-paste(opfile, "betalamdaFit_LL13_v2.Rdata", sep="")
  load(outputfile)
  
  outputfile<-paste(opfile, "betalamdaFit_LL14_v2.Rdata", sep="")
  load(outputfile)
  
  outputfile<-paste(opfile, "betalamdaFit_LL134_v2.Rdata", sep="")
  load(outputfile)
  
  outputfile<-paste(opfile, "betalamdaFit_LL124_v2.Rdata", sep="")
  load(outputfile)
  
  outputfile<-paste(opfile, "betalamdaFit_LL1234_v2.Rdata", sep="")
  load(outputfile)
  
  outputfile<-paste(opfile, "betalamdaFit_InfectedSnails_mx_v2.Rdata", sep="")
  load(outputfile)

  
  quartz()
  par(mfrow=c(2,2))
  plot( lamda_range, likelihood_range1, col=1, lwd=2, main="L1", type="l")
  plot( lamda_range, likelihood_range2, col=1, lwd=2, main="L2", type="l")
  plot( lamda_range, likelihood_range3, col=1, lwd=2, main="L3", type="l")
  plot( lamda_range, likelihood_range4, col=1, lwd=2, main="L4", type="l")
  
  quartz()
  par(mfrow=c(2,2))
  plot( lamda_range, likelihood_range124, col=1, lwd=2, main="L124", type="l")
  plot( lamda_range, likelihood_range3, col=1, lwd=2, main="L3", type="l")
  plot( lamda_range, likelihood_range14, col=1, lwd=2, main="L14", type="l")
  plot( lamda_range, likelihood_range1234, col=1, lwd=2, main="L1234", type="l")
  
  quartz()
  plot( lamda_range, InfectedSnails_mx, col=1, lwd=2, main="Infected Snails", type="l")
  
  likelihood_range1[which(likelihood_range1==0)]<-1e-320
  likelihood_range2[which(likelihood_range2==0)]<-1e-320
  likelihood_range3[which(likelihood_range3==0)]<-1e-320
  likelihood_range4[which(likelihood_range4==0)]<-1e-320
  likelihood_range124[which(likelihood_range124==0)]<-1e-320
  
  #low season
  negLL124<--log(likelihood_range124)
  negLL4<--log(likelihood_range4)
  negLL1234<--log(likelihood_range1234)
  negLL14<--log(likelihood_range14)
  
  ind_low1<-which(negLL124==min(negLL124))
  ind_low2<-which(negLL4==min(negLL4))
  ind_low3<-which(negLL1234==min(negLL1234))
  ind_low4<-which(negLL14==min(negLL14))
  
  
  bestLamda_low1<-lamda_range[ ind_low1]
  bestLamda_low2<-lamda_range[ ind_low2]
  bestLamda_low3<-lamda_range[ ind_low3]
  bestLamda_low4<-lamda_range[ ind_low4]
  
  bestLamda_low1<-lamda_range_1[j1]
  bestLamda_low2<-lamda_range_2[j2]
  
  params<-parameters_2pops_mda_Chris1
  params["beta"]<-bestBeta
  params["lamda"]<-bestLamda_low1
  params["phi_P"]<-1
  Ro_low1<-get_Ro_mesocosm_withPrawns_newModel(params, eqbm.b )
  
  params<-parameters_2pops_mda_Chris1
  params["beta"]<-bestBeta
  params["lamda"]<-bestLamda_low2
  params["phi_P"]<-1
  Ro_low2<-get_Ro_mesocosm_withPrawns_newModel(params, eqbm.b )
  
  lamda_timeweighted<-(bestLamda_low1*(560/760)) + (bestLamda_low2* (200/760))
  
  params<-parameters_2pops_mda_Chris1
  params["beta"]<-bestBeta
  params["lamda"]<-lamda_timeweighted
  params["phi_P"]<-1
  Ro_low3<-get_Ro_mesocosm_withPrawns_newModel(params, eqbm.b )
  
  
#   quartz()
#   par(mfrow=c(2,1))
#   plot( lamda_range, negLL124, col=1, lwd=2, main=paste("data points 1,2,4 , Ro = ",round(Ro_low1,4), sep=""), type="l")
#   abline(v=bestLamda_low1, col=3, lwd=2, lty=2)
#   plot( lamda_range, negLL4, col=1, lwd=2, main=paste("data point 4 , Ro = ",round(Ro_low2,4), sep=""), type="l")
#   abline(v=bestLamda_low2, col=3, lwd=2,lty=2)
#   
#   quartz()
#   plot( lamda_range, negLL1234, col=1, lwd=2, main=paste("data points 1,2,3,4 , Ro = ",round(Ro_low3,4), sep=""), type="l")
#   abline(v=bestLamda_low3, col=3, lwd=2, lty=2)
  
 
  #### high season
  negLL3<--log(likelihood_range3)
  ind_high<-which(negLL3==min(negLL3))
  bestLamda_high<-lamda_range[ind_high]
  
  params<-parameters_2pops_mda_Chris1
  params["beta"]<-bestBeta
  params["lamda"]<-bestLamda_high
  params["phi_P"]<-1
  Ro_high<-get_Ro_mesocosm_withPrawns_newModel(params, eqbm.b )
  
  
  quartz()
  plot( lamda_range, negLL3, col=1, lwd=2, main=paste("data point 3 , Ro = ",round(Ro_high,4), sep=""), type="l")
  abline(v=bestLamda_high, col=3, lwd=2, lty=2)
  
  
#   maxInd_3<-which(InfectedSnails_mx_filterOutput_3==max(InfectedSnails_mx_filterOutput_3), arr.ind=TRUE)
#   InfectedSnails_mx_filterOutput_3[which(InfectedSnails_mx_filterOutput_3==0)]<-1e-320
#   negLL_3<--log(InfectedSnails_mx_filterOutput_3)
#   beta_3<-beta_range[maxInd_3[1]]
#   lamda_3<-lamda_range[maxInd_3[2]]
#   
#   maxInd_4<-which(InfectedSnails_mx_filterOutput_4==max(InfectedSnails_mx_filterOutput_4), arr.ind=TRUE)
#   InfectedSnails_mx_filterOutput_4[which(InfectedSnails_mx_filterOutput_4==0)]<-1e-320
#   negLL<--log(InfectedSnails_mx_filterOutput_4)
#   beta_4<-beta_range[maxInd_4[1]]
#   lamda_4<-lamda_range[maxInd_4[2]]
#   
#   params<-parameters_2pops_mda_Chris
#   params["beta"]<-beta_3
#   params["lamda"]<-lamda_3
#   params["phi_P"]<-1
#   get_Ro_mesocosm_withPrawns_newModel(params, eqbm.b )
#   
#   params<-parameters_2pops_mda_Chris
#   params["beta"]<-beta_4
#   params["lamda"]<-lamda_4
#   params["phi_P"]<-1
#   get_Ro_mesocosm_withPrawns_newModel(params, eqbm.b )
#   
# #### run to equilibrium
# time<-seq(from=0, to=30*365, by=1)
# nstart=c(S=3892,E=3750,I=1200, Wt=20, Wu=20, P=p) 
# output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda_seas,params))
# output_eqbm<-output[dim(output)[1],]
# cov<-params["cov"]
# W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
# 
#           quartz()
#           plot(output$time, output$Wt, col=1, lwd=2, type="l", ylim=c(0, max(output$Wt)), main="ind4")
#           
#           quartz()
#           par(mfrow=c(2,2), oma=c(0,0,2,0))
#           plot(output$time, output$S, type="l", col=1, main="S")
#           plot(output$time, output$E, type="l", col=2, main="E")
#           plot(output$time, output$I, type="l", col=3, main="I")
#           plot(output$time, output$Wt, type="l", col=4, main="Wt")
#           title(get_Ro_mesocosm_withPrawns(params, output_eqbm), outer=TRUE)
#           output_eqbm$I/(output_eqbm$S+output_eqbm$E+output_eqbm$I)*100
# 
# 
# get_Ro_mesocosm_withPrawns_newModel(params, output_eqbm)
# 
# #
# 
# params<-parameters
# get_Ro_mesocosm_withPrawns()
#  chi_sq_CV<-5.991
#  UL95<-min(negLL1234)+chi_sq_CV
#  UL95_mx<-negLL1234
#  indToKeep<-which(UL95_mx <= UL95)
#  UL95_mx[-indToKeep]<-0
#  quartz()
#  image(t(UL95_mx), col=heat.colors(12), main="min to 95CI level are non red areas of negLL plot")
#  
#     
#   ll1234_sl<-likelihood_range1234[indInfected] #shortlist
#   ll12434_sl_max<-max(ll1234_sl)
#   min(ll1234_sl)
#   
#   likelihood_range1[,1]
#  which(likelihood_range1==max(likelihood_range1), arr.ind = TRUE  ) 
#   
#   quartz()
# image(1:dim(likelihood_range1)[2], 1:dim(likelihood_range1)[1], t(likelihood_range1), col=heat.colors(12))
#   
#     
#   quartz()
#   par(mfrow=c(2,2), oma=c(0,0,2,0))
#   plot(lamda_range, W_baseline_range, type="l", lwd=2,  xlab="lamda", ylab="W", main="Baseline")
#   plot(lamda_range, W_postmda_July12_range, type="l", lwd=2, main="First Follow Up",  xlab="lamda", ylab="W")
#   plot(lamda_range, W_postmda_Feb13_range, type="l", lwd=2, main="Second Follow Up",  xlab="lamda", ylab="W")
#   plot(lamda_range, W_postmda_Sep13_range, type="l", lwd=2, main="Third Follow Up",  xlab="lamda", ylab="W")
#   title("Variation of Worm Burden with Lamda", outer="TRUE")
# #   
#   quartz()
#   par(mfrow=c(2,2), oma=c(0,0,2,0))
#   plot(lamda_range, likelihood_range1, type="l", lwd=2, main="Baseline")
#   plot(lamda_range, likelihood_range2, type="l", lwd=2, main="First follow up")
#   plot(lamda_range, likelihood_range3, type="l", lwd=2, main="Second follow up")
#   plot(lamda_range, likelihood_range4, type="l", lwd=2, main="Third follow up")
#   title("Variation of Likelihood with Lamda", outer="TRUE")
#   
#   likelihood_range1[which(likelihood_range1==0)]<-1e-320
#   likelihood_range2[which(likelihood_range2==0)]<-1e-320
#   likelihood_range3[which(likelihood_range3==0)]<-1e-320
#   likelihood_range4[which(likelihood_range4==0)]<-1e-320
#   
#   quartz()
#   par(mfrow=c(2,2), oma=c(0,0,2,0))
#   plot(lamda_range, -log(likelihood_range1), type="l", lwd=2, main="Baseline")
#   plot(lamda_range, -log(likelihood_range2), type="l", lwd=2, main="First follow up")
#   plot(lamda_range, -log(likelihood_range3), type="l", lwd=2, main="Second follow up")
#   plot(lamda_range, -log(likelihood_range4), type="l", lwd=2, main="Third follow up")
#   title("Variation of Neg Log Likelihood with Lamda", outer="TRUE")
#   
#   negLogLL1<- -log(likelihood_range1)
#   negLogLL2<- -log(likelihood_range2)
#   negLogLL3<- -log(likelihood_range3)
#   negLogLL4<- -log(likelihood_range4)
#   
#   ind1<-which( negLogLL1==min(negLogLL1))
#   ind2<-which( negLogLL2==min(negLogLL2))
#   ind3<-which( negLogLL3==min(negLogLL3))
#   ind4<-which( negLogLL4==min(negLogLL4))
#   
#   params<-parameters_2pops_mda
#   params["lamda"]<-lamda_range[ind1]
#   Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
#   Ro
#   
#   params<-parameters_2pops_mda
#   params["lamda"]<-lamda_range[ind2]
#   Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
#   Ro
#   
#   params<-parameters_2pops_mda
#   params["lamda"]<-lamda_range[ind3]
#   Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
#   Ro
#   
#   params<-parameters_2pops_mda
#   params["lamda"]<-lamda_range[ind4]
#   Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
#   Ro
#   
#   
#   
#   quartz()
#   #par(mfrow=c(2,2), oma=c(0,0,2,0))
#   plot(lamda_range[(ind1-20):(ind1+20)], -log(likelihood_range1)[(ind1-20):(ind1+20)], type="l", lwd=2, main="Baseline")
#   abline(h=negLogLL1[ind1], col=3, lty=2)
#   abline(h=negLogLL1[ind1]+1.92, col=3, lty=1, lwd=2)
#   locator(n=2, type="p")
# 
#   quartz()
#   plot(lamda_range[(ind2-20):(ind2+20)], -log(likelihood_range2)[(ind2-20):(ind2+20)], type="l", lwd=2, main="First follow up")
#   abline(h=negLogLL2[ind2], col=3, lty=2)
#   abline(h=negLogLL2[ind2]+1.92, col=3, lty=1, lwd=2)
#   locator(n=2, type="p")
#   
#   
#   quartz()
#   plot(lamda_range[(ind3-60):(ind3+60)], -log(likelihood_range3)[(ind3-60):(ind3+60)], type="l", lwd=2, main="Second follow up")
#   abline(h=negLogLL3[ind3], col=3, lty=2)
#   abline(h=negLogLL3[ind3]+1.92, col=3, lty=1, lwd=2)
#   locator(n=2, type="p")
#   
#   quartz()
#   plot(lamda_range[(ind4-60):(ind4+60)], -log(likelihood_range4)[(ind4-60):(ind4+60)], type="l", lwd=2, main="Third follow up")
#   abline(h=negLogLL4[ind4], col=3, lty=2)
#   abline(h=negLogLL4[ind4]+1.92, col=3, lty=1, lwd=2)
#   locator(n=2, type="p")
#   
# ### combinations  
#   quartz()
#   par(mfrow=c(2,2), oma=c(0,0,2,0))
#   plot(lamda_range, likelihood_range13, type="l", lwd=2, main="Baseline+Second follow up")
#   plot(lamda_range, likelihood_range14, type="l", lwd=2, main="Baseline+Third follow up")
#   plot(lamda_range, likelihood_range134, type="l", lwd=2, main="Baseline+Second+Third follow up")
#   plot(lamda_range, likelihood_range1234, type="l", lwd=2, main="All 4 points")
#   title("Variation of Likelihood with Lamda", outer="TRUE")
#   
#   likelihood_range13[which(likelihood_range13==0)]<-1e-320
#   likelihood_range14[which(likelihood_range14==0)]<-1e-320
#   likelihood_range134[which(likelihood_range134==0)]<-1e-320
#   likelihood_range1234[which(likelihood_range1234==0)]<-1e-320
#   
#   negLogLL13<- -log(likelihood_range134)
#   negLogLL14<- -log(likelihood_range14)
#   negLogLL134<- -log(likelihood_range134)
#   negLogLL1234<- -log(likelihood_range1234)
#   
#   ind13<-which( negLogLL13==min(negLogLL13))
#   ind14<-which( negLogLL14==min(negLogLL14))
#   ind134<-which( negLogLL134==min(negLogLL134))
#   ind1234<-which( negLogLL1234==min(negLogLL1234))
#   
#   params<-parameters_2pops_mda
#   params["lamda"]<-lamda_range[ind13]
#   Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
#   Ro
#   params<-parameters_2pops_mda
#   params["lamda"]<-lamda_range[ind14]
#   Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
#   Ro
#   params<-parameters_2pops_mda
#   params["lamda"]<-lamda_range[ind134]
#   Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
#   Ro
#   params<-parameters_2pops_mda
#   params["lamda"]<-lamda_range[ind1234]
#   Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
#   Ro
#   
#   
#   
#   
#   quartz()
#   plot(lamda_range[(ind13-10):(ind13+10)], negLogLL13[(ind13-10):(ind13+10)], type="l", lwd=2, main="Baseline+Second follow up")
#   abline(h=min(negLogLL13), col=3, lty=2)
#   abline(h=min(negLogLL13)+1.96, col=3, lty=3)
#   locator(n=2, type="p")
#   
#   quartz()
#   plot(lamda_range[(ind14-20):(ind14+20)],negLogLL14[(ind14-20):(ind14+20)], type="l", lwd=2, main="Baseline+Third follow up")
#   abline(h=min(negLogLL14), col=3, lty=2)
#   abline(h=min(negLogLL14)+1.96, col=3, lty=3)
#   locator(n=2, type="p")
#   
#   quartz()
#   plot(lamda_range[(ind134-20):(ind134+20)], negLogLL134[(ind134-20):(ind134+20)], type="l", lwd=2, main="Baseline+Second+Third follow up")
#   abline(h=min(negLogLL134), col=3, lty=2)
#   abline(h=min(negLogLL134)+1.96, col=3, lty=3)
#   locator(n=2, type="p")
#   
#   quartz()
#   plot(lamda_range[(ind1234-20):(ind1234+20)], negLogLL1234[(ind1234-20):(ind1234+20)], type="l", lwd=2, main="All 4 points")
#   abline(h=min(negLogLL1234), col=3, lty=2)
#   abline(h=min(negLogLL1234)+1.96, col=3, lty=3)
#   locator(n=2, type="p")
#   
#   params<-parameters_2pops_mda
#   params["lamda"]<-lamda_range[ind3]
#   Ro<-get_Ro_mesocosm_withPrawns(params, eqbm)
#   Ro
#  
#  
#   
#   
# for(b in 1:3){
#   ### lets get the ouput curves for the 3 beta values
#   params<-parameters_2pops_mda
#   params["phi_P"]<-1/CC_P
#   params["beta"]<-W_bestBeta[b]
#   params["phi_N"]<-phi_N_new
#   params["mda"]<-0 #no mda
#   #### run to equilibrium
#   nstart=c(S=3892,E=3750,I=1200, Wt=20, Wu=20, P=p) 
#   time<-seq(from=0, to=30*365, by=1)
#   output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
#   output_eqbm<-output[dim(output)[1],]
#   cov<-params["cov"]
#   W_eqbm<- as.numeric( (cov*output_eqbm$Wt) + ((1-cov)*output_eqbm$Wu) )
#   
#   #apply mda
#   nstart1<-c(S=output_eqbm$S,E=output_eqbm$E,I=output_eqbm$I, Wt=output_eqbm$Wt, Wu=output_eqbm$Wu, P=p)                   
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
# params<-parameters_2pops_mda
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
