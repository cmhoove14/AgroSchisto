#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License <http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir and Justin Remais. This work was supported in part by the National Institutes of Health/National Science Foundation Ecology of Infectious Disease program funded by the Fogarty International Center (grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).Per the terms of this license, if you are making derivative use of this work, you must identify that your work is a derivative work, give credit to the original work, provide a link to the license, and indicate changes that were made.###############

#Read in function
require(deSolve)
require(tidyverse)
require(rootSolve)

source("~/ElimFeas_StochMod/PostReview/Refs/schisto_mods_pdd_nopdd_savio.R")

#Some useful functions #######
phi_Wk <- function(W, k) {
  if(W <= 0){
    val = 1
  } 
  #print(c(W,k))
  else{
    func <- function(x) {
    a <- ( W / (W + k) )
    b <- ((1-a)^(1+k))/(2*pi)
    return(( b*( 1-cos(x) ) / (( 1 + a*cos(x) )^(1+k)) ))
    }
    
  val = integrate(func, 0, 2*pi, subdivisions = 10000,
                  rel.tol = 1e-10, stop.on.error = FALSE)$value

  }  
    return(1-val)
}

st.er <- function(x) {
  sd(x)/sqrt(length(x))
} #Function to calculate standard error of the mean

Prevalence <- function(W, k) {
  p=1 - (1/(1+W/k)^(k))*(1+W/(1+W/k)) # fraction of humans with at least 2 parasites
  return(p)
}

#Get distributions for Epi datapoints given estimated W and k 
Prob_negbin<-function(i,k,m){
  p<-(gamma(k+i)/(gamma(i+1) * gamma(k))) * (1 + (m/k))^(-k-i) * ((m/k)^i)
  p
}
# i number of parasites per person
# k clumping parameter
# m mean burden
# refer wikipedia page - https://en.wikipedia.org/wiki/Negative_binomial_distribution
# formula refer - http://influentialpoints.com/Training/negative_binomial_distribution-principles-properties-assumptions.htm

Prob_gaussian<-function(y, mu, sd){
  p<-(1/(sqrt(2*pi)*sd) ) * exp( -((y-mu)^2)/(2*sd^2) )
  p
}

#Model parameters ####################
params = pars_Chris1

time=seq(0, 365*200, 30) #50 years to ensure equilibrium is reached

nstart=c(S=4000,E=2000,I=500, Wt=10, Wu=10)

#Epi data ###########
N <- 129 #number of people in community
EggToWormConvert <- 3.6 #convert egg burden to worm burden based on estimate of eggs/mated female from Cheever study
covrg <- 0.43 #43% coverage estimated in this community

params["cov"] <- covrg

#Baseline estimates
  W_baseline<-6.5/EggToWormConvert
  W_baseline_k<-0.08
  W_baseline_sd<-sqrt( (W_baseline)+((W_baseline^2)/W_baseline_k ) )
  W_baseline_se<-W_baseline_sd/sqrt(N)

#5 month epi data point estimates 
  W_5month_field_mean<-1.5/EggToWormConvert
  W_5month_field_k<-0.02
  W_5month_field_sd<-sqrt( (W_5month_field_mean)+((W_5month_field_mean^2)/W_5month_field_k ) )
  W_5month_field_se<-W_5month_field_sd/sqrt(N)

#Feb 13 epi data estimates (in the middle of high transmission season) 
  W_Feb13_field_mean<-161/EggToWormConvert
  W_Feb13_field_k<-0.21
  W_Feb13_field_sd<-sqrt( (W_Feb13_field_mean)+((W_Feb13_field_mean^2)/W_Feb13_field_k ) )
  W_Feb13_field_se<-W_Feb13_field_sd/sqrt(N)

#Sept 13 epi data estimates (in the middle of high transmission season)
  W_Sep13_field_mean<-17.6/EggToWormConvert
  W_Sep13_field_k<-0.29
  W_Sep13_field_sd<-sqrt( (W_Sep13_field_mean)+((W_Sep13_field_mean^2)/W_Sep13_field_k ) )
  W_Sep13_field_se<-W_Sep13_field_sd/sqrt(N)

#Timepoints of epi datapoint collection 
  timepoints<-c(1, 156, 375, 594)
  traj_time <- c(1:(max(timepoints) +2))    #time vector to run over 
  
#Estimates of measured worm burden at epi timepoints 
  measured_ws <- c(W_baseline, W_5month_field_mean, W_Feb13_field_mean, W_Sep13_field_mean)
  w_errors <- c(W_baseline_se, W_5month_field_se, W_Feb13_field_se, W_Sep13_field_se)
  measured_ks <- c(W_baseline_k, W_5month_field_k, W_Feb13_field_k, W_Sep13_field_k)

#timepoints at which transmission shifts from low to high
  seasons <- c(1, 134, 263, 502)

#Events data frame containing timepoints at which MDA occurs
  mdas <- data.frame(var = rep('Wt', 4),
                     time = c(2, 22, 157, 376),
                     value = rep(0.06, 4),     #94% efficacy
                     method = rep('mult', 4))

# Differential equation model with dynamic k and lambda ############
schisto_mod_pdd_add_ndds_dyna=function(t, n, parameters) { 
  with(as.list(parameters),{
  #set k based on time in simulation  
  k <- ifelse(t < timepoints[1], W_baseline_k, 
           ifelse(timepoints[1] <= t & t < timepoints[2], W_5month_field_k,
                  ifelse(timepoints[2] <= t & t < timepoints[3], W_Feb13_field_k,
                         ifelse(timepoints[3] <= t, W_Sep13_field_k, NA))))
    
  #set lambda (snail to human transmission) based on time in simulation  
  lamda <- ifelse(t < seasons[1], parameters[["lamda2"]], 
               ifelse(seasons[1] <= t & t < seasons[2], parameters[["lamda1"]],
                      ifelse(seasons[2] <= t & t < seasons[3], parameters[["lamda2"]],
                             ifelse(seasons[3] <= t, parameters[["lamda1"]], NA))))
    
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    N=S+E+I
    
    W=(cov*Wt) + ((1-cov)*Wu) #weighting treated and untreated populations
    #print(c(t, W, k, lamda))
    #Miracidia production; function of adult female worms alive in the system (W*H*0.5) assuming 1:1 sex ratio,
  #DD functions
    
    f = f_Wgk(W, gam, k)  
    R = R_Wv(W, v)
    phi = phi_Wk(W = W, k = k)  #Mating probability
    
    #print(c("Crowding" = f, "Immunity" = R, "Mating_prob" = phi))

    M=((0.5*W*H)*phi*f)#*m*u_H*(v*vq)
    

    dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S - beta*M*S #Susceptible snails
    
    dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails
    
    dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails
    
    #worm burden in human
    dWtdt= (lamda*I*R) - ((mu_W+mu_H)*Wt)
    dWudt= (lamda*I*R) - ((mu_W+mu_H)*Wu)
    
    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt)))
  }) 
} 

## Likelihood function for a single data point:
  pointLogLike <- function(obs_mu = measured_ws, #Measured worm burdens
                           obs_sd = w_errors,       #Measured worm burden uncertainty
                           modOut) {                #Model output W
    log(Prob_gaussian(y=modOut, mu=obs_mu, sd=obs_sd))
  }
  
#Function to simulate the transmission, given a triplet.(beta, lamda1, lamda2) #############
Tx_simulate<-function(lamdas){
  print(lamdas)
  params["lamda1"]<-lamdas[1]
  params["lamda2"]<-lamdas[2]
  params["lamda"]<-( (375/594) * lamdas[1]) +( (219/594) * lamdas[2]) 

    run_to_eqbm <- ode(nstart, time, schisto_mod_pdd_add_ndds, params)[max(time)/30,]
    
  	traj <- data.frame(ode(run_to_eqbm[c(2:6)], traj_time, schisto_mod_pdd_add_ndds_dyna,  
                           params, events = list(data = mdas)))
  	
  	observed <- traj$Wt[timepoints + 2] * 0.5 * mapply(phi_Wk, traj$Wt[timepoints + 2], measured_ks)
  	
  	point_LL <- mapply(pointLogLike, measured_ws, w_errors, observed)
  	
  	point_LL[!is.finite(point_LL)] <- NA
  	
  	logLike <- sum(point_LL, na.rm = TRUE)
  	
  	print(-logLike)
  	return(-logLike)

}

#Parameter ranges and optimization ####################
lamda1s<-seq(1e-6, 1e-3, length.out = 10)
lamda2s<-seq(5e-6, 5e-3, length.out = 10)

#Fit using Optim
op_min <- optim(par=c(lamda1s[2], lamda2s[2]), 
                Tx_simulate, method="Nelder-Mead")

#Extract best fitting parameters
  bestLamda1<-op_min$par[1]
  bestLamda2<-op_min$par[2]

  params["lamda1"]<-bestLamda1
  params["lamda2"]<-bestLamda2
  params["lamda"]<-( (375/594) * bestLamda1) +( (219/594) * bestLamda2) 
  
#Now do a grid search on a 3D grid of beta, lamda1 and lamda2 around the best fit parameter values ########
  lamda1_range<-seq(from=bestLamda1*0.1, to=bestLamda1*1.9, length.out = 50)
  lamda2_range<-seq(from=bestLamda2*0.1, to=bestLamda2*1.9, length.out = 50)
  
par_grid <- expand.grid(lamda1_range, lamda2_range)  

#Convert to list for use in map
par_grid_map <- as.list(data.frame(t(par_grid)))

negLL <- map(par_grid_map, Tx_simulate)

#combine with parameters
fin_pars <- cbind(par_grid, negLL = unlist(negLL))
fin_pars$lamda_twa <- (375/594) * fin_pars$Var1 +(219/594) * fin_pars$Var2 

#Value of negLL that serves as cutoff for 95%CI
boundary <- Tx_simulate(c(bestLamda1,bestLamda2)) + 7.815 #chi-square with 3 degree of freedom, 95% CI

#Keep obs within 95%CI
fin_pars95ci <- fin_pars %>% 
  filter(negLL <= boundary) %>% 
  arrange(lamda_twa)

save(fin_pars95ci, file = "~/ElimFeas_StochMod/PostReview/Outputs/model_fit_profile_likelihood_parameters.Rdata")