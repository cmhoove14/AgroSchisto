#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License <http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir and Justin Remais. This work was supported in part by the National Institutes of Health/National Science Foundation Ecology of Infectious Disease program funded by the Fogarty International Center (grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).Per the terms of this license, if you are making derivative use of this work, you must identify that your work is a derivative work, give credit to the original work, provide a link to the license, and indicate changes that were made.###############

#Read in function
require(deSolve)
require(tidyverse)
require(rootSolve)
require(coda)

source("~/ElimFeas_StochMod/PostReview/Refs/schisto_mods_pdd_nopdd_savio.R")

#Some useful functions #######
phi_Wk <- function(W, k) {
#  print(c(W,k))
  func <- function(x) {
    a <- ( W / (W + k) )
    b <- ((1-a)^(1+k))/(2*pi)
    return(( b*( 1-cos(x) ) / (( 1 + a*cos(x) )^(1+k)) ))
  }
  val = integrate(func, 0, 2*pi, subdivisions = 10000,
                  rel.tol = 1e-10, stop.on.error = FALSE)$value
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

epi_plot <- as.data.frame(cbind("time" = timepoints,  "Worm_burden" = measured_ws, "se" = w_errors))

epi_plot %>%
  ggplot(aes(x = timepoints, y = measured_ws)) + geom_point() + theme_bw() + ylim(0,50) + 
    geom_errorbar(ymin = measured_ws - w_errors, ymax = measured_ws + w_errors, width = 5)

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
  lamda <- ifelse(t < seasons[1], parameters[["lamda1"]], 
               ifelse(seasons[1] <= t & t < seasons[2], parameters[["lamda2"]],
                      ifelse(seasons[2] <= t & t < seasons[3], parameters[["lamda1"]],
                             ifelse(seasons[3] <= t, parameters[["lamda2"]], NA))))
    
    
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

#MCMC setup ######
#Prior function (don't really need this since drawing from uniform distributions)
  logPrior <- function(param_set){
    
    logPriorBeta <- log(dunif(param_set[["beta"]], min = 1e-8, max = 1e-2))
    logPriorLamda1 <- log(dunif(param_set[["lamda"]], min = 1e-8, max = 1e-2))
    logPriorLamda2 <- log(dunif(param_set[["lamda"]], min = 1e-8, max = 1e-2))
    
    return(logPriorBeta + logPriorLamda1 + logPriorLamda2)
    
  }

## Likelihood function for a single data point:
  pointLogLike <- function(obs_mu = measured_ws, #Measured worm burdens
                           obs_sd = w_errors,       #Measured worm burden uncertainty
                           modOut) {                #Model output W
    log(Prob_gaussian(y=modOut, mu=obs_mu, sd=obs_sd))
  }

## Likelihood function for all 4 data points:
  trajLogLike <- function(param_set) {
    run_to_eqbm <- ode(nstart, time, schisto_mod_pdd_add_ndds, param_set)[max(time)/30,]
    
  	traj <- data.frame(ode(run_to_eqbm[c(2:6)], traj_time, schisto_mod_pdd_add_ndds_dyna,  
                           param_set, events = list(data = mdas)))
  	
  	observed <- traj$Wt[timepoints + 2] * 0.5 * mapply(phi_Wk, traj$Wt[timepoints + 2], measured_ks)
  	
  	logLike <- 0
  	for (i in 1:length(observed)) {
  		logLike <- logLike + pointLogLike(obs_mu = measured_ws[i], obs_sd = w_errors[i], modOut = observed[i])
  	}
  	return(logLike)
  }
  
##Plot simulation 
  plottraj <- function(param_set) {
    run_to_eqbm <- ode(nstart, time, schisto_mod_pdd_add_ndds, param_set)[max(time)/30,]
    
  	traj <- data.frame(ode(run_to_eqbm[c(2:6)], traj_time, schisto_mod_pdd_add_ndds_dyna,  
                           param_set, events = list(data = mdas)))
  	
  	traj_full <- rbind(run_to_eqbm, traj)
  	traj_full$time[1] <- 0
  	
  	traj_full$W <- (covrg*traj_full$Wt) + ((1-covrg)*traj_full$Wu)
  	
  	traj_full %>%
  	  select(time, Wt, Wu) %>%
  	  mutate(W = covrg*Wt + (1-covrg)*Wu) %>%
  	  gather(key = "Treatment", value = "Worm_Burden", Wt, Wu, W) %>%
  	  ggplot(aes(x = time, y = Worm_Burden, lty = Treatment)) + theme_bw() + geom_line() 
  	  
  }	

## Posterior function:
  logPosterior <- function(param_set) {
    ## Calculate the log prior (logPrior) for the vector of model
    ## parameters (theta).
  	logPrior <- logPrior(param_set)
    
    ## Calculate the log likelihood (logLike) of the data given theta, the
    ## prevalence data (PrevData), and the initial values of the state
    ## variables (initState).
  	logLike <- trajLogLike(param_set)
    
    logPosterior <- logPrior + logLike
  	return(logPosterior)
  }

#Test with an example parameter set #############
  par_test <- params
  par_test["beta"] <- 1e-6
  par_test["lamda1"] <- 1e-6
  par_test["lamda2"] <- 1e-4

  #logPrior(par_test)
  #pointLogLike(modOut = 1.8)
  #trajLogLike(par_test)
  #plottraj(par_test)
  #ogPosterior(par_test)
  
## Metropolis-Hastings algorithm for searching parameter space
mcmcMH <- function(posterior, initParams, proposal_up, proposal_lo, numIterations) {
  
  # Evaluate the function "posterior" at initial parameter set, and assign to a
  # variable called posteriorThetaCurrent.
  posteriorThetaCurrent <- posterior(initParams)
  
  # Initialise variables to store the current value of parameters, the
  # vector of sample values, and the number of accepted proposals.
  thetaCurrent <- initParams
  samples <- c(initParams[["beta"]], initParams[["lamda1"]], initParams[["lamda2"]])
  accepted <- 0
  
  # Run the MCMC algorithm for numIterations iterations.
  for (i in 1:numIterations) {
    
    # Draw a new parameter set from a Gaussian proposal distribution and
    # assign this to a variable called thetaProposed.
    thetaProposed <- thetaCurrent
    
    thetaProposed[["beta"]] <- runif(1, min = proposal_lo[1], max = proposal_up[1])
    thetaProposed[["lamda1"]] <- runif(1, min = proposal_lo[2], max = proposal_up[2])
    thetaProposed[["lamda2"]] <- runif(1, min = proposal_lo[3], max = proposal_up[3])
    
    # Evaluate the log) posterior function at the proposed parameter set
    # value and assign to a variable called 
    # posteriorThetaProposed.
    posteriorThetaProposed <- posterior(thetaProposed)
    
    # Compute the Metropolis-Hastings (log) acceptance
    # probability and assign to a variable called
    # logAcceptance.
    logAcceptance <- posteriorThetaProposed - posteriorThetaCurrent
    
    # Draw a random number uniformly-distributed between 0 and 1
    # using "runif" and assign to a variable called randNum.
    randNum <- runif(n = 1, min = 0, max = 1)
    
    # Use the random number and the acceptance probability to 
    # determine if thetaProposed will be accepted.
    if (randNum < exp(logAcceptance)) {
      
      # If accepted, change the current value of theta to the
      # proposed value of theta.
      thetaCurrent <- thetaProposed
      
      # And update the current value of the posterior 
      # function.
      posteriorThetaCurrent <- posteriorThetaProposed
      
      # And update number of accepted proposals.
      accepted <- accepted + 1
    }
    
    # Add the current theta to the vector of samples.
    samples <- c(samples, thetaCurrent[["beta"]], thetaCurrent[["lamda1"]], thetaCurrent[["lamda2"]])
    
    # Print the current state of chain and acceptance rate.
    cat("iteration:", i, "chain:", thetaCurrent[["beta"]], thetaCurrent[["lamda1"]], thetaCurrent[["lamda2"]],
        "acceptance rate:", accepted / i, "\n")
  }
  return(samples)
}
 
# Running the MCMC algorithm to vary the parameters beta, lamda1, and lamda2 parameter sets
#mcmcTrace <- mcmcMH(posterior = logPosterior, # posterior distribution
#                    initParams = par_test, # intial parameter guess
#                    proposal_lo = rep(1e-7, 3), #upper bound on uniform distribution for each parameter
#                    proposal_up = rep(1e-3, 3), #lower bound on uniform distribution for each parameter
#                    numIterations = 2000) # number of iterations

#trace <- matrix(mcmcTrace, ncol = 2, byrow = T)
  
# Use the package "coda" to convert the trace into this format:
#trace <- mcmc(trace)
#plot(trace)
#summary(trace)

# Running the MCMC algorithm to vary the parameters beta, lamda1, and lamda2 parameter sets
mcmcTrace <- mcmcMH(posterior = logPosterior, # posterior distribution
                    initParams = par_test, # intial parameter guess
                    proposal_lo = rep(1e-7, 3), #upper bound on uniform distribution for each parameter
                    proposal_up = rep(1e-3, 3), #lower bound on uniform distribution for each parameter
                    numIterations = 10000) # number of iterations

trace <- matrix(mcmcTrace, ncol = 3, byrow = T)

# Use the package "coda" to convert the trace into this format:
trace <- mcmc(trace)

## Remove the first 1000 iterations to allow for burn-in:
traceBurn <- trace[-(1:1000),]
traceBurn <- mcmc(traceBurn)
plot(traceBurn)
summary(traceBurn)

## Check for autocorrelation:
autocorr.plot(traceBurn)

## Subsampling to account for autocorrelation:
subsample <- 10
traceBurnbeta <- traceBurn[,1]
traceBurnlamda1 <- traceBurn[,2]
traceBurnlamda2 <- traceBurn[,3]
traceBurnAndThinbeta <- traceBurnbeta[seq(1,length(traceBurnbeta),subsample)]
traceBurnAndThinlamda1 <- traceBurnAndThinlamda1[seq(1,length(traceBurnlamda1),subsample)]
traceBurnAndThinlamda2 <- traceBurnAndThinlamda2[seq(1,length(traceBurnlamda2),subsample)]
traceBurnAndThin <- cbind(traceBurnAndThinbeta, traceBurnAndThinlamda1, traceBurnAndThinlamda2)
traceBurnAndThin <- mcmc(traceBurnAndThin)
plot(traceBurnAndThin)
summary(traceBurnAndThin) 

## Check for autocorrelation:
autocorr.plot(traceBurnAndThin)
