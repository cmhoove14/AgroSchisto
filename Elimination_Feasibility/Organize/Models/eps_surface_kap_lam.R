#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

##################################################################################################

require(deSolve)
require(graphics)
require(parallel)
require(rootSolve)

no.cores = detectCores() - 1

source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/Reff_fn_lib_CH.R")
source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/lib_schistoModels_DDandNODD_CH.R")

#Use simplified r0 expression to determine range of lamda to test #########
parameters_2pops_mda_Chris1["beta"]<-1.63e-6      #from best fit params with low transmission
cov = 0.8
parameters_2pops_mda_Chris1["cov"] = cov          # 80% coverage
params<-parameters_2pops_mda_Chris1


get_Ro_lamda<-function(lamda) #variable params to be sampled
{ 
  f_N<-params["f_N"]
  phi_N<-params["phi_N"]
  z<-params["z"]
  mu_N<-params["mu_N"]
  sigma<-params["sigma"]
  mu_I<-params["mu_I"]
  mu_W<-params["mu_W"]
  m<-params["m"]
  H<-params["H"]
  mu_H<-params["mu_H"]
  beta<-params["beta"]
  
  
  #Equilibrium estimate of N given snail params
  N_eq = max(uniroot.all(f = function(y){(f_N)*(1-y/phi_N) - mu_N
  }, 
  c(0, as.numeric(phi_N))))
  
  if(N_eq < 0){
    N_eq = 0
  }
  
  Ro_est <- sqrt((0.5*beta*H*N_eq*lamda*sigma)/((mu_W+mu_H)*(mu_N+sigma)*(mu_N+mu_I)))
  
  Ro_est
  
}

init.lam = seq(1e-5,4e-4, length.out = 100)

windows()
plot(init.lam, sapply(init.lam, get_Ro_lamda, simplify = T), pch = 16, cex = 0.7)

#we'll go with a range of lambdas from 1e-4 to 3e-4 corresponding to r0 ~1-2

#parameter ranges, values and vectors #########
  sim.range = 50                                    #number of values within test range to simulate
  lam.range = seq(1e-4, 3e-4, length.out = sim.range)  #transmission intensity
  kap.range = seq(0, 2, length.out = sim.range)            #Pos. density dependence

#temp vectors to fill
  w.pre = as.numeric()
  w.pos = as.numeric()
  bbr = as.numeric()

#events data frame for MDA introduction (will be the same every time) #######
  eff = 0.99                                        # 99% efficacy
  mda.years = c(1:20) # annual MDA for 20 years
  
  mda.events = data.frame(var = rep('Wt', length(mda.years)),
                          time = c(mda.years*365),
                          value = rep((1 - eff), length(mda.years)),
                          method = rep('mult', length(mda.years)))
  
  teq<-seq(from=0, to=200*365, by=10)
  nstart=c(S=3892,E=3750,I=1200, Wt=50, Wu=50)
  time = c(1:(365*22))
  
#Function to estimate epsilon given transmisison intensity (lambda) and PDD (kappa) #########
eps.est = function(lam, kap){
  
  params["lamda"] = lam 
  params['k'] = kap

  output<-as.data.frame(ode(nstart,teq,schisto_halstead_2pops_mda,params))
  
  eqbm2 = c(S = output[dim(output)[1], 2], 
            E = output[dim(output)[1], 3], 
            I = output[dim(output)[1], 4], 
            Wt = output[dim(output)[1],5], 
            Wu = output[dim(output)[1],6])
  
  mda.sim = as.matrix(ode(eqbm2, time, schisto_halstead_2pops_mda, params,
                      events = list(data = mda.events)))
  
  for(y in mda.years){
    w.pre[y] = cov*mda.sim[(y*365), 5] + (1-cov) * mda.sim[(y*365), 6]
    w.pos[y] = cov*mda.sim[(y*365+1), 5] + (1-cov) * mda.sim[(y*365+1), 6]
  }
  
  for(b in mda.years[-20]){
    bbr[b] = (1/w.pos[b])*(w.pre[b+1] - w.pos[b])
  }
  
  eps = lm(bbr ~ c(1:19))$coefficients[2]
  
  return(as.numeric(eps))
    
}

eps.fill = matrix(ncol = sim.range, nrow = sim.range)    #matrix to fill with eps estimates

clusteps = makeCluster(no.cores)
clusterExport(cl = clusteps, 
              varlist = c('params','sim.range', 'lam.range', 'kap.range', 'ode', 'mda.years',
                          'mda.events', 'nstart', 'teq', 'schisto_halstead_2pops_mda',
                          'time', 'cov', 'w.pre', 'w.pos', 'bbr'))  

for(l in 1:sim.range){
  eps.fill[l,] = parSapply(clusteps, lam.range, eps.est, kap = kap.range[l], simplify = T)
  print(l)
} #lambda varies across a row; kappa varies down a column

stopCluster(clusteps)

windows()

persp(y = lam.range, ylim = range(lam.range), x = kap.range, xlim = range(kap.range),
      z = eps.fill, ticktype = 'detailed', nticks = 4, 
      xlab = 'Pos. Density Dependence',
      ylab = 'Transmission Intensity',
      zlab = 'Elimination Feasibility Estimator',
      phi = 0, theta = 90, shade = 0.4, col = 'lightblue')
