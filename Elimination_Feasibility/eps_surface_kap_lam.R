#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#This code adapted from Arathi's "schisto_halstead_2pops_mda_bouncebackRate_cloudParam.R"

# This code seeks to test the bounce back rate concept. Given the worm burden from a model, can we determine if
# it has positive DD and hence possible to eliminate? Or has no DD and is harder to eliminate.
# We first start with a deterministic model
# We will then fit it to the Lampsar II data to get a collection of parameter values that fit max likelihood bounds
# Each parameter value will give a unique W for pre and post mda in this model and a unique bounce back rate curve
# If there are 100 parameter values for example, how many have Reff<1 after 10 rounds?
##################################################################################################

require(deSolve)
require(graphics)

source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/Reff_fn_lib_CH.R")
source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/lib_schistoModels_DDandNODD_CH.R")

#parameter ranges
  sim.range = 20                                    #number of values within test range to simulate
  lam.range = seq(1.2e-4, 1.7e-4, length.out = sim.range)  #transmission intensity
  kap.range = seq(0, 2, length.out = sim.range)            #Pos. density dependence

#temp vectors to fill
  w.pre = as.numeric()
  w.pos = as.numeric()
  bbr = as.numeric()

#parameter assignment
  parameters_2pops_mda_Chris1["beta"]<-1.63e-6      #from best fit params with low transmission
  cov = 0.8
  parameters_2pops_mda_Chris1["cov"] = cov          # 80% coverage
  eff = 0.99                                        # 99% efficacy
  time = c(1:(365*22)) #1 year run up, 20 years annual mda followed by 1 extra year

#events data frame for MDA introduction (will be the same every time)
  mda.years = c(1:20) # annual MDA for 20 years
  
  mda.events = data.frame(var = rep('Wt', length(mda.years)),
                          time = c(mda.years*365),
                          value = rep((1 - eff), length(mda.years)),
                          method = rep('mult', length(mda.years)))
  
#Function to estimate epsilon given transmisison intensity (lambda) and PDD (kappa)
eps.est = function(lam, kap){
  params<-parameters_2pops_mda_Chris1
  
  params["lamda"] = lam 
  params['k'] = kap

  teq<-seq(from=0, to=200*365, by=1)
  nstart=c(S=3892,E=3750,I=1200, Wt=50, Wu=50)
  
  output<-as.data.frame(ode(nstart,teq,schisto_halstead_2pops_mda,params))
  
  baseline<-output[dim(output)[1],]
  
  eqbm = baseline[,c(2:6)]
  
  eqbm2 = c(S=eqbm$S, 
            E=eqbm$E, 
            I=eqbm$I, 
            Wt=eqbm$Wt, 
            Wu=eqbm$Wu)
  
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

for(l in 1:sim.range){
  eps.fill[l,] = sapply(lam.range, eps.est, kap = kap.range[l], simplify = T)
  print(l)
} #lambda varies across a row; kappa varies down a column

windows()

persp(y = lam.range, ylim = range(lam.range), x = kap.range, xlim = range(kap.range),
      z = eps.fill, ticktype = 'detailed', nticks = 4, 
      xlab = 'Pos. Density Dependence',
      ylab = 'Transmission Intensity',
      zlab = 'Elimination Feasibility Estimator',
      phi = 30, theta = 30, shade = 0.4, col = 'lightblue')
