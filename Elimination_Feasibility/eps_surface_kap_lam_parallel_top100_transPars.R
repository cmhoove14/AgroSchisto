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

#Try simulations with different transmission intensities
library(parallel)
require(deSolve)
  no.cores = detectCores() - 1

outputfile<-"C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/shortlist_51.csv"
  best100<-as.data.frame(read.csv(file=outputfile, header=TRUE, sep=","))
    sort_ind<-order(best100$negLL, decreasing=T)
    best100<-best100[sort_ind,]
      best100 = best100[c(1:100),] #Get rid of other parameter estimates; keep 100 parameter sets with max(negLL)

beta.vec = rep(best100$beta, times = 100)
lam.vec = rep(best100$lamda.twa, times = 100)
r0.vec = rep(best100$R0, times = 100)
kap.vec = rep(seq(0,5, length.out = 100), each = 100)
eps.vec = 0
  vecs = cbind(beta.vec, lam.vec, r0.vec, kap.vec, eps.vec)  

eps.est.beta = function(lam, beta, kap){
  params<-parameters_2pops_mda_Chris1
  
  params["lamda"] = lam 
  params["beta"] = beta
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


clusteps = makeCluster(no.cores)
  clusterExport(cl = clusteps, varlist = c('schisto_halstead_2pops_mda','parameters_2pops_mda_Chris1','ode'))  

  vecs[,5] = clusterMap(clusteps, 
                        fun = eps.est.beta, 
                        lam = lam.vec, 
                        beta = beta.vec, 
                        kap = kap.vec, 
                        SIMPLIFY = TRUE)

stopCluster(clusteps)
