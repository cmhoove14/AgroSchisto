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

source("~/ElimFeas_StochMod/PostReview/Refs/schisto_mods_pdd_nopdd_savio.R")

load("~/ElimFeas_StochMod/PostReview/Outputs/model_fit_profile_likelihood_parameters.Rdata")
load("~/ElimFeas_StochMod/PostReview/Refs/best_fit_params.Rdata")

#Use simplified r0 expression to determine range of lamda to test #########
covrg <- 0.8
params[["cov"]] <- covrg

#parameter ranges, values and vectors #########
  sim.range = 50                                    #number of values within test range to simulate
  lam.range = seq(min(fin_pars95ci$lamda_twa), max(fin_pars95ci$lamda_twa), length.out = sim.range)  #transmission intensity
  kap.range = seq(0, 2, length.out = sim.range)            #Pos. density dependence

#temp vectors to fill
  w.pre = as.numeric()
  w.pos = as.numeric()
  bbr = as.numeric()

#events data frame for MDA introduction (will be the same every time) #######
  eff = 0.94                                        # 94% efficacy
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

  output<-as.data.frame(ode(nstart,teq,schisto_mod_pdd_add_ndds,params))
  
  eqbm2 = c(S = output[dim(output)[1], 2], 
            E = output[dim(output)[1], 3], 
            I = output[dim(output)[1], 4], 
            Wt = output[dim(output)[1],5], 
            Wu = output[dim(output)[1],6])
  
  mda.sim = as.matrix(ode(eqbm2, time, schisto_mod_pdd_add_ndds, params,
                      events = list(data = mda.events)))
  
  for(y in mda.years){
    w.pre[y] = covrg*mda.sim[(y*365), 5] + (1-covrg) * mda.sim[(y*365), 6]
    w.pos[y] = covrg*mda.sim[(y*365+1), 5] + (1-covrg) * mda.sim[(y*365+1), 6]
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
                          'mda.events', 'nstart', 'teq', 'schisto_mod_pdd_add_ndds',
                          'time', 'covrg', 'w.pre', 'w.pos', 'bbr', 'f_Wgk', 'R_Wv', 'phi_Wk'))  

for(l in 1:sim.range){
  eps.fill[l,] = parSapply(clusteps, lam.range, eps.est, kap = kap.range[l], simplify = T)
  print(l)
} #lambda varies across a row; kappa varies down a column

stopCluster(clusteps)

save(eps.fill, file = "~/ElimFeas_StochMod/PostReview/Outputs/eps_surface_lamRange_kRange.Rdata")