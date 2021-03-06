#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


#Load packages and other files ##########
library(adaptivetau)
library(parallel)

n.cores = detectCores() - 1

min.sim = 1
max.sim = 500

source("~/ElimFeas_StochMod/lib_schistoModels_DDandNODD_CH.R")


#Set parameters #######
cov = 0.8
eff = 0.99
  mda.years = c(2:21)

params = as.list(parameters_2pops_mda_Chris1)
  params$cov = cov
  params$beta = 1.63e-06


#adaptivetau model ############
transitions = list(
  c(S = 1),             #New snail born
  c(S = -1),            #Susceptible snail dies
  c(S = -1, E = 1),     #Susceptible snail becomes exposed
  c(E = -1),            #Exposed snail dies
  c(E = -1, I = 1),     #Exposed snail becomes Infected
  c(I = -1),            #Infected snail dies
  c(Wt = 1, Wu = 1),    #Infected snail emits cercaria that produces an adult worm
  c(Wt = -1),           #Adult worm in the treated population dies
  c(Wu = -1))           #Adult worm in the untreated population dies

sfx <- function(x, p, t) {
  S = x['S']
  E = x['E']
  I = x['I']
  N = S + I + E
  Wt = x['Wt']
  Wu = x['Wu']
  W = cov*Wt+(1-cov)*Wu

#model    
  return(c(p$f_N * (1-N/p$phi_N) * (S + E),   #Snail birth
           p$mu_N * S,        #Susceptible snail death
           p$beta * 0.5 *  W * p$H * S * get_phi(W = W, k = p$k),  #Snail exposure
           p$mu_N * E,       #Exposed snail dies
           p$sigma * E,      #Exposed snail becomes infected
           (p$mu_N + p$mu_I) * I,   #Infected snail dies
           p$lamda * I,        #infected snail produces adult worm
           (p$mu_W + p$mu_H) * Wt,    #Adult worm in treated population dies
           (p$mu_W + p$mu_H) * Wu))    #Adult worm in untreated population dies
}

years = 61

#days where MDA is applied
year.days = as.numeric()
for(i in 1:20){
  year.days[i] = 365*i + (i-1)
}

#objects to fill
  fill = list()           #list to fill with simulations
  pe1 = as.numeric()      #elimination binary vector to fill for each sim
  w.pre = as.numeric()    #vector to fill with w_pre values
  w.pos = as.numeric()    #vector to fill with w_pos values
  bbr = as.numeric()      #vector to fill with bbr values from w_pre and w_pos  
  eps = as.numeric()      #vector to fill with epsilon (elim. feas. estimator) vals from slope of bbr

#function to simulate transmission over 61 years (1 year transmission spin up, 20 yrs MDA, 40 yrs recovery)

stoch.sim.noise = function(init, k, lam, sim){
  params['k'] = k
  params['lamda'] = lam
  
  init1 = setNames(as.numeric(round(init)), c('S', 'E', 'I', 'Wt', 'Wu'))
  
  set.seed(sim)
  
#simulate 1 year of transmission  
  fill[[1]] = ssa.adaptivetau(init1, transitions, 
                              sfx, params, tf=365) 
  
#simulate 20 years of MDA  
  for(m in 2:21){    
    init = setNames(as.numeric(fill[[m-1]][dim(fill[[m-1]])[1],c(2:6)]), 
                    colnames(fill[[m-1]])[c(2:6)]) #reset initial states
    
    init[4] = round(init[4]* (1-eff))  #apply MDA
    
    fill[[m]] = ssa.adaptivetau(init, transitions, 
                                sfx, params, tf=365) #stochastic sim for a year
    
    fill[[m]][,1] = fill[[m]][,1] + (365*(m-1)+(m-1))    #adjust time
  }
  
#simulate 40 years no MDA  
  for(f in 22:years){
    
    if(sum(as.numeric(fill[[f-1]][dim(fill[[f-1]])[1],c(3:6)])) > 0){
      
      init = setNames(as.numeric(fill[[f-1]][dim(fill[[f-1]])[1],c(2:6)]), 
                      colnames(fill[[f-1]])[c(2:6)]) #reset initial states
    
       fill[[f]] = ssa.adaptivetau(init, transitions, 
                                sfx, params, tf=365) #stochastic sim for a year
       
        #init[4] = round(init[4]* (1-eff))  #NO MDA

       fill[[f]][,1] = fill[[f]][,1] + (365*(f-1)+(f-1))    #adjust time
    } else {
      
      fill[[f]] = matrix(fill[[f-1]][dim(fill[[f-1]])[1],], ncol = 6, byrow = F)
      
    }
    
  }
  
  matfin = do.call(rbind,fill) #convert list to matrix
  
#add mean worm burden based on coverage partition between treated and untreated 
  matfin = cbind(matfin, Wm = cov*matfin[,5] + (1-cov)*matfin[,6])  
  
  if(sum(matfin[max(which(!is.na(matfin[,7]))),][3:7]) == 0){ 
    pe1 = 1   #if no exposed, infected snails and no adult worms, elimination = 1
  } else {
    pe1 = 0   #else, elimination = 0
  }
  
  
  for(i in 1:length(year.days)){
     w.pre[i] = mean(rnbinom(params$H, 
                       mu = matfin[ , 7][matfin[ , 1] %in% year.days[i]],
                       size = k)) #w.pre values sampled from 300 individuals
  
    w.pos[i] = mean(rnbinom(params$H, 
                         mu = matfin[ , 7][matfin[ , 1] %in% (year.days[i] + 1)],
                         size = k)) #w.pos values sampled from 300 individuals
  }
 
  
  bbr = (1/w.pos[c(1:19)])*(w.pre[c(2:20)] - w.pos[c(1:19)])  #bbr values
    bbr[is.infinite(bbr)] = NA
    
  #Estimate epsilon for each sim    
    eps = lm(bbr ~ c(1:19))$coefficients[2]
  
  return(c(k, lam, sim, pe1, eps))
  
} #end function

#Run simulations #######
#Necesssary parameter values: transmission, PDD, initial state variables  
par.sims = 50
  lam.range = seq(1e-4, 3e-4, length.out = par.sims)  #transmission intensity range
  kap.range = seq(0, 2, length.out = par.sims)        #Pos. density dependence range

#get matrix of parameter values with equilibrium state variables    
par.mat = read.csv('~/ElimFeas_StochMod/eq_vals_for_trans_pars.csv')
  par.mat[c(1:50), 2] = 0.01    #can't sample from neg. binom with k = 0, so replace
  par.mat = par.mat[c(min.sim:max.sim),]
  
stoch.sims = 1000  #number of simulations for each parameter set

#Final values array to fill with p(e), eps, eps.sd
fill.arr = array(data = NA, dim = c(max.sim - min.sim + 1, 4, stoch.sims))

#Make cluster ######
clust = makeCluster(n.cores)
clusterExport(cl = clust, 
              varlist = c('params','lam.range', 'kap.range', 'years', 'stoch.sim.noise', 'par.mat',
                          'mda.years', 'year.days', 'par.sims', 'stoch.sims',
                          'par.mat', 'fill', 'transitions', 'sfx', 'eff',
                          'fill.arr', 'ssa.adaptivetau', 'cov', 'w.pre', 'w.pos',
                          'get_phi', 'fx', 'pe1', 'bbr', 'eps', 'min.sim', 'max.sim'))  


#run all simulations######
for(s in c(min.sim:max.sim)){

#Run sims for parameter set  
  fill.arr[s, , ] = 
  parSapply(clust, c(1:stoch.sims), stoch.sim.noise, init = par.mat[s,c(3:7)], 
                                                lam = par.mat[s,1], 
                                                k = par.mat[s,2], simplify = T)[c(1,2,4,5),]
  
  print(s)
}  

stopCluster(clust)

arr.name = paste('~/ElimFeas_StochMod/wNoise_array', min.sim, '_', max.sim, '.Rdata', sep = '')
save(fill.arr, file = arr.name)
