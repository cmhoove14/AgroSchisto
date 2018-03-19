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
library(deSolve)

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
  W = covrg*Wt+(1-covrg)*Wu
  
  return(c(p$f_N * (1-N/p$C) * (S + E),   #Snail birth
           p$mu_N * S,        #Susceptible snail death
           p$beta * 0.5 * W * p$H * S * phi_Wk(W = W, k = p$k) * f_Wgk(W, p$gam, p$k),  #Snail exposure
           p$mu_N * E,       #Exposed snail dies
           p$sigma * E,      #Exposed snail becomes infected
           (p$mu_N + p$mu_I) * I,   #Infected snail dies
           p$lamda * I * R_Wv(W, p$v),        #infected snail produces adult worm
           (p$mu_W + p$mu_H) * Wt,    #Adult worm in treated population dies
           (p$mu_W + p$mu_H) * Wu))    #Adult worm in untreated population dies
}

#Full stochastic simulation model returning P(e) and Eps
stoch.sim = function(init, k, lam, sim){
  params['k'] = k
  params['lamda'] = lam
  
  init1 = setNames(as.numeric(round(init)), c('S', 'E', 'I', 'Wt', 'Wu'))
  
  set.seed(sim)
  
  fill[[1]] = ssa.adaptivetau(init1, transitions, 
                              sfx, params, tf=365)    #simulate 1 year of transmission
  
  for(m in 2:21){    #simulate 20 years of MDA
    init = setNames(as.numeric(fill[[m-1]][dim(fill[[m-1]])[1],c(2:6)]), 
                    colnames(fill[[m-1]])[c(2:6)]) #reset initial states
    
    init[4] = round(init[4]* (1-eff))  #apply MDA
    
    fill[[m]] = ssa.adaptivetau(init, transitions, 
                                sfx, params, tf=365) #stochastic sim for a year
    
    fill[[m]][,1] = fill[[m]][,1] + (365*(m-1)+(m-1))    #adjust time
  }
  
  for(f in 22:years){
    
    if(sum(as.numeric(fill[[f-1]][dim(fill[[f-1]])[1],c(3:6)])) > 0){  # Stop simulation if elimination has occurred
    
    init = setNames(as.numeric(fill[[f-1]][dim(fill[[f-1]])[1],c(2:6)]), 
                    colnames(fill[[f-1]])[c(2:6)]) #reset initial states
    
    fill[[f]] = ssa.adaptivetau(init, transitions, 
                                sfx, params, tf=365) #stochastic sim for a year
    
    fill[[f]][,1] = fill[[f]][,1] + (365*(f-1)+(f-1))    #adjust time
  } else {
      
      fill[[f]] = matrix(fill[[f-1]][dim(fill[[f-1]])[1],], ncol = 6, byrow = F)
      
    }
}
  
  matfin = do.call(rbind,fill)
  
  matfin = cbind(matfin, Wm = covrg*matfin[,5] + (1-covrg)*matfin[,6])
  
  if(sum(matfin[max(which(!is.na(matfin[,7]))),][3:7]) == 0){ 
    pe1 = 1   #if no exposed, infected snails and no adult worms, elimination = 1
  } else {
    pe1 = 0   #else, elimination = 0
  }
  
  w.pre = matfin[ , 7][matfin[ , 1] %in% year.days]     #w.pre values
  prev.pre = Prevalence(w.pre, params['k']$k)

  w.pos = matfin[ , 7][matfin[ , 1] %in% (year.days+1)] #w.pos values
  prev.pos = Prevalence(w.pos, params['k']$k)

  
  bbr.w = (1/w.pos[c(1:19)])*(w.pre[c(2:20)] - w.pos[c(1:19)])
    bbr.w[which(is.infinite(bbr.w))] <- NA

  bbr.prev = (1/prev.pos[c(1:19)])*(prev.pre[c(2:20)] - prev.pos[c(1:19)])
    bbr.prev[which(is.infinite(bbr.prev))] <- NA


  #Estimate epsilon for each sim  
  if(sum(!is.na(bbr.w)) >=3){
    
    eps.w = lm(bbr.w ~ c(1:19), na.action = na.exclude)$coefficients[2]
    eps.prev = lm(bbr.prev ~ c(1:19), na.action = na.exclude)$coefficients[2]
    
  } else {
    eps.w = NA
    eps.prev = NA
  }
    
  
  return(c(k, lam, sim, pe1, eps.w, eps.prev))
  
}

#Full stochastic simulation model with observation noise added to estimation of worm burden returning P(e) and Eps
stoch.sim.noise = function(init, k, lam, sim){
  params['k'] = k
  params['lamda'] = lam
  
  init1 = setNames(as.numeric(round(init)), c('S', 'E', 'I', 'Wt', 'Wu'))
  
  set.seed(sim)
  
  fill[[1]] = ssa.adaptivetau(init1, transitions, 
                              sfx, params, tf=365)    #simulate 1 year of transmission
  
  for(m in 2:21){    #simulate 20 years of MDA
    init = setNames(as.numeric(fill[[m-1]][dim(fill[[m-1]])[1],c(2:6)]), 
                    colnames(fill[[m-1]])[c(2:6)]) #reset initial states
    
    init[4] = round(init[4]* (1-eff))  #apply MDA
    
    fill[[m]] = ssa.adaptivetau(init, transitions, 
                                sfx, params, tf=365) #stochastic sim for a year
    
    fill[[m]][,1] = fill[[m]][,1] + (365*(m-1)+(m-1))    #adjust time
  }
  
  for(f in 22:years){
    
    if(sum(as.numeric(fill[[f-1]][dim(fill[[f-1]])[1],c(3:6)])) > 0){  # Stop simulation if elimination has occurred
    
    init = setNames(as.numeric(fill[[f-1]][dim(fill[[f-1]])[1],c(2:6)]), 
                    colnames(fill[[f-1]])[c(2:6)]) #reset initial states
    
    fill[[f]] = ssa.adaptivetau(init, transitions, 
                                sfx, params, tf=365) #stochastic sim for a year
    
    fill[[f]][,1] = fill[[f]][,1] + (365*(f-1)+(f-1))    #adjust time
  } else {
      
      fill[[f]] = matrix(fill[[f-1]][dim(fill[[f-1]])[1],], ncol = 6, byrow = F)
      
    }
}
  
  matfin = do.call(rbind,fill)
  
  matfin = cbind(matfin, Wm = covrg*matfin[,5] + (1-covrg)*matfin[,6])
  
  if(sum(matfin[max(which(!is.na(matfin[,7]))),][3:7]) == 0){ 
    pe1 = 1   #if no exposed, infected snails and no adult worms, elimination = 1
  } else {
    pe1 = 0   #else, elimination = 0
  }
  
#Estimates of W_pre and W_pos with observation noise resulting from resampling  
for(i in 1:length(year.days)){
    w.pre[i] = mean(rnbinom(params$H, 
                            mu = matfin[ , 7][matfin[ , 1] %in% year.days[i]],
                            size = k)) #w.pre values sampled from 300 individuals
    
    prev.pre[i] = Prevalence(w.pre[i], params['k']$k)
  
    w.pos[i] = mean(rnbinom(params$H, 
                            mu = matfin[ , 7][matfin[ , 1] %in% (year.days[i] + 1)],
                            size = k)) #w.pos values sampled from 300 individuals
    
    prev.pos[i] = Prevalence(w.pos[i], params['k']$k)

  }
  
  bbr.w = (1/w.pos[c(1:19)])*(w.pre[c(2:20)] - w.pos[c(1:19)])
    bbr.w[which(is.infinite(bbr.w))] <- NA
  
  bbr.prev = (1/prev.pos[c(1:19)])*(prev.pre[c(2:20)] - prev.pos[c(1:19)])
    bbr.prev[which(is.infinite(bbr.prev))] <- NA

  #Estimate epsilon for each sim    
  if(sum(!is.na(bbr.w)) >=3){
    
    eps.w = lm(bbr.w ~ c(1:19), na.action = na.exclude)$coefficients[2]
    eps.prev = lm(bbr.prev ~ c(1:19), na.action = na.exclude)$coefficients[2]
    
  } else {
    eps.w = NA
    eps.prev = NA
  }

  return(c(k, lam, sim, pe1, eps.w, eps.prev))
  
}
