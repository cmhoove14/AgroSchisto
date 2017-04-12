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

source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/lib_schistoModels_DDandNODD_CH.R")
source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/Reff_fn_lib_CH.R")

outputfile<-"C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/shortlist_51.csv"
shortlist_first100<-as.data.frame(read.csv(file=outputfile, header=TRUE, sep=","))
sort_ind<-order(shortlist_first100$R0, decreasing=FALSE)
shortlist_first100<-shortlist_first100[sort_ind,]
  shortlist_first100 = shortlist_first100[c(1:100),] #Get rid of other parameter estimates

#transmission parameter ranges to sample
par.sims = 10                                         #number of values within vector
  lam.range = seq(1e-4, 3e-4, length.out = par.sims)  #transmission intensity range
  kap.range = seq(0, 2, length.out = par.sims)        #Pos. density dependence range
  
#Set parameters #######
cov = 0.8
eff = 0.99
mda.years = c(2:21)

params = as.list(parameters_2pops_mda_Chris1)
params$beta = shortlist_first100$beta[1]
params$lamda = shortlist_first100$lamda.twa[1]
  params$cov = cov

start = c(S = 5000, # susceptible humans
          E = 2000, # infected humans
          I = 500, # infected humans
          Wt = 72,
          Wu = 72) # recovered (and immune) humans

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

#Test #######
lam = lam.range[par.sims/2]
kap = kap.range[par.sims/2]         #test parameter values

params['k'] = kap
params['lamda'] = lam

start2 = as.numeric(as.data.frame(ode(start,seq(0,200*365,30),
                    schisto_halstead_2pops_mda,params))[length(seq(0,200*365,30)), c(2:6)])

stoch.sims = 10
years = 61

year.days = as.numeric()
for(i in 1:20){
  year.days[i] = 365*i + (i-1)
}

plotters = c(1:stoch.sims)
#objects to fill
fill = list()
pe1 = as.numeric()
w.pre = as.numeric()   
w.pos = as.numeric()
bbr = as.numeric()
eps = as.numeric()     

#function to simulate transmission over 61 years (1 year transmission spin up, 20 yrs MDA, 40 yrs recovery)

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
    init = setNames(as.numeric(fill[[f-1]][dim(fill[[f-1]])[1],c(2:6)]), 
                    colnames(fill[[f-1]])[c(2:6)]) #reset initial states
    
    #init[4] = round(init[4]* (1-eff))  #NO MDA
    
    fill[[f]] = ssa.adaptivetau(init, transitions, 
                                sfx, params, tf=365) #stochastic sim for a year
    
    fill[[f]][,1] = fill[[f]][,1] + (365*(f-1)+(f-1))    #adjust time
  }
  
  matfin = do.call(rbind,fill)
  
  matfin = cbind(matfin, Wm = cov*matfin[,5] + (1-cov)*matfin[,6])
  
  if(sim %in% plotters){
    lines(matfin[,1], matfin[,7], col = sim)
  }
    
  
  if(sum(matfin[max(which(!is.na(matfin[,7]))),][3:7]) == 0){ 
    pe1 = 1   #if no exposed, infected snails and no adult worms, elimination = 1
  } else {
    pe1 = 0   #else, elimination = 0
  }
  
  w.pre = matfin[ , 7][matfin[ , 1] %in% year.days]     #w.pre values
  w.pos = matfin[ , 7][matfin[ , 1] %in% (year.days+1)] #w.pos values
  
  bbr = (1/w.pos[c(1:19)])*(w.pre[c(2:20)] - w.pos[c(1:19)])
  
  #Estimate epsilon for each sim    
    eps = lm(bbr ~ c(1:19))$coefficients[2]
  
  return(c(k, lam, sim, pe1, eps))
  
}

#Test
plot(0, 0, pch = 1, cex = 0.001, xlab = 'time', ylab = 'W',
     xlim = c(0, 365*61), ylim = c(0, start2[5]))

test = sapply(c(1:stoch.sims), stoch.sim, init = start2, k = kap, lam = lam, simplify = T)


#Run simulations #######
#Necesssary parameter values: transmission, PDD, initial state variables  
par.sims = 2
  lam.range = seq(1e-4, 3e-4, length.out = par.sims)  #transmission intensity range
  kap.range = seq(0, 2, length.out = par.sims)        #Pos. density dependence range

#get matrix of parameter values with equilibrium state variables    
#par.mat = read.csv('Elimination_Feasibility/eq_vals_for_trans_pars.csv')
par.mat = matrix(ncol = 7, nrow = par.sims^2)         #initialize matrix

#fill first two matrix columns with transmission intensity and PDD parameter
par.mat[,c(1:2)] = cbind(rep(lam.range, times = par.sims), rep(kap.range, each = par.sims)) 

#fill next five columns with initial state variable values given trans. intensity & PDD
for(i in 1:nrow(par.mat)){
  params['lamda'] = par.mat[i,1]
  params['k'] = par.mat[i,2]
  
  par.mat[i,c(3:7)] = as.numeric(as.data.frame(ode(start,seq(0,200*365,30),
                                                   schisto_halstead_2pops_mda,params))[length(seq(0,200*365,30)), c(2:6)])
  
  print(i)
}


stoch.sims = 1000  #number of simulations for ach parameter set
plot.path = 'Elimination_Feasibility/plots/simages/'

#Final values array to fill with p(e), eps, eps.sd
fill.arr = array(data = NA, dim = c(par.sims, par.sims, 2, stoch.sims))   

#run all simulations######
for(s in 1:nrow(par.mat)){

#plot annotations  
  plotters = round(runif(20, min = 1, max = stoch.sims))        #random sim numbers to plot
  plot.name = paste('sim.l', which(par.mat[s,1] == lam.range), 
                    '.k', which(par.mat[s,2] == kap.range), '.png', sep = '')
  plot.title = paste('lambda = ', par.mat[s,1], 
                     ' , kappa = ', par.mat[s,2], sep = '')
  
  png(filename = paste(plot.path, plot.name, sep = ''), width = 600, height = 450)  
  plot(0, 0, pch = 1, cex = 0.001, xlab = 'time', ylab = 'W',
       xlim = c(0, 365*61), ylim = c(0, par.mat[s,7]), main = plot.title) #initialize blank plot
  

#Run sims for parameter set  
  fill.arr[which(par.mat[s,1] == lam.range), which(par.mat[s,2] == kap.range), c(1:2), ] = 
                                           sapply(c(1:stoch.sims), stoch.sim, init = par.mat[s,c(3:7)], 
                                                                             lam = par.mat[s,1], 
                                                                             k = par.mat[s,2], simplify = T)[c(4,5),]
  
  dev.off()  

  print(s)
}  

#post process ########
#get p(e) for each parameter set
pe = as.numeric()
eps.mean = as.numeric()
eps.sd = as.numeric()
for(s in 1:nrow(par.mat)){
  pe[s] = sum(fill.arr[which(par.mat[s,1] == lam.range), which(par.mat[s,2] == kap.range), 1, ]) / stoch.sims
  eps.mean[s] = mean(fill.arr[which(par.mat[s,1] == lam.range), which(par.mat[s,2] == kap.range), 2, ])
  eps.sd[s] = sd(fill.arr[which(par.mat[s,1] == lam.range), which(par.mat[s,2] == kap.range), 2, ])
}

plot(eps.mean, pe, pch = 16, cex = 0.6,
     xlab = expression(epsilon), ylab = expression(italic('P(e)')))