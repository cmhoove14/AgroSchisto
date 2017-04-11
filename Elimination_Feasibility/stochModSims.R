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
#library(parallel)

#n.cores = detectCores() - 1

source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/lib_schistoModels_DDandNODD_CH.R")
source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/Reff_fn_lib_CH.R")

outputfile<-"C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/shortlist_51.csv"
shortlist_first100<-as.data.frame(read.csv(file=outputfile, header=TRUE, sep=","))
sort_ind<-order(shortlist_first100$R0, decreasing=FALSE)
shortlist_first100<-shortlist_first100[sort_ind,]
  shortlist_first100 = shortlist_first100[c(1:100),] #Get rid of other parameter estimates

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


#function to simulate transmission over 61 years (1 year transmission spin up, 20 yrs MDA, 40 yrs recovery)
stoch.sim = function(init, k, lam, sim){
  params['k'] = k
  params['lamda'] = lam
  
  eq = as.data.frame(ode(init,seq(0,200*365,30),
                         schisto_halstead_2pops_mda,params))[length(seq(0,200*365,30)), c(2:6)]
  
  init1 = setNames(as.numeric(round(eq)), colnames(eq))
  
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
  
  fill.test[c(1:nrow(matfin)), , sim] = cbind(matfin, Wm = cov*matfin[,5] + (1-cov)*matfin[,6])
  
  assign('fill.test', fill.test, envir = .GlobalEnv)
  
}

#Test #######
par.sims = 5                            #number of values within test parameter range to simulate
lam = 3e-4
kap = 0.4

stoch.sims = 10
years = 61

fill = list()
fill.test = array(data = NA, dim = c(years*365*2, 7, stoch.sims))    #array to fill with all simulations

sapply(c(1:stoch.sims), stoch.sim, init = start, k = kap, lam = lam, simplify = T)

  plot(fill.test[ , 1, 1], fill.test[ , 7, 1], type='l', lwd = 2)

  for(i in 2:stoch.sims){
    lines(fill.test[ , 1, i], fill.test[ , 7, i], col = i, lwd = 2)
  }

#Run simulations #######
#Necesssary parameter values  
par.sims = 10
  lam.range = seq(1e-4, 3e-4, length.out = par.sims)  #transmission intensity
  kap.range = seq(0, 2, length.out = par.sims)            #Pos. density dependence
par.mat = cbind(rep(lam.range, times = par.sims), rep(kap.range, each = par.sims))  

year.days = as.numeric()
for(i in 1:20){
  year.days[i] = 365*i + (i-1)
}

stoch.sims = 10  
plot.path = 'Elimination_Feasibility/plots/simages/'

#Final values array to fill
fin.vals = array(data = NA, dim = c(par.sims, par.sims, 3))   #across parameter ranges, save p(e), eps, eps.sd
#Make cluster (attempt)
#clust = makeCluster(n.cores)
#clusterExport(cl = clust, 
#              varlist = c('params','lam.range', 'kap.range', 'ode', 'years', 'stoch.sim',
#                          'start', 'mda.years', 'year.days', 'par.sims', 'stoch.sims',
#                          'plot.path', 'par.mat', 'fill', 'transitions', 'sfx', 'eff',
#                          'fin.vals', 'schisto_halstead_2pops_mda', 'ssa.adaptivetau', 'cov',
#                          'get_phi', 'fx', 'fill.test', 'pe1', 'w.prepos', 'eps.fill'))  


for(s in 1:nrow(par.mat)){
#Things to fill in each parameter set simulation
  fill.test = array(data = NA, dim = c(years*365*3, 7, stoch.sims))    #array to fill with all simulations
  pe1 = as.numeric()        #Vector for number of chains that go extinct
  
  w.prepos = array(data = NA, dim = c(length(year.days), stoch.sims, 3))   #array for W_pre, W_pos, and BBR
  eps.fill = as.numeric()     #vector of epsilon estimates for each parameter set
  
  
#plot annotations  
  plot.name = paste('sim.l', which(par.mat[s,1] == lam.range), 
                    '.k', which(par.mat[s,2] == kap.range), '.png', sep = '')
  plot.title = paste('lambda = ', par.mat[s,1], 
                     ' , kappa = ', par.mat[s,2], sep = '')

#Run sims for parameter set  
  #parSapply(clust, c(1:stoch.sims), stoch.sim, init = start, lam = par.mat[s,1], k = par.mat[s,2], simplify = T)
  sapply(c(1:stoch.sims), stoch.sim, init = start, lam = par.mat[s,1], k = par.mat[s,2], simplify = T)
  
#Save plot of sample simulations    
  png(filename = paste(plot.path, plot.name, sep = ''), width = 600, height = 450)  
  plot(fill.test[ , 1, 1], fill.test[ , 7, 1], type='l', lwd = 2,
         ylab = 'W', xlab = 'time', main = plot.title)
  
  for(i in seq(2, stoch.sims, stoch.sims/10)){
    lines(fill.test[ , 1, i], fill.test[ , 7, i], col = i, lwd = 2)
  }
  dev.off()  
  
#Get probability of elimination (P(e)) as number of chains that lead to extinction out of all chains
  for(p in 1:stoch.sims){
    if(sum(fill.test[max(which(!is.na(fill.test[ , 7, p]))), , p][3:7]) == 0){ 
      pe1[p] = 1   #if no exposed, infected snails and no adult worms, elimination = 1
    } else {
      pe1[p] = 0   #else, elimination = 0
    }
    
  }
  
  fin.vals[which(par.mat[s,1] == lam.range), which(par.mat[s,2] == kap.range), 1] = sum(pe1) / stoch.sims

#extract W_pre and W_pos values from simulations    
  for(w in 1:stoch.sims){
    w.prepos[, w, 1] = fill.test[ , 7, w][fill.test[ , 1, w] %in% year.days]     #w.pre values
    w.prepos[, w, 2] = fill.test[ , 7, w][fill.test[ , 1, w] %in% (year.days+1)] #w.pos values
  }
  
#Get BBR vals across time  
  for(b in 1:stoch.sims){
    w.prepos[, b, 3] = c((1/w.prepos[c(1:19), b, 2])*(w.prepos[c(2:20), b, 1] - w.prepos[c(1:19), b, 2]), 0)
  }  

#Estimate epsilon for each sim    
  for(e in 1:stoch.sims){
    eps.fill[e] = lm(w.prepos[c(1:19), e, 3] ~ c(1:19))$coefficients[2]
  }

#Save mean and sd of epsilon estimates for parameter vals across simulations    
  fin.vals[which(par.mat[s,1] == lam.range), which(par.mat[s,2] == kap.range), 2] = mean(eps.fill)
  fin.vals[which(par.mat[s,1] == lam.range), which(par.mat[s,2] == kap.range), 3] = sd(eps.fill)

#Save array to global environment    
  assign('fin.vals', fin.vals, envir = .GlobalEnv)
  
  print(s)  
} #End of simulation loop

plot(c(fin.vals[ , , 2]), c(fin.vals[ , , 1]), pch = 16, cex = 0.6,
     xlab = expression(epsilon), ylab = expression(italic('P(e)')))

#Piece meal testing of components in function above ##########
sapply(c(1:stoch.sims), stoch.sim, init = start, k = kap.range[1], lam = lam.range[1], simplify = T)

plot(fill.test[ , 1, 1], fill.test[ , 7, 1], type='l', lwd = 2)

  for(i in 2:stoch.sims){
    lines(fill.test[ , 1, i], fill.test[ , 7, i], col = i, lwd = 2)
  }

#Get probability of elimination (P(e)) as number of chains that lead to extinction out of all chains
pe1 = as.numeric()

for(p in 1:stoch.sims){
  if(sum(fill.test[max(which(!is.na(fill.test[ , 7, p]))), , p][3:7]) == 0){ 
    pe1[p] = 1   #if no exposed, infected snails and no adult worms, elimination = 1
  } else {
    pe1[p] = 0   #else, elimination = 0
  }
  
}

sum(pe1) / stoch.sims


#Get W-pre and pos values, estimate BBR profile, and calculate epsilon
year.days = as.numeric()
for(i in 1:20){
  year.days[i] = 365*i + (i-1)
}

w.prepos = array(data = NA, dim = c(length(year.days), stoch.sims, 3))

for(w in 1:stoch.sims){
  w.prepos[, w, 1] = fill.test[ , 7, w][fill.test[ , 1, w] %in% year.days]     #w.pre values
  w.prepos[, w, 2] = fill.test[ , 7, w][fill.test[ , 1, w] %in% (year.days+1)] #w.pos values
}
for(b in 1:stoch.sims){
  w.prepos[, b, 3] = c((1/w.prepos[c(1:19), b, 2])*(w.prepos[c(2:20), b, 1] - w.prepos[c(1:19), b, 2]), 0)
}  

plot(c(1:19), w.prepos[c(1:19), 1, 3], pch = 16, ylim = c(-2, 2))
  for(r in 2:stoch.sims){
    points(c(1:19), w.prepos[c(1:19), r, 3], pch = 16, col = r)
  }

eps.k1l1 = as.numeric()

for(e in 1:stoch.sims){
  eps.k1l1[e] = lm(w.prepos[c(1:19), e, 3] ~ c(1:19))$coefficients[2]
}