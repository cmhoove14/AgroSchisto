#Load packages and other files ##########
library(adaptivetau)
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

#parameter values to sim over & objects to fill in simulation #########
par.sims = 5                                    #number of values within test parameter range to simulate
lam.range = seq(1e-4, 3e-4, length.out = par.sims)  #transmission intensity
kap.range = seq(0, 2, length.out = par.sims)            #Pos. density dependence

stoch.sims = 10   #number of simulations for each parameter set

fill = list()     #list to fill in each simulation

years = 61

fill.k1l1 = array(data = NA, dim = c(years*365*2, 7, stoch.sims))    #array to fill with all simulations

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

  fill.k1l1[c(1:nrow(matfin)), , sim] = cbind(matfin, Wm = cov*matfin[,5] + (1-cov)*matfin[,6])
  
  assign('fill.k1l1', fill.k1l1, envir = .GlobalEnv)
  
}

#test ######
stoch.sim(init = start, k = kap.range[1], lam = lam.range[1], sim = 1)

  sim1 = fill.k1l1[ , , 1]
  
plot(sim1[,1], sim1[,7], type = 'l', lwd = 2)

stoch.sim(init = start, k = kap.range[1], lam = lam.range[1], sim = 2)

  sim2 = fill.k1l1[ , , 2]

lines(sim2[,1], sim2[,7], lwd = 2, col=2)

#Run simulations #######
sapply(c(1:stoch.sims), stoch.sim, init = start, k = kap.range[1], lam = lam.range[1], simplify = T)

plot(fill.k1l1[ , 1, 1], fill.k1l1[ , 7, 1], type='l', lwd = 2)

  for(i in 2:stoch.sims){
    lines(fill.k1l1[ , 1, i], fill.k1l1[ , 7, i], col = i, lwd = 2)
  }

#Get W-pre and pos values to estimate BBR
year.days = as.numeric()
for(i in 1:20){
  year.days[i] = 365*i + (i-1)
}

w.prepos = array(data = NA, dim = c(length(year.days), stoch.sims, 3))

for(w in 1:stoch.sims){
  w.prepos[, w, 1] = fill.k1l1[ , 7, w][fill.k1l1[ , 1, w] %in% year.days]     #w.pre values
  w.prepos[, w, 2] = fill.k1l1[ , 7, w][fill.k1l1[ , 1, w] %in% (year.days+1)] #w.pos values
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