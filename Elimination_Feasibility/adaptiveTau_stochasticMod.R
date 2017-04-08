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
source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/lib_schistoModels_DDandNODD_CH.R")
source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/Reff_fn_lib_CH.R")

outputfile<-"C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/shortlist_51.csv"
shortlist_first100<-as.data.frame(read.csv(file=outputfile, header=TRUE, sep=","))
sort_ind<-order(shortlist_first100$R0, decreasing=FALSE)
shortlist_first100<-shortlist_first100[sort_ind,]
  shortlist_first100 = shortlist_first100[c(1:100),] #Get rid of other parameter estimates

cov = 0.8
eff = 0.99
mda.years = c(2:21)
  
params = as.list(parameters_2pops_mda_Chris1)
  params$beta = shortlist_first100$beta[1]
  params$lamda = shortlist_first100$lamda.twa[1]
  params$cov = cov
  

init1 = c(S = 5000, # susceptible humans
          E = 2000, # infected humans
          I = 500, # infected humans
          Wt = 72,
          Wu = 72) # recovered (and immune) humans

eq = as.data.frame(ode(init1,seq(0,200*365,30),
      schisto_halstead_2pops_mda,params))[length(seq(0,200*365,30)), c(2:6)]

init1 = setNames(as.numeric(round(eq)), colnames(eq))

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

fill = list()

r1=as.data.frame(ssa.adaptivetau(init1, transitions, 
                                sfx, params, tf=365))

fill[[1]] = ssa.adaptivetau(init1, transitions, 
                               sfx, params, tf=365)

for(m in 2:max(mda.years)){
  init = setNames(as.numeric(fill[[m-1]][dim(fill[[m-1]])[1],c(2:6)]), 
                  colnames(fill[[m-1]])[c(2:6)]) #reset initial states
  
  init[4] = round(init[4]* (1-eff))  #apply MDA
  
  fill[[m]] = ssa.adaptivetau(init, transitions, 
                              sfx, params, tf=365) #stochastic sim for a year
  
  fill[[m]][,1] = fill[[m]][,1] + (365*(m-1)+1)    #adjust time
}

matfin = do.call(rbind,fill)

Wm = cov*matfin[,5] + (1-cov)*matfin[,6]

matfin = cbind(matfin, Wm)