## library of schisto models used to study Reff and bounceback rate #####
#########################################################################

#mating function that implements PDD#######
phi_Wk <- function(W, k) {
  if( k<= 0){
    return(1)
  } else {
    func <- function(x) {
    a <- ( W / (W + k) )
    b <- ((1-a)^(1+k))/(2*pi)
    return(( b*( 1-cos(x) ) / (( 1 + a*cos(x) )^(1+k)) ))
  }
  val = integrate(func, 0, 2*pi, subdivisions = 10000,
                  rel.tol = 1e-10, stop.on.error = FALSE)$value
  return(1-val)
  }
  
}
#Model with PDD ####################
schisto_mod_pdd=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    N=S+E+I
    
    W=(cov*Wt) + ((1-cov)*Wu) #weighting treated and untreated populations
    #Miracidia production; function of adult female worms alive in the system (W*H*0.5) assuming 1:1 sex ratio,
    #mating probability function (gamma) from Anderson and May,
    
    phi = phi_Wk(W = W, k = k)  #Mating probability

    #per-reproductive female egg production (m),
    #egg viability or the fraction of eggs that successfully hatch into viable miracidia (v),
    
    M=((0.5*W*H)*phi)#*m*u_H*(v*vq)
    

    dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S - beta*M*S #Susceptible snails
    
    dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails
    
    dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails
    
    #worm burden in human
    dWtdt= (lamda*I) - ((mu_W+mu_H)*Wt)
    dWudt= (lamda*I) - ((mu_W+mu_H)*Wu)
    
    
    
    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt)))
  }) 
} 

#Model with no pdd ####################
schisto_mod_nopdd=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    N=S+E+I
    
    W=(cov*Wt) + ((1-cov)*Wu) #weighting treated and untreated populations
    #Miracidia production; function of adult female worms alive in the system (W*H*0.5) assuming 1:1 sex ratio,
    #mating probability function (gamma) from Anderson and May,
    
    #per-reproductive female egg production (m),
    #egg viability or the fraction of eggs that successfully hatch into viable miracidia (v),
    
    M=(0.5*W*H) 
    
    #miracidial mortality and infectivity (perhaps influenced by agrochemicals) affects beta
    
    dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S  - beta*M*S #Susceptible snails
    
    dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails
    
    dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails
    
    #worm burden in human
    dWtdt= (lamda*I) - ((mu_W+mu_H)*Wt)
    dWudt= (lamda*I) - ((mu_W+mu_H)*Wu)
    
    
    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt)))
  }) 
} 

#Parameters #####################
pars_Chris1=c( #Excluding beta, Phi_Nq, f_Nq, and muPq which will be read into the R0 function
  ##standard snail parameters 
  f_N=0.10, # recruitment rate (from sokolow et al)
  C=10000, # carrying capacity from sokolow et al
  z=0.5, #Proportion of exposed snails that reproduce from sokolow et al
  mu_N=1/60, #Mortality rate from Sokolow et al
  sigma=1/40, #Transition rate from exposed to infected from sokolow et al
  mu_I=1/10 - 1/60, #additional snail death due to infection from sokolow et al
  
  #Adult Worm, Miracidia and Circariae Parameters
  #lamda=1.5e-5, #probability of snail shedding a cercariae that infects a human host and survives to reproduction
  mu_W=1/(3.3*365), # death rate of adult worms

  #Human parameters
  H=300, #number of humans
  mu_H=1/(60*365), #Assumes 60 year lifespan
  k=0.08, #clumping parameter of the negative binomial distribution
  u_H=1200, #mL urine per human/day (approximate, ranges from 800 - 2000)
  
  #Transmission parameters
  lamda=1.2e-4, #1.5e-5 #snail-to-man transmission: p(infected snail sheds cercariae that infects human and reaches adulthood)
  beta=1.6e-6,  #2.5e-5 #man-to-snail transmission: p(mated female worm produces a miracidia that infects a snail)
  
  cov = 0.8
)




#Model with PDD and additional NDDs: crowding of adult worms and adaptive immunity ##########
##Parasite Fecundity reduced due to crowding at high densities (NDD)
  f_Wgk <- function(W, gamma, k) {
    if(k <= 0){
      return(1)
    } else {
      return((1 + ((W*(1-(exp(-gamma))))/k))^(-k-1))
    }
    
  }
  
  pars_Chris1["gam"] <- 0.001

#Host acquired immunity which reduces transmission at high worm burdens    
  R_Wv <- function(W,v){
    exp(1-v*W-exp(-v*W))
  }
  
  pars_Chris1["v"] <- 0.0028

#Model with PDD and added NDDs  
schisto_mod_pdd_add_ndds=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    N=S+E+I
    
    W=(cov*Wt) + ((1-cov)*Wu) #weighting treated and untreated populations
    #Miracidia production; function of adult female worms alive in the system (W*H*0.5) assuming 1:1 sex ratio,
  #DD functions
    
    f = f_Wgk(W, gam, k)  
    R = R_Wv(W, v)
    phi = phi_Wk(W = W, k = k)  #Mating probability

    M=((0.5*W*H)*phi*f)#*m*u_H*(v*vq)
    

    dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S - beta*M*S #Susceptible snails
    
    dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails
    
    dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails
    
    #worm burden in human
    dWtdt= (lamda*I*R) - ((mu_W+mu_H)*Wt)
    dWudt= (lamda*I*R) - ((mu_W+mu_H)*Wu)
    
    
    
    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt)))
  }) 
} 

#Model with no PDD and added NDDs
schisto_mod_nopdd_add_ndds=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    N=S+E+I
    
    W=(cov*Wt) + ((1-cov)*Wu) #weighting treated and untreated populations
    #Miracidia production; function of adult female worms alive in the system (W*H*0.5) assuming 1:1 sex ratio,
    #mating probability function (gamma) from Anderson and May,
    
    f = f_Wgk(W, gam, k)  
    R = R_Wv(W, v)
    # phi = phi_Wk(W = W, k = k)  #Mating probability

    #per-reproductive female egg production (m),
    #egg viability or the fraction of eggs that successfully hatch into viable miracidia (v),
    
    M=((0.5*W*H)*f)#*m*u_H*(v*vq)
    

    dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S - beta*M*S #Susceptible snails
    
    dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails
    
    dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails
    
    #worm burden in human
    dWtdt= (lamda*I*R) - ((mu_W+mu_H)*Wt)
    dWudt= (lamda*I*R) - ((mu_W+mu_H)*Wu)
    
    
    
    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt)))
  }) 
} 
