## library of schisto models used to study Reff and bounceback rate #####
#########################################################################

#Model structure and equations ####################
schisto_halstead_2pops_mda=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    P=n[6]
    N=S+E+I
    
    W=(cov*Wt) + ((1-cov)*Wu) #weighting treated and untreated populations
    #Miracidia production; function of adult female worms alive in the system (W*H*0.5) assuming 1:1 sex ratio,
    #mating probability function (gamma) from Anderson and May,
    
    fx<-function(x, mean.worm, clump){
      alpha<-(mean.worm)/(clump+mean.worm)
      (1-cos(x)) / ( (1+(alpha*cos(x)))^(1+clump) )
    }
    phi_Wk<-function(W, k){
      alpha<-W/(W+k)
      1-( (1-alpha)^(k+1) * (integrate(fx, 0, 2*pi, W, k)$value) /(2*pi)  )
    }
    #     a<-W/(W+k)
    #     gamma = 1- ( (1-a)^(k+1) *(integrate(fx, 0, 2*pi, W, k)$value) / 2*pi   )
    gamma= phi_Wk(W,k)
    
    #per-reproductive female egg production (m),
    #egg viability or the fraction of eggs that successfully hatch into viable miracidia (v),
    
    M=((0.5*W*H)*gamma)#*m*u_H*(v*vq)
    
    #miracidial mortality and infectivity (perhaps influenced by agrochemicals) affects beta
    
    pred= (alpha*P)/(1+(alpha*N*Th)) #death rate of snails due to predators (Prawns)
    
    dSdt= f_N*(1-(N/(phi_N*phi_Nq)))*(S+E) - 
      mu_N*S - pred*S - beta*M*S #Susceptible snails
    
    dEdt= beta*M*S - (mu_N+pred+sigma)*E #Exposed snails
    
    dIdt= sigma*E - (mu_N+pred+mu_I)*I #Infected snails
    
    #worm burden in human
    dWtdt= (lamda*I) - ((mu_W+mu_H)*Wt) - (eff*Wt*mda )
    dWudt= (lamda*I) - ((mu_W+mu_H)*Wu)
    
    
    dPdt= f_P*(1-(P/phi_P))*P-(mu_P+muPq)*P #prawn population
    
    
    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt, dPdt)))
  }) 
} 

#Model structure and equations ####################
schisto_halstead_2pops_mda_noPhi=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    P=n[6]
    N=S+E+I
    
    W=(cov*Wt) + ((1-cov)*Wu) #weighting treated and untreated populations
    #Miracidia production; function of adult female worms alive in the system (W*H*0.5) assuming 1:1 sex ratio,
    #mating probability function (gamma) from Anderson and May,
    
    fx<-function(x, mean.worm, clump){
      alpha<-(mean.worm)/(clump+mean.worm)
      (1-cos(x)) / ( (1+(alpha*cos(x)))^(1+clump) )
    }
    phi_Wk<-function(W, k){
      alpha<-W/(W+k)
      1-( (1-alpha)^(k+1) * (integrate(fx, 0, 2*pi, W, k)$value) /(2*pi)  )
    }
    #     a<-W/(W+k)
    #     gamma = 1- ( (1-a)^(k+1) *(integrate(fx, 0, 2*pi, W, k)$value) / 2*pi   )
    gamma= phi_Wk(W,k)
    
    #per-reproductive female egg production (m),
    #egg viability or the fraction of eggs that successfully hatch into viable miracidia (v),
    
    M=(0.5*W*H) #*gamma)#*m*u_H*(v*vq)
    
    #miracidial mortality and infectivity (perhaps influenced by agrochemicals) affects beta
    
    pred= (alpha*P)/(1+(alpha*N*Th)) #death rate of snails due to predators (Prawns)
    
    dSdt= f_N*(1-(N/(phi_N*phi_Nq)))*(S+E) - 
      mu_N*S - pred*S - beta*M*S #Susceptible snails
    
    dEdt= beta*M*S - (mu_N+pred+sigma)*E #Exposed snails
    
    dIdt= sigma*E - (mu_N+pred+mu_I)*I #Infected snails
    
    #worm burden in human
    dWtdt= (lamda*I) - ((mu_W+mu_H)*Wt) - (eff*Wt*mda )
    dWudt= (lamda*I) - ((mu_W+mu_H)*Wu)
    
    
    dPdt= f_P*(1-(P/phi_P))*P-(mu_P+muPq)*P #prawn population
    
    
    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt, dPdt)))
  }) 
} 


#List parameters and values #####################
parameters_2pops_mda=c(
  ##standard snail parameters 
  f_N=0.1, # recruitment rate: snails/snail/day
  phi_N=10000, # carrying capacity: max snail population (corresponds to ~50/m^2)
  z=0.5, #Proportion of exposed snails that reproduce: density dependent, but assumed constant here
  mu_N=1/60, #Mortality rate from Anderson and May (from Chu 1966): deaths/snail/day
  sigma=1/40, #Transition rate from exposed to infected; ~latent period
  mu_I=1/10 - 1/60, #additional snail death due to infection
  ## snail parameters impacted by agrochemicals
  f_Nq=1, #Not affected in mesocosm
  phi_Nq=1, #Scalar of snail carrying capacity by chemical concentration INFORMED BY BOTTOM UP EFFECTS IN MESOCOSM
  #mu_Nq=0, #Chem concentrations too low to affect snails in mesocosom experiments
  
  
  #prawn parameters
  alpha=0.003, #attack rate
  Th=0.067,#~Prawn predation limit for procambarius clarkii: ~15 snails/day
  f_P=0.117,#prawn birth rate: prawns/prawn/day
  phi_P=120,  #prawn carrying capacity ~0.6/m^2
  mu_P= 1, # kill all prawns. 0.038095238, #observed daily 24 hr prawn mortality rate in chlorP-free tanks of mesocosm#0.00137, previous rate
  
  #prawn parameter impacted by agrochemicals
  muPq=0, #agrochemical dependent prawn mortality rate (added to baseline mortality rate) INFORMED BY MESOCOSM AND HALSTEAD 2015
  #alpha_q=1, #Scalar of predation rate due to sub-lethal toxicity; not considered in mesocosm
  
  #Worm parameters
  mu_W=1/(3.3*365), # death rate of adult worms
  m=.36, #eggs/mL urine/female adult worm
  v=0.084, #Egg viability controlling schistosome egg->infective miracidia
  vq=1, #agrochemical-caused reduction in egg viability
  
  #Human parameters
  H=300, #number of humans
  mu_H=1/(60*365), #Assumes ~60 year lifespan
  k=0.08, #clumping parameter of the negative binomial distribution
  u_H=1200, #mL urine per human/day (approximate, ranges from 800 - 2000)
  
  #Transmission parameters
  lamda=1.9e-04, #2.5e-5, #snail-to-man transmission: p(infected snail sheds cercariae that infects human and reaches adulthood)
  beta=1.63e-06, #1.75e-5, #man-to-snail transmission: p(mated female worm produces a miracidia that infects a snail)
  
  
  #treatment parameters
  cov=0.43, #coverage of treatment across the population, Lampsar I = 100/1000 =0.1 %, Lampsar II = 129/300 = 43%
  eff=0.95, # efficiency of the drug
  mda=0 # flag to indicate if mda is applied or not
)

fx<-function(x, mean.worm, clump){
  alpha<-(mean.worm)/(clump+mean.worm)
  (1-cos(x)) / ( (1+(alpha*cos(x)))^(1+clump) )
}
phi_Wk<-function(W, k){
  alpha<-W/(W+k)
  1-( (1-alpha)^(k+1) * (integrate(fx, 0, 2*pi, W, k)$value) /(2*pi)  )
}


parameters_2pops_mda_Chris1=c( #Excluding beta, Phi_Nq, f_Nq, and muPq which will be read into the R0 function
  ##standard snail parameters 
  f_N=0.10, # recruitment rate (from sokolow et al)
  phi_N=10000, # carrying capacity from sokolow et al
  z=0.5, #Proportion of exposed snails that reproduce from sokolow et al
  mu_N=1/60, #Mortality rate from Sokolow et al
  sigma=1/40, #Transition rate from exposed to infected from sokolow et al
  mu_I=1/10 - 1/60, #additional snail death due to infection from sokolow et al
  ## snail parameters impacted by agrochemicals
  f_Nq=1, #Not affected in mesocosm
  phi_Nq=1, #Scalar of snail carrying capacity by chemical concentration INFORMED BY BOTTOM UP EFFECTS IN MESOCOSM
  #mu_Nq=0, #Chem concentrations too low to affect snails in mesocosom experiments
  
  #prawn parameters
  alpha=0.003, #attack rate
  Th=0.067,#~Prawn predation limit
  f_P=0.234/2, #prawn birth rate from Cervantes-Santiago Aquaculture 2010 paper (/2 for 1:1 female-male ratio)
  phi_P=120,  #prawn carrying capacity
  mu_P= 0.03883984, #observed daily 24 hr prawn mortality rate in chlorP-free tanks of mesocosm#0.00137, previous rate
  muPq = 0,
  
  #Adult Worm, Miracidia and Circariae Parameters
  #lamda=1.5e-5, #probability of snail shedding a cercariae that infects a human host and survives to reproduction
  mu_W=1/(3.3*365), # death rate of adult worms
  m=0.5, #miracidial shedding rate per reproductive female divided by miracidial mortality; from sokolow et al
  
  #Human parameters
  H=300, #number of humans
  mu_H=1/(60*365), #Assumes 60 year lifespan
  k=0.08, #clumping parameter of the negative binomial distribution
  u_H=1200, #mL urine per human/day (approximate, ranges from 800 - 2000)
  
  #Transmission parameters
  lamda=1.5e-5, #snail-to-man transmission: p(infected snail sheds cercariae that infects human and reaches adulthood)
  beta=2.5e-5, #man-to-snail transmission: p(mated female worm produces a miracidia that infects a snail)
  
  #Treatment parameters
  eff=0.99, #Efficacy of PZQ treatment (fraction of worms cleared), as per Sanna's PNAS simulation
  cov=0.43, #coverage of PZQ treatment in the human population 129/300 in Lampsar II
  mda=0 #binary to indicate if MDA is applied
)


