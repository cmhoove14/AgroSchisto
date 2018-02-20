#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
require(rootSolve)
require(deSolve)

#R0(q) code to investigate influence of agrochemicals on schisto transmission
kmat = matrix(NA,2,2)
  kmat[1,1] = 0
  kmat[2,2] = 0

area = 1
#Model parameters ##############

parameters=c(
  # Location parameters
  A = area,          # Area of site of interest, m^2
  H = 1.5*area,      # Human population at site of interest (based on 300 people at 200m^2 water contact site from Sokolow PNAS)
  Om = 1 / sqrt(area),# degree of overlap between water contamination, snail, and human habitats
  
  # Snail reproductive parameters
  f_N = 0.1,        # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection - more like a recruitment rate); 
                     #   from Woolhouse & Chandiwana et al. 1990 
  phi_N = 50*area,   # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
  z = 0.7,           # Fraction of exposed snails that reproduce
  
  # Snail mortality parameters
  mu_N = 0.017,        # Natural mortality rate of large snails (deaths/snail/day; mean lifespan ~ 90 days (~1 month))
  mu_I = 0.083,        # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
  
  # Predator pop dynamic parameters (from pred tweaking code)
  f_P = 0.117,         #Predator intrinsic recruitment rate assuming 30000 eggs per year and 1% survival from eggs to repro maturity
                      #   from Adha-Ar et al 2016
  phi_P = 0.125*area,   #Predator CC (0.125/m^2 from Sokolow et al realized predator density, discounted here for natural pop)
  mu_P = 0.038,       #Predator mortality rate from Halstead et al paper
  
  # Predation parameters
  alpha = 0.02, # Predator attack rate at high prawn/snail weight ratio per Sokolow 2014 Acta Tropica; reduced by larger area
  Th = 0.3,        # Predator handling time at high prawn/snail weight ratio per Sokolow 2014 Acta Tropica
  nn = 1,            # exponent of the Holling's type III functional response
  
  # miracidia parameters
  m = 432,           # Miracidial shedding rate per adult female worm assuming 0.36 eggs/mL urine and 1200 mL urine per person per day
  v = 0.084,         # Egg viability of S. haematobium (i.e. miracidia/egg) from Halstead et al
  pi_M = 6.22,       # mean Miracidia-hrs per day in agrochemical-free water
  
  # cercariae parameters
  theta = 109/7,       # mean daily cercarial shedding rate of patently infected snails @25C (cercariae/I-snail/day); Pfluger 1984
  pi_C = 14.21,        # mean cercariae-hrs in agrochemical-free water
  
  # transmission parameters
  beta = 1.340292e-05 / sqrt(area),       # Human-to-snail infection probability in reference area (exposed snails/susceptible snail/miracidi-hr/day); 
                        #     divided by 24 to account for hourly scale of miracidial survival; average of best-fit beta values
                        #     from fitting procedure in Halstead et al, weighted by likelihood
  sigma = 1/40,         # Latent period for exposed snails (infectious snails/exposed snail/day))
  lamda = 4.722579e-05 / sqrt(area), # Snail-to-human infection probability in reference area (adult worms/cercariae-hr); 
                        #     divided by 24 to account for hourly scale of miracidial survival; average of best-fit twa lambda values
                        #     from fitting procedure in Halstead et al, weighted by likelihood
  k=0.08,               # Clumping parameter of negative binomial distribution of worms in humans
  
  # Schisto mortality parameters
  mu_W = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  mu_H = 1/(60*365),   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
  
  #treatment parameters
  cov=0.43, #coverage of treatment across the population, Lampsar I = 100/1000 =0.1 %, Lampsar II = 129/300 = 43%
  eff=0.95 # efficiency of the drug
)

#Null (default) functional responses ##############
  nil1<-function(In){ #For proportional responses
    return(1)
  }
  nil0<-function(In){ #For additive responses
    return(0)
  }

#R0(q) functions ###############
#r0 of insecticide function
r0.In = function(In = 0,
                f.f_Nq = nil1, 
                f.mu_Pq = nil0,
                f.phi_Nq = nil1, 
                f.mu_Nq = nil0, 
                f.alpha_q = nil1,
                f.theta_q = nil1, 
                f.pi_Mq = nil1, 
                f.pi_Cq = nil1, 
                f.v_q = nil1)
  
{ area = parameters['A']
  sigma = parameters['sigma']
  lamda = parameters['lamda']
  omega = parameters['Om']
  theta = parameters['theta']
  pi_C = parameters['pi_C']
  beta = parameters['beta']
  H = parameters['H']
  m = parameters['m']
  v = parameters['v']
  pi_M = parameters['pi_M']
  mu_N = parameters['mu_N']
  mu_I = parameters['mu_I']
  mu_H = parameters['mu_H']
  mu_W = parameters['mu_W']
  mu_P = parameters['mu_P']
  f_P = parameters['f_P']
  phi_P = parameters['phi_P']
  alpha = parameters['alpha']
  Th = parameters['Th']
  f_N = parameters['f_N']
  phi_N = parameters['phi_N']
  
  f_Nq = f_N * f.f_Nq(In)
  muPq = mu_P + f.mu_Pq(In)
  phi_Nq = phi_N * f.phi_Nq(In)
  mu_Nq = mu_N + f.mu_Nq(In)
  alpha_q = alpha * f.alpha_q(In)
  theta_q = theta * f.theta_q(In)
  pi_Mq = pi_M * f.pi_Mq(In)
  pi_Cq = pi_C * f.pi_Cq(In)
  v_q = v * f.v_q(In)

  #Equilibrium estimate of P given prawn predator parameters and q
  P.eq = phi_P*(1 - muPq/f_P)         
  
  if(P.eq<0){P.eq = 0}
  
  #Equilibrium estimate of N given snail parameters
  N.eq = max(uniroot.all(f = function(N){(f_Nq)*(1 - N/phi_Nq) - mu_Nq - (P.eq*alpha_q*(1/area))/(1+alpha*Th*(N/area))},
                         c(0, as.numeric(phi_Nq))))
  
  if(N.eq<0){N.eq = 0}
  
  #Equilibrium predation rate estimate
  psi.eq = (alpha_q*(N.eq/area))/(1+alpha*Th*(N.eq/area))
  
#R_0 of q estimate
  kmat[1,2] = (beta*omega*N.eq*H*m*v_q*pi_Mq)/(2*mu_H+mu_W)
  kmat[2,1] = (sigma*lamda*omega*theta_q*pi_Cq)/((mu_Nq+P.eq*psi.eq+sigma)*(mu_Nq+mu_I+P.eq*psi.eq))
  
  r0 = max(eigen(kmat)$values)
  
  #r0 = sqrt((theta_q*H*lamda*pi_Cq*omega*sigma*m*beta*N.eq*v_q*pi_Mq*omega) / 
  #            (2*(mu_Nq + mu_I + P.eq*psi.eq)*(mu_Nq + P.eq*psi.eq + sigma)*(mu_H + mu_W)))
  
  
  return(c(N.eq, P.eq, r0))

}
#r0 of herbicide function
r0.He = function(He = 0,
                 f.f_Nq = nil1, 
                 f.mu_Pq = nil0,
                 f.phi_Nq = nil1, 
                 f.mu_Nq = nil0, 
                 f.alpha_q = nil1,
                 f.theta_q = nil1, 
                 f.pi_Mq = nil1, 
                 f.pi_Cq = nil1, 
                 f.v_q = nil1)
  
{ area = parameters['A']
  sigma = parameters['sigma']
  lamda = parameters['lamda']
  omega = parameters['Om']
  theta = parameters['theta']
  pi_C = parameters['pi_C']
  beta = parameters['beta']
  H = parameters['H']
  m = parameters['m']
  v = parameters['v']
  pi_M = parameters['pi_M']
  mu_N = parameters['mu_N']
  mu_I = parameters['mu_I']
  mu_H = parameters['mu_H']
  mu_W = parameters['mu_W']
  mu_P = parameters['mu_P']
  f_P = parameters['f_P']
  phi_P = parameters['phi_P']
  alpha = parameters['alpha']
  Th = parameters['Th']
  f_N = parameters['f_N']
  phi_N = parameters['phi_N']
  
  f_Nq = f_N * f.f_Nq(He)
  muPq = mu_P + f.mu_Pq(He)
  phi_Nq = phi_N * f.phi_Nq(He)
  mu_Nq = mu_N + f.mu_Nq(He)
  alpha_q = alpha * f.alpha_q(He)
  theta_q = theta * f.theta_q(He)
  pi_Mq = pi_M * f.pi_Mq(He)
  pi_Cq = pi_C * f.pi_Cq(He)
  v_q = v * f.v_q(He)

#Equilibrium estimate of P given prawn predator parameters and q
  P.eq = phi_P*(1 - muPq/f_P)         

    if(P.eq<0){P.eq = 0}

  #Equilibrium estimate of N given snail parameters
  N.eq = max(uniroot.all(f = function(N){(f_Nq)*(1 - N/phi_Nq) - mu_Nq - (P.eq*alpha_q*(1/area))/(1+alpha*Th*(N/area))},
                         c(0, as.numeric(phi_Nq))))
  
  if(N.eq<0){N.eq = 0}
  
  #Equilibrium predation rate estimate
  psi.eq = (alpha_q*(N.eq/area))/(1+alpha*Th*(N.eq/area))
  
#R_0 of q estimate
  kmat[1,2] = (beta*omega*N.eq*H*m*v_q*pi_Mq)/(2*mu_H+mu_W)
  kmat[2,1] = (sigma*lamda*omega*theta_q*pi_Cq)/((mu_Nq+P.eq*psi.eq+sigma)*(mu_Nq+mu_I+P.eq*psi.eq))
  
  r0 = max(eigen(kmat)$values)
  #r0 = sqrt((theta_q*H*lamda*pi_Cq*omega*sigma*m*beta*N.eq*v_q*pi_Mq*omega) / 
  #           (2*(mu_Nq + mu_I + P.eq*psi.eq)*(mu_Nq + P.eq*psi.eq + sigma)*(mu_H + mu_W)))

return(c(N.eq, P.eq, r0))

}
#r0 of fertilizer function  
r0.Fe = function(Fe = 0,
                 f.f_Nq = nil1, 
                 f.mu_Pq = nil0,
                 f.phi_Nq = nil1, 
                 f.mu_Nq = nil0, 
                 f.alpha_q = nil1,
                 f.theta_q = nil1, 
                 f.pi_Mq = nil1, 
                 f.pi_Cq = nil1, 
                 f.v_q = nil1)
  
{ area = parameters['A']
  sigma = parameters['sigma']
  lamda = parameters['lamda']
  omega = parameters['Om']
  theta = parameters['theta']
  pi_C = parameters['pi_C']
  beta = parameters['beta']
  H = parameters['H']
  m = parameters['m']
  v = parameters['v']
  pi_M = parameters['pi_M']
  mu_N = parameters['mu_N']
  mu_I = parameters['mu_I']
  mu_H = parameters['mu_H']
  mu_W = parameters['mu_W']
  mu_P = parameters['mu_P']
  f_P = parameters['f_P']
  phi_P = parameters['phi_P']
  alpha = parameters['alpha']
  Th = parameters['Th']
  f_N = parameters['f_N']
  phi_N = parameters['phi_N']
  
  f_Nq = f_N * f.f_Nq(Fe)
  muPq = mu_P + f.mu_Pq(Fe)
  phi_Nq = phi_N * f.phi_Nq(Fe)
  mu_Nq = mu_N + f.mu_Nq(Fe)
  alpha_q = alpha * f.alpha_q(Fe)
  theta_q = theta * f.theta_q(Fe)
  pi_Mq = pi_M * f.pi_Mq(Fe)
  pi_Cq = pi_C * f.pi_Cq(Fe)
  v_q = v * f.v_q(Fe)

#Equilibrium estimate of P given prawn predator parameters and q
  P.eq = phi_P*(1 - muPq/f_P)         

    if(P.eq<0){P.eq = 0}

  #Equilibrium estimate of N given snail parameters
  N.eq = max(uniroot.all(f = function(N){(f_Nq)*(1 - N/phi_Nq) - mu_Nq - (P.eq*alpha_q*(1/area))/(1+alpha*Th*(N/area))},
                         c(0, as.numeric(phi_Nq))))
  
  if(N.eq<0){N.eq = 0}
  
  #Equilibrium predation rate estimate
  psi.eq = (alpha_q*(N.eq/area))/(1+alpha*Th*(N.eq/area))
  
#R_0 of q estimate
  kmat[1,2] = (beta*omega*N.eq*H*m*v_q*pi_Mq)/(2*mu_H+mu_W)
  kmat[2,1] = (sigma*lamda*omega*theta_q*pi_Cq)/((mu_Nq+P.eq*psi.eq+sigma)*(mu_Nq+mu_I+P.eq*psi.eq))
  
  r0 = max(eigen(kmat)$values)
  #r0 = sqrt((theta_q*H*lamda*pi_Cq*omega*sigma*m*beta*N.eq*v_q*pi_Mq*omega) / 
  #            (2*(mu_Nq + mu_I + P.eq*psi.eq)*(mu_Nq + P.eq*psi.eq + sigma)*(mu_H + mu_W)))

return(c(N.eq, P.eq, r0))

}
#r0 with fixed parameter values
r0.fix = function(f_Nqx = 1, 
                  mu_Pqx = 0,
                  phi_Nqx = 1, 
                  mu_Nqx = 0, 
                  alpha_qx = 1,
                  theta_qx = 1, 
                  pi_Mqx = 1, 
                  pi_Cqx = 1, 
                  v_qx = 1)
  
{ area = parameters['A']
sigma = parameters['sigma']
lamda = parameters['lamda']
omega = parameters['Om']
theta = parameters['theta']
pi_C = parameters['pi_C']
beta = parameters['beta']
H = parameters['H']
m = parameters['m']
v = parameters['v']
pi_M = parameters['pi_M']
mu_N = parameters['mu_N']
mu_I = parameters['mu_I']
mu_H = parameters['mu_H']
mu_W = parameters['mu_W']
mu_P = parameters['mu_P']
f_P = parameters['f_P']
phi_P = parameters['phi_P']
alpha = parameters['alpha']
Th = parameters['Th']
f_N = parameters['f_N']
phi_N = parameters['phi_N']

f_Nq = f_N * f_Nqx
muPq = mu_P + mu_Pqx
phi_Nq = phi_N * phi_Nqx
mu_Nq = mu_N + mu_Nqx
alpha_q = alpha * alpha_qx
theta_q = theta * theta_qx
pi_Mq = pi_M * pi_Mqx
pi_Cq = pi_C * pi_Cqx
v_q = v * v_qx

#Equilibrium estimate of P given prawn predator parameters and q
P.eq = phi_P*(1 - muPq/f_P)         

  if(P.eq<0){P.eq = 0}

#Equilibrium estimate of N given snail parameters
  N.eq = max(uniroot.all(f = function(N){(f_Nq)*(1 - N/phi_Nq) - mu_Nq - (P.eq*alpha_q*(1/area))/(1+alpha*Th*(N/area))},
                         c(0, as.numeric(phi_Nq))))
  
  if(N.eq<0){N.eq = 0}

#Equilibrium predation rate estimate
  psi.eq = (alpha_q*(N.eq/area))/(1+alpha*Th*(N.eq/area))

#R_0 of q estimate
  kmat[1,2] = (beta*omega*N.eq*H*m*v_q*pi_Mq)/(2*mu_H+mu_W)
  kmat[2,1] = (sigma*lamda*omega*theta_q*pi_Cq)/((mu_Nq+P.eq*psi.eq+sigma)*(mu_Nq+mu_I+P.eq*psi.eq))
  
  r0 = max(eigen(kmat)$values)
  #r0 = sqrt((theta_q*H*lamda*pi_Cq*omega*sigma*m*beta*N.eq*v_q*pi_Mq*omega) / 
  #            (2*(mu_Nq + mu_I + 0.1*P.eq*psi.eq)*(mu_Nq + 0.9*P.eq*psi.eq + sigma)*(mu_H + mu_W)))
  
return(c(N.eq, P.eq, r0))

}

#r0 with variable beta & lambda
r0.bl = function(beta, lambda)
  
{ area = parameters['A']
sigma = parameters['sigma']
lamda = lambda
omega = parameters['Om']
theta = parameters['theta']
pi_C = parameters['pi_C']
beta = beta
H = parameters['H']
m = parameters['m']
v = parameters['v']
pi_M = parameters['pi_M']
mu_N = parameters['mu_N']
mu_I = parameters['mu_I']
mu_H = parameters['mu_H']
mu_W = parameters['mu_W']
mu_P = parameters['mu_P']
f_P = parameters['f_P']
phi_P = parameters['phi_P']
alpha = parameters['alpha']
Th = parameters['Th']
f_N = parameters['f_N']
phi_N = parameters['phi_N']

#Equilibrium estimate of P given prawn predator parameters and q
P.eq = phi_P*(1 - mu_P/f_P)         

if(P.eq<0){P.eq = 0}

#Equilibrium estimate of N given snail parameters
N.eq = max(uniroot.all(f = function(N){(f_N)*(1 - N/phi_N) - mu_N - (P.eq*alpha*(1/area))/(1+alpha*Th*(N/area))},
                       c(0, as.numeric(phi_N))))

if(N.eq<0){N.eq = 0}

#Equilibrium predation rate estimate
psi.eq = (alpha*(N.eq/area))/(1+alpha*Th*(N.eq/area))

#R_0 of q estimate
  kmat[1,2] = (beta*omega*N.eq*H*m*v*pi_M)/(2*mu_H+mu_W)
  kmat[2,1] = (sigma*lamda*omega*theta*pi_C)/((mu_N+P.eq*psi.eq+sigma)*(mu_N+mu_I+P.eq*psi.eq))

  r0 = max(eigen(kmat)$values)
#r0 = sqrt((theta*H*lamda*pi_C*omega*sigma*m*beta*N.eq*v*pi_M*omega) / 
#            (2*(mu_N + mu_I + P.eq*psi.eq)*(mu_N + P.eq*psi.eq + sigma)*(mu_H + mu_W)))

return(c(N.eq, P.eq, r0))

}
#r0 for agrochemical combinations
r0.Ag = function(In = 0,
                 He = 0,
                 Fe = 0,
                 
                 f.in.f_Nq = nil1, 
                 f.he.f_Nq = nil1, 
                 f.fe.f_Nq = nil1, 
                 
                 f.in.mu_Pq = nil0,
                 f.he.mu_Pq = nil0,
                 f.fe.mu_Pq = nil0,
                 
                 f.in.phi_Nq = nil1, 
                 f.he.phi_Nq = nil1, 
                 f.fe.phi_Nq = nil1, 
                 
                 f.in.mu_Nq = nil0, 
                 f.he.mu_Nq = nil0, 
                 f.fe.mu_Nq = nil0,
                 
                 f.in.alpha_q = nil1,
                 f.he.alpha_q = nil1,
                 f.fe.alpha_q = nil1,
                 
                 f.in.theta_q = nil1, 
                 f.he.theta_q = nil1, 
                 f.fe.theta_q = nil1, 
                 
                 f.in.pi_Mq = nil1, 
                 f.he.pi_Mq = nil1, 
                 f.fe.pi_Mq = nil1, 
                 
                 f.in.pi_Cq = nil1, 
                 f.he.pi_Cq = nil1, 
                 f.fe.pi_Cq = nil1, 
                 
                 f.in.v_q = nil1,
                 f.he.v_q = nil1,
                 f.fe.v_q = nil1)

{ area = parameters['A']
sigma = parameters['sigma']
lamda = parameters['lamda']
omega = parameters['Om']
theta = parameters['theta']
pi_C = parameters['pi_C']
beta = parameters['beta']
H = parameters['H']
m = parameters['m']
v = parameters['v']
pi_M = parameters['pi_M']
mu_N = parameters['mu_N']
mu_I = parameters['mu_I']
mu_H = parameters['mu_H']
mu_W = parameters['mu_W']
mu_P = parameters['mu_P']
f_P = parameters['f_P']
phi_P = parameters['phi_P']
alpha = parameters['alpha']
Th = parameters['Th']
f_N = parameters['f_N']
phi_N = parameters['phi_N']

f_Nq = f_N * (1-sum(1-c(f.in.f_Nq(In), f.he.f_Nq(He), f.fe.f_Nq(Fe))))
  if(f_Nq < 0) f_Nq = 0
muPq = mu_P + (f.in.mu_Pq(In) + f.he.mu_Pq(He) + f.fe.mu_Pq(Fe))
phi_Nq = phi_N * (1-sum(1-c(f.in.phi_Nq(In), f.he.phi_Nq(He), f.fe.phi_Nq(Fe))))
  if(phi_Nq < 0) phi_Nq = 0
mu_Nq = mu_N + (f.in.mu_Nq(In) + f.he.mu_Nq(He) + f.fe.mu_Nq(Fe))
alpha_q = alpha * (1-sum(1-c(f.in.alpha_q(In), f.he.alpha_q(He), f.fe.alpha_q(Fe))))
  if(alpha_q < 0) alpha_q = 0
theta_q = theta * (1-sum(1-c(f.in.theta_q(In), f.he.theta_q(He), f.fe.theta_q(Fe))))
  if(theta_q < 0) theta_q = 0
pi_Mq = 1-sum(1-c(f.in.pi_Mq(In), f.he.pi_Mq(He), f.fe.pi_Mq(Fe)))
  if(pi_Mq < 0) pi_Mq = 0
pi_Cq = 1-sum(1-c(f.in.pi_Cq(In), f.he.pi_Cq(He), f.fe.pi_Cq(Fe)))
  if(pi_Cq < 0) pi_Cq = 0
v_q = v * (1-sum(1-c(f.in.v_q(In), f.he.v_q(He), f.fe.v_q(Fe))))
  if(v_q < 0) v_q = 0
#Equilibrium estimate of P given prawn predator parameters and q
P.eq = phi_P*(1 - muPq/f_P)         

if(P.eq<0){P.eq = 0}

#Equilibrium estimate of N given snail parameters
N.eq = max(uniroot.all(f = function(N){(f_Nq)*(1 - N/phi_Nq) - mu_Nq - (P.eq*alpha_q*(1/area))/(1+alpha*Th*(N/area))},
                       c(0, as.numeric(phi_Nq))))

if(N.eq<0){N.eq = 0}

#Equilibrium predation rate estimate
psi.eq = (alpha_q*(N.eq/area))/(1+alpha*Th*(N.eq/area))

#R_0 of q estimate
  kmat[1,2] = (beta*omega*N.eq*H*m*v_q*pi_Mq)/(2*mu_H+mu_W)
  kmat[2,1] = (sigma*lamda*omega*theta_q*pi_Cq)/((mu_Nq+P.eq*psi.eq+sigma)*(mu_Nq+mu_I+P.eq*psi.eq))
  
  r0 = max(eigen(kmat)$values)
#r0 = sqrt((theta_q*H*lamda*pi_Cq*omega*sigma*m*beta*N.eq*v_q*pi_Mq*omega) / 
#            (2*(mu_Nq + mu_I + P.eq*psi.eq)*(mu_Nq + P.eq*psi.eq + sigma)*(mu_H + mu_W)))

return(c(N.eq, P.eq, r0))

}

#Do some tests #############  
#Influence of variable predator densities on r0 #########
pred.dens = seq(0,1,0.01)
  r0.pd = as.numeric()
  n.pd = as.numeric()

for(i in 1:length(pred.dens)){
  parameters['phi_P'] = pred.dens[i]*area
  r0.pd[i] = r0.In(In = 0)[3]
  n.pd[i] = r0.In(In = 0)[1]
}
 
  plot(pred.dens, r0.pd, type = 'l', lwd = 2, xlab = 'pred CC', ylab = 'R0')
  plot(pred.dens, n.pd, type = 'l', lwd = 2, col = 3, xlab = 'pred CC', ylab = 'N*')
  
  parameters['phi_P'] = 0.125*area
  
r0.In(In = 0)  #R0 = 2.00: coexistence of snail population, pred population, and disease @ realized pred density of ~0.1/m^2

#influence of various mortality rates on r0 ###############
mup.dens = seq(0,1,0.01) - parameters['mu_P']
r0.mup = as.numeric()

for(i in 1:length(mup.dens)){
  parameters['mu_P'] = mup.dens[i]+0.038
  r0.mup[i] = r0.In(In = 0)[3]
}

plot(mup.dens, r0.mup, type = 'l', lwd = 2, xlab = 'pred mortality rate', ylab = 'R0')

parameters['mu_P'] = 0.038 #Reset mortality rate to original

#influence of various snail carrying capacities on r0 ###############
phin.dens = seq(10, 200, 10)*area
r0.phin = as.numeric()

for(i in 1:length(phin.dens)){
  parameters['phi_N'] = phin.dens[i]
  r0.phin[i] = r0.In(In = 0)[3]
}

plot(phin.dens, r0.phin, type = 'l', lwd = 2, xlab = 'snail carrying capacity', ylab = 'R0')

parameters['phi_N'] = 50*area #Reset carrying capacity to original 

#influence of various snail reproductive rates on r0 ###############
fn.dens = seq(parameters['f_N'], 0, length.out = 50)
  r0.fn = as.numeric()

  for(i in 1:length(fn.dens)){
    parameters['f_N'] = fn.dens[i]
    r0.fn[i] = r0.In(In = 0)[3]
  }

plot(fn.dens, r0.fn, type = 'l', lwd = 2, xlab = 'snail reproductive rate', ylab = 'R0')

parameters['f_N'] = 0.1 #Reset reproductive rate to original

#influence of various snail mortality rates on r0 ###############
mun.dens = seq(0.0001, 1, length.out = 100)
r0.mun = as.numeric()

for(i in 1:length(mun.dens)){
  parameters['mu_N'] = mun.dens[i]
  r0.mun[i] = r0.In(In = 0)[3]
}

plot(mun.dens, r0.mun, type = 'l', lwd = 2, xlab = 'snail mortality rate', ylab = 'R0')

parameters['mu_N'] = 0.017 #Reset reproductive rate to original


#Influence of variable beta and lambda on r0 #################
beta.vec = seq(1e-9, 1e-5, length.out = 1000)
  r0.beta = as.numeric()
  
for(i in 1:length(beta.vec)){
  r0.beta[i] = r0.bl(beta = beta.vec[i], lambda = 1e-6)[3] 
}
  
  plot(beta.vec, r0.beta, type = 'l', lwd = 2, xlab = 'transmission parameter', ylab = 'R0')
  
lambda.vec = seq(1e-9, 1e-5, length.out = 1000)
  r0.lambda = as.numeric()
  
for(i in 1:length(lambda.vec)){
  r0.lambda[i] = r0.bl(lambda = lambda.vec[i], beta = 1e-6)[3] 
}
  
  lines(lambda.vec, r0.lambda, type = 'l', lwd = 2, col = 2, lty = 2)

# Re read Model parameters to make sure everything is set properly ##############
  parameters=c(
    # Location parameters
    A = area,          # Area of site of interest, m^2
    H = 1.5*area,      # Human population at site of interest (based on 300 people at 200m^2 water contact site from Sokolow PNAS)
    Om = 1 / sqrt(area),# degree of overlap between water contamination, snail, and human habitats
    
    # Snail reproductive parameters
    f_N = 0.1,        # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection - more like a recruitment rate); 
    #   from Woolhouse & Chandiwana et al. 1990 
    phi_N = 50*area,   # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
    z = 0.7,           # Fraction of exposed snails that reproduce
    
    # Snail mortality parameters
    mu_N = 0.017,        # Natural mortality rate of large snails (deaths/snail/day; mean lifespan ~ 90 days (~1 month))
    mu_I = 0.083,        # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
    
    # Predator pop dynamic parameters (from pred tweaking code)
    f_P = 0.117,         #Predator intrinsic recruitment rate assuming 30000 eggs per year and 1% survival from eggs to repro maturity
    #   from Adha-Ar et al 2016
    phi_P = 0.15*area,   #Predator CC (0.125/m^2 from Sokolow et al realized predator density, discounted here for natural pop)
    mu_P = 0.038,       #Predator mortality rate from Halstead et al paper
    
    # Predation parameters
    alpha = 0.003/sqrt(area), # Predator attack rate at high prawn/snail weight ratio per Sokolow 2014 Acta Tropica; reduced by larger area
    Th = 0.067,        # Predator handling time at high prawn/snail weight ratio per Sokolow 2014 Acta Tropica
    nn = 1,            # exponent of the Holling's type III functional response
    
    # miracidia parameters
    m = 432,           # Miracidial shedding rate per adult female worm assuming 0.36 eggs/mL urine and 1200 mL urine per person per day
    v = 0.084,         # Egg viability of S. haematobium (i.e. miracidia/egg) from Halstead et al
    pi_M = 6.22,       # mean Miracidia-hrs per day in agrochemical-free water
    
    # cercariae parameters
    theta = 109/7,       # mean daily cercarial shedding rate of patently infected snails @25C (cercariae/I-snail/day); Pfluger 1984
    pi_C = 14.21,        # mean cercariae-hrs in agrochemical-free water
    
    # transmission parameters
    beta = 1.340292e-05 / sqrt(area),       # Human-to-snail infection probability in reference area (exposed snails/susceptible snail/miracidi-hr/day); 
    #     divided by 24 to account for hourly scale of miracidial survival; average of best-fit beta values
    #     from fitting procedure in Halstead et al, weighted by likelihood
    sigma = 1/40,         # Latent period for exposed snails (infectious snails/exposed snail/day))
    lamda = 4.722579e-05 / sqrt(area), # Snail-to-human infection probability in reference area (adult worms/cercariae-hr); 
    #     divided by 24 to account for hourly scale of miracidial survival; average of best-fit twa lambda values
    #     from fitting procedure in Halstead et al, weighted by likelihood
    k=0.08,               # Clumping parameter of negative binomial distribution of worms in humans
    
    # Schisto mortality parameters
    mu_W = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
    mu_H = 1/(60*365),   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
    
    #treatment parameters
    cov=0.43, #coverage of treatment across the population, Lampsar I = 100/1000 =0.1 %, Lampsar II = 129/300 = 43%
    eff=0.95 # efficiency of the drug
  )  
  
  r0.In(In = 0)