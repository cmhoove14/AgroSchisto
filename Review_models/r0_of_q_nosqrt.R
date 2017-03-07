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

#Model parameters ##############
area = 200
parameters=c(
  # Location parameters
  A = area,          # Area of site of interest, m^2
  H = 1.5*area,      # Human population at site of interest (based on 300 people at 200m^2 water contact site from Sokolow PNAS)
  Om = 1,            # degree of overlap between water contamination, snail, and human habitats
  
  # Snail reproductive parameters
  f_N = 0.60,        # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection - more like a recruitment rate); 
  #   from Woolhouse & Chandiwana et al. 1990 
  phi_N = 50*area,   # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
  z = 0.7,           # Fraction of exposed snails that reproduce
  
  # Snail mortality parameters
  mu_N = 0.03,        # Natural mortality rate of large snails (deaths/snail/day; mean lifespan ~ 33 days (~1 month))
  mu_I = 1/10,        # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
  
  # Predator pop dynamic parameters (from pred tweaking code)
  f_P = 0.02,          #Predator intrinsic recruitment rate
  phi_P = 0.1*area,    #Predator carrying capacity (0.1/m^2)
  mu_P = 0.0026,       #Predator mortality rate
  
  # Predation parameters
  alpha = 0.003,     # predator attack rate on snails
  Th = 1/14,         # Predator handling time
  nn = 1,            # exponent of the Holling's type III functional response
  
  # miracidia parameters
  m = 432,           # Miracidial shedding rate per adult female worm assuming 0.36 eggs/mL urine and 1200 mL urine per person per day
  v = 0.084,         # Egg viability of S. haematobium (i.e. miracidia/egg)
  pi_M = 1,          # Miracidial infectivity parameter
  
  # cercariae parameters
  theta = 109,       # cercarial shedding rate in snails (cercariae/I-snail/day); Pfluger 1984
  pi_C = 1,          # cercarial infectivity parameter
  
  # transmission parameters
  beta = 2e-4/area,  # Human-to-snail infection probability in reference area (infected snails/miracidia/snail/day)
  sigma = 1/40,      # Latent period for exposed snails (infectious snails/exposed snail/day))
  lamda = 2e-4/area, # Snail-to-human infection probability per cercaria
  k=0.2,             # Clumping parameter of negative binomial distribution of worms in humans
  
  # Schisto mortality parameters
  mu_W = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  mu_H = 1/(60*365)   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
)

#Null (default) functional responses ##############
  nil1<-function(In){ #For proportional responses
    return(1)
  }
  nil0<-function(In){ #For additive responses
    return(0)
  }

#R0(q) functions ###############
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
  pi_Mq = f.pi_Mq(In)
  pi_Cq = f.pi_Cq(In)
  v_q = v * f.v_q(In)

  #Equilibrium estimate of P given prawn predator parameters and q
  P.eq = phi_P*(1 - muPq/f_P)         
  
  if(P.eq<0){P.eq = 0}
  
  #Equilibrium estimate of N given snail parameters
  N.eq = max(uniroot.all(f = function(N){(f_Nq)*(1 - N/phi_Nq) - mu_Nq - (P.eq*alpha_q)/(1+alpha*Th*N)},
                         c(0, as.numeric(phi_Nq))))
  
  if(N.eq<0){N.eq = 0}
  
  #Equilibrium predation rate estimate
  psi.eq = alpha_q/(1+alpha*Th*(N.eq/area))
  
  #R_0 of q estimate
  e.hat = (theta_q*H*lamda*pi_Cq*omega) / (mu_Nq + mu_I + P.eq*psi.eq)
  i.hat = sigma / (mu_Nq + P.eq*psi.eq + sigma)
  w.hat = (m*beta*N.eq*v_q*pi_Mq*omega) / (mu_H + mu_W)
  
  r0 = e.hat * i.hat * w.hat
  
  return(c(N.eq, P.eq, r0))

}

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
  pi_Mq = f.pi_Mq(He)
  pi_Cq = f.pi_Cq(He)
  v_q = v * f.v_q(He)

#Equilibrium estimate of P given prawn predator parameters and q
  P.eq = phi_P*(1 - muPq/f_P)         

    if(P.eq<0){P.eq = 0}

#Equilibrium estimate of N given snail parameters
  N.eq = max(uniroot.all(f = function(N){(f_Nq)*(1 - N/phi_Nq) - mu_Nq - (P.eq*alpha_q)/(1+alpha*Th*N)},
                         c(0, as.numeric(phi_Nq))))

    if(N.eq<0){N.eq = 0}

#Equilibrium predation rate estimate
  psi.eq = alpha_q/(1+alpha*Th*(N.eq/area))

#R_0 of q estimate
  e.hat = (theta_q*H*lamda*pi_Cq*omega) / (mu_Nq + mu_I + P.eq*psi.eq)
  i.hat = sigma / (mu_Nq + P.eq*psi.eq + sigma)
  w.hat = (m*beta*N.eq*v_q*pi_Mq*omega) / (mu_H + mu_W)

    r0 = e.hat * i.hat * w.hat

return(c(N.eq, P.eq, r0))

}
  
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
  pi_Mq = f.pi_Mq(Fe)
  pi_Cq = f.pi_Cq(Fe)
  v_q = v * f.v_q(Fe)

#Equilibrium estimate of P given prawn predator parameters and q
  P.eq = phi_P*(1 - muPq/f_P)         

    if(P.eq<0){P.eq = 0}

#Equilibrium estimate of N given snail parameters
  N.eq = max(uniroot.all(f = function(N){(f_Nq)*(1 - N/phi_Nq) - mu_Nq - (P.eq*alpha_q)/(1+alpha*Th*N)},
                         c(0, as.numeric(phi_Nq))))

    if(N.eq<0){N.eq = 0}

#Equilibrium predation rate estimate
  psi.eq = alpha_q/(1+alpha*Th*(N.eq/area))

#R_0 of q estimate
  e.hat = (theta_q*H*lamda*pi_Cq*omega) / (mu_Nq + mu_I + P.eq*psi.eq)
  i.hat = sigma / (mu_Nq + P.eq*psi.eq + sigma)
  w.hat = (m*beta*N.eq*v_q*pi_Mq*omega) / (mu_H + mu_W)

  r0 = e.hat * i.hat * w.hat

return(c(N.eq, P.eq, r0))

}

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
pi_Mq = pi_Mqx
pi_Cq = pi_Cqx
v_q = v * v_qx

#Equilibrium estimate of P given prawn predator parameters and q
P.eq = phi_P*(1 - muPq/f_P)         

  if(P.eq<0){P.eq = 0}

#Equilibrium estimate of N given snail parameters
  N.eq = max(uniroot.all(f = function(N){(f_Nq)*(1 - N/phi_Nq) - mu_Nq - (P.eq*alpha_q)/(1+alpha*Th*N)},
                         c(0, as.numeric(phi_Nq))))

  if(N.eq<0){N.eq = 0}

#Equilibrium predation rate estimate
  psi.eq = alpha_q/(1+alpha*Th*(N.eq/area))

#R_0 of q estimate
  e.hat = (theta_q*H*lamda*pi_Cq*omega) / (mu_Nq + mu_I + P.eq*psi.eq)
  i.hat = sigma / (mu_Nq + P.eq*psi.eq + sigma)
  w.hat = (m*beta*N.eq*v_q*pi_Mq*omega) / (mu_H + mu_W)

r0 = e.hat * i.hat * w.hat

return(c(N.eq, P.eq, r0))

}

#Do some tests #############  
pred.dens = seq(0,1,0.01)
  r0.pd = as.numeric()

for(i in 1:length(pred.dens)){
  parameters['phi_P'] = pred.dens[i]*area
  r0.pd[i] = r0.In(In = 0)[3]
}
 
  plot(pred.dens, r0.pd, type = 'l', lwd = 2, xlab = 'pred density', ylab = 'R0')
  
  parameters['phi_P'] = 0.1*area
  
r0.In(In = 0)  #R0 = 4.01: coexistence of snail population, pred population, and disease @ pred density of 0.1/m^2