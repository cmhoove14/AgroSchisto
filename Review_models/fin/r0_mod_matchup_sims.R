
library(foreach)
library(doParallel)
library(parallel)
library(fBasics)
library(ggplot2)

source('~/RemaisWork/Schisto/R Codes/ag_schist/Review_models/fin/r0_of_q_fin.R')
source('~/RemaisWork/Schisto/R Codes/ag_schist/Review_models/fin/mod_q_fin.R')

t.pars = expand.grid(lamda = seq(1e-6, 5e-5, length.out = 20), beta = seq(1e-5, 5e-4, length.out = 20))

t.pars$r0 = 0
t.pars$W = 0

tfx = function(b,l){
  r0 = r0.bl(beta = b, lambda = l)[3]
  
  ac.pars['beta'] = b
  ac.pars['lamda'] = l
  
  W = ode(ac.start, ac.time, agrochem_mod, ac.pars)[length(ac.time), 5]
  
  return(c(r0, W))
}

clus1 = makeCluster(detectCores()-1)
clusterExport(clus1, c('uniroot.all', 'kmat', 'tfx', 'r0.bl', 'parameters', 
                       'ac.pars', 'ode', 'agrochem_mod'))

for(i in 1:nrow(t.pars)){
  t.pars[i,3:4] = tfx(b = t.pars[i,2], l = t.pars[i,1])
  
  if(i %% 100 == 0) print(i)
}

#t.pars[,3:4]<-clusterMap(clus1, tfx, b = t.pars$beta, l = t.pars$lamda)

stopCluster(clus1)

ggplot(t.pars, aes(x = r0, y = W)) +
  theme_bw() +
  geom_line()

#Let's see if distribution of predation is the reason that r0 and mean worm burden aren't mathing up ###########

r0.bl2 = function(beta, lambda,i.frac, e.frac)
  
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
r0 = sqrt((theta*H*lamda*pi_C*omega*sigma*m*beta*N.eq*v*pi_M*omega) / 
            ((mu_N + mu_I + (i.frac*P.eq*psi.eq))*(mu_N + (e.frac*P.eq*psi.eq) + sigma)*2*(mu_H + mu_W)))

return(c(N.eq, P.eq, r0))

}

t.pars$r02 = 0
t.pars$W2 = 0

tfx2 = function(b,l){
  ac.pars['beta'] = b
  ac.pars['lamda'] = l
  
  sv = ode(ac.start, ac.time, agrochem_mod, ac.pars)[length(ac.time), ]
  W = sv[5]
  
  i.f = sv[4] / sum(sv[2:4])
  e.f = sv[3] / sum(sv[2:4])
  
  r0 = r0.bl2(beta = b, lambda = l, i.frac = i.f, e.frac = e.f)[3]
  
  return(c(r0, W))
}

clus1 = makeCluster(detectCores()-1)
clusterExport(clus1, c('uniroot.all', 'kmat', 'tfx', 'r0.bl', 'parameters', 
                       'ac.pars', 'ode', 'agrochem_mod'))

for(i in 1:nrow(t.pars)){
  t.pars[i,5:6] = tfx2(b = t.pars[i,2], l = t.pars[i,1])
  
  if(i %% 100 == 0) print(i)
}

stopCluster(clus1)

ggplot(t.pars, aes(x = r02, y = W2)) +
  theme_bw() +
  geom_line()

ggplot(t.pars, aes(x = r0, y = r02)) +
  theme_bw() +
  geom_line()

hist(t.pars$r02 - t.pars$r0, breaks = 30)

t.pars.reduce<-t.pars[which(t.pars$r02 < 5 & t.pars$r02 > 1 & t.pars$W2 < 80),]

ggplot(t.pars, aes(x = W, y = W2)) +
  theme_bw() +
  geom_line()

#Check area scaling ###########

area = 200

#Model parameters

parameters=c(
  # Location parameters
  A = area,          # Area of site of interest, m^2
  H = 1.5*area,      # Human population at site of interest (based on 300 people at 200m^2 water contact site from Sokolow PNAS)
  Om = 1 / sqrt(area),# degree of overlap between water contamination, snail, and human habitats
  
  # Snail reproductive parameters
  f_N = 0.10,        # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection - more like a recruitment rate); 
  #   from Woolhouse & Chandiwana et al. 1990 
  phi_N = 50*area,   # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
  z = 0.7,           # Fraction of exposed snails that reproduce
  
  # Snail mortality parameters
  mu_N = 1/90,        # Natural mortality rate of large snails (deaths/snail/day; mean lifespan ~ 90 days (~1 month))
  mu_I = 1/10-1/90,        # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
  
  # Predator pop dynamic parameters (from pred tweaking code)
  f_P = 0.117,         #Predator intrinsic recruitment rate assuming 30000 eggs per year and 1% survival from eggs to repro maturity
  #   from Adha-Ar et al 2016
  phi_P = 0.125*area,   #Predator CC (0.125/m^2 from Sokolow et al realized predator density, discounted here for natural pop)
  mu_P = 0.0038,       #Predator mortality rate from Halstead et al paper
  
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


t.pars$r0.a200 = 0
t.pars$W.a200 = 0
t.pars$r02.a200 = 0
t.pars$W2.a200 = 0

clus1 = makeCluster(detectCores()-1)
clusterExport(clus1, c('uniroot.all', 'kmat', 'tfx', 'r0.bl', 'parameters', 
                       'ac.pars', 'ode', 'agrochem_mod'))

for(i in 1:nrow(t.pars)){
  t.pars[i,7:8] = tfx(b = t.pars[i,2]/sqrt(area), l = t.pars[i,1]/sqrt(area))
  t.pars[i,9:10] = tfx2(b = t.pars[i,2]/sqrt(area), l = t.pars[i,1]/sqrt(area))
  
  if(i %% 100 == 0) print(i)
}

#t.pars[,3:4]<-clusterMap(clus1, tfx, b = t.pars$beta, l = t.pars$lamda)

stopCluster(clus1)