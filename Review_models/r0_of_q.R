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

#R0(q) code to investigate influence of chlorpyrifos

#Model parameters ##############
area = 200
parameters=c(
  # Location parameters
  A = area,          # Area of site of interest, m^2
  H = 1.5*area,      # Human population at site of interest (based on 300 people at 200m^2 water contact site from Sokolow PNAS)
  Om = 1,            # degree of overlap between water contamination, snail, and human habitats
  
  # Snail reproductive parameters
  f_N = 0.10,        # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection - more like a recruitment rate); 
  #   from Sokolow et al. 2015 
  phi_N = 50*area,        # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
  z = 0.7,           # Fraction of exposed snails that reproduce
  
  # Snail mortality parameters
  mu_N = 1/60,        # Natural mortality rate of large snails (deaths/snail/day; assume mean lifespan = 60 days)
  mu_I = 1/10,        # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
  
  # Predator pop dynamic parameters (from pred tweaking code)
  f_P = 0.02,          #Predator intrinsic recruitment rate
  phi_P = 0.3*area,         #Predator carrying capacity (0.5/m^2)
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
  beta = 2e-4/area,       # Human-to-snail infection probability in reference area (infected snails/miracidia/snail/day)
  sigma = 1/40,      # Latent period for exposed snails (infectious snails/exposed snail/day))
  lamda = 2e-4/area,      # Snail-to-human infection probability per cercaria
  k=0.2,             # Clumping parameter of negative binomial distribution of worms in humans
  
  # Schisto mortality parameters
  mu_W = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  mu_H = 1/(60*365)   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
)

#Chlorpyrifos functional responses ##############
  source('T_variant_Agro/ChlorP_Tox.R')
  nil1<-function(In){
    return(1)
  }
  nil0<-function(In){
    return(0)
  }

#R0(q) function ###############
r0.In = function(conc = 0, 
                 f.f_Nq = nil1, 
                 f.mu_Pq = nil0,
                 f.phi_Nq = nil1, 
                 f.mu_Nq = nil0, 
                 f.alpha_q = nil1,
                 f.theta_q = nil1, 
                 f.pi_Mq = nil1, 
                 f.pi_Cq = nil1, 
                 f.v_q = nil1)
{In = conc
  
  area = parameters['A']
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
  N.eq = max(uniroot.all(f = function(y){
    (f_Nq)*(1 - (y)/phi_Nq) - mu_Nq - (P.eq*alpha_q)/(1+alpha*Th*(y))
    }, c(0, as.numeric(phi_Nq))))
  
  if(N.eq<0){N.eq = 0}
  
  #Equilibrium predation rate estimate
  psi.eq = alpha_q/(1+alpha*Th*N.eq)
  
  #R_0 of q estimate
  r0 <- sqrt((sigma*lamda*omega^2*theta_q*pi_Cq*beta*N.eq*H*m*v_q*pi_Mq) / 
              (2*(mu_Nq + P.eq*psi.eq + sigma)*(mu_Nq + P.eq*psi.eq + mu_I)*(mu_H + mu_W)))
  
  return(c(N.eq, P.eq, r0))

}
  
pred.dens = seq(0,1,0.01)
  r0.pd = as.numeric()

for(i in 1:length(pred.dens)){
  parameters['phi_P'] = pred.dens[i]*area
  r0.pd[i] = r0.In(conc = 0)[3]
  print(r0.In(conc = 0))
}
 
  plot(pred.dens, r0.pd, type = 'l', lwd = 2, xlab = 'pred density', ylab = 'R0')
  
parameters['phi_P'] = 0.1*area  
  
r0.In(conc = 0)  #R0 = 3.38: coexistence of snail population, pred population, and disease @ pred density of 0.1/m^2

r0.In(conc = 5, f.mu_Pq = muPq_ch_sat09) #R0 = 5.39: Pred population eliminated, max r0 in pred-free environment

r0.conc = matrix(ncol = 101, nrow = 12)

for(i in 1:100){
  r0.conc[1,i] = r0.In(conc = (i-1))[3]                       #Baseline r0 with no agrochemical effects included
  r0.conc[2,i] = r0.In(conc = (i-1),
                       f.f_Nq = f_Nq_chlor_ibrahim92)[3]      #Toxic effects on snail reproduction
  r0.conc[3,i] = r0.In(conc = (i-1),
                       f.mu_Pq = muPq_ch_sat09)[3]            #Mortality effects on predator population from Sat'09
  r0.conc[4,i] = r0.In(conc = (i-1),
                       f.alpha_q = alpha_q_ch_sat09)[3]       #Feeding reduction effects on predator population
  r0.conc[5,i] = r0.In(conc = (i-1),
                       f.mu_Nq = muNq_chlor_ibrahim92)[3]     #Mortality effects on snail population
  r0.conc[6,i] = r0.In(conc = (i-1),
                       f.mu_Nq = muNq_chlor_ibrahim92,
                       f.f_Nq = f_Nq_chlor_ibrahim92)[3]      #Both repro and mortality effects in snails
  r0.conc[7,i] = r0.In(conc = (i-1),
                       f.mu_Pq = muPq_ch_sat09,
                       f.alpha_q = alpha_q_ch_sat09)[3]       #Both attack rate and mortality effects in preds
  r0.conc[8,i] = r0.In(conc = (i-1),
                       f.mu_Pq = muPq_ch_Halstead)[3]         #Mortality effects on predator population from Halstead'15
  r0.conc[9,i] = r0.In(conc = (i-1),
                       f.pi_Mq = pi_Mq_ch_has11)[3]           #Mortality effects on miracidia from Hasheesh11
  r0.conc[10,i] = r0.In(conc = (i-1),
                       f.pi_Cq = pi_Cq_ch_has11)[3]           #Mortality effects on cercariae from Hasheesh11
  r0.conc[11,i] = r0.In(conc = (i-1),
                        f.pi_Mq = pi_Mq_ch_has11,
                        f.pi_Cq = pi_Cq_ch_has11)[3]          #Mortality effects on miracidia & cercariae from Hasheesh11
  r0.conc[12,i] = r0.In(conc = (i-1),
                        f.f_Nq = f_Nq_chlor_ibrahim92,
                        f.mu_Pq = muPq_ch_Halstead,
                        f.mu_Nq = muNq_chlor_ibrahim92,
                        f.pi_Mq = pi_Mq_ch_has11,
                        f.pi_Cq = pi_Cq_ch_has11)[3]          #Combined effects
}

plot(x=c(0:100), y=r0.conc[1,], type = 'l', lwd = 3, ylim = c(0,6), col = 'grey40',
     ylab = expression('R'[0]), xlab = 'Chlorpyrifos Concentration (ppb)',
     main = expression(paste('R'[0],'(q) simulations with Chlorpyrifos', sep='')))
  legend('bottomleft', legend = expression('baseline R'[0]), lwd=3, col = 'grey40', cex=0.75)
 #Snail effects
  lines(x=c(0:100), y=r0.conc[2,], lwd=2, col='lightblue')
  lines(x=c(0:100), y=r0.conc[5,], lwd=2, col='navy')
  lines(x=c(0:100), y=r0.conc[6,], lwd=2, col='blue', lty=2)
    legend('bottomleft', legend = c(expression('baseline R'[0]),
                                    expression(paste('f'[N],'(q)', sep='')),
                                    expression(paste(mu[N],'(q)', sep='')),
                                    'combined'), 
           lwd=c(3,2,2,2), lty = c(1,1,1,2), 
           col = c('grey40', 'lightblue', 'navy', 'blue'), cex=0.75)
  
 #Predator effects 
  lines(x=c(0:100), y=r0.conc[3,], lwd=2, col='olivedrab1')
  lines(x=c(0:100), y=r0.conc[4,], lwd=2, col='olivedrab4')
  lines(x=c(0:100), y=r0.conc[7,], lwd=2, col='palegreen', lty=2)
  lines(x=c(0:100), y=r0.conc[8,], lwd=2, col='green')
    legend('bottomleft', legend = c(expression('baseline R'[0]),
                                    expression(paste(mu[P],'(q) Halstead', sep='')),
                                    expression(paste(alpha,'(q) Sat.', sep='')),
                                    expression(paste(mu[P],'(q) Sat.', sep='')),
                                    'combined'), 
           lwd=c(3,2,2,2,2), lty = c(1,1,1,1,2), 
           col = c('grey40', 'green', 'olivedrab4', 'olivedrab1', 'palegreen'), cex=0.75)
 #Larval effects
  lines(x=c(0:100), y=r0.conc[9,], lwd=2, col='indianred4')
  lines(x=c(0:100), y=r0.conc[10,], lwd=2, col='lightcoral')
  lines(x=c(0:100), y=r0.conc[11,], lwd=2, col='red', lty=2)
    legend('bottomleft', legend = c(expression('baseline R'[0]),
                                    expression(paste(pi[M],'(q)', sep='')),
                                    expression(paste(pi[C],'(q)', sep='')),
                                    'combined'), 
           lwd=c(3,2,2,2), lty = c(1,1,1,2), 
           col = c('grey40', 'indianred4', 'lightcoral', 'red'), cex=0.75)
 #All effects combined
  lines(x=c(0:100), y=r0.conc[12,], lwd=3, col='black', lty=3)
    legend('bottomleft', legend = c(expression('baseline R'[0]),
                                    expression(paste('R'[0],'(q)', sep=''))), 
           lwd=c(3,3), lty = c(1,3), col = c('grey40', 'black'), cex=0.75)
    