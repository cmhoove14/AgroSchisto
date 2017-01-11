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

#Model parameters ##############
area = 1
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
  phi_P = 0.5*area,         #Predator carrying capacity (0.5/m^2)
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
  beta = 1e-5,       # Human-to-snail infection probability in reference area (infected snails/miracidia/snail/day)
  sigma = 1/40,      # Latent period for exposed snails (infectious snails/exposed snail/day))
  lamda = 2.8e-4,      # Snail-to-human infection probability per cercaria
  k=0.2,             # Clumping parameter of negative binomial distribution of worms in humans
  
  # Schisto mortality parameters
  mu_W = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  mu_H = 1/(60*365)   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
)

require(sensitivity)

#Get parameter sets to sample from #######################
sims = 50
#parameters ranges
f_N.range<-seq(0.05, 0.45, length.out = sims)
phi_N.range<-seq(10, 70, length.out = sims)
z.range<-seq(0.2, 1.0, length.out = sims)
mu_N.range<-seq(1/20, 1/100, length.out = sims)
mu_I.range<-seq(1/2, 1/20, length.out = sims)
f_P.range<-seq(0.001, 0.2, length.out = sims)
phi_P.range<-seq(0.1, 5.0, length.out = sims)
mu_P.range<-seq(1/100, 1/(365*2), length.out = sims)
alpha.range<-seq(0.00001, 0.3, length.out = sims)
Th.range<-seq(1/30, 1/5, length.out = sims)
nn.range<-seq(1, 3, length.out = sims)
m.range<-seq(30, 600, length.out = sims)
v.range<-seq(0.02, .20, length.out = sims)
pi_M.range<-seq(0, 1.0, length.out = sims)
theta.range<-seq(10, 500, length.out = sims)
pi_C.range<-seq(0, 1.0, length.out = sims)
beta.range<-seq(1e-7, 1e-3, length.out = sims)
sigma.range<-seq(1/10, 1/60, length.out = sims)
lamda.range<-seq(1e-6, 1e-3, length.out = sims)
k.range<-seq(0.01, 1.0, length.out = sims)
mu_W.range<-seq(1/365, 1/(365*10), length.out = sims)
mu_H.range<-seq(1/(365*10), 1/(365*80), length.out = sims)

#parameter names
paranges<-cbind("f_N" = f_N.range,
                "phi_N" = phi_N.range,
                "z" = z.range,
                "mu_N" = mu_N.range,
                "mu_I" = mu_I.range,
                "f_P" = f_P.range,
                "phi_P" = phi_P.range,
                "mu_P" = mu_P.range,
                "alpha" = alpha.range,
                "Th" = Th.range,
                "nn" = nn.range,
                "m" = m.range,
                "v" = v.range,
                "pi_M" = pi_M.range,
                "theta" = theta.range,
                "pi_C" = pi_C.range,
                "beta" = beta.range,
                "sigma" = sigma.range,
                "lamda" = lamda.range,
                "k" = k.range,
                "mu_W" = mu_W.range,
                "mu_H" = mu_H.range)

constantparams<-matrix(ncol = length(parameters), nrow = sims)

for(i in 1:length(parameters)){
  constantparams[,i] = rep(parameters[i],sims)
}

colnames(constantparams)<-names(parameters)
vars<-colnames(paranges)
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
  
  area = prams['A']
  sigma = prams['sigma']
  lamda = prams['lamda']
  omega = prams['Om']
  theta = prams['theta']
  pi_C = prams['pi_C']
  beta = prams['beta']
  H = prams['H']
  m = prams['m']
  v = prams['v']
  pi_M = prams['pi_M']
  mu_N = prams['mu_N']
  mu_I = prams['mu_I']
  mu_H = prams['mu_H']
  mu_W = prams['mu_W']
  mu_P = prams['mu_P']
  f_P = prams['f_P']
  phi_P = prams['phi_P']
  alpha = prams['alpha']
  Th = prams['Th']
  f_N = prams['f_N']
  phi_N = prams['phi_N']
  
  f_Nq = f_N
  muPq = mu_P
  phi_Nq = phi_N
  mu_Nq = mu_N
  alpha_q = alpha
  theta_q = theta
  pi_Mq = pi_M
  pi_Cq = pi_C
  v_q = v

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

nil1<-function(In){
  return(1)
}
nil0<-function(In){
  return(0)
}

#First check scatter plots of outcomes across each parameter range while holding other parameters equal ############
outputfillr0<-matrix(ncol = length(vars), nrow = sims)
outputfillN.eq<-matrix(ncol = length(vars), nrow = sims)

for(j in 1:length(vars)){
  for(i in 1:sims){
    print(c(j, i))
    
    other<-constantparams[, -which(colnames(constantparams) %in% vars[j])]
    test<-paranges[, which(colnames(paranges) %in% vars[j])]
    
    parametersuse<-cbind(other, test) 
    colnames(parametersuse)[dim(parametersuse)[2]]<-vars[j]
    
    prams<-parametersuse[i,]
    
    outputfillr0[i,j] = r0.In(conc=0)[3]
    outputfillN.eq[i,j] = r0.In(conc=0)[1]
  }
  
  mypath <- file.path("C:","Users","chris_hoover","Documents","RemaisWork","Schisto","R Codes",
                      "ag_schist","Review_models","Sensitivity_Plots", "r0",
                      paste("r0_only-sens_", vars[j], ".jpg", sep = ""))
  
  jpeg(file=mypath, width = 750, height = 625, units = "px")
  
  par(mfrow = c(1,1), mar = c(4,3.75,1,0.4)+0.1)
  
  plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfillr0[,j], 
       xlab = vars[j], ylab = 'r0', type = 'l', lwd=2,
       pch = 16, cex = 0.75, col = 'navy', 
       ylim = c(0,5))
  
  #plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfillN.eq[,j], 
   #    xlab = vars[j], ylab = 'N.eq', type = 'l', lwd=2,
    #   pch = 16, cex = 0.75, col = 'olivedrab', 
     #  ylim = c(0,60))
  
  dev.off()
  
}