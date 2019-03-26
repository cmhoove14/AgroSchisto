#Agrochem model parameters

area = 200

init_pars=c(
  # Location parameters
  A = area,          # Area of site of interest, m^2
  H = 1.5*area,      # Human population at site of interest (based on 300 people at 200m^2 water contact site from Sokolow PNAS)
  #Om = 1 / sqrt(area),# degree of overlap between water contamination, snail, and human habitats
  
  # Snail reproductive parameters
  f_N = 0.1,        # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection - more like a recruitment rate); 
                     #   from Woolhouse & Chandiwana et al. 1990 
  K_N = 50*area,   # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
  z = 0.7,           # Fraction of exposed snails that reproduce
  
  # Snail mortality parameters
  mu_N = 1/60,        # Natural mortality rate of large snails (deaths/snail/day; mean lifespan ~ 60 days (~2 months))
  mu_I = 1/10-1/60,   # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
  
  # Predator pop dynamic parameters (from pred tweaking code)
  f_P = 0.117,  #Predator intrinsic recruitment rate assuming 30000 eggs per year and 1% survival from eggs to repro maturity
                  #   from Adha-Ar et al 2016
  K_P = 0.125*area,  #Predator CC (0.125/m^2 from Sokolow et al realized predator density, discounted here for natural pop)
  mu_P = 0.038,       #Predator mortality rate from Halstead et al paper
  
  # Predation parameters
  alpha = 2,       # Predator attack rate at high prawn/snail weight ratio per Sokolow 2014 Acta Tropica
  eps = 0.1,       #Attack rate penalty estimated as 90% reduction from lab conditions  
  Th = 0.1,        # Predator handling time at high prawn/snail weight ratio per Sokolow 2014 Acta Tropica
  nn = 1,          # exponent of the Holling's type III functional response
  
  # miracidia parameters
  m = 432,           # Miracidial shedding rate per adult female worm assuming 0.36 eggs/mL urine and 1200 mL urine per person per day
  v = 0.084,         # Egg viability of S. haematobium (i.e. miracidia/egg) from Halstead et al
  pi_M = 6.22,       # mean Miracidia-hrs per day in agrochemical-free water
  
  # cercariae parameters
  theta = 109/7,       # mean daily cercarial shedding rate of patently infected snails @25C (cercariae/I-snail/day); Pfluger 1984
  pi_C = 14.21,        # mean cercariae-hrs in agrochemical-free water
  
  # transmission parameters
  beta = 1.6e-5/24,       # Human-to-snail infection probability in reference area (exposed snails/susceptible snail/miracidi-hr/day); 
                        #     divided by 24 to account for hourly scale of miracidial survival; from fitting procedure in Halstead et al
  sigma = 1/40,         # Latent period for exposed snails (infectious snails/exposed snail/day))
  lambda = 3.7e-6/24, # Snail-to-human infection probability in reference area (adult worms/cercariae-hr); 
                        #     divided by 24 to account for hourly scale of miracidial survival; from fitting procedure in Halstead et al
  kappa = 0.08,               # Clumping parameter of negative binomial distribution of worms in humans
  
  # Schisto mortality parameters
  mu_W = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  mu_H = 1/(60*365),   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
  
  #treatment parameters
  cvrg=0.43, #coverage of treatment across the population, Lampsar I = 100/1000 =0.1 %, Lampsar II = 129/300 = 43%
  eff=0.94, # efficacy of the drug
  
  #agrochemical parameters set to null values
  Knq = 1,
  fnq = 1,
  munq = 0,
  vq = 1,
  pimq = 1,
  thetaq = 1,
  picq = 1,
  mupq = 0,
  psiq = 1

)
