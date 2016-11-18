require(deSolve)

snails = function(t, n, parameters){
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    
    N = S + E + I
    
    #mating function
    fx<-function(x, mean.worm = W, clump = k){
      (1-cos(x))/((1+((mean.worm/(clump+mean.worm))*cos(x)))^1+clump)
    }
    gamma = 1 - ((1-(W/(W+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value)       
    
    mu_Nq = mu_N #*f(Q)
    phi_Nq = phi_N
    v_q = v
    pi_Mq = pi_M
    pi_Cq = pi_C
    theta_q = theta
    pred = 0
    
    #Schistosome larval concentration equations
    Wf = 0.5*W*H*gamma
    
    M = Wf*m*v_q*pi_Mq 
    
    C = theta_q*I*pi_Cq 
    
    dSdt= f_N*(1-N/(phi_Nq*A))*(S+z*E) - mu_Nq*S - pred*(S/A)^nn - beta*M*Om*S       #Susceptible snails
    
    dEdt= beta*M*S - mu_Nq*E - pred*(E/A)^nn - sigma*E                           #Exposed snails
    
    dIdt= sigma*E - mu_Nq*I - mu_I*I - pred*(I/A)^nn                             #Infected snails
    
    dWdt= lamda*Om*C - (mu_W+mu_H)*W                                      #mean worm burden in human population
    
    return(list(c(dSdt, dEdt, dIdt, dWdt)))
  })
}

area=200 #m^2
snstart = c(S=30*area, 
            E=0, 
            I=0, 
            W=5)

time = seq(0,365*30,1)

sparms=c(
  # Location parameters
  A = area,          # Area of site of interest, m^2
  H = 300,           # Human population at site of interest
  Om = 1,            # degree of overlap between water contamination, snail, and human habitats
  
  # Snail reproductive parameters
  f_N = 0.10,        # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection - more like a recruitment rate); 
                     #   from Sokolow et al. 2015 
  phi_N = 50,        # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
  z = 0.7,           # Fraction of exposed snails that reproduce
  
  # Snail mortality parameters
  mu_N = 1/60,        # Natural mortality rate of large snails (deaths/snail/day; assume mean lifespan = 60 days)
  mu_I = 1/10,        # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
  
  nn = 2,            # exponent of the Holling's type III functional response
  
  # miracidia parameters
  m = 432,           # Miracidial shedding rate per adult female worm assuming 0.36 eggs/mL urine and 1200 mL urine per person per day
  v = 0.084,         # Egg viability of S. haematobium (i.e. miracidia/egg)
  pi_M = 1,          # Miracidial infectivity parameter
  
  # cercariae parameters
  theta = 109,       # cercarial shedding rate in snails (cercariae/I-snail/day); Pfluger 1984
  pi_C = 1,          # cercarial infectivity parameter
  
  # transmission parameters
  beta = 1e-7,       # Human-to-snail infection probability in reference area (infected snails/miracidia/snail/day)
  sigma = 1/40,      # Latent period for exposed snails (infectious snails/exposed snail/day))
  lamda = 1e-6,      # Snail-to-human infection probability per cercaria
  k=0.2,             # Clumping parameter of negative binomial distribution of worms in humans
  
  # Schisto mortality parameters
  mu_W = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  mu_H = 1/(60*365)   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
)

soutput = as.data.frame(ode(snstart, time, snails, sparms))
  soutput$N = (soutput$S + soutput$E + soutput$I)
  
  plot(x = soutput$time, y = soutput$N, lwd=2, xlab = 'time', ylab = 'snail dynamics', 
       type = 'l', ylim = c(0, max(soutput$N)))
    lines(soutput$time, soutput$S, lwd=2, col = 'blue')
    lines(soutput$time, soutput$E, lwd=2, col = 'orange')
    lines(soutput$time, soutput$I, lwd=2, col = 'red')
  
  plot(x = soutput$time, y = soutput$W, lwd=2, xlab = 'time', ylab = 'mean worm burden', 
       type = 'l', col = 'purple', ylim = c(0, max(soutput$W)))
  
snstart2 = c(S=soutput[max(soutput$time),c(2)], 
             E=soutput[max(soutput$time),c(3)], 
             I=soutput[max(soutput$time),c(4)], 
             W=soutput[max(soutput$time),c(5)])  
  
#Halt reproduction and assess snail pop #############
  sparms['f_N'] = 0
  time2 = seq(0,365,1)
  
soutput2 = as.data.frame(ode(snstart2, time2, snails, sparms))
  soutput2$N = (soutput2$S + soutput2$E + soutput2$I)
  
  plot(x = soutput2$time, y = soutput2$N, lwd=2, xlab = 'time', ylab = 'snail dynamics', 
       type = 'l', ylim = c(0, max(soutput2$N)))
    lines(soutput2$time, soutput2$S, lwd=2, col = 'blue')
    lines(soutput2$time, soutput2$E, lwd=2, col = 'orange')
    lines(soutput2$time, soutput2$I, lwd=2, col = 'red')
    
# Calculate PDF f(t) = -dS(t)/dt (i.e. numerically approximate derivative of survival curve) ###########
    vect = soutput2$N
    deriv = numeric(length(vect)-1)
    for (i in 1:(length(vect)-1)) {
      deriv[i] = vect[i] - vect[i+1]
    }
    
    deriv2 = deriv/soutput2$N[1]
    deriv2 = deriv2[-(length(deriv2))]
    sum(deriv2) 
    
    x = 0:(length(deriv2)-1)
    
    plot(x = x, y = deriv2, type = 'l')
    
    # Calculate life expectancy
    sum(x*deriv2)   #Mean life expectancy is about a month, that's what we want    