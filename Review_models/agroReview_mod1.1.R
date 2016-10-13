#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Schisto model with time-variant agrochemical concentration
  #To do list
  #(1) Incorporate agrochemical functions
  #(2) Fit prawn pop dynamics to data
  #(3) Fit prawn attack rate on snails to PNAS data?

require(deSolve)

#Agrochemical response functions ###############
  f_phi_Nq = function(Q1){
    phi_Nq = phi_N + phi_N * (atra.slope*log(Q1))
    return(phi_Nq)
  }

#Model structure and equations #############
mod1 = function(t, n, parameters) {
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    P=n[5]
    Q1=n[6]
    
    
#Dynamic variables ####################
    #Total number of snails
      N = S+E+I  
    #mating function
      fx<-function(x, mean.worm = W, clump = k){
        (1-cos(x))/((1+((mean.worm/(clump+mean.worm))*cos(x)))^1+clump)
      }
      gamma = 1 - ((1-(W/(W+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value)       
    
    #per capita predation of predators on snails with variable functional response based on parameter nn  
      pred= (alpha*P)/(1+(alpha*((N/A)^nn)*Th))   
      
    #Agrochem parameters
      phi_Nq = f_phi_Nq(Q1)
      mu_Nq = mu_N #* f2(Q1)
      mu_Pq = mu_P #+ exp(-3.074e+02*exp(-Q1))/10
      v_q = v #* f4(Q1)
      theta_q = theta #* f5(Q1)
      pi_Mq = pi_M #* f6(Q1)
      pi_Cq = pi_C #* f7(Q1)
    
    #Schistosome larval concentration equations
    Wf = 0.5*W*H*gamma
    
    M = Wf*m*v_q*pi_Mq 
    
    C = theta_q*I*pi_Cq 
    
    #print(c(paste('M = ', round(M, digits = 2), sep = ''),
    #        paste('C = ', round(C, digits = 2), sep = ''),
    #        paste('gam = ', round(gamma, digits = 2), sep = ''),
    #        paste('pred = ', pred, sep = ''),
    #paste('N = ', round(N, digits = 2), sep = '')))
    
#differential equations ###################
    
    dSdt= f_N*(1-N/(phi_Nq*A))*(S+z*E) - mu_Nq*S - pred*(S/A)^nn - beta*M*Om*S       #Susceptible snails
    
    dEdt= beta*M*S - mu_Nq*E - pred*(E/A)^nn - sigma*E                           #Exposed snails
    
    dIdt= sigma*E - mu_Nq*I - mu_I*I - pred*(I/A)^nn                             #Infected snails
    
    dWdt= lamda*Om*C - (mu_W+mu_H)*W                                      #mean worm burden in human population
    
    dPdt= f_P*(1-P/(phi_P*A))*P - mu_Pq*P                                    #prawn population (number individuals)
    
    dQ1dt= -k_Q1*Q1                                                    #Agrochemical concentration
 
    
    return(list(c(dSdt, dEdt, dIdt, dWdt, dPdt, dQ1dt)))
  }) 
}

#Initial values and parameters #####################
area=1 #m^2
nstart1 = c(S=5*area, 
           E=0, 
           I=0, 
           W=5, 
           P=0, 
           Q1=0)
yrs=30
time = seq(0,365*yrs,1)

parameters=c(
  # Location parameters
  A = area,          # Area of site of interest, m^2
  H = 1.5*area,      # Human population at site of interest (based on 300 people at 200m^2 water contact site from Sokolow PNAS)
  Om = 1,            # degree of overlap between water contamination, snail, and human habitats
  
  # Snail reproductive parameters
  f_N = 0.10,        # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection - more like a recruitment rate); 
                     #   from Sokolow et al. 2015 
  phi_N = 50,        # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
  z = 0.7,           # Fraction of exposed snails that reproduce

  # Snail mortality parameters
  mu_N = 1/60,        # Natural mortality rate of large snails (deaths/snail/day; assume mean lifespan = 60 days)
  mu_I = 1/10,        # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
  
  # Predator pop dynamic parameters (from pred tweaking code)
  f_P = 0.02,          #Predator intrinsic recruitment rate
  phi_P = 0.5,         #Predator carrying capacity
  mu_P = 0.0026,       #Predator mortality rate
  
  # Predation parameters
  alpha = 0.003,     # predator attack rate on snails
  Th = 1/14,         # Predator handling time
  nn = 2,            # exponent of the Holling's type III functional response
  
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
  
  #Agrochemical parameters
  k_Q1 = 0.0126,     # half-life of agrochemical (Chlorpyrifos here) in water according to Cornell database
  
  # Schisto mortality parameters
  mu_W = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  mu_H = 1/(60*365)   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
)

#Model run with no predators, no agrochemicals and plot ####################
output1 = as.data.frame(ode(nstart1, time, mod1, parameters))
  output1$N = (output1$S + output1$E + output1$I)

plot(x = output1$time, y = output1$N, lwd=2, xlab = 'time', ylab = 'snail dynamics', 
     type = 'l', ylim = c(0, max(output1$N)))
  lines(output1$time, output1$S, lwd=2, col = 'blue')
  lines(output1$time, output1$E, lwd=2, col = 'orange')
  lines(output1$time, output1$I, lwd=2, col = 'red')
  
plot(output1$time, output1$W, lwd=2, col = 'purple', type = 'l', 
     xlab = 'time', ylab = 'mean worm burden', ylim = c(0,max(output1$W)))  

#Model run adding in predators ##################
nstart2 = c(S=output1[max(output1$time),c(2)], 
              E=output1[max(output1$time),c(3)], 
              I=output1[max(output1$time),c(4)], 
              W=output1[max(output1$time),c(5)], 
              P=0.2*area, 
              Q1=0)
yrs=20
time = seq(0,365*yrs,1)
  
output2 = as.data.frame(ode(nstart2, time, mod1, parameters))
  output2$N = (output2$S + output2$E + output2$I)
  
plot(x = output2$time, y = output2$N, lwd=2, xlab = 'time', ylab = 'snail dynamics', 
       type = 'l', ylim = c(0, max(output2$N)))
  lines(output2$time, output2$S, lwd=2, col = 'blue')
  lines(output2$time, output2$E, lwd=2, col = 'orange')
  lines(output2$time, output2$I, lwd=2, col = 'red')
  
plot(output2$time, output2$P, lwd=2, type = 'l', col = 'green', 
       xlab = 'time', ylab = 'prawn dynamics', ylim = c(0, max(output2$P)))
  
plot(output2$time, output2$W, lwd=2, col = 'purple', type = 'l', xlab = 'time', 
     ylab = 'mean worm burden', ylim = c(0,max(output2$W)))  
  
#Model run adding in example agrochemical, chlorpyrifos (Q1) ##################
nstart3 = c(S=output2[max(output2$time),c(2)], 
            E=output2[max(output2$time),c(3)], 
            I=output2[max(output2$time),c(4)], 
            W=output2[max(output2$time),c(5)], 
            P=output2[max(output2$time),c(6)], 
            Q1=20)
yrs=2
time = seq(0,365*yrs,1)
  
output3 = as.data.frame(ode(nstart3, time, mod1, parameters))
  output3$N = output3$S + output3$E + output3$I

  plot(x = output3$time, y = output3$N, lwd=2, xlab = 'time', ylab = 'snail dynamics', 
       type = 'l', ylim = c(0, max(output3$N)))
  lines(output3$time, output3$S, lwd=2, col = 'blue')
  lines(output3$time, output3$E, lwd=2, col = 'orange')
  lines(output3$time, output3$I, lwd=2, col = 'red')
  
plot(output3$time, output3$P, lwd=2, type = 'l', col = 'green', 
       xlab = 'time', ylab = 'preds', ylim = c(0, max(output3$P)))
  lines(output3$time, output3$Q1, lwd=2, col = 'red', lty = 2)
  legend('topleft', cex = 0.75, lty = 2, col = 'red', lwd = 2, legend = 'log(Q1)')

plot(output3$time, output3$W, lwd=2, col = 'purple', type = 'l', xlab = 'time', 
     ylab = 'mean worm burden', ylim = c(0,max(output3$W)))  