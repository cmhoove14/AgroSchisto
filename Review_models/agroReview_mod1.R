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

require(deSolve)

mod1 = function(t, n, parameters) {
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    P=n[5]
    Q1=n[6]
    
    #agrochemical parameters 
    f_red = exp(-0.0028479*Q1) #from negative exponential fit to Ibrahim 1991 data
    
    mu_Pq = (exp(-3.074e+02*exp(-Q1)))/10 #from gompertz curve fit to Halstead 2015 10-day mortality
    
    pred_red = 1 #NO DATA FOR CHLORPYRIFOS; see Yuan 2004 for Malathion & Paraquat toxicity
    
    i_Cq = 1 #NO DATA FOR CHLORPYRIFOS; see Rohr 2008 for Malathion, Carbaryl, Atrazine toxicity
    
    i_Mq = 1 #NO DATA FOUND
    
    v_Mq = 1 #NO DATA FOUND
    
    #Other agrochemical variables
    phi_Nq = 1
    
    #Dynamic variables
    N = S+E+I  #Total number of snails
    
    fx<-function(x, mean.worm = W, clump = k){
      (1-cos(x))/((1+((mean.worm/(clump+mean.worm))*cos(x)))^1+clump)
    }
    gamma = 1 - ((1-(W/(W+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value)           #mating function
    
    # Mean and total prawn biomass, converting from length (mm) to weight (g)
    Bm.p = (a.p*(L/10)^b.p)/10  
    Bm.t = P*Bm.p
    
    #Predation
    Bm.n = a.s*0.8^b.s  # 8mm size class
    
    Bm.r = Bm.p / Bm.n
    
    alpha_star = ar*Bm.r
    
    alpha = alpha_star/A                                                            #alpha corrected for natural setting and low density
    
    Th = 1/(th*Bm.r)                                                                #handling time 
    
    pred= (alpha*P)/(1+(alpha*((N/A)^nn)*(Th*pred_red)))                            #per capita predation of predators on snails; a Holling's type III functional response
    
    print(pred)
    
    #Schistosome larval concentration equations
    M = 0.5*W*H*gamma*m*(v_M*v_Mq)*i_Mq #- (mu_M + mu_Mq)*M - beta*M*N 
    C= theta*I*i_Cq #- (mu_C + mu_Cq)*C - lamda*H*C
    
    #differential equations
    
    dSdt= (f_N*f_red)*(1-(N/(phi_N*phi_Nq*A)))*(S+E) - mu_N*S - pred*((S/A)^nn) - beta*M*S  #Susceptible snails
    
    dEdt= beta*M*S - mu_N*E - pred*((E/A)^nn) - sigma*E                                     #Exposed snails
    
    dIdt= sigma*E - (mu_N + mu_I)*I - pred*((I/A)^nn)                                       #Infected snails
    
    dWdt= lamda*C - (mu_W+mu_H)*W                                                       #mean worm burden in human population
    
    dPdt= -P*(mu_P*L^d1 + mu_Pq*L^d2 + Bm.t/phi_P)                                      #prawn population (number individuals)
    
    dLdt= k_P/(1+gam*Bm.t)*(linf - L)                                                   #mean prawn length (mm)
    
    dQ1dt= -k_Q1*Q1                                                                     #Agrochemical concentration (insecticide)
 
    
    return(list(c(dSdt, dEdt, dIdt, dWdt, dPdt, dLdt, dQ1dt)))
  }) 
}

#Initial values and parameters #####################
area=200
nstart1 = c(S=15*area, 
           E=10*area, 
           I=5*area, 
           W=5, 
           P=0,#2*area, 
           L=25, 
           Q1=0)
yrs=30
time = seq(0,365*yrs,1)

parameters=c(
  # Location parameters
  A = area,          # Area of site of interest, m^2
  H = 300,           # Human population at site of interest
  
  # Snail reproductive parameters
  f_N = 0.10,        # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection - more like a recruitment rate); 
                     #   from Sokolow et al. 2015 
  phi_N = 50,        # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015

  # Snail mortality parameters
  mu_N = 1/60,        # Natural mortality rate of large snails (deaths/snail/day; assume mean lifespan = 50 days)
  mu_I = 1/10,        # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
  
  # Prawn growth parameters
  a.p = 0.096868,    # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, growout phase)
  b.p = 3.2944,      # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, growout phase)
  k_P = 0.00339726,  # Growth coefficient, from Nwosu & Wolfi 2006 (M. vollenhovenii); alternate value for M. rosenbergii, from Sampaio & Wagner 1996: 0.0104333333
  linf = 206,        # Max length (mm), from Nwosu & Wolfi 2006 (M. vollenhovenii)
  gam = 1e-6,        # Density-dependent growth parameter (based on biomass per hectare); not yet fit
  
  # Prawn mortality parameters
  mu_P = 0.006136986,# Natural prawn mortality rate, from Nwosu & Wolfi 2006 (M. vollenhovenii)
  d1 = -0.25,        # Exponential parameter relating size with mortality; no source
  d2 = -0.1,         # Exponential parameter relating size with insecticide mortality; no source
  phi_P = 1000*area, # Density-dependent mortality parameter (based on biomass per hectare); not yet fit
  
  # Predation parameters
  a.s = 0.187178454, # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
  b.s = 2.536764792, # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
  ar = 0.037192,     # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
  th = 0.40450,      # Coefficient for relationship between biomass ratio and handling time, fitted to data from Sokolow et al. 2014
  nn = 2,            # exponent of the Holling's type III functional response
  
  
  # miracidia parameters
  m = 432,           # Miracidial shedding rate per adult female worm assuming 0.36 eggs/mL urine and 1200 mL urine per person per day
  v_M = 0.084,       # Egg viability of S. haematobium
  
  # cercariae parameters
  theta = 109,       # cercarial shedding rate in snails (cercariae/I-snail/day); Pfluger 1984

  # transmission parameters
  beta = 1e-7,       # Human-to-snail infection probability in reference area (infected snails/miracidia/snail/day)
  sigma = 1/40,      # Latent period for exposed snails (infectious snails/exposed snail/day); adjusted from Sokolow et al. 2015 (original value: 1/50)
  lamda = 1e-6,      # Snail-to-human infection probability scaled to 1 m^2 (composite including mortality, infection, survival to patency)
  k=0.2,             # Clumping parameter of negative binomial distribution of worms in humans
  
  #Agrochemical parameters
  k_Q1 = 0.0126,     # half-life of chlorpyrifos in water according to Cornell database
  
  # Schisto mortality parameters
  mu_W = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  mu_H = 1/(60*365)   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
)

#Model run and plot ####################
output = as.data.frame(ode(nstart1, time, agro_time, parameters))
  output$N = output$S + output$E + output$I
  output$P.B = output$P * (parameters['a.p']*(output$L/10)^parameters['b.p'])/10
  eqbm<-output[max(output$time),c(2:8)]
  
par(mfrow = c(3,1), mar = c(3.0,4,1,2.0))
plot(x = output$time, y = output$N, lwd=2, xlab = 'time', ylab = 'snail dynamics', 
     type = 'l', ylim = c(0, max(output$N)))
  lines(output$time, output$S, lwd=2, col = 'blue')
  lines(output$time, output$E, lwd=2, col = 'orange')
  lines(output$time, output$I, lwd=2, col = 'red')
  
plot(output$time, output$P, lwd=2, type = 'l', col = 'green', 
     xlab = 'time', ylab = 'prawn', ylim = c(0, max(output$P.B/10)))
  lines(output$time, output$P.B/10, lwd=2, col = 'darkgreen')
  lines(output$time, output$Q1,lwd=2, col = 'yellow')
  
plot(output$time, output$W, lwd=2, col = 'purple', type = 'l', xlab = 'time', ylab = 'mean worm burden', ylim = c(0,100))  

#Model run adding in predators with starting values=above eqbm ##################
nstart2 = c(S=output[max(output$time),c(2)], 
              E=output[max(output$time),c(3)], 
              I=output[max(output$time),c(4)], 
              W=output[max(output$time),c(5)], 
              P=2*area, 
              L=25, 
              Q1=0)
yrs=2
time = seq(0,365*yrs,1)
  
output2 = as.data.frame(ode(nstart2, time, agro_time, parameters))
  output2$N = output2$S + output2$E + output2$I
  output2$P.B = output2$P * (parameters['a.p']*(output2$L/10)^parameters['b.p'])/10
  output2[max(output2$time),]
  
par(mfrow = c(3,1), mar = c(3.0,4,1,2.0))
plot(x = output2$time, y = output2$N, lwd=2, xlab = 'time', ylab = 'snail dynamics', 
       type = 'l', ylim = c(0, max(output2$N)))
  lines(output2$time, output2$S, lwd=2, col = 'blue')
  lines(output2$time, output2$E, lwd=2, col = 'orange')
  lines(output2$time, output2$I, lwd=2, col = 'red')
  
plot(output2$time, output2$P, lwd=2, type = 'l', col = 'green', 
       xlab = 'time', ylab = 'prawn dynamics', ylim = c(0, max(output2$P.B/10)))
  lines(output2$time, output2$P.B/10, lwd=2, col = 'darkgreen')
  lines(output2$time, output2$Q1,lwd=2, col = 'yellow')
  
plot(output2$time, output2$W, lwd=2, col = 'purple', type = 'l', xlab = 'time', ylab = 'mean worm burden', ylim = c(0,100))  
  
#Model run adding in chlorpyrifos(Q1) ##################
nstart3 = c(S=output[max(output$time),c(2)], 
              E=output[max(output$time),c(3)], 
              I=output[max(output$time),c(4)], 
              W=output[max(output$time),c(5)], 
              P=2*area, 
              L=25, 
              Q1=20)
yrs=2
time = seq(0,365*yrs,1)
  
output3 = as.data.frame(ode(nstart3, time, agro_time, parameters))
  output3$N = output3$S + output3$E + output3$I
  output3$P.B = output3$P * (parameters['a.p']*(output3$L/10)^parameters['b.p'])/10
  output3[max(output3$time),]
  
par(mfrow = c(3,1), mar = c(3.0,4,1,2.0))
  plot(x = output3$time, y = output3$N, lwd=2, xlab = 'time', ylab = 'snail dynamics', 
       type = 'l', ylim = c(0, max(output3$N)))
  lines(output3$time, output3$S, lwd=2, col = 'blue')
  lines(output3$time, output3$E, lwd=2, col = 'orange')
  lines(output3$time, output3$I, lwd=2, col = 'red')
  
plot(output3$time, output3$P, lwd=2, type = 'l', col = 'green', 
       xlab = 'time', ylab = 'prawn', ylim = c(0, max(output3$P.B/10)))
  lines(output3$time, output3$P.B/10, lwd=2, col = 'darkgreen')
  lines(output3$time, output3$Q1,lwd=2, col = 'yellow')
  
plot(output3$time, output3$W, lwd=2, col = 'purple', type = 'l', xlab = 'time', ylab = 'mean worm burden', ylim = c(0,100))  

#Insecticide info for modeling concentration over time ################
#calculate decay rates (k) from hydrolysis half lives for each chemical in Halstead 2015 (from table S1 and pmep.cce.cornell.edu)
  mal.k = -log(0.5)/6.2         #from table s1 and in agreement of "less than 1 week in raw river water" from Cornell
  chlor.k = -log(0.5)/25.5      #from table S1; within the range of reported half life from cornell
  terb.k = -log(0.5)/6.5        #from table s1; in agreement with Cornell estimate of 5.5 days at pH of 7; degrades into formaldehyde
  lamcy.k = -log(0.5)/1         #very fast; 0 in table S1; "Not expected to be prevalent in surface waters" according to cornell website
  esfen.k = -log(0.5)/10        #cornell 4-15 days half life in water
  perm.k = -log(0.5)/2          #cornell "half life of less than 2.5 days"

#median EECs from Halstead 2015
  med.mal = 0.778
  med.chlor = 5.810
  med.terb = 1.435
  med.lamcy = 0.649
  med.esfen = 0.311
  med.perm = 1.420
  
#suggested application intervals for each agrochemical from Halstead 2015
  mal.days = c(1,6)             #two applications 5 days apart
  chlor.days = c(1,11,21)       #three applications 10 days apart
  terb.days = c(1)              #single application
  lamcy.days = seq(1, 61, by=4) #16 applications 4 days apart
  esfen.days = seq(1, 21, by=5) #5 applications 5 days apart
  perm.days = seq(1, 36, by=5)  #8 applications 5 days apart
  
#Model structure to dynamically vary agrochemical application, MDA introduction, and prawn introduction ####################
  nstart.eqbm = c(S=output[max(output$time),c(2)], 
                  E=output[max(output$time),c(3)], 
                  I=output[max(output$time),c(4)], 
                  W=output[max(output$time),c(5)], 
                  P=0, 
                  L=0, 
                  Q1=0)

agroc_prawn_mda<-function(agroc.add, chem, mda.add, coverage, efficacy, prawn.add, prawn.n, prawn.l){
  #first run for a year
  time1 = c(1:365)
  output.all<-as.data.frame(ode(nstart.eqbm, time1, agro_time, parameters))
}