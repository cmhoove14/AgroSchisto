# Differential equation model with dynamic k and lambda ############
agrochem_mda_mod_fit=function(t, n, parameters) { 
  with(as.list(parameters),{
  #set k based on time in simulation  
  k <- ifelse(t < timepoints[1], W_baseline_k, 
           ifelse(timepoints[1] <= t & t < timepoints[2], W_5month_field_k,
                  ifelse(timepoints[2] <= t & t < timepoints[3], W_Feb13_field_k,
                         ifelse(timepoints[3] <= t, W_Sep13_field_k, NA))))
    
  #set lambda (snail to human transmission) based on time in simulation  
  lambda <- ifelse(t < seasons[1], parameters[["lambda2"]], 
               ifelse(seasons[1] <= t & t < seasons[2], parameters[["lambda1"]],
                      ifelse(seasons[2] <= t & t < seasons[3], parameters[["lambda2"]],
                             ifelse(seasons[3] <= t, parameters[["lambda1"]], NA))))
    
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    P=n[6]
    
#Dynamic variables
    N = S+E+I  #Total number of snails
    
    W = cvrg*Wt + (1-cvrg)*Wu
    
    phi = phi_Wk(W, kappa)           

    pred_S = ((alpha*eps)*(S/A)^nn) / (1 + sum((alpha*eps)*Th*(S/A)^nn,
                                               (alpha*eps)*Th*(E/A)^nn,
                                               (alpha*eps)*Th*(I/A)^nn))
    
    pred_E = ((alpha*eps)*(E/A)^nn) / (1 + sum((alpha*eps)*Th*(S/A)^nn,
                                               (alpha*eps)*Th*(E/A)^nn,
                                               (alpha*eps)*Th*(I/A)^nn))
    
    pred_I = ((alpha*eps)*(I/A)^nn) / (1 + sum((alpha*eps)*Th*(S/A)^nn,
                                               (alpha*eps)*Th*(E/A)^nn,
                                               (alpha*eps)*Th*(I/A)^nn))    
    
#agrochemical parameters adjusted
    Knadj = K_N * Knq
    fnadj = f_N * fnq
    munadj = mu_N + munq
    vadj = v * vq
    pimadj = pi_M * pimq
    thetaadj = theta * thetaq
    picadj = pi_C * picq
    mupadj = mu_P + mupq
    psiadj = psiq
    
    pred_Sadj = pred_S * psiadj
    pred_Eadj = pred_E * psiadj
    pred_Iadj = pred_I * psiadj
    
#Schistosome larval concentration equations
    M = 0.5*W*H*phi*m*vadj*pimadj #- (mu_M + mu_Mq)*M - beta*M*N 
    C = thetaadj*I*picadj #- (mu_C + mu_Cq)*C - lambda*H*C
    
    #differential equations
    
    dSdt= fnadj*(1-(N/Knadj))*(S+E) - munadj*S - pred_Sadj*P - beta*M*S  #Susceptible snails
    
    dEdt= beta*M*S - munadj*E - pred_Eadj*P - sigma*E                                     #Exposed snails
    
    dIdt= sigma*E - (munadj + mu_I)*I - pred_Iadj*P                                       #Infected snails
    
    dWtdt= lambda*C - (mu_W+mu_H)*Wt                                      #mean worm burden in treated human population
    dWudt= lambda*C - (mu_W+mu_H)*Wu                                  #mean worm burden in untreated human population

    dPdt= f_P*(1-(P/K_P))*P - mupadj*P                                      #prawn population (number individuals)

    
    return(list(c(dSdt, dEdt, dIdt, dWtdt, dWudt, dPdt)))
  }) 
} 

# Differential equation model with dynamic k and seasonality in snail reproduction ############
agrochem_mda_seasonal_mod_fit=function(t, n, parameters) { 
  with(as.list(parameters),{
  #set k based on time in simulation  
  k <- ifelse(t < timepoints[1], W_baseline_k, 
           ifelse(timepoints[1] <= t & t < timepoints[2], W_5month_field_k,
                  ifelse(timepoints[2] <= t & t < timepoints[3], W_Feb13_field_k,
                         ifelse(timepoints[3] <= t, W_Sep13_field_k, NA))))
    
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    P=n[6]
    
#Dynamic variables
    N = S+E+I  #Total number of snails
    
    W = cvrg*Wt + (1-cvrg)*Wu
    
    phi = phi_Wk(W, kappa)           

    pred_S = ((alpha*eps)*(S/A)^nn) / (1 + sum((alpha*eps)*Th*(S/A)^nn,
                                               (alpha*eps)*Th*(E/A)^nn,
                                               (alpha*eps)*Th*(I/A)^nn))
    
    pred_E = ((alpha*eps)*(E/A)^nn) / (1 + sum((alpha*eps)*Th*(S/A)^nn,
                                               (alpha*eps)*Th*(E/A)^nn,
                                               (alpha*eps)*Th*(I/A)^nn))
    
    pred_I = ((alpha*eps)*(I/A)^nn) / (1 + sum((alpha*eps)*Th*(S/A)^nn,
                                               (alpha*eps)*Th*(E/A)^nn,
                                               (alpha*eps)*Th*(I/A)^nn))    
    
#agrochemical parameters adjusted
    Knadj = K_N * Knq
    fnadj = f_N * fnq
    munadj = mu_N + munq
    vadj = v * vq
    pimadj = pi_M * pimq
    thetaadj = theta * thetaq
    picadj = pi_C * picq
    mupadj = mu_P + mupq
    psiadj = psiq
    
    pred_Sadj = pred_S * psiadj
    pred_Eadj = pred_E * psiadj
    pred_Iadj = pred_I * psiadj
    
#Schistosome larval concentration equations
    M = 0.5*W*H*phi*m*vadj*pimadj #- (mu_M + mu_Mq)*M - beta*M*N 
    C = thetaadj*I*picadj #- (mu_C + mu_Cq)*C - lambda*H*C
    
    #differential equations
    
    dSdt= fnadj*(1-(N/Knadj))*(S+E) - munadj*S - pred_Sadj*P - beta*M*S  #Susceptible snails
    
    dEdt= beta*M*S - munadj*E - pred_Eadj*P - sigma*E                                     #Exposed snails
    
    dIdt= sigma*E - (munadj + mu_I)*I - pred_Iadj*P                                       #Infected snails
    
    dWtdt= lambda*((cos(2*pi*(t+180)/365) + 2)/2)*C - (mu_W+mu_H)*Wt                                      #mean worm burden in treated human population
    dWudt= lambda*((cos(2*pi*(t+180)/365) + 2)/2)*C - (mu_W+mu_H)*Wu                                  #mean worm burden in untreated human population

    dPdt= f_P*(1-(P/K_P))*P - mupadj*P                                      #prawn population (number individuals)

    
    return(list(c(dSdt, dEdt, dIdt, dWtdt, dWudt, dPdt)))
  }) 
} 

## Likelihood function for a single data point:
  pointLogLike <- function(obs_mu = measured_ws, #Measured worm burdens
                           obs_sd = w_errors,       #Measured worm burden uncertainty
                           modOut) {                #Model output W
    log(Prob_gaussian(y=modOut, mu=obs_mu, sd=obs_sd))
  }
 
#Get distributions for Epi datapoints given estimated W and k 
Prob_negbin<-function(i,k,m){
  p<-(gamma(k+i)/(gamma(i+1) * gamma(k))) * (1 + (m/k))^(-k-i) * ((m/k)^i)
  p
}
# i number of parasites per person
# k clumping parameter
# m mean burden
# refer wikipedia page - https://en.wikipedia.org/wiki/Negative_binomial_distribution
# formula refer - http://influentialpoints.com/Training/negative_binomial_distribution-principles-properties-assumptions.htm

Prob_gaussian<-function(y, mu, sd){
  p<-(1/(sqrt(2*pi)*sd) ) * exp( -((y-mu)^2)/(2*sd^2) )
  p
}
 