#base model code that contains no agrochemical influences or predators
base_mod = function(t, n, parameters) {
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]

#Dynamic variables
    N = S+E+I  #Total number of snails
    
    phi = phi_Wk(W, kappa)           

#Schistosome larval concentration equations
    M = 0.5*W*H*phi*m*v*pi_M #- (mu_M + mu_Mq)*M - beta*M*N 
    C = theta*I*pi_C #- (mu_C + mu_Cq)*C - lambda*H*C
    
#differential equations
    
    dSdt= f_N*(1-(N/K_N))*(S+E) - mu_N*S - beta*M*S  #Susceptible snails
    
    dEdt= beta*M*S - mu_N*E - sigma*E                                     #Exposed snails
    
    dIdt= sigma*E - (mu_N + mu_I)*I                                      #Infected snails
    
    dWdt= lambda*C - (mu_W+mu_H)*W                                      #mean worm burden in human population
    
    
    return(list(c(dSdt, dEdt, dIdt, dWdt)))

  })
}  

#predator only code that builds on base model to add dynamic predator population
pred_mod = function(t, n, parameters) {
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    P=n[5]
    
#Dynamic variables
    N = S+E+I  #Total number of snails
    
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
    
#    pred_S = (alpha*S)/(1+(alpha*Th*N^nn))  #per capita predation of predators on susceptible snails
#    pred_E = (alpha*E)/(1+(alpha*Th*N^nn))  #per capita predation of predators on exposed snails
#    pred_I = (alpha*I)/(1+(alpha*Th*N^nn))  #per capita predation of predators on infected snails
    
#    pred_S = (alpha*(S/A))/(1+(alpha*Th*(N/A)^nn))  #per capita predation of predators on susceptible snails
#    pred_E = (alpha*(E/A))/(1+(alpha*Th*(N/A)^nn))  #per capita predation of predators on exposed snails
#    pred_I = (alpha*(I/A))/(1+(alpha*Th*(N/A)^nn))  #per capita predation of predators on infected snails
    
#Schistosome larval concentration equations
    M = 0.5*W*H*phi*m*v*pi_M #- (mu_M + mu_Mq)*M - beta*M*N 
    C = theta*I*pi_C #- (mu_C + mu_Cq)*C - lambda*H*C
    
    #differential equations
    
    dSdt= f_N*(1-(N/K_N))*(S+E) - mu_N*S - pred_S*P - beta*M*S  #Susceptible snails
    
    dEdt= beta*M*S - mu_N*E - pred_E*P - sigma*E                                     #Exposed snails
    
    dIdt= sigma*E - (mu_N + mu_I)*I - pred_I*P                                       #Infected snails
    
    dWdt= lambda*C - (mu_W+mu_H)*W                                      #mean worm burden in human population
    
    dPdt= f_P*(1-(P/K_P))*P - mu_P*P                                      #prawn population (number individuals)

    
    return(list(c(dSdt, dEdt, dIdt, dWdt, dPdt)))
  }) 
}

#agrochemical model code to investigate influence of agrochemicals on schisto transmission
agrochem_mod = function(t, n, parameters) {
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    P=n[5]
    
#Dynamic variables
    N = S+E+I  #Total number of snails
    
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
    
    dWdt= lambda*C - (mu_W+mu_H)*W                                      #mean worm burden in human population
    
    dPdt= f_P*(1-(P/K_P))*P - mupadj*P                                      #prawn population (number individuals)

    
    return(list(c(dSdt, dEdt, dIdt, dWdt, dPdt)))
  }) 
}

#agrochemical + mda model code to investigate influence of agrochemicals, predators, mda and their interactions on schisto transmission
agrochem_mda_mod = function(t, n, parameters) {
  with(as.list(parameters),{
    
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

# mda model code to investigate influence of predators, mda and their interactions on schisto transmission
mda_mod = function(t, n, parameters) {
  
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    P=n[6]
    
#Dynamic variables
    N = S+E+I  #Total number of snails
    
    W = parameters["cvrg"]*Wt + (1-parameters["cvrg"])*Wu
    
    phi = phi_Wk(W, parameters["kappa"])           

    pred_S = ((parameters["alpha"]*parameters["eps"])*(S/parameters["A"])^parameters["nn"]) / 
      (1 + sum((parameters["alpha"]*parameters["eps"])*parameters["Th"]*(S/parameters["A"])^parameters["nn"],
               (parameters["alpha"]*parameters["eps"])*parameters["Th"]*(E/parameters["A"])^parameters["nn"],
               (parameters["alpha"]*parameters["eps"])*parameters["Th"]*(I/parameters["A"])^parameters["nn"]))
    
    pred_E = ((parameters["alpha"]*parameters["eps"])*(E/parameters["A"])^parameters["nn"]) / 
      (1 + sum((parameters["alpha"]*parameters["eps"])*parameters["Th"]*(S/parameters["A"])^parameters["nn"],
               (parameters["alpha"]*parameters["eps"])*parameters["Th"]*(E/parameters["A"])^parameters["nn"],
               (parameters["alpha"]*parameters["eps"])*parameters["Th"]*(I/parameters["A"])^parameters["nn"]))
    
    pred_I = ((parameters["alpha"]*parameters["eps"])*(I/parameters["A"])^parameters["nn"]) / 
      (1 + sum((parameters["alpha"]*parameters["eps"])*parameters["Th"]*(S/parameters["A"])^parameters["nn"],
               (parameters["alpha"]*parameters["eps"])*parameters["Th"]*(E/parameters["A"])^parameters["nn"],
               (parameters["alpha"]*parameters["eps"])*parameters["Th"]*(I/parameters["A"])^parameters["nn"]))    
    
#Schistosome larval concentration equations
    M = 0.5*W*parameters["H"]*phi*parameters["m"]*parameters["v"]*parameters["pi_M"] #- (mu_M + mu_Mq)*M - beta*M*N 
    C = parameters["theta"]*I*parameters["pi_C"] #- (mu_C + mu_Cq)*C - lambda*H*C
    
    #differential equations
    
    dSdt= parameters["f_N"]*(1-(N/parameters["K_N"]))*(S+E) - parameters["mu_N"]*S - pred_S*P - parameters["beta"]*M*S  #Susceptible snails
    
    dEdt= parameters["beta"]*M*S - parameters["mu_N"]*E - pred_E*P - parameters["sigma"]*E                                     #Exposed snails
    
    dIdt= parameters["sigma"]*E - (parameters["mu_N"] + parameters["mu_I"])*I - pred_I*P                                       #Infected snails
    
    dWtdt= parameters["lambda"]*C - (parameters["mu_W"]+parameters["mu_H"])*Wt  #mean worm burden in treated human population
    dWudt= parameters["lambda"]*C - (parameters["mu_W"]+parameters["mu_H"])*Wu  #mean worm burden in untreated human population

    dPdt= parameters["f_P"]*(1-(P/parameters["K_P"]))*P - parameters["mu_P"]*P  #prawn population (number individuals)

    
    return(list(c(dSdt, dEdt, dIdt, dWtdt, dWudt, dPdt)))
}

# mda model code to investigate influence of predators, mda and their interactions on schisto transmission
mda_pred_free_mod = function(t, n, parameters) {
  
    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]

#Dynamic variables
    N = S+E+I  #Total number of snails
    
    W = parameters["cvrg"]*Wt + (1-parameters["cvrg"])*Wu
    
    phi = phi_Wk(W, parameters["kappa"])           

#Schistosome larval concentration equations
    M = 0.5*W*parameters["H"]*phi*parameters["m"]*parameters["v"]*parameters["pi_M"] #- (mu_M + mu_Mq)*M - beta*M*N 
    C = parameters["theta"]*I*parameters["pi_C"] #- (mu_C + mu_Cq)*C - lambda*H*C
    
    #differential equations
    
    dSdt= parameters["f_N"]*(1-(N/parameters["K_N"]))*(S+E) - parameters["mu_N"]*S - parameters["beta"]*M*S  #Susceptible snails
    
    dEdt= parameters["beta"]*M*S - parameters["mu_N"]*E - parameters["sigma"]*E                                     #Exposed snails
    
    dIdt= parameters["sigma"]*E - (parameters["mu_N"] + parameters["mu_I"])*I                                    #Infected snails
    
    dWtdt= parameters["lambda"]*C - (parameters["mu_W"]+parameters["mu_H"])*Wt  #mean worm burden in treated human population
    dWudt= parameters["lambda"]*C - (parameters["mu_W"]+parameters["mu_H"])*Wu  #mean worm burden in untreated human population

    
    return(list(c(dSdt, dEdt, dIdt, dWtdt, dWudt)))
}
