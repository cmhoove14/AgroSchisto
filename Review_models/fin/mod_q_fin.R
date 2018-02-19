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

source('~/RemaisWork/Schisto/R Codes/ag_schist/Review_models/fin/r0_of_q_fin.R')
fx<-function(x, mean.worm, clump){
  alpha<-(mean.worm)/(clump+mean.worm)
  (1-cos(x)) / ( (1+(alpha*cos(x)))^(1+clump) )
}
mateprob<-function(W, k){
  alpha<-W/(W+k)
  1-( (1-alpha)^(k+1) * (integrate(fx, 0, 2*pi, W, k, stop.on.error = F, subdivisions = 1000)$value) /(2*pi)  )
}


#full model code to investigate influence of agrochemicals on schisto transmission

agrochem_mod = function(t, n, parameters) {
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    P=n[5]
    
#Dynamic variables
    N = S+E+I                                                                       #Total number of snails
    
    gamma = mateprob(W = W, k = k)           
    
    pred_S = (alpha*(S/A))/(1+(alpha*Th*(N/A)^nn))  #per capita predation of predators on susceptible snails
    pred_E = (alpha*(E/A))/(1+(alpha*Th*(N/A)^nn))  #per capita predation of predators on susceptible snails
    pred_I = (alpha*(I/A))/(1+(alpha*Th*(N/A)^nn))  #per capita predation of predators on susceptible snails
    
#agrochemical parameters adjusted
    phinadj = phi_N * phinq
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
    M = 0.5*W*H*gamma*m*vadj*pimadj #- (mu_M + mu_Mq)*M - beta*M*N 
    C = thetaadj*I*picadj #- (mu_C + mu_Cq)*C - lamda*H*C
    
    #differential equations
    
    dSdt= fnadj*(1-(N/phinadj))*(S+E) - mu_N*S - pred_Sadj*P - beta*Om*M*S  #Susceptible snails
    
    dEdt= beta*M*S - mu_N*E - pred_Eadj*P - sigma*E                                     #Exposed snails
    
    dIdt= sigma*E - (mu_N + mu_I)*I - pred_Iadj*P                                       #Infected snails
    
    dWdt= lamda*Om*C - (mu_W+mu_H)*W                                                       #mean worm burden in human population
    
    dPdt= f_P*(1-(P/phi_P))*P - mupadj*P                                      #prawn population (number individuals)

    
    return(list(c(dSdt, dEdt, dIdt, dWdt, dPdt)))
  }) 
}

ac.pars = parameters

ac.pars['phinq'] = 1
ac.pars['fnq'] = 1
ac.pars['munq'] = 0
ac.pars['vq'] = 1
ac.pars['pimq'] = 1
ac.pars['thetaq'] = 1
ac.pars['picq'] = 1
ac.pars['mupq'] = 0
ac.pars['psiq'] = 1

#ac.pars['beta'] = 1.340292e-05/area
#ac.pars['lamda'] = 4.722579e-05/area

ac.start = c(S = 20*area,
             E = 5*area,
             I = 1*area,
             W = 30,
             P = 0.1*area)

ac.time = seq(0, (365*50), 5)

ac.sim = as.data.frame(ode(ac.start, ac.time, agrochem_mod, ac.pars))
  ac.sim$N = ac.sim$S + ac.sim$E + ac.sim$I
  ac.sim$Wf = ac.sim$W * sapply(ac.sim$W, mateprob, k=ac.pars['k']) * 0.5
  
  ac.eqbm = ac.sim[dim(ac.sim)[1],]
    
ggplot(ac.sim, aes(x = time)) +
  theme_bw() +
  geom_line(aes(y = N), size = 1.25, col = 'black') +
  geom_line(aes(y = S), size = 1.25, col = 'green') +
  geom_line(aes(y = E), size = 1.25, col = 'orange') +
  geom_line(aes(y = I), size = 1.25, col = 'red') 

ggplot(ac.sim, aes(x = time)) +
  theme_bw() +
  geom_line(aes(y = P), size = 1.25, col = 'blue') 

ggplot(ac.sim, aes(x = time)) +
  theme_bw() +
  geom_line(aes(y = W), size = 1.25, col = 'purple') 

#Get probability of having high worm burden, with high worm burden defined as 100
  rho = 1 - pnbinom(100, mu = ac.eqbm$Wf, size = ac.pars['k'])
