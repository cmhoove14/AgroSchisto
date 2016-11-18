source('Review_models/agroReview_mod1.1.R')

#NEED TO FINISH MERGING R0q code and older get_Ro_beta_lamda code

R0q<-function(phi_Nq = parameters['phi_N'],
              muNq = parameters['mu_N'],
              muPq = parameters['mu_P'],
              f_Nq = parameters['f_N'],
              psi_q = 1){
  
  P_eq<-(1-((muPq)/f_P))*phi_P #Equilibrium estimate of P given prawn predator parameters
    if(P_eq<0){P_eq=0}
  
  
}

get_Ro_beta_lamda<-function(muPq = 0, phi_Nq = 1, beta, lamda, f_Nq = 1) #variable parameters to be sampled
{ 
  f_N<-parameters["f_N"]
  phi_N<-parameters["phi_N"]
  z<-parameters["z"]
  mu_N<-parameters["mu_N"]
  sigma<-parameters["sigma"]
  mu_I<-parameters["mu_I"]
  alpha<-parameters["alpha"]
  nn<-parameters["nn"]
  Th<-parameters["Th"]
  f_P<-parameters["f_P"]
  phi_P<-parameters["phi_P"]
  mu_P<-parameters["mu_P"]
  mu_W<-parameters["mu_W"]
  m<-parameters["m"]
  H<-parameters["H"]
  mu_H<-parameters["mu_H"]
  
  P_eq<-(1-((muPq+mu_P)/f_P))*phi_P #Equilibrium estimate of P given prawn predator parameters
  
  if(P_eq<0){P_eq=0}
  
  #Equilibrium estimate of N given snail parameters
  N_eq = max(uniroot.all(f = function(y){(f_N*f_Nq)*(1-y/(phi_N*phi_Nq)) - 
      mu_N - 
      (P_eq*alpha*(y/200)^(nn-1))/(1+alpha*Th*(y/200)^nn)}, 
      c(0, as.numeric(phi_N*phi_Nq))))
  
  if(N_eq < 0){
    N_eq = 0
  }
  
  print(N_eq)
  
  pred<-(alpha*P_eq*(N_eq/200)^(nn-1))/(1+(alpha*(N_eq/200)^nn*Th))#death rate of snails due to predators given equilibrium estimates of P and N
  
  Ro_est <- sqrt((0.5*beta*H*N_eq*lamda*sigma)/((mu_W+mu_H)*(mu_N+(pred/5)+sigma)*(mu_N+(pred/10)+mu_I)))
  
  return(c(N_eq,P_eq,Ro_est ))
  
} #End R0 function  