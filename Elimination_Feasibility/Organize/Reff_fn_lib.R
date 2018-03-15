##### BASIC Function lib for Reff #####

getReff<-function(parameters, W, k){ #no prawns, no psi=0, Set mu_P<-1
  f_N<-parameters["f_N"]
  phi_N<-parameters["phi_N"]
  z<-parameters["z"]
  mu_N<-parameters["mu_N"]
  sigma<-parameters["sigma"]
  mu_I<-parameters["mu_I"]
  alpha<-parameters["alpha"]
  Th<-parameters["Th"]
  f_P<-parameters["f_P"]
  phi_P<-parameters["phi_P"]
  mu_P<-parameters["mu_P"]
  mu_W<-parameters["mu_W"]
  m<-parameters["m"]
  H<-parameters["H"]
  mu_H<-parameters["mu_H"]
  muPq<-parameters["muPq"]
  phi_Nq<-parameters["phi_Nq"]
  beta<-parameters["beta"] 
  lamda<-parameters["lamda"]
  f_Nq<-parameters["f_Nq"]
  
  k1<-sigma/(mu_N+mu_I)
  k2<-beta*0.5*W*H*phi_Wk(W,k)/(mu_N+sigma)
  
  Num_1<-lamda*k1*k2*phi_N*phi_Nq
  Num_2<- ( (1+k2)*f_N ) - ( mu_N + (beta*0.5*H*W*phi_Wk(W,k) ) )
  Den<-( 1+k2+(k1*k2) ) * ( 1+k2 ) * ( mu_W+mu_H )* W* f_N
  
  Reff<-as.numeric(Num_1*Num_2/Den)
  
  k1_nophi<-sigma/(mu_N+mu_I)
  k2_nophi<-beta*0.5*W*H/(mu_N+sigma)
  
  Num_1_nophi<-lamda*k1_nophi*k2_nophi*phi_N*phi_Nq
  Num_2_nophi<- ( (1+k2_nophi)*f_N ) - ( mu_N + (beta*0.5*H*W ) )
  Den_nophi<-( 1+k2_nophi+(k1_nophi*k2_nophi) ) * ( 1+k2_nophi ) * ( mu_W+mu_H )* W* f_N
  
  Reff_noMatingFn<-as.numeric(Num_1_nophi*Num_2_nophi/Den_nophi)
  
  list(Reff=Reff, Reff_noPhi=Reff_noMatingFn)
}

bouncebackRate<-function(W2, t2, W1, t1){
  
  bbr<-((W2-W1)/(t2-t1))*(1/W1)
  bbr
}

### Obtain the negative binomial distribution for worm burden #####
#The ecological model for the neg bin prob 
# i number of parasites per person
# k clumping parameter
# m mean burden
# refer wikipedia page - https://en.wikipedia.org/wiki/Negative_binomial_distribution
# formula refer - http://influentialpoints.com/Training/negative_binomial_distribution-principles-properties-assumptions.htm

Prob_negbin<-function(i,k,m){
  p<-(gamma(k+i)/(gamma(i+1) * gamma(k))) * (1 + (m/k))^(-k-i) * ((m/k)^i)
  p
}


Prob_gaussian<-function(y, mu, sd){
  p<-(1/(sqrt(2*pi)*sd) ) * exp( -((y-mu)^2)/(2*sd^2) )
  p
}
####################################################################################
getReff_withP<-function(parameters,W, k, phi_P, N_star){ #with finite number of prawns
  f_N<-parameters["f_N"]
  phi_N<-parameters["phi_N"]
  z<-parameters["z"]
  mu_N<-parameters["mu_N"]
  sigma<-parameters["sigma"]
  mu_I<-parameters["mu_I"]
  alpha<-parameters["alpha"]
  Th<-parameters["Th"]
  f_P<-parameters["f_P"]
  # phi_P<-parameters["phi_P"]
  mu_P<-parameters["mu_P"]
  mu_W<-parameters["mu_W"]
  m<-parameters["m"]
  H<-parameters["H"]
  mu_H<-parameters["mu_H"]
  muPq<-parameters["muPq"]
  phi_Nq<-parameters["phi_Nq"]
  beta<-parameters["beta"] 
  lamda<-parameters["lamda"]
  f_Nq<-parameters["f_Nq"]
  
  
  
  P_star<- (1-  ((mu_P+muPq)/f_P) )*phi_P
  if(P_star<0) P_star<-0
  psi_star<- (alpha*P_star)/(1+(alpha*N_star*Th) )
  
  k1<-sigma/(mu_N+psi_star+mu_I)
  k2<-(beta*0.5*W*H*phi_Wk(W,k))/(mu_N+psi_star+sigma)
  
  Num1<-lamda*k1*k2*phi_N*phi_Nq
  Num2<-f_N*(1+k2) - mu_N - psi_star - (beta*0.5*W*H*phi_Wk(W,k))
  Den<-(mu_W+mu_H)*W*f_N*(1+k2)*(1+k2+(k1*k2))
  
  Reff<-Num1*Num2/Den
  Reff
  
  k1_nophi<-sigma/(mu_N+psi_star+mu_I)
  k2_nophi<-(beta*0.5*W*H)/(mu_N+psi_star+sigma)
  
  Num1_nophi<-lamda*k1_nophi*k2_nophi*phi_N*phi_Nq
  Num2_nophi<-f_N*(1+k2_nophi) - mu_N - psi_star - (beta*0.5*W*H)
  Den_nophi<-(mu_W+mu_H)*W*f_N*(1+k2_nophi)*(1+k2_nophi+(k1_nophi*k2_nophi))
  
  Reff_nophi<-Num1_nophi*Num2_nophi/Den_nophi
  Reff_nophi
  
  list(Reff=Reff, Reff_nophi=Reff_nophi)
}
