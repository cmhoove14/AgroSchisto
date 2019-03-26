#Predator consumption rate equation (holling's)
pred_consumption_rate_mod <- function(target, S, E, I, area, alpha, eps, Th, n, q_red){
  q_red * ((alpha*eps)*(target/area)^n) / (1 + sum((alpha*eps)*Th*(S/area)^n,
                                                    (alpha*eps)*Th*(E/area)^n,
                                                    (alpha*eps)*Th*(I/area)^n))
  }

#Equilibrium snail population from parameters
get_N.eq <- function(f_Nq, K_Nq, mu_Nq, P.eq, area, alpha, eps, Th, n, q_red){
  max(uniroot.all(f = function(N){
    f_Nq*(1 - N/K_Nq)*N - mu_Nq*N - (P.eq * q_red * ((alpha*eps)*(N/area)^n) / (1 + (alpha*eps)*Th*(N/area)^n))
    },
    c(0, as.numeric(K_Nq))))
}

#Mating probability given worm burden, W, and clumping parameter, k
phi_Wk <- function(W, k) {
  if(W <= 0){
    val = 1
  } 
  #print(c(W,k))
  else{
    func <- function(x) {
    a <- ( W / (W + k) )
    b <- ((1-a)^(1+k))/(2*pi)
    return(( b*( 1-cos(x) ) / (( 1 + a*cos(x) )^(1+k)) ))
    }
    
  val = integrate(func, 0, 2*pi, subdivisions = 10000,
                  rel.tol = 1e-10, stop.on.error = FALSE)$value

  }  
    return(1-val)
}

#Standard error given vector 
st.er <- function(x) {
  sd(x)/sqrt(length(x))
} #Function to calculate standard error of the mean

#Prevalence from worm burden and clumping parameter 
Prevalence <- function(W, k) {
  p=1 - (1/(1+W/k)^(k))*(1+W/(1+W/k)) # fraction of humans with at least 2 parasites
  return(p)
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
