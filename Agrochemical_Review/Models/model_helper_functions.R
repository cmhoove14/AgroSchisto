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

#Function to estimate prevalence given mean worm burden and clumping parameter assuming neg. binomial dist'n
get_prev <- function(clump, 
                     burden){
  pnbinom(1, size = clump, mu = burden, lower.tail = FALSE)
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


#Function to estimate schisto DALYs
est_dalys <- function(burden,     # mean worm burden of NB
                      clump,      # clumping parameter of NB
                      weights_lo, # disability weight for low infection intensity 
                      weights_hi, # disability weight for high infection intensity >=50 eggs/10mL
                      epmL,       # estimate of eggs/10mL per mated female worm
                      pop,        # human population
                      daily = TRUE){  # estimating DALYs on a daily basis?
  eggs_mL <- rnbinom(pop, size = clump, mu = burden) * phi_Wk(burden, clump) * 0.5 * epmL # estimate of egg burden converted from worm burden
  
  n_lo <- length(eggs_mL[eggs_mL > 0 & eggs_mL < 50])
  n_hi <- length(eggs_mL[eggs_mL >= 50])
  
  if(daily){
    
    dalys_lo <- n_lo * weights_lo/365
    dalys_hi <- n_hi * weights_hi/365
    
  } else {
    
    dalys_lo <- n_lo * weights_lo
    dalys_hi <- n_hi * weights_hi

  }
  
  return(dalys_hi + dalys_lo)
}

#Function to return summary of a vector as median (IQR)
get_sum <- function(vec){
  paste0(round(median(vec, na.rm = T), 2), 
         " (", round(quantile(vec, 0.25, na.rm = T), 2), 
         " - ", round(quantile(vec, 0.75, na.rm = T), 2), ")"  )
} 
