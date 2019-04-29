#Null (default) functional responses ##############
  nil1<-function(In){ #For proportional responses
    return(1)
  }
  
  nil0<-function(In){ #For additive responses
    return(0)
  }

r0.Ag = function(In = 0,
                 He = 0,
                 Fe = 0,
                 parameters,
                 
                 f.in.f_Nq = nil1, 
                 f.he.f_Nq = nil1, 
                 f.fe.f_Nq = nil1, 
                 
                 f.in.mu_Pq = nil0,
                 f.he.mu_Pq = nil0,
                 f.fe.mu_Pq = nil0,
                 
                 f.in.K_Nq = nil1, 
                 f.he.K_Nq = nil1, 
                 f.fe.K_Nq = nil1, 
                 
                 f.in.mu_Nq = nil0, 
                 f.he.mu_Nq = nil0, 
                 f.fe.mu_Nq = nil0,
                 
                 f.in.alpha_q = nil1,
                 f.he.alpha_q = nil1,
                 f.fe.alpha_q = nil1,
                 
                 f.in.theta_q = nil1, 
                 f.he.theta_q = nil1, 
                 f.fe.theta_q = nil1, 
                 
                 f.in.pi_Mq = nil1, 
                 f.he.pi_Mq = nil1, 
                 f.fe.pi_Mq = nil1, 
                 
                 f.in.pi_Cq = nil1, 
                 f.he.pi_Cq = nil1, 
                 f.fe.pi_Cq = nil1, 
                 
                 f.in.v_q = nil1,
                 f.he.v_q = nil1,
                 f.fe.v_q = nil1)

{with(as.list(parameters),{
# get agrochemical-influenced parameters based on agrochemical concentrations and assuming purely additive effects _q#######

  #Snail birth rate  
    f_Nq = f_N * (1-sum(1-c(f.in.f_Nq(In), f.he.f_Nq(He), f.fe.f_Nq(Fe))))
      if(f_Nq < 0) f_Nq = 0
    
  #Predator mortality rate    
    muPq = mu_P + (f.in.mu_Pq(In) + f.he.mu_Pq(He) + f.fe.mu_Pq(Fe))
  
  #Snail carrying capacity
    K_Nq = K_N * (1-sum(1-c(f.in.K_Nq(In), f.he.K_Nq(He), f.fe.K_Nq(Fe))))
      if(K_Nq < 0) K_Nq = 0
  
  #Snail mortality rate
    mu_Nq = mu_N + (f.in.mu_Nq(In) + f.he.mu_Nq(He) + f.fe.mu_Nq(Fe))
  
  #Predator consumption rate
    pred_red = 1-sum(1-c(f.in.alpha_q(In), f.he.alpha_q(He), f.fe.alpha_q(Fe)))
      if(pred_red < 0) pred_red = 0
  
  #Cercarial shedding rate
    theta_q = theta * (1-sum(1-c(f.in.theta_q(In), f.he.theta_q(He), f.fe.theta_q(Fe))))
      if(theta_q < 0) theta_q = 0
  
  #Miracidial survival
    pi_Mq = pi_M * (1-sum(1-c(f.in.pi_Mq(In), f.he.pi_Mq(He), f.fe.pi_Mq(Fe))))
      if(pi_Mq < 0) pi_Mq = 0
  
  #Cercarial survival
    pi_Cq = pi_C*(1-sum(1-c(f.in.pi_Cq(In), f.he.pi_Cq(He), f.fe.pi_Cq(Fe))))
      if(pi_Cq < 0) pi_Cq = 0
    
  #Schisto egg viability
    v_q = v * (1-sum(1-c(f.in.v_q(In), f.he.v_q(He), f.fe.v_q(Fe))))
      if(v_q < 0) v_q = 0

#Equilibrium estimates of pred pop, snail pop, and predation rate #########
  P.eq = K_P*(1 - muPq/f_P)         
  
  if(P.eq<0){P.eq = 0}

#Equilibrium estimate of N given snail parameters
  N.eq = get_N.eq(f_Nq, K_Nq, mu_Nq, P.eq, area, alpha, eps, Th, n=nn, q_red=pred_red)
  
  if(N.eq<0){N.eq = 0}

#Equilibrium predation rate estimate
  psi.eq = pred_consumption_rate_mod(target = N.eq, S = N.eq, E = 0, I = 0, 
                                     area, alpha, eps, Th, n=nn, q_red = pred_red)

#R_0 of q estimate ##########
  #Transmission matrix
  T_mat <- matrix(c(0,0,0.5*beta*H*m*v_q*pi_Mq*N.eq,
                    0,0,0,
                    0,lambda*theta_q*pi_Cq,0), nrow = 3)
  
  #transition matrix
  S_mat <- matrix(c(-(mu_Nq+P.eq*((pred_red*alpha*eps)/(area+alpha*eps*Th*N.eq))+sigma),0,0,
                    sigma, -(mu_Nq+P.eq*((pred_red*alpha*eps)/(area+alpha*eps*Th*N.eq))+mu_I),0,
                    0,0,-(mu_H+mu_W)), nrow = 3)
  
  #Next generation matrix with large domain (K_L = T*-S^-1), dominant eigenvalue of which is R0
  K_L <- T_mat %*% -solve(S_mat)
  
  r0 <- max(eigen(K_L)$values)
  
return(c("N_eq" = N.eq, 
         "P_eq" = P.eq, 
         "R0" = r0))

  })
}

r0.Ag.pars <- function(A, H, f_N, K_N, z, mu_N, mu_I,
                       f_P, K_P, mu_P, alpha, eps, Th, nn, 
                       m, v, pi_M, theta, pi_C,
                       beta, sigma, lambda, kappa, 
                       mu_W, mu_H, cvrg, eff, 
                       Knq, fnq, munq, vq, pimq, thetaq, picq, mupq, psiq){
  
# get agrochemical-influenced parameters #######

  #Snail birth rate  
    f_Nq = f_N * fnq
    
  #Predator mortality rate    
    muPq = mu_P + mupq
  
  #Snail carrying capacity
    K_Nq = K_N * Knq
  
  #Snail mortality rate
    mu_Nq = mu_N + munq
  
  #Predator consumption rate
    pred_red = psiq
  
  #Cercarial shedding rate
    theta_q = theta * thetaq
  
  #Miracidial survival
    pi_Mq = pi_M * pimq
  
  #Cercarial survival
    pi_Cq = pi_C * picq
    
  #Schisto egg viability
    v_q = v * vq

#Equilibrium estimates of pred pop, snail pop, and predation rate #########
  P.eq = K_P*(1 - muPq/f_P)         
  
  if(P.eq<0){P.eq = 0}

#Equilibrium estimate of N given snail parameters
  N.eq = get_N.eq(f_Nq, K_Nq, mu_Nq, P.eq, area, alpha, eps, Th, n=nn, q_red=pred_red)
  
  if(N.eq<0){N.eq = 0}

#Equilibrium predation rate estimate
  psi.eq = pred_consumption_rate_mod(target = N.eq, S = N.eq, E = 0, I = 0, 
                                     area, alpha, eps, Th, n=nn, q_red = pred_red)

#R_0 of q estimate ##########
  #Transmission matrix
  T_mat <- matrix(c(0,0,0.5*beta*H*m*v_q*pi_Mq*N.eq,
                    0,0,0,
                    0,lambda*theta_q*pi_Cq,0), nrow = 3)
  
  #transition matrix
  S_mat <- matrix(c(-(mu_Nq+P.eq*((pred_red*alpha*eps)/(area+alpha*eps*Th*N.eq))+sigma),0,0,
                    sigma, -(mu_Nq+P.eq*((pred_red*alpha*eps)/(area+alpha*eps*Th*N.eq))+mu_I),0,
                    0,0,-(mu_H+mu_W)), nrow = 3)
  
  #Next generation matrix with large domain (K_L = T*-S^-1), dominant eigenvalue of which is R0
  K_L <- T_mat %*% -solve(S_mat)
  
  r0 <- max(eigen(K_L)$values)
  
return(data.frame("N_eq" = N.eq, 
                  "P_eq" = P.eq, 
                  "R0" = r0))

}

r0.Base = function(parameters)

{with(as.list(parameters),{
#Equilibrium estimate of N given snail parameters
  N.eq = K_N*(1 - mu_N/f_N)

#R_0 estimate ##########
  #Transmission matrix
  T_mat <- matrix(c(0,0,0.5*beta*H*m*v*pi_M*N.eq,
                    0,0,0,
                    0,lambda*theta*pi_C,0), nrow = 3)
  
  #transition matrix
  S_mat <- matrix(c(-(mu_N+sigma),0,0,
                    sigma, -(mu_N+mu_I),0,
                    0,0,-(mu_H+mu_W)), nrow = 3)
  
  #Next generation matrix with large domain (K_L = T*-S^-1), dominant eigenvalue of which is R0
  K_L <- T_mat %*% -solve(S_mat)
  
  r0 <- max(eigen(K_L)$values)
  
return(c("N_eq" = N.eq,
       "P_eq" = NA,
       "R0_base" = r0))

  })
}

r0.pred = function(parameters)

{with(as.list(parameters),{
#equilibrium estimate of P given parameters
  P.eq = K_P*(1 - mu_P/f_P)         
 
#Equilibrium estimate of N given snail parameters
  N.eq = get_N.eq(f_N, K_N, mu_N, P.eq, area, alpha, eps, Th, n=nn, q_red=1)

#R_0 estimate ##########
  #Transmission matrix
  T_mat <- matrix(c(0,0,0.5*beta*H*m*v*pi_M*N.eq,
                    0,0,0,
                    0,lambda*theta*pi_C,0), nrow = 3)
  
  #transition matrix
  S_mat <- matrix(c(-(mu_N+P.eq*((alpha*eps)/(area+alpha*eps*Th*N.eq))+sigma),0,0,
                    sigma, -(mu_N+P.eq*((alpha*eps)/(area+alpha*eps*Th*N.eq))+mu_I),0,
                    0,0,-(mu_H+mu_W)), nrow = 3)
  
  #Next generation matrix with large domain (K_L = T*-S^-1), dominant eigenvalue of which is R0
  K_L <- T_mat %*% -solve(S_mat)
  
  r0 <- max(eigen(K_L)$values)
  
return(c("N_eq" = N.eq,
         "P_eq" = P.eq,
         "R0_pred" = r0))

  })
}
