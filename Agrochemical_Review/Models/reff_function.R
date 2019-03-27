N_star <- function(parameters, W, dd = T){
  with(as.list(parameters),{
    
    P.eq = phi_P*(1 - mu_P/f_P)
    
    if(dd){
      N_est <- max(uniroot.all(f = function(N){  
        f_N*(1-N/phi_N)*(1+(beta*W*H*m*phi_Wk(W, k = k)*v*pi_M)/(2*(mu_N+sigma+(P.eq*alpha*eps)/(A+alpha*eps*Th*N)))) - 
          (mu_N + 0.5*beta*W*H*m*phi_Wk(W, k = k)*v*pi_M + (P.eq*alpha*eps)/(A+alpha*eps*Th*N))
    }, c(0, as.numeric(phi_N))))

    } else {
      N_est <- max(uniroot.all(f = function(N){  
        f_N*(1-N/phi_N)*(1+(beta*W*H*m*v*pi_M)/(2*(mu_N+sigma+(P.eq*alpha*eps)/(A+alpha*eps*Th*N)))) - 
          (mu_N + 0.5*beta*W*H*m*v*pi_M + (P.eq*alpha*eps)/(A+alpha*eps*Th*N))
    }, c(0, as.numeric(phi_N))))

    }
    
  return(N_est)  
  })
}

S_star <- function(parameters, W, dd = T){
  with(as.list(parameters),{
    
    P.eq = phi_P*(1 - mu_P/f_P)
    
  N_use <- N_star(parameters = parameters, W = W, dd = dd)
  
  if(dd){
  S_est = N_use / (1 + 
                      (beta*W*H*m*phi_Wk(W, k = k)*v*pi_M)/(2*(mu_N+sigma+(P.eq*alpha*eps)/(A+alpha*eps*Th*N_use))) + #T2
                      sigma / (mu_N + mu_I + (P.eq*alpha*eps)/(A+alpha*eps*Th*N_use)) * #T1
                      (beta*W*H*m*phi_Wk(W, k = k)*v*pi_M)/(2*(mu_N+sigma+(P.eq*alpha*eps)/(A+alpha*eps*Th*N_use)))) #T2
    
  } else {
  S_est = N_use / (1 + 
                      (beta*W*H*m*v*pi_M)/(2*(mu_N+sigma+(P.eq*alpha*eps)/(A+alpha*eps*Th*N_use))) + #T2
                      sigma / (mu_N + mu_I + (P.eq*alpha*eps)/(A+alpha*eps*Th*N_use)) * #T1
                      (beta*W*H*m*v*pi_M)/(2*(mu_N+sigma+(P.eq*alpha*eps)/(A+alpha*eps*Th*N_use)))) #T2
    
  }
    return(S_est)
  
  })
  
}

Reff <- function(parameters, W, dd = T){
  with(as.list(parameters),{
    
    P.eq = phi_P*(1 - mu_P/f_P)
    
    N_use <- N_star(parameters = parameters, W = W, dd = dd)
    
    if(dd){
    (lambda * theta * pi_C * 
      sigma / (mu_N + mu_I + (P.eq*alpha*eps)/(A+alpha*eps*Th*N_use)) * #T1
      (beta*W*H*m*phi_Wk(W, k = k)*v*pi_M)/(2*(mu_N+sigma+(P.eq*alpha*eps)/(A+alpha*eps*Th*N_use))) * 
      S_star(parameters = parameters, W = W, dd = dd)) / 
      ((mu_W + mu_H)*W)
      
    } else {
    (lambda * theta * pi_C * 
      sigma / (mu_N + mu_I + (P.eq*alpha*eps)/(A+alpha*eps*Th*N_use)) * #T1
      (beta*W*H*m*v*pi_M)/(2*(mu_N+sigma+(P.eq*alpha*eps)/(A+alpha*eps*Th*N_use))) * 
      S_star(parameters = parameters, W = W, dd = dd)) / 
      ((mu_W + mu_H)*W)
      
    }
    
  })
  
}

test_Ws <- c(seq(1e-6, 1e-3, length.out = 100), seq(1e-3, 1, length.out = 300), 1:100)

plot(log(test_Ws),sapply(test_Ws, Reff, parameters = init_pars), type = "l", lwd = 2)
  lines(log(test_Ws),sqrt(sapply(test_Ws, Reff, parameters = init_pars, dd = F)), lty = 2, lwd = 2)
  