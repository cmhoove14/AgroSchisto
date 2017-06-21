require(deSolve)

#Pop dynamics ##########
prawns = function(t, n, parameters){
  with(as.list(parameters),{
    
    mu_Pq = mu_P #*f(Q)
    
    P=n[1]
    
    dPdt= f_P*(1-P/(phi_P*A))*P - mu_Pq*P
    
    return(list(dPdt))
  })
}

area = 200
nstrt = c(P=1)
time = seq(0, 365*5, 1)

prpar = c(f_P = 1,
          phi_P = 0.125,
          mu_P = 0.038,
          A = area)

proutput = as.data.frame(ode(nstrt, time, prawns, prpar))

plot(x = proutput$time, y = proutput$P, lwd=2, xlab = 'time', ylab = 'pred pop', 
     type = 'l', ylim = c(0, max(proutput$P+5)))

nstrt2 = c(P = max(proutput$P))

#Halt reproduction and see how long population persists
prpar2 = c(f_P = 0,
           phi_P = 0.1,
           mu_P = 0.038,
           A = area)

proutput2 = as.data.frame(ode(nstrt2, time, prawns, prpar2))

plot(x = proutput2$time, y = proutput2$P, lwd=2, xlab = 'time', ylab = 'pred pop', 
     type = 'l', ylim = c(0, max(proutput2$P+5)))

# Calculate PDF f(t) = -dP(t)/dt (i.e. numerically approximate derivative of survival curve)
vect = proutput2$P
  deriv = numeric(length(vect)-1)
for (i in 1:(length(vect)-1)) {
  deriv[i] = vect[i] - vect[i+1]
}

deriv2 = deriv/proutput2$P[1]
deriv2 = deriv2[-(length(deriv2))]
sum(deriv2) 

x = 0:(length(deriv2)-1)

# Calculate life expectancy
sum(x*deriv2)   #Estimate mean of exponentially distrubted lifespans

#Check stability of the model in relation to pred mortality
mup.seq = seq(1e-3, 0.3, length.out = 10)
mup.arr = array(data = NA, dim = c(length(mup.seq), length(time), 2))

for(i in 1:length(mup.seq)){
  prpar['mu_P'] = mup.seq[i]
  mup.arr[i, , ] = ode(nstrt, time, prawns, prpar)
}

plot(mup.arr[1, ,1], mup.arr[1, ,2], type = 'l', xlab = 'time', ylab = 'P',
     ylim = c(0, max(mup.arr[, , 2])+2), xlim = c(0,2500))
  for(j in 2:length(mup.seq)){
    lines(mup.arr[j, ,1], mup.arr[j, ,2], col = j)
  } #population persists up to mortality rates of ~0.2
    legend('topright', legend = round(mup.seq,3), lty = 1, col = c(1:length(mup.seq)),
           bty = 'n', cex = 0.6, title = 'Pred Mortality Rate')
    
#Check out predator function in relation to handling time ############
pred.dy = function(phi_P = 0.1, muP = 0.038, f_P = 0.82,
                   f_N = 0.6, phi_N = 50, mu_N = 0.03,
                   alpha = 0.6, Th = 0.2, area = 200){
  #Equilibrium estimate of P given prawn predator parameters and q
  P.eq = phi_P*area*(1 - muP/f_P)         
  attack = alpha/area/30
  
  #Equilibrium estimate of N given snail parameters
  N.eq = max(uniroot.all(f = function(N){(f_N)*(1 - N/(phi_N*area)) - mu_N - (P.eq*attack*(N/area))/(1+attack*Th*(N/area))},
                         c(0, as.numeric(phi_N*area))))
  
  #Equilibrium predation rate estimate
  psi.eq = (P.eq*attack*(N.eq/area))/(1+attack*Th*(N.eq/area))
  
  return(c( N.eq, P.eq, psi.eq, psi.eq/P.eq))
} 
  pred.dy()
  
  #Th relation to N.eq  
  plot(seq(0,50,0.1), sapply(seq(0,50,0.1), pred.dy, phi_P = 0.1, muP = 0.038, f_P = 0.82,
                             f_N = 0.6, phi_N = 50, mu_N = 0.03,
                             alpha = 0.6, area = 200)[1,], 
       type = 'l', xlab = 'Th', ylab = 'N.eq')
  #alpha relation to N.eq
  plot(seq(0,10,0.1), sapply(seq(0,10,0.1), pred.dy, phi_P = 0.1, muP = 0.038, f_P = 0.82,
                             f_N = 0.6, phi_N = 50, mu_N = 0.03,
                             Th = 0.2, area = 200)[1,], 
       type = 'l', xlab = 'alpha', ylab = 'N.eq')
  #phi_N relation to N.eq
  plot(seq(10,100,1), sapply(seq(10,100,1), pred.dy, phi_P = 0.1, muP = 0.038, f_P = 0.82,
                             f_N = 0.6, mu_N = 0.03,
                             alpha = 0.6, Th = 0.2, area = 200)[1,], 
       type = 'l', xlab = 'phi_N', ylab = 'N.eq')
    