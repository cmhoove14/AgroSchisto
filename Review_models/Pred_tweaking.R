require(deSolve)

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

prpar = c(f_P = 0.02,
          phi_P = 0.5,
          mu_P = 0.0026,
          A = area)

proutput = as.data.frame(ode(nstrt, time, prawns, prpar))

plot(x = proutput$time, y = proutput$P, lwd=2, xlab = 'time', ylab = 'pred pop', 
     type = 'l', ylim = c(0, max(proutput$P+5)))

nstrt2 = c(P = max(proutput$P))

prpar2 = c(f_P = 0,
           phi_P = 0.5,
           mu_P = 0.0026,
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

plot(x = x, y = deriv2, type = 'l')

# Calculate life expectancy
sum(x*deriv2)   #Mean life expectancy is about a year, that's what we want