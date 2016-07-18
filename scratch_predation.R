Ws = c(1:100)
k = 0.2
m.f = as.numeric()

fx<-function(x, mean.worm = W, clump = k){
  (1-cos(x))/((1+((mean.worm/(clump+mean.worm))*cos(x)))^(1+clump))
}

for(i in 1:length(Ws)){
  W = Ws[i]
  m.f[i] = 1 - ((1-(W/(W+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value)
}

plot(Ws, m.f, type='l', lwd = 2)

#Pred functional response ####################
dens = seq(0,10000,1)

psis1 = as.numeric()
psis0.5 = as.numeric()
psis1.5 = as.numeric()
psis2 = as.numeric()
psis2.5 = as.numeric()
psis3 = as.numeric()
psis4 = as.numeric()

alpha = 0.003
Th.use = 0.067

#get.psi = function(N, n){
#  psi = (alpha*(N^n)) / (1+alpha*Th.use*(N^n))
#  psi
#}

#get.psi(N = 5000/200, n =2)

for(i in 1:length(dens)){
  N = dens[i]
  psis1[i] = (alpha*(N^1)) / (1+alpha*Th.use*(N^1))
}

for(i in 1:length(dens)){
  N = dens[i]
  psis0.5[i] = (alpha*((N/200)^0.5)) / (1+alpha*Th.use*((N/200)^0.5))
}

for(i in 1:length(dens)){
  N = dens[i]
  psis1.5[i] = (alpha*((N/200)^1.5)) / (1+alpha*Th.use*((N/200)^1.5))
}

for(i in 1:length(dens)){
  N = dens[i]
  psis2[i] = (alpha*((N/200)^2)) / (1+alpha*Th.use*((N/200)^2))
}

for(i in 1:length(dens)){
  N = dens[i]
  psis2.5[i] = (alpha*((N/200)^2.5)) / (1+alpha*Th.use*((N/200)^2.5))
}

for(i in 1:length(dens)){
  N = dens[i]
  psis3[i] = (alpha*((N/200)^3)) / (1+alpha*Th.use*((N/200)^3))
}

for(i in 1:length(dens)){
  N = dens[i]
  psis4[i] = (alpha*((N/200)^4)) / (1+alpha*Th.use*((N/200)^4))
}


plot(dens, psis1, type='l', lwd=2, ylim = c(0, 15), xlab = 'snail density (N/m^2)', ylab = 'predation rate')
  lines(dens, psis0.5, lwd = 2, col = 'orange')
  lines(dens, psis1.5, lwd = 2, col = 'blue')
  lines(dens, psis2, lwd = 2, col = 'green')
  lines(dens, psis2.5, lwd = 2, col = 'red')
  lines(dens, psis3, lwd = 2, col = 'purple')
  lines(dens, psis4, lwd = 2, col = 'pink')
  legend('topleft', lwd = 2, col = c('black', 'orange', 'blue', 'green', 'red', 'purple', 'pink'), 
         legend = c(1,0.5,1.5,2,2.5,3,4), title = 'n =', cex=0.6)
  
plot(dens, psis1, type='l', lwd=2, ylim = c(0,alpha+alpha/2), #xlim = c(0,1.5),
     xlab = 'snail density (N/m^2)', ylab = 'predation rate')
  lines(dens, psis0.5, lwd = 2, col = 'orange')
  lines(dens, psis1.5, lwd = 2, col = 'blue')
  lines(dens, psis2, lwd = 2, col = 'green')
  lines(dens, psis2.5, lwd = 2, col = 'red')
  lines(dens, psis3, lwd = 2, col = 'purple')
  lines(dens, psis4, lwd = 2, col = 'pink')
  legend('topleft', lwd = 2, col = c('black', 'orange', 'blue', 'green', 'red', 'purple', 'pink'), 
         legend = c(1,0.5,1.5,2,2.5,3,4), title = 'n =', cex=0.6) 
  

#Predation function testing############

f_N<-parameters["f_N"]
phi_N<-parameters["phi_N"]
z<-parameters["z"]
mu_N<-parameters["mu_N"]
sigma<-parameters["sigma"]
mu_I<-parameters["mu_I"]
alpha<-parameters["alpha"]
Th<-parameters["Th"]
nn<-parameters['nn']
f_P<-parameters["f_P"]
phi_P<-parameters["phi_P"]
mu_P<-parameters["mu_P"]
mu_W<-parameters["mu_W"]
m<-parameters["m"]
H<-parameters["H"]
mu_H<-parameters["mu_H"]

eq_N<-function(y){
  
  f_N*(1-y/phi_N) - 
    mu_N - 
    P_eq*(alpha*((y)^(nn-1))/(1+alpha*Th*((y)^nn)))
  
} 
nn=2
alpha=0.003
P_eq = 60
curve(eq_N, -10000, 10000)
  abline(h = 0, lty=2)
max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))

tests<-data.frame(muPq = seq(0, p.dead, p.dead/100),
                  P_eq = 0,
                  N_eq4 = 0,
                  N_eq3 = 0,
                  N_eq2 = 0,
                  N_eq2.5 = 0,
                  N_eq1.5 = 0,
                  N_eq1 = 0)

for(i in 1:nrow(tests)){
  tests[i,2] = (1-((tests[i,1]+mu_P)/f_P))*phi_P
  P_eq = tests[i,2]
  nn=4
  tests[i,3] = max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))
  nn=3
  tests[i,4] = max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))
  nn=2
  tests[i,5] = max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))
  nn=2.5
  tests[i,6] = max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))
  nn=1.5
  tests[i,7] = max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))
  nn=1
  tests[i,8] = max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))
  if(tests[i,8] < 0){
    tests[i,8] = 0
  }
}

plot(tests$P_eq, tests$N_eq1, type = 'l', ylim = c(0,max(tests$N_eq1)), lwd = 2,
     xlab = 'Equilibrium predator population (P*)', ylab = 'Equilibrium snail population (N*)')
  lines(tests$P_eq, tests$N_eq3, col = 'purple', lwd = 2, lty=3)
  lines(tests$P_eq, tests$N_eq2, col = 'green', lwd=2, lty = 3)
  lines(tests$P_eq, tests$N_eq4, col = 'pink', lwd=2, lty = 3)
  lines(tests$P_eq, tests$N_eq2.5, col = 'red', lwd=2, lty = 3)
  lines(tests$P_eq, tests$N_eq1.5, col = 'blue', lwd=2, lty = 3)
  legend('topright', lwd = 2, col = c('black', 'blue', 'green', 'red', 'purple', 'pink'), 
         legend = c(1,1.5,2,2.5,3,4), title = 'n =', cex=0.8) 
  
plot(tests$P_eq, tests$N_eq1, type = 'l', ylim = c(0,250), lwd = 2, 
       xlab = 'Equilibrium predator population (P*)', ylab = 'Equilibrium snail population (N*)')
  lines(tests$P_eq, tests$N_eq3, col = 'purple', lwd = 2, lty=3)
  lines(tests$P_eq, tests$N_eq2, col = 'green', lwd=2, lty = 3)
  lines(tests$P_eq, tests$N_eq4, col = 'pink', lwd=2, lty = 3)
  lines(tests$P_eq, tests$N_eq2.5, col = 'red', lwd=2, lty = 3)
  lines(tests$P_eq, tests$N_eq1.5, col = 'blue', lwd=2, lty = 3)
  legend('topright', lwd = 2, col = c('black', 'blue', 'green', 'red', 'purple', 'pink'), 
         legend = c(1,1.5,2,2.5,3,4), title = 'n =', cex=0.8)   

#R0 testing ###############  
  
  get_Ro_beta_lamda_nn<-function(muPq = 0, phi_Nq = 1, beta, lamda, f_Nq = 1, nn=1) #variable parameters to be manipulated
  { 
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
    
    P_eq<-(1-((muPq+mu_P)/f_P))*phi_P #Equilibrium estimate of P given prawn predator parameters
    
    if(P_eq<0){P_eq=0}
    
    #Equilibrium estimate of N given snail parameters
    N_eq = max(uniroot.all(f = function(y){f_N*(1-y/phi_N) - 
                               mu_N - 
                               (P_eq*alpha*(y/500)^(nn-1))/(1+alpha*Th*(y/500)^nn)}, 
                           c(0, as.numeric(phi_N))))
    
    if(N_eq < 0){
      N_eq = 0
    }

    pred<-(alpha*P_eq*(N_eq/200)^(nn-1))/(1+(alpha*(N_eq/200)^nn*Th))#death rate of snails due to predators given equilibrium estimates of P and N
    
    Ro_est <- sqrt((0.5*beta*H*N_eq*lamda*sigma)/((mu_W+mu_H)*(mu_N+0.3*pred+sigma)*(mu_N+0.3*pred+mu_I)))
    
    return(c(N_eq,P_eq,pred,Ro_est ))
    
  } #End R0 function  
  
  get_Ro_beta_lamda_nn(muPq = 0.08,lamda = lamda.use, beta = beta.use, nn=1)
  
  p.dead = parameters["f_P"] - parameters["mu_P"] #muPq value at which predator population is eliminated 
  
tests2<-data.frame(muPq = seq(0, p.dead, p.dead/100),
                    P_eq = 0,
                    Ro4 = 0,
                    Ro3 = 0,
                    Ro2.5 = 0,
                    Ro2 = 0,
                    Ro1.5 = 0,
                    Ro1 = 0)  

for(i in 1:nrow(tests2)){
  tests2[i,2] = (1-((tests2[i,1]+mu_P)/f_P))*phi_P
  tests2[i,3] = get_Ro_beta_lamda_nn(muPq = tests2[i,1], phi_Nq = 1, beta=beta.use, lamda=lamda.use, f_Nq = 1, nn=4)[4]
  tests2[i,4] = get_Ro_beta_lamda_nn(muPq = tests2[i,1], phi_Nq = 1, beta=beta.use, lamda=lamda.use, f_Nq = 1, nn=3)[4]
  tests2[i,5] = get_Ro_beta_lamda_nn(muPq = tests2[i,1], phi_Nq = 1, beta=beta.use, lamda=lamda.use, f_Nq = 1, nn=2.5)[4]
  tests2[i,6] = get_Ro_beta_lamda_nn(muPq = tests2[i,1], phi_Nq = 1, beta=beta.use, lamda=lamda.use, f_Nq = 1, nn=2)[4]
  tests2[i,7] = get_Ro_beta_lamda_nn(muPq = tests2[i,1], phi_Nq = 1, beta=beta.use, lamda=lamda.use, f_Nq = 1, nn=1.5)[4]
  tests2[i,8] = get_Ro_beta_lamda_nn(muPq = tests2[i,1], phi_Nq = 1, beta=beta.use, lamda=lamda.use, f_Nq = 1, nn=1)[4]
}

plot(tests2$P_eq, tests2$Ro1, type = 'l', ylim = c(0,max(tests2$Ro1)), lwd = 2,
     xlab = 'Equilibrium predator population (P*)', ylab = 'R0 estimate')
  lines(tests2$P_eq, tests2$Ro3, col = 'purple', lwd = 2, lty=3)
  lines(tests2$P_eq, tests2$Ro2, col = 'green', lwd=2, lty = 3)
  lines(tests2$P_eq, tests2$Ro2.5, col = 'red', lwd=2, lty = 3)
  lines(tests2$P_eq, tests2$Ro1.5, col = 'blue', lwd=2, lty = 3)
  lines(tests2$P_eq, tests2$Ro4, col = 'pink', lwd=2, lty = 3)
  legend('topright', lwd = 2, col = c('black', 'blue', 'green', 'red', 'purple', 'pink'), 
         legend = c(1,1.5,2,2.5,3,4), title = 'n =', cex=0.8)  
  
plot(tests2$P_eq, tests2$Ro1, type = 'l', ylim = c(0,0.75), lwd = 2,
       xlab = 'Equilibrium predator population (P*)', ylab = 'R0 estimate')
  lines(tests2$P_eq, tests2$Ro3, col = 'purple', lwd = 2, lty=3)
  lines(tests2$P_eq, tests2$Ro2, col = 'green', lwd=2, lty = 3)
  lines(tests2$P_eq, tests2$Ro2.5, col = 'red', lwd=2, lty = 3)
  lines(tests2$P_eq, tests2$Ro1.5, col = 'blue', lwd=2, lty = 3)
  lines(tests2$P_eq, tests2$Ro4, col = 'pink', lwd=2, lty = 3)
  legend('topright', lwd = 2, col = c('black', 'blue', 'green', 'red', 'purple', 'pink'), 
         legend = c(1,1.5,2,2.5,3,4), title = 'n =', cex=0.8)   
  
  
  
  
  
  
  
#With prey refuge #################
  mm = 0.1 #proportion of snails that are protected by refugia
  
#Pred functional response ####################
  dens = seq(0,50,0.1)
  
  psis1 = as.numeric()
  psis0.5 = as.numeric()
  psis1.5 = as.numeric()
  psis2 = as.numeric()
  psis2.5 = as.numeric()
  psis3 = as.numeric()
  psis4 = as.numeric()
  
  alpha = 0.003
  Th.use = 0.067
  
  #get.psi = function(N, n){
  #  psi = (alpha*(N^n)) / (1+alpha*Th.use*(N^n))
  #  psi
  #}
  
  #get.psi(N = 5000/200, n =2)
  
  for(i in 1:length(dens)){
    N = dens[i]
    psis1[i] = (alpha*((1-mm)^1*N^1)) / (1+alpha*Th.use*((1-mm)^1*N^1)^1)
  }
  
  for(i in 1:length(dens)){
    N = dens[i]
    psis0.5[i] = (alpha*((1-mm)^0.5*N^0.5)) / (1+alpha*Th.use*((1-mm)^0.5*N^0.5))
  }
  
  for(i in 1:length(dens)){
    N = dens[i]
    psis1.5[i] = (alpha*((1-mm)^1.5*N^1.5)) / (1+alpha*Th.use*((1-mm)^1.5*N^1.5))
  }
  
  for(i in 1:length(dens)){
    N = dens[i]
    psis2[i] = (alpha*((1-mm)^2*N^2)) / (1+alpha*Th.use*((1-mm)^2*N^2))
  }
  
  for(i in 1:length(dens)){
    N = dens[i]
    psis2.5[i] = (alpha*((1-mm)^2.5*N^2.5)) / (1+alpha*Th.use*((1-mm)^2.5*N^2.5))
  }
  
  for(i in 1:length(dens)){
    N = dens[i]
    psis3[i] = (alpha*((1-mm)^3*N^3)) / (1+alpha*Th.use*((1-mm)^3*N^3))
  }
  
  for(i in 1:length(dens)){
    N = dens[i]
    psis4[i] = (alpha*((1-mm)^4*N^4)) / (1+alpha*Th.use*((1-mm)^4*N^4))
  }
  
  
  plot(dens, psis1, type='l', lwd=2, ylim = c(0, 15), xlab = 'snail pop (N)', ylab = 'predation rate')
  lines(dens, psis0.5, lwd = 2, col = 'orange')
  lines(dens, psis1.5, lwd = 2, col = 'blue')
  lines(dens, psis2, lwd = 2, col = 'green')
  lines(dens, psis2.5, lwd = 2, col = 'red')
  lines(dens, psis3, lwd = 2, col = 'purple')
  lines(dens, psis4, lwd = 2, col = 'pink')
  legend('topleft', lwd = 2, col = c('black', 'orange', 'blue', 'green', 'red', 'purple', 'pink'), 
         legend = c(1,0.5,1.5,2,2.5,3,4), title = 'n =', cex=0.8)
  
  plot(dens, psis1, type='l', lwd=2, ylim = c(0,alpha+alpha/2), xlim = c(0,1.5),
       xlab = 'snail pop (N)', ylab = 'predation rate')
  lines(dens, psis0.5, lwd = 2, col = 'orange')
  lines(dens, psis1.5, lwd = 2, col = 'blue')
  lines(dens, psis2, lwd = 2, col = 'green')
  lines(dens, psis2.5, lwd = 2, col = 'red')
  lines(dens, psis3, lwd = 2, col = 'purple')
  lines(dens, psis4, lwd = 2, col = 'pink')
  legend('topleft', lwd = 2, col = c('black', 'orange', 'blue', 'green', 'red', 'purple', 'pink'), 
         legend = c(1,0.5,1.5,2,2.5,3,4), title = 'n =', cex=0.8) 
  
  
#Predation function testing############
  
  f_N<-parameters["f_N"]
  phi_N<-parameters["phi_N"]
  z<-parameters["z"]
  mu_N<-parameters["mu_N"]
  sigma<-parameters["sigma"]
  mu_I<-parameters["mu_I"]
  alpha<-parameters["alpha"]
  Th<-parameters["Th"]
  nn<-parameters['nn']
  f_P<-parameters["f_P"]
  phi_P<-parameters["phi_P"]
  mu_P<-parameters["mu_P"]
  mu_W<-parameters["mu_W"]
  m<-parameters["m"]
  H<-parameters["H"]
  mu_H<-parameters["mu_H"]
  
  eq_N<-function(y){
    
    f_N*(1-y/phi_N) - 
      mu_N - 
      P_eq*((alpha*(1-mm)^(nn-1)*(y/200)^(nn-1))/(1+alpha*Th*(1-mm)^nn*(y/200)^nn))
    
  } 
  
  tests<-data.frame(muPq = seq(0, p.dead, p.dead/100),
                    P_eq = 0,
                    N_eq4 = 0,
                    N_eq3 = 0,
                    N_eq2 = 0,
                    N_eq2.5 = 0,
                    N_eq1.5 = 0,
                    N_eq1 = 0)
  
  for(i in 1:nrow(tests)){
    tests[i,2] = (1-((tests[i,1]+mu_P)/f_P))*phi_P
    P_eq = tests[i,2]
    nn=4
    tests[i,3] = max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))
    nn=3
    tests[i,4] = max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))
    nn=2
    tests[i,5] = max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))
    nn=2.5
    tests[i,6] = max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))
    nn=1.5
    tests[i,7] = max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))
    nn=1
    tests[i,8] = max(uniroot.all(eq_N, c(0, as.numeric(phi_N))))
    if(tests[i,8] < 0){
      tests[i,8] = 0
    }
  }
  
  plot(tests$P_eq, tests$N_eq1, type = 'l', ylim = c(0,max(tests$N_eq1)), lwd = 2,
       xlab = 'Equilibrium predator population (P*)', ylab = 'Equilibrium snail population (N*)')
  lines(tests$P_eq, tests$N_eq3, col = 'purple', lwd = 2, lty=3)
  lines(tests$P_eq, tests$N_eq2, col = 'green', lwd=2, lty = 3)
  lines(tests$P_eq, tests$N_eq4, col = 'pink', lwd=2, lty = 3)
  lines(tests$P_eq, tests$N_eq2.5, col = 'red', lwd=2, lty = 3)
  lines(tests$P_eq, tests$N_eq1.5, col = 'blue', lwd=2, lty = 3)
  legend('topright', lwd = 2, col = c('black', 'blue', 'green', 'red', 'purple', 'pink'), 
         legend = c(1,1.5,2,2.5,3,4), title = 'n =', cex=0.8) 
  
  plot(tests$P_eq, tests$N_eq1, type = 'l', ylim = c(0,250), lwd = 2, 
       xlab = 'Equilibrium predator population (P*)', ylab = 'Equilibrium snail population (N*)')
  lines(tests$P_eq, tests$N_eq3, col = 'purple', lwd = 2, lty=3)
  lines(tests$P_eq, tests$N_eq2, col = 'green', lwd=2, lty = 3)
  lines(tests$P_eq, tests$N_eq4, col = 'pink', lwd=2, lty = 3)
  lines(tests$P_eq, tests$N_eq2.5, col = 'red', lwd=2, lty = 3)
  lines(tests$P_eq, tests$N_eq1.5, col = 'blue', lwd=2, lty = 3)
  legend('topright', lwd = 2, col = c('black', 'blue', 'green', 'red', 'purple', 'pink'), 
         legend = c(1,1.5,2,2.5,3,4), title = 'n =', cex=0.8)   
  
  #R0 testing ###############  
  
  get_Ro_beta_lamda_nn<-function(muPq = 0, phi_Nq = 1, beta, lamda, f_Nq = 1, nn=1) #variable parameters to be manipulated
  { 
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
    
    P_eq<-(1-((muPq+mu_P)/f_P))*phi_P #Equilibrium estimate of P given prawn predator parameters
    
    if(P_eq<0){P_eq=0}
    
    #Equilibrium estimate of N given snail parameters
    N_eq = max(uniroot.all(f = function(y){f_N*(1-y/phi_N) - 
        mu_N - 
        (P_eq*alpha*(1-mm)^(nn-1)*(y/200)^(nn-1))/(1+alpha*Th*(1-mm)^nn*(y/200)^nn)}, 
        c(0, as.numeric(phi_N))))
    
    if(N_eq < 0){
      N_eq = 0
    }
    
    pred<-(alpha*P_eq*(1-mm)^(nn-1)*(N_eq/200)^(nn-1))/(1+(alpha*(1-mm)^nn*(N_eq/200)^nn*Th))#death rate of snails due to predators given equilibrium estimates of P and N
    
    Ro_est <- sqrt((0.5*beta*H*N_eq*lamda*sigma)/((mu_W+mu_H)*(mu_N+0.3*pred+sigma)*(mu_N+0.3*pred+mu_I)))
    
    return(c(N_eq,P_eq,pred,Ro_est ))
    
  } #End R0 function  
  
  get_Ro_beta_lamda_nn(muPq = 0.08,lamda = lamda.use, beta = beta.use, nn=1)
  
  p.dead = parameters["f_P"] - parameters["mu_P"] #muPq value at which predator population is eliminated 
  
  tests2<-data.frame(muPq = seq(0, p.dead, p.dead/100),
                     P_eq = 0,
                     Ro4 = 0,
                     Ro3 = 0,
                     Ro2.5 = 0,
                     Ro2 = 0,
                     Ro1.5 = 0,
                     Ro1 = 0)  
  
  for(i in 1:nrow(tests2)){
    tests2[i,2] = (1-((tests2[i,1]+mu_P)/f_P))*phi_P
    tests2[i,3] = get_Ro_beta_lamda_nn(muPq = tests2[i,1], phi_Nq = 1, beta=beta.use, lamda=lamda.use, f_Nq = 1, nn=4)[4]
    tests2[i,4] = get_Ro_beta_lamda_nn(muPq = tests2[i,1], phi_Nq = 1, beta=beta.use, lamda=lamda.use, f_Nq = 1, nn=3)[4]
    tests2[i,5] = get_Ro_beta_lamda_nn(muPq = tests2[i,1], phi_Nq = 1, beta=beta.use, lamda=lamda.use, f_Nq = 1, nn=2.5)[4]
    tests2[i,6] = get_Ro_beta_lamda_nn(muPq = tests2[i,1], phi_Nq = 1, beta=beta.use, lamda=lamda.use, f_Nq = 1, nn=2)[4]
    tests2[i,7] = get_Ro_beta_lamda_nn(muPq = tests2[i,1], phi_Nq = 1, beta=beta.use, lamda=lamda.use, f_Nq = 1, nn=1.5)[4]
    tests2[i,8] = get_Ro_beta_lamda_nn(muPq = tests2[i,1], phi_Nq = 1, beta=beta.use, lamda=lamda.use, f_Nq = 1, nn=1)[4]
  }
  
  plot(tests2$P_eq, tests2$Ro1, type = 'l', ylim = c(0,max(tests2$Ro1)), lwd = 2,
       xlab = 'Equilibrium predator population (P*)', ylab = 'R0 estimate')
  lines(tests2$P_eq, tests2$Ro3, col = 'purple', lwd = 2, lty=3)
  lines(tests2$P_eq, tests2$Ro2, col = 'green', lwd=2, lty = 3)
  lines(tests2$P_eq, tests2$Ro2.5, col = 'red', lwd=2, lty = 3)
  lines(tests2$P_eq, tests2$Ro1.5, col = 'blue', lwd=2, lty = 3)
  lines(tests2$P_eq, tests2$Ro4, col = 'pink', lwd=2, lty = 3)
  legend('topright', lwd = 2, col = c('black', 'blue', 'green', 'red', 'purple', 'pink'), 
         legend = c(1,1.5,2,2.5,3,4), title = 'n =', cex=0.8)  
  
  plot(tests2$P_eq, tests2$Ro1, type = 'l', ylim = c(0,0.75), lwd = 2,
       xlab = 'Equilibrium predator population (P*)', ylab = 'R0 estimate')
  lines(tests2$P_eq, tests2$Ro3, col = 'purple', lwd = 2, lty=3)
  lines(tests2$P_eq, tests2$Ro2, col = 'green', lwd=2, lty = 3)
  lines(tests2$P_eq, tests2$Ro2.5, col = 'red', lwd=2, lty = 3)
  lines(tests2$P_eq, tests2$Ro1.5, col = 'blue', lwd=2, lty = 3)
  lines(tests2$P_eq, tests2$Ro4, col = 'pink', lwd=2, lty = 3)
  legend('topright', lwd = 2, col = c('black', 'blue', 'green', 'red', 'purple', 'pink'), 
         legend = c(1,1.5,2,2.5,3,4), title = 'n =', cex=0.8)   
  
  