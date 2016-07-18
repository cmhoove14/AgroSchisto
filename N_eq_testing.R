require(ggplot2)

#Parameter values ########################
parameters=c( #Excluding beta, Phi_Nq, f_Nq, and muPq which will be read into the R0 function
  ##standard snail parameters 
  f_N=0.10, # recruitment rate (from sokolow et al)
  phi_N=10000, # carrying capacity from sokolow et al
  z=0.5, #Proportion of exposed snails that reproduce from sokolow et al
  mu_N=1/60, #Mortality rate from Sokolow et al
  sigma=1/40, #Transition rate from exposed to infected from sokolow et al
  mu_I=1/10, #additional snail death due to infection from sokolow et al
  
  #prawn parameters
  alpha=0.003, #attack rate
  Th=0.067,#~Prawn predation limit
  f_P=0.234/2, #prawn birth rate from Cervantes-Santiago Aquaculture 2010 paper (/2 for 1:1 female-male ratio)
  phi_P=120,  #prawn carrying capacity
  mu_P= 0.03809524, #observed daily 24 hr prawn mortality rate in chlorP-free tanks of mesocosm#0.00137, previous rate
  
  #Adult Worm, Miracidia and Circariae Parameters
  lamda=1.5e-5, #probability of snail shedding a cercariae that infects a human host and survives to reproduction
  mu_W=1/(3.3*365), # death rate of adult worms
  #m=0.5, #miracidial shedding rate per reproductive female divided by miracidial mortality; from sokolow et al
  
  #Human parameters
  H=300, #number of humans
  mu_H=1/(60*365) #Assumes 60 year lifespan
)

beta_0= 7.446512e-05
beta_up=9.5484e-05 #NEED TO UPDATE FROM ARATHI< THIS IS JUST A PLACEHOLDER
beta_lo=7.67848e-06

#R0 function #########################

get_Ro<-function(muPq = 0, phi_Nq = 1, beta, f_Nq = 1) #variable parameters to be manipulated
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
  lamda<-parameters["lamda"]
  mu_W<-parameters["mu_W"]
  m<-parameters["m"]
  H<-parameters["H"]
  mu_H<-parameters["mu_H"]
  
  P_eq<-(1-((muPq+mu_P)/f_P))*phi_P #Equilibrium estimate of P given prawn predator parameters
  if(P_eq<0){
    P_eq=0
  }
  #Equilibrium estimate of N given snail parameters
  #Shorthand values to use in N_eq expression
  #a= -(alpha*Th*f_N*f_Nq*phi_N*phi_Nq^-1)
  a = -(alpha*Th*f_N*f_Nq)/(phi_N*phi_Nq)
  #b= f_N*f_Nq*alpha*Th - mu_N*alpha*Th - f_N*f_Nq*phi_N*phi_Nq^-1
  b = (f_N*f_Nq)/(phi_N*phi_Nq) + alpha*Th*f_N*f_Nq - alpha*Th*mu_N
  #c= f_N*f_Nq - mu_N - alpha*P_eq
  c = f_N*f_Nq - alpha*P_eq - mu_N
  
  if((b^2-4*a*c)<0){ #If prawn population sufficient to eliminate snails, N_eq=0
    N_eq1=0
  } else {
    N_eq1 <- (-b + sqrt(b^2-4*a*c)) / (2*a) #Function to solve quadratic expression for N_eq
  }
  
  if((b^2-4*a*c)<0){ #If prawn population sufficient to eliminate snails, N_eq=0
    N_eq2=0
  } else {
    N_eq2 <- (-b - sqrt(b^2-4*a*c)) / (2*a) #Function to solve quadratic expression for N_eq
  }
  
  if(N_eq1 > N_eq2){ #If prawn population sufficient to eliminate snails, N_eq=0
    N_eq=N_eq1
  } else {
    N_eq = N_eq2 #Function to solve quadratic expression for N_eq
  }
  
  
  pred<-(alpha*P_eq)/(1+(alpha*N_eq*Th))#death rate of snails due to predators given equilibrium estimates of P and N
  
  T1<-0.5*beta*H*N_eq
  T2<-lamda*sigma
  T3<- (mu_W+mu_H)*(mu_N+pred+sigma)*(mu_N+pred+mu_I)
  
  Ro_est <- sqrt((T1*T2)/T3)
  
  print(N_eq)
  print(P_eq)
  Ro_est 
  
}  


get_Ro(muPq=0.06, f_Nq =1.3, beta = beta_0) 

mu_Pqs<-seq(from=0, to=0.09, by=0.001) 
R0s.muPs<-rep(0, 91)
  for(i in 1:91){
    R0s.muPs[i]=get_Ro(muPq=mu_Pqs[i], phi_Nq =1, beta = beta_0)
  }

plot(x=mu_Pqs, y=R0s.muPs, xlab = "Prawn mortality rate (muPq)", ylab = "R0 estimate")
  
#Weird testing ##############
neeqs<-seq(0, 12000, by=100)

#N_eq profile when P_eq=83
  ys.83<-rep(0, 121)
  
  for(i in 1:121){
    ys.83[i]=4.425e-9*neeqs[i]^2+2.95e-5*neeqs[i]+0.1475-0.003*83
  }

    plot(x=neeqs, y=ys.83, type='l')

#N_eq profile when P_eq=50
  ys.50<-rep(0, 121)
  
  for(i in 1:121){
    ys.50[i]=4.425e-9*neeqs[i]^2+2.95e-5*neeqs[i]+0.1475-0.003*50
  }

    lines(x=neeqs, y=ys.50, col='red')
    
#N_eq profile when P_eq=25
  ys.25<-rep(0, 121)
    
    for(i in 1:121){
      ys.25[i]=4.425e-9*neeqs[i]^2+2.95e-5*neeqs[i]+0.1475-0.003*25
    }
    
    lines(x=neeqs, y=ys.25, col='blue')

#Inside R0 function to check out equilibrium estimates of snail population from quadratic ##############
get_N.eq<-function(muPq = 0, f_Nq = 1, phi_Nq = 1)
{   #HAVE TO SET muPq and phi_Nq in function call
  f_N<-parameters["f_N"]
  phi_N<-parameters["phi_N"]
  z<-parameters["z"]
  mu_N<-parameters["mu_N"]
  beta<-parameters["beta"]
  sigma<-parameters["sigma"]
  mu_I<-parameters["mu_I"]
  alpha<-parameters["alpha"]
  Th<-parameters["Th"]
  f_P<-parameters["f_P"]
  phi_P<-parameters["phi_P"]
  mu_P<-parameters["mu_P"]
  lamda<-parameters["lamda"]
  mu_W<-parameters["mu_W"]
  m<-parameters["m"]
  H<-parameters["H"]
  mu_H<-parameters["mu_H"]
  
  P_eq<-(1-((muPq+mu_P)/f_P))*phi_P #Equilibrium estimate of P given prawn predator parameters
  if(P_eq<0){
    P_eq=0
  }
  #Equilibrium estimate of N given snail parameters
  #Shorthand values to use in N_eq expression
  #a= -(alpha*Th*f_N*f_Nq*phi_N*phi_Nq^-1)
  a = -(alpha*Th*f_N*f_Nq)/(phi_N*phi_Nq)
  #b= f_N*f_Nq*alpha*Th - mu_N*alpha*Th - f_N*f_Nq*phi_N*phi_Nq^-1
  b = (f_N*f_Nq)/(phi_N*phi_Nq) + alpha*Th*f_N*f_Nq - alpha*Th*mu_N
  #c= f_N*f_Nq - mu_N - alpha*P_eq
  c = f_N*f_Nq - alpha*P_eq - mu_N
  
    if((b^2-4*a*c)<0){ #If prawn population sufficient to eliminate snails, N_eq=0
      N_eq1=0
    } else {
      N_eq1 <- (-b + sqrt(b^2-4*a*c)) / (2*a) #Function to solve quadratic expression for N_eq
    }
  
    if((b^2-4*a*c)<0){ #If prawn population sufficient to eliminate snails, N_eq=0
      N_eq2=0
    } else {
      N_eq2 <- (-b - sqrt(b^2-4*a*c)) / (2*a) #Function to solve quadratic expression for N_eq
    }
  
  print(c(P_eq,N_eq1, N_eq2))
  
}
    
get_N.eq(muPq = 0.06, phi_Nq =1.61)
get_N.eq(muPq = 0.73)


#Check out influence of snail fertility rate scalars across muPq profile ######################
  fnqs<-seq(1, 2, by=0.025) 
  mups<-seq(0, 0.1, by=0.0025)
  
  tests<-data.frame("f_Nq" = rep(fnqs, 41),
                    "muPq" = rep(mups, each=41),
                    "P_eq" = rep(0, length(fnqs)*41),
                    "N_eq1" = rep(0, length(fnqs)*41),
                    "N_eq2" = rep(0, length(fnqs)*41),
                    "R0_B1" = rep(0, length(fnqs)*41),
                    "R0_B2" = rep(0, length(fnqs)*41),
                    "R0_B3" = rep(0, length(fnqs)*41))

  for(i in 1:nrow(tests)){
    tests[i,3] = get_N.eq(muPq = tests[i,2], f_Nq = tests[i,1])[1]
    tests[i,4] = get_N.eq(muPq = tests[i,2], f_Nq = tests[i,1])[2]
    tests[i,5] = get_N.eq(muPq = tests[i,2], f_Nq = tests[i,1])[3]
    tests[i,6] = get_Ro(muPq = tests[i,2], f_Nq = tests[i,1], beta = 1.1128e-6*(10/3))
    tests[i,7] = get_Ro(muPq = tests[i,2], f_Nq = tests[i,1], beta = 7.038e-6*(10/3))
    tests[i,8] = get_Ro(muPq = tests[i,2], f_Nq = tests[i,1], beta = 4.8372e-5)
  }
  
  plot(x=tests$muPq[tests$f_Nq==1] ,y=tests$N_eq2[tests$f_Nq==1], type='l', lwd=2, col = 'red', ylim = c(0,12000),
       xlab = expression(paste('Predator mortality (', mu[P][q], ')', sep='')),
       ylab = expression(paste('Predicted equilibrium snail population')))
    lines(x=tests$muPq[tests$f_Nq==1.19] , y=tests$N_eq2[tests$f_Nq==1.19], col = 'pink', lwd=2)
    lines(x=tests$muPq[tests$f_Nq==1.74] , y=tests$N_eq2[tests$f_Nq==1.74], col = 'orange', lwd=2)
    lines(x=tests$muPq[tests$f_Nq==1.62] , y=tests$N_eq2[tests$f_Nq==1.62], col = 'blue', lwd=2)
    legend('bottomright', legend = c('Ref', 'Fe', 'At', 'At:Fe'), lwd=2, col=c('red', 'pink', 'orange', 'blue'))
    
    
#Plot heat map of R0 values between pred mortality and snailbirth rate 
  tests.gg1<-subset(tests, R0_B1 <= 1.025 & R0_B1 >= 0.985)
  tests.gg2<-subset(tests, R0_B2 <= 1.015 & R0_B2 >= 0.985)
  tests.gg3<-subset(tests, R0_B3 <= 1.03 & R0_B3 >= 0.97)

 
    
    ggplot(tests, aes(x=muPq, y=f_Nq, fill=R0_B1))+
      theme_bw()+
      theme(axis.title=element_text(size=14),
            axis.text=element_text(size=15))+
      scale_y_continuous(breaks=c(1, 1.2, 1.4, 1.6, 1.8, 2.0))+
      scale_x_continuous(breaks=c(0, 0.025, 0.05, 0.075, 0.1))+
      geom_tile(color='white', size=0.1)+
      scale_fill_continuous(low='green', high='red') +
      geom_line(data = tests.gg1, aes(x=muPq, y=f_Nq), col = 'white', size=1.5)+
      geom_line(data = tests.gg2, aes(x=muPq, y=f_Nq), col = 'grey50', size=1)+
      geom_line(data = tests.gg3, aes(x=muPq, y=f_Nq), col = 'black', size=1)

    
#Check out influence of snail carrying capacity scalars across muPq profile ##################### 
    phinqs<-seq(1, 2, by=0.025) 
    mups<-seq(0, 0.1, by=0.0025)
    
    tests2<-data.frame("phi_Nq" = rep(phinqs, 41),
                      "muPq" = rep(mups, each=41),
                      "P_eq" = rep(0, length(phinqs)*41),
                      "N_eq1" = rep(0, length(phinqs)*41),
                      "N_eq2" = rep(0, length(phinqs)*41),
                      "R0_B1" = rep(0, length(phinqs)*41),
                      "R0_B2" = rep(0, length(phinqs)*41),
                      "R0_B3" = rep(0, length(phinqs)*41))
    
    for(i in 1:nrow(tests2)){
      tests2[i,3] = get_N.eq(muPq = tests2[i,2], phi_Nq = tests2[i,1])[1]
      tests2[i,4] = get_N.eq(muPq = tests2[i,2], phi_Nq = tests2[i,1])[2]
      tests2[i,5] = get_N.eq(muPq = tests2[i,2], phi_Nq = tests2[i,1])[3]
      tests2[i,6] = get_Ro(muPq = tests2[i,2], phi_Nq = tests2[i,1], beta = 1.00e-5)
      tests2[i,7] = get_Ro(muPq = tests2[i,2], phi_Nq = tests2[i,1], beta = 2.4376e-6)
      tests2[i,8] = get_Ro(muPq = tests2[i,2], phi_Nq = tests2[i,1], beta = 2.8372e-5)
    }
    
    
  #Plot heat map of R0 values between pred mortality and snail carrying capacity  
    tests2.gg1<-subset(tests2, R0_B1 <= 1.01 & R0_B1 >= 0.99)
    tests2.gg2<-subset(tests2, R0_B2 <= 1.005 & R0_B2 >= 0.995)
    tests2.gg3<-subset(tests2, R0_B3 <= 1.015 & R0_B3 >= 0.985)
    
    
    
    ggplot(tests2, aes(x=muPq, y=phi_Nq, fill=R0_B3))+
      theme_bw()+
      theme(axis.title=element_text(size=14),
            axis.text=element_text(size=15))+
      scale_y_continuous(breaks=c(1, 1.2, 1.4, 1.6, 1.8, 2.0))+
      scale_x_continuous(breaks=c(0, 0.025, 0.05, 0.075, 0.1))+
      geom_tile(color='white', size=0.1)+
      scale_fill_continuous(low='green', high='red') +
      geom_line(data = tests2.gg1, aes(x=muPq, y=phi_Nq), col = 'black', size=1)+
      geom_line(data = tests2.gg2, aes(x=muPq, y=phi_Nq), col = 'black', size=1)+
      geom_line(data = tests2.gg3, aes(x=muPq, y=phi_Nq), col = 'black', size=1)
  
    plot(x=sort(unique(tests2$muPq), decreasing = TRUE), 
         y = sort(c(unique(tests2$P_eq),rep(0,8))), 
         type = 'l', lwd=2, col='red', 
         xlab = expression(paste(mu[p][q], sep='')),
         ylab  = expression(paste('P'^'*', sep='')))
    axis(side = 3,at = c(0,0.025,0.05,0.075,0.10))
    axis(side = 4,at = c(0,54,63,74.5), col='green')
    abline(a = 53.9, b=0, lty=2, lwd=2)
    abline(a = 63.25, b=0, lty=2, lwd=2, col='lightgreen')
    abline(a = 74.5, b=0, lty=2, lwd=2 , col='green')
    
    
  tests2$phi_Nq<-factor(tests2$phi_Nq, levels = sort(unique(tests2$phi_Nq)))
    
    
    plot(x=tests2$muPq[tests2$phi_Nq==1] ,y=tests2$N_eq2[tests2$phi_Nq==1], type='l', lwd=2, col = 'red', ylim = c(0,20000),
         xlab = expression(paste('Predator mortality (', mu[P][q], ')', sep='')),
         ylab = expression(paste('Predicted equilibrium snail population')))
    lines(x=tests2$muPq[tests2$phi_Nq=='1.175'] , y=tests2$N_eq2[tests2$phi_Nq=='1.175'], col = 'pink', lwd=2)
    lines(x=tests2$muPq[tests2$phi_Nq=='1.6'] , y=tests2$N_eq2[tests2$phi_Nq=='1.6'], col = 'orange', lwd=2)
    lines(x=tests2$muPq[tests2$phi_Nq=='1.5'] , y=tests2$N_eq2[tests2$phi_Nq=='1.5'], col = 'blue', lwd=2)
    legend('bottomright', legend = c('Ref', 'Fe', 'At', 'At:Fe'), lwd=2, col=c('red', 'pink', 'orange', 'blue'))
      
    
  plot(tests$muPq[tests$f_Nq==1], tests$N_eq1[tests$f_Nq==1], pch=16, type='l', lwd=2,
       xlab = expression(paste('Predator mortality (', mu[P][q], ')', sep='')), 
       ylab = "Equilibrium N", ylim=c(-5000, 20000))
    lines(tests$muPq[tests$f_Nq==1], tests$N_eq2[tests$f_Nq==1], pch=16, col='red', lwd=2)
  
  #Something weird going on @muPq = 0.018; investigate
    tests2<-data.frame("phiNq" = rep(1, 1001),
                      "muPq" =seq(0.0175, 0.0185, by=0.000001),
                      "P_eq" = rep(0, 1001),
                      "N_eq1" = rep(0, 1001),
                      "N_eq2" = rep(0, 1001),
                      "R0" = rep(0, 1001))
    
    for(i in 1:nrow(tests2)){
      tests2[i,3] = get_N.eq(muPq = tests2[i,2], phi_Nq = tests2[i,1])[1]
      tests2[i,4] = get_N.eq(muPq = tests2[i,2], phi_Nq = tests2[i,1])[2]
      tests2[i,5] = get_N.eq(muPq = tests2[i,2], phi_Nq = tests2[i,1])[3]
      tests2[i,6] = get_Ro(muPq = tests2[i,2], phi_Nq = tests2[i,1], beta = beta_0)
    }
    
    plot(tests2$muPq[tests2$phiNq==1], tests2$N_eq1[tests2$phiNq==1], pch=16, 
         xlab = expression(paste('Predator mortality (', mu[P][q], ')', sep='')), 
         ylab = "Equilibrium N", ylim=c(-5000, 20000))
    points(tests2$muPq[tests2$phiNq==1], tests2$N_eq2[tests2$phiNq==1], pch=16, col='red', lwd=2)
    
    ggplot(tests, aes(x=P_eq, y=N_eq2, fill=R0))+
      theme_bw()+
      #scale_y_continuous(breaks=c(1, 1.2, 1.4, 1.6, 1.8, 2.0))+
      #scale_x_continuous(breaks=c(0, 0.025, 0.05, 0.075, 0.1))+
      geom_tile(color='white', size=0.1)+
      #coord_equal()+
      scale_fill_continuous(low='green', high='red')

#R0 function where equilibrium snails and predator populations are set    
get_Ro.eq<-function(N_eq, P_eq)
    {   
      f_N<-parameters["f_N"]
      phi_N<-parameters["phi_N"]
      z<-parameters["z"]
      mu_N<-parameters["mu_N"]
      beta<-parameters["beta"]
      sigma<-parameters["sigma"]
      mu_I<-parameters["mu_I"]
      alpha<-parameters["alpha"]
      Th<-parameters["Th"]
      f_P<-parameters["f_P"]
      phi_P<-parameters["phi_P"]
      mu_P<-parameters["mu_P"]
      lamda<-parameters["lamda"]
      mu_W<-parameters["mu_W"]
      m<-parameters["m"]
      H<-parameters["H"]
      mu_H<-parameters["mu_H"]
      
      pred<-(alpha*P_eq)/(1+(alpha*N_eq*Th))#death rate of snails due to predators given equilibrium estimates of P and N
      
      T1<-0.5*beta*m*H*N_eq
      T2<-lamda*sigma
      T3<- (mu_W+mu_H)*(mu_N+pred+sigma)*(mu_N+pred+mu_I)
      
      Ro_est <- sqrt((T1*T2)/T3)

      Ro_est 
      
}

get_Ro.eq(N_eq = 10000, P_eq = 5)
    
#Get 10-day chlorP data to incorporate #################
ecotox_10.<-data.frame('chem'=rep("Chlorpyrifos", 30),
                      'dose'=c(rep(0,5), rep(0.64,5), rep(3.2,5), 
                               rep(6.4,5), rep(32,5), rep(64,5)),
                      'dead'=c(rep(0,5), rep(0,5), rep(0,5), 
                               rep(0,2), rep(1,13)))

  ecotox10..<-glm(dead ~ dose, family=binomial(link="probit"),data=ecotox_10.)
  summary(ecotox10..)
  
  #Extrapolate response to constant gradient of Chlorpyrifos concentration
  p.ecotox10<-data.frame(dose=seq(from=0, to=100, by=0.01))
  p.ecotox10[, c('mortality', 'st.er')]<-predict(ecotox10.., p.ecotox10, 
                                                 type = "response", se.fit=TRUE)
  
  for(i in 1:nrow(p.ecotox10)){
    p.ecotox10[i,4] = get_N.eq(muPq = p.ecotox10[i,2]/10, phi_Nq = 1)[1]
    p.ecotox10[i,5] = get_N.eq(muPq = p.ecotox10[i,2]/10, phi_Nq = 1)[2]
    p.ecotox10[i,6] = get_N.eq(muPq = p.ecotox10[i,2]/10, phi_Nq = 1)[3]
    p.ecotox10[i,7] = get_Ro(muPq = p.ecotox10[i,2]/10, phi_Nq = 1, beta = beta_0)
    p.ecotox10[i,8] = get_N.eq(muPq = p.ecotox10[i,2]/10, phi_Nq = 1.16)[1]
    p.ecotox10[i,9] = get_N.eq(muPq = p.ecotox10[i,2]/10, phi_Nq = 1.16)[2]
    p.ecotox10[i,10] = get_N.eq(muPq = p.ecotox10[i,2]/10, phi_Nq = 1.16)[3]
    p.ecotox10[i,11] = get_Ro(muPq = p.ecotox10[i,2]/10, phi_Nq = 1.16, beta = beta_0)
    p.ecotox10[i,12] = get_N.eq(muPq = p.ecotox10[i,2]/10, phi_Nq = 1.61)[1]
    p.ecotox10[i,13] = get_N.eq(muPq = p.ecotox10[i,2]/10, phi_Nq = 1.61)[2]
    p.ecotox10[i,14] = get_N.eq(muPq = p.ecotox10[i,2]/10, phi_Nq = 1.61)[3]
    p.ecotox10[i,15] = get_Ro(muPq = p.ecotox10[i,2]/10, phi_Nq = 1.61, beta = beta_0)
    p.ecotox10[i,16] = get_N.eq(muPq = p.ecotox10[i,2]/10, phi_Nq = 1.51)[1]
    p.ecotox10[i,17] = get_N.eq(muPq = p.ecotox10[i,2]/10, phi_Nq = 1.51)[2]
    p.ecotox10[i,18] = get_N.eq(muPq = p.ecotox10[i,2]/10, phi_Nq = 1.51)[3]
    p.ecotox10[i,19] = get_Ro(muPq = p.ecotox10[i,2]/10, phi_Nq = 1.51, beta = beta_0)
  }
  colnames(p.ecotox10)[c(4:19)]<-c('Ch_P_eq', 'Ch_N_eq1', 'Ch_N_eq2', 'Ch_R0',
                                   'Ch_Fe_P_eq', 'Ch_Fe_N_eq1', 'Ch_Fe_N_eq2', 'Ch_Fe_R0',
                                   'Ch_At_P_eq', 'Ch_At_N_eq1', 'Ch_At_N_eq2', 'Ch_At_R0',
                                   'Ch_At_Fe_P_eq', 'Ch_At_Fe_N_eq1', 'Ch_At_Fe_N_eq2', 'Ch_At_Fe_R0')
  
  plot(x=p.ecotox10$dose ,y=p.ecotox10$Ch_N_eq2, type='l', 
       lwd=2, col = 'red', ylim = c(0,20000), xlim=c(0,10),
       xlab = expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep='')),
       ylab = expression(paste('Predicted equilibrium snail population')))
  lines(x=p.ecotox10$dose , y=p.ecotox10$Ch_Fe_N_eq2, col = 'pink', lwd=2)
  lines(x=p.ecotox10$dose , y=p.ecotox10$Ch_At_N_eq2, col = 'orange', lwd=2)
  lines(x=p.ecotox10$dose , y=p.ecotox10$Ch_At_Fe_N_eq2, col = 'blue', lwd=2)
  legend('bottomright', legend = c('Ref', 'Fe', 'At', 'At:Fe'), lwd=2, 
         col=c('red', 'pink', 'orange', 'blue'))
  
#Get 4-day chlorP data to incorporate #############
ecotox_4.<-data.frame('chem'=rep("Chlorpyrifos", 30),
                         'dose'=c(rep(0,5), rep(0.64,5), rep(3.2,5), 
                                  rep(6.4,5), rep(32,5), rep(64,5)),
                         'dead'=c(rep(0,5), rep(0,5), rep(0,5), 
                                  rep(0,5), 0, rep(1,9)))
  
ecotox_4..<-glm(dead ~ dose, family=binomial(link="probit"),data=ecotox_4.)
  summary(ecotox_4..)
  
#Extrapolate response to constant gradient of Chlorpyrifos concentration
  p.ecotox4<-data.frame(dose=seq(from=0, to=100, by=0.01))
  p.ecotox4[, c('mortality', 'st.er')]<-predict(ecotox_4.., p.ecotox4, 
                                                 type = "response", se.fit=TRUE)
  
  for(i in 1:nrow(p.ecotox4)){
    p.ecotox4[i,4] = get_N.eq(muPq = p.ecotox4[i,2]/4, phi_Nq = 1)[1]
    p.ecotox4[i,5] = get_N.eq(muPq = p.ecotox4[i,2]/4, phi_Nq = 1)[2]
    p.ecotox4[i,6] = get_N.eq(muPq = p.ecotox4[i,2]/4, phi_Nq = 1)[3]
    p.ecotox4[i,7] = get_Ro(muPq = p.ecotox4[i,2]/4, phi_Nq = 1, beta = beta_0)
    p.ecotox4[i,8] = get_N.eq(muPq = p.ecotox4[i,2]/4, phi_Nq = 1.16)[1]
    p.ecotox4[i,9] = get_N.eq(muPq = p.ecotox4[i,2]/4, phi_Nq = 1.16)[2]
    p.ecotox4[i,10] = get_N.eq(muPq = p.ecotox4[i,2]/4, phi_Nq = 1.16)[3]
    p.ecotox4[i,11] = get_Ro(muPq = p.ecotox4[i,2]/4, phi_Nq = 1.16, beta = beta_0)
    p.ecotox4[i,12] = get_N.eq(muPq = p.ecotox4[i,2]/4, phi_Nq = 1.61)[1]
    p.ecotox4[i,13] = get_N.eq(muPq = p.ecotox4[i,2]/4, phi_Nq = 1.61)[2]
    p.ecotox4[i,14] = get_N.eq(muPq = p.ecotox4[i,2]/4, phi_Nq = 1.61)[3]
    p.ecotox4[i,15] = get_Ro(muPq = p.ecotox4[i,2]/4, phi_Nq = 1.61, beta = beta_0)
    p.ecotox4[i,16] = get_N.eq(muPq = p.ecotox4[i,2]/4, phi_Nq = 1.51)[1]
    p.ecotox4[i,17] = get_N.eq(muPq = p.ecotox4[i,2]/4, phi_Nq = 1.51)[2]
    p.ecotox4[i,18] = get_N.eq(muPq = p.ecotox4[i,2]/4, phi_Nq = 1.51)[3]
    p.ecotox4[i,19] = get_Ro(muPq = p.ecotox4[i,2]/4, phi_Nq = 1.51, beta = beta_0)
  }
  colnames(p.ecotox4)[c(4:19)]<-c('Ch_P_eq', 'Ch_N_eq1', 'Ch_N_eq2', 'Ch_R0',
                                   'Ch_Fe_P_eq', 'Ch_Fe_N_eq1', 'Ch_Fe_N_eq2', 'Ch_Fe_R0',
                                   'Ch_At_P_eq', 'Ch_At_N_eq1', 'Ch_At_N_eq2', 'Ch_At_R0',
                                   'Ch_At_Fe_P_eq', 'Ch_At_Fe_N_eq1', 'Ch_At_Fe_N_eq2', 'Ch_At_Fe_R0')
  
  plot(x=p.ecotox4$dose ,y=p.ecotox4$Ch_N_eq2, type='l', 
       lwd=2, col = 'red', ylim = c(0,20000), xlim=c(0,30),
       xlab = expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep='')),
       ylab = expression(paste('Predicted equilibrium snail population')))
  lines(x=p.ecotox4$dose , y=p.ecotox4$Ch_Fe_N_eq2, col = 'pink', lwd=2)
  lines(x=p.ecotox4$dose , y=p.ecotox4$Ch_At_N_eq2, col = 'orange', lwd=2)
  lines(x=p.ecotox4$dose , y=p.ecotox4$Ch_At_Fe_N_eq2, col = 'blue', lwd=2)
  legend('topleft', legend = c('Ref', 'Fe', 'At', 'At:Fe'), lwd=2, 
         col=c('red', 'pink', 'orange', 'blue'))
  
  
  
#3d response of N_eq to variable atrazine and chlorP doses ###################
  