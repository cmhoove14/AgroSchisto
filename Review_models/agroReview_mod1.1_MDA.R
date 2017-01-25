source('Response_Fxs/griggs08_piC.R')
source('Response_Fxs/Baxter_Rohr_Atrazine2011.R')
source('Review_models/agroReview_mod1.1.R')
require(deSolve)

mod.mda = function(t, n, parameters) {
  with(as.list(parameters),{
    
    S=n[1]  #Susceptible snails
    E=n[2]  #Exposed snails
    I=n[3]  #Infected snails
    Wu=n[4] #Worm burden untreated
    Wt=n[5] #Worm burden treated
    P=n[6]  #Predator population
    He=n[7] #Herbicide concentration
    Fe=n[8] #Fertilizer concentration
    In=n[9] #Insecticide concentration
    
    W=cov*Wt + (1-cov)*Wu
    
    #Dynamic variables ####################
    #Total number of snails
    N = S+E+I  
    #mating function
    fx<-function(x, mean.worm = W, clump = k){
      (1-cos(x))/((1+((mean.worm/(clump+mean.worm))*cos(x)))^1+clump)
    }
    gamma = 1 - ((1-(W/(W+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value)       
    
    #Agrochem parameters (set to null values for now)
    phi_Nq = f_phi_Nq_at(He) #* (He_fx(He) + Fe_fx(Fe))
    mu_Nq = mu_N # 
    mu_Pq = mu_P
    v_q = v #* f4(Q1)
    theta_q = theta #* f5(Q1)
    pi_Mq = pi_M #* f6(Q1)
    pi_Cq = pi_C_atr_grg08(He) #* f7(Q1)
    pred_red = 1
    
    #per capita predation of predators on snails with variable functional response based on parameter nn  
    pred= ((alpha*P)/(1+(alpha*((N/A)^nn)*Th)))*pred_red   
    
    #Schistosome larval concentration equations
    Wf = 0.5*W*H*gamma
    
    M = Wf*m*v_q*pi_Mq 
    
    C = theta_q*I*pi_Cq 
    
    #print(c(paste('M = ', round(M, digits = 2), sep = ''),
    #        paste('C = ', round(C, digits = 2), sep = ''),
    #        paste('gam = ', round(gamma, digits = 2), sep = ''),
    #        paste('pred = ', pred, sep = ''),
    #paste('N = ', round(N, digits = 2), sep = '')))
    
    #differential equations ###################
    
    dSdt= f_N*(1-N/(phi_Nq*A))*(S+z*E) - mu_Nq*S - pred*(S/A)^nn - beta*M*Om*S       #Susceptible snails
    
    dEdt= beta*M*S - mu_Nq*E - pred*(E/A)^nn - sigma*E                           #Exposed snails
    
    dIdt= sigma*E - mu_Nq*I - mu_I*I - pred*(I/A)^nn                             #Infected snails
    
    dWudt= lamda*Om*C - (mu_W+mu_H)*Wu
    
    dWtdt= lamda*Om*C - (mu_W+mu_H)*Wt                                      #mean worm burden in human population
    
    dPdt= f_P*(1-P/(phi_P*A))*P - mu_Pq*P                                    #prawn population (number individuals)
    
    dHedt= -k_He*He                                                    #Herbicide concentration
    
    dFedt= -k_Fe*Fe                                                    #Fertilizer concentration
    
    dIndt= -k_In*In                                                    #Insecticide concentration
    
    return(list(c(dSdt, dEdt, dIdt, dWudt, dWtdt, dPdt, dHedt, dFedt, dIndt)))
  }) 
}

nstart.mda = c(S=output2[max(output2$time),c(2)], 
               E=output2[max(output2$time),c(3)], 
               I=output2[max(output2$time),c(4)], 
               Wu=output2[max(output2$time),c(5)], 
               Wt=output2[max(output2$time),c(5)], 
               P=output2[max(output2$time),c(6)], 
               He=0,
               Fe=0,
               In=0)

cov = 0.80
n.years = 10
  t.total = n.years*365
  t.mda = seq(0, t.total,1)
mda.year = 100
mda.time = round(seq(100, t.total-100, length.out = n.years))
mda.eff = 0.95

mda.events = data.frame(var=rep('Wt', times = n.years),
                        time = mda.time,
                        value = rep((1-mda.eff), times = n.years),
                        method = rep("mult", times = n.years))

output.mda = as.data.frame(ode(nstart.mda, t.mda, mod.mda, parameters,
                               events = list(data = mda.events)))
output.mda$N = output.mda$S + output.mda$E + output.mda$I
output.mda$W = output.mda$Wt*cov + output.mda$Wu*(1-cov)

plot(output.mda$time, output.mda$W, type = 'l', lwd = 2, col = 'purple', ylim = c(0,60),
     xlab = 'time (Days)', ylab = 'Mean Worm Burden')

nstart.mda.he = c(S=output2[max(output2$time),c(2)], 
               E=output2[max(output2$time),c(3)], 
               I=output2[max(output2$time),c(4)], 
               Wu=output2[max(output2$time),c(5)], 
               Wt=output2[max(output2$time),c(5)], 
               P=output2[max(output2$time),c(6)], 
               He=50,
               Fe=0,
               In=0)


output.mda.he = as.data.frame(ode(nstart.mda.he, t.mda, mod.mda, parameters,
                               events = list(data = mda.events)))
output.mda.he$N = output.mda.he$S + output.mda.he$E + output.mda.he$I
output.mda.he$W = output.mda.he$Wt*cov + output.mda.he$Wu*(1-cov)

  lines(output.mda.he$time, output.mda.he$W, type = 'l', lwd = 2, lty = 2, col = 'purple')

plot(output.mda.he$time, output.mda.he$He, type = 'l', lwd = 2, col = 'gold2', ylim = c(0,50),
     xlab = 'time (Days)', ylab = 'Atrazine concentration')
  
plot(output.mda$time, output.mda$N, type = 'l', lwd = 2, col = 'blue', ylim = c(0,50),
     xlab = 'time (Days)', ylab = 'Snail Density')
  lines(output.mda.he$time, output.mda.he$N, type = 'l', lwd = 2, lty = 2, col = 'blue')