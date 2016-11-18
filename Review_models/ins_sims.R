#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(deSolve)
source('Review_models/agroReview_mod1.1.R')
source('Review_models/Response_fxs.R')

mod.ins.test = function(t, n, parameters) {
  with(as.list(parameters),{
    
    S=n[1]  #Susceptible snails
    E=n[2]  #Exposed snails
    I=n[3]  #Infected snails
    W=n[4]  #Mean worm burden
    P=n[5]  #Predator population
    He=n[6] #Herbicide concentration
    Fe=n[7] #Fertilizer concentration
    In=n[8] #Insecticide concentration
    
    
    #Dynamic variables ####################
    #Total number of snails
    N = S+E+I  
    #mating function
    fx<-function(x, mean.worm = W, clump = k){
      (1-cos(x))/((1+((mean.worm/(clump+mean.worm))*cos(x)))^1+clump)
    }
    gamma = 1 - ((1-(W/(W+k)))^(1+k)/2*pi)*(integrate(fx, 0, 2*pi)$value)       
    
    #Agrochem parameters (set to null values for now)
    phi_Nq = phi_N #* (He_fx(He) + Fe_fx(Fe))
    mu_Nq = mu_N # 
    mu_Pq = In_fx(In)
      if(mu_Pq != parameters['mu_P']){
         print(In_fx(In))
      }
    v_q = v #* f4(Q1)
    theta_q = theta #* f5(Q1)
    pi_Mq = pi_M #* f6(Q1)
    pi_Cq = pi_C #* f7(Q1)
    pred_red = 1 #In_fx(In)[2]
    
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
    
    dWdt= lamda*Om*C - (mu_W+mu_H)*W                                      #mean worm burden in human population
    
    dPdt= f_P*(1-P/(phi_P*A))*P - mu_Pq*P                                    #prawn population (number individuals)
    
    dHedt= -k_He*He                                                    #Herbicide concentration
    
    dFedt= -k_Fe*Fe                                                    #Fertilizer concentration
    
    dIndt= -k_In*In                                                   #Insecticide concentration
    
    return(list(c(dSdt, dEdt, dIdt, dWdt, dPdt, dHedt, dFedt, dIndt)))
  }) 
}

In_fx.null = function(In){
  return(parameters['mu_P'])
}

#calculate decay rates (k) from hydrolysis half lives for each chemical in Halstead 2015 (from table S1 and pmep.cce.cornell.edu) ############
mal.k = -log(0.5)/6.2         #from table s1 and in agreement of "less than 1 week in raw river water" from Cornell
chlor.k = -log(0.5)/25.5      #from table S1; within the range of reported half life from cornell
terb.k = -log(0.5)/6.5        #from table s1; in agreement with Cornell estimate of 5.5 days at pH of 7; degrades into formaldehyde
lamcy.k = -log(0.5)/1         #very fast; 0 in table S1; "Not expected to be prevalent in surface waters" according to cornell website
esfen.k = -log(0.5)/10        #cornell 4-15 days half life in water
perm.k = -log(0.5)/2          #cornell "half life of less than 2.5 days"

#median EECs from Halstead 2015
med.mal = 0.778
med.chlor = 5.810
med.terb = 1.435
med.lamcy = 0.649
med.esfen = 0.311
med.perm = 1.420

#max EECs from Halstead 2015
max.mal = 18.4
max.chlor = 28.0
max.terb = 36.6
max.lamcy = 1.77
max.esfen = 1.03
max.perm = 5.98

#Model structure to dynamically vary insecticide application ###################

ins.test<-function(chem, chem.k, chem.med){
  
  parameters['k_In'] = chem.k
  time = seq(1, 365, by=1)

  nstart = c(S = as.numeric(p.eqbm[2]),
             E = as.numeric(p.eqbm[3]),
             I = as.numeric(p.eqbm[4]),
             W = as.numeric(p.eqbm[5]),
             P = as.numeric(p.eqbm[6]),
             He = as.numeric(p.eqbm[7]),
             Fe = as.numeric(p.eqbm[8]),
             In = as.numeric(p.eqbm[9]))
  
  events.df<-data.frame(var= 'In',
                        time = 50,
                        value = chem.med,
                        method = 'add')
  
  output.all2 = as.data.frame(ode(nstart, time, mod.ins.test, parameters,
                                 events = list(data = events.df)))
  output.all2$N = output.all2$S + output.all2$E + output.all2$I

par(mfrow = c(4,1), mar = c(2.5,4,1.5,1.0), cex.lab = 1.2)
  plot(output.all$time, output.all$In,lwd=2, col = 'red', xlab = 'time', ylab = 'concentration (ppb)', 
       type = 'l', main = chem, ylim = c(0, 30))
  
  plot(output.all$time, output.all$P, lwd=2, type = 'l', col = 'darkgreen', 
       xlab = 'time', ylab = 'predators (P)', ylim = c(0, parameters['phi_P']))
  
  plot(x = output.all$time, y = output.all$N, lwd=2, xlab = 'time', ylab = 'snail dynamics', 
       type = 'l', ylim = c(0, max(output.all$N)+10))
      lines(output.all$time, output.all$S, lwd=2, col = 'blue')
      lines(output.all$time, output.all$E, lwd=2, col = 'orange')
      lines(output.all$time, output.all$I, lwd=2, col = 'red')
      legend('topleft', legend = c('N', 'S', 'E', 'I'), lty = 1, 
             col = c('black', 'blue', 'orange', 'red'), cex = 0.8,
             ncol = 2, bty = 'n')
      
  plot(output.all$time, output.all$W, lwd=2, col = 'purple', type = 'l', 
       xlab = 'time', ylab = 'mean worm burden', ylim = c(20,23))  
  
}

#Simulations with different agrochemicals, prawn addition timing ###################
#Run with no agrochemical introducation
In_fx<-function(In){
  return(In_fx.null(In))
}
ins.test(chem = 'None',
         chem.k = chlor.k,
         chem.med = 0)

#Run with chlorpyrifos
In_fx<-function(In){
  return(f_muPq_ch_Halstead(In))
}
ins.test(chem = 'Chlorpyrifos',
         chem.k = chlor.k,
         chem.med = max.chlor)

#Run with esfenvalerate
In_fx<-function(In){
  return(f_muPq_esfen_Halstead(In))
}
ins.test(chem = 'Esfenvalerate',
         chem.k = esfen.k,
         chem.med = max.esfen)

#Run with permethrin
In_fx<-function(In){
  return(f_muPq_perm_Halstead(In))
}
ins.test(chem = 'Permethrin',
         chem.k = perm.k,
         chem.med = max.perm)

#Run with lambda-cyhalothrn
In_fx<-function(In){
  return(f_muPq_lamcy_Halstead(In))
}
ins.test(chem = 'L-Cyhalothrin',
         chem.k = lamcy.k,
         chem.med = max.lamcy)

#Run with terbufos
In_fx<-function(In){
  return(f_muPq_terb_Halstead(In))
}
ins.test(chem = 'Terbufos',
         chem.k = terb.k,
         chem.med = max.terb)

#Run with malathion
In_fx<-function(In){
  return(f_muPq_mal_Halstead(In))
}
ins.test(chem = 'Malathion',
         chem.k = mal.k,
         chem.med = max.mal)

par(mfrow = c(1,1), mar = c(4.5,4,1.5,1.0))
plot(output.all$time, output.all$W, lwd=2, lty = 2, col = 'purple', type = 'l', 
     xlab = 'time (days)', ylab = 'mean worm burden', ylim = c(20,23), xlim = c(0, 365))  

  lines(output.all2$time, output.all2$W, lwd=2, col = 'purple', type = 'l')
  legend('topleft', legend = c('Esfenvalerate', 'Chlorpyrifos'), lty = c(1,2), lwd = 2, 
         col = c('purple'), ncol = 2, bty = 'n')