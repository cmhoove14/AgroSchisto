parameters=c(#Excluding beta, Phi_Nq, f_Nq, and muPq which will be read into the R0 function
  ##standard snail parameters 
  f_N=0.10, # recruitment rate (from sokolow et al)
  phi_N=10000, # carrying capacity from sokolow et al
  z=0.5, #Proportion of exposed snails that reproduce from sokolow et al
  mu_N=1/60, #Mortality rate from Sokolow et al
  sigma=1/40, #Transition rate from exposed to infected from sokolow et al
  mu_I=1/10 - 1/60, #additional snail death due to infection from sokolow et al
  
  #prawn parameters
  alpha=0.003, #attack rate
  Th=0.067,#~Prawn predation limit
  f_P=0.234/2, #prawn birth rate from Cervantes-Santiago Aquaculture 2010 paper (/2 for 1:1 female-male ratio)
  phi_P=120,  #prawn carrying capacity
  mu_P= muP, #observed daily 24 hr prawn mortality rate in chlorP-free tanks of mesocosm#0.00137, previous rate
  
  #Adult Worm, Miracidia and Circariae Parameters
  mu_W=1/(3.3*365), # death rate of adult worms
  m=0.5, #miracidial shedding rate per reproductive female divided by miracidial mortality; from sokolow et al
  
  #Human parameters
  H=300, #number of humans
  mu_H=1/(60*365) #Assumes 60 year lifespan
)


get_Ro_beta_lamda<-function(muPq = 0, phi_Nq = 1, beta, lamda, f_Nq = 1) #variable parameters to be manipulated
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
  #Shorthand values to use in N_eq expression
  N_eq = max(uniroot.all(f = function(y){(f_N*f_Nq)*(1-y/(phi_N*phi_Nq)) - 
      mu_N - 
      (P_eq*alpha)/(1+alpha*Th*y)}, 
      c(0, as.numeric(phi_N*phi_Nq))))
  
  if(N_eq < 0){
    N_eq = 0
  }
  
  pred<-(alpha*P_eq)/(1+(alpha*N_eq*Th))#death rate of snails due to predators given equilibrium estimates of P and N
  
  Ro_est <- sqrt((0.5*beta*H*N_eq*lamda*sigma)/((mu_W+mu_H)*(mu_N+pred+sigma)*(mu_N+pred+mu_I)))
  
  return(c(N_eq,P_eq,Ro_est, pred))
  
} #End R0 function  


r0.atra.chlor<-data.frame("Atra" = rep(c(0:100), times = 65),
                          "logatra" = rep(log(c(0:100)+1), times = 65),
                          "dose" = rep(c(0:64), each = 101),
                          "phi_Nq" = 0,
                          "rate" = 0,
                          "R0" =0)


r0.atra.chlor$rate<-(predict(ecotox10_mod, r0.atra.chlor, 
                             type = "response", se.fit=TRUE)$fit)/10 #fill mortality rate data from model

r0.atra.chlor$rate = r0.atra.chlor$rate - r0.atra.chlor[1,5] 

r0.atra.chlor$phi_Nq<-(predict(atra_mod, r0.atra.chlor, 
                               type = "response", se.fit=TRUE)$fit) #fill bottom-up effect data from model

r0.atra.chlor$phi_Nq<-r0.atra.chlor$phi_Nq + (1 - r0.atra.chlor[1,4]) #Normalize to 1

for(i in 1:nrow(r0.atra.chlor)){
  r0.atra.chlor[i,6] = get_Ro_beta_lamda(muPq = r0.atra.chlor[i,5],
                                         beta = beta.use,
                                         lamda = lamda.use,
                                         phi_N = r0.atra.chlor[i,4])[3]
  r0.atra.chlor[i,7] = get_Ro_beta_lamda(muPq = r0.atra.chlor[i,5],
                                         beta = beta.use,
                                         lamda = lamda.use,
                                         phi_N = r0.atra.chlor[i,4])[2]
  r0.atra.chlor[i,8] = get_Ro_beta_lamda(muPq = r0.atra.chlor[i,5],
                                         beta = beta.use,
                                         lamda = lamda.use,
                                         phi_N = r0.atra.chlor[i,4])[1]
  r0.atra.chlor[i,9] = get_Ro_beta_lamda(muPq = r0.atra.chlor[i,5],
                                         beta = beta.use,
                                         lamda = lamda.use,
                                         phi_N = r0.atra.chlor[i,4])[4]
}

colnames(r0.atra.chlor)[c(7:9)]<-c('P_eq', 'N_eq', 'pred')


plot(x = r0.atra.chlor$Atra[r0.atra.chlor$dose == 12], y = r0.atra.chlor$R0[r0.atra.chlor$dose == 12],
     xlab = 'atrazine concentration', ylab = 'R0', type = 'l', lwd = 2, ylim = c(0,5))   
for(i in 13:17){
  lines(x = r0.atra.chlor$Atra[r0.atra.chlor$dose == i], y = r0.atra.chlor$R0[r0.atra.chlor$dose == i],
        lwd = 2, col = 2+(i-13))
}

plot(x = r0.atra.chlor$pred[r0.atra.chlor$dose == 12], y = r0.atra.chlor$N_eq[r0.atra.chlor$dose == 12], 
     type='l', lwd=2, xlab ='', ylab ='')



#Step 3.3 combine chlorpyrifos and atrazine d/r data in heat map ##############
r0.atra.chlor<-data.frame("Atra" = 0,
                          "logatra" = rep(seq(from=0, to=log(101), by=0.1), times = 65),
                          "dose" = rep(seq(0,64,1), each = 47),
                          "logdose" = rep(log(seq(0,64,1)+1), each = 47),
                          "phi_Nq" = 0,
                          "rate" = 0,
                          "R0" =0)

r0.atra.chlor$Atra<-exp(r0.atra.chlor$logatra)

r0.atra.chlor$rate<-(predict(ecotox10_mod, r0.atra.chlor, 
                             type = "response", se.fit=TRUE)$fit)/10 #fill mortality rate data from model

r0.atra.chlor$rate = r0.atra.chlor$rate - r0.atra.chlor[1,6] #normalize to 0; no excess death with no concentration

r0.atra.chlor$phi_Nq<-(predict(atra_mod, r0.atra.chlor, 
                               type = "response", se.fit=TRUE)$fit) #fill bottom-up effect data from model

r0.atra.chlor$phi_Nq<-r0.atra.chlor$phi_Nq + (1 - r0.atra.chlor[1,5]) #Normalize to 1

for(i in 1:nrow(r0.atra.chlor)){
  r0.atra.chlor[i,7] = get_Ro_beta_lamda(muPq = r0.atra.chlor[i,6],
                                         beta = beta.use,
                                         lamda = lamda.use,
                                         phi_N = r0.atra.chlor[i,5])[3]
  r0.atra.chlor[i,8] = get_Ro_beta_lamda(muPq = r0.atra.chlor[i,6],
                                         beta = beta.use,
                                         lamda = lamda.use,
                                         phi_N = r0.atra.chlor[i,5])[2]
  r0.atra.chlor[i,9] = get_Ro_beta_lamda(muPq = r0.atra.chlor[i,6],
                                         beta = beta.use,
                                         lamda = lamda.use,
                                         phi_N = r0.atra.chlor[i,5])[1]
}

colnames(r0.atra.chlor)[c(8,9)]<-c('P_eq', 'N_eq')

plot(x = r0.atra.chlor$logatra[r0.atra.chlor$dose == 64],
     y = r0.atra.chlor$R0[r0.atra.chlor$dose == 64],
     xlab = 'atrazine concentration', ylab = 'R0',
     type = 'l', lwd = 2)

ggplot(r0.atra.chlor, aes(x=logatra, y=dose, fill=R0))+
  theme_bw()+
  #scale_fill_brewer(type = 'div', palette = 'RdYlGn', direction = -1)+
  scale_fill_distiller(palette = "Spectral")+
  scale_x_continuous(breaks = log(c(0,1,10,100)+1), 
                     limits = log(c(0,100)+1),
                     labels = c('0','1','10','100'))+
  scale_y_continuous(breaks = c(0,20,40,60), limits = c(0,64))+
  geom_raster(interpolate = TRUE)+
  coord_equal(ratio = 1/20)+
  #geom_vline(xintercept = c(unique(r0.atra.chlor$logatra)), col = 'black', alpha = 0.5)+
  geom_hline(yintercept = c(0,20,40,60), col = 'lightgrey', alpha = 0.25)+
  labs(y=expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep = '')), 
       x=expression(paste('log+1 Atrazine concentration (', mu, 'g/L)', sep = '')))+
  theme(axis.ticks=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=14),
        legend.title=element_text(size=15), legend.text=element_text(size=12))+
  geom_text(label='D', x=-3, y=65, size=6, alpha=.50)    
