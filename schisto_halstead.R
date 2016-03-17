#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#function used to model halstead et al mesocosm results impact on human infection

require(deSolve)
require(ggplot2)

#Model structure and equations ####################
  #Represents a simplified version only including agroC effects found significant in mesocosms
schisto_halstead=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    P=n[5]
    N=S+E+I
    
    #Miracidia population
    M=Eg*V #Eg = egg input in mesocosm = 2276 eggs/tank/week = 325/day ; V=8.4%
    
    #Circariae population
    #C=I*Theta*Rc
    
    #Rate of predation 
    #TOOK OUT alpha_q to simplify because of high dose of ChlorP 
    pred= (alpha*P)/(1+(alpha*N*Th)) #death rate of snails due to predators (Prawns)
    
    
  #snail compartment models
  #CHANGED to multiply by phi_Nq so that scalar of snail carrying capacity derived from 
    #mesocosm results can be input directly
  #TOOK OUT f_Nq and mu_Nq because they were not detectable in mesocosm experiments
    dSdt= (f_N)*((1-((phi_N*(1/phi_Nq))*N))*(S+(z*E))) - 
      ((mu_N)*S) - (pred*S) - (beta*(M)*S)
    dEdt= beta*(M)*S - ((mu_N+pred+sigma)*E)
    dIdt= (sigma*E) - ((mu_N+pred+mu_I)*I)
    
    #worm burden in human
    dWdt= (lamda*I) - ((mu_W+mu_H)*W)
    
    #prawn population
    dPdt= (f_P*(1-(phi_P*P))*P)-((mu_P+mu_Pq)*P)
    
    return(list(c(dSdt,dEdt,dIdt,dWdt, dPdt)))
  }) 
} 

#List parameters and values #####################
parameters=c(
  ##standard snail parameters 
    f_N=0.16, # birth rate
    phi_N=(1-1/(80*0.16))/10000, # carrying capacity
    z=0.5, #Proportion of exposed snails that reproduce
    mu_N=1/80, #Mortality rate
    
    beta=0.0000011128, #Best fit data to first two points in Lampsar 2
    
    sigma=0.02, #Transition rate from exposed to infected
    mu_I=0.05, #additional snail death due to infection
  
  ## snail parameters impacted by agrochemicals
    #f_Nq=1, #Not affected in mesocosm
    phi_Nq=1, #Scalar of snail carrying capacity by chemical concentration INFORMED BY BOTTOM UP EFFECTS IN MESOCOSM
    #mu_Nq=0, #Chem concentrations too low to affect in mesocosom experiments
  
  #prawn parameters
    alpha=0.003, #attack rate
    Th=0.1,#~Prawn predation limit
    f_P=0.128,#prawn birth rate
    phi_P=1/50,  #prawn carrying capacity
    mu_P= 0.00137, #prawn mortality
  
  #prawn parameter impacted by agrochemicals
    mu_Pq=0, #can be a function or fixed (0.02 for eg)
    #alpha_q=1, #Scalar of predation rate due to sub-lethal toxicity; not considered in mesocosm
  
  #Adult Worm, Miracidia and Circariae Parameters
  lamda=0.00004, #probability of cercaria surviving to reproduction after infecting a human host
  mu_W=1/(3*365), # death rate of adult worms
  #m=0.5, #Miracidial production by adult female worms NOT NEEDED BECAUSE OF CONSTANT EGG INPUT
  Eg=325*40, #constant daily input (from weekly input reported in halstead et al) scaled up to 200 m2 surface area
  V=.084, #Egg viability controlling schistosome egg->infective miracidia
  #Rm=1, #Scalar to miracidial input due to direct toxicity to larval schistosomes; NOT DETECTED IN MESOCOSM
  #Rc=1, #probability of cercaria infecting human
  #Theta=, #per capita cercarial shedding rate of infected snails
  
  #NOT CONSIDERED IN MESOCOSM: parameters impacted by agrochemicals 
  #V_q=1,#Egg viability (produced eggs that hatch to miracidia)
  #R_q=1,#Larval mortality
  #Theta_q=1,#Scalar to infected snail cercarial shedding rate based on agrochemical exposure
  
  #Human parameters
  H=1000, #number of humans
  mu_H=1/(60*365) #Assumes 60 year lifespan
)

#Run model ###############################

#Set initial values
yrs=30
time=seq(0,365*yrs,1)
p=0

nstart=c(S=7000,E=3750,I=1200, W=20, P=p)
#Run model and save equilibrium values for 5 state variables
output=as.data.frame(ode(nstart,time,schisto_halstead,parameters)) 

eqbm<-output[365*yrs,]
eqbm

#plot time series of estimated values to make sure state variables reached equilibrium ######
  plot(output[,1], output[,2], type='l', xlab="time",ylab="System Variables", 
       ylim=c(0,max( output[,2]+output[,3]+output[,4] )), col=1, lty=1, lwd=2)
    lines(output[,1],output[,3],col=2, lty=2, lwd=2)
    lines(output[,1],output[,4],col=3, lty=2, lwd=2)
    lines(output[,1],output[,2]+output[,3]+output[,4],col=4, lty=2, lwd=2)
    lines(output[,1],output[,5],col=5, lty=2)
    lines(output[,1],output[,6],col=6, lty=2, lwd=2)
    
  #Plot worm burden only  
    plot(output[,1],output[,5],col=5, ylab='Mean worm burden (W)', xlab='time')  

#R_0 function with variable prawn mortality (mu_Pq) and snail carrying capacity (phi_Nq) #######    
  
  #muPq data from halstead et al 2015
  
  chlor<-data.frame(dose=c(0,0.64,3.2,6.4,32,64), 
                    mortality=c(0,0,0,0,0.8, 1.00))
  
  chlor$mu_P<- -0.25*log(1-chlor[,2])
  chlor[6,3]=9.0109
  
  #Probit analysis of data (probit analysis most often used for toxicological 
  #data investigating mortality responses to toxic exposures)
  tox_prawn<-glm(mortality ~ dose, family=binomial(link="probit"),data=chlor)
  summary(tox_prawn)
  
  #Extrapolate response to constant gradient of Chlorpyrifos concentration
  response_prawn<-data.frame(dose=seq(from=0, to=500, by=1))
  response_prawn[, c("mortality", "se")]<-predict(tox_prawn, response_prawn, 
                                                  type = "response", se.fit=TRUE)
  
  #Convert %mortality to mortality rate
  mu_agro_p<- -0.25*log(1-response_prawn[,2])
  mu_agro_p.sd<- -0.25*log(1-response_prawn[,3])*sqrt(5)
  
  #End chlorP response ################
    
get_Ro<-function(muPq, phi_Nq, parameters, eqbm)
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
  Eg<-parameters["Eg"]
  V<-parameters["V"]
  H<-parameters["H"]
  mu_H<-parameters["mu_H"]
  
  N<-eqbm$S+eqbm$E+eqbm$I
  P<-eqbm$P
  pred<- as.numeric((alpha*P)/(1+(alpha*N*Th)) )#death rate of snails due to predators (Prawns)
  
  T1<-beta*Eg*V
  T2<-(( 1-((mu_N+pred)/(f_N)))/(phi_N*(1/phi_Nq))) / (mu_N + sigma+pred)
  T3<- sigma/(mu_N+mu_I+pred)
  T4<-lamda/(mu_W+mu_H)
  
  Ro_miraCirc <- T1*T2*T3*T4
  
  Ro_miraCirc 
  
}

#Run with 10000 generated data points from predator mortality and carrying capacity parameters ###############

#phi_Nq variability derived from mesocosom data in Halstead_et_al_bottom_up_predict code
  phi_base<-rnorm(n=10000, mean=1, sd=(0.1488876*sqrt(10)))
  phi_fert<-rnorm(n=10000, mean=1.159642, sd=(0.1330912*sqrt(5)))
  phi_atra<-rnorm(n=10000, mean=1.614304, sd=(0.1223979*sqrt(5)))
  phi_atfe<-rnorm(n=10000, mean=1.509579, sd=(0.1711861*sqrt(5)))
#muPq data derived from halstead et al 2015 data to prawn mortality  
  muPq10k<-rnorm(n=10000, mean=mu_agro_p[64], sd=mu_agro_p.sd[64])
  
  
R0_base<-rep(0, length(phi_base))
for(i in 1:length(phi_base)){
  parameters["phi_Nq"]=phi_base[i]
  parameters["mu_Pq"]=muPq10k[i]
  output=as.data.frame(ode(nstart,time,schisto_halstead,parameters)) 
  eqbm<-output[365*yrs,]
  
  R0_base[i]<-get_Ro(muPq = muPq10k[i], phi_Nq = phi_base[i], 
                parameters=parameters, eqbm=eqbm)
}    
  
R0_base<-get_Ro(muPq = muPq10k, phi_Nq = phi_base, 
                parameters=parameters, eqbm=eqbm)
R0_fert<-get_Ro(muPq = muPq10k, phi_Nq = phi_fert, 
                parameters=parameters, eqbm=eqbm)
R0_atra<-get_Ro(muPq = muPq10k, phi_Nq = phi_atra, 
                parameters=parameters, eqbm=eqbm)
R0_atfe<-get_Ro(muPq = muPq10k, phi_Nq = phi_atfe, 
                parameters=parameters, eqbm=eqbm)

bot_up_r0<-data.frame("Treatment"=c('ChlorP Only', 'ChlorP+Fert', 'ChlorP+Atra', 'All three'),
                      "meanR0"=c(mean(R0_base), mean(R0_fert),mean(R0_atra),mean(R0_atfe)),
                      "st.erR0"=c(st.er(R0_base), st.er(R0_fert),st.er(R0_atra),st.er(R0_atfe)))

bot_up_r0$Treatment<-factor(bot_up_r0$Treatment, levels = c('ChlorP Only', 'ChlorP+Fert', 
                                                            'ChlorP+Atra', 'All three'))

ggplot(bot_up_r0, aes(x=Treatment, y=meanR0, fill=Treatment))+
  theme_bw()+
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=15))+
  scale_fill_manual(values=c('red', 'green', 'yellow', 'blue')) +
  ylab("Mean +/- SEM predicted R-0")+
  geom_bar(position=position_dodge(), stat="identity", width = .7) +
  geom_errorbar(aes(ymin=meanR0-st.erR0,
                    ymax=meanR0+st.erR0),
                width=.2, position=position_dodge(.7))
