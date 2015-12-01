#First version of agro model with all agro-imopacted functions incorporated.
require(deSolve)
######################################################################


#this function has 5 state variables, the 4 from Sokolow, 
#P, the prawn polulation and 
#the miracidia(M) and circariae(C) populations as a function of agrochemical concentration q
schisto_prawnsMiraCircFunc_agro=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    P=n[5]
    N=S+E+I
    
      
    #Agrochemical impacted parameters
    
    #1. Rate of successful hatching of schistosome eggs
    V<-function(V_q){
      V_o * V_q
    }
    
    #2. Infectivity of miracidia or circariae, as a fraction, includes rate of infection/mortality.
    R<-function(R_q){
      R_o * R_q
    }
    
    #3. Rate of shedding circariae by infected snails
    Theta<-function(Theta_q){
      Theta_o * Theta_q
    }
    
    #4. Attack rate of prawns
    alpha<-function(alpha_q){
      alpha_o*alpha_q
    }
    
    #Miracidia population
    M=0.5*W*H*m*V(V_q)*R(R_q)
    
    #Circariae population
    C=I*Theta(Theta_q)*R(R_q)
    
    #Rate of predation
    pred= (alpha(alpha_q)*P)/(1+(alpha(alpha_q)*N*Th)) #death rate of snails due to predators (Prawns)
    
    
    #snail compartment models
    dSdt= (f_N*f_Nq)*((1-((phi_N+phi_Nq)*N))*(S+(z*E))) - 
      ((mu_N+mu_Nq)*S) - (pred*S) - (beta*(M)*S)
    dEdt= beta*(M)*S - ((mu_N+mu_Nq+pred+sigma)*E)
    dIdt= (sigma*E) - ((mu_N+mu_Nq+pred+mu_I)*I)
    
    #worm burden in human
    dWdt= (lamda*C) - ((mu_W+mu_H)*W)
    
    #prawn population
    dPdt= (f_P*(1-(phi_P*P))*P)-((mu_P+mu_Pq)*P)
    
    return(list(c(dSdt,dEdt,dIdt,dWdt, dPdt)))
  }) 
} 

#Set time to run and number of starting prawns. P>=1 in order for prawn population model to populate
yrs=30
time=seq(0,365*yrs,1)
p=0


parameters=c(
  ##standard snail parameters 
  f_N=0.16, # birth rate
  phi_N=(1-1/(80*0.16))/10000, # carrying capacity
  z=0.5, #Proportion of exposed snails that reproduce
  mu_N=1/80, #Mortality rate
  beta=0.000004, #Infection rate by miracidia
  sigma=0.02, #Transition rate from exposed to infected
  mu_I=0.05, #additional snail death due to infection
  
  ## snail parameters impacted by agrochemicals
 
  f_Nq=1, #Scalar of snail birth rate by chemical concentration
  phi_Nq=0, #Scalar of snail carrying capacity by chemical concentration
  mu_Nq=0, #Extra mortality rate from chemical toxicity
  
  #agrochemical concentration
  #q=0,
  
  #prawn parameters
  
  alpha_o=0.003, #attack rate
  Th=0.1,#~Prawn predation limit
  f_P=0.128,#prawn birth rate
  phi_P=1/50,  #prawn carrying capacity
  mu_P= 0.00137, #prawn mortality
  
  #prawn parameter impacted by agrochemicals
  mu_Pq=0, #can be a function or fixed (0.02 for eg)
  alpha_q=1, #Scalar of predation rate due to sub-lethal toxicity
  
  #Adult Worm, Miracidia and Circariae Parameters
  lamda=1, #fraction of circariae that infect humans
  mu_W=1/(3*365), # death rate of adult worms
  m=0.5, #Miracidial production by adult female worms
  V_o=1, #Scalar to miracidial input due to agrochemical toxicity to schistosome eggs
  R_o=1, #Scalar to miracidial/cercarial input due to direct toxicity to larval schistosomes
  Theta_o=0.00004, #miracidia shedding rate divided by mortality
  
  #parameters impacted by agrochemicals
  V_q=1,#Egg viability (produced eggs that hatch to miracidia)
  R_q=1,#Larval mortality
  Theta_q=1,#Scalar to infected snail cercarial shedding rate based on agrochemical exposure
  
  #Human parameters
  H=1000, #number of humans
  mu_H=1/(60*365)
)

#end of parameters ###############################

#Set initial values
nstart=c(S=7000,E=3750,I=1200, W=20, P=p)
#Run model and save equilibrium values for 5 state variables
output=as.data.frame(ode(nstart,time,schisto_prawnsMiraCircFunc_agro,parameters)) 

eqbm<-output[365*yrs,]
eqbm

#Plot changes over time period #############################################
quartz()
par(mfrow=c(1,1))
plot(output[,1], output[,2], type='l', xlab="time",ylab="System Variables",
     main=paste("Prawns start=",p, "N=", output[365*yrs,2]+output[365*yrs,3]+output[365*yrs,4], sep=" "), sub = paste("inv of prawn carrying cap = ", parameters["phi_P"], sep=""),
     ylim=c(0,max( output[,2]+output[,3]+output[,4] )), col=1, lty=1, lwd=2)
lines(output[,1],output[,3],col=2, lty=2, lwd=2)
lines(output[,1],output[,4],col=3, lty=3, lwd=2)
lines(output[,1],output[,2]+output[,3]+output[,4],col=4, lty=4, lwd=2)
lines(output[,1],output[,5],col=5, lty=5)
lines(output[,1],output[,6],col=6, lty=6, lwd=2)

text(x=365*yrs, y=eqbm[2:6], labels=as.character(round(eqbm[2:6],2)), cex=0.5 )
legend("topright",c("S","E","I","N", "W", "P"),col=c(1:6),lty=1:6,cex=0.7, lwd=2)

#Calculate approximate Ro ###########################
#From Sokolow et al 2015: 
#Ro = (betha*msr*eta*H*((1-mu/f)/phi)/(mu+sig))*(sig/(mu+alfa))*(lambda/(mua+dh)); 
# need to recalculate for functional response therefore function returns Sokolow based Ro
## as well as Ro that incorporates Ag parameters 
get_Ro<-function(parameters, eqbm)
{
  f_N<-parameters["f_N"]
  phi_N<-parameters["phi_N"]
  z<-parameters["z"]
  mu_N<-parameters["mu_N"]
  beta<-parameters["beta"]
  sigma<-parameters["sigma"]
  mu_I<-parameters["mu_I"]
  f_Nq<-parameters["f_Nq"]
  phi_Nq<-parameters["phi_Nq"]
  mu_Nq<-parameters["mu_Nq"]
  #q<-parameters["q"]
  alpha_o<-parameters["alpha_o"]
  Th<-parameters["Th"]
  f_P<-parameters["f_P"]
  phi_P<-parameters["phi_P"]
  mu_P<-parameters["mu_P"]
  muPq<-parameters["mu_Pq"]
  alpha_q<-parameters["alpha_q"]
  lamda<-parameters["lamda"]
  mu_W<-parameters["mu_W"]
  m<-parameters["m"]
  V_o<-parameters["V_o"]
  R_o<-parameters["R_o"]
  Theta_o<-parameters["Theta_o"]
  V_q<-parameters["V_q"]
  R_q<-parameters["R_q"]
  Theta_q<-parameters["Theta_q"]
  H<-parameters["H"]
  mu_H<-parameters["mu_H"]
  
  
  Ro_sok = (beta*m*0.5*H*((1-mu_N/f_N)/phi_N)/(mu_N+sigma))*(sigma/(mu_N+mu_I))*
    (Theta_o/(mu_W+mu_H)); # need to recalculate for functional response
  print(paste("Ro using Sokolow model = ",Ro_sok, sep=" " ))
   
  #Agrochemical impacted parameters
  
  #1. Rate of successful hatching of schistosome eggs
  V<-function(V_q){
    V_o * V_q
  }
  
  #2. Infectivity of miracidia or circariae, as a fraction
  R<-function(R_q){
    R_o * R_q
  }
  
  #3. Rate of shedding circariae by infected snails
  Theta<-function(Theta_q){
    Theta_o * Theta_q
  }
  
  #4. Attack rate of prawns
  alpha<-function(alpha_q){
    alpha_o*alpha_q
  }
  N<-eqbm$S+eqbm$E+eqbm$I
  P<-eqbm$P
 pred<- as.numeric((alpha(alpha_q)*P)/(1+(alpha(alpha_q)*N*Th)) )#death rate of snails due to predators (Prawns)
  
  
  T1<-beta*0.5*H*m*V(V_q)*R(R_q)
  T2<- ( ( 1-((mu_N+mu_Nq+pred)/(f_N+f_Nq)) )/(phi_N+phi_Nq) ) / ( mu_N + mu_Nq+ sigma+pred)
  T3<- sigma/(mu_N+mu_Nq+mu_I+pred)
  T4<-lamda*Theta(Theta_q)*R(R_q)/(mu_W+mu_H)
  Ro_miraCirc <- T1*T2*T3*T4
    
    
    
    
    #(beta*((m_egg*p_egg2mira)/mu_M)*0.5*H*((1-mu_N/f_N)/phi)/(mu_N+sigma))*(sigma/(mu_N+mu_I))* ((beta_H*lamda_C*lamda_I)/(mu_C+(beta_H*H) ) )*(1/(mu_W+mu_H))
   print(paste("Ro with P=",P, "using our modified model = ",Ro_miraCirc, sep=" " ))
   
   Ro_miraCirc 
   
}

#Run Ro with base parameters and no agroC inputs #########################
get_Ro(parameters=parameters, eqbm=eqbm)
r<-as.numeric(get_Ro(parameters, eqbm))
r

# No agrochemicals, only Prawn carrying capacity variation. ##################

Ro_range<-numeric()
CC_range<-1:60
phi_P_range<-1/CC_range

for(phi_P in phi_P_range){
  nstart=c(S=7000,E=3750,I=1200, W=20, P=p)
  parameters["phi_P"]<-phi_P
  output=as.data.frame(ode(nstart,time,schisto_prawnsMiraCircFunc_agro,parameters)) #params_withextra
  #output[100:500,]
  eqbm<-output[365*yrs,]
  eqbm
  r<-as.numeric(get_Ro(parameters, eqbm))
  Ro_range<-c(Ro_range,r)
}

#plot this
quartz()
plot(CC_range, Ro_range, type='l', col=3, lwd=2)
abline(h=1, col=3, lty=2)


# Read in Chlorpyrifos data, modeling and functions for example agroC ###################
require(aod)
require(ggplot2)

#Probit analysis of toxicity data to prawns from Halstead et al 2015
chlor<-data.frame(dose=c(0,0.64,3.176863,6.331487,29.814438,58.029675), 
                  mortality=c(0,0,0,0,0.81,1))

tox_prawn<-glm(mortality ~ dose, family=binomial(link="probit"),data=chlor)
summary(tox_prawn) #mu_Pq

response_prawn<-data.frame(dose=seq(from=0, to=500, by=1))
response_prawn[, c("mortality", "se")]<-predict(tox_prawn, response_prawn, 
                                                type = "response", se.fit=TRUE)

#Plot to visualize fit to data (looks good, but few data points)
plot(response_prawn$dose, response_prawn$mortality, type='l') #mu_Pq
points(chlor$dose, chlor$mortality, pch=18, col=2)

#Reductions in predator attack rate (alpha) caused by sublethal doses of insecticide
##Actually modeled from a different organophosphate insecticide, Malathion, so ok to disclude this
prawn_alpha<-data.frame(dose=c(0, 10,50,100), a_red=c(0, 0.16,0.31,0.44))
plot(prawn_alpha$dose, prawn_alpha$a_red, type='l')

alpha<-lm(a_red ~ dose +0, data=prawn_alpha)
summary(alpha)

res_prawn_a<-data.frame(dose=seq(from=0, to=500, by=1))
res_prawn_a[, c("a_red", "se")]<-predict(alpha, res_prawn_a, 
                                         type = "response", se.fit=TRUE)

res_prawn_a$a_red[res_prawn_a$a_red>=1.00]<- 1.00 #alpha_Nq

#Snail mortality rate increase (i.e. mu_N_q) taken from Ibrahim 1992
##Modeling derived in excel
model_mu_q<-vector(mode="numeric", length=501)
for(i in seq(from=0, to=501, by=1)){
  model_mu_q[i]=0.0125*2.71828^(0.0028*i)-0.0125 #Exponential model derived from raw data in excel
}

snail_mu<-data.frame(dose=seq(from=0, to=500, by=1), mu_q=model_mu_q)

plot(snail_mu$dose, snail_mu$mu_q, type="l")
#Assess snail reductions in fecundity resulting from Chlorpyrifos exposure
#Dose=ppb chlorpyrifos, response=proportion reduction of control daily reproduction rate
snail_f<-data.frame(dose=c(0,125,250,500), f_red=c(0,(1-.73),(1-.48),(1-.37)))
plot(snail_f$dose, snail_f$f_red, type="l")

snail<-lm(f_red ~ dose + 0, data=snail_f)
response_snail<-data.frame(dose=seq(from=0, to=500, by=1))
response_snail[, c("f_red", "se")]<-predict(snail, response_snail, 
                                            type = "response", se.fit=TRUE) 

mu_chem<-rep(0, times=dim(response_snail)[1])

mu_chem<- (-0.25*log(1-agroChem_data[,2]) )#mu_Pq

#Combine all chlorpyrifos effects in one data frame
agroChem_data<-data.frame(dose=response_prawn$dose, prawn_mort=response_prawn$mortality, 
                          snail_birth=response_snail$f_red, alpha_red=res_prawn_a$a_red, 
                          snail_death=snail_mu$mu_q)#dose, 1-mu_Pq, 1-f_Nq, 1-alpha_q, mu_Nq.


#Plotting functions, etc. below ###################################


# effect of agrochemicals compiled in a matrix
Ro_matrix<-matrix(0, nrow=dim(agroChem_data)[1], ncol=length(phi_P_range))
#create a dataframe with the diff agrochem impacted parameters, keeping other ones fixed.
#for(i in 1:dim(agroChem_data)[1]){
for(i in 499:dim(agroChem_data)[1]){
  agroData<-agroChem_data[i,]
  agro_parameters<-parameters
  agro_parameters["mu_Pq"]<-as.numeric(mu_chem[i])
  agro_parameters["f_Nq"]<-as.numeric( 1- agroData[3])
  agro_parameters["alpha_q"]<-as.numeric(1-agroData[4])
  agro_parameters["mu_Nq"]<-as.numeric( agroData[5])
  
  print(i)
  for(phi_P in phi_P_range){
    j<-which(phi_P_range==phi_P)
    nstart=c(S=7000,E=3750,I=1200, W=20, P=p)
    agro_parameters["phi_P"]<-phi_P
    output=as.data.frame(ode(nstart,time,schisto_prawnsMiraCircFunc_agro,agro_parameters)) #params_withextra
    #output[100:500,]
    eqbm<-output[365*yrs,]
    eqbm
    r<-as.numeric(get_Ro(agro_parameters, eqbm))
    Ro_matrix[i,j]<-r
    #Ro_range<-c(Ro_range,r)
  }#end of phi_P range
  
  #save the file
  wd<-"/Users/Arathi/Documents/2015/Monash /Data/"
  outputfile<-paste(wd,"Ro_matrix.Rdata", sep="")
  save(Ro_matrix,file=outputfile)
  
}
#save the file
wd<-"/Users/Arathi/Documents/2015/Monash /Data/"
outputfile<-paste(wd,"Ro_matrix.Rdata", sep="")
#save(Ro_matrix,file=outputfile)
load(file=outputfile)

#2D plot of Carrying capacity vs Ro
quartz()
plot(CC_range, Ro_range, lwd=2, col=1, type='l', xlab="Prawn Carrying Capacity", ylab="Ro")
lines(CC_range, Ro_matrix[501,], lwd=2, col=2)
legend(x=35,y=4, legend=c("No Agrochemicals", "With Agrochemicals"), col=c(1,2), lwd=2)

#2D plot of chloropyrifos concentration versus Ro at a given carrying capacity (say 1, 20 or 50).
CC_fixed1<-1
CC_fixed2<-20
CC_fixed3<-40
CC_fixed4<-50
quartz()
plot(agroChem_data[,1], Ro_matrix[,CC_fixed1], col=1, type='l', xlab="Concentration of Chloropyrifos (ppm)", ylab="Ro", lwd=2 )
lines(agroChem_data[,1], Ro_matrix[,CC_fixed2], col=2, lwd=2)
lines(agroChem_data[,1], Ro_matrix[,CC_fixed3], col=3, lwd=2)
lines(agroChem_data[,1], Ro_matrix[,CC_fixed4], col=4, lwd=2)
legend(x=400, y=4, legend=c("0 ppm", "20 ppm", "40 ppm", "50 ppm"), col=1:4, lwd=2)

#3D line plot
library(plot3D)
#create the 3D data
data_3d <-numeric()
for(i in 1:dim(Ro_matrix)[1]){
  for(j in 1:dim(Ro_matrix)[2]){
    x<-agroChem_data[i,1]
    y<-CC_range[j]
    z<-Ro_matrix[i,j]
    row<-c(x,y,z)
    data_3d<-rbind(data_3d, row)
    
  }
}

dose0<-which(data_3d[,1]==0)
dose20<-which(data_3d[,1]==20)
dose500<-which(data_3d[,1]==500)
ignore<-which(data_3d[,3]==0) # not yet computed these entries
#save the 3D data
wd<-"/Users/Arathi/Documents/2015/Monash /Data/"
outputfile<-paste(wd,"data_3d.Rdata", sep="")
#save(data_3d,file=outputfile)
load(outputfile)
quartz()
lines3D(x=data_3d[-ignore,1], y=data_3d[-ignore,2], z=data_3d[-ignore,3], col=round(data_3d[-ignore,3]), xlab="dose(ppm)", ylab="carrying capacity of prawns", zlab="Ro")

lines3D(x=data_3d[dose20,1], y=data_3d[dose20,2], z=data_3d[dose20,3], col=round(data_3d[dose20,3]), xlab="dose(ppm)", ylab="carrying capacity of prawns", zlab="Ro")

scatter3D(x=agroChem_data[t_range,1], y=agroChem_output[t_range,1]+agroChem_output[t_range,2]+agroChem_output[t_range,3], z=agroChem_output[t_range,5], type='l', col=2, xlab="C", ylab="N", zlab="P")    

