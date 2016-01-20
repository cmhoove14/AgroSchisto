# This code generates the sensitivity analysis of our schisto model.
# it uses LHC PRCC technique (Sally Blower 1991,1994) to identify the 
# agro impacted parameters that most greatly impact the outputs Ro and W.

######The model adapted from Sokolow 2015, to incorporate dynamic prawn population 
##### and agrochem impact.
##### Agrochem impacts the following parameters:fNq, phi_Nq, mu_Nq, alpha_q, mu_Pq, Theta_q, R_q, V_q.
require(deSolve)

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
    
    
    #snail compartment model
    dSdt= (f_N*f_Nq)*((1-((phi_N+phi_Nq)*N))*(S+(z*E))) - ((mu_N+mu_Nq)*S) - (pred*S) - (beta*(M)*S)
    dEdt= beta*(M)*S - ((mu_N+mu_Nq+pred+sigma)*E)
    dIdt= (sigma*E) - ((mu_N+mu_Nq+pred+mu_I)*I)
    
    #worm burden in human
    dWdt= (lamda*C) - ((mu_W+mu_H)*W)
    
    #prawn population
    dPdt= (f_P*(1-(phi_P*P))*P)-((mu_P+mu_Pq)*P)
    
    return(list(c(dSdt,dEdt,dIdt,dWdt, dPdt)))
  }) 
} 

yrs=30
time=seq(0,365*yrs,1)
p=1
CC_P<-1 #prawn carrying capacity

parameters=c(
  ##standard snail parameters 
  f_N=0.16,
  phi_N=(1-1/(80*0.16))/10000, 
  z=0.5, 
  mu_N=1/80, 
  beta=0.000004,
  sigma=0.02, 
  mu_I=0.05, #additional snail death due to infection
  
  ## snail parameters impacted by agrochemicals
  
  f_Nq=1,
  phi_Nq=0,
  mu_Nq=0,
  
  #agrochemical concentration
  #q=0,
  
  #prawn parameters
  
  alpha_o=0.003, #attack rate
  Th=0.1,
  f_P=0.128,#prawn birth rate
  phi_P=1/CC_P,  #carrying capacity
  mu_P= 0.00137, #prawn mortality
  
  #prawn parameter impacted by agrochemicals
  mu_Pq=0, #can be a function or fixed (0.02 for eg)
  alpha_q=1,
  
  #Adult Worm, Miracidia and Circariae Parameters
  lamda=1, #fraction of circariae that infect humans
  mu_W=1/(3*365), # death rate of worm
  m=0.5,
  V_o=1,
  R_o=1,
  Theta_o=0.00004,
  #parameters impacted by agrochemicals
  V_q=1,
  R_q=1,
  Theta_q=1,
  
  #Human parameters
  H=1000, #number of humans
  mu_H=1/(60*365)
)

#end of parameters


nstart=c(S=7000,E=3750,I=1200, W=20, P=p)
output=as.data.frame(ode(nstart,time,schisto_prawnsMiraCircFunc_agro,parameters)) 

eqbm<-output[365*yrs,]
eqbm

#############################################
#Plot changes over time period
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
  
  
  #   Ro_sok = (beta*m*0.5*H*((1-mu_N/f_N)/phi_N)/(mu_N+sigma))*(sigma/(mu_N+mu_I))*(Theta_o/(mu_W+mu_H)); # need to recalculate for functional response
  #   print(paste("Ro using Sokolow model = ",Ro_sok, sep=" " ))
  #    
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
  N<-eqbm["S"]+eqbm["E"]+eqbm["I"]
  P<-eqbm["P"]
  pred<- as.numeric((alpha(alpha_q)*P)/(1+(alpha(alpha_q)*N*Th)) )#death rate of snails due to predators (Prawns)
  
  
  T1<-beta*0.5*H*m*V(V_q)*R(R_q)*N
  T2<- sigma*lamda*Theta(Theta_q)*R(R_q)
  T3<-(mu_W+mu_H)*( mu_N + mu_Nq+ sigma+pred)*(mu_N+mu_Nq+mu_I+pred)
  
  #   print(T1*T2)
  #   print(T3)
  #   print(T1*T2/T3)
  Ro_miraCirc <- sqrt( (T1*T2)/T3 ) 
  
  
  
  
  #(beta*((m_egg*p_egg2mira)/mu_M)*0.5*H*((1-mu_N/f_N)/phi)/(mu_N+sigma))*(sigma/(mu_N+mu_I))* ((beta_H*lamda_C*lamda_I)/(mu_C+(beta_H*H) ) )*(1/(mu_W+mu_H))
  #  print(paste("Ro with P=",P, "using our modified model = ",Ro_miraCirc, sep=" " ))
  
  Ro_miraCirc 
  
}


get_Ro(parameters=parameters, eqbm=eqbm)
r<-as.numeric(get_Ro(parameters, eqbm))
r


######################################################################
#### SENSITIVITY ANALYSIS ############################################
######################################################################

CC_P<-1 #prawn carrying capacity

#test the variation in output variables W and Ro with change in agrochem parameters
# Agrochem impacts the following parameters:f_Nq, phi_Nq, mu_Nq, alpha_q, mu_Pq, Theta_q, R_q, V_q.
epsilon<-0.1*parameters["phi_N"]
f_Nq<-seq(from=0.001, to=1, by=0.001)#1
phi_Nq<-seq(from=(-parameters["phi_N"]+epsilon)/1000, to=(-parameters["phi_N"]+epsilon), by=(-parameters["phi_N"]+epsilon)/1000)
mu_Nq<-seq(from=0.00006, to=0.06, by=0.00006)
alpha_q<-seq(from=0.001, to=1, by=0.001)
mu_Pq<-seq(from=0.01, to=10, by=0.01)
Theta_q<-seq(from=0.001, to=1, by=0.001)
R_q<-seq(from=0.001, to=1, by=0.001)
V_q<-seq(from=0.001, to=1, by=0.001)#8



#create a matrix of indices for the LHC
LHC_indices<-matrix(0, nrow=length(f_Nq), ncol=8)

for(j in 1:dim(LHC_indices)[2]){
  LHC_indices[,j]<-sample(1:length(f_Nq), size=length(f_Nq), replace=FALSE)
}#end of j loop

params00<-cbind(f_Nq=f_Nq[LHC_indices[,1]], phi_Nq=phi_Nq[LHC_indices[,2]], mu_Nq=mu_Nq[LHC_indices[,3]], alpha_q=alpha_q[LHC_indices[,4]], mu_Pq=mu_Pq[LHC_indices[,5]], Theta_q=Theta_q[LHC_indices[,6]], R_q=R_q[LHC_indices[,7]], V_q=V_q[LHC_indices[,8]])
#params00<-expand.grid( f_Nq=f_Nq, phi_Nq=phi_Nq, mu_Nq=mu_Nq, alpha_q=alpha_q, mu_Pq=mu_Pq, Theta_q=Theta_q, R_q=R_q, V_q=V_q)

sims = dim(params00)[1]

constantparams<-cbind(
  f_N=rep(0.16, sims),
  phi_N=rep((1-1/(80*0.16))/10000, sims), 
  z=rep(0.5, sims), 
  mu_N=rep(1/80, sims), 
  beta=rep(0.000004, sims),
  sigma=rep(0.02, sims), 
  mu_I=rep(0.05,sims), #additional snail death due to infection
  
  ## snail parameters impacted by agrochemicals
  
  #f_Nq=rep(1, sims),
  #phi_Nq=rep(0, sims),
  #mu_Nq=rep(0, sims),
  
  #agrochemical concentration
  #q=0,
  
  #prawn parameters
  
  alpha_o=rep(0.003, sims), #attack rate
  Th=rep(0.1, sims),
  f_P=rep(0.128, sims),#prawn birth rate
  phi_P=rep(1/CC_P, sims),  #carrying capacity
  mu_P= rep(0.00137, sims), #prawn mortality
  
  #prawn parameter impacted by agrochemicals
  #mu_Pq=rep(0, sims), #can be a function or fixed (0.02 for eg)
  #alpha_q=rep(1, sims),
  
  #Adult Worm, Miracidia and Circariae Parameters
  lamda=rep(1, sims), #fraction of circariae that infect humans
  mu_W=rep(1/(3*365), sims), # death rate of worm
  m=rep(0.5, sims),
  V_o=rep(1, sims),
  R_o=rep(1, sims),
  Theta_o=rep(0.00004, sims), #miracidia shedding rate divided by mortality
  
  #parameters impacted by agrochemicals
  #V_q=rep(1, sims),
  #R_q=rep(1, sims),
  #Theta_q=rep(1, sims),
  
  #Human parameters
  H=rep(1000, sims), #number of humans
  mu_H=rep(1/(60*365), sims)
  
  
)
params1 <- cbind(params00,constantparams)
head(params1)

yrs=50
time=seq(0,365*yrs,1)

p<-1
op<-numeric()

for(i in 1:sims){
  parameters<-params1[i,]
  nstart=c(S=7000,E=3750,I=1200, W=20, P=p)
  output=as.data.frame(ode(nstart,time,schisto_prawnsMiraCircFunc_agro,parameters)) 
  eqbm<-output[365*yrs,]
  eqbm
  
  r<-as.numeric(get_Ro(parameters, eqbm))
  opRow<-cbind(eqbm,r)
  op<-rbind(op,opRow)
  print(c(i,r, eqbm$W))
}
wd<-"/Users/Arathi/Documents/2015/Monash /Data/"
outputfile<-paste(wd,"varyAgroParam_CCP", CC_P,".Rdata", sep="")
save(params00,file=outputfile)

outputfile<-paste(wd,"Output_varyAgroParam_CCP", CC_P,".Rdata", sep="")
save(op,file=outputfile)

op_mx<-cbind(params00,op[,5], op[,7])
NA_index<-which(is.na(op[,5])==TRUE)
if(length(NA_index)>0) op_mx<-op_mx[-NA_index,]

NA_index<-which(is.na(op[,7])==TRUE)
if(length(NA_index)>0) op_mx<-op_mx[-NA_index,]

####Get the basic stats of the variables and outputs
op_stats_mx<-matrix(0, nrow=5, ncol=dim(op_mx)[2])
rownames(op_stats_mx)<-c("Mean", "Median", "SD", "Q1", "Q3")

for(j in 1:dim(op_stats_mx)[2]){
  colValues<-op_mx[,j]
  qnt<-quantile(colValues, probs = c(0.25, 0.5, 0.75), na.rm = FALSE,
                names = FALSE, type = 1)
  op_stats_mx[1,j]<-mean(colValues)
  op_stats_mx[2,j]<-qnt[2]
  op_stats_mx[3,j]<-sqrt(var(colValues))
  op_stats_mx[4,j]<-qnt[1]
  op_stats_mx[5,j]<-qnt[3]
}

SA_W_R<-matrix(0, nrow=4, ncol=dim(params00)[2])

for(j in 1:dim(SA_W_R)[2]){
  c1<-cor.test(op_mx[,9], op_mx[,j],
               alternative = "two.sided",
               method = "spearman",
               exact = NULL, conf.level = 0.95, continuity = TRUE)
  
  SA_W_R[1,j]<-c1$estimate
  SA_W_R[2,j]<-c1$p.value
  
  c2<-cor.test(op_mx[,10], op_mx[,j],
               alternative = "two.sided",
               method = "spearman",
               exact = NULL, conf.level = 0.95, continuity = TRUE)
  SA_W_R[3,j]<-c2$estimate
  SA_W_R[4,j]<-c2$p.value
}#end of j loop

colnames(SA_W_R)<- colnames(params00)

order_ind1<-order(abs(SA_W_R[1,]), decreasing=TRUE)
order_ind2<-order(abs(SA_W_R[3,]), decreasing=TRUE)

quartz()
par(mfrow=c(1,1))
barplot(SA_W_R[1,], horiz=TRUE, col="cyan", xlab="Sensitivity to W, mean worm burden", xlim=c(-1,1), cex.names=0.5)
quartz()
par(mfrow=c(1,1))
barplot(SA_W_R[3,], horiz=TRUE, col="green", xlab="Sensitivity to Ro", xlim=c(-1,1), cex.names=0.5)


##########################################################


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

#Plot to visualize fir to data
plot(response_prawn$dose, response_prawn$mortality, type='l') #mu_Pq
points(chlor$dose, chlor$mortality, pch=18, col=2)

#Reductions in predator attack rate (alpha) caused by sublethal doses of insecticide
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

#Combine all chlorpyrifos effects in one data frame
agroChem_data<-data.frame(dose=response_prawn$dose, prawn_mort=response_prawn$mortality, 
                          snail_birth=response_snail$f_red, alpha_red=res_prawn_a$a_red, 
                          snail_death=snail_mu$mu_q)#dose, 1-mu_Pq, 1-f_Nq, 1-alpha_q, mu_Nq.

mu_chem<-rep(0, times=dim(response_snail)[1])

mu_chem<- (-0.25*log(1-agroChem_data[,2]) )#mu_Pq



##########################################################
CC_P_range<-1:60
phi_P_range<-1/CC_P_range

# effect of agrochemicals.
Ro_matrix<-matrix(0, nrow=dim(agroChem_data)[1], ncol=length(phi_P_range))
#create a dataframe with the diff agrochem impacted parameters, keeping other ones fixed.
for(i in 1:dim(agroChem_data)[1]){
  #for(i in 419:498){
  agroData<-agroChem_data[i,]
  agro_parameters<-parameters
  agro_parameters["mu_Pq"]<-as.numeric(mu_chem[i])
  agro_parameters["f_Nq"]<-as.numeric( 1- agroData[3])
  agro_parameters["alpha_q"]<-as.numeric(1-agroData[4])
  agro_parameters["mu_Nq"]<-as.numeric( agroData[5])
  
  print(i)
  for(phi_P in phi_P_range[c(1,16,17,48,60)]){
    
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
  outputfile<-paste(wd,"Ro_matrix_1.Rdata", sep="")
  save(Ro_matrix,file=outputfile)
  
}
#save the file
wd<-"/Users/Arathi/Documents/2015/Monash /Data/"
outputfile<-paste(wd,"Ro_matrix_1.Rdata", sep="")
#save(Ro_matrix,file=outputfile)
load(file=outputfile)

#2D plot of Carrying capacity vs Ro
AgroConc<-c(seq(from=0, to=50, by=10), seq(from=100, to=500, by=100) )
#AgroConc<-seq(from=20, to=30, by=1)
quartz()
plot(CC_range, Ro_matrix[AgroConc[1]+1,], lwd=2, col=1, type='l', xlab="Prawn Carrying Capacity", ylab="Ro", lty=1)
for(i in 2:length(AgroConc)){
  lines(CC_range, Ro_matrix[AgroConc[i]+1,], lwd=2, col=i, lty=i)
  
}
legend(x=35,y=4, legend=as.character(AgroConc), col=1:length(AgroConc), lwd=2, lty=1:length(AgroConc))

#2D plot of chloropyrifos concentration versus Ro at a given carrying capacity (say 1, 20 or 50).
#CC_fixed<-seq(from=1, to=51, by=5)
CC_fixed<-c(1,16,17,48,60)
row_range<-1:500
quartz()
plot(agroChem_data[row_range,1], Ro_matrix[row_range,CC_fixed[1]], col=1, type='l', xlab="Concentration of Chloropyrifos (ppm)", ylab="Ro", lwd=2, ylim=c(0,5) )
for(i in 2:length(CC_fixed)){
  lines(agroChem_data[row_range,1], Ro_matrix[row_range,CC_fixed[i]], col=i, lwd=2)
}
legend(x=300, y=4, legend=as.character(CC_fixed), col=1:length(CC_fixed), lwd=2)

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




