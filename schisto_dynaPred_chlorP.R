# code to compute effect of a fixed concentration of agrochemical at eqbm. Get 3 d plot of concentration vs N and P.

#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.##########

require(deSolve)
require(aod)
require(ggplot2)
require(nlstools)
require(nls2)
require(graphics)

#Chlorpyrifos toxicity/functional responses 
#Probit analysis of toxicity data to prawns from Halstead et al 2015 ######################
chlor<-data.frame(dose=c(0,0.64,3.176863,6.331487,29.814438,58.029675), 
                  mortality=c(0,0,0,0,0.81, 1.00))

chlor$mu_P<- -0.25*log(1-chlor[,2])
chlor[6,3]=9.0109

#Probit analysis of data (probit analysis most often used for toxicological data investigating mortality responses to toxic exposures)
tox_prawn<-glm(mortality ~ dose, family=binomial(link="probit"),data=chlor)
summary(tox_prawn)

#Extrapolate response to constant gradient of Chlorpyrifos concentration
response_prawn<-data.frame(dose=seq(from=0, to=500, by=1))
response_prawn[, c("mortality", "se")]<-predict(tox_prawn, response_prawn, 
                                                type = "response", se.fit=TRUE)

#Convert %mortality to mortality rate
mu_agro_p<- -0.25*log(1-response_prawn[,2])

#Plot to visualize fit to data
plot(response_prawn$dose, response_prawn$mortality, type='l', 
     xlab=expression(paste("Chlorpyrifos (", mu, "g/L)")))
  points(chlor$dose, chlor$mortality, pch=18, col=2, cex=1.5)
  
  
#Reductions in predator attack rate (alpha) caused by sublethal doses of insecticide ###########
  #Data taken from Yuan et al 2004: prawn exposure to Malathion, a different organophosphate insecticide
prawn_alpha<-data.frame(dose=c(0, 10,50,100), a_red=c(0, 0.16,0.31,0.44))
  plot(prawn_alpha$dose, prawn_alpha$a_red, ylim=c(0,1), xlim=c(0,500), pch=18, col="red")
    x=c(1:500)
    lines(x, 0.1*log(x), type='l') #Rough log function appears to fit pretty well
#Test a linear regression to data  
alpha<-lm(a_red ~ dose +0, data=prawn_alpha)
  summary(alpha)

#Create data frame and fill with model predictions    
res_prawn_a<-data.frame(dose=seq(from=0, to=500, by=1))
  res_prawn_a[, c("a_red", "a_red_lwr", "a_red_upr")]<-predict(alpha, res_prawn_a, interval="confidence") 
  
res_prawn_a$a_red[res_prawn_a$a_red>=1.00]<- 1.00 # Proportional reduction therefore !>1.00

#Test a non-linear fit to data
alpha2<-nls(a_red ~ b*log((dose+1)), data=prawn_alpha, start=list(b=0.1), trace=TRUE)
  summary(alpha2)

#Add to predictions data frame  
res_prawn_a[, c("a_red2")]<-predict(alpha2, res_prawn_a)

#Plot to compare two proposed models to observed data
plot(prawn_alpha$dose, 1-prawn_alpha$a_red, ylim=c(0,1), xlim=c(0,500), col="red", 
     pch=18, cex=1.4, ylab="alpha_red", xlab=expression(paste("Chlorpyrifos (", mu, "g/L)")))
  lines(res_prawn_a$dose, 1-res_prawn_a$a_red)
    lines(res_prawn_a$dose, 1-res_prawn_a$a_red_lwr, lty=3)
    lines(res_prawn_a$dose, 1-res_prawn_a$a_red_upr, lty=3)
  lines(res_prawn_a$dose, 1-res_prawn_a$a_red2, col="blue")
  legend("topright", legend=c("linear model", "log model"), col=c(1,4), lwd=1)

  
  
#Snail mortality rate increase (i.e. mu_N_q) taken from Ibrahim 1992 #####################
mu_n_obs<-data.frame(dose=c(0,125,250,500), mu_n_q=c(0.0125,0.01526559, 0.01673515, 0.06557157))

#Model enhanced mortality as a function of dose. AgroC-free(i.e. control) mortality=0.0125  
  mu_n_mod<-nls(mu_n_q ~ 0.0125*exp(b*dose)-0.0125,data = mu_n_obs, start=list(b=0.002), trace=TRUE)
    summary(mu_n_mod)
#Fill data frame with modeled response
  model_mu_q<-data.frame(dose=seq(from=0, to=500, by=1))
    model_mu_q[,"mu_n_q"]<-predict(mu_n_mod, model_mu_q, 
                        type = "response", se.fit=TRUE)   
    
#Plot fit to data
  par(mfrow=c(1,1))
    plot(model_mu_q$mu_n_q, type='l', ylab="mu_n_q", xlim=c(0,500), ylim=c(0,0.1), 
         xlab=expression(paste("Chlorpyrifos (", mu, "g/L)")))
      points(mu_n_obs$dose, mu_n_obs$mu_n_q-0.0125, pch=18, col="red", cex=1.3)

snail_mu<-data.frame(dose=seq(from=0, to=500, by=1), mu_q=model_mu_q$mu_n_q)
#Assess snail reductions in fecundity resulting from Chlorpyrifos exposure from Ibrahim 1992############################
#Dose=ppb chlorpyrifos, response=proportion reduction of control daily reproduction rate
snail_f<-data.frame(dose=c(0,125,250,500), f_red=c(0,(1-.73),(1-.48),(1-.37)))

#Linear model fit to data
snail<-lm(f_red ~ dose + 0, data=snail_f)
  summary(snail)
#Create and fill model results in data frame  
response_snail<-data.frame(dose=seq(from=0, to=500, by=1))
response_snail[, c("f_red", "f_red_lwr", "f_red_upr")]<-predict(snail, response_snail, interval="confidence") 

#Logarithmic fit to data
snail2<-nls(f_red ~ b*log(dose+1), data=snail_f, start=list(b=0.1), trace=TRUE)
  summary(snail2)
#Add to data frame
  response_snail[,"f_red2"]<-predict(snail2, response_snail)
  
#Plot two models and assess fit to data
  plot(snail_f$dose, 1-snail_f$f_red, pch=18, col="red", cex=1.3, ylim=c(0,1),xlim=c(0,500),
     xlab=expression(paste("Chlorpyrifos (", mu, "g/L)")), ylab="f_n_red")
  lines(response_snail$dose, 1-response_snail$f_red)
    lines(response_snail$dose, 1-response_snail$f_red_lwr, lty=3) #lower 95%Ci
    lines(response_snail$dose, 1-response_snail$f_red_upr, lty=3) #upper 95%Ci
  lines(response_snail$dose, 1-response_snail$f_red2, col="blue") #log fit
  legend("topright", legend=c("linear model", "log model"), col=c(1,4), lwd=1)
  
#Combine all chlorpyrifos effects in one data frame #######################
agroChem_data<-data.frame(dose=response_prawn$dose, prawn_mort=mu_agro_p, 
                          snail_birth=response_snail$f_red, alpha_red=res_prawn_a$a_red, 
                          snail_death=snail_mu$mu_q)
#Combine all chlorpyrifos effects in one data frame with log functions instead of linear ones#######################  
agroChem_data<-data.frame(dose=response_prawn$dose, prawn_mort=mu_agro_p, 
                          snail_birth=response_snail$f_red2, alpha_red=res_prawn_a$a_red2, 
                          snail_death=snail_mu$mu_q)  

#Model with ChlorP impacts included #####################################################################


schisto_dynaPred_chlor=function(t, n, parameters) { 
  with(as.list(parameters),{
    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]
    P=n[5]
    N=S+E+I
    pred= (alpha*alpha_q*P)/(1+(alpha*alpha_q)*N*Th)
    inf= 0.5*beta*W*H*m*S
    dSdt= (f_N*f_red)*((1-(phi_N*N))*(S+(z*E))) - (mu_N*S) - (mu_N_q*S) - (pred*S) - inf
    dEdt= inf - ((mu_N+mu_N_q+pred+sigma)*E)
    dIdt= (sigma*E) - ((mu_N+mu_N_q+mu_I)*I) - (pred*I)
    dWdt= (lamda*I) - ((mu_W+mu_H)*W)
    dPdt= f_p*(1-(phi_p*P))*P-(mu_p + mu_p_q)*P
    return(list(c(dSdt,dEdt,dIdt,dWdt,dPdt)))
  }) 
} 
# Read in parameter and starting values ###################################
time=seq(0,365*30,1)
p=1 #p should be >=1 in order to populate prawn population

parameters=c(f_N=0.16, #Instantaneous snail mortality rate
             phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
             z=0.5, #Fraction of exposed snails that reproduce
             mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
             mu_N_q=0, #Enhanced snail mortality from agroC exposure
             alpha=0.003, #per capita prawn attack rate on snails
             Th=0.1, #Prawn predation max
             beta=0.000004, #Snail infection probability
             H=1000, # human population
             sigma=0.02, # exposed snails becoming infected
             mu_I=0.05, #Enhanced mortality of infected snails
             lamda=0.00004, #Human infection probability from infected snails
             mu_H=0.00004566, #Worm mortality from host deaths
             mu_W=0.00091324, #Adult worm natural mortality
             m=0.5, #Miracidial production per reproductive female
             alpha_q=1, #Reduction in prawn attack rate from AgroC
             f_p=0.128, #Prawn per capita daily reproduction rate
             phi_p=0.02, #Carrying capacity of prawn population
             mu_p= 1/(2*365), #Natural mortality rate of prawns (assuming lifespan=2 years)
             mu_p_q= 0, #Additional prawn mortality rate caused by agroC exposure
             f_red=1 #Fraction of snail reproductive rate remaining given exposure to doe of agrochemical
)


nstart=c(S=2435,E=4445,I=1422, W=59.3, P=p) #Prawn-free equilibrium values
output=as.data.frame(ode(nstart,time,schisto_dynaPred_chlor,parameters)) 
head(output)
eqbm<-output[max(time),2:6]
eqbm

#Plot results of model to observe changes over time #############################
t_range<-time
par(mfrow=c(1,1))
plot(output[t_range,1], output[t_range,2], type='l', xlab="conc",ylab="Number of Snails",
     ylim=c(0,max(output[t_range,2]+output[t_range,3]+output[t_range,4] )))
lines(output[t_range,1],output[t_range,3],col="red")
lines(output[t_range,1],output[t_range,4],col="blue")
lines(output[t_range,1],output[t_range,2]+output[t_range,3]+output[t_range,4],col="green")
lines(output[t_range,1],output[t_range,5],col="brown")
lines(output[t_range,1],output[t_range,6],col="purple")

legend("topright",paste(c("S","E","I","N","W","P"), as.integer(c(eqbm[,1], eqbm[,2], eqbm[,3],
                                                      (eqbm[,1]+eqbm[,2]+eqbm[,3]), eqbm[,4], eqbm[,5]))),
       col=c("black","red","blue","green", "brown", "purple"), lty=1,cex=0.7)

###Plot four parameters affected by Chlorpyrifos for poster ##################

par(mfrow=c(4,1),cex.lab=2, bg="white", fg="black",cex.axis=1.5, col.axis="black", col.lab="black")
mat<-matrix(c(1,1,1,1,
              2,2,2,2,
              3,3,3,3,
              4,4,4,4,4,4), nrow=9, byrow=TRUE)
layout(mat)

#PLOT1-Prawn mortality rate
par(mai=c(0.15,0.85,0.05,0.3))
plot(agroChem_data[,1], 1-agroChem_data[,3], type='l', col=6, lwd=3,ylim=c(0,1),
     xaxt='n', xlab='', ylab = "")
  mtext(expression(italic('f'[N*','~q])),side=2,cex=2, line=2, col=6)
  points(snail_f$dose, 1-snail_f$f_red, pch=18, cex=2.3, col="black")
  text(500,0.75, labels="A",col="black", cex=2.5)
  legend("bottomleft", pch=18, cex=2, legend="Observed data pts", 
         col="black", bty='n')
  axis(1,labels=FALSE)
  
#PLOT2-Prawn attack rate fraction 
par(mai=c(0.15,0.85,0.15,0.3))
plot(agroChem_data[,1], agroChem_data[,5], type='l', col=6, lwd=3,ylim=c(0,0.06),
     xaxt='n',xlab='',ylab = "")
  mtext(expression(italic(mu[N*','~q])),side=2,cex=2, line=2.2, col=6)
  points(mu_n_obs$dose, mu_n_obs$mu_n_q, pch=18, cex=2.3, col="black")
  text(500,0.02, labels="B",col="black", cex=2.5)
  axis(1,labels=FALSE)
  
#PLOT3-Snail reproduction rate fraction 
par(mai=c(0.2,0.85,0.15,0.3))
plot(agroChem_data[,1], mu_chem, type='l', col=6,lwd=3, ylim=c(0,10), 
     xaxt='n',xlab='',ylab = "")
  mtext(expression(italic(mu[P*','~q])),side=2,cex=2, line=2.3, col=6)
  points(x=chlor$dose, y=chlor$mu_P, pch=18, cex=2.3, col="black")
  text(500,5, labels="C", col="black", cex=2.5)
  axis(1,labels=FALSE)

#PLOT4-Snail mortality rate increase   
par(mai=c(0.75,0.85,0.05,0.3))
plot(agroChem_data[,1], 1-agroChem_data[,4], type='l', col=6, lwd=3,ylim=c(0,1),
     xlab=expression(paste("Chlorpyrifos (", mu, "g/L)")), ylab = "")
  mtext(expression(italic(alpha[q])),side=2,cex=2, line=2.3, col=6)
  points(prawn_alpha$dose, 1-prawn_alpha$a_red, pch=18, cex=2.3, col="black")
  text(500,0.5, labels="D",col="black", cex=2.5)

#NOTE: No reason to rerun below model (time-consuming); csv files of outputs are available
#ALSO: Be sure to take note of which agroChem_data file is being used: linear models or log models (see above)  
#Model equilibrium values of state variables across range of chlorpyrifos values###################
  #This is based on an assumption of a constant concentration of the agrochemical having a constant agrochemical
  #effect across the entire model run (30 yrs) very unrealistic, but this is just to get an idea of the potential effects
agroChem_output<-numeric()
for(i in 1:dim(agroChem_data)[1]){
  p_mu<- agroChem_data[i,2]
  f_r<- 1-agroChem_data[i,3]
  alph_r<- 1-agroChem_data[i,4]
  n_mu<-agroChem_data[i,5]
  
  p=1
  parameters=c(f_N=0.16, #Instantaneous snail mortality rate
               phi_N= (1-1/(80*0.16))/10000, #~Inverse of snail carrying capacity
               z=0.5, #Fraction of exposed snails that reproduce
               mu_N=1/80, #Snail mortality rate; assumed lifespan=12 weeks
               mu_N_q=n_mu, #Enhanced snail mortality from agroC exposure
               alpha=0.003, #per capita prawn attack rate on snails
               Th=0.1, #Prawn predation max
               beta=0.000004, #Snail infection probability
               H=1000, # human population
               sigma=0.02, # exposed snails becoming infected
               mu_I=0.05, #Enhanced mortality of infected snails
               lamda=0.00004, #Human infection probability from infected snails
               mu_H=0.00004566, #Worm mortality from host deaths
               mu_W=0.00091324, #Adult worm natural mortality
               m=0.5, #Miracidial production per reproductive female
               alpha_q=alph_r, #Reduction in prawn attack rate from AgroC
               f_p=0.128, #Prawn per capita daily reproduction rate
               phi_p=0.02, #Carrying capacity of prawn population
               mu_p= 1/(2*365), #Natural mortality rate of prawns (assuming lifespan=2 years)
               mu_p_q= p_mu, #Additional prawn mortality rate caused by agroC exposure
               f_red=f_r #Fraction of snail reproductive rate remaining given exposure to doe of agrochemical
               )
   
  nstart=c(S=2435,E=4445,I=1422, W=59.3, P=p) #VALUES USED HERE ARE PRAWN-FREE EQUILIBRIUM
  nstart[5]<-p
  output=as.data.frame(ode(nstart,time,schisto_dynaPred_chlor,parameters)) 
  eqbm_val<-output[length(time),2:6]
 agroChem_output<-rbind(agroChem_output, eqbm_val) 
  print(i)
}

#Plot equilibrium values of state variables across range of chlorpyrifos concentrations ##################
t_range<-1:501
par(mfrow=c(1,1))
plot(agroChem_data[t_range,1], agroChem_output[t_range,1], type='l',ylab="State variables eqbm values (t=30 yrs)",
     xlab=expression(paste("Chlorpyrifos (", mu, "g/L)")),
     ylim=c(0,max(agroChem_output[t_range,1]+agroChem_output[t_range,2]+agroChem_output[t_range,3])))
lines(agroChem_data[t_range,1],agroChem_output[t_range,2],col="red")
lines(agroChem_data[t_range,1],agroChem_output[t_range,3],col="blue")
lines(agroChem_data[t_range,1],agroChem_output[t_range,1]+agroChem_output[t_range,2]+agroChem_output[t_range,3],
      col="green")
lines(agroChem_data[t_range,1],agroChem_output[t_range,4],col="brown")
lines(agroChem_data[t_range,1],agroChem_output[t_range,5],col="purple")

legend("topright",c("S","E","I","N", "W", "P"),col=c("black","red","blue","green", "brown", "purple"),
       lty=1,cex=0.7)

#Combine results across chlorpyrifos concentrations in one data frame
chlor_effects<-data.frame(chlor=agroChem_data[c(1:501),1], 
                          N_snails=(agroChem_output[c(1:501),1]+agroChem_output[c(1:501),2]+agroChem_output[c(1:501),3]),
                          I_snails=agroChem_output[c(1:501),3],
                          S_snails=agroChem_output[c(1:501),1],
                          W_humans=agroChem_output[c(1:501),4],
                          N_preds=agroChem_output[c(1:501),5])

#Plot preds across concentrations
plot(chlor_effects$chlor, chlor_effects$N_preds, type='l', xlab="chlor_conc", ylab="preds")

#Save files for later use ################
#Linear models used
write.csv(chlor_effects, 
          "E:/Remais_Work_7_14_15/Senegal_Schisto/AgroData/Data/chlor_tox.csv", 
          row.names = FALSE)

#Log models used
write.csv(chlor_effects, 
          "E:/Remais_Work_7_14_15/Senegal_Schisto/AgroData/Data/chlor_tox_log.csv", 
          row.names = FALSE)

#Call saved data to use in comparisons ########################
chlor_data<-read.csv("E:/Remais_Work_7_14_15/Senegal_Schisto/AgroData/Data/chlor_tox.csv")
chlor_data2<-read.csv("E:/Remais_Work_7_14_15/Senegal_Schisto/AgroData/Data/chlor_tox_log.csv")

#Plot to compare results of using liog vs linear response functions for alpha_red and f_n_red #################
plot(chlor_data$chlor, chlor_data$N_snails, type = 'l', col= "green", lwd=2,
     xlab = "Chlorpyrifos (ppb)", ylab = "State variables eqbm values")
  lines(chlor_data$chlor, chlor_data$I_snails, col= "blue", lwd=2)
  lines(chlor_data$chlor, chlor_data$E_snails, col= "red", lwd=2)
  lines(chlor_data$chlor, chlor_data2$I_snails, col= "blue", lwd=2,lty=2)
  lines(chlor_data$chlor, chlor_data2$E_snails, col= "red", lwd=2,lty=2)
  lines(chlor_data$chlor, chlor_data2$N_snails, col= "green", lwd=2,lty=2)
  legend("topright", legend=c("linear models", "log models"), lty=c(1,2))
  legend(x=437, y=7118, legend=c("N", "E", "I"), col=c("green", "red", "blue"), lwd=1.2)


#Plot both snail population and predator population across chlor concentration ###########################
par(mar=c(5, 4, 4, 6) + 0.1)
plot(chlor_data$chlor, chlor_data$N_snails, axes=F, ylim=c(0,8000), xlab="",
     ylab="",type='l',col="blue", main="")
  axis(4, ylim=c(0,8000),col="blue", las=1)
  mtext("Preds", side=2, line=2.5, col="red")
  
par(new=TRUE)

plot(chlor_data$chlor, chlor_data$N_preds, type = 'l', col= "red",
     xlab = "", ylab = "", ylim=c(0,50), axes=FALSE)
  axis(2, ylim=c(0,50), col="red", las=1)
  mtext("Snails", side=4, col="blue", line=4)
  
  axis(1, xlim=c(0,500), col="black")
  mtext("Chlorpyrifos (ppb)", side=1, col= "black", line=2.5)
