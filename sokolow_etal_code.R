#Susanne Sokolow Shistosomiasi model
#version 1.7
#
#(R code edited by Giulio A. De Leo)
#
#

#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

require(deSolve)


params <- c(

 P = 0,			# prawn density
 H = 1000, 		#human population abundance

#  snail dynamics
 f 	= 0.16, 	# instantaneous intrinsic fertility rate of snails times probability of survival to juvenile(day^-1)
 phi= (1-1/(80*0.16))/10000,	# density-dependent parameter for fertility (assuming cannibalism)
 mu	= 1/80, 	# snail mortality (mean life expectancy: 12 weeks)
 ze = .5, 		# fraction of exposed snails that reproduce

# epidemiological parameters for infected snails
 betha 	= (4*10^-6),	# per capita infection rate
 k		= 0.25,		# clumping parameter of teh negative binomila distribution
 msr	= 0.5,		# miracidia shedding rate per reproductive female (devided by miracidia mortality rate)
 sig	= 1/50,		# rate from exposed to infected
 alfa 	= 1/20,		# mortality rate of infected snail (if negative means that castrated snail increase their life expectancy - tyet it could be zero)

## predation rates
 q	= 0.003,		# predation rate on susceptible snails
 Th = 0.1,


## parameters for human host
 lambda	= 0.00004,		# cercaria shadding rate divided by cercaria mortality
 mua  	= 1/(3*365), # adult worm natural mortality =1/3 years
 dh		= 1/(60*365)# mortality/fertility rate of human population

) # end parameters

params

#######################################################################
#######################################################################


times <- seq(from=0,to=365*30,by=1) # time horizon seq(from=0,to=365*30,by=1) 
xstart <- c(M=20,S=3892, E=3750, I=1200) # initial conditions c(M=63,S=3892, E=3750, I=1200) which is prawn free equilibrium with these parameters

#######################################################################
#######################################################################

## MDA function

mda = function(t) {
    ifelse(t %% 365 == 0, 1, 0)}

###########################################################################
const.prawn.pop.model <- function (t, x, params) {

## first extract the state variables
M <- x[1]	# mean parasite burden in humans 
S <- x[2]	# density of susceptible snails
E <- x[3]	# density of infected but not yet infective snails
I <- x[4]	# density of infected snails



## now extract the parameters
P 	<- params["P"]	# prawn density
H	<- params["H"]		#human population abundance

#  snail dynamics
f 	<- params["f"]  	# instantaneous intrinsic fertility rate of snails
phi <- params["phi"] 	# density-dependent parameter for fertility (assuming cannibalism)
mu	<- params["mu"]		# snail mortality
ze 	<- params["ze"]		# fraction of exposed snails that reproduce

# epidemiological parameters for infected snails
betha <- params["betha"]	# per capita infection rate
k	<- params["k"]		# clumping parameter of teh negative binomila distribution
msr	<- params["msr"]	# miracidia shading rate per reproductive female
sig	<- params["sig"]	# rate from exposed to infected
alfa<- params["alfa"]	# mortality rate of infected snail (if negative means that castrated snail increas their life expectancy - tyet it could be zero)

## predation rates
q	<- params["q"]		# predation rate on susceptible smails

Th  <- params["Th"]		# handling time parameter from functional response


## parameters for the adult worms dynamics in the human hosts
lambda<-params["lambda"]	# cercarie shadding rate divided by cercarie mortality
mua	<- params["mua"]	# adult warm natural mortality =1/4 
dh	<- params["dh"]		# human mortality/fertility rate =1/60 

##compute the number of reproductive females (rf) assuming:
#  - 1:1 sex ratio 
#  - negative binomial distribution with clumping parameter k
# "rf" is equal to one half (for the 1:1 sex ratio, as the parasite is dioecious) the total number of parasites (M*H) 
# minus "single-worm" infections. Single worm infection are equal to the the number of human hosts with exactly one (1) parasite
# which can be computed analytically (see recursive formula 3.105, pag 87, Hilborn and Mangel, The ecological Detective!) 

# Accordingly, we could compute "eta", the fraction of reproductive females over the total parasite abundance, as follows:
    eta <- 0.5  # *(1 - 1/(1+M/k)^(k+1)) 
# yet, at invasion at least, we assume that an individual harboring two parasites (one male and one female) 
# is introduced in a population of H susceptible individuals and "1/phi" susceptible snails (assuming no prawns). 
# To compute threshold invasion, i.e. Ro in this case, I think that it is correct to assume that eta=0.5
# infact, at the onset of the infection, with only two parasite, the distribution cannot be certainly binomial.
# In conclusion, the formulation reported above is corrected as long as we assume that, 
# for instance, the mean parasite burden M is equal to or grater than 1, with M*H>>1 the total number of parasites



rF 	<- eta*M*H 


# total snail population
N = S+E+I

## now code the model equations
dSdt <- f*(1 - phi*N)*(S+ze*E) -          mu*S  - betha*msr*rF*S - (q*S*P/(1+q*N*Th))
dEdt <- betha*msr*rF*S         -  (mu + sig)*E - (q*P*E/(1+q*N*Th))
dIdt <- sig*E                  - (mu + alfa)*I - (q*P*I/(1+q*N*Th))
dMdt <- lambda*I               -  (mua + dh)*M 







## combine results into a single vector
dxdt <- c(dMdt, dSdt, dEdt, dIdt)
## return result as a list!
list(dxdt)
}
##########################################################################################

## now extract the parameters
P 	<- params["P"]	# prawn density
H	<- params["H"]		#human population abundance

#  snail dynamics
f 	<- params["f"]  	# instantaneous intrinsic fertility rate of snails
phi <- params["phi"] 	# density-dependent parameter for fertility (assuming cannibalism)
mu	<- params["mu"]		# snail mortality
ze 	<- params["ze"]		# fraction of exposed snails that reproduce

# epidemiological parameters for infected snails
betha <- params["betha"]	# per capita infection rate
k	<- params["k"]		# clumping parameter of teh negative binomila distribution
msr	<- params["msr"]	# miracidia shading rate per reproductive female
sig	<- params["sig"]	# rate from exposed to infected
alfa<- params["alfa"]	# mortality rate of infected snail (if negative means that castrated snail increas their life expectancy - tyet it could be zero)

## predation rates
q	<- params["q"]		# predation rate on susceptible smails
Th  <- params["Th"]		# handling time parameter from functional response


## parameters for the adult worms dynamics in the human hosts
lambda<-params["lambda"]	# cercarie shadding rate divided by cercarie mortality
mua	<- params["mua"]	# adult warm natural mortality =1/4 
dh	<- params["dh"]		# human mortality/fertility rate =1/60 

xstart
M <- xstart["M"]
S <- xstart["S"]
E <- xstart["E"]
I <- xstart["I"]

eta <- 0.5  #*(1 - 1/(1+M/k)^(k+1)) 
rF 	<- eta*M*H 
N = S+E+I

dSdt <- f*(1 - phi*N)*(S+ze*E) -          mu*S - q*S*P/(1+q*N*Th) - betha*msr*rF*S
dEdt <- betha*msr*rF*S         -  (mu + sig)*E - q*P*E/(1+q*N*Th)
dIdt <- sig*E                  - (mu + alfa)*I - q*P*I/(1+q*N*Th)
dMdt <- lambda*I               -  (mua + dh)*M 





Ro = (betha*msr*eta*H*((1-mu/f)/phi)/(mu+sig+))*(sig/(mu+alfa))*(lambda/(mua+dh)); # need to recalculate for functional response
y$Ro = (y$betha*y$msr*0.5*1000*((1-y$mu/y$f)/y$phi)/(y$mu+y$sig))*(y$sig/(y$mu+y$alfa))*(y$lambda/(y$mua+y$dh));
y3$Ro = (y3$betha*y3$msr*0.5*1000*((1-y3$mu/y3$f)/y3$phi)/(y3$mu+y3$sig))*(y3$sig/(y3$mu+y3$alfa))*(y3$lambda/(y3$mua+y3$dh));

##without prawns: 
Ro = (betha*msr*eta*H*((1-mu/f)/phi)/(mu+sig))*(sig/(mu+alfa))*(lambda/(mua+dh)); # need to recalculate for functional response
cat("this is Ro:", Ro, "\n")
1/Ro

betha*msr*eta*H*((1-mu/f)/phi)/(mu+sig)
sig/(mu+alfa)
lambda/(mua+dh)

0.5*(1 - 1/(1+5/(k))^(k+1))

(1-mu/f)

#################################################################
# now integrate the system
out <- as.data.frame(
 ode(
   func		= const.prawn.pop.model,
   y		= xstart,
   times	= times,
   parms	= params,
  )
)

## or using lsoda:

out <- as.list(data.frame(lsoda(xstart,times, const.prawn.pop.model, params)))

###################################################################
# plot the results over the entire simulation time
par(mfrow=c(3,1))
time_year = out$time/365
plot(    S~time_year, data=out,type='l',col='green', lwd = 2, ylim=c(0,max(c(S,E,I))), ylab = "Susceptible snails")
  points(E~time_year, data=out,type='l', ylab= "Exposed snails", lwd = 2, col = "blue")
  points(I~time_year, data=out,type='l', ylab= "Infected snails", lwd = 2, col = "red")
plot(M~time_year, data=out,type='l',col='green', lwd = 2,  ylab = "adult worm per host")
plot(I~M, data=out,type='l',col='green', lwd = 2,  ylab = "phase plan")
###################################################################
# now plot the results over the first year (365 days)
str(out)
tail(out)
par(mfrow=c(3,1))
plot(  S~time, data=out,type='l',col='green', lwd = 2, ylim=c(0,max(c(S,E,I))), ylab = "Susceptible snails", xlim=c(0,365))
points(E~time, data=out,type='l', ylab= "Exposed snails", lwd = 2, col = "blue", xlim=c(0,365))
points(I~time, data=out,type='l', ylab= "Infected snails", lwd = 2, col = "red", xlim=c(0,365))

plot(M~time, data=out,type='l',col='green', lwd = 2,  ylab = "adult worm per host", xlim=c(0,365))

plot(I~M, data=out,type='l',col='green', lwd = 2,  ylab = "phase plan",)



N=out$S + out$E + out$I
tail(N)
buf2=as.matrix(cbind(out$S/N, out$E/N, out$I/N), nrow=3)
str(buf2)
summary(buf2)
tail(buf2)

#prevalence of infected humans
MM=10; MM
k
#create a function for prevalence to use with Sapply() on the M vector
Prevalence <- function(x) {
	k=as.numeric(params["k"])
	p=1 - (1/(1+x/k)^(k))*(1+x/(1+x/k)) # fraction of humans with at least 2 parasites
	return(p)
}
1/((1+MM/k)^k)

# prevalence of humans infected with at least 2 parasites
1/(1+MM/k)^(k+1) 

#######################################################################
#######################################################################


times <- seq(from=0,to=365*30,by=1) # time horizon seq(from=0,to=365*30,by=1) 
xstart <- c(M=20, U=20, S=3892, E=3750, I=1200) # initial conditions c(M=63,S=3892, E=3750, I=1200) which is prawn free equilibrium with these parameters

#######################################################################
#######################################################################
#####################################################################################
const.prawn.pop.model.2pops <- function (t, x, params, co) {

## first extract the state variables
M <- x[1]	# mean parasite burden in the portion of the population treated  
S <- x[2]	# density of susceptible snails
E <- x[3]	# density of infected but not yet infective snails
I <- x[4]	# density of infected snails
U <- x[1]	# mean parasite burden in the portion of the Untreated population   



## now extract the parameters
P 	<- params["P"]		# prawn density
H	<- params["H"]		#human population abundance
co  <- params["co"]		# population coverage

#  snail dynamics
f 	<- params["f"]  	# instantaneous intrinsic fertility rate of snails
phi <- params["phi"] 	# density-dependent parameter for fertility (assuming cannibalism)
mu	<- params["mu"]		# snail mortality
ze 	<- params["ze"]		# fraction of exposed snails that reproduce

# epidemiological parameters for infected snails
betha <- params["betha"]	# per capita infection rate
k	<- params["k"]		# clumping parameter of teh negative binomila distribution
msr	<- params["msr"]	# miracidia shading rate per reproductive female
sig	<- params["sig"]	# rate from exposed to infected
alfa<- params["alfa"]	# mortality rate of infected snail (if negative means that castrated snail increas their life expectancy - tyet it could be zero)

## predation rates
q	<- params["q"]		# predation rate on susceptible smails

Th  <- params["Th"]		# handling time parameter from functional response


## parameters for the adult worms dynamics in the human hosts
lambda<-params["lambda"]	# cercarie shadding rate divided by cercarie mortality
mua	<- params["mua"]	# adult warm natural mortality =1/4 
dh	<- params["dh"]		# human mortality/fertility rate =1/60 


eta <- 0.5  # *(1 - 1/(1+M/k)^(k+1)) 

rF <- eta*(M*co + U*(1-co))*H # need to wait for the size of the population treated and untreated

# total snail population
N = S+E+I

## now code the model equations
dSdt <- f*(1 - phi*N)*(S+ze*E)           - mu*S - betha*msr*rF*S - (q*S*P/(1+q*N*Th))
dEdt <- betha*msr*rF*S           - (mu + sig)*E - (q*P*E/(1+q*N*Th))
dIdt <- sig*E                  - (mu + alfa)*I - (q*P*I/(1+q*N*Th))
dMdt <- lambda*I                - (mua + dh)*M - log(1/0.1)*M*mda(t)
dUdt <- lambda*I                - (mua + dh)*U 

## combine results into a single vector
dxdt <- c(dMdt, dSdt, dEdt, dIdt, dUdt)
## return result as a list!
list(dxdt)
}




#####################################################################################################
#Simulation loop to explore sensitivity to prawn predation rate using ODE:
#####################################################################################################

## load libraries
library(ggplot2)
library(grid)
library(lattice)

## first create a data frame of the variable parameters:
P <-seq(from=0,to=200,by=10) #varying prawn density from 0 to 200/enclosure
q<-seq(from=0,to=0.03,by=.005)# varying the attack rate on susceptible snails from 0.03 to .005
Thinv<-seq(from=2, to=10, by=2) 
Th<-1/Thinv # varying snails killed per day from 1 to 10
msr<-seq(from=0.1,to=2,by=.1)#varying miracidial shedding rate
siginv<-seq(from=30,to=100,by=10)
sig<-1/siginv# varying latency period from 30 to 100 days
muinv<-seq(from=25,to=100,by=25)
mu=1/muinv# varying mean snail lifespan from 30 to 360 days
f=seq(from=0.05, to=0.6, by=0.05) # sensitivity on fecundity rate in snails
params00<- expand.grid(P=P, Th=Th)

sims = dim(params0)[1] # gives the length of each column of params0
sims = dim(params00)[1] # gives the length of each column of params0

## set c, the fraction treated
co=0.8

## next create a data frame of constant parameters
constantparams<- cbind(
 H = rep(1000,sims),  #human population abundance

#  snail dynamics
 f= rep(.16,sims), # instantaneous intrinsic fertility rate of snails (day^-1)
 phi= rep((1-1/(80*0.16))/10000,sims),	# density-dependent parameter for fertility 
 mu	= rep(1/80,sims), 	# snail mortality (mean life expectancy: 2 months)
 ze = rep(.5,sims), 		# fraction of exposed snails that reproduce

# epidemiological parameters for infected snails
 betha 	= rep((4*10^-6),sims),	# per capita infection rate
 k		= rep(0.2,sims),		# clumping parameter of teh negative binomila distribution
 msr	= rep(.5, sims),		# miracidia shedding rate per reproductive female (devided by miracidia mortality rate)
 sig	= rep(1/50,sims),		# rate from exposed to infected
 alfa 	= rep(1/20,sims),		# mortality rate of infected snail (if negative means that castrated snail increase their life expectancy - tyet it could be zero)


## parameters for the adult worms dynamics in the human hosts
 lambda	= rep(.00005,sims),		# cercaria shadding rate divided by cercaria mortality
 mua  	= rep(1/(3*365),sims),# adult worm natural mortality =1/3 years
 dh		= rep(1/(60*365),sims),# mortality/fertility rate of human population

## predation rates
 q	= rep(0.0003,sims)  #,	# predation rate on snails,
 ##Th = rep(0.1,sims)  # ave 5 snails per day
 
) # end constant parameters

## combine variable and constant parameters into one:
params1 <- cbind(params00,constantparams)
head(params1)

simprawn2 <- function (params1) 
### simprawn() explores the sensitivity of 
###   const.prawn.pop.model.2pops() to parameter combinations 
###   related to prawn density and consumption of snails.
###   It takes a dataframe of parameter combinations and 
###   outputs a dataframe in which each row represents a 
###   different parameter combination together with its 
###   resulting equilibrium simulation output from the program 
###   const.prawn.pop.model().
###
### usage: first, create the dataframe params0 by hand or 
###   using code such as
###   params0 = expand.grid(m=c(0.1,0.2,0.3), nu=c(0.1,0.2))
###   simprawnout = simprawn(params0)
###   Now run stats on, plot etc, simprawnout
###
### version 0.1  15 Sept 2013
### Sanna Sokolow sokolow@stanford.edu and Giulio DeLeo 
### deleo@stanford.edu
###
{ 	
	## first we create a new dataframe whose first columns 
	##   are parameter values and whose last ones are 
	##   simulation results
	sims = dim(params1)[1] # gives number of rows in params1
	m = numeric(sims) # creates a vector of zero's of length sims
	s = numeric(sims) # creates a vector of zero's of length sims
	e = numeric(sims) # creates a vector of zero's of length sims
	I = numeric(sims) # creates a vector of zero's of length sims
	

	## each parameter combination is run through 
	## const.prawn.pop.model:
	for(i in 1:sims)
	{
		# simout contains the entire output corresponding 
		# to one row of parameters. We then pull out 
		# interesting data from it.
		
		paramsi<-unlist(as.vector(params1[i,]))
		# now integrate the system
		out <- as.data.frame(
 		ode(
   		func		= const.prawn.pop.model,
   			y		= xstart,
   			times	= times,
   			parms	= paramsi
  )
)
		m[i]=as.numeric(tail(out$M, n=1))
		s[i]=as.numeric(tail(out$S, n=1))
		e[i]=as.numeric(tail(out$E, n=1))
		I[i]=as.numeric(tail(out$I, n=1))
		
	}
	simprawnout= cbind(params1,m,s,e,I)
	return (simprawnout)
}

y=simprawn2(params=params1)
y$z=y$e+y$s+y$I
y$pi=y$I/y$z
y$prev=Prevalence(y$m)

wireframe(prev ~ P*(1/Th), data = y,
  xlab = "P", ylab = "consumption rate at saturation", zlab= list(label="proportion infected people", rot=90),ylim=c(0,max(1/Th)),
  main = "Sensitivity to prawn stocking and prawn mean size",
  drape = TRUE,
  colorkey = TRUE,
  screen = list(z = -120, x = -60),scales=list(arrows=FALSE)
)

y$Ro = (((y$betha*y$msr*(y$f-y$mu))/(y$f*y$phi))/(y$mu+y$sig))*(y$sig/(y$mu+y$alfa))*((y$lambda*0.5*1000)/(y$mua+y$dh))

par(mfrow=c(1,2))
yp=subset(y, Th==0.1)
plot(yp$P,yp$z, type='l', ylim=c(0,max(yp$z)), xlab="number of prawns per enclosure", ylab="number of snails")
lines(yp$P,yp$I, col="red")

yp2=subset(y, Th==1/6)
lines(yp2$P,yp2$z, type='l', lty="dotted")

lines(yp2$P,yp2$I, col="red", lty="dotted")

yp3=subset(y, Th==1/4)
lines(yp3$P,yp3$z, type='l', lty="dotted")

lines(yp3$P,yp3$I, col="red", lty="dotted")

plot(yp$P, yp$prev*100, type='l', ylim=c(0,100), xlab="number of prawns per enclosure", ylab="human burden/prevalence")
lines(yp$P,yp$m, col="red")
lines(yp2$P,yp2$prev*100, type='l', lty="dotted")

lines(yp2$P,yp2$m, col="red", lty="dotted")
lines(yp3$P,yp3$prev*100, type='l', lty="dotted")

lines(yp3$P,yp3$m, col="red", lty="dotted")

y0=subset(y, P==0)
y25=subset(y,P==25)
y50=subset(y,P==50)
y75=subset(y, P==75)
y100=subset(y,P==100)
y200=subset(y,P==200)

par(mfrow=c(1,2))
plot(y25$Ro,100-100*(y25$m/y0$m), type ='l', xlab="prawn free Ro", ylab="percent reduction in prevalence or burden", ylim=c(0,100))
lines(y50$Ro,100-100*(y50$m/y0$m))
lines(y75$Ro,100-100*(y75$m/y0$m))
lines(y25$Ro,100-100*(y25$prev/y0$prev), lty='dashed')
lines(y50$Ro,100-100*(y50$prev/y0$prev), lty='dashed')
lines(y75$Ro,100-100*(y75$prev/y0$prev), lty='dashed')

plot(y25$Ro,100-(100*(y25$I/(y0$I))), type ='l', lty="dotted", xlab="prawn free Ro", ylab="percent reduction in snails", ylim=c(0,100))
lines(y50$Ro,100-(100*(y50$I/(y0$I))), type ='l', lty="dotted")
lines(y75$Ro,100-(100*(y75$I/(y0$I))), type ='l', lty="dotted")
lines(y25$Ro,100-(100*(y25$z/(y0$z))), type ='l')
lines(y50$Ro,100-(100*(y50$z/(y0$z))), type ='l')
lines(y75$Ro,100-(100*(y75$z/(y0$z))), type ='l')



lines(y25$Ro,100-(100*(y25$e/(y0$e))), type ='l')
lines(y50$Ro,100-(100*(y50$e/(y0$e))), type ='l')
lines(y75$Ro,100-(100*(y75$e/(y0$e))), type ='l')

plot(y25$Ro,y25$m/y0$m, type ='l', xlab="prawn free Ro", ylab="relative prevalence or burden", ylim=c(0,1))
lines(y50$Ro,y50$m/y0$m)
lines(y100$Ro,y100$m/y0$m)
lines(y25$Ro,y25$prev/y0$prev, lty='dashed')
lines(y50$Ro,y50$prev/y0$prev, lty='dashed')
lines(y100$Ro,y100$prev/y0$prev, lty='dashed')

plot(y50$Ro, y50$I/y0$I, type='l', xlab="prawn free Ro", ylab="infected and total snail relative abundance", ylim=c(0,1)) 
lines(y100$Ro, y100$I/y0$I)
lines(y25$Ro, y25$I/y0$I)
lines(y50$Ro,y50$z/y0$z, lty="dashed")
lines(y100$Ro, y100$z/y0$z, lty="dashed")
lines(y25$Ro, y25$z/y0$z, lty="dashed")

plot(y50$Ro,y50$m/y0$m, type ='l', xlab="prawn free Ro", ylab="relative prevalence or burden", ylim=c(0,1))
lines(y100$Ro,y100$m/y0$m)
lines(y200$Ro,y200$m/y0$m)
lines(y50$Ro,y50$prev/y0$prev, lty='dashed')
lines(y100$Ro,y100$prev/y0$prev, lty='dashed')
lines(y200$Ro,y200$prev/y0$prev, lty='dashed')

plot(y50$Ro, y50$I/y0$I, type='l', xlab="prawn free Ro", ylab="infected and total snail relative abundance", ylim=c(0,1)) 
lines(y100$Ro, y100$I/y0$I)
lines(y200$Ro, y200$I/y0$I)
lines(y50$Ro,y50$z/y0$z, lty="dashed")
lines(y100$Ro, y100$z/y0$z, lty="dashed")
lines(y200$Ro, y200$z/y0$z, lty="dashed")

y$i=y$I
y2=subset(y, 1/Th==2)
y4=subset(y,1/Th==4)
y6=subset(y,1/Th==6)
y8=subset(y, 1/Th==8)
y10=subset(y, 1/Th==10)

plot1=qplot(P,m, data=y10, geom='line')+theme_classic()
plot2=qplot(P,i, geom='line', data=y10)+theme_classic()
plot3=qplot(P,prev, geom='line', data=y10)+theme_classic()
plot4=qplot(P,z, geom='line', data = y10)+theme_classic()
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(1, 2))
print(plot3, vp = vplayout(2, 1))
print(plot4, vp = vplayout(2, 2))

plot1=qplot(P,m, data=y6, geom='line')+theme_classic()
plot2=qplot(P,i, geom='line', data=y6)+theme_classic()
plot3=qplot(P,prev, geom='line', data=y6)+theme_classic()
plot4=qplot(P,z, geom='line', data = y6)+theme_classic()
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(1, 2))
print(plot3, vp = vplayout(2, 1))
print(plot4, vp = vplayout(2, 2))

yp=subset(y, msr==0.8)
plot(yp$P,yp$z)

##3d graphics
library(lattice)
wireframe(m ~ P*Th, data = y,
  xlab = "X Coordinate (feet)", ylab = "Y Coordinate (feet)",
  main = "Surface elevation data",
  drape = TRUE,
  colorkey = TRUE,
  screen = list(z = -60, x = -60)
)

attach(y)
library(rgl)

plot3d(Th, P, m, col="red", size=3)

plot1=qplot(y$P,y$m, geom='line')
plot2=qplot(y$P,y$I, geom='line')
plot3=qplot(y$P,y$e, geom='line')
plot4=qplot(y$P,y$s, geom='line')

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(1, 2))
print(plot3, vp = vplayout(2, 1))
print(plot4, vp = vplayout(2, 2))





# end simprawn sensitivity analysis


##################################################################################
const.prawn.pop.model.discrete.2pops <- function (t, x, params, co) {

### const.prawn.pop.model simulates the status of a coral metapopulation in which 
###   M is the mean parasite burden in humans, 
###   S gives the density of susceptible snails,
###   E gives the density of exposed but not yet infective 
###      snails, 
###   I gives the density of infected snails. 
### 
###    version 0.1 15 Sept 2013
###     SSokolow sokolow@stanford.edu G DeLeo deleo@stanford.edu
### 

   times=times
   tmax = length(times)
   M = numeric(tmax+1) #the history of M
   S = numeric(tmax+1) # history of S
   E = numeric(tmax+1) # history of E 
   I = numeric(tmax+1) # history of I 
   U=numeric(tmax+1)  # history of U
   M <- x[1]	# mean parasite burden in humans 
S <- x[2]	# density of susceptible snails
E <- x[3]	# density of infected but not yet infective snails
I <- x[4]	# density of infected snails
U<-x[1]  # mean parasite burden in proportion of the untreated population
      ## main simulation loop 
   for(t in 1:tmax) 
   { 
       ## first extract the parameters
P 	<- params["P"]	# prawn density
H	<- params["H"]		#human population abundance
co<-co # population coverage as an input to the simulation

#  snail dynamics
f 	<- params["f"]  	# instantaneous intrinsic fertility rate of snails
phi <- params["phi"] 	# density-dependent parameter for fertility (assuming cannibalism)
mu	<- params["mu"]		# snail mortality
ze 	<- params["ze"]		# fraction of exposed snails that reproduce

# epidemiological parameters for infected snails
betha <- params["betha"]	# per capita infection rate
k	<- params["k"]		# clumping parameter of teh negative binomila distribution
msr	<- params["msr"]	# miracidia shading rate per reproductive female
sig	<- params["sig"]	# rate from exposed to infected
alfa<- params["alfa"]	# mortality rate of infected snail (if negative means that castrated snail increas their life expectancy - tyet it could be zero)

## predation rates
q	<- params["q"]		# predation rate on snails
Th <- params["Th"] # handling time
## parameters for the adult worms dynamics in the human hosts
lambda <-params["lambda"]	# cercarie shadding rate divided by cercarie mortality
mua	<- params["mua"]	# adult warm natural mortality =1/4 
dh	<- params["dh"]		# human mortality/fertility rate =1/60 
## now code the model equations

S[t+1] <- as.numeric(S[t] + f*(1 - phi*S[t]- phi*E[t] - phi*I[t])*(S[t]+ze*E[t]) -          mu*S[t] - (q*S[t]*P/(1+q*Th*(S[t]+E[t]+I[t]))) - betha*msr*0.5*(M[t]*co + U[t]*(1-co))*S[t]*H)
E[t+1] <- as.numeric(E[t] + betha*msr*0.5*(M[t]*co + U[t]*(1-co))*S[t]*H         -  (mu + sig)*E[t] - (q*P*E[t])/(1+q*Th*(S[t]+E[t]+I[t])))
I[t+1] <- as.numeric(I[t] + sig*E[t]                  - (mu + alfa)*I[t] - (q*P*I[t]/(1+q*Th*(S[t]+E[t]+I[t]))))
M[t+1] <- as.numeric(M[t] + lambda*I[t]               -  (mua + dh)*M[t] - 0.99*mda(t)*M[t])
U[t+1] <- as.numeric(U[t] + lambda*I[t]  -  (mua + dh)*U[t])

} 
   return(data.frame(cbind(t=times, M=M[-1], S=S[-1], E=E[-1], I=I[-1], U=U[-1])))
} 
####################################################################







x=const.prawn.pop.model.discrete(t=times, x=xstart, params=params)
par(mfrow=c(3,1))
time_year = x$t/365
plot(    S~time_year, data=x,type='l',col='green', lwd = 2, ylim=c(0,max(c(S,E,I))), ylab = "Susceptible snails")
points(E~time_year, data=x,type='l', ylab= "Exposed snails", lwd = 2, col = "blue")
points(I~time_year, data=x,type='l', ylab= "Infected snails", lwd = 2, col = "red")
plot(M~time_year, data=x,type='l',col='green', lwd = 2,  ylab = "adult worm per host")
plot(I~M, data=x,type='l',col='green', lwd = 2,  ylab = "phase plan")

## Comparing with and without prawns

## MDA function q 1 yrs
co<-0.10

mda = function(t) {
   ifelse((t == 365*4.4)|(t == 365*5+540)|(t == 365*5+690) | (t == 365*5+900), 1, 0)}
## ifelse(t %% 365 == 0, 1, 0)}
## load params with 0 prawns
params <- c(
 P = 0,			# prawn density
 H = 1000, 		#human population abundance

#  snail dynamics
 f 	= 0.16, 	# instantaneous intrinsic fertility rate of snails times probability of survival to juvenile(day^-1)
 phi= (1-1/(80*0.16))/10000,	# density-dependent parameter for fertility (assuming cannibalism)
 mu	= 1/80, 	# snail mortality (mean life expectancy: 12 weeks)
 ze = .5, 		# fraction of exposed snails that reproduce

# epidemiological parameters for infected snails
 betha 	= (4*10^-6),	# per capita infection rate
 k		= 0.25,		# clumping parameter of teh negative binomila distribution
 msr	= 0.5,		# miracidia shedding rate per reproductive female (devided by miracidia mortality rate)
 sig	= 1/50,		# rate from exposed to infected
 alfa 	= 1/20,		# mortality rate of infected snail (if negative means that castrated snail increase their life expectancy - tyet it could be zero)

## predation rates
 q	= 0.003,		# predation rate on susceptible snails
 Th = 0.2,
## parameters for human host
 lambda	= 0.00004,		# cercaria shadding rate divided by cercaria mortality
 mua  	= 1/(3*365), # adult worm natural mortality =1/3 years
 dh		= 1/(60*365)# mortality/fertility rate of human population

) # end parameters

x=const.prawn.pop.model.discrete(t=times, x=xstart, params=params)
x4=const.prawn.pop.model.discrete.2pops(t=times, x=xstart, params=params, co=co)
Ro = (((params["betha"]*params["msr"]*(params["f"]-params["mu"]))/(params["f"]*params["phi"]))/(params["mu"]+params["sig"]))*(params["sig"]/(params["mu"]+params["alfa"]))*((params["lambda"]*0.5*1000)/(params["mua"]+params["dh"]))
Ro

## reload with prawns - no mda

params["P"]<-25
x2=const.prawn.pop.model.discrete.nomda(t=times, x=xstart, params=params)
time_year = x$t/365



## reload params with prawns + mda



x3=const.prawn.pop.model.discrete(t=times, x=xstart, params=params)
x5=const.prawn.pop.model.discrete.2pops(t=times, x=xstart, params=params, co=co)


Prevalence <- function(x) {
	k=as.numeric(params["k"])
	p=1 - (1/(1+x/k)^(k))*(1+x/(1+x/k)) # fraction of humans with at least 2 parasites
	return(p)
}

Prevalence_heavy <- function(x) {
	k=as.numeric(params["k"])
	p= (1+x)^-21*(x)^20 # fraction of humans with at least 2 parasites
	return(p)
}

x$prevheavy=Prevalence_heavy(x$M)
plot(x$prevheavy)

x$prev=Prevalence(x$M)
x2$prev=Prevalence(x2$M)
x3$prev=Prevalence(x3$M)
x4$prev=Prevalence(x4$M)
x4$H= (x4$M*co) + (x4$U*(1-co))
x4$prev2=Prevalence(x4$H)
x4$prev3=Prevalence(x4$U)
x5$prev=Prevalence(x5$M)
x5$H= (x5$M*co) + (x5$U*(1-co))
x5$prev2=Prevalence(x5$H)
x5$prev3=Prevalence(x5$U)

par(mfrow=c(2,3))
time_year = x$t/365
plot(M~time_year, data=x4, xlim=c(0,20),type='l', lwd = 2,  cex.lab=1.5, cex.axis=1.5, ylab = "Mean adult worms per human host")
lines(H~time_year, data=x4, lty="dotted")
lines(U~time_year, data=x4, lty="dashed")
plot(M~time_year, data=x2,xlim=c(0,20),type='l', cex.lab=1.5, cex.axis=1.5, lwd = 2,  ylab = "Mean adult worms per human host")
plot(M~time_year, data=x5,xlim=c(0,20),type='l', cex.lab=1.5, cex.axis=1.5, lwd = 2,  ylab = "Mean adult worms per human host")
lines(H~time_year, data=x5, lty="dotted")
lines(U~time_year, data=x5, lty="dashed")

plot(prev*100~time_year, data=x4, ylim=c(0,100), xlim=c(0,20),type='l',cex.lab=1.5, cex.axis=1.5,  lwd = 2,  ylab = "Prevalence")
lines(prev2*100~time_year, data=x4, lty="dotted")
lines(prev3*100~time_year, data=x4, lty="dashed")

plot(prev*100~time_year, data=x2, ylim=c(0,100), xlim=c(0,20),type='l', cex.lab=1.5, cex.axis=1.5, lwd = 2,  ylab = "Prevalence")
plot(prev*100~time_year, data=x5, ylim=c(0,100), xlim=c(0,20),type='l', cex.lab=1.5, cex.axis=1.5, lwd = 2,  ylab = "Prevalence")
lines(prev2*100~time_year, data=x5, lty="dotted")
lines(prev3*100~time_year, data=x5, lty="dashed")

## effect of prazi @ co% coverage for comparison with data (set co to 0.1, set treatment frequency to q 180 days)
par(mfrow=c(1,3))
time_year = x$t/365
plot(prev*100~time_year, data=x4, ylim=c(0,100), xlim=c(0,3),type='l',cex.lab=1.5, cex.axis=1.5,  lwd = 2,  ylab = "Prevalence")
lines(prev*100~time_year, data=x5, lty="dotted")
lines(((S+E+I)/1000)+60)~time_year, data=x5, lty="dotted")

plot(M~time_year, data=x4, xlim=c(0,3),type='l', lwd = 2,  cex.lab=1.5, cex.axis=1.5, ylab = "Mean adult worms per human host")
lines(M~time_year, data=x5, lty="dotted")

plot(I~time_year, ylim=c(0,max(I)), data=x4, xlim=c(0,3),type='l', lwd = 2,  cex.lab=1.5, cex.axis=1.5, ylab = "snails")
 lines(I~time_year, data=x5, lty="dotted")


1-x5$M[890]/x4$M[890]
1-x5$prev[890]/x4$prev[890]

## With seasonality######################################################################################################
const.prawn.pop.model.discrete <- function (t, x, params) {

### const.prawn.pop.model simulates the status of a coral metapopulation in which 
###   M is the mean parasite burden in humans, 
###   S gives the density of susceptible snails,
###   E gives the density of exposed but not yet infective 
###      snails, 
###   I gives the density of infected snails. 
### 
###    version 0.1 15 Sept 2013
###     SSokolow sokolow@stanford.edu G DeLeo deleo@stanford.edu
### 

   times=times
   tmax = length(times)
   M = numeric(tmax+1) #the history of M
   S = numeric(tmax+1) # history of S
   E = numeric(tmax+1) # history of E 
   I = numeric(tmax+1) # history of I 
   M <- x[1]	# mean parasite burden in humans 
S <- x[2]	# density of susceptible snails
E <- x[3]	# density of infected but not yet infective snails
I <- x[4]	# density of infected snails

      ## main simulation loop 
   for(t in 1:tmax) 
   { 
       ## first extract the parameters
P 	<- params["P"]	# prawn density
H	<- params["H"]		#human population abundance

#  snail dynamics
f 	<- params["f"]  	# instantaneous intrinsic fertility rate of snails
phi <- params["phi"] 	# density-dependent parameter for fertility (assuming cannibalism)
mu	<- params["mu"]		# snail mortality
ze 	<- params["ze"]		# fraction of exposed snails that reproduce

# epidemiological parameters for infected snails
betha <- params["betha"]	# per capita infection rate
k	<- params["k"]		# clumping parameter of teh negative binomila distribution
msr	<- params["msr"]	# miracidia shading rate per reproductive female
sig	<- params["sig"]	# rate from exposed to infected
alfa<- params["alfa"]	# mortality rate of infected snail (if negative means that castrated snail increas their life expectancy - tyet it could be zero)

## predation rates
q	<- params["q"]		# predation rate on snails
Th <- params["Th"] # handling time
## parameters for the adult worms dynamics in the human hosts
lambda <-params["lambda"]	# cercarie shadding rate divided by cercarie mortality
mua	<- params["mua"]	# adult warm natural mortality =1/4 
dh	<- params["dh"]		# human mortality/fertility rate =1/60 
## now code the model equations

S[t+1] <- as.numeric(S[t] + f*((cos(2*pi*t/365) + 2)/2)*(1 - phi*S[t]- phi*E[t] - phi*I[t])*(S[t]+ze*E[t]) -          mu*S[t] - (q*S[t]*P/(1+q*Th*(S[t]+E[t]+I[t]))) - betha*msr*0.5*M[t]*S[t]*H)
E[t+1] <- as.numeric(E[t] + betha*msr*0.5*M[t]*S[t]*H         -  (mu + sig)*E[t] - (q*P*E[t])/(1+q*Th*(S[t]+E[t]+I[t])))
I[t+1] <- as.numeric(I[t] + sig*E[t]                  - (mu + alfa)*I[t] - (q*P*I[t]/(1+q*Th*(S[t]+E[t]+I[t]))))
M[t+1] <- as.numeric(M[t] + lambda*I[t]               -  (mua + dh)*M[t] - 0.99*mda(t)*M[t])

} 
   return(data.frame(cbind(t=times, M=M[-1], S=S[-1], E=E[-1], I=I[-1])))
} 

## Now with no mda:
###########################################################################
const.prawn.pop.model.discrete.nomda <- function (t, x, params) {

### const.prawn.pop.model simulates the status of a coral metapopulation in which 
###   M is the mean parasite burden in humans, 
###   S gives the density of susceptible snails,
###   E gives the density of exposed but not yet infective 
###      snails, 
###   I gives the density of infected snails. 
### 
###    version 0.1 15 Sept 2013
###     SSokolow sokolow@stanford.edu G DeLeo deleo@stanford.edu
### 

   times=times
   tmax = length(times)
   M = numeric(tmax+1) #the history of M
   S = numeric(tmax+1) # history of S
   E = numeric(tmax+1) # history of E 
   I = numeric(tmax+1) # history of I 
   M <- x[1]	# mean parasite burden in humans 
S <- x[2]	# density of susceptible snails
E <- x[3]	# density of infected but not yet infective snails
I <- x[4]	# density of infected snails

      ## main simulation loop 
   for(t in 1:tmax) 
   { 
       ## first extract the parameters
P 	<- params["P"]	# prawn density
H	<- params["H"]		#human population abundance

#  snail dynamics
f 	<- params["f"]  	# instantaneous intrinsic fertility rate of snails
phi <- params["phi"] 	# density-dependent parameter for fertility (assuming cannibalism)
mu	<- params["mu"]		# snail mortality
ze 	<- params["ze"]		# fraction of exposed snails that reproduce

# epidemiological parameters for infected snails
betha <- params["betha"]	# per capita infection rate
k	<- params["k"]		# clumping parameter of teh negative binomila distribution
msr	<- params["msr"]	# miracidia shading rate per reproductive female
sig	<- params["sig"]	# rate from exposed to infected
alfa<- params["alfa"]	# mortality rate of infected snail (if negative means that castrated snail increas their life expectancy - tyet it could be zero)

## predation rates
q	<- params["q"]		# predation rate on snails
Th <- params["Th"] # handling time
## parameters for the adult worms dynamics in the human hosts
lambda <-params["lambda"]	# cercarie shadding rate divided by cercarie mortality
mua	<- params["mua"]	# adult warm natural mortality =1/4 
dh	<- params["dh"]		# human mortality/fertility rate =1/60 
## now code the model equations

S[t+1] <- as.numeric(S[t] + f*((cos(2*pi*t/365) + 2)/2)*(1 - phi*S[t]- phi*E[t] - phi*I[t])*(S[t]+ze*E[t]) -          mu*S[t] - (q*S[t]*P/(1+q*Th*(S[t]+E[t]+I[t]))) - betha*msr*0.5*M[t]*S[t]*H)
E[t+1] <- as.numeric(E[t] + betha*msr*0.5*M[t]*S[t]*H         -  (mu + sig)*E[t] - (q*P*E[t])/(1+q*Th*(S[t]+E[t]+I[t])))
I[t+1] <- as.numeric(I[t] + sig*E[t]                  - (mu + alfa)*I[t] - (q*P*I[t]/(1+q*Th*(S[t]+E[t]+I[t]))))
M[t+1] <- as.numeric(M[t] + lambda*I[t]               -  (mua + dh)*M[t])

} 
   return(data.frame(cbind(t=times, M=M[-1], S=S[-1], E=E[-1], I=I[-1])))
} 
##################################################################################
const.prawn.pop.model.discrete.2pops <- function (t, x, params, co) {

### const.prawn.pop.model simulates the status of a coral metapopulation in which 
###   M is the mean parasite burden in humans, 
###   S gives the density of susceptible snails,
###   E gives the density of exposed but not yet infective 
###      snails, 
###   I gives the density of infected snails. 
### 
###    version 0.1 15 Sept 2013
###     SSokolow sokolow@stanford.edu G DeLeo deleo@stanford.edu
### 

   times=times
   tmax = length(times)
   M = numeric(tmax+1) #the history of M
   S = numeric(tmax+1) # history of S
   E = numeric(tmax+1) # history of E 
   I = numeric(tmax+1) # history of I 
   U=numeric(tmax+1)  # history of U
   M <- x[1]	# mean parasite burden in humans 
S <- x[2]	# density of susceptible snails
E <- x[3]	# density of infected but not yet infective snails
I <- x[4]	# density of infected snails
U<-x[1]  # mean parasite burden in proportion of the untreated population
      ## main simulation loop 
   for(t in 1:tmax) 
   { 
       ## first extract the parameters
P 	<- params["P"]	# prawn density
H	<- params["H"]		#human population abundance
co<-co # population coverage as an input to the simulation

#  snail dynamics
f 	<- params["f"]  	# instantaneous intrinsic fertility rate of snails
phi <- params["phi"] 	# density-dependent parameter for fertility (assuming cannibalism)
mu	<- params["mu"]		# snail mortality
ze 	<- params["ze"]		# fraction of exposed snails that reproduce

# epidemiological parameters for infected snails
betha <- params["betha"]	# per capita infection rate
k	<- params["k"]		# clumping parameter of teh negative binomila distribution
msr	<- params["msr"]	# miracidia shading rate per reproductive female
sig	<- params["sig"]	# rate from exposed to infected
alfa<- params["alfa"]	# mortality rate of infected snail (if negative means that castrated snail increas their life expectancy - tyet it could be zero)

## predation rates
q	<- params["q"]		# predation rate on snails
Th <- params["Th"] # handling time
## parameters for the adult worms dynamics in the human hosts
lambda <-params["lambda"]	# cercarie shadding rate divided by cercarie mortality
mua	<- params["mua"]	# adult warm natural mortality =1/4 
dh	<- params["dh"]		# human mortality/fertility rate =1/60 
## now code the model equations

S[t+1] <- as.numeric(S[t] + f*((cos(2*pi*t/365) + 2)/2)*(1 - phi*S[t]- phi*E[t] - phi*I[t])*(S[t]+ze*E[t]) -          mu*S[t] - (q*S[t]*P/(1+q*Th*(S[t]+E[t]+I[t]))) - betha*msr*0.5*(M[t]*co + U[t]*(1-co))*S[t]*H)
E[t+1] <- as.numeric(E[t] + betha*msr*0.5*(M[t]*co + U[t]*(1-co))*S[t]*H         -  (mu + sig)*E[t] - (q*P*E[t])/(1+q*Th*(S[t]+E[t]+I[t])))
I[t+1] <- as.numeric(I[t] + sig*E[t]                  - (mu + alfa)*I[t] - (q*P*I[t]/(1+q*Th*(S[t]+E[t]+I[t]))))
M[t+1] <- as.numeric(M[t] + lambda*I[t]               -  (mua + dh)*M[t] - 0.99*mda(t)*M[t])
U[t+1] <- as.numeric(U[t] + lambda*I[t]  -  (mua + dh)*U[t])

} 
   return(data.frame(cbind(t=times, M=M[-1], S=S[-1], E=E[-1], I=I[-1], U=U[-1])))
} 
####################################################################
x=const.prawn.pop.model.discrete(t=times, x=xstart, params=params)
par(mfrow=c(3,1))
time_year = x$t/365
plot(    S~time_year, data=x,type='l',col='green', lwd = 2, ylim=c(0,max(c(S,E,I))), ylab = "Susceptible snails")
points(E~time_year, data=x,type='l', ylab= "Exposed snails", lwd = 2, col = "blue")
points(I~time_year, data=x,type='l', ylab= "Infected snails", lwd = 2, col = "red")
plot(M~time_year, data=x,type='l',col='green', lwd = 2,  ylab = "adult worm per host")
plot(I~M, data=x,type='l',col='green', lwd = 2,  ylab = "phase plan")

## Comparing with and without prawns

## MDA function q 1 yrs
co<-0.10

mda = function(t) {
   ifelse((t == (365*4)+150-630)|(t == (365*4)+150)|(t == (365*4)+170)|(t == (365*4)+300)|(t == (365*4)+510) | (t == (365*4)+690), 1, 0)}## ifelse(t %% 365 == 0, 1, 0)}
## load params with 0 prawns
params <- c(
 P = 0,			# prawn density
 H = 1000, 		#human population abundance

#  snail dynamics
 f 	= 0.16, 	# instantaneous intrinsic fertility rate of snails times probability of survival to juvenile(day^-1)
 phi= (1-1/(80*0.16))/10000,	# density-dependent parameter for fertility (assuming cannibalism)
 mu	= 1/80, 	# snail mortality (mean life expectancy: 12 weeks)
 ze = .5, 		# fraction of exposed snails that reproduce

# epidemiological parameters for infected snails
 betha 	= (4*10^-6),	# per capita infection rate
 k		= 0.25,		# clumping parameter of teh negative binomila distribution
 msr	= 0.5,		# miracidia shedding rate per reproductive female (devided by miracidia mortality rate)
 sig	= 1/50,		# rate from exposed to infected
 alfa 	= 1/20,		# mortality rate of infected snail (if negative means that castrated snail increase their life expectancy - tyet it could be zero)

## predation rates
 q	= 0.003,		# predation rate on susceptible snails
 Th = 0.1,
## parameters for human host
 lambda	= 0.00004,		# cercaria shadding rate divided by cercaria mortality
 mua  	= 1/(3*365), # adult worm natural mortality =1/3 years
 dh		= 1/(60*365)# mortality/fertility rate of human population

) # end parameters

x=const.prawn.pop.model.discrete(t=times, x=xstart, params=params)
x4=const.prawn.pop.model.discrete.2pops(t=times, x=xstart, params=params, co=co)
Ro = (((params["betha"]*params["msr"]*(params["f"]-params["mu"]))/(params["f"]*params["phi"]))/(params["mu"]+params["sig"]))*(params["sig"]/(params["mu"]+params["alfa"]))*((params["lambda"]*0.5*1000)/(params["mua"]+params["dh"]))
Ro

## reload with prawns - no mda

params["P"]<-25
x2=const.prawn.pop.model.discrete.nomda(t=times, x=xstart, params=params)
time_year = x$t/365
## reload params with prawns + mda
x3=const.prawn.pop.model.discrete(t=times, x=xstart, params=params)
x5=const.prawn.pop.model.discrete.2pops(t=times, x=xstart, params=params, co=co)
Prevalence <- function(x) {
	k=as.numeric(params["k"])
	p=1 - (1/(1+x/k)^(k))*(1+x/(1+x/k)) # fraction of humans with at least 2 parasites
	return(p)
}

Prevalence_heavy <- function(x) {
	k=as.numeric(params["k"])
	p= (1+x)^-21*(x)^20 # fraction of humans with at least 2 parasites
	return(p)
}

x$prevheavy=Prevalence_heavy(x$M)
plot(x$prevheavy)

x$prev=Prevalence(x$M)
x2$prev=Prevalence(x2$M)
x3$prev=Prevalence(x3$M)
x4$prev=Prevalence(x4$M)
x4$H= (x4$M*co) + (x4$U*(1-co))
x4$prev2=Prevalence(x4$H)
x4$prev3=Prevalence(x4$U)
x5$prev=Prevalence(x5$M)
x5$H= (x5$M*co) + (x5$U*(1-co))
x5$prev2=Prevalence(x5$H)
x5$prev3=Prevalence(x5$U)

par(mfrow=c(2,3))
time_year = x$t/365
plot(M~time_year, data=x4, xlim=c(0,20),type='l', lwd = 2,  cex.lab=1.5, cex.axis=1.5, ylab = "Mean adult worms per human host")
lines(H~time_year, data=x4, lty="dotted")
lines(U~time_year, data=x4, lty="dashed")
plot(M~time_year, data=x2,xlim=c(0,20),type='l', cex.lab=1.5, cex.axis=1.5, lwd = 2,  ylab = "Mean adult worms per human host")
plot(M~time_year, data=x5,xlim=c(0,20),type='l', cex.lab=1.5, cex.axis=1.5, lwd = 2,  ylab = "Mean adult worms per human host")
lines(H~time_year, data=x5, lty="dotted")
lines(U~time_year, data=x5, lty="dashed")

plot(prev*100~time_year, data=x4, ylim=c(0,100), xlim=c(0,20),type='l',cex.lab=1.5, cex.axis=1.5,  lwd = 2,  ylab = "Prevalence")
lines(prev2*100~time_year, data=x4, lty="dotted")
lines(prev3*100~time_year, data=x4, lty="dashed")

plot(prev*100~time_year, data=x2, ylim=c(0,100), xlim=c(0,20),type='l', cex.lab=1.5, cex.axis=1.5, lwd = 2,  ylab = "Prevalence")
plot(prev*100~time_year, data=x5, ylim=c(0,100), xlim=c(0,20),type='l', cex.lab=1.5, cex.axis=1.5, lwd = 2,  ylab = "Prevalence")
lines(prev2*100~time_year, data=x5, lty="dotted")
lines(prev3*100~time_year, data=x5, lty="dashed")

## effect of prazi @ co% coverage for comparison with data (set co to 0.1, set treatment frequency to q 180 days)
par(mfrow=c(1,3))
time_year = x$t/365
plot(prev*100~time_year, data=x4, ylim=c(0,100), xlim=c(0,3),type='l',cex.lab=1.5, cex.axis=1.5,  lwd = 2,  ylab = "Prevalence")
lines(prev*100~time_year, data=x5, lty="dotted")

plot(M~time_year, data=x4, xlim=c(0,3),type='l', lwd = 2,  cex.lab=1.5, cex.axis=1.5, ylab = "Mean adult worms per human host")
lines(M~time_year, data=x5, lty="dotted")

plot(S+E+I~time_year, ylim=c(0,max(S+E+I)), data=x5, xlim=c(0,3),type='l', lwd = 2,  cex.lab=1.5, cex.axis=1.5, ylab = "snails")
 lines(I~time_year, data=x5, lty="dotted")

effect of prazi @ co% coverage for comparison with data (set co to 0.1, set treatment frequency to q 180 days)
par(mfrow=c(1,3))
time_year = x$t/365
plot(prev*100~time_year, data=x4, ylim=c(0,100), xlim=c(4,10),type='l',cex.lab=1.5, cex.axis=1.5,  lwd = 2,  ylab = "Prevalence")
lines(prev*100~time_year, data=x5, lty="dotted")
lines(((S+E+I)/1000)+80~time_year, data=x5, lty="dotted")

plot(M~time_year, data=x4, xlim=c(4,10),type='l', lwd = 2,  cex.lab=1.5, cex.axis=1.5, ylab = "Mean adult worms per human host")
lines(M~time_year, data=x5, lty="dotted")

plot(I~time_year, ylim=c(0,max(I)), data=x4, xlim=c(4,10),type='l', lwd = 2,  cex.lab=1.5, cex.axis=1.5, ylab = "snails")
 lines(I~time_year, data=x5, lty="dotted")

1-x5$M[2148]/x4$M[2148]
1-x5$prev[2148]/x4$prev[2148]

1-x5$prev[1967]/x4$prev[1967]
(x4$M[1967]-x5$M[1967])/x4$M[1967]

1-x5$M[1759]/x4$M[1759]
1-x5$prev[1759]/x4$prev[1759]

std <- function(x) sd(x)/sqrt(length(x))



######################################################################
######################################################################
## Monte Carlo sensitivity analysis - by brute force
######################################################################

MCprawn <- function (params) 
### MCprawn() explores the sensitivity of 
###   const.prawn.pop.model() to parameter combinations 
###   using a monte carlo method.
###   It takes a dataframe of parameter combinations and 
###   outputs a dataframe in which each row represents a 
###   different parameter combination together with its 
###   resulting simulation output from the program 
###   const.prawn.pop.model().
###
### usage: first, create the dataframe params0 by hand or 
###   using code such as
###   params0 = expand.grid(m=c(0.1,0.2,0.3), nu=c(0.1,0.2))
###   simprawnout = simprawn(params0)
###   Now run stats on, plot etc, simprawnout
###
### version 0.1  15 April 2014
### Sanna Sokolow sokolow@stanford.edu and Giulio DeLeo 
### deleo@stanford.edu
###
{ 	
	## first we create a new dataframe whose first columns 
	##   are parameter values and whose last ones are 
	##   simulation results
	sims = length(params[,1]) # gives number of rows in params1
	m=matrix(nrow=length(times),ncol=sims)
	s=matrix(nrow=length(times),ncol=sims)
	e=matrix(nrow=length(times),ncol=sims)
	I=matrix(nrow=length(times),ncol=sims)
	
	## each parameter combination is run through 
	## const.prawn.pop.model:
	for(i in 1:sims)
	{
		# simout contains the entire output corresponding 
		# to one row of parameters. We then pull out 
		# interesting data from it.
		
		paramsi<-params[i,]
		# now integrate the system
		out <- as.matrix(
 		ode(
   		func		= const.prawn.pop.model,
   			y		= xstart,
   			times	= times,
   			parms	= paramsi
  ))

		m[,i]=out[,"M"]
		s[,i]=out[,"S"]
		e[,i]=out[,"E"]
		I[,i]=out[,"I"]		
	}
	z=seq(1:sims)
	colnames(m)=colnames(m)=c(sapply(z, function(z) paste("M",z, sep="")))
	colnames(s)=colnames(s)=c(sapply(z, function(z) paste("S",z, sep="")))
	colnames(e)=colnames(e)=c(sapply(z, function(z) paste("E",z, sep="")))
	colnames(I)=colnames(I)=c(sapply(z, function(z) paste("I",z, sep="")))	
	m=cbind(times,m)
	s=cbind(times,s)
	e=cbind(times,e)
	I=cbind(times,I)
	
	simprawnout= list(params,m,s,e,I)
	return (simprawnout)
}
############################# using sapply - doesn't work!


MCprawn.discrete <- function (params) 
### MCprawn() explores the sensitivity of 
###   const.prawn.pop.model() to parameter combinations 
###   using a monte carlo method.
###   It takes a dataframe of parameter combinations and 
###   outputs a dataframe in which each row represents a 
###   different parameter combination together with its 
###   resulting simulation output from the program 
###   const.prawn.pop.model().
###
### usage: first, create the dataframe params0 by hand or 
###   using code such as
###   params0 = expand.grid(m=c(0.1,0.2,0.3), nu=c(0.1,0.2))
###   simprawnout = simprawn(params0)
###   Now run stats on, plot etc, simprawnout
###
### version 0.1  15 April 2014
### Sanna Sokolow sokolow@stanford.edu and Giulio DeLeo 
### deleo@stanford.edu
###
{ 	
	## first we create a new dataframe whose first columns 
	##   are parameter values and whose last ones are 
	##   simulation results
	sims = length(params[,1]) # gives number of rows in params1
	m=matrix(nrow=length(times),ncol=sims)
	s=matrix(nrow=length(times),ncol=sims)
	e=matrix(nrow=length(times),ncol=sims)
	I=matrix(nrow=length(times),ncol=sims)
	u=matrix(nrow=length(times),ncol=sims)
	
	## each parameter combination is run through 
	## const.prawn.pop.model:
	for(i in 1:sims)
	{
		# simout contains the entire output corresponding 
		# to one row of parameters. We then pull out 
		# interesting data from it.
		
		paramsi<-params[i,]
		out=const.prawn.pop.model.discrete.2pops(t=times, x=xstart, params=paramsi, co=co)
	
		m[,i]=out[,"M"]
		s[,i]=out[,"S"]
		e[,i]=out[,"E"]
		I[,i]=out[,"I"]		
		u[,i]=out[,"U"]
	}
	z=seq(1:sims)
	colnames(m)=colnames(m)=c(sapply(z, function(z) paste("M",z, sep="")))
	colnames(s)=colnames(s)=c(sapply(z, function(z) paste("S",z, sep="")))
	colnames(e)=colnames(e)=c(sapply(z, function(z) paste("E",z, sep="")))
	colnames(I)=colnames(I)=c(sapply(z, function(z) paste("I",z, sep="")))	
	colnames(u)=colnames(u)=c(sapply(z, function(z) paste("U",z, sep="")))	
	m=cbind(times,m)
	s=cbind(times,s)
	e=cbind(times,e)
	I=cbind(times,I)
	u=cbind(times, u)
	
	MCout= list(params,m,s,e,I,u)
	return (MCout)
}

############################# using sapply - doesn't work!
MCprawn <- function (params) 
### MCprawn() explores the sensitivity of 
###   const.prawn.pop.model() to parameter combinations 
###   using a monte carlo method.
###   It takes a dataframe of parameter combinations and 
###   outputs a dataframe in which each row represents a 
###   different parameter combination together with its 
###   resulting simulation output from the program 
###   const.prawn.pop.model().
###
### usage: first, create the dataframe params0 by hand or 
###   using code such as
###   params0 = expand.grid(m=c(0.1,0.2,0.3), nu=c(0.1,0.2))
###   simprawnout = simprawn(params0)
###   Now run stats on, plot etc, simprawnout
###
### version 0.1  15 April 2014
### Sanna Sokolow sokolow@stanford.edu and Giulio DeLeo 
### deleo@stanford.edu
###
{ 	
	## first we create a new dataframe whose first columns 
	##   are parameter values and whose last ones are 
	##   simulation results
	sims = length(params[,1])# gives number of rows in params1
	m=matrix(nrow=length(times),ncol=sims)
	s=matrix(nrow=length(times),ncol=sims)
	e=matrix(nrow=length(times),ncol=sims)
	I=matrix(nrow=length(times),ncol=sims)
	
	## each parameter combination is run through 
	## const.prawn.pop.model:
	sims = length(params[,1])
	v=(1:sims)
	sapply(v, function (v)
	{
		# simout contains the entire output corresponding 
		# to one row of parameters. We then pull out 
		# interesting data from it.
		
		paramsi<-params[v,]
		# now integrate the system
		out <- as.matrix(
 		ode(
   		func		= const.prawn.pop.model,
   			y		= xstart,
   			times	= times,
   			parms	= paramsi
))
		m[,v]<-out[,"M"]
		s[,v]=out[,"S"]
		e[,v]=out[,"E"]
		I[,v]=out[,"I"]
	})
		
	simprawnout= list(params,m,s,e,I)
	return (simprawnout)
}

colnames(m)=seq(1:sims)
	colnames(s)=seq(1:sims)
	colnames(e)=seq(1:sims)
	colnames(I)=seq(1:sims)
	m=cbind(times,m)
	s=cbind(times,s)
	e=cbind(times,e)
	I=cbind(times,I)

#########################################################################################
## run the function and extract relevant data

x=MCprawn(params=params1)

##plot all trajectories of M (2nd element in the list after params)
sims= length(params1[,1])
y=subset(x[[2]], select=-times)
z=seq(1:sims)
plot(times,y[,2], type='l', ylab="M", ylim=c(0,max(y)))
sapply(z, function(z) lines(times, y[,z]))

#using ggplot
library(ggplot2)
y=x[[2]]
q=qplot(data=y, times, M2, geom='line', ylab="M")
a=1:sims
z=sapply(a, function(a) paste("M",a, sep=""))
sapply(z, function(z) q=q+geom_line(aes(times, z)))


## creating subsets at each level of P
y=subset(x[[2]], select=-times)
yp=cbind(x[[1]],t(y))
y3=as.data.frame(yp)
y3$Ro = y3$betha*y3$msr*eta*y3$H*((1-(y3$mu/y3$f))/y3$phi)*(1/(y3$mu+y3$sig))*(y3$sig/(y3$mu+y3$alfa))*(y3$lambda/(y3$mua+y3$dh))

ypsub=subset(y3, P==60)
ypsub=ypsub[,17:10967]
ypsub=t(ypsub)

## order columns with respect to equilibrium value of M


ypsub=ypsub[,order(ypsub[nrow(ypsub),])]


## plot median, 95th and 5th percentiles for M
plot(times, ypsub[,7], type='l', lty=4, ylab="M", ylim=c(0,max(ypsub)))
lines(times, ypsub[,55], lty=1)
lines(times, ypsub[,110], lty=4)

## plotting correlations
par(mfrow=c(1,4))
plot(yp[,1],yp[,10967])  ## P and equil M
plot(yp[,2],yp[,10967]) ## H and equil M
plot(yp[,3],yp[,10967]) ## f and equil M
plot(yp[,4],yp[,10967])## phi and equil M
plot(yp[,5],yp[,10967])## mu and equil M
plot(yp[,6],yp[,10967])## ze and equil M
plot(yp[,7],yp[,10967])## k and equil M
plot(yp[,8],yp[,10967])## msr and equil M
plot(yp[,9],yp[,10967])## sig and equil M
plot(yp[,10],yp[,10967])## alfa and equil M
plot(yp[,11],yp[,10967])## q and equil M
plot(yp[,12],yp[,10967])## Th and equil M
plot(yp[,13],yp[,10967])## lambda and equil M
plot(yp[,14],yp[,10967])## mua and equil M

##plot all trajectories of S (3rd element in the list after params)
y=x[[3]]
z=seq(1:sims)
plot(times,y[,1], type='l')
sapply(z, function(z) lines(times, y[,z]))

## order columns with respect to equilibrium value of S
y=subset(y,select=-times)
y=y[,order(y[nrow(y),])]
y=cbind(times,y)
tail(y)

## plot median, 95th and 5th percentiles for S
plot(times, y[,951], type='l', lty=4, ylab="S", ylim=c(min(c(y[,951],y[,51],y[,501])), max(c(y[,951],y[,51],y[,501]))))
lines(times, y[,501], lty=1)
lines(times, y[,51], lty=4)

## create data frame with equilibrium values of M,S,E,I and parameters and plot parameters versus equil values
a=x[[1]]
y=x[[2]]
y=subset(y,select=-times)
b=cbind(a, M=y[nrow(y),])
z=x[[3]]
z=subset(z,select=-times)
b=cbind(b, S=z[nrow(z),])
zz=x[[4]]
zz=subset(zz,select=-times)
b=cbind(b, E=zz[nrow(zz),])
zzz=x[[5]]
zzz=subset(zzz,select=-times)
b=cbind(b, I=zzz[nrow(zzz),])

b=as.data.frame(b)
b$Ro = b$betha*b$msr*eta*b$H*((1-(b$mu/b$f))/b$phi)*(1/(b$mu+b$sig))*(b$sig/(b$mu+b$alfa))*(b$lambda/(b$mua+y3$dh))

plot(b$P, b[,length(b-3)], type='l', ylab="M")

plot(b$P, b[,length(b-2)], type='l', ylab="S")

plot(b$P, b[,length(b-1)], type='l', ylab="E")

plot(b$P, b[,length(b)], type='l', ylab="I")

## create data frame with time-to-quasi-elimination and parameters and plot parameters versus time to elimination
a=x[[1]]
y=x[[5]] ## infected snails
y=y[10:10951,]
y=as.data.frame(y)
z=seq(1,sims)
q=sapply(z, function(z) min(which(y[,z]<1))) ## first row (first day in simulation) when no. infected snails <1
b=cbind(a, q)
plot(b$P,b[,length(b)]/365, type='l')

## create a matrix for input parameters
## assume a given size prawn
## start with large prawns as in field experiment

## calculate params:
Carrycapacity=rnorm(1000, mean=10000, sd = 1000)


P=sample(seq(0, 300, by=1), 1000, replace=T) ## number prawns regular interval
H=rnorm(1000, mean=1000, sd=100)  ## human population size random f=rnorm(1000, mean=0.6, sd=0.06) ## fecundity random rnorm 0.16 with 10% CV
phi=1/Carrycapacity ## ~inverse carrying capacity random
f=rnorm(1000, mean=0.16, sd= 0.16/10) ## fecundity
mu=rnorm(1000, mean=0.0125, sd = 0.0125/10) ## mu random positive
betha=rnorm(1000, mean=0.000004, sd=0.000004/10) ## transmission rate man-snail random
k=runif(1000, 0.2, 0.7) ## clumping parameter random uniform in empirically derived bounds
msr=rnorm(1000, mean=0.8, sd = 0.8/10) ## msr random 
ze=rnorm(1000, mean=0.5, sd= 0.5/10) ## fraction of exposed that reproduce
sig=rnorm(1000, mean=0.02, sd = 0.02/10) ## sig random norm
alfa=rnorm(1000, mean=0.05, sd = 0.05/10) ## alfa random norm
q=runif(1000, 0.0003, 0.003) ## q random unif from 10X lower to equal to that found in lab
maxconsump=rnorm(1000, mean=7.9, sd=1.02) ## random normal based on empirical lab studies on satiation
Th=1/maxconsump
lambda=rnorm(1000, mean=0.00005, sd = 0.00005/10) ## lambda random norm
mua=rnorm(1000, mean=1/(3*365), sd=1/(3*365*10)) ##mua random norm
dh=1/(3*365) ##dh constant


params1=as.matrix(cbind(P,H,f,phi,mu,ze,betha,k,msr,sig,alfa,q,Th, lambda,mua,dh))

## next create a data frame of constant parameters
constantparams<- cbind(
 #H = rep(1000,sims),  #human population abundance

#  snail dynamics
 f 	= rep(0.16,sims), # instantaneous intrinsic fertility rate of snails (day^-1)
 phi= rep((1-1/(80*0.16))/10000,sims),	# density-dependent parameter for fertility 
 mu	= rep(1/80,sims), 	# snail mortality (mean life expectancy: 50days)
 ze = rep(.5,sims), 		# fraction of exposed snails that reproduce

# epidemiological parameters for infected snails
 betha 	= rep((4*10^-6),sims),	# per capita infection rate
 k		= rep(1,sims),		# clumping parameter of teh negative binomila distribution
 msr	= rep(0.5,sims),		# miracidia shedding rate per reproductive female (devided by miracidia mortality rate)
 sig	= rep(1/50,sims),		# rate from exposed to infected
 alfa 	= rep(1/20,sims),		# mortality rate of infected snail (if negative means that castrated snail increase their life expectancy - tyet it could be zero)


## parameters for the adult worms dynamics in the human hosts
 lambda	= rep(.00005,sims),		# cercaria shadding rate divided by cercaria mortality
 mua  	= rep(1/(3*365),sims),# adult worm natural mortality =1/3 years
 dh		= rep(1/(60*365),sims),# mortality/fertility rate of human population

## predation rates
 q	= rep(0.003,sims)	,	# predation rate on snails
 Th = rep(0.1, sims)  # handling time parameter
 
) # end constant parameters


#######################################################################################
## New simulation framework - for each row in the parameters, run all leves of prawn density to get PDE and PDSE (prawn density for elimination and prawn density for snail eradication)
sims=1000
## Creat matrix with sims combinations of the parameters other than prawns
H=rnorm(sims, mean=1000, sd=100)  ## human population size random 
Carrycapacity=rnorm(sims, mean=10000, sd = 1000)
phi=1/Carrycapacity ## ~inverse carrying capacity random
f=rnorm(sims, mean=0.16, sd= 0.16/10) ## fecundity random rnorm 0.16 with 10% CVmu=rnorm(sims, mean=0.0125, sd = 0.0125/10) ## mu random positive
betha=rnorm(sims, mean=0.000004, sd=0.000004/10) ## transmission rate man-snail random
k=runif(sims, 0.2, 0.7) ## clumping parameter random uniform in empirically derived bounds
msr=rnorm(sims, mean=0.8, sd = 0.8/10) ## msr random 
ze=rnorm(sims, mean=0.5, sd= 0.5/10) ## fraction of exposed that reproduce
sig=rnorm(sims, mean=0.02, sd = 0.02/10) ## sig random norm
alfa=rnorm(sims, mean=0.05, sd = 0.05/10) ## alfa random norm
q=runif(sims, 0.0003, 0.003) ## q random unif from 10X lower to equal to that found in lab
maxconsump=rnorm(sims, mean=7.9, sd=1.02) ## random normal based on empirical lab studies on satiation
Th=1/maxconsump
lambda=rnorm(sims, mean=0.00005, sd = 0.00005/10) ## lambda random norm
mua=rnorm(sims, mean=1/(3*365), sd=1/(3*365*10)) ##mua random norm
dh=1/(3*365) ##dh constant


params2=as.matrix(cbind(H,f,phi,mu,ze,betha,k,msr,sig,alfa,q,Th, lambda,mua,dh))


######################################################################
######################################################################
## Monte Carlo sensitivity analysis - all levels of prawns of 1000 different combinations of parameters
######################################################################

MC2prawn <- function (params=params2, sims=sims) 
### MCprawn() explores the sensitivity of 
###   const.prawn.pop.model() to parameter combinations 
###   using a monte carlo method.
###   It takes a dataframe of parameter combinations and 
###   outputs a dataframe in which each row represents a 
###   different parameter combination together with its 
###   resulting simulation output from the program 
###   const.prawn.pop.model().
###
### usage: first, create the dataframe params0 by hand or 
###   using code such as
###   params0 = expand.grid(m=c(0.1,0.2,0.3), nu=c(0.1,0.2))
###   simprawnout = simprawn(params0)
###   Now run stats on, plot etc, simprawnout
###
### version 0.1  15 April 2014
### Sanna Sokolow sokolow@stanford.edu and Giulio DeLeo 
### deleo@stanford.edu
###
{ 	
	## first we create vectors to store simulation results	
	PDE=numeric(sims)
	PDSE=numeric(sims)
    analytic_pdse=numeric(sims)
    
    stop=500
	
	## each parameter combination is run through 
	## const.prawn.pop.model:
	for(i in 1:sims)
	{
		maxprawn=as.integer((params[i,"f"]-params[i,"mu"])/params[i,"q"])
		analytic_pdse[i] = maxprawn
		
		PD4E = FALSE
		
		for(r in 1:stop)
		{
			P=r	
			# simout contains the entire output corresponding 
			# to one row of parameters. We then pull out 
			# interesting data from it.
		
			paramsi<-c(params[i,], P=P)
			# now integrate the system
			out <- as.matrix(
 			ode(
   			func		= const.prawn.pop.model, 	y = xstart,
   				times	= times,
   				parms	= paramsi		
   				))
   		
		if(out[10951,"I"]	< 1 &  PD4E == FALSE) 
		   {
			PDE[i]=P
			PD4E = TRUE
			}

		if((out[10951,"S"]+out[10951,"I"]+out[10951,"E"])	< 1) 
		   {
			PDSE[i]=P
			break
			
			}	
				
		PDSE[i]=maxprawn
		}
		}
	
	simprawnout= cbind(params,PDE, PDSE)
	return (simprawnout)
}

######################################################################
######################################################################
## Monte Carlo sensitivity analysis - all levels of prawns of 1000 different combinations of parameters
######################################################################

MC3prawn <- function (params=params2, sims=sims) 
### MCprawn() explores the sensitivity of 
###   const.prawn.pop.model() to parameter combinations 
###   using a monte carlo method.
###   It takes a dataframe of parameter combinations and 
###   outputs a dataframe in which each row represents a 
###   different parameter combination together with its 
###   resulting simulation output from the program 
###   const.prawn.pop.model().
###
### usage: first, create the dataframe params0 by hand or 
###   using code such as
###   params0 = expand.grid(m=c(0.1,0.2,0.3), nu=c(0.1,0.2))
###   simprawnout = simprawn(params0)
###   Now run stats on, plot etc, simprawnout
###
### version 0.1  15 April 2014
### Sanna Sokolow sokolow@stanford.edu and Giulio DeLeo 
### deleo@stanford.edu
###
{ 	
	## first we create vectors to store simulation results	
	PDE=numeric(sims)
	PDSE=numeric(sims)
    pde=numeric(0)	
	
	## each parameter combination is run through 
	## const.prawn.pop.model:
	for(i in 1:sims)
	{
		
		maxprawn=as.integer((params[i,"f"]-params[i,"mu"])/params[i,"q"])
		stop=500
		for(r in 1:stop)
		{
			P=r
			
			{
			# simout contains the entire output corresponding 
			# to one row of parameters. We then pull out 
			# interesting data from it.
		
			paramsi<-c(params[i,], P=P)
			# now integrate the system
			out <- as.matrix(
 			ode(
   			func		= const.prawn.pop.model,
   				y		= xstart,
   				times	= times,
   				parms	= paramsi		
   				))
   		
		if((out[10951,"S"]+out[10951,"I"]+out[10951,"E"])	< 1) 
		{
			
			break
			
			}
			PDE[i]=P
		}		
		PDSE[i]=maxprawn
		}
		}
	
	simprawnout= cbind(params,PDE, PDSE)
	return (simprawnout)
}

