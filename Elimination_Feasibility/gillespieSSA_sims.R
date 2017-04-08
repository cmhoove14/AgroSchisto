#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Load packages and other files ##########

require(deSolve)
require(graphics)
require(GillespieSSA)

source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/Reff_fn_lib_CH.R")
source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/lib_schistoModels_DDandNODD_CH.R")

outputfile<-"C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/shortlist_51.csv"
shortlist_first100<-as.data.frame(read.csv(file=outputfile, header=TRUE, sep=","))
sort_ind<-order(shortlist_first100$R0, decreasing=FALSE)
shortlist_first100<-shortlist_first100[sort_ind,]
  shortlist_first100 = shortlist_first100[c(1:100),] #Get rid of other parameter estimates

phi.params = read.csv('Elimination_Feasibility/phi_parameters_simplified.csv')  

#Other initial dependencies ########
  params = parameters_2pops_mda_Chris1
  cov = 0.8
  params['cov'] = cov
  k.range = phi.params$k
  
  x0.init = c(S=3892,E=3750,I=1200, Wt=50, Wu=50)
  
#SSA function ########
  getSSA<-function(params, x0, time, k){
    params['k'] = k
    params['phi_a'] = phi.params$a[phi.params$k == k]
    params['phi_b'] = phi.params$b[phi.params$k == k]
    
    a1<-"f_N*(1-((S+E+I)/phi_N) )*(S+E)"
    a2<-"0.5*beta*(cov*Wt+((1-cov)*Wu))*H*S*(phi_a-(phi_a/((cov*Wt+(1-cov)*Wu)+1)^phi_b))"
    a3<-"mu_N*E"
    a4<-"(sigma*E)"
    a5<-"( (mu_N+mu_I)*I)"
    a6<-"(lamda*I)" 
    a7<-"( (mu_W+mu_H)*Wt)"
    a8<-"( (mu_W+mu_H)*Wu)"
    a9<-"(eff*Wt*mda )"
    a10<-"(mu_N*S)"
    
    a<-c(a1,a2,a3,a4,a5,a6,a7, a8, a9,a10)
    
    nu <- matrix(c(1,-1, 0, 0, 0, 0, 0, 0, 0,-1,
                   0, 1,-1,-1, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 1,-1, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 1,-1, 0,-1, 0,
                   0, 0, 0, 0, 0, 1, 0,-1, 0, 0),nrow=5,byrow=TRUE)
    
    out <- (ssa(x0,                              # initial state vector
                a,                               # propensity vector (i.e. model)
                nu,                              # state-change matrix
                parms = params,                  # parameter values
                tf=time,                         # final run time
                simName="basic model",           
                censusInterval=0, 
                ignoreNegativeState=TRUE, 
                method="ETL", 
                verbose=F, 
                tau=0.3))$data
    
    
    out
  }
  
getSSA(params = params, x0 = x0.init, time = c(1:(365*2)), k = k.range[5])  