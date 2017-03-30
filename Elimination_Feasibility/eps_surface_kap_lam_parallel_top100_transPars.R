#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

##################################################################################################
#Packages and source other scripts for model and parameter fitting values #####
require(deSolve)
require(graphics)
require(parallel)

source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/Reff_fn_lib_CH.R")
source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/lib_schistoModels_DDandNODD_CH.R")

no.cores = detectCores() - 1

#Get best fit parameter values ###### 
outputfile<-"C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/shortlist_51.csv"
  best100<-as.data.frame(read.csv(file=outputfile, header=TRUE, sep=","))
    sort_ind<-order(best100$negLL, decreasing=T)
    best100<-best100[sort_ind,]
      best100 = best100[c(1:100),] #Get rid of other parameter estimates; keep 100 parameter sets with max(negLL)

#Vectors for simulation and matrix to fill with epsilon values ############      
beta.vec = rep(best100$beta, times = 100)
lam.vec = rep(best100$lamda.twa, times = 100)
r0.vec = rep(best100$R0, times = 100)
kap.vec = rep(seq(0,3, length.out = 100), each = 100)
eps.vec = 0
  vecs = cbind(beta.vec, lam.vec, r0.vec, kap.vec, eps.vec) 
  
  teq<-seq(from=0, to=200*365, by=10)            #time sequence to run to eqbm
  nstart=c(S=3892,E=3750,I=1200, Wt=50, Wu=50)  #nstart vector for initialization
  time = c(1:(365*22))                          #time series to simulate MDA
  cov = 0.8
  parameters_2pops_mda_Chris1["cov"] = cov      # 80% coverage
  eff = 0.99                                    # 99% efficacy
  params = parameters_2pops_mda_Chris1
  
  mda.years = c(1:20) # annual MDA for 20 years
  mda.events = data.frame(var = rep('Wt', length(mda.years)),
                          time = c(mda.years*365),
                          value = rep((1 - eff), length(mda.years)),
                          method = rep('mult', length(mda.years)))
  
#bbr function ########
  pre.days = mda.years * 365
  pos.days = mda.years * 365 +1    
      
  bbr.est = function(input){ #input is a data.frame output from ode solver
    
    bbr = (1/((cov*input$Wt[input$time %in% pos.days] + 
             (1-cov) * input$Wu[input$time %in% pos.days])[c(1:(length(pre.days)-1))])) *
          (((cov*input$Wt[input$time %in% pre.days] + 
             (1-cov) * input$Wu[input$time %in% pre.days])[c(2:length(pre.days))]) - 
          ((cov*input$Wt[input$time %in% pos.days] + 
             (1-cov) * input$Wu[input$time %in% pos.days])[c(1:(length(pre.days)-1))]))
    
    bbr
  }  

#Epsilon estimator given transmission parameters and clumping parameter ######   
eps.est.beta = function(lam, beta, kap){
#Update parameter estimates for simulation
  params["lamda"] = lam 
  params["beta"] = beta
  params['k'] = kap
#run to equilibrium  
  output<-as.data.frame(ode(y = nstart,
                            times = teq,
                            func = schisto_halstead_2pops_mda,
                            parms = params))
  
  eqbm2 = c(S=output[dim(output)[1],2], 
            E=output[dim(output)[1],3], 
            I=output[dim(output)[1],4], 
            Wt=output[dim(output)[1],5], 
            Wu=output[dim(output)[1],6])
  
#This simulates mda across mda years, calculates bbr each year, and estimates epsilon
  #from the generated data.frame
  as.numeric(lm(bbr.est(as.data.frame(ode(y = eqbm2, 
                                          times = time, 
                                          func = schisto_halstead_2pops_mda, 
                                          parms = params,
                                          events = list(data = mda.events)))) ~ 
                  c(1:(length(mda.years)-1)))$coefficients[2])
  
}

clusteps = makeCluster(no.cores)
  clusterExport(cl = clusteps, 
                varlist = c('params','lam.vec', 'beta.vec', 'kap.vec', 'ode', 'vecs',  'mda.years',
                            'mda.events', 'nstart', 'teq', 'schisto_halstead_2pops_mda',
                            'bbr.est', 'pre.days', 'pos.days', 'time', 'cov'))  
#Below used to test function as it was being tweaked to minimize computation time
#  system.time(clusterMap(clusteps, 
#                         fun = eps.est.beta, 
#                         lam = lam.vec, 
#                         beta = beta.vec, 
#                         kap = kap.vec, 
#                         SIMPLIFY = TRUE))
  
  
  vecs[,5] = clusterMap(clusteps, 
                        fun = eps.est.beta, 
                        lam = lam.vec, 
                        beta = beta.vec, 
                        kap = kap.vec, 
                        SIMPLIFY = TRUE)

stopCluster(clusteps)

write.csv(vecs, 'Elimination_Feasibility/eps_estimates_100bestfit_kap0-3.csv',
          row.names = FALSE)