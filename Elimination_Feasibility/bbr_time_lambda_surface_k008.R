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

source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/Reff_fn_lib_CH.R")
source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/lib_schistoModels_DDandNODD_CH.R")

#parameter vectors and values ######
params = parameters_2pops_mda_Chris1
  params['beta'] = 1.63e-6
  cov = 0.8
  params['cov'] = cov

lam.range = seq(1e-4, 3e-4, length.out = 50)  #Range of transmission intensities to test

teq<-seq(from=0, to=200*365, by=10)            #time sequence to run to eqbm
nstart=c(S=3892,E=3750,I=1200, Wt=50, Wu=50)  #nstart vector for initialization
time = c(1:(365*22))                          #time series to simulate MDA
eff = 0.99                                    # 99% efficacy

mda.years = c(1:20) # annual MDA for 20 years
mda.events = data.frame(var = rep('Wt', length(mda.years)),
                        time = c(mda.years*365),
                        value = rep((1 - eff), length(mda.years)),
                        method = rep('mult', length(mda.years)))

#BBR funtion ######
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
#BBR profile simulator #####
bbr.prof = function(lam){
  #Update parameter estimates for simulation
  params["lamda"] = lam 
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
  bbr.est(as.data.frame(ode(y = eqbm2, 
                            times = time, 
                            func = schisto_halstead_2pops_mda, 
                            parms = params,
                            events = list(data = mda.events))))
  
}

#simulate bbr profile across range of transmission intensities with k = 0.08 (best fit value)
  bbrm = matrix(ncol = length(lam.range), nrow = length(mda.years)-1)
  for(l in 1:length(lam.range)){
    bbrm[, l] = bbr.prof(lam = lam.range[l])
    print(l)
  }
  
colnames(bbrm) = round(lam.range, digits = 6)
rownames(bbrm) = mda.years[-20]

write.csv(bbrm, 'Elimination_Feasibility/bbr_surface_time_lambda_k008.csv', row.names = F)