#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#This code adapted from Arathi's "schisto_halstead_2pops_mda_bouncebackRate_cloudParam.R"

# This code seeks to test the bounce back rate concept. Given the worm burden from a model, can we determine if
# it has positive DD and hence possible to eliminate? Or has no DD and is harder to eliminate.
# We first start with a deterministic model
# We will then fit it to the Lampsar II data to get a collection of parameter values that fit max likelihood bounds
# Each parameter value will give a unique W for pre and post mda in this model and a unique bounce back rate curve
# If there are 100 parameter values for example, how many have Reff<1 after 10 rounds?
##################################################################################################

require(deSolve)
require(graphics)
require(ggplot2)

source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/Reff_fn_lib_CH.R")
source("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/lib_schistoModels_DDandNODD_CH.R")

outputfile<-"C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Monash/shortlist_51.csv"
shortlist_first100<-as.data.frame(read.csv(file=outputfile, header=TRUE, sep=","))
sort_ind<-order(shortlist_first100$R0, decreasing=FALSE)
shortlist_first100<-shortlist_first100[sort_ind,]
  shortlist_first100 = shortlist_first100[c(1:100),] #Get rid of other parameter estimates
  
eqbm = data.frame(eqbm.s = 0, eqbm.e = 0, eqbm.i = 0,
                       eqbm.wt = 0, eqbm.wu = 0)
eqbm.pdd = data.frame(eqbm.pdd.s = 0, eqbm.pdd.e = 0, eqbm.pdd.i = 0,
                           eqbm.pdd.wt = 0, eqbm.pdd.wu = 0)
#Fill shortlist parameter vector with equilibrium estimates
  k.fit = 0.08 #clumping parameter derived from fit to epi data
  cov = 0.8 #80% coverage assumed
  
for(i in 1:nrow(shortlist_first100)){
    params<-parameters_2pops_mda_Chris1
    
    params["beta"]<-shortlist_first100$beta[i]
    params["lamda"]<-shortlist_first100$lamda.twa[i] 
    params['cov']<-cov
    
    time<-seq(from=0, to=200*365, by=1)
    nstart=c(S=3892,E=3750,I=1200, Wt=50, Wu=50)
    
    params['k'] = 0
      output<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))

    params['k'] = k.fit  
      output.PDD<-as.data.frame(ode(nstart,time,schisto_halstead_2pops_mda,params))
    
    baseline.PDD<-output.PDD[dim(output.PDD)[1],]
    baseline<-output[dim(output)[1],]
    
    eqbm[i,] = baseline[,c(2:6)]
    eqbm.pdd[i,] = baseline.PDD[,c(2:6)]
        
    print(i)

  }

# next create 100 models with these beta, lamda combinations and run with DD, get 100 BBrate curves
cov = 0.8  # 80% coverage
eff = 0.99 # 99% efficacy
mda.years = c(1:20) # annual MDA for 20 years
time = c(1:(365*61)) #1 year run up, 20 years annual mda followed by 40 years rebound with no MDA

#1 annual round of MDA
mda.events = data.frame(var = rep('Wt', length(mda.years)),
                         time = c(mda.years*365),
                         value = rep((1 - eff), length(mda.years)),
                         method = rep('mult', length(mda.years)))

#double dose annual MDA separated by 3 weeks
mda.events2 = data.frame(var = rep('Wt', length(mda.years*2)),
                         time = sort(c(c(mda.years*365), c(mda.years*365)+21)),
                         value = rep((1 - eff), length(mda.years*2)),
                         method = rep('mult', length(mda.years*2)))

#Objects to fill in loop over 100 parameter sets
#array with all 100 60-year runs with and w/o PDD
  det.runs.k008.mda1 = array(data = NA, dim = c(length(time), ncol(eqbm)+1, nrow(shortlist_first100), 2))
#Vector of W-pre values
  w.pre.k008.mda1 = array(data = NA, dim = c(nrow(shortlist_first100), length(mda.years)))
  w.pos.k008.mda1 = array(data = NA, dim = c(nrow(shortlist_first100), length(mda.years)))
  w.pre.k008.mda1.pdd = array(data = NA, dim = c(nrow(shortlist_first100), length(mda.years)))
  w.pos.k008.mda1.pdd = array(data = NA, dim = c(nrow(shortlist_first100), length(mda.years)))


for(y in 1:nrow(shortlist_first100)){
  params<-parameters_2pops_mda_Chris1
  
    params["beta"]<-shortlist_first100$beta[y]
    params["lamda"]<-shortlist_first100$lamda.twa[y] 
    params["cov"]<-cov
    params['k']<-0
    
  nstart = c(S=eqbm$eqbm.s[y], 
             E=eqbm$eqbm.e[y], 
             I=eqbm$eqbm.i[y], 
             Wt=eqbm$eqbm.wt[y], 
             Wu=eqbm$eqbm.wu[y])
  
  det.runs.k008.mda1[, , y, 1] = ode(nstart, time, schisto_halstead_2pops_mda, params,
                           events = list(data = mda.events))
    for(i in mda.years){
      w.pre.k008.mda1[y,i] = cov*det.runs.k008.mda1[(i*365), 5, y, 1] + (1-cov) * det.runs.k008.mda1[(i*365), 6, y, 1]
      w.pos.k008.mda1[y,i] = cov*det.runs.k008.mda1[(i*365+1), 5, y, 1] + (1-cov) * det.runs.k008.mda1[(i*365+1), 6, y, 1]
    } 
  
  params['k']<-k.fit
  
  nstart.pdd = c(S=eqbm.pdd$eqbm.pdd.s[y], 
                 E=eqbm.pdd$eqbm.pdd.e[y], 
                 I=eqbm.pdd$eqbm.pdd.i[y], 
                 Wt=eqbm.pdd$eqbm.pdd.wt[y], 
                 Wu=eqbm.pdd$eqbm.pdd.wu[y])
  
  det.runs.k008.mda1[, , y, 2] = ode(nstart.pdd, time, schisto_halstead_2pops_mda, params,
                           events = list(data = mda.events))
  
  for(i in mda.years){
    w.pre.k008.mda1.pdd[y,i] = cov*det.runs.k008.mda1[(i*365), 5, y, 2] + (1-cov) * det.runs.k008.mda1[(i*365), 6, y, 2]
    w.pos.k008.mda1.pdd[y,i] = cov*det.runs.k008.mda1[(i*365+1), 5, y, 2] + (1-cov) * det.runs.k008.mda1[(i*365+1), 6, y, 2]
  } 
  print(y)
}
  
#Some post processing
  mean.w.k008.mda1 = matrix(ncol = 3, nrow = length(time))
    mean.w.k008.mda1[,1] = rowMeans(det.runs.k008.mda1[ , 5, , 1])
    mean.w.k008.mda1[,2] = rowMeans(det.runs.k008.mda1[ , 6, , 1])
    mean.w.k008.mda1[,3] = cov * mean.w.k008.mda1[,1] + (1-cov) * mean.w.k008.mda1[,2]
    
  mean.w.pdd.k008.mda1 = matrix(ncol = 3, nrow = length(time))
    mean.w.pdd.k008.mda1[,1] = rowMeans(det.runs.k008.mda1[ , 5, , 2])
    mean.w.pdd.k008.mda1[,2] = rowMeans(det.runs.k008.mda1[ , 6, , 2])
    mean.w.pdd.k008.mda1[,3] = cov * mean.w.pdd.k008.mda1[,1] + (1-cov) * mean.w.pdd.k008.mda1[,2]
   
  plot(time/365, mean.w.pdd.k008.mda1[,3], type = 'l', lwd = 2, xlim = c(0,61), ylim = c(0,100),
       xlab = 'time (yrs)', ylab = 'W')
    lines(time/365, mean.w.pdd.k008.mda1[,1], lty=2)
    lines(time/365, mean.w.pdd.k008.mda1[,2], lty=3)
    
  plot(time/365, mean.w.k008.mda1[,3], type = 'l', lwd = 2, xlim = c(0,61), ylim = c(0,100),
       xlab = 'time (yrs)', ylab = 'W')
    lines(time/365, mean.w.k008.mda1[,1], lty=2)
    lines(time/365, mean.w.k008.mda1[,2], lty=3)
    