#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(deSolve)
require(sensitivity)
source('Review_models/agroReview_mod1.R')

#Get parameter sets to sample from #######################
  sims = 50
#parameters ranges
  f_N.range<-seq(0.05, 0.45, length.out = sims)
  phi_N.range<-seq(10, 70, length.out = sims)
  z.range<-seq(0.2, 1.0, length.out = sims)
  mu_N.range<-seq(1/20, 1/100, length.out = sims)
  mu_I.range<-seq(1/2, 1/20, length.out = sims)
  f_P.range<-seq(0.001, 0.2, length.out = sims)
  phi_P.range<-seq(0.1, 5.0, length.out = sims)
  mu_P.range<-seq(1/100, 1/(365*2), length.out = sims)
  alpha.range<-seq(0.00001, 0.3, length.out = sims)
  Th.range<-seq(1/30, 1/5, length.out = sims)
  nn.range<-seq(1, 3, length.out = sims)
  m.range<-seq(30, 600, length.out = sims)
  v.range<-seq(0.02, .20, length.out = sims)
  pi_M.range<-seq(0, 1.0, length.out = sims)
  theta.range<-seq(10, 500, length.out = sims)
  pi_C.range<-seq(0, 1.0, length.out = sims)
  beta.range<-seq(1e-7, 1e-3, length.out = sims)
  sigma.range<-seq(1/10, 1/60, length.out = sims)
  lamda.range<-seq(1e-6, 1e-3, length.out = sims)
  k.range<-seq(0.01, 1.0, length.out = sims)
  mu_W.range<-seq(1/365, 1/(365*10), length.out = sims)
  mu_H.range<-seq(1/(365*10), 1/(365*80), length.out = sims)
  
#parameter names
  paranges<-cbind("f_N" = f_N.range,
                  "phi_N" = phi_N.range,
                  "z" = z.range,
                  "mu_N" = mu_N.range,
                  "mu_I" = mu_I.range,
                  "f_P" = f_P.range,
                  "phi_P" = phi_P.range,
                  "mu_P" = mu_P.range,
                  "alpha" = alpha.range,
                  "Th" = Th.range,
                  "nn" = nn.range,
                  "m" = m.range,
                  "v" = v.range,
                  "pi_M" = pi_M.range,
                  "theta" = theta.range,
                  "pi_C" = pi_C.range,
                  "beta" = beta.range,
                  "sigma" = sigma.range,
                  "lamda" = lamda.range,
                  "k" = k.range,
                  "mu_W" = mu_W.range,
                  "mu_H" = mu_H.range)

constantparams<-matrix(ncol = length(parameters), nrow = sims)

for(i in 1:length(parameters)){
  constantparams[,i] = rep(parameters[i],sims)
}

colnames(constantparams)<-names(parameters)
  vars<-colnames(paranges)
  
#First check scatter plots of outcomes across each parameter range while holding other parameters equal ############
  outputfill1<-matrix(ncol = length(vars), nrow = sims)
  outputfill2<-matrix(ncol = length(vars), nrow = sims)
  outputfill3<-matrix(ncol = length(vars), nrow = sims)
  
  senstart = nstart3
  senstart['Q1'] = 0 #Assessing model sensitivities without agrochemicals
  time = seq(0, 365*25, 1)
  
  for(j in 10:length(vars)){
    for(i in 1:sims){
      print(c(j, i))
      
      other<-constantparams[, -which(colnames(constantparams) %in% vars[j])]
      test<-paranges[, which(colnames(paranges) %in% vars[j])]
      
      parametersuse<-cbind(other, test) 
      colnames(parametersuse)[dim(parametersuse)[2]]<-vars[j]
      
      prams<-parametersuse[i,]
      
      outputest = as.data.frame(ode(senstart, time, mod1, prams))
      
      outputfill1[i,j] = sum(outputest$S[dim(outputest)[1]] + outputest$E[dim(outputest)[1]] + outputest$I[dim(outputest)[1]])
      outputfill2[i,j] = outputest$I[dim(outputest)[1]]
      outputfill3[i,j] = outputest$W[dim(outputest)[1]]
      
    }
    par(mfrow = c(3,1), mar = c(4,3.75,1,0.4)+0.1)
    
    plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill1[,j], 
         xlab = vars[j], ylab = 'snail pop',
         pch = 16, cex = 0.75, col = 'blue', 
         ylim = c(0,75))
    
    plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill2[,j], 
         xlab = vars[j], ylab = 'infecteds',
         pch = 16, cex = 0.75, col = 'red', 
         ylim = c(0,20))
    
    plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill3[,j], 
         xlab = vars[j], ylab = 'mean worm burden',
         pch = 16, cex = 0.75, col = 'purple', 
         ylim = c(0,300))
  }