#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#get p(e) for each parameter set
pe = as.numeric()
eps.mean = as.numeric()
eps.sd = as.numeric()
for(s in 1:nrow(par.mat)){
  pe[s] = sum(fill.arr[which(par.mat[s,1] == lam.range), which(par.mat[s,2] == kap.range), 1, ]) / stoch.sims
  eps.mean[s] = mean(fill.arr[which(par.mat[s,1] == lam.range), which(par.mat[s,2] == kap.range), 2, ])
  eps.sd[s] = sd(fill.arr[which(par.mat[s,1] == lam.range), which(par.mat[s,2] == kap.range), 2, ])
}

plot(eps.mean, pe, pch = 16, cex = 0.6,
     xlab = expression(epsilon), ylab = expression(italic('P(e)')))