#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Load datasets from simulations with no observation noise #########
library(tidyverse)

load("Elimination_Feasibility/Organize/Models/Outputs_Refs/model_fit_profile_likelihood_parameters.Rdata")
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/eqbm_stoch_sims2.Rdata")  

par.sims = 50
  lam.range = seq(min(fin_pars95ci$lamda_twa), max(fin_pars95ci$lamda_twa), length.out = par.sims)  #transmission intensity range
  kap.range = seq(0, 2, length.out = par.sims)        #Pos. density dependence range

  par.mat = expand.grid(kap = kap.range, lam = lam.range)
  
  par.mat <- cbind(par.mat, stoch_sims_eqbm_vals_fin)
  
  par.mat$kap.sim = NA     #add column for kappa used in the simulation
  par.mat$lam.sim = NA     #add column for lambda used in the simulation 
  par.mat$eps = NA         #add column for elim. feas estimator (eps) derived from worm burden data
  par.mat$eps.sd = NA      #add column for elim. feas estimator st. dev 
  par.mat$eps.prev = NA    #add column for elim. feas estimator (eps) derived from prevalence data
  par.mat$eps.prev.sd = NA #add column for elim. feas estimator (eps) derived from prevalence data st. dev
  par.mat$pe = NA          #add column for proba elimination (P(e))
  
#Use array loads from sim runs on savio to fill the data frame #########
#In each array dimensions are as follows:
# [here, , ] is kappa and lamda parameters used in the given set of simulations
# [ , here, ] is the 5 outputs of interest: 1: kappa parameter used, 2: lamda parameter used, 3: resulting p(e), 4: epsilon estimate from worm burden data, 5: epsilon estimate from prevalence data
# [, , here] is the result from all 1000 runs of the stochastic simulation with the given set of parameters
load('Elimination_Feasibility/Organize/Models/Outputs_Refs/stoch_sims1_500.Rdata')
  inarr = fill.arr
  chunk = 1
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$kap.sim[m+lo-1] = mean(inarr[m, 1, ], na.rm = TRUE)
    par.mat$lam.sim[m+lo-1] = mean(inarr[m, 2, ], na.rm = TRUE)
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ], na.rm = TRUE)
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ], na.rm = TRUE)
    par.mat$eps.prev[m+lo-1] = mean(inarr[m, 5, ], na.rm = TRUE)
    par.mat$eps.prev.sd[m+lo-1] = sd(inarr[m, 5, ], na.rm = TRUE)
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }
  
load('Elimination_Feasibility/Organize/Models/Outputs_Refs/stoch_sims501_1000.Rdata')
  inarr = fill.arr
  chunk = 2
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$kap.sim[m+lo-1] = mean(inarr[m, 1, ], na.rm = TRUE)
    par.mat$lam.sim[m+lo-1] = mean(inarr[m, 2, ], na.rm = TRUE)
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$eps.prev[m+lo-1] = mean(inarr[m, 5, ], na.rm = TRUE)
    par.mat$eps.prev.sd[m+lo-1] = sd(inarr[m, 5, ], na.rm = TRUE)
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }
  
load('Elimination_Feasibility/Organize/Models/Outputs_Refs/stoch_sims1001_1500.Rdata')
  inarr = fill.arr
  chunk = 3
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$kap.sim[m+lo-1] = mean(inarr[m, 1, ], na.rm = TRUE)
    par.mat$lam.sim[m+lo-1] = mean(inarr[m, 2, ], na.rm = TRUE)
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$eps.prev[m+lo-1] = mean(inarr[m, 5, ], na.rm = TRUE)
    par.mat$eps.prev.sd[m+lo-1] = sd(inarr[m, 5, ], na.rm = TRUE)
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }
  
load('Elimination_Feasibility/Organize/Models/Outputs_Refs/stoch_sims1501_2000.Rdata')
  inarr = fill.arr
  chunk = 4
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$kap.sim[m+lo-1] = mean(inarr[m, 1, ], na.rm = TRUE)
    par.mat$lam.sim[m+lo-1] = mean(inarr[m, 2, ], na.rm = TRUE)
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$eps.prev[m+lo-1] = mean(inarr[m, 5, ], na.rm = TRUE)
    par.mat$eps.prev.sd[m+lo-1] = sd(inarr[m, 5, ], na.rm = TRUE)
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }
  
load('Elimination_Feasibility/Organize/Models/Outputs_Refs/stoch_sims2001_2500.Rdata')
  inarr = fill.arr
  chunk = 5
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$kap.sim[m+lo-1] = mean(inarr[m, 1, ], na.rm = TRUE)
    par.mat$lam.sim[m+lo-1] = mean(inarr[m, 2, ], na.rm = TRUE)
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$eps.prev[m+lo-1] = mean(inarr[m, 5, ], na.rm = TRUE)
    par.mat$eps.prev.sd[m+lo-1] = sd(inarr[m, 5, ], na.rm = TRUE)
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }
  
#Quick sanity check to make sure that parameters used in simulations all match up
  sum(par.mat$kap == par.mat$k) == 2500 
  sum(par.mat$kap == par.mat$kap.sim) == 2500

par.mat2 = par.mat 
#plots of and analysis of pe(e) / eps #######
#plot P(e) across eps
par.mat2 %>% 
  ggplot(aes(x = eps, y = pe)) + geom_point(pch = 18, cex = 0.8) + 
    labs(x = expression(epsilon), 
         y = expression(italic('P(e)'))) +
    theme_bw()

#Save plot as tiff
tiff("Elimination_Feasibility/plots/PLoS_Figs_PostReview/Fig5.tiff", 
     height = 13.2, width = 13.2, units = 'cm', compression = "lzw", res = 300) 
par.mat2 %>% 
  ggplot(aes(x = eps, y = pe)) + geom_point(pch = 18, cex = 0.8) + 
    labs(x = expression(epsilon), 
         y = expression(italic('P(e)'))) +
    theme_bw()
dev.off()

par.mat2 %>% 
  ggplot(aes(x = eps, y = pe, col = lambda)) + geom_point(pch = 18, cex = 0.8) + 
    labs(x = expression(epsilon), 
         y = expression(italic('P(e)')),
         col = expression(lambda)) +
    scale_color_gradient(low = 'wheat1', high = 'red') +
    theme_bw()

par.mat2 %>% 
  ggplot(aes(x = eps, y = pe, col = kap)) + geom_point(pch = 18, cex = 0.8) + 
    labs(x = expression(epsilon), 
         y = expression(italic('P(e)')),
         col = expression(kappa)) +
    scale_color_gradient(low = 'wheat1', high = 'red') +
    theme_bw()

par.mat2 %>% 
  ggplot(aes(x = eps.prev, y = pe)) + geom_point(pch = 18, cex = 0.8) + 
    labs(x = expression(epsilon), 
         y = expression(italic('P(e)'))) +
    theme_bw()
  
par.mat2 %>% 
  ggplot(aes(x = eps.prev, y = pe, col = lambda)) + geom_point(pch = 18, cex = 0.8) + 
    labs(x = expression(epsilon), 
         y = expression(italic('P(e)')),
         col = expression(lambda)) +
    scale_color_gradient(low = 'wheat1', high = 'red') +
    theme_bw()

#prevalence based
plot(par.mat2$eps.prev, par.mat2$pe, pch = 18, cex = 0.6, ylim = c(0,1),
     xlab = expression(epsilon), ylab = expression(italic('P(e)')))

#ggplot to show val of lambda / kappa in addition to scatterP
par.mat2 %>% 
  ggplot(aes(x = eps, y = pe, col = lambda)) +
    theme_bw() +
    geom_point(shape = 18, size = 1.2) +
    scale_color_gradient(low = 'wheat1', high = 'red') +
    labs(col = expression(lambda),
         x = expression(epsilon),
         y = expression(italic(P(e))))

par.mat2 %>% 
  ggplot(aes(x = eps.prev, y = pe, col = lambda)) +
    theme_bw() +
    geom_point(shape = 18, size = 1.2) +
    scale_color_gradient(low = 'wheat1', high = 'red') +
    labs(col = expression(lambda),
         x = expression(epsilon),
         y = expression(italic(P(e))))

#only plot for vals where p(e) != 1 or 0
par.mat0_1 <- par.mat2 %>% 
  filter(pe < 1 & pe >0 & kap != 0)

plot(par.mat0_1$eps, par.mat0_1$pe, pch = 18, cex = 0.6, ylim = c(0,1),
     xlab = expression(epsilon), ylab = '',
     cex.lab = 1.4)
  mtext(side = 2, text = expression(italic('P(e)')), line = 2.4, cex = 1.2)

#logit transform p(e)
plot(par.mat0_1$eps, log(par.mat0_1$pe / (1 - par.mat0_1$pe)), 
     pch = 18, cex = 0.6, #ylim = c(0,1),
     xlab = expression(epsilon), ylab = expression(italic('P(e)')))

#regression of eps to pe ########
#with manual logit transform
pe.mod = lm(log(pe / (1 - pe)) ~ eps, data = par.mat0_1)
  summary(pe.mod)
  abline(a = pe.mod$coefficients[1], b = pe.mod$coefficients[2],
         lty = 2, col = 2, lwd = 2)
#weighted
pe.mod.weighted = lm(log(pe / (1 - pe)) ~ eps, data = par.mat0_1, weights = eps.sd)
  summary(pe.mod.weighted)
  #abline(a = pe.mod.weighted$coefficients[1], b = pe.mod.weighted$coefficients[2],
  #       lty = 2, col = 3, lwd = 2)
  
#Doesn't look like adding weights changes much, so stick with unweighted regression  
  
#including transmission intensity and pdd as regression coefficients
pe.mod.all =  lm(log(pe / (1 - pe)) ~ eps + kap + lam, data = par.mat0_1) 
  summary(pe.mod.all)
AIC(pe.mod, pe.mod.all)  #looks like this actually produces a better model... 
  
#plot regression prediction with actual data #########
  plot(par.mat0_1$eps, par.mat0_1$pe, pch = 18, cex = 0.6, ylim = c(0,1),
       xlab = expression(epsilon), ylab = expression(italic('P(e)')))
  
  back.trans = function(eps){
    p = exp(pe.mod$coefficients[1] + eps*pe.mod$coefficients[2])
    pe.pred = p / (1+p)
    return(pe.pred)
  }
    lines(seq(-0.05, 0.02, by = 0.001), sapply(seq(-0.05, 0.02, by = 0.001), back.trans, simplify = T),
          lty = 2,col = 2, lwd = 2)
#other plots #######    
#plot relationship between eps and its standard deviation
plot(par.mat2$eps, par.mat2$eps.sd, pch = 18, cex = 0.6, #ylim = c(0,1),
     ylab = expression(paste(epsilon, ' St. Dev.', sep='')), 
     xlab = expression(epsilon))

#p(e) across transmission intensity
plot(par.mat2$lam, par.mat2$pe, pch = 18, cex = 0.6, ylim = c(0,1),
     xlab = expression(lambda), ylab = expression(italic('P(e)')))

#p(e) across pos. dens. dep.
plot(par.mat2$kap, par.mat2$pe, pch = 18, cex = 0.6, ylim = c(0,1),
     xlab = expression(kappa), ylab = expression(italic('P(e)')))

#surface of p(e) across kappa and lambda
lams = unique(par.mat2$lam)
kaps = unique(par.mat2$kap)

lkmat = matrix(par.mat2$pe, ncol=50, nrow=50, byrow = T)

persp(y = lams, ylim = range(lams), x = kaps, xlim = range(kaps),
      z = lkmat, ticktype = 'detailed', nticks = 4, 
      xlab = 'Pos. Density Dependence',
      ylab = 'Transmission Intensity',
      zlab = 'Probability of Elimination',
      phi = 20, theta = 135, shade = 0.4, col = 'lightblue')


#Use array loads from sim runs on savio to create new data frame with each row as a simulation #########
mat.all = data.frame()

load('Elimination_Feasibility/Savio/Results/fill_array1_500.Rdata')
  inarr = fill.arr
  chunk = 1
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(j in 1:500){
    mat.all = rbind(mat.all, t(inarr[j, , ]))
  }
  
load('Elimination_Feasibility/Savio/Results/fill_array501_1000.Rdata')
  inarr = fill.arr
  chunk = 2
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(j in 1:500){
    mat.all = rbind(mat.all, t(inarr[j, , ]))
  }
  
load('Elimination_Feasibility/Savio/Results/fill_array1001_1500.Rdata')
  inarr = fill.arr
  chunk = 3
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(j in 1:500){
    mat.all = rbind(mat.all, t(inarr[j, , ]))
  }
  
load('Elimination_Feasibility/Savio/Results/fill_array1501_2000.Rdata')
  inarr = fill.arr
  chunk = 4
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(j in 1:500){
    mat.all = rbind(mat.all, t(inarr[j, , ]))
  }
  
load('Elimination_Feasibility/Savio/Results/fill_array2001_2500.Rdata')
  inarr = fill.arr
  chunk = 5
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(j in 1:500){
    mat.all = rbind(mat.all, t(inarr[j, , ]))
  }
  
colnames(mat.all) = c('kappa', 'lambda', 'pe', 'eps')

#logistic regression with glm 
pe.glm = glm(pe ~ eps + kappa + lambda, data = mat.all, 
             family = 'binomial'(link = 'logit'))
summary(pe.glm)

