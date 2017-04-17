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
par.mat = read.csv('Elimination_Feasibility/eq_vals_for_trans_pars.csv')   #load parameter sets data frame with eq
  par.mat$eps = NA       #add column for elim. feas estimator (eps)
  par.mat$eps.sd = NA    #add column for elim. feas estimator st. dev 
  par.mat$pe = NA        #add column for proba elimination (P(e))

#Use array loads from sim runs on savio to fill the data frame
load('Elimination_Feasibility/Savio/fill_array1_500.Rdata')
  inarr = fill.arr
  chunk = 1
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }
  
load('Elimination_Feasibility/Savio/fill_array501_1000.Rdata')
  inarr = fill.arr
  chunk = 2
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }
load('Elimination_Feasibility/Savio/fill_array1001_1500.Rdata')
  inarr = fill.arr
  chunk = 3
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }
load('Elimination_Feasibility/Savio/fill_array1501_2000.Rdata')
  inarr = fill.arr
  chunk = 4
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }
load('Elimination_Feasibility/Savio/fill_array2001_2500.Rdata')
  inarr = fill.arr
  chunk = 5
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }

par.mat2 = par.mat  
#plot P(e) across eps
plot(par.mat2$eps, par.mat2$pe, pch = 18, cex = 0.6, ylim = c(0,1),
     xlab = expression(epsilon), ylab = expression(italic('P(e)')))
  #for(e in 1:nrow(par.mat2)){
  #    segments(y0 = par.mat2$pe[e], y1 = par.mat2$pe[e],
  #             x0 = par.mat2$eps[e] - par.mat2$eps.sd[e],
  #             x1 = par.mat2$eps[e] + par.mat2$eps.sd[e])
  #  } #plot error bars, but makes the plot unreadable

#only plot for vals where p(e) != 1 or 0
par.mat0_1 = subset(par.mat2, pe != 0 & pe != 1 & V2 != 0)

plot(par.mat0_1$eps, par.mat0_1$pe, pch = 18, cex = 0.6, ylim = c(0,1),
     xlab = expression(epsilon), ylab = expression(italic('P(e)')))

#logit transform p(e)
plot(par.mat0_1$eps, log(par.mat0_1$pe / (1 - par.mat0_1$pe)), 
     pch = 18, cex = 0.6, #ylim = c(0,1),
     xlab = expression(epsilon), ylab = expression(italic('P(e)')))

#regression
pe.mod = lm(log(par.mat0_1$pe / (1 - par.mat0_1$pe)) ~ par.mat0_1$eps, data = par.mat0_1)
  summary(pe.mod)
  abline(a = pe.mod$coefficients[1], b = pe.mod$coefficients[2],
         lty = 2, col = 2)
#plot relationship between eps and its standard deviation
plot(par.mat2$eps, par.mat2$eps.sd, pch = 18, cex = 0.6, #ylim = c(0,1),
     ylab = expression(paste(epsilon, ' St. Dev.', sep='')), 
     xlab = expression(epsilon))

#p(e) across transmission intensity
plot(par.mat2$V1, par.mat2$pe, pch = 18, cex = 0.6, ylim = c(0,1),
     xlab = expression(lambda), ylab = expression(italic('P(e)')))

#p(e) across pos. dens. dep.
plot(par.mat2$V2, par.mat2$pe, pch = 18, cex = 0.6, ylim = c(0,1),
     xlab = expression(kappa), ylab = expression(italic('P(e)')))

#surface of p(e) across kappa and lambda
lams = unique(par.mat2$V1)
kaps = unique(par.mat2$V2)

lkmat = matrix(par.mat2$pe, ncol=50, nrow=50)

persp(x = lams, xlim = range(lams), y = kaps, ylim = range(kaps),
      z = lkmat, ticktype = 'detailed', nticks = 4, 
      ylab = 'Pos. Density Dependence',
      xlab = 'Transmission Intensity',
      zlab = 'Probability of Elimination',
      phi = 20, theta = 85, shade = 0.4, col = 'lightblue')

#Same routine for simulations with observation noise #########
par.mat = read.csv('Elimination_Feasibility/eq_vals_for_trans_pars.csv')   #load parameter sets data frame with eq
  par.mat$eps = NA       #add column for elim. feas estimator (eps)
  par.mat$eps.sd = NA    #add column for elim. feas estimator st. dev 
  par.mat$pe = NA        #add column for proba elimination (P(e))

load('Elimination_Feasibility/Savio/wNoise_array1_500.Rdata')
  inarr = fill.arr
  chunk = 1
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }

load('Elimination_Feasibility/Savio/wNoise_array501_1000.Rdata')
inarr = fill.arr
chunk = 2
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }

load('Elimination_Feasibility/Savio/wNoise_array1001_1500.Rdata')
inarr = fill.arr
chunk = 3
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }

load('Elimination_Feasibility/Savio/wNoise_array1501_2000.Rdata')
inarr = fill.arr
chunk = 4
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }

load('Elimination_Feasibility/Savio/wNoise_array2001_2500.Rdata')
inarr = fill.arr
chunk = 5
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    #print(m)
  }

#plot P(e) across eps
plot(par.mat$eps, par.mat$pe, pch = 18, cex = 0.6, ylim = c(0,1),
     xlab = expression(epsilon), ylab = expression(italic('P(e)')))
  