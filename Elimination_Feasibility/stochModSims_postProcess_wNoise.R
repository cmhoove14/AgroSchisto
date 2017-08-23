#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Load datasets from simulations with observation noise #########
library(ggplot2)
library(aod)
par.mat = read.csv('Elimination_Feasibility/eq_vals_for_trans_pars.csv')   #load parameter sets data frame with eq
  par.mat$eps = NA       #add column for elim. feas estimator (eps)
  par.mat$eps.sd = NA    #add column for elim. feas estimator st. dev 
  par.mat$pe = NA        #add column for prob elimination (P(e))
  par.mat$eps.med = NA   #add column for median value of epsilon
  
#Same routine for simulations with observation noise #########
load('Elimination_Feasibility/Savio/Results/wNoise_array1_500.Rdata')
  inarr = fill.arr
  chunk = 1
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    par.mat$eps.med[m+lo-1] = median(inarr[m, 4, ])
    #print(m)
  }

load('Elimination_Feasibility/Savio/Results/wNoise_array501_1000.Rdata')
inarr = fill.arr
chunk = 2
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    par.mat$eps.med[m+lo-1] = median(inarr[m, 4, ])
    #print(m)
  }

load('Elimination_Feasibility/Savio/Results/wNoise_array1001_1500.Rdata')
inarr = fill.arr
chunk = 3
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    par.mat$eps.med[m+lo-1] = median(inarr[m, 4, ])
    #print(m)
  }

load('Elimination_Feasibility/Savio/Results/wNoise_array1501_2000.Rdata')
inarr = fill.arr
chunk = 4
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    par.mat$eps.med[m+lo-1] = median(inarr[m, 4, ])
    #print(m)
  }

load('Elimination_Feasibility/Savio/Results/wNoise_array2001_2500.Rdata')
inarr = fill.arr
chunk = 5
  lo = chunk*dim(inarr)[1] - (dim(inarr)[1]-1)
  hi = chunk*dim(inarr)[1]
  
  for(m in 1:500){
    par.mat$eps[m+lo-1] = mean(inarr[m, 4, ])
    par.mat$eps.sd[m+lo-1] = sd(inarr[m, 4, ])
    par.mat$pe[m+lo-1] = sum(inarr[m, 3, ]) / dim(inarr)[3]  #P(e) = sims that end in elimination / total sims
    par.mat$eps.med[m+lo-1] = median(inarr[m, 4, ])
    #print(m)
  }

colnames(par.mat)[c(1:7)] = c('lambda', 'kappa', 'S', 'E', 'I', 'Wt', 'Wu')  

#Some summaries of P(e) across the transmission parameters ##########
mean(par.mat$pe[par.mat$lambda == min(par.mat$lambda)])
mean(par.mat$pe[par.mat$lambda == max(par.mat$lambda)])

mean(par.mat$pe[par.mat$kappa == min(par.mat$kappa)])
mean(par.mat$pe[par.mat$kappa == max(par.mat$kappa)])

mean(par.mat$pe[round(par.mat$eps, digits = 3) == 0])

lams = unique(par.mat$lambda)
kaps = unique(par.mat$kappa)

lkmat.pe = matrix(par.mat$pe, ncol=50, nrow=50, byrow = T)
            colnames(lkmat.pe) = round(unique(lams), digits = 6)
            rownames(lkmat.pe) = round(unique(kaps), digits = 3)
  
  heat.pe = data.matrix(lkmat.pe)
  
lkmat.eps = matrix(par.mat$eps, ncol=50, nrow=50, byrow = T)
  colnames(lkmat.eps) = round(unique(lams), digits = 6)
  rownames(lkmat.eps) = round(unique(kaps), digits = 3)
  
  heat.eps = data.matrix(lkmat.eps)
  
heatmap(heat.pe, Rowv=NA, Colv=NA, col = heat.colors(256), scale="column")
heatmap(heat.eps, Rowv=NA, Colv=NA, col = heat.colors(256), scale="column")


#plots of and analysis of pe(e) / eps with obs noise #######
plot(par.mat$eps, par.mat$pe, pch = 18, cex = 0.6, ylim = c(0,1),
     xlab = expression(epsilon), ylab = expression(italic('P(e)')))
  points(par.mat$eps[par.mat$kappa == 0], par.mat$pe[par.mat$kappa == 0], 
         pch = 17, cex = 0.8, col = 2)
  legend('bottomleft', legend = c(expression(paste(kappa, ' = 0', sep = '')),
                                expression(paste(kappa, ' > 0', sep = ''))),
         pch = c(17,18), col = c(2,1), bty = 'n', cex = 0.8)
  
#scatterP for publication according to PLos guidelines
tiff("Elimination_Feasibility/plots/PLoS_Figs/Fig4.tiff", height = 11, width = 19.05, units = 'cm', 
      compression = "lzw", res = 300) 

plot(par.mat$eps, par.mat$pe, pch = 18, cex = 0.6, ylim = c(0,1),
     xlab = '', ylab = '', cex.lab = 1.25, mar = c(4.6,6.6,1.1,0.6))
  mtext(side = 2, text = expression(italic('P(e)')), line = 2.4, cex = 1.5)
  mtext(side = 1, text = expression(italic(epsilon)), line = 2.4, cex = 1.5)
  
dev.off()


#st. dev of epsilon
plot(par.mat$eps, par.mat$eps.sd, pch = 18, cex = 0.6, ylim = range(par.mat$eps.sd),# col = 4,
     ylab = expression(paste('st. dev (', epsilon, ')', sep = '')), 
     xlab = expression(epsilon))
  points(par.mat$eps[par.mat$kappa == 0], par.mat$eps.sd[par.mat$kappa == 0], 
         pch = 17, cex = 0.8, col = 2)
  legend('left', legend = c(expression(paste(kappa, ' = 0', sep = '')),
                              expression(paste(kappa, ' > 0', sep = ''))),
         pch = c(17,18), col = c(2,1), bty = 'n', cex = 0.8)
  
#so what about median rather than mean of epsilon vals?
plot(par.mat$eps.med, par.mat$pe, pch = 18, cex = 0.6, ylim = c(0,1),
     xlab = expression(epsilon), ylab = expression(italic('P(e)')))
  points(par.mat$eps.med[par.mat$kappa == 0], par.mat$pe[par.mat$kappa == 0], 
         pch = 17, cex = 0.8, col = 2)
  legend('bottom', legend = c(expression(paste(kappa, ' = 0', sep = '')),
                              expression(paste(kappa, ' > 0', sep = ''))),
         pch = c(17,18), col = c(2,1), bty = 'n', cex = 0.8)
  
#ggplot to show val of lambda / kappa in addition to scatterP
  ggplot(data = par.mat, aes(x = eps, y = pe, col = lambda)) +  
    theme_bw() +
    geom_point(shape = 18, size = 1.2) +
    scale_color_gradient(low = 'wheat1', high = 'red') +
    labs(col = expression(lambda),
         x = expression(epsilon),
         y = expression(italic(P(e))))
  
#only plot for vals where p(e) != 1 or 0
  par.mat01 = subset(par.mat, pe != 0 & pe != 1 & kappa != 0)
  
  plot(par.mat01$eps, par.mat01$pe, pch = 18, cex = 0.6, ylim = c(0,1),
       xlab = expression(epsilon), ylab = '',
       cex.lab = 1.4)
    mtext(side = 2, text = expression(italic('P(e)')), line = 2.4, cex = 1.2)
  
  
#logit transform p(e)
  plot(par.mat01$eps, log(par.mat01$pe / (1 - par.mat01$pe)), 
       pch = 18, cex = 0.6, #ylim = c(0,1),
       xlab = expression(epsilon), ylab = expression(italic('P(e)')))
  
#regression of eps to pe ########
  #with manual logit transform
  pe.mod = lm(log(pe / (1 - pe)) ~ eps, data = par.mat01)
    summary(pe.mod)
  abline(a = pe.mod$coefficients[1], b = pe.mod$coefficients[2],
         lty = 2, col = 2, lwd = 2)
  
  #weighted based on inverse of st. dev of epsilon
  pe.mod.weight = lm(log(pe / (1 - pe)) ~ eps, data = par.mat01, weights = eps.sd^-1)
    summary(pe.mod.weight)
    
  #doesn't change results  
    
#regression of eps to lambda and kappa ########
  plot(par.mat01$lambda, log(par.mat01$pe / (1 - par.mat01$pe)), 
       pch = 18, cex = 0.6, #ylim = c(0,1),
       xlab = expression(lambda), ylab = expression(italic('P(e)')))
  pe.mod.lam = lm(log(pe / (1 - pe)) ~ lambda, data = par.mat01)
    summary(pe.mod.lam)
    abline(a = pe.mod.lam$coefficients[1], b = pe.mod.lam$coefficients[2],
           lty = 2, col = 2, lwd = 2)
    
  plot(par.mat01$kappa, log(par.mat01$pe / (1 - par.mat01$pe)), 
       pch = 18, cex = 0.6, #ylim = c(0,1),
       xlab = expression(kappa), ylab = expression(italic('P(e)')))
  pe.mod.kap = lm(log(pe / (1 - pe)) ~ kappa, data = par.mat01)
    summary(pe.mod.kap)
    abline(a = pe.mod.kap$coefficients[1], b = pe.mod.kap$coefficients[2],
           lty = 2, col = 2, lwd = 2)
    

#plot residuals of logit(p(e)) ~ lambda across epsilon
  plot(par.mat01$eps, resid(pe.mod.lam), pch = 18, cex = 0.6) 
    lm.res.eps = lm(resid(pe.mod.lam) ~ par.mat01$eps)
      summary(lm.res.eps)
    abline(a = lm.res.eps$coefficients[1], b = lm.res.eps$coefficients[2],
           lty = 2, col = 2, lwd = 2)
    
#plot residuals of logit(p(e)) ~ kappa across epsilon and lambda
  plot(par.mat01$eps, resid(pe.mod.kap), pch = 18, cex = 0.6) 
    lm.res.kap = lm(resid(pe.mod.kap) ~ par.mat01$eps)
    summary(lm.res.kap)
    abline(a = lm.res.kap$coefficients[1], b = lm.res.kap$coefficients[2],
           lty = 2, col = 2, lwd = 2)
    
  plot(par.mat01$lambda, resid(pe.mod.kap), pch = 18, cex = 0.6) 
    lm.res.lam = lm(resid(pe.mod.kap) ~ par.mat01$lambda)
    summary(lm.res.lam)
    abline(a = lm.res.lam$coefficients[1], b = lm.res.lam$coefficients[2],
           lty = 2, col = 2, lwd = 2)
    
#definitely cross-correlation between lambda and epsilon... check with plot
  plot(par.mat01$lambda, par.mat01$eps, 
       pch = 18, cex = 0.6, #ylim = c(0,1),
       ylab = expression(epsilon), xlab = expression(lambda))
  eps.lam.mod = lm(eps ~ lambda, data = par.mat01)
      summary(eps.lam.mod)
    abline(a = eps.lam.mod$coefficients[1], b = eps.lam.mod$coefficients[2],
           lty = 2, col = 2, lwd = 2)
    
#Formal nested regression modeling ##########
  mod.lek = lm(log(pe / (1 - pe)) ~ eps + lambda * kappa, data = par.mat01)
  
  anova(mod.lek)
  