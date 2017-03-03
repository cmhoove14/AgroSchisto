#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Run this and skip to plotting functions unless manipulating parameters, simulations, etc. #############
  load('Review_models/r0_of_Atrazine_ws.RData')
  require(rootSolve)

#Compare different piC functions ##########
  plot(lowess(conc.atr, conc.atr.means.r0[,1] - r0.He(0)[3], f = 0.25), 
       ylim = c(-3,0), xlim = c(0, 501), type = 'l', lwd=2,
       xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(atrazine)')))
    for(i in 2:4){
      lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3], f = 0.25), col = i, lwd = 2)
    }
    for(i in 1:4){
      lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] + conc.atr.sds.r0[,i], f = 0.25), 
            col = i, lty=2)
      lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] - conc.atr.sds.r0[,i], f = 0.25), 
            col = i, lty=2)
    }
  legend('bottomleft', title = expression(paste(pi[C], ' Sources')), 
         legend = c('Meta', 'Griggs 2008', 'Koprivnikar 2006', 'Rohr 2008', 'St. Dev.'),
         col = c(1,2,3,4,1), cex = 0.7, lwd=c(rep(2,4),1), lty=c(rep(1,4),2))  
  title(expression(paste(pi[C], ' Influence on ', R[0], ' from different sources')), cex = 0.8)
  
#Rohr piC, Bakry muN, and Baxter/Rohr phiN ##########
plot(lowess(conc.atr, conc.atr.means.r0[,4] - r0.He(0)[3], f = 0.1), col = 4, ylim = c(-3,3), xlim = c(0, 500), type = 'l', lwd=2,
     xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(atrazine)')))
  for(i in 5:6){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3], f = 0.1), lwd=2, col = i)
  }
  for(i in 4:6){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] + conc.atr.sds.r0[,i], f = 0.1), col = i, lty=2)
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] - conc.atr.sds.r0[,i], f = 0.1), col = i, lty=2)
  }
  #Add point estimates from studies
  points(r0.atr.fix.n0$atr[r0.atr.fix.n0$par == 'phi_N'], 
         r0.atr.fix.n0$r0[r0.atr.fix.n0$par == 'phi_N'] - r0.He(0)[3], pch = 15)  
  points(r0.atr.fix.n0$atr[r0.atr.fix.n0$study == 'Rohr08'], 
         r0.atr.fix.n0$r0[r0.atr.fix.n0$study == 'Rohr08'] - r0.He(0)[3], pch = 16)  
  points(r0.atr.fix.n0$atr[r0.atr.fix.n0$par == 'muN'], 
         r0.atr.fix.n0$r0[r0.atr.fix.n0$par == 'muN'] - r0.He(0)[3], pch = 17)  
  legend('bottomleft', title = 'Parameters', 
         legend = c(expression(paste(pi[C], ' Rohr08')), 
                    expression(paste(mu[N], ' Bakry12')),
                    expression(paste(Phi[N], ' Baxter11')),
                    'St. Dev.'),
         col = c(4,5,6,1), cex = 0.7,  lwd=c(rep(2,3),1), lty=c(rep(1,3),2), bty='n') 
  legend(100, -1.75, title = 'Observed Points', 
         legend = c(expression(paste(pi[C], ' Rohr08')), 
                    expression(paste(mu[N], ' Bakry12')),
                    expression(paste(Phi[N], ' Baxter11'))),
         cex = 0.7, pch = c(16,17,15), bty='n') 
  title(expression(paste('Individual atrazine parameters influence on ', R[0])), cex = 0.8)

#Pairwise combinations (with Rohr08 piC parameter) influence on r0 ##########
plot(lowess(conc.atr, conc.atr.means.r0[,9] - r0.He(0)[3], f = 0.1), ylim = c(-3,3), xlim = c(0, 500), type = 'l', lwd=2,
     xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(atrazine)')))
  for(i in c(10,11,13)){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3], f = 0.1), lwd=2, col = i-8)
  }
  for(i in c(9:11,13)){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] + conc.atr.sds.r0[,i], f = 0.1), col = i-8, lty=2)
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] - conc.atr.sds.r0[,i], f = 0.1), col = i-8, lty=2)
  }
  legend('bottomleft', title = 'Parameters', 
         legend = c(expression(paste(pi[C], ' & ', Phi[N])), 
                    expression(paste(pi[C], ' & ', mu[N])),
                    expression(paste(Phi[N], ' & ', mu[N])),
                    expression(paste(pi[C], ', ', mu[N], ' & ', Phi[N])),
                    'St. Dev.'),
         col = c(9:11,13)-8, cex = 0.7,  lwd=c(rep(2,4),1), lty=c(rep(1,4),2), bty='n') 
  title(expression(paste('Parameter combinations influence on ', R[0])), cex = 0.8)  
  
#Pairwise combinations (with Rohr08 piC parameter) influence on N eq ##########
plot(lowess(conc.atr, ((conc.atr.means.Neq[,9] - r0.He(0)[1])/parameters['A']), f = 0.1), ylim = c(-50,100), xlim = c(0, 500), type = 'l', lwd=2,
     xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, 'Snail Density')))
  for(i in c(10,11,13)){
    lines(lowess(conc.atr, (conc.atr.means.Neq[,i] - r0.He(0)[1])/parameters['A'], f = 0.1), lwd=2, col = i-8)
  }
  for(i in c(9:11,13)){
    lines(lowess(conc.atr, (conc.atr.means.Neq[,i] - r0.He(0)[1] + conc.atr.sds.Neq[,i])/parameters['A'], 
                 f = 0.1), col = i-8, lty=2)
    lines(lowess(conc.atr, (conc.atr.means.Neq[,i] - r0.He(0)[1] - conc.atr.sds.Neq[,i])/parameters['A'], 
                 f = 0.1), col = i-8, lty=2)
  }
  legend('bottomleft', title = 'Parameters', 
         legend = c(expression(paste(pi[C], ' & ', Phi[N])), 
                    expression(paste(pi[C], ' & ', mu[N])),
                    expression(paste(Phi[N], ' & ', mu[N])),
                    expression(paste(pi[C], ', ', mu[N], ' & ', Phi[N])),
                    'St. Dev.'),
         col = c(9:11,13)-8, cex = 0.7,  lwd=c(rep(2,4),1), lty=c(rep(1,4),2), bty='n') 
  title(expression(paste('Parameter combinations influence on equilibrium snail density, ', N[eq])), cex = 0.8)  

#individual phiN and muN functions with combined muN/phiN function #########
plot(lowess(conc.atr, conc.atr.means.r0[,5] - r0.He(0)[3], f = 0.1), ylim = c(-3,3), xlim = c(0, 500), 
     type = 'l', lwd=2, col = 5,
     xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(atrazine)')))
  for(i in c(6,11)){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3], f = 0.1), lwd=2, col = i)
  }
  for(i in c(5,6,11)){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] + conc.atr.sds.r0[,i], f = 0.1), col = i, lty=2)
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] - conc.atr.sds.r0[,i], f = 0.1), col = i, lty=2)
  }
  legend('bottomleft', title = 'Parameters', 
         legend = c(expression(mu[N]),
                    expression(Phi[N]), 
                    expression(paste(Phi[N], ' & ', mu[N])),
                    'St. Dev.'),
         col = c(5,6,11,1), cex = 0.7, lwd=c(rep(2,3),1), lty=c(rep(1,3),2), bty='n') 
  title(expression(paste(mu[N],' & ', Phi[N], ' influence on ', R[0])), cex = 0.8)  
#3d plot of direct effect (muN) and bottom up effect (phi_N) influence on R0 #################
  muN.range = seq(min(conc.atr.means.par[, 5]), max(conc.atr.means.par[, 5]), 
                  length.out = 101)
  phiN.range = seq(min(conc.atr.means.par[, 6]), max(conc.atr.means.par[, 6]), 
                   length.out = 101)
  
  m.munphin = matrix(nrow = 101,
                     ncol = 101)
  for(i in 1:length(muN.range)){
    for(j in 1:length(phiN.range)){
      m.munphin[i,j] = r0.fix(mu_Nqx = muN.range[i],
                              phi_Nqx = phiN.range[j])[3]
    }
  }
  
  persp(x = muN.range + parameters['mu_N'], xlim = range(muN.range) + parameters['mu_N'], 
        y = phiN.range * parameters['phi_N']/parameters['A'], ylim = range(phiN.range) * parameters['phi_N']/parameters['A'],
        z = m.munphin - r0.He(0)[3], ticktype = 'detailed', nticks = 4, zlim = c(-3,3),
        xlab = 'Snail mortality rate', ylab = 'Snail carrying capacity', zlab = 'Change in R0',
        phi = 30, theta = 30, shade = 0.4, col = 'lightblue')
  
#loess regression lines fit to single simulations
plot(c(0,max(conc.atr)), c(rep(r0.He(0)[3],2) - r0.He(0)[3]), type='l', lwd = 2, ylim = c(-2.5,2.5), xlim = c(0, max(conc.atr)),
     xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  lines(l1.atr.df$atr, l1.atr.df$pred.2 - r0.He(0)[3], lty = 2, col=2)
    lines(l1.atr.df$atr, l1.atr.df$pred.2 + l1.atr.df$se.2*qnorm(1-.05/2) - r0.He(0)[3], lty = 3, col=2)
    lines(l1.atr.df$atr, l1.atr.df$pred.2 - l1.atr.df$se.2*qnorm(1-.05/2) - r0.He(0)[3], lty = 3, col=2)
#Don't show lines from individual piC experiments to reduce clutter 
  #lines(l2.atr.df$atr, l2.atr.df$pred.2, lty = 2, col=5)
  #  lines(l2.atr.df$atr, l2.atr.df$pred.2 + l2.atr.df$se.2*qnorm(1-.05/2), lty = 3, col=5)
  #  lines(l2.atr.df$atr, l2.atr.df$pred.2 - l2.atr.df$se.2*qnorm(1-.05/2), lty = 3, col=5)
  #lines(l3.atr.df$atr, l3.atr.df$pred.2, lty = 2, col=7)
  #  lines(l3.atr.df$atr, l3.atr.df$pred.2 + l3.atr.df$se.2*qnorm(1-.05/2), lty = 3, col=7)
  #  lines(l3.atr.df$atr, l3.atr.df$pred.2 - l3.atr.df$se.2*qnorm(1-.05/2), lty = 3, col=7)
  lines(l4.atr.df$atr, l4.atr.df$pred.2 - r0.He(0)[3], lty = 2, col=4)
    lines(l4.atr.df$atr, l4.atr.df$pred.2 + l4.atr.df$se.2*qnorm(1-.05/2) - r0.He(0)[3], lty = 3, col=4)
    lines(l4.atr.df$atr, l4.atr.df$pred.2 - l4.atr.df$se.2*qnorm(1-.05/2) - r0.He(0)[3], lty = 3, col=4)
  lines(l5.atr.df$atr, l5.atr.df$pred.2 - r0.He(0)[3], lty = 2, col=6)
    lines(l5.atr.df$atr, l5.atr.df$pred.2 + l5.atr.df$se.2*qnorm(1-.05/2) - r0.He(0)[3], lty = 3, col=6)
    lines(l5.atr.df$atr, l5.atr.df$pred.2 - l5.atr.df$se.2*qnorm(1-.05/2) - r0.He(0)[3], lty = 3, col=6)
  lines(l6.atr.df$atr, l6.atr.df$pred.2 - r0.He(0)[3], lty = 2, col=3, lwd=2)
    lines(l6.atr.df$atr, l6.atr.df$pred.2 + l6.atr.df$se.2*qnorm(1-.05/2) - r0.He(0)[3], lty = 3, col=3, lwd=2)
    lines(l6.atr.df$atr, l6.atr.df$pred.2 - l6.atr.df$se.2*qnorm(1-.05/2) - r0.He(0)[3], lty = 3, col=3, lwd=2)
#Point estimates from studies
  points(r0.atr.fix.n0$atr[r0.atr.fix.n0$par == 'phi_N'], r0.atr.fix.n0$r0[r0.atr.fix.n0$par == 'phi_N'] - r0.He(0)[3],
         pch = 17, cex = 0.8, col = 6)  
  points(r0.atr.fix.n0$atr[r0.atr.fix.n0$par == 'piC'], r0.atr.fix.n0$r0[r0.atr.fix.n0$par == 'piC'] - r0.He(0)[3],
         pch = 17, cex = 0.8, col = 2)  
  points(r0.atr.fix.n0$atr[r0.atr.fix.n0$par == 'muN'], r0.atr.fix.n0$r0[r0.atr.fix.n0$par == 'muN'] - r0.He(0)[3],
         pch = 17, cex = 0.8, col = 4)  
  legend('bottomleft', title = 'Parameters', legend = c('Baseline',
                                                       expression(paste(pi[C], ' meta')),
                                                       expression(paste(mu[N], ' Bakry 2012')),
                                                       expression(paste(Phi[N], ' Baxter 2011')),
                                                       'Combined functions'),
         lty = c(1, rep(2,4)), col = c(1,2,4,6,3), cex = 0.8)  
  title(expression(paste('Atrazine influence on ', R[0], sep ='')))
    
#loess regression lines fit to all simulations
plot(c(0,max(conc.atr)), c(rep(r0.He(0)[3],2)), type='l', lwd = 2, ylim = c(0,7), xlim = c(0, max(conc.atr)),
     xlab = 'Atrazine (ppb)', ylab = expression(paste(R[0], '(q)')))
  lines(l1.atr.df$atr, l1.atr.df$pred, lty = 2, col=2)
    lines(l1.atr.df$atr, l1.atr.df$pred + l1.atr.df$se*qnorm(1-.05/2), lty = 3, col=2)
    lines(l1.atr.df$atr, l1.atr.df$pred - l1.atr.df$se*qnorm(1-.05/2), lty = 3, col=2)
  lines(l2.atr.df$atr, l2.atr.df$pred, lty = 2, col=3)
    lines(l2.atr.df$atr, l2.atr.df$pred + l2.atr.df$se*qnorm(1-.05/2), lty = 3, col=3)
    lines(l2.atr.df$atr, l2.atr.df$pred - l2.atr.df$se*qnorm(1-.05/2), lty = 3, col=3)
  lines(l3.atr.df$atr, l3.atr.df$pred, lty = 2, col=4)
    lines(l3.atr.df$atr, l3.atr.df$pred + l1.atr.df$se*qnorm(1-.05/2), lty = 3, col=4)
    lines(l3.atr.df$atr, l3.atr.df$pred - l1.atr.df$se*qnorm(1-.05/2), lty = 3, col=4)
  lines(l4.atr.df$atr, l4.atr.df$pred, lty = 2, col=5)
    lines(l4.atr.df$atr, l4.atr.df$pred + l1.atr.df$se*qnorm(1-.05/2), lty = 3, col=5)
    lines(l4.atr.df$atr, l4.atr.df$pred - l1.atr.df$se*qnorm(1-.05/2), lty = 3, col=5)
  lines(l5.atr.df$atr, l5.atr.df$pred, lty = 2, col=6)
    lines(l5.atr.df$atr, l5.atr.df$pred + l1.atr.df$se*qnorm(1-.05/2), lty = 3, col=6)
    lines(l5.atr.df$atr, l5.atr.df$pred - l1.atr.df$se*qnorm(1-.05/2), lty = 3, col=6)
  lines(l6.atr.df$atr, l6.atr.df$pred, lty = 2, col=7)
    lines(l6.atr.df$atr, l6.atr.df$pred + l1.atr.df$se*qnorm(1-.05/2), lty = 3, col=7)
    lines(l6.atr.df$atr, l6.atr.df$pred - l1.atr.df$se*qnorm(1-.05/2), lty = 3, col=7)