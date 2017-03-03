#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Load results of simulations #############
  load('Review_models/r0_of_malathion2_ws.RData')
  require(rootSolve)

#Plot results ###################
#mean of all sims as points  
plot(lowess(conc.mal, conc.mal.means.r0[,1] - r0.In(0)[3], f = 0.25), lty=2, lwd=2, type = 'l', col = 2, 
     xlab = 'Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')),
     ylim = c(-3,3), xlim = c(0, max(conc.mal)))
  for(i in 2:(length(parfx)+2)){
    lines(lowess(conc.mal, conc.mal.means.r0[,i] - r0.In(0)[3], f = 0.25), lty=2, lwd=2, col = i+1)
  }
  legend('topright', lty=2, col = c(2,3,6,7,8), cex = 0.7,
         legend = c('piC', 'piM', 'muP', 'muN', 'piC & piM & muP'))
  
plot(conc.mal, conc.mal.means.r0[,1], pch = 17, cex=0.5, col = 2, xlab = 'Malathion (ppb)', ylab = expression(paste(R[0], '(q)')),
     ylim = c(0,4), xlim = c(0, 1000))
  for(i in 2:(length(parfx)+2)){
    points(conc.mal, conc.mal.means.r0[,i], pch = 17, cex=0.5, col = i+1)
  }  
  legend('bottomleft', pch = 17, col = c(4,5,7), cex = 0.7,
         legend = c('fN', 'muN', 'fN & muN'))

#loess regression lines of single sims
plot(c(0,max(conc.mal)), c(rep(r0.In(0)[3],2) - r0.In(0)[3]), type='l', lwd = 2, ylim = c(-4,3), xlim = c(0, max(conc.mal)),
     xlab = 'Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  lines(l1.mal.df$mal, l1.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=2)
    lines(l1.mal.df$mal, l1.mal.df$pred.2 + l1.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=2)
    lines(l1.mal.df$mal, l1.mal.df$pred.2 - l1.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=2)
  lines(l2.mal.df$mal, l2.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=5)
    lines(l2.mal.df$mal, l2.mal.df$pred.2 + l2.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=5)
    lines(l2.mal.df$mal, l2.mal.df$pred.2 - l2.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=5)
  lines(l3.mal.df$mal, l3.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=7)
    lines(l3.mal.df$mal, l3.mal.df$pred.2 + l3.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=7)
    lines(l3.mal.df$mal, l3.mal.df$pred.2 - l3.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=7)
  lines(l4.mal.df$mal, l4.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=4)
    lines(l4.mal.df$mal, l4.mal.df$pred.2 + l4.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=4)
    lines(l4.mal.df$mal, l4.mal.df$pred.2 - l4.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=4)
  lines(l5.mal.df$mal, l5.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=6)
    lines(l5.mal.df$mal, l5.mal.df$pred.2 + l5.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=6)
    lines(l5.mal.df$mal, l5.mal.df$pred.2 - l5.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=6)
  lines(l6.mal.df$mal, l6.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=8)
    lines(l6.mal.df$mal, l6.mal.df$pred.2 + l6.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=8)
    lines(l6.mal.df$mal, l6.mal.df$pred.2 - l6.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=8)
  lines(l7.mal.df$mal, l7.mal.df$pred.2 - r0.In(0)[3], lty = 2, col=3)
    lines(l7.mal.df$mal, l7.mal.df$pred.2 + l7.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=3)
    lines(l7.mal.df$mal, l7.mal.df$pred.2 - l7.mal.df$se.2*qnorm(1-.05/2) - r0.In(0)[3], lty = 3, col=3)
    
      #Point estimates from studies
  points(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'phi_N'], r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'phi_N'] - r0.In(0)[3],
         pch = 17, cex = 0.8, col = 6)  
  points(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'piC'], r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'piC'] - r0.In(0)[3],
         pch = 17, cex = 0.8, col = 2)  
  points(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'muN'], r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'muN'] - r0.In(0)[3],
         pch = 17, cex = 0.8, col = 4)  
  legend('bottomleft', title = 'Parameters', legend = c('Baseline',
                                                        expression(paste(pi[C], ' meta')),
                                                        expression(paste(mu[N], ' Bakry 2012')),
                                                        expression(paste(Phi[N], ' Baxter 2011')),
                                                        'Combined functions'),
         lty = c(1, rep(2,4)), col = c(1,2,4,6,3), cex = 0.8)  
  title(expression(paste('Malathion influence on ', R[0], sep ='')))
  

  legend('topright', legend = c(expression(paste(pi[C], '(q)')),
                                expression(paste(pi[M], '(q)')),
                                expression(paste(f[N], '(q)')),
                                expression(paste(mu[N], '(q)')),
                                expression(paste(mu[P], '(q)'))), pch = 17, col = c(2:6), cex = 0.7)