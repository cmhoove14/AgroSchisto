#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

  load('Review_models/r0_of_Atrazine_ws.RData')
  require(rootSolve)

#Compare different piC functions ##########
  plot(lowess(conc.atr, conc.atr.means.r0[,1] - r0.He(0)[3], f = 0.25), 
       ylim = c(-2,0), xlim = c(0, 300), type = 'l', lwd=2,
       xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(atrazine)')))
    for(i in 2:4){
      lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3], f = 0.25), col = i, lwd = 2)
    }
    for(i in 1:4){
      lines(lowess(conc.atr, (conc.atr.means.r0[,i] + conc.atr.sds.r0[,i]) - r0.He(0)[3], f = 0.25), 
            col = i, lty=2)
      lines(lowess(conc.atr, (conc.atr.means.r0[,i] - conc.atr.sds.r0[,i]) - r0.He(0)[3], f = 0.25), 
            col = i, lty=2)
    }
  legend('bottomleft', title = expression(paste(pi[C], ' Sources')), 
         legend = c('Meta', 'Griggs 2008', 'Koprivnikar 2006', 'Rohr 2008', 'St. Dev.'),
         col = c(1,2,3,4,1), cex = 0.7, lwd=c(rep(2,4),1), lty=c(rep(1,4),2))  
  title(expression(paste(pi[C], ' Influence on ', R[0], ' from different sources')), cex = 0.8)
  
#Rohr piC, Bakry muN, and Baxter/Rohr phiN ##########
plot(lowess(conc.atr, conc.atr.means.r0[,4] - r0.He(0)[3], f = 0.01), 
     ylim = c(-3,3), xlim = c(0, 500), type = 'l', lwd=2,
     xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(atrazine)')))
  for(i in 5:6){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3], f = 0.01), lwd=2, col = i-3)
  }
  for(i in 4:6){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] + conc.atr.sds.r0[,i], f = 0.01), 
          col = i-3, lty=2)
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] - conc.atr.sds.r0[,i], f = 0.01), 
          col = i-3, lty=2)
  }
  #Add point estimates from studies
  points(r0.atr.fix.n0$atr[r0.atr.fix.n0$par == 'phi_N'], 
         r0.atr.fix.n0$r0[r0.atr.fix.n0$par == 'phi_N'] - r0.He(0)[3], pch = 17, col = 4)  
  points(r0.atr.fix.n0$atr[r0.atr.fix.n0$study == 'Rohr08'], 
         r0.atr.fix.n0$r0[r0.atr.fix.n0$study == 'Rohr08'] - r0.He(0)[3], pch = 17, col = 8)  
  points(r0.atr.fix.n0$atr[r0.atr.fix.n0$par == 'muN'], 
         r0.atr.fix.n0$r0[r0.atr.fix.n0$par == 'muN'] - r0.He(0)[3], pch = 17, col = 6)  
  legend('bottomleft', title = 'Parameters', 
         legend = c(expression(paste(pi[C], ' Rohr08')), 
                    expression(paste(mu[N], ' Bakry12')),
                    expression(paste(Phi[N], ' Baxter11')),
                    'St. Dev.'),
         col = c(1,2,3,1), cex = 0.7,  lwd=c(rep(2,3),1), lty=c(rep(1,3),2), bty='n') 
  legend(100, -1.75, title = 'Observed Points', 
         legend = c(expression(paste(pi[C], ' Rohr08')), 
                    expression(paste(mu[N], ' Bakry12')),
                    expression(paste(Phi[N], ' Baxter11'))),
         cex = 0.7, pch = 17, col = c(8,6,4), bty='n') 
  title(expression(paste('Individual atrazine parameters influence on ', R[0])), cex = 0.8)

#Pairwise combinations (with Rohr08 piC parameter) influence on r0 ##########
plot(lowess(conc.atr, conc.atr.means.r0[,9] - r0.He(0)[3], f = 0.01), 
     ylim = c(-3,3), xlim = c(0, 500), type = 'l', lwd=2, col = 3,
     xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(atrazine)')))
  for(i in c(10,11,13)){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3], f = 0.01), lwd=2, col = i-6)
  }
  for(i in c(9:11,13)){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] + conc.atr.sds.r0[,i], f = 0.01), col = i-6, lty=2)
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] - conc.atr.sds.r0[,i], f = 0.01), col = i-6, lty=2)
  }
  legend('bottomleft', title = 'Parameters', 
         legend = c(expression(paste(pi[C], ' & ', Phi[N])), 
                    expression(paste(pi[C], ' & ', mu[N])),
                    expression(paste(Phi[N], ' & ', mu[N])),
                    expression(paste(pi[C], ', ', mu[N], ' & ', Phi[N])),
                    'St. Dev.'),
         col = c(9:11,13)-6, cex = 0.7,  lwd=c(rep(2,4),1), lty=c(rep(1,4),2), bty='n') 
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
plot(lowess(conc.atr, conc.atr.means.r0[,5] - r0.He(0)[3], f = 0.01), ylim = c(-3,3), xlim = c(0, 500), 
     type = 'l', lwd=2, col = 5,
     xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(atrazine)')))
  for(i in c(6,11)){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3], f = 0.01), lwd=2, col = i)
  }
  for(i in c(5,6,11)){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] + conc.atr.sds.r0[,i], f = 0.01), col = i, lty=2)
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] - conc.atr.sds.r0[,i], f = 0.01), col = i, lty=2)
  }
  legend('bottomleft', title = 'Parameters', 
         legend = c(expression(mu[N]),
                    expression(Phi[N]), 
                    expression(paste(Phi[N], ' & ', mu[N])),
                    'St. Dev.'),
         col = c(5,6,11,1), cex = 0.7, lwd=c(rep(2,3),1), lty=c(rep(1,3),2), bty='n') 
  title(expression(paste(mu[N],' & ', Phi[N], ' influence on ', R[0])), cex = 0.8)  