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
       ylim = c(-2,0), xlim = c(0, 2000), type = 'l', lwd=2,
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
         col = c(1,2,3,4,1), cex = 0.7, lwd=c(rep(2,4),1), lty=c(rep(1,4),2), bty = 'n')  
  title(expression(paste(pi[C], ' Influence on ', R[0], ' from different sources')), cex = 0.8)

#Rohr 2008 provides the best estimate as it is stable over wide range of concentrations and provides solid data    
  
#Rohr piC, Bakry muN, Baxter/Rohr phiN, and Omran&Salama muN ##########
plot(lowess(conc.atr, conc.atr.means.r0[,4] - r0.He(0)[3], f = 0.01), 
     ylim = c(-4,4), xlim = c(0, 2000), type = 'l', lwd=2,
     xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(atrazine)')))
  for(i in 5:7){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3], f = 0.01), lwd=2, col = i-3)
  }
  for(i in 4:7){
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
  points(r0.atr.fix.n0$atr[r0.atr.fix.n0$study == 'Bakry12'], 
         r0.atr.fix.n0$r0[r0.atr.fix.n0$study == 'Bakry12'] - r0.He(0)[3], pch = 17, col = 6)  
  points(r0.atr.fix.n0$atr[r0.atr.fix.n0$study == 'Omran&Salama'], 
         r0.atr.fix.n0$r0[r0.atr.fix.n0$study == 'Omran&Salama'] - r0.He(0)[3], pch = 17, col = 7)  
  legend('bottomright', title = 'Modeled Functions', 
         legend = c(expression(paste(pi[C], ' Rohr08')), 
                    expression(paste(mu[N], ' Bakry12')),
                    expression(paste(Phi[N], ' Baxter11')),
                    expression(paste(mu[N], ' Omran&Salama')),
                    'St. Dev.'),
         col = c(1,2,3,4,1), cex = 0.7,  lwd=c(rep(2,4),1), lty=c(rep(1,4),2), bty='n') 
  legend('right', title = 'Observed Points', 
         legend = c(expression(paste(pi[C], ' Rohr08')), 
                    expression(paste(mu[N], ' Bakry12')),
                    expression(paste(Phi[N], ' Baxter11'))),
         cex = 0.7, pch = 17, col = c(8,6,4), bty='n') 
  title(expression(paste('Individual atrazine parameters influence on ', R[0])), cex = 0.8)

  
  
#Pairwise combinations (with Rohr08 piC parameter) influence on r0 ##########
plot(lowess(conc.atr, conc.atr.means.r0[,10] - r0.He(0)[3], f = 0.01), 
     ylim = c(-4,4), xlim = c(0, 2000), type = 'l', lwd=2, col = 4,
     xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(atrazine)')))
  for(i in c(11:14)){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3], f = 0.01), lwd=2, col = i-6)
  }
  for(i in c(10:14)){
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] + conc.atr.sds.r0[,i], f = 0.01), col = i-6, lty=2)
    lines(lowess(conc.atr, conc.atr.means.r0[,i] - r0.He(0)[3] - conc.atr.sds.r0[,i], f = 0.01), col = i-6, lty=2)
  }
  legend('bottomleft', title = 'Parameters', 
         legend = c(expression(paste(pi[C], ' & ', Phi[N])), 
                    expression(paste(pi[C], ' & Bakry12 ', mu[N])),
                    expression(paste(Phi[N], ' & Bakry12 ', mu[N])),
                    expression(paste(pi[C], '& Omran&Salama ', mu[N])),
                    expression(paste(Phi[N], '& Omran&Salama ', mu[N])),
                    'St. Dev.'),
         col = c(10:14)-6, cex = 0.7,  lwd=c(rep(2,5),1), lty=c(rep(1,5),2), bty='n') 
  title(expression(paste('Parameter combinations influence on ', R[0])), cex = 0.8)  
  
#Pairwise combinations (with Rohr08 piC parameter) influence on N eq ##########
  plot(lowess(conc.atr, conc.atr.means.Neq[,10] - r0.He(0)[1], f = 0.01), 
       xlim = c(0, 500), ylim = c(-10000, 20000), type = 'l', lwd=2, col = 4,
       xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(atrazine)')))
  for(i in c(11:14)){
    lines(lowess(conc.atr, conc.atr.means.Neq[,i] - r0.He(0)[1], f = 0.01), lwd=2, col = i-6)
  }
  for(i in c(10:14)){
    lines(lowess(conc.atr, conc.atr.means.Neq[,i] - r0.He(0)[1] + conc.atr.sds.Neq[,i], f = 0.01), col = i-6, lty=2)
    lines(lowess(conc.atr, conc.atr.means.Neq[,i] - r0.He(0)[1] - conc.atr.sds.Neq[,i], f = 0.01), col = i-6, lty=2)
  }
  legend('bottomleft', title = 'Parameters', 
         legend = c(expression(paste(pi[C], ' & ', Phi[N])), 
                    expression(paste(pi[C], ' & Bakry12 ', mu[N])),
                    expression(paste(Phi[N], ' & Bakry12 ', mu[N])),
                    expression(paste(pi[C], '& Omran&Salama ', mu[N])),
                    expression(paste(Phi[N], '& Omran&Salama ', mu[N])),
                    'St. Dev.'),
         col = c(10:14)-6, cex = 0.7,  lwd=c(rep(2,5),1), lty=c(rep(1,5),2), bty='n') 
  title(expression(paste('Parameter combinations influence on equilibrium snail density, ', N[eq])), cex = 0.8)  

#Plot total r0 of q function ###############
  #with reduction in snail reproduction
  plot(lowess(conc.atr, conc.atr.means.r0[,17] - r0.He(0)[3], f = 0.01), 
       lwd=2, type = 'l', col = 3, ylim = c(-4,4), xlim = c(0, 2000),
       xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  lines(lowess(conc.atr, 
               (conc.atr.means.r0[,17] + conc.atr.sds.r0[,17] - r0.He(0)[3]), f = 0.01), 
        col = 3, lty=2)
  lines(lowess(conc.atr, 
               (conc.atr.means.r0[,17] - conc.atr.sds.r0[,17] - r0.He(0)[3]), f = 0.01), 
        col = 3, lty=2)
  title('Atrazine - all response functions incorporated')
  
  #Restrict to 200ppb
  plot(lowess(conc.atr, conc.atr.means.r0[,17] - r0.He(0)[3], f = 0.01), 
       lwd=2, type = 'l', col = 3, ylim = c(-4,4), xlim = c(0, 200),
       xlab = 'Atrazine (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  lines(lowess(conc.atr, 
               (conc.atr.means.r0[,17] + conc.atr.sds.r0[,17] - r0.He(0)[3]), f = 0.01), 
        col = 3, lty=2)
  lines(lowess(conc.atr, 
               (conc.atr.means.r0[,17] - conc.atr.sds.r0[,17] - r0.He(0)[3]), f = 0.01), 
        col = 3, lty=2)
  title('Atrazine - all response functions incorporated')
  