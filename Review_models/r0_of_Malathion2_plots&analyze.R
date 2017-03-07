#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

  load('Review_models/r0_of_malathion2_ws.RData') #conc range = 0 - 2000 by 2s, 2000 - 2e5 by 200s
  require(rootSolve)

#Plot results of individual parameter functions ###################
#mean of all r0 sims   
plot(lowess(conc.mal, conc.mal.means.r0[,1] - r0.In(0)[3], f = 0.01), 
     lwd=2, type = 'l', ylim = c(-3,3), xlim = c(0, 150000),
     xlab = 'Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  for(i in 2:length(parfx)){
    lines(lowess(conc.mal, conc.mal.means.r0[,i] - r0.In(0)[3], f = 0.01), lwd=2, col = i)
  }
  for(i in 1:length(parfx)){
    lines(lowess(conc.mal, (conc.mal.means.r0[,i] + conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i, lty=2)
    lines(lowess(conc.mal, (conc.mal.means.r0[,i] - conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i, lty=2)
  }
  legend(135000, 2, lwd=2, col = c(1:5), cex = 0.7, bty = 'n',
         legend = c(expression(paste(pi[C], sep = '')), 
                    expression(paste(pi[M], sep = '')), 
                    expression(paste('f'[N], sep = '')), 
                    expression(paste(mu[P], sep = '')), 
                    expression(paste(mu[N], sep = ''))))
  
#zoom   
plot(lowess(conc.mal, conc.mal.means.r0[,1] - r0.In(0)[3], f = 0.01), 
     lwd=2, type = 'l', ylim = c(-4,0), xlim = c(0, 2000),
     xlab = 'Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  for(i in 3:4){
    lines(lowess(conc.mal[c(0:1000)], (conc.mal.means.r0[c(0:1000),i] - r0.In(0)[3]), f = 0.01), 
          lwd=2, col = i)
  }
  for(i in 3:4){
    lines(lowess(conc.mal[c(0:1000)], 
                 (conc.mal.means.r0[c(0:1000),i] + conc.mal.sds.r0[c(0:1000),i]) - r0.In(0)[3], f = 0.01), 
          col = i, lty=2)
    lines(lowess(conc.mal[c(0:1000)], 
                 (conc.mal.means.r0[c(0:1000),i] - conc.mal.sds.r0[c(0:1000),i]) - r0.In(0)[3], f = 0.01), 
          col = i, lty=2)
  }
    legend('right', lwd=2, col = c(1,3,4), cex = 0.7, bty = 'n',
           legend = c(expression(paste(pi[C], sep = '')), 
                      expression(paste('f'[N], sep = '')), 
                      expression(paste(mu[N], sep = ''))))
  
#log scale   
plot(lowess(log(conc.mal+1), conc.mal.means.r0[,1] - r0.In(0)[3], f = 0.01), 
     lwd=2, type = 'l', ylim = c(-4,4), xlim = log(c(0, 200000)+1),
     xlab = 'ln+1 Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  for(i in 2:length(parfx)){
    lines(lowess(log(conc.mal+1), conc.mal.means.r0[,i] - r0.In(0)[3], f = 0.01), lwd=2, col = i)
  }
  for(i in 1:length(parfx)){
    lines(lowess(log(conc.mal+1), 
                 (conc.mal.means.r0[,i] + conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i, lty=2)
    lines(lowess(log(conc.mal+1), 
                 (conc.mal.means.r0[,i] - conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i, lty=2)
  }
legend('topleft', lwd=2, col = c(1:5), cex = 0.7, bty = 'n',
       legend = c(expression(paste(pi[C], sep = '')), 
                  expression(paste(pi[M], sep = '')), 
                  expression(paste('f'[N], sep = '')), 
                  expression(paste(mu[N], sep = '')), 
                  expression(paste(mu[P], sep = ''))))

#Add Point estimates from studies
  points(log(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'muN']+1), 
         r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'muN'] - r0.In(0)[3],
         pch = 17, cex = 0.8, col = 6)  
  points(log(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'piC']+1), 
         r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'piC'] - r0.In(0)[3],
         pch = 17, cex = 0.8, col = 'gold')  
  points(log(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'piM']+1), 
         r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'piM'] - r0.In(0)[3],
         pch = 17, cex = 0.8, col = 8)  
  points(log(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'muP']+1), 
         r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'muP'] - r0.In(0)[3],
         pch = 17, cex = 0.8, col = 'orange')
  legend('bottomleft', title = 'Observed Points', 
         legend = c(expression(paste(mu[N], ' Bakry 2011')),
                    expression(paste(pi[C], ' Tchounwou 1992')),
                    expression(paste(pi[M], ' Tchounwou 1991')),
                    expression(paste(mu[P], ' Halstead 2015'))),
         pch = 17, col = c(6,'gold', 8, 'orange'), cex = 0.7, bty = 'n')  
  title('Malathion')
  
#Plot results of parameter combinations #########
#Larval toxicity  
  plot(lowess(log(conc.mal+1), conc.mal.means.r0[,1] - r0.In(0)[3], f = 0.01), 
       lwd=2, type = 'l', ylim = c(-4,0), xlim = log(c(0, 200000)+1),
       xlab = 'ln+1 Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  for(i in c(2,6)){
    lines(lowess(log(conc.mal+1), conc.mal.means.r0[,i] - r0.In(0)[3], f = 0.01), lwd=2, col = i)
  }
  for(i in c(1,2,6)){
    lines(lowess(log(conc.mal+1), 
                 (conc.mal.means.r0[,i] + conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i, lty=2)
    lines(lowess(log(conc.mal+1), 
                 (conc.mal.means.r0[,i] - conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i, lty=2)
  }
  legend('bottomleft', lwd=2, col = c(1,2,6), cex = 0.7, bty = 'n',
         legend = c(expression(paste(pi[C], sep = '')), 
                    expression(paste(pi[M], sep = '')), 
                    expression(paste(pi[C], ' & ', pi[M], sep = ' '))))
  title('Malathion larval effects')
  
#Snail toxicity
  plot(lowess(log(conc.mal+1), conc.mal.means.r0[,3] - r0.In(0)[3], f = 0.01), 
       lwd=2, type = 'l', col = 3, ylim = c(-4,0), xlim = log(c(0, 200000)+1),
       xlab = 'ln+1 Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  for(i in c(4,7)){
    lines(lowess(log(conc.mal+1), conc.mal.means.r0[,i] - r0.In(0)[3], f = 0.01), lwd=2, col = i)
  }
  for(i in c(3,4,7)){
    lines(lowess(log(conc.mal+1), 
                 (conc.mal.means.r0[,i] + conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i, lty=2)
    lines(lowess(log(conc.mal+1), 
                 (conc.mal.means.r0[,i] - conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i, lty=2)
  }
  legend('bottomleft', lwd=2, col = c(3,4,7), cex = 0.7, bty = 'n',
         legend = c(expression(paste('f'[N], sep = '')), 
                    expression(paste(mu[N], sep = '')), 
                    expression(paste('f'[N], ' & ', mu[N], sep = ''))))
  title('Malathion snail effects')
  
#Predator mortality with larval and snail effects
  plot(lowess(log(conc.mal+1), conc.mal.means.r0[,5] - r0.In(0)[3], f = 0.01), 
       lwd=2, type = 'l', col = 5, ylim = c(-4,4), xlim = log(c(0, 200000)+1),
       xlab = 'ln+1 Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  lines(lowess(log(conc.mal+1), 
               (conc.mal.means.r0[,5] + conc.mal.sds.r0[,5] - r0.In(0)[3]), f = 0.01), 
        col = 5, lty=2)
  lines(lowess(log(conc.mal+1), 
               (conc.mal.means.r0[,5] - conc.mal.sds.r0[,5] - r0.In(0)[3]), f = 0.01), 
        col = 5, lty=2)
  for(i in c(8,9)){
    lines(lowess(log(conc.mal+1), conc.mal.means.r0[,i] - r0.In(0)[3], f = 0.01), lwd=2, col = i-2)
  }
  for(i in c(8,9)){
    lines(lowess(log(conc.mal+1), 
                 (conc.mal.means.r0[,i] + conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i-2, lty=2)
    lines(lowess(log(conc.mal+1), 
                 (conc.mal.means.r0[,i] - conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i-2, lty=2)
  }
  legend('topleft', lwd=2, col = c(5,6,7), cex = 0.7, bty = 'n',
         legend = c(expression(paste(mu[P], sep = '')),
                    expression(paste(mu[P], ' & ', pi[C], ' & ', pi[C])),
                    expression(paste(mu[P], ' & f'[N], ' & ', mu[N]))))
  
#Plot total r0 of q function ###############
  #with reduction in snail reproduction
  plot(lowess(log(conc.mal+1), conc.mal.means.r0[,10] - r0.In(0)[3], f = 0.01), 
       lwd=2, type = 'l', col = 2, ylim = c(-4,4), xlim = log(c(0, 200000)+1),
       xlab = 'ln+1 Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  lines(lowess(log(conc.mal+1), 
               (conc.mal.means.r0[,10] + conc.mal.sds.r0[,10] - r0.In(0)[3]), f = 0.01), 
        col = 2, lty=2)
  lines(lowess(log(conc.mal+1), 
               (conc.mal.means.r0[,10] - conc.mal.sds.r0[,10] - r0.In(0)[3]), f = 0.01), 
        col = 2, lty=2)
  title('Malathion - all possible response functions')
  
  #W/o reduction in snail reproduction 
  #with reduction in snail reproduction
  plot(lowess(log(conc.mal+1), conc.mal.means.r0[,11] - r0.In(0)[3], f = 0.01), 
       lwd=2, type = 'l', col = 2, ylim = c(-4,4), xlim = log(c(0, 200000)+1),
       xlab = 'ln+1 Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  lines(lowess(log(conc.mal+1), 
               (conc.mal.means.r0[,11] + conc.mal.sds.r0[,11] - r0.In(0)[3]), f = 0.01), 
        col = 2, lty=2)
  lines(lowess(log(conc.mal+1), 
               (conc.mal.means.r0[,11] - conc.mal.sds.r0[,11] - r0.In(0)[3]), f = 0.01), 
        col = 2, lty=2)
  title('Malathion - no reduction in snail reproduction')
  