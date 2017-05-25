#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

load('Review_models/Savio/Malathion/r0_malathion_savio.Rdata')

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
  legend(100000, 2.4, lwd=2, col = c(1:7), cex = 0.6, bty = 'n',
         legend = c(expression(paste(pi[C], ' - Tchounwou 1992', sep = '')), 
                    expression(paste(pi[M], ' - Tchounwou 1991a', sep = '')), 
                    expression(paste('f'[N], ' - Bakry 2011', sep = '')), 
                    expression(paste(mu[N], ' - Bakry 2011', sep = '')),
                    expression(paste(mu[P], ' - Halstead 2015', sep = '')), 
                    expression(paste('f'[N], ' - Tchounwou 1991b', sep = '')), 
                    expression(paste(mu[N], ' - Tchounwou 1991b', sep = ''))))
  
#zoom   
plot(lowess(conc.mal, conc.mal.means.r0[,1] - r0.In(0)[3], f = 0.01), 
     lwd=2, type = 'l', ylim = c(-4,0), xlim = c(0, 1000),
     xlab = 'Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  for(i in 2:length(parfx)){
    lines(lowess(conc.mal[c(0:1000)], (conc.mal.means.r0[c(0:1000),i] - r0.In(0)[3]), f = 0.01), 
          lwd=2, col = i)
  }
  for(i in 1:length(parfx)){
    lines(lowess(conc.mal[c(0:1000)], 
                 (conc.mal.means.r0[c(0:1000),i] + conc.mal.sds.r0[c(0:1000),i]) - r0.In(0)[3], f = 0.01), 
          col = i, lty=2)
    lines(lowess(conc.mal[c(0:1000)], 
                 (conc.mal.means.r0[c(0:1000),i] - conc.mal.sds.r0[c(0:1000),i]) - r0.In(0)[3], f = 0.01), 
          col = i, lty=2)
  }
legend('bottomleft', lwd=2, col = c(1:7), cex = 0.6, bty = 'n',
       legend = c(expression(paste(pi[C], ' - Tchounwou 1992', sep = '')), 
                  expression(paste(pi[M], ' - Tchounwou 1991a', sep = '')), 
                  expression(paste('f'[N], ' - Bakry 2011', sep = '')), 
                  expression(paste(mu[N], ' - Bakry 2011', sep = '')),
                  expression(paste(mu[P], ' - Halstead 2015', sep = '')), 
                  expression(paste('f'[N], ' - Tchounwou 1991b', sep = '')), 
                  expression(paste(mu[N], ' - Tchounwou 1991b', sep = ''))))
  
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
legend('topleft', lwd=2, col = c(1:7), cex = 0.6, bty = 'n',
       legend = c(expression(paste(pi[C], ' - Tchounwou 1992', sep = '')), 
                  expression(paste(pi[M], ' - Tchounwou 1991a', sep = '')), 
                  expression(paste('f'[N], ' - Bakry 2011', sep = '')), 
                  expression(paste(mu[N], ' - Bakry 2011', sep = '')),
                  expression(paste(mu[P], ' - Halstead 2015', sep = '')), 
                  expression(paste('f'[N], ' - Tchounwou 1991b', sep = '')), 
                  expression(paste(mu[N], ' - Tchounwou 1991b', sep = ''))))

#Get rid of Bakry 2011 functions (only two reliable data points on non-schisto snail species) ##########

conc.mal.means.r02 = conc.mal.means.r0[, -c(3,4)]
conc.mal.sds.r02 = conc.mal.sds.r0[, -c(3,4)]

#mean of r0 sims w/o Bakry study
plot(lowess(conc.mal, conc.mal.means.r02[,1] - r0.In(0)[3], f = 0.01), 
     lwd=2, type = 'l', ylim = c(-3,3), xlim = c(0, 150000),
     xlab = 'Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
for(i in 2:(length(parfx)-2)){
  lines(lowess(conc.mal, conc.mal.means.r02[,i] - r0.In(0)[3], f = 0.01), lwd=2, col = i)
}
for(i in 1:(length(parfx)-2)){
  lines(lowess(conc.mal, (conc.mal.means.r02[,i] + conc.mal.sds.r02[,i] - r0.In(0)[3]), f = 0.01), 
        col = i, lty=2)
  lines(lowess(conc.mal, (conc.mal.means.r02[,i] - conc.mal.sds.r02[,i] - r0.In(0)[3]), f = 0.01), 
        col = i, lty=2)
}
legend('bottomleft', lwd=2, col = c(1:5), cex = 0.5, bty = 'n', title = 'Modeled Functions',
       legend = c(expression(paste(pi[C], ' - Tchounwou 1992', sep = '')), 
                  expression(paste(pi[M], ' - Tchounwou 1991a', sep = '')), 
                  expression(paste(mu[P], ' - Halstead 2015', sep = '')), 
                  expression(paste(mu[N], ' - Tchounwou 1991b', sep = '')), 
                  expression(paste('f'[N], ' - Tchounwou 1991b', sep = ''))))

#Add Point estimates from studies
  points(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'muN' & r0.mal.fix.n0$study == 'Tch91b'], 
         r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'muN' & r0.mal.fix.n0$study == 'Tch91b'] - r0.In(0)[3],
         pch = 17, col = 6)  
  points(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'fN' & r0.mal.fix.n0$study == 'Tch91b'], 
         r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'fN' & r0.mal.fix.n0$study == 'Tch91b'] - r0.In(0)[3],
         pch = 17, col = 2)
  points(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'piC'], 
         r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'piC'] - r0.In(0)[3],
         pch = 17, col = 'gold')  
  points(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'piM'], 
         r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'piM'] - r0.In(0)[3],
         pch = 17)  
  points(r0.mal.fix.n0$mal[r0.mal.fix.n0$par == 'muP'], 
         r0.mal.fix.n0$r0[r0.mal.fix.n0$par == 'muP'] - r0.In(0)[3],
         pch = 17, col = 'orange')
  legend(27500, -1.8, title = 'Observed Points', 
         legend = c(expression(paste(pi[C], ' - Tchounwou 1992', sep = '')), 
                    expression(paste(pi[M], ' - Tchounwou 1991a', sep = '')), 
                    expression(paste(mu[P], ' - Halstead 2015', sep = '')), 
                    expression(paste(mu[N], ' - Tchounwou 1991b', sep = '')), 
                    expression(paste('f'[N], ' - Tchounwou 1991b', sep = ''))),
         pch = 17, col = c('gold', 1, 'orange', 6, 2), cex = 0.5, bty = 'n')  
  title('Malathion component effects')
  
#Plot results of parameter combinations #########
#Larval toxicity  
  plot(lowess(conc.mal, conc.mal.means.r0[,1] - r0.In(0)[3], f = 0.01), 
       lwd=2, type = 'l', ylim = c(-4,0), xlim = c(0, 150000),
       xlab = 'Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  for(i in c(2,8)){
    lines(lowess(conc.mal, conc.mal.means.r0[,i] - r0.In(0)[3], f = 0.01), lwd=2, col = i)
  }
  for(i in c(1,2,8)){
    lines(lowess(conc.mal, 
                 (conc.mal.means.r0[,i] + conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i, lty=2)
    lines(lowess(conc.mal, 
                 (conc.mal.means.r0[,i] - conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i, lty=2)
  }
  legend('bottomleft', lwd=2, col = c(1,2,8), cex = 0.7, bty = 'n',
         legend = c(expression(paste(pi[C], sep = '')), 
                    expression(paste(pi[M], sep = '')), 
                    expression(paste(pi[C], ' & ', pi[M], sep = ' '))))
  title('Malathion combined larval effects')
  
#Snail toxicity
  plot(lowess(conc.mal, conc.mal.means.r0[,6] - r0.In(0)[3], f = 0.01), 
       lwd=2, type = 'l', col = 4, ylim = c(-4,0), xlim = c(0, 150000),
       xlab = 'Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  for(i in c(7,10)){
    lines(lowess(conc.mal, conc.mal.means.r0[,i] - r0.In(0)[3], f = 0.01), lwd=2, col = i-2)
  }
  for(i in c(6,7,10)){
    lines(lowess(conc.mal, 
                 (conc.mal.means.r0[,i] + conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i-2, lty=2)
    lines(lowess(conc.mal, 
                 (conc.mal.means.r0[,i] - conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i-2, lty=2)
  }
  legend('bottomleft', lwd=2, col = c(4,5,8), cex = 0.7, bty = 'n',
         legend = c(expression(paste(mu[N], sep = '')), 
                    expression(paste('f'[N], sep = '')), 
                    expression(paste('f'[N], ' & ', mu[N], sep = ''))))
  title('Malathion combined snail effects')
  
#Predator mortality with larval and snail effects
  plot(lowess(conc.mal, conc.mal.means.r0[,5] - r0.In(0)[3], f = 0.01), 
       lwd=2, type = 'l', col = 3, ylim = c(-4,4), xlim = c(0, 150000),
       xlab = 'Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  lines(lowess(conc.mal, 
               (conc.mal.means.r0[,5] + conc.mal.sds.r0[,5] - r0.In(0)[3]), f = 0.01), 
        col = 3, lty=2)
  lines(lowess(conc.mal, 
               (conc.mal.means.r0[,5] - conc.mal.sds.r0[,5] - r0.In(0)[3]), f = 0.01), 
        col = 3, lty=2)
  for(i in c(11,12)){
    lines(lowess(conc.mal, conc.mal.means.r0[,i] - r0.In(0)[3], f = 0.01), lwd=2, col = i-6)
  }
  for(i in c(11,12)){
    lines(lowess(conc.mal, 
                 (conc.mal.means.r0[,i] + conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i-6, lty=2)
    lines(lowess(conc.mal, 
                 (conc.mal.means.r0[,i] - conc.mal.sds.r0[,i] - r0.In(0)[3]), f = 0.01), 
          col = i-6, lty=2)
  }
  legend('bottomleft', lwd=2, col = c(3,5,6), cex = 0.7, bty = 'n',
         legend = c(expression(paste(mu[P], sep = '')),
                    expression(paste(mu[P], ' & ', pi[C], ' & ', pi[M])),
                    expression(paste(mu[P], ' & f'[N], ' & ', mu[N]))))
  title('Malathion combined predator & snail/larval effects')
  
#Plot total r0 of q function ###############
  #with reduction in snail reproduction
  plot(lowess(conc.mal, conc.mal.means.r0[,13] - r0.In(0)[3], f = 0.01), 
       lwd=2, type = 'l', col = 2, ylim = c(-4,4), xlim = c(0, 200000),
       xlab = 'Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  lines(lowess(conc.mal, 
               (conc.mal.means.r0[,13] + conc.mal.sds.r0[,13] - r0.In(0)[3]), f = 0.01), 
        col = 2, lty=2)
  lines(lowess(conc.mal, 
               (conc.mal.means.r0[,13] - conc.mal.sds.r0[,13] - r0.In(0)[3]), f = 0.01), 
        col = 2, lty=2)
  title('Malathion - all response functions incorporated')
  
  #Restrict to 25ppm
  plot(lowess(conc.mal, conc.mal.means.r0[,13] - r0.In(0)[3], f = 0.01), 
       lwd=2, type = 'l', col = 2, ylim = c(-4,4), xlim = c(0, 25000),
       xlab = 'Malathion (ppb)', ylab = expression(paste(Delta, R[0], '(q)')))
  lines(lowess(conc.mal, 
               (conc.mal.means.r0[,13] + conc.mal.sds.r0[,13] - r0.In(0)[3]), f = 0.01), 
        col = 2, lty=2)
  lines(lowess(conc.mal, 
               (conc.mal.means.r0[,13] - conc.mal.sds.r0[,13] - r0.In(0)[3]), f = 0.01), 
        col = 2, lty=2)
  title('Malathion - all response functions incorporated')