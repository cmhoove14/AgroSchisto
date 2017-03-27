#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############


require(deSolve)
require(graphics)
require(ggplot2)

#Simulations with 100 transmission parameter sets, k=2, double round of MDA annually #####
load("Elimination_Feasibility/det_sim_k200_mda1&2.RData")

#Some post-processing of data
lowest = as.data.frame(det.runs.k2.mda1[ , , 1, 1]  )
  colnames(lowest) = c('time', 'S', 'E', 'I', 'Wt', 'Wu')
  
lowest$W = cov*lowest$Wt + (1-cov)*lowest$Wu
  plot(lowest$time/365, lowest$W, type = 'l', lwd = 2, ylim = c(0,80))
    lines(lowest$time/365, lowest$Wt, type = 'l', lty = 2)
    lines(lowest$time/365, lowest$Wu, type = 'l', lty = 3)
    #points(c(1:20), w.pre.k2.mda1[1,], pch = 16, col=3)
    #points(c(1:20), w.pos.k2.mda1[1,], pch = 16, col=4)
    
lowest.pdd = as.data.frame(det.runs.k2.mda1[ , , 1, 2])
  colnames(lowest.pdd) = c('time', 'S', 'E', 'I', 'Wt', 'Wu')
  
lowest.pdd$W = cov*lowest.pdd$Wt + (1-cov)*lowest.pdd$Wu
  plot(lowest.pdd$time/365, lowest.pdd$W, type = 'l', lwd = 2, ylim = c(0,80))
    lines(lowest.pdd$time/365, lowest.pdd$Wt, type = 'l', lty = 2)
    lines(lowest.pdd$time/365, lowest.pdd$Wu, type = 'l', lty = 3)

#Find BBR for PDD-free and PDD models ###########
#PDD-free    
bbr = (1/w.pos.k2.mda1[,c(1:19)])*(w.pre.k2.mda1[,c(2:20)] - w.pos.k2.mda1[,c(1:19)])
  bbr.mean = colMeans(bbr)
  bbr.sd = apply(bbr, 2, sd)

#PDD     
bbr.pdd = (1/w.pos.k2.mda1.pdd[,c(1:19)])*(w.pre.k2.mda1.pdd[,c(2:20)] - w.pos.k2.mda1.pdd[,c(1:19)])
  bbr.pdd.mean = colMeans(bbr.pdd)
  bbr.pdd.sd = apply(bbr.pdd, 2, sd)  
 
#plot   
    plot(c(1:19), bbr.mean, pch = 16, xlab = 'time (yrs)', ylab = 'BBR', ylim = c(-0.35, 0.35))
      for(i in 1:length(bbr.sd)){
        segments(x0 = i, x1 = i, 
                 y0 = bbr.mean[i] + bbr.sd[i], y1 = bbr.mean[i] - bbr.sd[i])
      }
    
    points(c(1:19), bbr.pdd.mean, pch = 16, col = 2)
      for(i in 1:length(bbr.pdd.sd)){
        segments(x0 = i, x1 = i, 
                 y0 = bbr.pdd.mean[i] + bbr.pdd.sd[i], y1 = bbr.pdd.mean[i] - bbr.pdd.sd[i],
                 col = 2)
      }
    
#Seems like something strange going on with quick drop in first year, is BBR different in treated population? ###########
#get w-treated pre and post treatment
#PDD-free
  wt.pos.k2.mda1 = det.runs.k2.mda1[seq(366, 22265, 365), 5, , 1][c(1:19),]
  wt.pre.k2.mda1 = det.runs.k2.mda1[seq(365, 22265, 365), 5, , 1][c(2:20),]

  bbr.wt = (1/wt.pos.k2.mda1) * (wt.pre.k2.mda1 - wt.pos.k2.mda1)
    
    bbr.mean.wt = rowMeans(bbr.wt)
    bbr.sd.wt = apply(bbr.wt, 1, sd)
    
    
#PDD 
  wt.pos.k2.mda1.pdd = det.runs.k2.mda1[seq(366, 22265, 365), 5, , 2][c(1:19),]
  wt.pre.k2.mda1.pdd = det.runs.k2.mda1[seq(365, 22265, 365), 5, , 2][c(2:20),]
    
  bbr.wt.pdd = (1/wt.pos.k2.mda1.pdd) * (wt.pre.k2.mda1.pdd - wt.pos.k2.mda1.pdd)
    
    bbr.mean.wt.pdd = rowMeans(bbr.wt.pdd)
    bbr.sd.wt.pdd = apply(bbr.wt.pdd, 1, sd)

    
#plot   
  plot(c(1:19), bbr.mean.wt, pch = 16, xlab = 'time (yrs)', ylab = 'BBR (Wt)', ylim = c(0,70), cex = 0.7)
    for(i in 1:length(bbr.sd.wt)){
      segments(x0 = i, x1 = i, 
               y0 = bbr.mean.wt[i] + bbr.sd.wt[i], y1 = bbr.mean.wt[i] - bbr.sd.wt[i])
    }
    
  points(c(1:19), bbr.mean.wt.pdd, pch = 16, col = 2, cex = 0.7)
    for(i in 1:length(bbr.pdd.sd)){
      segments(x0 = i, x1 = i, 
               y0 = bbr.pdd.mean[i] + bbr.pdd.sd[i], y1 = bbr.pdd.mean[i] - bbr.pdd.sd[i],
               col = 2)
    }
    
#Estimate slope of BBR line (aka elimination feasibility estimator) for PDD-free and PDD models #######
#PDD-free model
  bbr.df = data.frame('time' = c(1:19),
                      'bbr' = bbr.mean,
                      'bbr.sd' = bbr.sd,
                      'bbr.slp' = 0,
                      'bbr.slp.025' = 0,
                      'bbr.slp.975' = 0)  

    for(m in 3:length(bbr.mean)){ #only calculate slope with 3+ data points
      mod = lm(bbr[1:m] ~ time[1:m], data = bbr.df, weights = bbr.sd[1:m]^-1)
        bbr.df[m, 4] = mod$coefficients[2]
        bbr.df[m, c(5:6)] = confint(mod, level = 0.95)[2,]
    }
  
#PDD model  
  bbr.pdd.df = data.frame('time' = c(1:19),
                          'bbr' = bbr.pdd.mean,
                          'bbr.sd' = bbr.pdd.sd,
                          'bbr.slp' = 0,
                          'bbr.slp.025' = 0,
                          'bbr.slp.975' = 0)  
    
    for(m in 3:length(bbr.pdd.mean)){
      mod.pdd = lm(bbr[1:m] ~ time[1:m], data = bbr.pdd.df, weights = bbr.sd[1:m]^-1)
        bbr.pdd.df[m, 4] = mod.pdd$coefficients[2]
        bbr.pdd.df[m, c(5:6)] = confint(mod.pdd, level = 0.95)[2,]
    }  

#plot  
  plot(bbr.df$time, bbr.df$bbr.slp, type = 'l', lwd = 2
       , xlab = 'time (yrs)', ylim = c(-0.3, 0.2), xlim=c(2,20), 
       ylab = expression(paste('Elimination Feasibility Estimator ( ', epsilon, ')',sep = '')))
    lines(bbr.pdd.df$time, bbr.pdd.df$bbr.slp, lwd = 2, col=2)
  #95% confints
    lines(bbr.df$time, bbr.df$bbr.slp.025, lty = 2)
    lines(bbr.df$time, bbr.df$bbr.slp.975, lty = 2)
    lines(bbr.pdd.df$time, bbr.pdd.df$bbr.slp.025, lty = 2, col=2)
    lines(bbr.pdd.df$time, bbr.pdd.df$bbr.slp.975, lty = 2, col=2)
 
#Test for differences between PDD and PDD-free ########
  bbr.diff = rbind(bbr.df, bbr.pdd.df)
    bbr.diff$pdd = c(rep('no', 19), rep('yes', 19))
    p.diff = as.numeric()
    
  #include presence/absence of pdd as a regression term in the model; store p-value associated with it  
    
    for(m in 3:length(bbr.pdd.mean)){
      test.seq = c(1:m)
      test.df = subset(bbr.diff, time %in% test.seq)
        mod.diff = lm(bbr ~ time + pdd, data = test.df, weights = bbr.sd^-1)
        p.diff[m-2] = summary(mod.diff)$coefficients[3,4]
    }  
    lines(c(3:nrow(bbr.df)), p.diff, lty=3, col=3, lwd=2)
    
  #Add final legend  
    legend('bottomright', bty='n', cex = 0.8, lty = c(1,1,2,2,3), lwd = c(2,2,1,1,3), col = c(1,2,1,2,3),
           legend = c('PDD-free model', 'PDD model', '95% CI PDD-free', '95% CI PDD', 'p-value')) 
    
#Simulations with 100 transmission parameter sets, k=0.08, single round of MDA annually ####
  load("Elimination_Feasibility/det_sim_k008_mda1.RData")
    
#Some post-processing of data; lowest transmisison intensity plotting
  lowest = as.data.frame(det.runs.k008.mda1[ , , 1, 1]  )
    colnames(lowest) = c('time', 'S', 'E', 'I', 'Wt', 'Wu')
    
  lowest$W = cov*lowest$Wt + (1-cov)*lowest$Wu
    plot(lowest$time/365, lowest$W, type = 'l', lwd = 2, ylim = c(0,80))
      lines(lowest$time/365, lowest$Wt, type = 'l', lty = 2)
      lines(lowest$time/365, lowest$Wu, type = 'l', lty = 3)
    
  lowest.pdd = as.data.frame(det.runs.k008.mda1[ , , 1, 2])
    colnames(lowest.pdd) = c('time', 'S', 'E', 'I', 'Wt', 'Wu')
    
  lowest.pdd$W = cov*lowest.pdd$Wt + (1-cov)*lowest.pdd$Wu
    plot(lowest.pdd$time/365, lowest.pdd$W, type = 'l', lwd = 2, ylim = c(0,80))
      lines(lowest.pdd$time/365, lowest.pdd$Wt, type = 'l', lty = 2)
      lines(lowest.pdd$time/365, lowest.pdd$Wu, type = 'l', lty = 3)
  
lowest$W[lowest$time == 366]      
      
#Find BBR for PDD-free and PDD models #########
  #PDD-free    
  bbr = (1/w.pos.k008.mda1[,c(1:19)])*(w.pre.k008.mda1[,c(2:20)] - w.pos.k008.mda1[,c(1:19)])
    bbr.mean = colMeans(bbr)
    bbr.sd = apply(bbr, 2, sd)
      
  #PDD     
  bbr.pdd = (1/w.pos.k008.mda1.pdd[,c(1:19)])*(w.pre.k008.mda1.pdd[,c(2:20)] - w.pos.k008.mda1.pdd[,c(1:19)])
    bbr.pdd.mean = colMeans(bbr.pdd)
    bbr.pdd.sd = apply(bbr.pdd, 2, sd)  
      
#plot   
  plot(c(1:19), bbr.mean, pch = 16, xlab = 'time (yrs)', ylab = 'BBR', ylim = c(-0.35, 0.35))
    for(i in 1:length(bbr.sd)){
        segments(x0 = i, x1 = i, 
                 y0 = bbr.mean[i] + bbr.sd[i], y1 = bbr.mean[i] - bbr.sd[i])
    }
      
    points(c(1:19), bbr.pdd.mean, pch = 16, col = 2)
      for(i in 1:length(bbr.pdd.sd)){
        segments(x0 = i, x1 = i, 
                 y0 = bbr.pdd.mean[i] + bbr.pdd.sd[i], y1 = bbr.pdd.mean[i] - bbr.pdd.sd[i],
                 col = 2)
      }

#Estimate slope of BBR line (aka elimination feasibility estimator) for PDD-free and PDD models ###########
#PDD-free model
  bbr.df = data.frame('time' = c(1:19),
                      'bbr' = bbr.mean,
                      'bbr.sd' = bbr.sd,
                      'bbr.slp' = 0,
                      'bbr.slp.025' = 0,
                      'bbr.slp.975' = 0)  

    for(m in 3:length(bbr.mean)){ #only calculate slope with 3+ data points
      mod = lm(bbr[1:m] ~ time[1:m], data = bbr.df, weights = bbr.sd[1:m]^-1)
      bbr.df[m, 4] = mod$coefficients[2]
      bbr.df[m, c(5:6)] = confint(mod, level = 0.95)[2,]
    }

#PDD model  
  bbr.pdd.df = data.frame('time' = c(1:19),
                          'bbr' = bbr.pdd.mean,
                          'bbr.sd' = bbr.pdd.sd,
                          'bbr.slp' = 0,
                          'bbr.slp.025' = 0,
                          'bbr.slp.975' = 0)  

  for(m in 3:length(bbr.pdd.mean)){
    mod.pdd = lm(bbr[1:m] ~ time[1:m], data = bbr.pdd.df, weights = bbr.sd[1:m]^-1)
    bbr.pdd.df[m, 4] = mod.pdd$coefficients[2]
    bbr.pdd.df[m, c(5:6)] = confint(mod.pdd, level = 0.95)[2,]
  }  

#plot  
plot(bbr.df$time, bbr.df$bbr.slp, type = 'l', lwd = 2, 
     xlab = 'time (yrs)', ylim = c(-0.3, 0.2), xlim=c(2,20), 
     ylab = expression(paste('Elimination Feasibility Estimator ( ', epsilon, ')',sep = '')))
  lines(bbr.pdd.df$time, bbr.pdd.df$bbr.slp, lwd = 2, col=2)
  
  #95% confints
    lines(bbr.df$time, bbr.df$bbr.slp.025, lty = 2)
    lines(bbr.df$time, bbr.df$bbr.slp.975, lty = 2)
    lines(bbr.pdd.df$time, bbr.pdd.df$bbr.slp.025, lty = 2, col=2)
    lines(bbr.pdd.df$time, bbr.pdd.df$bbr.slp.975, lty = 2, col=2)

#Test for differences between PDD and PDD-free ###########
#include presence/absence of pdd as a regression term in the model; store p-value associated with it  

  for(m in 3:length(bbr.pdd.mean)){
    test.seq = c(1:m)
    test.df = subset(bbr.diff, time %in% test.seq)
    mod.diff = lm(bbr ~ time + pdd, data = test.df, weights = bbr.sd^-1)
    p.diff[m-2] = summary(mod.diff)$coefficients[3,4]
  }  
    
    lines(c(3:nrow(bbr.df)), p.diff, lty=3, col=3, lwd=2)

#Add final legend  
  legend('bottomright', bty='n', cex = 0.8, lty = c(1,1,2,2,3), lwd = c(2,2,1,1,3), col = c(1,2,1,2,3),
          legend = c('PDD-free model', 'PDD model', '95% CI PDD-free', '95% CI PDD', 'p-value')) 
  
#get epsilon for each sim & test if bbr slope (epsilon) significantly different from 0 over time ###########
  eps.fill = matrix(ncol = length(3:ncol(bbr)),
                    nrow = nrow(bbr))
  eps.pdd.fill = matrix(ncol = length(3:ncol(bbr)),
                        nrow = nrow(bbr))

#get epsilon (slope of BBR profile) from year 3-20 for all 100 sims    
  for(m in 1:ncol(eps.fill)){
    for(i in 1:nrow(bbr)){
      eps.fill[i,m] = lm(bbr[i,c(1:(m+2))] ~ c(1:(m+2)))$coefficients[2]
      eps.pdd.fill[i,m] = lm(bbr.pdd[i,c(1:(m+2))] ~ c(1:(m+2)))$coefficients[2]
    }
    
  } #end function
  
  #descriptive stats for epsilon estimates
    eps.mean = colMeans(eps.fill)
    eps.sd = apply(eps.fill, 2, sd)
    
    eps.pdd.mean = colMeans(eps.pdd.fill)
    eps.pdd.sd = apply(eps.pdd.fill, 2, sd)
    
#vectors to fill with p-values from t tests    
  t.fill.g0 = as.numeric()        #is pdd-free model slope greater than 0?
  t.fill.l0 = as.numeric()        #is pdd-free model slope less than 0?
  t.fill.pdd.l0 = as.numeric()    #is pdd model slope less than 0/
  t.fill.pdd.diff = as.numeric()  #are pdd and pdd-free slopes different?

  for(t in 1:ncol(eps.fill)){
    t.fill.g0[t] = t.test(x = eps.fill[,t], mu = 0, alternative = 'greater')$p.value
    t.fill.l0[t] = t.test(x = eps.fill[,t], mu = 0, alternative = 'less')$p.value
    t.fill.pdd.l0[t] = t.test(x = eps.pdd.fill[,t], mu = 0, alternative = 'less')$p.value
    t.fill.pdd.diff[t] = t.test(x = eps.fill[,t], y = eps.pdd.fill[,t], alternative = 'two.sided')$p.value
  }

#plot epsilon over time 
  
  
  plot(c(3:ncol(bbr)), eps.mean, pch = 16, cex = 0.8, xlim = c(0,20), ylim = c(-0.03, 0.01), xlab = 'time (yrs)',
        ylab = expression(paste('Elimination Feasibility Estimator ( ', epsilon, ')',sep = '')))
    points(c(3:ncol(bbr)), eps.pdd.mean, pch = 16, cex = 0.8, col=2)
    
    #Add error bars as st dev
    for(s in 1:length(eps.sd)){
      segments(x0 = s+2, x1 = s+2, 
               y0 = eps.mean[s] + eps.sd[s], y1 = eps.mean[s] - eps.sd[s])
      segments(x0 = s+2, x1 = s+2, 
               y0 = eps.pdd.mean[s] + eps.pdd.sd[s], y1 = eps.pdd.mean[s] - eps.pdd.sd[s], col = 2)
    }

#Create df for using ggplot
  eps.gg = data.frame(eps = c(eps.mean, eps.pdd.mean),
                      eps.sd = c(eps.sd, eps.pdd.sd),
                      time = rep(c(3:ncol(bbr)), 2),
                      Model = c(rep('PDD-free', length(eps.mean)), rep('PDD', length(eps.mean))))  
  
eps.ggp = ggplot(eps.gg, aes(x = time, y = eps, fill = Model)) +
            theme_bw() +
            theme(plot.margin = unit(c(1.5,0.5,1,1), 'cm')) +
            xlim(0,20) +
            ylim(-0.035,0.01) +
            labs(x = 'time (yrs)',
                 y = expression(paste('Elimination Feasibility Estimator (', epsilon, ')',sep = ''))) +
            geom_bar(position = 'dodge', stat = 'identity', width = 0.8) +
            scale_fill_manual(values = alpha(c('red', 'black'), .6)) +
            geom_errorbar(aes(x = time, ymin = (eps - eps.sd), ymax = (eps + eps.sd)), 
                          width = 0.25, position = position_dodge(width = 0.8)) +
            annotate('text', x = c(3:ncol(bbr)), y = 0.008, label = c(rep('a', 3), rep('b', 14)), size = 3)
eps.ggp
    
#manuscript plot version
  opar<-par()

windows(width = 11)    
  par(mfrow = c(2,2), pin = c(11,7), mar = c(3,4,3,1)+0.1)
    plot(lowest$time/365, lowest$W, type = 'l', lwd = 2, ylim = c(0,80),
         xlab = 'Time (years)', ylab = expression(italic(W)))
      #lines(lowest$time/365, lowest$Wt, type = 'l', lty = 2)
      #lines(lowest$time/365, lowest$Wu, type = 'l', lty = 3)
  
    plot(lowest.pdd$time/365, lowest.pdd$W, type = 'l', lwd = 2, ylim = c(0,80), col=2,
         xlab = 'Time (years)', ylab = expression(italic(W)))
      #lines(lowest.pdd$time/365, lowest.pdd$Wt, type = 'l', lty = 2)
      #lines(lowest.pdd$time/365, lowest.pdd$Wu, type = 'l', lty = 3)
      
    plot(c(1:19), bbr.mean, pch = 16, cex = 0.8, ylim = c(-0.1, 0.35), xlim = c(0,20),
         xlab = 'Time (yrs)', ylab = expression(italic(BBR)))
      for(i in 1:length(bbr.sd)){
        segments(x0 = i, x1 = i, 
                 y0 = bbr.mean[i] + bbr.sd[i], y1 = bbr.mean[i] - bbr.sd[i])
      }
      
      points(c(1:19), bbr.pdd.mean, pch = 16, cex = 0.8, col = 2)
      for(i in 1:length(bbr.pdd.sd)){
        segments(x0 = i, x1 = i, 
                 y0 = bbr.pdd.mean[i] + bbr.pdd.sd[i], y1 = bbr.pdd.mean[i] - bbr.pdd.sd[i],
                 col = 2)
      }
      
      legend('bottomright', legend = c('PDD-free model', 'PDD-model'), 
             pch = 16, cex = 0.8, col = c(1,2), bty='n')
      
      plot(c(3:ncol(bbr)), eps.mean, pch = 16, cex = 0.8, xlim = c(0,20), ylim = c(-0.03, 0.01), xlab = 'time (yrs)',
           ylab = expression(paste('Elimination Feasibility Estimator ( ', epsilon, ')',sep = '')))
        points(c(3:ncol(bbr)), eps.pdd.mean, pch = 16, cex = 0.8, col=2)
      
      #Add error bars as st dev
        for(s in 1:length(eps.sd)){
          segments(x0 = s+2, x1 = s+2, 
                   y0 = eps.mean[s] + eps.sd[s], y1 = eps.mean[s] - eps.sd[s])
          segments(x0 = s+2, x1 = s+2, 
                   y0 = eps.pdd.mean[s] + eps.pdd.sd[s], y1 = eps.pdd.mean[s] - eps.pdd.sd[s], col = 2)
        }
  

#BELOW IS DEPRECATED #########
#Estimate slope of BBR line (aka elimination feasibility estimator) for PDD-free and PDD models
      #PDD-free model
      bbr.df = data.frame('time' = c(1:19),
                          'bbr' = bbr.mean,
                          'bbr.sd' = bbr.sd,
                          'bbr.slp' = 0,
                          'bbr.slp.025' = 0,
                          'bbr.slp.975' = 0)  
      
      for(m in 3:length(bbr.mean)){ #only calculate slope with 3+ data points
        mod = lm(bbr[1:m] ~ time[1:m], data = bbr.df, weights = bbr.sd[1:m]^-1)
        bbr.df[m, 4] = mod$coefficients[2]
        bbr.df[m, c(5:6)] = confint(mod, level = 0.95)[2,]
      }
      
      #PDD model  
      bbr.pdd.df = data.frame('time' = c(1:19),
                              'bbr' = bbr.pdd.mean,
                              'bbr.sd' = bbr.pdd.sd,
                              'bbr.slp' = 0,
                              'bbr.slp.025' = 0,
                              'bbr.slp.975' = 0)  
      
      for(m in 3:length(bbr.pdd.mean)){
        mod.pdd = lm(bbr[1:m] ~ time[1:m], data = bbr.pdd.df, weights = bbr.sd[1:m]^-1)
        bbr.pdd.df[m, 4] = mod.pdd$coefficients[2]
        bbr.pdd.df[m, c(5:6)] = confint(mod.pdd, level = 0.95)[2,]
      }  
      
      #plot  
      plot(bbr.df$time, bbr.df$bbr.slp, type = 'l', lwd = 2, 
           xlab = 'time (yrs)', ylim = c(-0.3, 0.2), xlim=c(2,20), 
           ylab = expression(paste('Elimination Feasibility Estimator ( ', epsilon, ')',sep = '')))
      lines(bbr.pdd.df$time, bbr.pdd.df$bbr.slp, lwd = 2, col=2)
      
      #95% confints
      lines(bbr.df$time, bbr.df$bbr.slp.025, lty = 2)
      lines(bbr.df$time, bbr.df$bbr.slp.975, lty = 2)
      lines(bbr.pdd.df$time, bbr.pdd.df$bbr.slp.025, lty = 2, col=2)
      lines(bbr.pdd.df$time, bbr.pdd.df$bbr.slp.975, lty = 2, col=2)
      
      #Test for differences between PDD and PDD-free
      #include presence/absence of pdd as a regression term in the model; store p-value associated with it  
      
      for(m in 3:length(bbr.pdd.mean)){
        test.seq = c(1:m)
        test.df = subset(bbr.diff, time %in% test.seq)
        mod.diff = lm(bbr ~ time + pdd, data = test.df, weights = bbr.sd^-1)
        p.diff[m-2] = summary(mod.diff)$coefficients[3,4]
      }  
      
      lines(c(3:nrow(bbr.df)), p.diff, lty=3, col=3, lwd=2)
      
      #Add final legend  
      legend('bottomright', bty='n', cex = 0.8, lty = c(1,1,2,2,3), lwd = c(2,2,1,1,3), col = c(1,2,1,2,3),
             legend = c('PDD-free model', 'PDD model', '95% CI PDD-free', '95% CI PDD', 'p-value')) 
      