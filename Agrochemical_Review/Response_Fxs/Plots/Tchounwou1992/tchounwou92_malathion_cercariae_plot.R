#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/tchounwou92_malathion_cercariae_fit.R")
#Tchounwou Data plotted ############
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1992/tchounwou92_malathion_cercariae_data.png")

plot(cerc.mal$time_hrs[cerc.mal$conc==0], cerc.mal$surv[cerc.mal$conc==0]/100, pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tch92.piC.mal, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)

  for(i in 2:length(unique(cerc.mal$conc[cerc.mal$chem == 'malathion']))){
    points(cerc.mal$time_hrs[cerc.mal$conc==unique(cerc.mal$conc)[i]], 
           cerc.mal$surv[cerc.mal$conc==unique(cerc.mal$conc)[i]]/100, pch=16,
           col = i)
    lines(time, predict(tch92.piC.mal, 
                        data.frame(time_hrs=time, conc = unique(cerc.mal$conc)[i])), lty = 2, col = i)
  }
legend('topright', title = 'Mal(ppm)', legend = c(0,50,100,150,200,250), pch = c(17,rep(16,4)),
       col = c(1:6), cex=0.7, bty = 'n')

dev.off()

#Data frame of LL.2 parameters across maltahion concentrations ##################    
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1992/tchounwou92_malathion_cercariae_dose_pars.png")

  plot(parms.df$mal, parms.df$e, pch = 16, xlab = 'malathion (ppm)', ylab = 'LL.2 Parameters',
       ylim = c(0, 15))
    points(parms.df$mal, parms.df$b, pch = 17, col=2)
    for(i in 1:length(parms.df$mal)){
      segments(x0 = parms.df$mal[i], y0 = parms.df$e[i] + parms.df$e.se[i],
               x1 = parms.df$mal[i], y1 = parms.df$e[i] - parms.df$e.se[i])
      segments(x0 = parms.df$mal[i]+1000, y0 = parms.df$b[i] + parms.df$b.se[i],
               x1 = parms.df$mal[i]+1000, y1 = parms.df$b[i] - parms.df$b.se[i], col=2)
    }
    legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
           cex = 0.8, bty = 'n', title = 'Observed points')
    
    lines(seq(0,300, 3),sapply(seq(0,300, 3), fx.e.mod)[1,], lty = 2)
    lines(seq(0,300, 3),sapply(seq(0,300, 3), fx.e.mod)[2,], lty = 3)
    lines(seq(0,300, 3),sapply(seq(0,300, 3), fx.e.mod)[3,], lty = 3)
    
    lines(seq(0,300, 3),sapply(seq(0,300, 3), fx.e.mod2)[1,], lty = 2, col = 3)
    lines(seq(0,300, 3),sapply(seq(0,300, 3), fx.e.mod2)[2,], lty = 3, col = 3)
    lines(seq(0,300, 3),sapply(seq(0,300, 3), fx.e.mod2)[3,], lty = 3, col = 3)
  
    lines(seq(0,300, 3),sapply(seq(0,300, 3), fx.b.mod)[1,], lty = 2, col = 2)
    lines(seq(0,300, 3),sapply(seq(0,300, 3), fx.b.mod)[2,], lty = 3, col = 2)
    lines(seq(0,300, 3),sapply(seq(0,300, 3), fx.b.mod)[3,], lty = 3, col = 2)
  
    
  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  

dev.off() 

#Qualitative model validation ###############
#Regenerate plot of observed data
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1992/tchounwou92_malathion_cercariae_validation.png")

plot(cerc.mal$time_hrs[cerc.mal$conc==0], cerc.mal$surv[cerc.mal$conc==0]/100, pch=16, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24), 
     main = expression(paste(pi[C], ' Model sim-val')))
  for(i in 2:length(unique(cerc.mal$conc[cerc.mal$chem == 'malathion']))){
    points(cerc.mal$time_hrs[cerc.mal$conc==unique(cerc.mal$conc[cerc.mal$chem == 'malathion'])[i]], 
           cerc.mal$surv[cerc.mal$conc==unique(cerc.mal$conc[cerc.mal$chem == 'malathion'])[i]]/100, pch=16,
           col = i)
  }
legend('topright', title = 'Mal(ppm)', legend = c(0,50,100,150,200,250), pch = c(17,rep(16,4)),
       col = c(1:6), cex=0.7, bty = 'n')

#plot model predictions
for(i in c(0,50,100,150,200,250)*1000){
  c = i/50000 + 1
  print(c)
    replicate(10, pred.fx.plot(In = i, clr = c))
  }
  
dev.off()
#plot model output compared to observed points
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1992/tchounwou92_malathion_piC_validation.png")

plot(parms.df$mal, tch92.mal.piC.aucs/tch92.mal.piC.aucs[1], pch = 16, ylim = c(0,1),
     xlab = 'Malathion (ppm)', ylab = expression(paste(pi[C], ' estimate', sep = ' ')))

  set.seed = 43093

  points(seq(0, 300, 1), sapply(seq(0, 300, 1)*1000, piC.tch92_mal_unc), pch = 5, col=4, cex = 0.5)
  
  legend("bottomleft", legend = c("Observed", "Simulated"), pch = c(16, 5), col = c(1,4), bty = 'n', cex = 0.7)

dev.off()  