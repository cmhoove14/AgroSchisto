#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/monte2016_glyphosate_larvae_fit.R")

#Monte 2016 cercariae Data plotted ############
png("Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_cercariae_data.png")

plot(monte_dat$time_hrs[monte_dat$conc == 0 &monte_dat$larvae == "cercariae"],
     monte_dat$surv[monte_dat$conc == 0 &monte_dat$larvae == "cercariae"]/100, pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
   lines(time, predict(tch92.piC.mal, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)

  for(i in 1:length(unique(cerc$conc))){
    points(cerc$time_hrs[cerc$conc==unique(cerc$conc)[i]], 
           cerc$surv[cerc$conc==unique(cerc$conc)[i]]/100, pch=16,
           col = i+1)
    lines(time, predict(mont16.piC.gly, 
                        data.frame(time_hrs=time, conc = unique(cerc$conc)[i])), lty = 2, col = i+1)
  }
legend('topright', title = 'Glyphosate (ppm)', legend = c(0,unique(cerc$conc)), pch = c(rep(17, 16,5)),
       col = c(1:6), cex=0.7, bty = 'n')

dev.off()

#Data frame of LL.2 parameters across glyphosate cercariae concentrations ##################    
png("Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_cercariae_dose_pars.png")

  plot(mont16_pars$gly, mont16_pars$e, pch = 16, xlab = 'Glyphosate (ppm)', ylab = 'LL.2 Parameters',
       ylim = c(0, 15))
    points(mont16_pars$gly, mont16_pars$b, pch = 17, col=2)
    for(i in 1:length(mont16_pars$gly)){
      segments(x0 = mont16_pars$gly[i], y0 = mont16_pars$e[i] + mont16_pars$e.se[i],
               x1 = mont16_pars$gly[i], y1 = mont16_pars$e[i] - mont16_pars$e.se[i])
      segments(x0 = mont16_pars$gly[i], y0 = mont16_pars$b[i] + mont16_pars$b.se[i],
               x1 = mont16_pars$gly[i], y1 = mont16_pars$b[i] - mont16_pars$b.se[i], col=2)
    }
    legend('topleft', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
           cex = 0.8, bty = 'n', title = 'Observed points')
    
    lines(seq(0,200, 2),sapply(seq(0,200, 2), fx.e.mod)[1,], lty = 2)
    lines(seq(0,200, 2),sapply(seq(0,200, 2), fx.e.mod)[2,], lty = 3)
    lines(seq(0,200, 2),sapply(seq(0,200, 2), fx.e.mod)[3,], lty = 3)
    
    lines(seq(0,200, 2),sapply(seq(0,200, 2), fx.e.mod2)[1,], lty = 2, col = 3)
    lines(seq(0,200, 2),sapply(seq(0,200, 2), fx.e.mod2)[2,], lty = 3, col = 3)
    lines(seq(0,200, 2),sapply(seq(0,200, 2), fx.e.mod2)[3,], lty = 3, col = 3)
  
    lines(seq(0,200, 2),sapply(seq(0,200, 2), fx.b.mod)[1,], lty = 2, col = 2)
    lines(seq(0,200, 2),sapply(seq(0,200, 2), fx.b.mod)[2,], lty = 3, col = 2)
    lines(seq(0,200, 2),sapply(seq(0,200, 2), fx.b.mod)[3,], lty = 3, col = 2)
  
    
  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  

dev.off() 

#Qualitative cercariae model validation ###############
#Regenerate plot of observed data
png("Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_cercariae_validation.png")

plot(monte_dat$time_hrs[monte_dat$conc == 0 &monte_dat$larvae == "cercariae"],
     monte_dat$surv[monte_dat$conc == 0 &monte_dat$larvae == "cercariae"]/100, pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))

  for(i in 1:length(unique(cerc$conc))){
    points(cerc$time_hrs[cerc$conc==unique(cerc$conc)[i]], 
           cerc$surv[cerc$conc==unique(cerc$conc)[i]]/100, pch=16,
           col = i+1)
  }

#plot model predictions
  replicate(10, monte16_pred_ts(He = 0*1000, clr = 1))
  replicate(10, monte16_pred_ts(He = unique(cerc$conc)[1]*1000, clr = 2))
  replicate(10, monte16_pred_ts(He = unique(cerc$conc)[2]*1000, clr = 3))
  replicate(10, monte16_pred_ts(He = unique(cerc$conc)[3]*1000, clr = 4))
  replicate(10, monte16_pred_ts(He = unique(cerc$conc)[4]*1000, clr = 5))
  replicate(10, monte16_pred_ts(He = unique(cerc$conc)[5]*1000, clr = 6))

dev.off()
#plot cercariae model output compared to observed AUC points #########
png("Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_cercariae_piC_validation.png")

plot(mont16_pars$gly[-6], mont16.gly.piC.aucs/mont16.gly.piC.aucs[1], pch = 16, ylim = c(0,1),
     xlab = 'Glyphosate (ppm)', ylab = expression(paste(pi[C], ' estimate', sep = ' ')))

  set.seed = 43093

  points(seq(0, 200, 1), sapply(seq(0, 200, 1)*1000, piC.mont16_gly_unc), pch = 5, col=4, cex = 0.5)
  
  legend("bottomleft", legend = c("Observed", "Simulated"), pch = c(16, 5), col = c(1,4), bty = 'n', cex = 0.7)

dev.off()  

#Monte 2016 miracidia Data plotted ############
png("Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_miracidia_data.png")

plot(monte_dat$time_hrs[monte_dat$conc == 0 & monte_dat$larvae == "miracidia"],
     monte_dat$surv[monte_dat$conc == 0 &monte_dat$larvae == "miracidia"]/100, pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
   lines(time, predict(mont16.piM.gly, 
                       data.frame(time_hrs=time, conc = 0)), lty = 2)

  for(i in 2:length(unique(mir$conc))){
    points(mir$time_hrs[mir$conc==unique(mir$conc)[i]], 
           mir$surv[mir$conc==unique(mir$conc)[i]]/100, pch=16,
           col = i)
    lines(time, predict(mont16.piM.gly, 
                        data.frame(time_hrs=time, conc = unique(mir$conc)[i])), lty = 2, col = i)
  }
legend('topright', title = 'Glyphosate (ppm)', legend = unique(mir$conc), pch = c(17, rep(16,5)),
       col = c(1:6), cex=0.7, bty = 'n')

dev.off()

#Data frame of LL.2 parameters across glyphosate miracidia concentrations ##################    
png("Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_miracidia_dose_pars.png")

  plot(mont16_pars_mir$gly, mont16_pars_mir$e, pch = 16, xlab = 'Glyphosate (ppm)', ylab = 'LL.2 Parameters',
       ylim = c(0, 15))
    points(mont16_pars_mir$gly, mont16_pars_mir$b, pch = 17, col=2)
    for(i in 1:length(mont16_pars_mir$gly)){
      segments(x0 = mont16_pars_mir$gly[i], y0 = mont16_pars_mir$e[i] + mont16_pars_mir$e.se[i],
               x1 = mont16_pars_mir$gly[i], y1 = mont16_pars_mir$e[i] - mont16_pars_mir$e.se[i])
      segments(x0 = mont16_pars_mir$gly[i], y0 = mont16_pars_mir$b[i] + mont16_pars_mir$b.se[i],
               x1 = mont16_pars_mir$gly[i], y1 = mont16_pars_mir$b[i] - mont16_pars_mir$b.se[i], col=2)
    }
    legend('topleft', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
           cex = 0.8, bty = 'n', title = 'Observed points')
    
    lines(seq(0,24, .2),sapply(seq(0,24, .2), fx.e.mod.mir)[1,], lty = 2)
    lines(seq(0,24, .2),sapply(seq(0,24, .2), fx.e.mod.mir)[2,], lty = 3)
    lines(seq(0,24, .2),sapply(seq(0,24, .2), fx.e.mod.mir)[3,], lty = 3)
    
    lines(seq(0,24, .2),sapply(seq(0,24, .2), fx.e.mod2.mir)[1,], lty = 2, col = 3)
    lines(seq(0,24, .2),sapply(seq(0,24, .2), fx.e.mod2.mir)[2,], lty = 3, col = 3)
    lines(seq(0,24, .2),sapply(seq(0,24, .2), fx.e.mod2.mir)[3,], lty = 3, col = 3)
  
    lines(seq(0,24, .2),sapply(seq(0,24, .2), fx.b.mod.mir)[1,], lty = 2, col = 2)
    lines(seq(0,24, .2),sapply(seq(0,24, .2), fx.b.mod.mir)[2,], lty = 3, col = 2)
    lines(seq(0,24, .2),sapply(seq(0,24, .2), fx.b.mod.mir)[3,], lty = 3, col = 2)
  
    
  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  

dev.off() 


#Qualitative miracidia model validation ###############
#Regenerate plot of observed data
png("Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_miracidia_validation.png")

plot(monte_dat$time_hrs[monte_dat$conc == 0 & monte_dat$larvae == "miracidia"],
     monte_dat$surv[monte_dat$conc == 0 &monte_dat$larvae == "miracidia"]/100, pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))

  for(i in 2:length(unique(mir$conc))){
    points(mir$time_hrs[mir$conc==unique(mir$conc)[i]], 
           mir$surv[mir$conc==unique(mir$conc)[i]]/100, pch=16,
           col = i)
  }

#plot model predictions
for(i in 1:length(unique(mir$conc))){
  replicate(10, monte16_pred_ts_mir(He = unique(mir$conc)[i]*1000, clr = i))
}

dev.off()

#plot miracidia model output compared to observed AUC points #########
png("Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_miracidia_piM_validation.png")

plot(mont16_pars_mir$gly, mont16.gly.piM.aucs/mont16.gly.piM.aucs[1], pch = 16, ylim = c(0,1),
     xlab = 'Glyphosate (ppm)', ylab = expression(paste(pi[M], ' estimate', sep = ' ')))

  set.seed = 43093

  points(seq(0, 24, .1), sapply(seq(0, 24, .1)*1000, piM.mont16_gly_unc), pch = 5, col=4, cex = 0.5)
  
  legend("bottomleft", legend = c("Observed", "Simulated"), pch = c(16, 5), col = c(1,4), bty = 'n', cex = 0.7)

dev.off()  
