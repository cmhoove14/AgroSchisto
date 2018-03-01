#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Ghaffar 2016 LARVAL (S.mansoni) data
source("Agrochemical_Review/Response_Fxs/Ghaffar2016_butralin_gylphosate_pendimethalin_cercariae_fit.R")

#butralin cercariae raw data ###############
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_raw_cercariae_mortality_data_butralin.png")
  plot(cer.ctrl$time_hrs, cer.ctrl$surv, ylim = c(0,1), xlim = c(0,24), pch=16, 
       xlab = 'time(hrs)', ylab = 'prop surviving')
    lines(time, L.3.fx(lc50 = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], t = time), lty = 2)
    
    for(i in 1:length(unique(cer.but$conc))){
      points(cer.but$time_hrs[cer.but$conc==unique(cer.but$conc)[i]], 
             cer.but$surv[cer.but$conc==unique(cer.but$conc)[i]], pch=17,
             col = i+1)
      lines(time, L.3.fx(lc50 = gaf.but.drc.cer$coefficients[i+5], slp = gaf.but.drc.cer$coefficients[i], t = time), lty = 2, col = i+1)
    }
  legend('topright', title = 'Butralin (ppb)', legend = but.vals, pch = c(16,rep(17,5)),
         col = c(1:6), cex=0.7, bty = 'n')
  title('Ghaffar2016 butralin toxicity to cercariae')
  
dev.off()  

#Fitted d-r function parameters ###########
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_fitted_d-r_parameters_and_functions_butralin.png")
  plot(but.cer$butralin, but.cer$e, pch = 16, xlab = 'butralin (ppb)', 
       ylab = 'LL.2 Parameters (cercariae)', ylim = c(0,10))
  
    points(but.cer$butralin, but.cer$b, pch = 17, col=2)
  
  for(i in 1:length(but.cer$butralin)){
    segments(x0 = but.cer$butralin[i], y0 = but.cer$e[i] + but.cer$e.se[i],
             x1 = but.cer$butralin[i], y1 = but.cer$e[i] - but.cer$e.se[i])
    segments(x0 = but.cer$butralin[i], y0 = but.cer$b[i] + but.cer$b.se[i],
             x1 = but.cer$butralin[i], y1 = but.cer$b[i] - but.cer$b.se[i], col=2)
  } 

    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.cer.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.cer.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.cer.pred, simplify = T)[3,], lty = 3)

    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.cer.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.cer.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.cer.pred2, simplify = T)[3,], lty = 3, col=3)

    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.cer.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.cer.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.cer.pred.b, simplify = T)[3,], lty = 3, col = 2)

legend('topleft', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
       legend = c('LC50 - Linear',
                  'LC50 - Exponential',
                  'Slp - linear',
                  '95% CI'), cex = 0.7)  
legend('top', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
       cex = 0.8, bty = 'n', title = 'Observed points')

dev.off()

#plot to test function ########
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_function_simulate_cercariae_butralin.png")

plot(but.cer$butralin, gaf.but.aucs.cer/gaf.but.aucs.cer[1],
     xlab = 'Butralin (ppm)', ylab = 'relative AUC (cercariae-hours)',
     pch = 16, ylim = c(0,1))

set.seed(43093)

  points(seq(0, 18000,100)/2000, sapply(seq(0, 18000,100)/2, piC.ghaf_butr.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0, 18000,100)/2000, sapply(seq(0, 18000,100)/2, piC.ghaf_butr.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)
  
  legend('bottomleft', pch = c(16,5,5), col = c(1,2,4), legend = c("observed", "log-linear lc50 fit", "linear lc50 fit"),
         cex = 0.7, bty = 'n')
  

dev.off()

#glyphosate cercariae raw data ###############
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_raw_cercariae_mortality_data_glyphosate.png")

  plot(cer$time_hrs[cer$conc==0], cer$surv[cer$conc==0], 
       pch=16, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
    lines(time, L.3.fx(lc50 = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], t = time), lty = 2)
  
  for(i in 1:length(unique(cer.gly$conc))){
    points(cer.gly$time_hrs[cer.gly$conc==unique(cer.gly$conc)[i]], 
           cer.gly$surv[cer.gly$conc==unique(cer.gly$conc)[i]], pch=17,
           col = i+1)
    lines(time, L.3.fx(lc50 = gaf.gly.drc.cer$coefficients[i+5], slp = gaf.gly.drc.cer$coefficients[i], t = time), lty = 2, col = i+1)
  }
  legend('topright', title = 'Glyphosate (ppm)', legend = gly.vals, pch = c(16,rep(17,5)),
         col = c(1:6), cex=0.7, bty = 'n')
  title('Ghaffar2016 Glyphosate toxicity to cercariae')
  
dev.off()

#Fitted d-r function parameters ########
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_fitted_d-r_parameters_and_functions_glyphosate.png")

plot(gly.cer$glyphosate, gly.cer$e, pch = 16, xlab = 'glyphosate (ppb)', 
     ylab = 'LL.2 Parameters (cercariae)', ylim = c(0,6))

  points(gly.cer$glyphosate, gly.cer$b, pch = 17, col=2)

for(i in 1:length(gly.cer$glyphosate)){
  segments(x0 = gly.cer$glyphosate[i], y0 = gly.cer$e[i] + gly.cer$e.se[i],
           x1 = gly.cer$glyphosate[i], y1 = gly.cer$e[i] - gly.cer$e.se[i])
  segments(x0 = gly.cer$glyphosate[i], y0 = gly.cer$b[i] + gly.cer$b.se[i],
           x1 = gly.cer$glyphosate[i], y1 = gly.cer$b[i] - gly.cer$b.se[i], col=2)
} 

    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.cer.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.cer.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.cer.pred, simplify = T)[3,], lty = 3)

    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.cer.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.cer.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.cer.pred2, simplify = T)[3,], lty = 3, col=3)

    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.cer.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.cer.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.cer.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('topleft', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('top', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

dev.off()

#plot to test function ########
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_function_simulate_cercariae_glyphosate.png")

plot(gly.cer$glyphosate, gaf.gly.aucs.cer/gaf.gly.aucs.cer[1],
     xlab = 'Glyphosate (ppm)', ylab = 'relative AUC (cercariae-hours)',
     pch = 16, ylim = c(0,1))

set.seed(43093)

  points(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3, piC.ghaf_gly.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3, piC.ghaf_gly.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)

  legend('bottomleft', pch = c(16,5,5), col = c(1,2,4), legend = c("observed", "log-linear lc50 fit", "linear lc50 fit"),
         cex = 0.7, bty = 'n')
  
dev.off()
#pendimethalin cercariae raw data ###############
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_raw_cercariae_mortality_data_pendimethalin.png")

plot(cer$time_hrs[cer$conc==0], cer$surv[cer$conc==0], 
     pch=16, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, L.3.fx(lc50 = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], t = time), lty = 2)
  
  for(i in 1:length(unique(cer.pen$conc))){
    points(cer.pen$time_hrs[cer.pen$conc==unique(cer.pen$conc)[i]], 
           cer.pen$surv[cer.pen$conc==unique(cer.pen$conc)[i]], pch=17,
           col = i+1)
    lines(time, L.3.fx(lc50 = gaf.pen.drc.cer$coefficients[i+5], slp = gaf.pen.drc.cer$coefficients[i], t = time), lty = 2, col = i+1)
  }
    legend('topright', title = 'Pendimethalin (ppm)', legend = pen.vals, pch = c(16,rep(17,5)),
           col = c(1:6), cex=0.7, bty = 'n')
    title('Ghaffar2016 Pendimethalin toxicity to cercariae')

dev.off()
#Get estimate of cercariae-hrs for each concentration    
gaf.pen.aucs.cer = as.numeric()
gaf.pen.aucs.cer[1] = integrate(f = L.3.fx, lc50 = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], 
                                lower=0, upper=24)[1]$value  

  for(j in 1:length(unique(cer.pen$conc))){
    gaf.pen.aucs.cer[j+1] = integrate(f = L.3.fx, lc50 = gaf.pen.drc.cer$coefficients[j+5], slp = gaf.pen.drc.cer$coefficients[j],
                                      lower=0, upper=24)[1]$value  
  }

#Fitted d-r function parameters ######## 
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_fitted_d-r_parameters_and_functions_pendimethalin.png")

plot(pen.cer$pendimethalin, pen.cer$e, pch = 16, xlab = 'Pendimethalin (ppb)', 
     ylab = 'LL.2 Parameters (cercariae)', ylim = c(0,6))

  points(pen.cer$pendimethalin, pen.cer$b, pch = 17, col=2)

  for(i in 1:length(pen.cer$pendimethalin)){
    segments(x0 = pen.cer$pendimethalin[i], y0 = pen.cer$e[i] + pen.cer$e.se[i],
             x1 = pen.cer$pendimethalin[i], y1 = pen.cer$e[i] - pen.cer$e.se[i])
    segments(x0 = pen.cer$pendimethalin[i], y0 = pen.cer$b[i] + pen.cer$b.se[i],
             x1 = pen.cer$pendimethalin[i], y1 = pen.cer$b[i] - pen.cer$b.se[i], col=2)
  } 

    lines(seq(0,4000,100)/1000, sapply(seq(0,4000,100)/1000, pen.cer.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,4000,100)/1000, sapply(seq(0,4000,100)/1000, pen.cer.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,4000,100)/1000, sapply(seq(0,4000,100)/1000, pen.cer.pred, simplify = T)[3,], lty = 3)

    lines(seq(0,4000,100)/1000, sapply(seq(0,4000,100)/1000, pen.cer.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,4000,100)/1000, sapply(seq(0,4000,100)/1000, pen.cer.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,4000,100)/1000, sapply(seq(0,4000,100)/1000, pen.cer.pred2, simplify = T)[3,], lty = 3, col=3)
  
    lines(seq(0,4000,100)/1000, sapply(seq(0,4000,100)/1000, pen.cer.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,4000,100)/1000, sapply(seq(0,4000,100)/1000, pen.cer.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,4000,100)/1000, sapply(seq(0,4000,100)/1000, pen.cer.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

dev.off()

#plot to test function ########
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_function_simulate_cercariae_pendimethalin.png")

plot(pen.cer$pendimethalin, gaf.pen.aucs.cer/gaf.pen.aucs.cer[1],
     xlab = 'pendimethalin (ppm)', ylab = 'relative AUC (cercariae-hours)',
     pch = 16, ylim = c(0,1))

set.seed(43093)

  points(seq(0,4000,20)/1000, sapply(seq(0,4000,20), piC.ghaf_pen.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,4000,20)/1000, sapply(seq(0,4000,20), piC.ghaf_pen.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)
  
  legend('bottomleft', pch = c(16,5,5), col = c(1,2,4), legend = c("observed", "log-linear lc50 fit", "linear lc50 fit"),
         cex = 0.7, bty = 'n')
  
dev.off()  