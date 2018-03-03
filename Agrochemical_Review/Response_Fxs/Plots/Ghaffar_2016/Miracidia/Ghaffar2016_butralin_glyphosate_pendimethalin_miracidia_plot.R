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
source("Agrochemical_Review/Response_Fxs/Ghaffar2016_butralin_glyphosate_pendimethalin_miracidia_fit.R")

#butralin ###############
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_raw_miracidia_mortality_data_butralin.png")

plot(mir.ctrl$time_hrs, mir.ctrl$surv, ylim = c(0,1), xlim = c(0,24), pch=16, 
     xlab = 'time(hrs)', ylab = 'prop surviving',
     main = 'Abdel-Ghaffar 2016: Butralin toxicity to miracidia')
  lines(time, L.3.fx(lc50 = gaf.ctrl.drc.mir$coefficients[2], slp = gaf.ctrl.drc.mir$coefficients[1], t =time), lty = 2)
  
  for(i in 1:length(unique(mir.but$conc))){
    points(mir.but$time_hrs[mir.but$conc==unique(mir.but$conc)[i]], 
           mir.but$surv[mir.but$conc==unique(mir.but$conc)[i]], pch=17,
           col = i+1)
    lines(time, L.3.fx(lc50 = gaf.but.drc.mir$coefficients[i+5], slp = gaf.but.drc.mir$coefficients[i], t =time), 
          lty = 2, col = i+1)
  }
  
  legend('topright', title = 'Butralin (ppb)', legend = but.vals, pch = c(16,rep(17,5)),
         col = c(1:6), cex=0.7, bty = 'n')
  
dev.off()  

#Fitted d-r function parameters ###########
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_fitted_d-r_parameters_and_functions_miracidia_butralin.png")

plot(but.mir$butralin, but.mir$e, pch = 16, xlab = 'butralin (ppb)', 
     ylab = 'L.2 Parameters (miracidia)', ylim = c(0,6))

  points(but.mir$butralin, but.mir$b, pch = 17, col=2)

  for(i in 1:length(but.mir$butralin)){
    segments(x0 = but.mir$butralin[i], y0 = but.mir$e[i] + but.mir$e.se[i],
             x1 = but.mir$butralin[i], y1 = but.mir$e[i] - but.mir$e.se[i])
    segments(x0 = but.mir$butralin[i], y0 = but.mir$b[i] + but.mir$b.se[i],
             x1 = but.mir$butralin[i], y1 = but.mir$b[i] - but.mir$b.se[i], col=2)
  } 
  
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred, simplify = T)[3,], lty = 3)

    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred2, simplify = T)[3,], lty = 3, col=3)

    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,9000,100)/1000, sapply(seq(0,9000,100)/1000, but.mir.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topleft', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

dev.off()


#plot to test function ########
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_function_simulate_miracidia_butralin.png")

plot(but.mir$butralin, gaf.but.aucs.mir/gaf.but.aucs.mir[1],
     xlab = 'Butralin (ppm)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))

  set.seed(43093)

  points(seq(0, 9000,100)/1000, sapply(seq(0, 9000, 100), piM.ghaf_butr.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0, 9000, 100)/1000, sapply(seq(0, 9000, 100), piM.ghaf_butr.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)
  
  legend('bottomleft', pch = c(16,5,5), col = c(1,2,4), legend = c("observed", "log-linear lc50 fit", "linear lc50 fit"),
         cex = 0.7, bty = 'n')
  
dev.off()  

#glyphosate ###############
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_raw_miracidia_mortality_data_glyphosate.png")

plot(mir.ctrl$time_hrs, mir.ctrl$surv, ylim = c(0,1), xlim = c(0,24), pch=16, 
     xlab = 'time(hrs)', ylab = 'prop surviving',
     main = 'Abdel-Ghaffar 2016: Glyphosate toxicity to miracidia')
  lines(time, L.3.fx(lc50 = gaf.ctrl.drc.mir$coefficients[2], slp = gaf.ctrl.drc.mir$coefficients[1], t =time), lty = 2)
  
  for(i in 1:length(unique(mir.gly$conc))){
    points(mir.gly$time_hrs[mir.gly$conc==unique(mir.gly$conc)[i]], 
           mir.gly$surv[mir.gly$conc==unique(mir.gly$conc)[i]], pch=17,
           col = i+1)
    lines(time, L.3.fx(lc50 = gaf.gly.drc.mir$coefficients[i+5], slp = gaf.gly.drc.mir$coefficients[i], t =time), 
          lty = 2, col = i+1)
  }

  legend('topright', title = 'Glyphosate (ppm)', legend = gly.vals, pch = c(16,rep(17,5)),
         col = c(1:6), cex=0.7, bty = 'n')
  
dev.off()  
  
#Get estimate of miracidia-hrs for each concentration    
  gaf.gly.aucs.mir = as.numeric()
  gaf.gly.aucs.mir[1] = integrate(f = L.3.fx, lc50 = gaf.ctrl.drc.mir$coefficients[2], slp = gaf.ctrl.drc.mir$coefficients[1], 
                                  lower=0, upper=24)[1]$value  
  
  for(j in 1:length(unique(mir.gly$conc))){
    gaf.gly.aucs.mir[j+1] = integrate(f = L.3.fx, lc50 = gaf.gly.drc.mir$coefficients[j+5], slp = gaf.gly.drc.mir$coefficients[j],
                                      lower=0, upper=24)[1]$value  
  }

#Fitted d-r function parameters ###########
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_fitted_d-r_parameters_and_functions_miracidia_glyphosate.png")
  
plot(gly.mir$glyphosate, gly.mir$e, pch = 16, xlab = 'glyphosate (ppm)', 
     ylab = 'LL.2 Parameters (miracidia)', ylim = c(0,6))

  points(gly.mir$glyphosate, gly.mir$b, pch = 17, col=2)

for(i in 1:length(gly.mir$glyphosate)){
  segments(x0 = gly.mir$glyphosate[i], y0 = gly.mir$e[i] + gly.mir$e.se[i],
           x1 = gly.mir$glyphosate[i], y1 = gly.mir$e[i] - gly.mir$e.se[i])
  segments(x0 = gly.mir$glyphosate[i], y0 = gly.mir$b[i] + gly.mir$b.se[i],
           x1 = gly.mir$glyphosate[i], y1 = gly.mir$b[i] - gly.mir$b.se[i], col=2)
} 

    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred, simplify = T)[3,], lty = 3)

    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred2, simplify = T)[3,], lty = 3, col=3)

    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,9000,100)*3/1000, sapply(seq(0,9000,100)*3/1000, gly.mir.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topleft', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

dev.off()

auc.gly.mir.lin0 = function(He){
    e0 = as.numeric(predict(gly.mir.lm.e, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(gly.mir.lm.b, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
    
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
    auc0
  }
  
piM.ghaf_gly.lin_unc = function(He){
  #if(He == 0) piM = 1 else{
    Heu = He/1000
    e = as.numeric(predict(gly.mir.lm.e, newdata = data.frame(glyphosate = Heu), se.fit = TRUE)[1:2])
    b = as.numeric(predict(gly.mir.lm.b, newdata = data.frame(glyphosate = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    #print(e.use)
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    #print(b.use)
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piM = auc/auc.gly.mir.lin0(0) #gaf.gly.aucs.mir[1]
  #  if(piM > 1) piM = 1
  #}
  return(piM)
}  #function to estimate AUC 

auc.gly.mir.exp0 = function(He){
  e0 = as.numeric(predict(gly.mir.lm.e2, newdata = data.frame(logglyphosate = log(He+1)), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(gly.mir.lm.b, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
  auc0
}


piM.ghaf_gly.exp_unc = function(He){
  #if(He == 0) piM = 1 else{
    Heu = He/1000
    e = as.numeric(predict(gly.mir.lm.e2, newdata = data.frame(logglyphosate = log(Heu+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(gly.mir.lm.b, newdata = data.frame(glyphosate = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    #print(e.use)
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, e[1], e[2])
    #print(b.use)
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piM = auc/auc.gly.mir.exp0(0) #gaf.gly.aucs.mir[1]
  #  if(piM > 1) piM = 1
  #}
  return(piM)
}  #function to estimate AUC 

keep.gaf.gly.mir = c('piM.ghaf_gly.lin_unc', 'gly.mir.lm.e', 'gly.mir.lm.b', 'L.3.fx', 'gaf.gly.aucs.mir',
                     'piM.ghaf_gly.exp_unc', 'gly.mir.lm.e2', 'auc.gly.mir.lin0', 'auc.gly.mir.exp0')

#plot to test function ########
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_function_simulate_miracidia_glyphosate.png")

plot(gly.mir$glyphosate, gaf.gly.aucs.mir/gaf.gly.aucs.mir[1],
     xlab = 'Glyphosate (ppb)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))

  set.seed(43093)

  points(seq(0,25500,100)/1000, sapply(seq(0,25500,100), piM.ghaf_gly.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,25500,100)/1000, sapply(seq(0,25500,100), piM.ghaf_gly.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)

  legend('bottomleft', pch = c(16,5,5), col = c(1,2,4), legend = c("observed", "log-linear lc50 fit", "linear lc50 fit"),
         cex = 0.7, bty = 'n')
  
dev.off()
#pendimethalin ###############
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_raw_miracidia_mortality_data_pendimethalin.png")

plot(mir.ctrl$time_hrs, mir.ctrl$surv, ylim = c(0,1), xlim = c(0,24), pch=16, 
     xlab = 'time(hrs)', ylab = 'prop surviving',
     main = 'Abdel-Ghaffar 2016: Pendimethalin toxicity to miracidia')
  lines(time, L.3.fx(lc50 = gaf.ctrl.drc.mir$coefficients[2], slp = gaf.ctrl.drc.mir$coefficients[1], t =time), lty = 2)

  for(i in 1:length(unique(mir.pen$conc))){
    points(mir.pen$time_hrs[mir.pen$conc==unique(mir.pen$conc)[i]], 
           mir.pen$surv[mir.pen$conc==unique(mir.pen$conc)[i]], pch=17,
           col = i+1)
    lines(time, L.3.fx(lc50 = gaf.pen.drc.mir$coefficients[i+5], slp = gaf.pen.drc.mir$coefficients[i], t =time), 
          lty = 2, col = i+1)
  }

  legend('topright', title = 'pendimethalin (ppm)', legend = pen.vals, pch = c(16,rep(17,5)),
         col = c(1:6), cex=0.7, bty = 'n')

dev.off()  
#Get estimate of miracidia-hrs for each concentration    
gaf.pen.aucs.mir = as.numeric()
  gaf.pen.aucs.mir[1] = integrate(f = L.3.fx, lc50 = gaf.ctrl.drc.mir$coefficients[2], slp = gaf.ctrl.drc.mir$coefficients[1], 
                                  lower=0, upper=24)[1]$value  
  
  for(j in 1:length(unique(mir.pen$conc))){
    gaf.pen.aucs.mir[j+1] = integrate(f = L.3.fx, lc50 = gaf.pen.drc.mir$coefficients[j+5], slp = gaf.pen.drc.mir$coefficients[j],
                                      lower=0, upper=24)[1]$value  
  }
  
#Fitted d-r function parameters ###########
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_fitted_d-r_parameters_and_functions_miracidia_pendimethalin.png")
  
plot(pen.mir$pendimethalin, pen.mir$e, pch = 16, xlab = 'Pendimethalin (ppb)', 
     ylab = 'L.2 Parameters (miracidia)', ylim = c(0,5))

  points(pen.mir$pendimethalin, pen.mir$b, pch = 17, col=2)

  for(i in 1:length(pen.mir$pendimethalin)){
    segments(x0 = pen.mir$pendimethalin[i], y0 = pen.mir$e[i] + pen.mir$e.se[i],
             x1 = pen.mir$pendimethalin[i], y1 = pen.mir$e[i] - pen.mir$e.se[i])
    segments(x0 = pen.mir$pendimethalin[i], y0 = pen.mir$b[i] + pen.mir$b.se[i],
             x1 = pen.mir$pendimethalin[i], y1 = pen.mir$b[i] - pen.mir$b.se[i], col=2)
  } 

    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred, simplify = T)[1,], lty = 2)
    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred, simplify = T)[2,], lty = 3)
    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred, simplify = T)[3,], lty = 3)

    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred2, simplify = T)[1,], lty = 2, col=3)
    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred2, simplify = T)[2,], lty = 3, col=3)
    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred2, simplify = T)[3,], lty = 3, col=3)

    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred.b, simplify = T)[1,], lty = 2, col = 2)
    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred.b, simplify = T)[2,], lty = 3, col = 2)
    lines(seq(0,4000,25)/1000, sapply(seq(0,4000,25)/1000, pen.mir.pred.b, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topleft', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

dev.off()  

#plot to test function ########
png("Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_function_simulate_miracidia_pendimethalin.png")

plot(pen.mir$pendimethalin, gaf.pen.aucs.mir/gaf.pen.aucs.mir[1],
     xlab = 'pendimethalin (ppm)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))
  
set.seed(43093)

  points(seq(0,4000,25)/1000, sapply(seq(0,4000,25), piM.ghaf_pen.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0,4000,25)/1000, sapply(seq(0,4000,25), piM.ghaf_pen.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)
  
  legend('bottomleft', pch = c(16,5,5), col = c(1,2,4), legend = c("observed", "log-linear lc50 fit", "linear lc50 fit"),
         cex = 0.7, bty = 'n')
  
dev.off()