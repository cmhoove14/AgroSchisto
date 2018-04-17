#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
#Data extraction and model fitting to Tantawy 2002 data
source("Agrochemical_Review/Response_Fxs/tantawy2002_butachlor_fpb_miracidia_fit.R")

#butachlor toxicity to miracidia  ###############
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Miracidia/tantawy2002_butachlor_piM_data.png")

plot(tant.mir.but.dat$time_hrs[tant.mir.but.dat$conc==0], tant.mir.but.dat$surv[tant.mir.but.dat$conc==0]/100, 
     pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tant.piM.but, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in 2:length(unique(tant.mir.but.dat$conc))){
    points(tant.mir.but.dat$time_hrs[tant.mir.but.dat$conc==unique(tant.mir.but.dat$conc)[i]], 
           tant.mir.but.dat$surv[tant.mir.but.dat$conc==unique(tant.mir.but.dat$conc)[i]]/100, pch=16,
           col = i)
    lines(time, predict(tant.piM.but, 
                        data.frame(time_hrs=time, conc = unique(tant.mir.but.dat$conc)[i])),
          lty = 2, col = i)
  }

title('butachlor toxicity to miracidia')
legend('topright', legend = c('control', .650,1.500,4.500,6.500), title = 'butachlor (ppm)',
       pch = c(17,16,16,16,16), col = c(1:5), cex=0.8, bty = 'n')  

dev.off()

#Compile LL.2 data for functional responses #############
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Miracidia/tantawy2002_butachlor_piM_parameter_models.png")

plot(mir.but$but, mir.but$e, pch = 16, xlab = 'butachlor (ppm)', ylab = 'LL.2 Parameters (miracidia)',
     ylim = c(0,20))

  points(mir.but$but, mir.but$b, pch = 17, col=2)

  for(i in 1:length(mir.but$but)){
    segments(x0 = mir.but$but[i], y0 = mir.but$e[i] + mir.but$e.se[i],
             x1 = mir.but$but[i], y1 = mir.but$e[i] - mir.but$e.se[i])
    segments(x0 = mir.but$but[i], y0 = mir.but$b[i] + mir.but$b.se[i],
             x1 = mir.but$but[i], y1 = mir.but$b[i] - mir.but$b.se[i], col=2)
  } 
  
#fit models to LL.2 parameters across concentration ########
    lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred.mir, simplify = T)[1,], lty = 2)
    lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred.mir, simplify = T)[2,], lty = 3)
    lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred.mir, simplify = T)[3,], lty = 3)

  lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred2.mir, simplify = T)[1,], lty = 2, col=3)
  lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred2.mir, simplify = T)[2,], lty = 3, col=3)
  lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred2.mir, simplify = T)[3,], lty = 3, col=3)

  lines(seq(0,7,.07), sapply(seq(0,7,.07), bl.pred.mir, simplify = T)[1,], lty = 2, col = 2)
  lines(seq(0,7,.07), sapply(seq(0,7,.07), bl.pred.mir, simplify = T)[2,], lty = 3, col = 2)
  lines(seq(0,7,.07), sapply(seq(0,7,.07), bl.pred.mir, simplify = T)[3,], lty = 3, col = 2)
  
  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

dev.off()
#plot sample outputs compared to observed points ##########
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Miracidia/tantawy2002_butachlor_piM_data_sim.png")

plot(mir.but$but*1000, tant.but.piM.aucs/tant.but.piM.aucs[1],
     xlab = 'Butachlor (ppb)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))

  set.seed(43093)
  
  points(seq(0, 6500,25), sapply(seq(0, 6500,25), piM.tant02_but.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0, 6500,25), sapply(seq(0, 6500,25), piM.tant02_but.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)
  
  legend("bottomleft", legend = c("observed", "linear LC50", "exponential LC50"),
         pch = c(16,5,5), col = c(1,4,2), cex = 0.6, bty = "n")

  dev.off()    

#fluazifop-p-butyl toxicity to miracidia ###########
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Miracidia/tantawy2002_fpb_piM_data.png")

plot(tant.mir.fpb.dat$time_hrs[tant.mir.fpb.dat$conc==0], tant.mir.fpb.dat$surv[tant.mir.fpb.dat$conc==0]/100, 
     pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tant.piM.fpb, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in 2:length(unique(tant.mir.fpb.dat$conc))){
    points(tant.mir.fpb.dat$time_hrs[tant.mir.fpb.dat$conc==unique(tant.mir.fpb.dat$conc)[i]], 
           tant.mir.fpb.dat$surv[tant.mir.fpb.dat$conc==unique(tant.mir.fpb.dat$conc)[i]]/100, pch=16,
           col = i)
    lines(time, predict(tant.piM.fpb, 
                        data.frame(time_hrs=time, conc = unique(tant.mir.fpb.dat$conc)[i])),
          lty = 2, col = i)
  }

  title('fluazifop-p-butyl toxicity to miracidia')
  legend('topright', legend = c('control', 1.760,4.500,9.000,17.600), 
         pch = c(17,16,16,16,16), col = c(1:5), cex=0.8, bty = 'n')  

dev.off()  

#Compile LL.2 data for functional responses #############
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Miracidia/tantawy2002_fpb_piM_parameter_models.png")
  
  plot(mir.fpb$fpb, mir.fpb$e, pch = 16, xlab = 'butachlor (ppm)', ylab = 'LL.2 Parameters (miracidia',
       ylim = c(0,20))
  
  points(mir.fpb$fpb, mir.fpb$b, pch = 17, col=2)
  
  for(i in 1:length(mir.fpb$fpb)){
    segments(x0 = mir.fpb$fpb[i], y0 = mir.fpb$e[i] + mir.fpb$e.se[i],
             x1 = mir.fpb$fpb[i], y1 = mir.fpb$e[i] - mir.fpb$e.se[i])
    segments(x0 = mir.fpb$fpb[i], y0 = mir.fpb$b[i] + mir.fpb$b.se[i],
             x1 = mir.fpb$fpb[i], y1 = mir.fpb$b[i] - mir.fpb$b.se[i], col=2)
  } 
  
#fit models to LL.2 parameters across concentration ########
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.mir.fpb, simplify = T)[1,], lty = 2)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.mir.fpb, simplify = T)[2,], lty = 3)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.mir.fpb, simplify = T)[3,], lty = 3)

  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred2.mir.fpb, simplify = T)[1,], lty = 2, col=3)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred2.mir.fpb, simplify = T)[2,], lty = 3, col=3)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred2.mir.fpb, simplify = T)[3,], lty = 3, col=3)

  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), bl.pred.mir.fpb, simplify = T)[1,], lty = 2, col = 2)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), bl.pred.mir.fpb, simplify = T)[2,], lty = 3, col = 2)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), bl.pred.mir.fpb, simplify = T)[3,], lty = 3, col = 2)

  legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
         legend = c('LC50 - Linear',
                    'LC50 - Exponential',
                    'Slp - linear',
                    '95% CI'), cex = 0.7)  
  legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
         cex = 0.8, bty = 'n', title = 'Observed points')

dev.off()

#plot sample output compared to observed AUCs ##########
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Miracidia/tantawy2002_fpb_piM_data_sim.png")

plot(mir.fpb$fpb, tant.fpb.piM.aucs/tant.fpb.piM.aucs[1],
     xlab = 'Fluazifop-p-butyl (ppb)', ylab = 'relative AUC (miracidia-hours)',
     pch = 16, ylim = c(0,1))
  
  set.seed(43093)

  points(seq(0, 18000,50)/1000, sapply(seq(0, 18000,50), piM.tant02_fpb.lin_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)
  points(seq(0, 18000,50)/1000, sapply(seq(0, 18000,50), piM.tant02_fpb.exp_unc, simplify = T),
         pch = 5, col = 2, cex = 0.5)
  
  legend("bottomleft", legend = c("observed", "linear LC50", "exponential LC50"),
         pch = c(16,5,5), col = c(1,4,2), cex = 0.6, bty = "n")

dev.off()