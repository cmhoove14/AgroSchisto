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
source("Agrochemical_Review/Response_Fxs/tantawy2002_butachlor_fpb_cercariae_fit.R")

#butachlor toxicity to cercariae  ###############
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Cercariae/tantawy2002_butachlor_piC_data.png")

plot(tant.cerc.but.dat$time_hrs[tant.cerc.but.dat$conc==0], tant.cerc.but.dat$surv[tant.cerc.but.dat$conc==0]/100, 
       pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tant.piC.but, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in 2:length(unique(tant.cerc.but.dat$conc))){
    points(tant.cerc.but.dat$time_hrs[tant.cerc.but.dat$conc==unique(tant.cerc.but.dat$conc)[i]], 
           tant.cerc.but.dat$surv[tant.cerc.but.dat$conc==unique(tant.cerc.but.dat$conc)[i]]/100, pch=16,
           col = i)
    lines(time, predict(tant.piC.but, 
                        data.frame(time_hrs=time, conc = unique(tant.cerc.but.dat$conc)[i])),
          lty = 2, col = i)
  }
  
  title('butachlor toxicity to cercariae')
  legend('topright', legend = c('control', .650,1.500,4.500,6.500), title = 'butachlor (ppm)',
         pch = c(17,16,16,16,16), col = c(1:5), cex=0.8, bty = 'n')  

dev.off()  

#compile butachlor data for function ###############
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Cercariae/tantawy2002_butachlor_piC_parameter_models.png")
    
plot(cerc.but$but, cerc.but$e, pch = 16, xlab = 'butachlor (ppm)', ylab = 'LL.2 Parameters',
     ylim = c(0,20))

  points(cerc.but$but, cerc.but$b, pch = 17, col=2)
  
    for(i in 1:length(cerc.but$but)){
      segments(x0 = cerc.but$but[i], y0 = cerc.but$e[i] + cerc.but$e.se[i],
               x1 = cerc.but$but[i], y1 = cerc.but$e[i] - cerc.but$e.se[i])
      segments(x0 = cerc.but$but[i], y0 = cerc.but$b[i] + cerc.but$b.se[i],
               x1 = cerc.but$but[i], y1 = cerc.but$b[i] - cerc.but$b.se[i], col=2)
    }    
  
#fit models to L.2 parameters across concentration ########
    lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred)[1,], lty = 2)
    lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred)[2,], lty = 3)
    lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred)[3,], lty = 3)

  lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred2)[1,], lty = 2, col=3)
  lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred2)[2,], lty = 3, col=3)
  lines(seq(0,7,.07), sapply(seq(0,7,.07), el.pred2)[3,], lty = 3, col=3)
    
    lines(seq(0,7,.07), sapply(seq(0,7,.07), bl.pred)[1,], lty = 2, col = 2)
    lines(seq(0,7,.07), sapply(seq(0,7,.07), bl.pred)[2,], lty = 3, col = 2)
    lines(seq(0,7,.07), sapply(seq(0,7,.07), bl.pred)[3,], lty = 3, col = 2)
    
    legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
           legend = c('LC50 - Linear',
                      'LC50 - Exponential',
                      'Slp - linear',
                      '95% CI'), cex = 0.7)  
    legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
           cex = 0.7, bty = 'n', title = 'Observed points')

dev.off()

#plot sample output compared to observed AUCs ##########
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Cercariae/tantawy2002_butachlor_piC_data_sim.png")
  plot(cerc.but$but, tant.but.pic.aucs / tant.but.pic.aucs[1],
        xlab = 'Butachlor (ppm)', ylab = 'relative AUC (cercariae-hours)',
       pch = 16, ylim = c(0,1))

  set.seed(43093)
  
    points(seq(0, 6500,50)/1000, sapply(seq(0, 6500,50), piC.tant02_but.lin_unc),
           pch = 5, col = 4, cex = 0.5)
    points(seq(0, 6500,50)/1000, sapply(seq(0, 6500,50), piC.tant02_but.exp_unc),
           pch = 5, col = 2, cex = 0.5)
    
  legend("bottomleft", legend = c("observed", "linear LC50", "exponential LC50"),
         pch = c(16,5,5), col = c(1,4,2), cex = 0.6, bty = "n")

dev.off()

#fluazifop-p-butyl toxicity to cercariae ###########
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Cercariae/tantawy2002_fpb_piC_data.png")

  plot(tant.cerc.fpb.dat$time_hrs[tant.cerc.fpb.dat$conc==0], tant.cerc.fpb.dat$surv[tant.cerc.fpb.dat$conc==0]/100, 
       pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tant.piC.fpb, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
    for(i in 2:length(unique(tant.cerc.fpb.dat$conc))){
      points(tant.cerc.fpb.dat$time_hrs[tant.cerc.fpb.dat$conc==unique(tant.cerc.fpb.dat$conc)[i]], 
             tant.cerc.fpb.dat$surv[tant.cerc.fpb.dat$conc==unique(tant.cerc.fpb.dat$conc)[i]]/100, pch=16,
             col = i)
      lines(time, predict(tant.piC.fpb, 
                          data.frame(time_hrs=time, conc = unique(tant.cerc.fpb.dat$conc)[i])),
            lty = 2, col = i)
    }
  
  title('fluazifop-p-butyl toxicity to cercariae')
  legend('topright', legend = c('control', 1.760,4.500,9.000,17.600), 
         pch = c(17,16,16,16,16), col = c(1:5), cex=0.8, bty = 'n')  

dev.off()  

#compile fluazifop-p-butyl data for function ###############
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Cercariae/tantawy2002_fpb_piC_parameter_models.png")
  
  plot(cerc.fpb$fpb, cerc.fpb$e, pch = 16, xlab = 'fluazifop-p-butyl (ppm)', 
       ylab = 'L.2 Parameters', ylim = c(0,20))
  
  points(cerc.fpb$fpb, cerc.fpb$b, pch = 17, col=2)
  
  for(i in 1:length(cerc.fpb$fpb)){
    segments(x0 = cerc.fpb$fpb[i], y0 = cerc.fpb$e[i] + cerc.fpb$e.se[i],
             x1 = cerc.fpb$fpb[i], y1 = cerc.fpb$e[i] - cerc.fpb$e.se[i])
    segments(x0 = cerc.fpb$fpb[i]+50, y0 = cerc.fpb$b[i] + cerc.fpb$b.se[i],
             x1 = cerc.fpb$fpb[i]+50, y1 = cerc.fpb$b[i] - cerc.fpb$b.se[i], col=2)
  }    
  
#fit models to L.2 parameters across concentration ########
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.fpb)[1,], lty = 2)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.fpb)[2,], lty = 3)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.fpb)[3,], lty = 3)
    
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.fpb2)[1,], lty = 2, col=3)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.fpb2)[2,], lty = 3, col=3)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), el.pred.fpb2)[3,], lty = 3, col=3)
    
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), bl.pred.fpb)[1,], lty = 2, col = 2)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), bl.pred.fpb)[2,], lty = 3, col = 2)
  lines(seq(0,18,0.1), sapply(seq(0,18,0.1), bl.pred.fpb)[3,], lty = 3, col = 2)
    
    legend('top', lty = c(2,2,2,3), col = c(1,3,2,1), bty = 'n', title = 'Fit models',
           legend = c('LC50 - Linear',
                      'LC50 - Exponential',
                      'Slp - linear',
                      '95% CI'), cex = 0.7)  
    legend('topright', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
           cex = 0.8, bty = 'n', title = 'Observed points')

dev.off()

#plot sample output compared to observed AUCs ##########
png("Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Cercariae/tantawy2002_fpb_piC_data_sim.png")
  plot(cerc.fpb$fpb, tant.fpb.pic.aucs/tant.fpb.pic.aucs[1],
        xlab = 'Fluazifop-p-butyl (ppb)', ylab = 'relative AUC (cercariae-hours)',
        pch = 16, ylim = c(0,1))

    set.seed(43093)
    
    points(seq(0, 18000,100)/1000, sapply(seq(0, 18000,100), piC.tant02_fpb.lin_unc),
           pch = 5, col = 4, cex = 0.5)
    points(seq(0, 18000,100)/1000, sapply(seq(0, 18000,100), piC.tant02_fpb.exp_unc),
           pch = 5, col = 2, cex = 0.5)
    
  legend("bottomleft", legend = c("observed", "linear LC50", "exponential LC50"),
         pch = c(16,5,5), col = c(1,4,2), cex = 0.6, bty = "n")

dev.off()