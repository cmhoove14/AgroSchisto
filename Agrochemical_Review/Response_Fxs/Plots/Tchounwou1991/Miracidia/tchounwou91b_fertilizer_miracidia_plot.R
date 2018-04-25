#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/tchounwou91b_fertilizer_miracidia_fit.R")

#Ammonium Sulphate concentration ############
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_ammonium_phosphate_miracidia_time_data.png")

plot(mir.amm$time_hrs[mir.amm$conc==0], mir.amm$alive[mir.amm$conc==0]/mir.amm$total[mir.amm$conc==0], 
     pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tch91.piM.amm, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in 2:length(unique(mir.amm$conc))){
    points(mir.amm$time_hrs[mir.amm$conc==unique(mir.amm$conc)[i]], 
           mir.amm$alive[mir.amm$conc==unique(mir.amm$conc)[i]] /
           mir.amm$total[mir.amm$conc==unique(mir.amm$conc)[i]], pch=16,
           col = i)
    lines(time, predict(tch91.piM.amm, 
                        data.frame(time_hrs=time, conc = unique(mir.amm$conc)[i])),
          lty = 2, col = i)
  }

title('ammonium phosphate toxicity to miracidia')
legend('topright', legend = c('control', unique(mir.amm$conc)[-1]), 
       pch = c(17,16,16,16,16,16), col = c(1:6), cex=0.8, bty = 'n')  

dev.off()

#Create data frame with parameter values and ammonium sulphate concentrations #######################
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_ammonium_phosphate_miracidia_d-r_pars.png")

plot(miramm.df$amm, miramm.df$e, pch = 16, xlab = 'amm. sulphate (ppm)', ylab = 'LL.2 Parameters',
     ylim = c(0, 7))
  points(miramm.df$amm, miramm.df$b, pch = 17, col=2)
    for(i in 1:length(miramm.df$amm)){
      segments(x0 = miramm.df$amm[i], y0 = miramm.df$e[i] + miramm.df$e.se[i],
               x1 = miramm.df$amm[i], y1 = miramm.df$e[i] - miramm.df$e.se[i])
      segments(x0 = miramm.df$amm[i], y0 = miramm.df$b[i] + miramm.df$b.se[i],
               x1 = miramm.df$amm[i], y1 = miramm.df$b[i] - miramm.df$b.se[i], col=2)
    }
#linear fit
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.e.amm.lin.pred)[1,],
          lty = 2)
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.e.amm.lin.pred)[2,],
          lty = 3)
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.e.amm.lin.pred)[3,],
          lty = 3)
    
#exponential fit
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.e.amm.exp.pred)[1,],
          lty = 2, col = 3)
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.e.amm.exp.pred)[2,],
          lty = 3, col = 3)
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.e.amm.exp.pred)[3,],
          lty = 3, col = 3)    
    
AIC(tch91.e.amm.lin, tch91.e.amm.exp) #exponential is a better fit
    
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.b.amm.pred)[1,],
          lty = 2, col = 2)
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.b.amm.pred)[2,],
          lty = 3, col = 2)
    lines(seq(0,max(miramm.df$amm)+1,50), sapply(seq(0,max(miramm.df$amm)+1,50), tch91.b.amm.pred)[3,],
          lty = 3, col = 2)

  dev.off()  
#Qualitative model validation ###############
#Regenerate plot of observed data
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_ammonium_phosphate_miracidia_validate.png")

plot(mir.amm$time_hrs[mir.amm$conc==0], mir.amm$surv[mir.amm$conc==0], pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 2:length(unique(mir.amm$conc))){
    points(mir.amm$time_hrs[mir.amm$conc==unique(mir.amm$conc)[i]], 
           mir.amm$surv[mir.amm$conc==unique(mir.amm$conc)[i]], pch=16,
           col = i)
  }
  title('ammonium phosphate toxicity to miracidia')
  legend('topright', legend = c('control', unique(mir.amm$conc)[-1]), 
         pch = c(17,16,16,16,16,16), col = c(1:6), cex=0.8, bty = 'n')

#plot model predictions
set.seed(43093)
for(i in c(unique(mir.amm$conc)*1000)){
  c = i/1000000 + 1
  print(c)
  replicate(10, predm.amm.plot(In = i, clr = c))
}

dev.off()  
#plot AUC predictions
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_ammonium_phosphate_miracidia_piM_validate.png")
plot(mir.auc.amm[,1], mir.auc.amm[,2], pch = 16, cex = 1.2)
set.seed(43093)
  points(seq(0,max(miramm.df$amm)+1,25), sapply(seq(0,max(miramm.df$amm)*1000,25000) , piM.tch91_amm_lin),
         pch = 17, cex = 0.5, col=4)
  points(seq(0,max(miramm.df$amm)+1,25), sapply(seq(0,max(miramm.df$amm)*1000,25000) , piM.tch91_amm_unc),
         pch = 17, cex = 0.5, col=2)
  
    legend('bottomleft', pch = c(16,17,17), col = c(1,4,2), legend = c("observed", "simulated-linear", "simulated-exponential"),
         cex = 0.8, bty = 'n')

dev.off() 

#Urea concentration ############
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_urea_miracidia_time_data.png")

plot(mir.ure$time_hrs[mir.ure$conc==0], mir.ure$alive[mir.ure$conc==0]/mir.ure$total[mir.ure$conc==0], 
     pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tch91.piM.ure, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
    for(i in 2:length(unique(mir.ure$conc))){
      points(mir.ure$time_hrs[mir.ure$conc==unique(mir.ure$conc)[i]], 
             mir.ure$alive[mir.ure$conc==unique(mir.ure$conc)[i]] /
               mir.ure$total[mir.ure$conc==unique(mir.ure$conc)[i]], pch=16,
             col = i)
      lines(time, predict(tch91.piM.ure, 
                          data.frame(time_hrs=time, conc = unique(mir.ure$conc)[i])),
            lty = 2, col = i)
    }

title('urea toxicity to miracidia')
legend('topright', legend = c('control', unique(mir.ure$conc)[-1]), 
       pch = c(17,16,16,16,16,16), col = c(1:6), cex=0.8, bty = 'n')  

dev.off()


#Create data frame with parameter values and urea concentrations #######################
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_urea_miracidia_piM_dose_pars.png")
plot(mirure.df$ure, mirure.df$e, pch = 16, xlab = 'urea (ppb)', ylab = 'LL.2 Parameters',
     ylim = c(0, 7))
  points(mirure.df$ure, mirure.df$b, pch = 17, col=2)
  for(i in 1:length(mirure.df$ure)){
    segments(x0 = mirure.df$ure[i], y0 = mirure.df$e[i] + mirure.df$e.se[i],
             x1 = mirure.df$ure[i], y1 = mirure.df$e[i] - mirure.df$e.se[i])
    segments(x0 = mirure.df$ure[i], y0 = mirure.df$b[i] + mirure.df$b.se[i],
             x1 = mirure.df$ure[i], y1 = mirure.df$b[i] - mirure.df$b.se[i], col=2)
  }
#linear fit
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.e.ure.lin.pred)[1,],
          lty = 2)
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.e.ure.lin.pred)[2,],
          lty = 3)
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.e.ure.lin.pred)[3,],
          lty = 3)
  
#exponential fit
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.e.ure.exp.pred)[1,],
          lty = 2, col = 3)
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.e.ure.exp.pred)[2,],
          lty = 3, col = 3)
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.e.ure.exp.pred)[3,],
          lty = 3, col = 3)    

AIC(tch91.e.ure.lin, tch91.e.ure.exp) #Linear is a better fit

    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.b.ure.pred)[1,],
          lty = 2, col = 2)
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.b.ure.pred)[2,],
          lty = 3, col = 2)
    lines(seq(0,max(mirure.df$ure)+1,50), sapply(seq(0,max(mirure.df$ure)+1,50), tch91.b.ure.pred)[3,],
          lty = 3, col = 2)

dev.off()

#Qualitative model validation ###############
#Regenerate plot of observed data
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_urea_miracidia_validate.png")

plot(mir.ure$time_hrs[mir.ure$conc==0], mir.ure$surv[mir.ure$conc==0], pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
for(i in 2:length(unique(mir.ure$conc))){
  points(mir.ure$time_hrs[mir.ure$conc==unique(mir.ure$conc)[i]], 
         mir.ure$surv[mir.ure$conc==unique(mir.ure$conc)[i]], pch=16,
         col = i)
}
title('urea toxicity to miracidia')
legend('topright', legend = c('control', unique(mir.ure$conc)[-1]), 
       pch = c(17,16,16,16,16,16), col = c(1:6), cex=0.8, bty = 'n')


#plot model predictions
set.seed(43093)
for(i in c(unique(mir.ure$conc)*1000)){
  c = i/2000000 + 1
  print(c)
  replicate(10, predm.ure.plot(In = i, clr = c))
}
dev.off()
#plot AUC predictions
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_urea_miracidia_piM_validate.png")

plot(mir.auc.ure[,1], mir.auc.ure[,2], pch = 16, cex = 1.2, ylab = "piM", xlab = "Urea concantration (ppm)")
  set.seed(43093)
  points(seq(0,max(mirure.df$ure)+1,25), sapply(seq(0,max(mirure.df$ure)*1000,25000) , piM.tch91_ure_unc),
         pch = 17, cex = 0.5, col=4)
  legend('bottomleft', pch = c(16,17), col = c(1,4), legend = c("observed", "simulated"),
         cex = 0.8, bty = 'n')

dev.off()