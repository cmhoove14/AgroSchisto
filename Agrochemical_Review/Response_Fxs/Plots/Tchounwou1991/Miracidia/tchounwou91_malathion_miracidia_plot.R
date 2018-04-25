#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/tchounwou91_malathion_miracidia_fit.R")

#Tchounwou Data plotted ############
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_malathion_miracidia_time_data.png")

plot(mir.mal$time_hrs[mir.mal$conc==0], mir.mal$alive[mir.mal$conc==0]/mir.mal$total[mir.mal$conc==0], 
     pch=17, xlab = 'time(hrs)', ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  lines(time, predict(tch91.piM.mal, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in 2:length(unique(mir.mal$conc))){
    points(mir.mal$time_hrs[mir.mal$conc==unique(mir.mal$conc)[i]], 
           mir.mal$alive[mir.mal$conc==unique(mir.mal$conc)[i]] /
           mir.mal$total[mir.mal$conc==unique(mir.mal$conc)[i]], pch=16,
           col = i)
    lines(time, predict(tch91.piM.mal, 
                        data.frame(time_hrs=time, conc = unique(mir.mal$conc)[i])),
          lty = 2, col = i)
  }

title('malathion toxicity to cercariae')
legend('topright', title = 'Mal(ppm)', legend = c('control', unique(mir.mal$conc)[-1]), 
       pch = c(17,16,16,16,16,16), col = c(1:6), cex=0.8, bty = 'n')  

dev.off()
#Create data frame with parameter values and malathion concentrations #######################
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_malathion_miracidia_dose_pars.png")

plot(mirp.df$mal, mirp.df$e, pch = 16, xlab = 'malathion (ppb)', ylab = 'LL.2 Parameters',
     ylim = c(0, 7))
  points(mirp.df$mal, mirp.df$b, pch = 17, col=2)
    for(i in 1:length(mirp.df$mal)){
      segments(x0 = mirp.df$mal[i], y0 = mirp.df$e[i] + mirp.df$e.se[i],
               x1 = mirp.df$mal[i], y1 = mirp.df$e[i] - mirp.df$e.se[i])
      segments(x0 = mirp.df$mal[i], y0 = mirp.df$b[i] + mirp.df$b.se[i],
               x1 = mirp.df$mal[i], y1 = mirp.df$b[i] - mirp.df$b.se[i], col=2)
    }

    lines(seq(0,150,0.5), sapply(seq(0,150,0.5), tch91.em.fx)[1,],
          lty = 2)
    lines(seq(0,150,0.5), sapply(seq(0,150,0.5), tch91.em.fx)[2,],
          lty = 3)
    lines(seq(0,150,0.5), sapply(seq(0,150,0.5), tch91.em.fx)[3,],
          lty = 3)

    lines(seq(0,150,0.5), sapply(seq(0,150,0.5), tch91.em.fx2)[1,],
          lty = 2, col = 3)
    lines(seq(0,150,0.5), sapply(seq(0,150,0.5), tch91.em.fx2)[2,],
          lty = 3, col = 3)
    lines(seq(0,150,0.5), sapply(seq(0,150,0.5), tch91.em.fx2)[3,],
          lty = 3, col = 3)
  
AIC(tch91.em.mod, tch91.em.mod2)  #linear is a better fit  
              
    lines(seq(0,150,0.5), sapply(seq(0,150,0.5), tch91.bm.fx)[1,],
          lty = 2, col = 2)
    lines(seq(0,150,0.5), sapply(seq(0,150,0.5), tch91.bm.fx)[2,],
          lty = 3, col = 2)
    lines(seq(0,150,0.5), sapply(seq(0,150,0.5), tch91.bm.fx)[3,],
          lty = 3, col = 2)

dev.off()    
#Qualitative model validation ###############
#Regenerate plot of observed data
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_malathion_miracidia_validate.png")

plot(mir.mal$time_hrs[mir.mal$conc==0], mir.mal$surv[mir.mal$conc==0], pch=17, xlab = 'time(hrs)',
     ylab = 'prop surviving', ylim = c(0,1), xlim = c(0,24))
  for(i in 2:length(unique(mir.mal$conc))){
    points(mir.mal$time_hrs[mir.mal$conc==unique(mir.mal$conc)[i]], 
           mir.mal$surv[mir.mal$conc==unique(mir.mal$conc)[i]], pch=16,
           col = i)
  }
  legend('topright', legend = c(0, 30, 60, 90, 120, 150), title = 'Malathion (ppm)',
         pch = c(17,rep(16,5)), col = c(1,2:6), cex = 0.8, bty = 'n')

#plot model predictions
for(i in c(0,30,60,90,120,150)*1000){
  c = i/30000 + 1
  print(c)
  replicate(10, predm.fx.plot(In = i, clr = c))
}

dev.off()
#plot model output compared to observed points
png("Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_malathion_piM_validate.png")
plot(mirp.df$mal, tch91.mal.pim.aucs/tch91.mal.pim.aucs[1], pch = 16, ylim = c(0,1),
     xlab = 'Malathion (ppm)', ylab = expression(paste(pi[M])),
     main = 'Sample Output of miracidial mortality function')
  points(seq(0,max(mirp.df$mal)*1000,500)/1000, sapply(seq(0,max(mirp.df$mal)*1000,500), piM.tch91_mal_unc), 
         pch = 5, col=4, cex = 0.5)
  
  legend('bottomleft', pch = c(16,5), col = c(1,4), legend = c("observed", "simulated"),
         cex = 0.8, bty = 'n')
  dev.off()