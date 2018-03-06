#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/griggs08_atrazine_metolachlor_cercariae_fit.R")

#Plot control survival curve, estimate time-dependent die-off function, and get auc from 0-24 hours ########
png("Agrochemical_Review/Response_Fxs/Plots/Griggs2008/Griggs_2008_raw_cercariae_mortality_data_atrazine.png")

plot(x = cerc.g0$time_hrs[cerc.g0$conc == 0], y = cerc.g0$prop_surv[cerc.g0$conc == 0],
       xlab = 'time (hrs)', ylab = 'prop alive', pch = 16, cex = 1.2, xlim = c(0,25), ylim = c(0,1))
    lines(time, predict(grg.mod, 
                        data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in unique(cerc.g0$conc)[c(2:3)]){
    points(cerc.g0$time_hrs[cerc.g0$conc == i], cerc.g0$prop_surv[cerc.g0$conc == i], 
           pch = 17, col = i+1)
    lines(time, predict(grg.mod, 
                        data.frame(time_hrs=time, conc = i)), lty = 2, col = i+1)
  } 

  title(main='Griggs2008 Atrazine-Cercarial mortality (E.trivolvis)')
    legend('bottomleft', legend = c('control', '15ppb', '100ppb'), 
           pch = c(16,17,17), col = c(1,8,5), cex=0.8, bty = 'n')
    
dev.off()    

#Fitted d-r function parameters ###########
modgdf= data.frame(atr = c(0:300),
                     logatr = log(c(0:300)+1),
                     pred.e = 0,
                     pred.e.se = 0,
                     pred.e2 = 0,
                     pred.e2.se = 0,
                     pred.b = 0,
                     pred.b.se = 0)
  
  modgdf[,3:4] = predict(eg.mod, newdata = modgdf, se.fit = TRUE)[1:2]
  modgdf[,5:6] = predict(eg.mod2, newdata = modgdf, se.fit = TRUE)[1:2]
  modgdf[,7:8] = predict(bg.mod, newdata = modgdf, se.fit = TRUE)[1:2]

png("Agrochemical_Review/Response_Fxs/Plots/Griggs2008/Griggs_2008_fitted_d-r_parameters_and_functions_cercariae_atrazine.png")
  
plot(grgc.df$atr, grgc.df$e, pch = 16, xlab = 'atrazine (ppb)', ylab = 'LL.2 Parameters',
     ylim = c(0,17), xlim = c(0,300))
  points(grgc.df$atr+3, grgc.df$b, pch = 17, col=2)
    for(i in 1:length(grgc.df$atr)){
      segments(x0 = grgc.df$atr[i], y0 = grgc.df$e[i] + grgc.df$e.se[i],
               x1 = grgc.df$atr[i], y1 = grgc.df$e[i] - grgc.df$e.se[i])
      segments(x0 = grgc.df$atr[i]+3, y0 = grgc.df$b[i] + grgc.df$b.se[i],
               x1 = grgc.df$atr[i]+3, y1 = grgc.df$b[i] - grgc.df$b.se[i], col=2)
    }
  
lines(modgdf$atr, modgdf$pred.e, lty = 2)
  lines(modgdf$atr, modgdf$pred.e + 1.96*modgdf$pred.e.se, lty = 3)
  lines(modgdf$atr, modgdf$pred.e - 1.96*modgdf$pred.e.se, lty = 3)

lines(modgdf$atr, modgdf$pred.e2, lty = 2, col=3)
  lines(modgdf$atr, modgdf$pred.e2 + 1.96*modgdf$pred.e2.se, lty = 3, col=3)
  lines(modgdf$atr, modgdf$pred.e2 - 1.96*modgdf$pred.e2.se, lty = 3, col=3)  
  
lines(modgdf$atr, modgdf$pred.b, lty = 2, col=2)
  lines(modgdf$atr, modgdf$pred.b + 1.96*modgdf$pred.b.se, lty = 3, col=2)
  lines(modgdf$atr, modgdf$pred.b - 1.96*modgdf$pred.b.se, lty = 3, col=2)
  
legend('topright', legend = c('slp', 'lc50'), pch = c(17, 16), col=c(2,1), cex = 0.7, bty='n')  
title('Griggs 08 cercarial survival parameters')  

dev.off()

#plot to test function ########
png("Agrochemical_Review/Response_Fxs/Plots/Griggs2008/Griggs_2008_function_simulate_cercariae_atrazine.png")

plot(c(0,15,100), grg08_atr_aucs/grg08_atr_aucs[1], ylim = c(0,1),pch = 16, cex = 1.2,
     xlab = 'atrazine (ppb)', ylab = expression(paste(pi[C], 'estimate')))
  
  points(c(0:100), sapply(c(0:100), piC.grg08_atr_unc), pch = 5, col = 2, cex = 0.5)
  points(c(0:100), sapply(c(0:100), piC.grg08_atr_unc2), pch = 5, cex = 0.5, col = 4)
  legend('bottomleft', legend = c('linear', 'exponential'), pch = 5, col = c(2,4),
         title = 'Function fit to lc50 parameter', cex = 0.7, bty='n')
  
dev.off()