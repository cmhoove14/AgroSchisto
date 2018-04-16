#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/koprivnikar2006_atrazine_cercariae_fit.R")

#Cercarial mortality (E. trivolvis) from Koprivnikar 2006 ##################
png("Agrochemical_Review/Response_Fxs/Plots/Koprivnikar2006/Koprivnikar_2006_function_raw_cercariae_mortality_atrazine.png")

plot(x = kop.c$time_hrs[kop.c$conc == 0], y = kop.c$surv[kop.c$conc == 0],
     xlab = 'time (hrs)', ylab = 'prop alive', pch = 16, cex = 1.2, xlim = c(0,25), ylim = c(0,1))
  lines(time, L.3.fx(lc50 = kop.mod$coefficients[4], slp = kop.mod$coefficients[1], t = time), lty = 2)

  for(i in c(2:3)){
    points(kop.c$time_hrs[kop.c$conc == unique(kop.c$conc)[i]], kop.c$surv[kop.c$conc == unique(kop.c$conc)[i]], 
           pch = 17, col = i)
    #lines(time, L.3.fx(lc50 = kop.mod$coefficients[4], slp = kop.mod$coefficients[1], t = time), lty = 2)

    lines(time, L.3.fx(lc50 = kop.mod$coefficients[3+i], slp = kop.mod$coefficients[i], t = time), lty = 2, col = i)
  } 

  title(main='Koprivnikar Atrazine-Cercarial mortality (E.trivolvis)')
    legend('bottomleft', legend = c('control', '20ppb', '200ppb'), 
           pch = c(16,17,17), col = c(1,2,3), cex=0.8, bty = 'n')

dev.off()    

#PLot parameter regressions and data
png("Agrochemical_Review/Response_Fxs/Plots/Koprivnikar2006/Koprivnikar_2006_fitted_d-r_parameters_and_functions_cercariae_atrazine.png")

plot(kopc.df$atrazine, kopc.df$e, pch = 16, xlab = 'atrazine (ppb)', ylab = 'LL.2 Parameters',
     ylim = c(0,15), xlim = c(0,300))
  points(kopc.df$atrazine+3, kopc.df$b, pch = 17, col=2)
    for(i in 1:length(kopc.df$atrazine)){
      segments(x0 = kopc.df$atrazine[i], y0 = kopc.df$e[i] + kopc.df$e.se[i],
               x1 = kopc.df$atrazine[i], y1 = kopc.df$e[i] - kopc.df$e.se[i])
      segments(x0 = kopc.df$atrazine[i]+3, y0 = kopc.df$b[i] + kopc.df$b.se[i],
               x1 = kopc.df$atrazine[i]+3, y1 = kopc.df$b[i] - kopc.df$b.se[i], col=2)
    }
  
    lines(c(1:300), sapply(c(1:300), kop.atr.e.lin.pred, simplify = T)[1,], lty = 2)
    lines(c(1:300), sapply(c(1:300), kop.atr.e.lin.pred, simplify = T)[2,], lty = 3)
    lines(c(1:300), sapply(c(1:300), kop.atr.e.lin.pred, simplify = T)[3,], lty = 3)

    lines(c(1:300), sapply(c(1:300), kop.atr.e.exp.pred, simplify = T)[1,], lty = 2, col=3)
    lines(c(1:300), sapply(c(1:300), kop.atr.e.exp.pred, simplify = T)[2,], lty = 3, col=3)
    lines(c(1:300), sapply(c(1:300), kop.atr.e.exp.pred, simplify = T)[3,], lty = 3, col=3)

    lines(c(1:300), sapply(c(1:300), kop.atr.b.lin.pred, simplify = T)[1,], lty = 2, col = 2)
    lines(c(1:300), sapply(c(1:300), kop.atr.b.lin.pred, simplify = T)[2,], lty = 3, col = 2)
    lines(c(1:300), sapply(c(1:300), kop.atr.b.lin.pred, simplify = T)[3,], lty = 3, col = 2)

    lines(c(1:300), sapply(c(1:300), kop.atr.b.exp.pred, simplify = T)[1,], lty = 2, col = 4)
    lines(c(1:300), sapply(c(1:300), kop.atr.b.exp.pred, simplify = T)[2,], lty = 3, col = 4)
    lines(c(1:300), sapply(c(1:300), kop.atr.b.exp.pred, simplify = T)[3,], lty = 3, col = 4)

legend('topright', lty = c(2,2,2,2,3), col = c(1,3,2,4,1), bty = 'n', title = 'Fit models',
       legend = c('LC50 - Linear',
                  'LC50 - Exponential',
                  'Slp - linear',
                  'Slp - exponential',
                  '95% CI'), cex = 0.7)  
legend('top', pch = c(16,17), col = c(1,2), legend = c('LC50', 'slp'),
       cex = 0.8, bty = 'n', title = 'Observed points')

    title('Koprivnikar cercarial survival parameters') 
 
dev.off()    

#plot to test function
    pred.fx.plot = function(He, clr){
      e = as.numeric(predict(kop.atr.e.lin, newdata = data.frame(atrazine = He), se.fit = TRUE)[1:2])
      b = as.numeric(predict(kop.atr.b.exp, newdata = data.frame(logatr = log(He+1)), se.fit = TRUE)[1:2])
      
      e.use = rnorm(1, e[1], e[2])
        while(e.use < 0) e.use = rnorm(1, e[1], e[2])
      b.use = rnorm(1, b[1], b[2])
        while(b.use < 0) b.use = rnorm(1, e[1], e[2])
      
      lines(time, L.3.fx(time, slp = b.use, lc50 = e.use), lty=2, col = clr)
    }   
    
png("Agrochemical_Review/Response_Fxs/Plots/Koprivnikar2006/Koprivnikar_2006_function_simulate_d-r_atrazine.png")

  plot(x = kop.c$time_hrs[kop.c$chem == 'control'], y = kop.c$surv[kop.c$chem == 'control'],
       xlab = 'time (hrs)', ylab = 'prop dead', pch = 16)
    points(x = kop.c$time_hrs[kop.c$conc==20], y = kop.c$surv[kop.c$conc==20], pch = 16, col = 2)
    points(x = kop.c$time_hrs[kop.c$conc==200], y = kop.c$surv[kop.c$conc==200], pch = 16, col = 3)
  
    for(i in c(0,20,200)){
      c = i/20 + 1
      print(c)
      replicate(10, pred.fx.plot(He = i, clr = c))
    }
    
dev.off()

#PLot to test auc function
png("Agrochemical_Review/Response_Fxs/Plots/Koprivnikar2006/Koprivnikar_2006_function_simulate_piC_atrazine.png")

plot(unique(kop.c$conc), kop.atr.aucs.cer/kop.atr.aucs.cer[1],
     xlab = 'atrazine (ppm)', ylab = 'relative AUC (cercariae-hours)',
     pch = 16, ylim = c(0,1), xlim = c(0,300))

set.seed(43093)

  points(c(1:300), sapply(c(1:300), piC_kop_atr_unc, simplify = T),
         pch = 5, col = 4, cex = 0.5)

  legend('bottomleft', pch = c(16,5), col = c(1,4), legend = c("observed", "simulated"),
         cex = 0.8, bty = 'n')
  
dev.off()  
