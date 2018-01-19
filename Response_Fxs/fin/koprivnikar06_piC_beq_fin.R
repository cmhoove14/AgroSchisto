#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(drc)

L.3.fx = function(t, lc50 = lc50, slp = slp){
  1 / (1+exp(slp*log(t / lc50)))
}

#Cercarial mortality (E. trivolvis) from Koprivnikar 2006 ##################
kop.c = read.csv('~/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/Koprivnikar2006.csv')
  kop.c[,5:8] = kop.c[,5:8]/100 #Convert survival measures to proportions
  time = c(0:25)
  kop.c$total = 10 #reported cercariae per treatment is 5-17, 10 per group represents conservative estimate
  
kop.cc = subset(kop.c, chem == 'control')  

kop.mod = drm(surv ~ time_hrs, conc, weights = surv.se^-1, data = kop.c, type = 'binomial', 
              fct = LL.3(names = c('b', 'd', 'e'), 
                         fixed = c(NA, 1, NA)))
  summary(kop.mod)
  plot(kop.mod)

plot(x = kop.c$time_hrs[kop.c$conc == 0], y = kop.c$surv[kop.c$conc == 0],
     xlab = 'time (hrs)', ylab = 'prop alive', pch = 16, cex = 1.2, xlim = c(0,25), ylim = c(0,1))
  lines(time, predict(kop.mod, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  
  for(i in unique(kop.c$conc)[c(2:3)]){
    points(kop.c$time_hrs[kop.c$conc == i], kop.c$surv[kop.c$conc == i], 
           pch = 17, col = i/20+1)
    lines(time, predict(kop.mod, 
                        data.frame(time_hrs=time, conc = i)), lty = 2, col = i/20+1)
  } 

  title(main='Koprivnikar Atrazine-Cercarial mortality (E.trivolvis)')
    legend('bottomleft', legend = c('control', '20ppb', '200ppb'), 
           pch = c(16,17,17), col = c(1,2,3), cex=0.8, bty = 'n')

#Create data frame with parameter values and atrazine concentrations #######################
kopc.df = data.frame(atr = c(0,20,200),
                     logatr = log(c(0,20,200)+1),
                     e = c(summary(kop.mod)$coefficients[c(4:6),1]),
                     e.se = c(summary(kop.mod)$coefficients[c(4:6),2]),
                     b = c(summary(kop.mod)$coefficients[c(1:3),1]),
                     b.se = c(summary(kop.mod)$coefficients[c(1:3),2]))

plot(kopc.df$atr, kopc.df$e, pch = 16, xlab = 'atrazine (ppb)', ylab = 'LL.2 Parameters',
     ylim = c(0,13), xlim = c(0,300))
  points(kopc.df$atr+3, kopc.df$b, pch = 17, col=2)
    for(i in 1:length(kopc.df$atr)){
      segments(x0 = kopc.df$atr[i], y0 = kopc.df$e[i] + kopc.df$e.se[i],
               x1 = kopc.df$atr[i], y1 = kopc.df$e[i] - kopc.df$e.se[i])
      segments(x0 = kopc.df$atr[i]+3, y0 = kopc.df$b[i] + kopc.df$b.se[i],
               x1 = kopc.df$atr[i]+3, y1 = kopc.df$b[i] - kopc.df$b.se[i], col=2)
    }
#parameters as function of atrazine  
ek.mod = lm(e ~ atr, weights = e.se^-1, data = kopc.df) 
ek.mod2 = lm(e ~ logatr, weights = e.se^-1, data = kopc.df)  
  AIC(ek.mod, ek.mod2) #linear fits better
bk.mod = lm(b ~ atr, weights = b.se^-1, data = kopc.df) 

modkdf= data.frame(atr = c(0:300),
                   logatr = log(c(0:300)+1),
                   pred.e = 0,
                   pred.e.se = 0,
                   pred.e2 = 0,
                   pred.e2.se = 0,
                   pred.b = 0,
                   pred.b.se = 0)

modkdf[,3:4] = predict(ek.mod, newdata = modkdf, se.fit = TRUE)[1:2]
modkdf[,5:6] = predict(ek.mod2, newdata = modkdf, se.fit = TRUE)[1:2]
modkdf[,7:8] = predict(bk.mod, newdata = modkdf, se.fit = TRUE)[1:2]

  lines(modkdf$atr, modkdf$pred.e, lty = 2)
    lines(modkdf$atr, modkdf$pred.e + 1.96*modkdf$pred.e.se, lty = 3)
    lines(modkdf$atr, modkdf$pred.e - 1.96*modkdf$pred.e.se, lty = 3)
  
  lines(modkdf$atr, modkdf$pred.e2, lty = 2, col=4)
    lines(modkdf$atr, modkdf$pred.e2 + 1.96*modkdf$pred.e2.se, lty = 3, col=4)
    lines(modkdf$atr, modkdf$pred.e2 - 1.96*modkdf$pred.e2.se, lty = 3, col=4)  
  
  lines(modkdf$atr, modkdf$pred.b, lty = 2, col=2)
    lines(modkdf$atr, modkdf$pred.b + 1.96*modkdf$pred.b.se, lty = 3, col=2)
    lines(modkdf$atr, modkdf$pred.b - 1.96*modkdf$pred.b.se, lty = 3, col=2)

  legend('topright', legend = c('slp', 'lc50'), pch = c(17, 16), col=c(2,1), cex = 0.7, bty='n')  
    title('Koprivnikar cercarial survival parameters') 
    
#Create function to generate d-r function with linear fit to lc50 parameter#####################
    auc.kop.lin.atr0 = function(He){
      e0 = as.numeric(predict(ek.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
      b0 = as.numeric(predict(bk.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
      
      e0.use = rnorm(1, e0[1], e0[2])
        while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
      b0.use = rnorm(1, b0[1], b0[2])
        while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
      auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                       stop.on.error = FALSE)[1]$value
      auc0
    }
    
    piC.kop_atr_unc = function(He){
      #if(He == 0) piC = 1 else{
      e = as.numeric(predict(ek.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
      b = as.numeric(predict(bk.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
      
      e.use = rnorm(1, e[1], e[2])
        while(e.use < 0) e.use = rnorm(1, e[1], e[2])
      b.use = rnorm(1, b[1], b[2])
        while(b.use < 0) b.use = rnorm(1, e[1], e[2])
      auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, stop.on.error = FALSE)[1]$value
      piC = auc/auc.kop.lin.atr0(0) 
      #  if(piC > 1) piC = 1
      #}
      return(piC)
    }  
    
    
#plot to test function
  plot(x = kop.c$time_hrs[kop.c$chem == 'control'], y = kop.c$surv[kop.c$chem == 'control'],
       xlab = 'time (hrs)', ylab = 'prop dead', pch = 16)
    points(x = kop.c$time_hrs[kop.c$conc==20], y = kop.c$surv[kop.c$conc==20], pch = 16, col = 2)
    points(x = kop.c$time_hrs[kop.c$conc==200], y = kop.c$surv[kop.c$conc==200], pch = 16, col = 3)
  
    
    pred.fx.plot = function(He, clr){
      e = as.numeric(predict(ek.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
      b = as.numeric(predict(bk.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
      
      e.use = rnorm(1, e[1], e[2])
        while(e.use < 0) e.use = rnorm(1, e[1], e[2])
      b.use = rnorm(1, b[1], b[2])
        while(b.use < 0) b.use = rnorm(1, b[1], b[2])
      
      lines(time, L.3.fx(time, slp = b.use, lc50 = e.use), lty=2, col = clr)
    }   
    
    for(i in c(0,20,200)){
      c = i/20 + 1
      print(c)
      replicate(10, pred.fx.plot(He = i, clr = c))
    }

keep.kop06.beq = c('auc.kop.lin.atr0', 'piC.kop_atr_unc', 'ek.mod', 'bk.mod')        