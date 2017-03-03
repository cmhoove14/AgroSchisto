atr = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Cercarial Mortality/rohr08_atr.csv')

plot(atr$log_conc, atr$surv, xlab = 'log+1 Atrazine (ppb)', ylab = 'Prop alive @14-18 hrs', 
     xlim = c(0,8), ylim = c(0.25,0.65), pch = 16)
  for(i in 1:length(atr$log_conc)){
    segments(x0 = atr$log_conc[i], y0 = 1-atr$hi[i],
             x1 = atr$log_conc[i], y1 = 1-atr$lo[i])
  }

atr.lin = lm(surv ~ log_conc, weights = st_err^-1, data = atr)

atr.drm = drm(alive/total ~ conc, total, data=atr, type = 'binomial',
              fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                         fixed = c(NA, 0, atr$surv[1], NA)))
  
  atr.test = data.frame(conc = c(0:2000),
                        log_conc = log(c(0:2000)+1),
                        pred.lin = 0,
                        pred.lin.hi = 0,
                        pred.lin.lo = 0,
                        pred.drm = 0,
                        pred.drm.hi = 0,
                        pred.drm.lo = 0)

    atr.test[,3:5] = predict(atr.lin, newdata = atr.test, interval = 'confidence', level = 0.95)
    
      lines(atr.test$log_conc, atr.test$pred.lin, lty=2)
      lines(atr.test$log_conc, atr.test$pred.lin.hi, lty=3)
      lines(atr.test$log_conc, atr.test$pred.lin.lo, lty=3)
      
    atr.test[,6:8] = predict(atr.drm, newdata = atr.test, interval = 'confidence', level = 0.95) 
    
      lines(atr.test$log_conc, atr.test$pred.drm, lty=2, col=2)
      lines(atr.test$log_conc, atr.test$pred.drm.hi, lty=3, col=2)
      lines(atr.test$log_conc, atr.test$pred.drm.lo, lty=3, col=2)
      
#Final function estimates relative mortality to He=0 and is interpreted as pi_C parameter      
  piC.atr.rohr08.uncertainty = function(He){
    rdrm(1, LL.2(), coef(atr.drm), He, yerror = 'rbinom', ypar = 100)$y / 100
  }    
  
  piC.atr.rohr08.lin = function(He){
    if(He == 0){piC = 1} else {
      ts = as.numeric(predict(atr.lin, newdata = data.frame(log_conc = log(He+1)), se.fit = TRUE)[1:2])
      piC = rnorm(1, ts[1], ts[2]) / atr.lin$coefficients[1]
    }

      piC
  } 
  
  
#Plot to see how it does
  plot(atr$conc, atr$surv / atr$surv[1], pch = 16, xlim = c(0,2000), ylim = c(0,1),
       ylab = 'rel prop alive @ 14-18 hrs', xlab = 'Atrazine (ppb)',
       main = 'Sample output of cercarial mortality')
    points(seq(0, 2000, 5), sapply(seq(0, 2000, 5), piC.atr.rohr08.lin), pch = 17, col=3, cex = 0.5) 
    legend('bottomleft', legend = c('Observed', 'Modeled'), pch = c(16,17), col = c(1,3), cex = 0.7)

#Vector of items to keep
  keep.atr.rohr08 = c('atr', 'piC.atr.rohr08.uncertainty', 'atr.drm', 'atr.lin', 'piC.atr.rohr08.lin')