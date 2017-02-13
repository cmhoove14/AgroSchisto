#Toxicity to Biomphalaria snails from Ibrahim 1992 ###################
# reduction in snail fecundity: from table 1, relative change in number of juveniles produced over the 8 week experiment
snail.repro = data.frame(dose = c(0,125,250,500),
                         juvs = c(4225, 3005, 1749, 1313))

chlor.fN.predict = drm(juvs ~ dose, data = snail.repro, type = 'continuous',
                       fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                                  fixed = c(NA, 0, 4225, NA)))
  summary(chlor.fN.predict)
  plot(chlor.fN.predict)
  
f_N_chlor_ibr92 = function(In){
    predict(chlor.fN.predict, data.frame(dose = In)) / predict(chlor.fN.predict, data.frame(dose = 0))
}

ibr92.fn.df = data.frame(dose = c(0:500),
                         Prediction = 0,
                         Lower = 0,
                         Upper = 0)

ibr92.fn.df[,2:4] <- predict(chlor.fN.predict, newdata = ibr92.fn.df, 
                             interval = 'confidence', level = 0.95)

plot(snail.repro$dose, snail.repro$juvs / snail.repro$juvs[1], pch = 16, ylim = c(0,1),
     ylab = 'relative produced over 8 weeks', xlab = 'ChlorP ppb', main='Ibrahim92 snail recruitment')
  
  lines(ibr92.fn.df$dose, ibr92.fn.df$Prediction / ibr92.fn.df$Prediction[1], col = 2, lty=2)
  lines(ibr92.fn.df$dose, ibr92.fn.df$Lower / ibr92.fn.df$Prediction[1], col = 2, lty=3)
  lines(ibr92.fn.df$dose, ibr92.fn.df$Upper / ibr92.fn.df$Prediction[1], col = 2, lty=3)
  
plot(snail.repro$dose, snail.repro$juvs, pch = 16, ylim = c(0,4500),
      ylab = 'juveniles produced over 8 weeks', xlab = 'ChlorP ppb', main='Ibrahim92 snail recruitment')
  
  lines(ibr92.fn.df$dose, ibr92.fn.df$Prediction, col = 2, lty=2)
  lines(ibr92.fn.df$dose, ibr92.fn.df$Lower, col = 2, lty=3)
  lines(ibr92.fn.df$dose, ibr92.fn.df$Upper, col = 2, lty=3)

#Snail mortality #############
  #-- relative change in survival over entire experiment period (12 weeks)
  snail.mort = data.frame(dose = c(0,125,250,500),
                          total = rep(30,4),
                          dead = c(5,7,8,30))
  
  snail.mort$mort = snail.mort$dead / snail.mort$total
  snail.mort$mean.daily.rate = snail.mort$mort / 84
  
    
  ibr_muNq<-drm(dead/total ~ dose, weights = total, data = snail.mort,
                type = 'binomial', fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"),
                                              fixed = c(NA, 0, 1, NA)))
  
    summary(ibr_muNq)
    plot(ibr_muNq, ylim = c(0,1))
    
  mu_N_chlor_ibr92 = function(In){
      predict(ibr_muNq, data.frame(dose = In))
  }
    
  ibr92.muN.df = data.frame(dose = c(0:500),
                             Prediction = 0,
                             Lower = 0,
                             Upper = 0)
    
    ibr92.muN.df[,2:4] <- predict(ibr_muNq, newdata = ibr92.muN.df, 
                                 interval = 'confidence', level = 0.95)
    
plot(snail.mort$dose, snail.mort$mort, ylim = c(0,1), pch = 16, 
       ylab = '12-week mortality rate', xlab = 'ChlorP ppb', main = 'Ibrahim92 12-week snail mortality')

    
    lines(ibr92.muN.df$dose, ibr92.muN.df$Prediction, col = 2, lty=2)
    lines(ibr92.muN.df$dose, ibr92.muN.df$Lower, col = 2, lty=3)
    lines(ibr92.muN.df$dose, ibr92.muN.df$Upper, col = 2, lty=3)