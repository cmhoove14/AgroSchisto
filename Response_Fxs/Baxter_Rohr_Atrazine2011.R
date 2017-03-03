require(ggplot2)
require(drc)

#Atrazine effect on phi_Nq (snail carrying capacity) from Baxter et al 2011 data with Rohr et al analysis####
atra.df<-data.frame('atra' = c(0,1,10,30,100),              #Raw atrazine concentration (ppb)
                    'logatra' = log(c(0,1,10,30,100)+1),    #Log atrazine concentration (ppb)
                    'growthrate' = c(0.119406, 0.153891, 0.19744, 0.118918, 0.27719), #Peak growth rate in Rohr reanalysis of Baxter data  
                    'st.err' = c(0.015949, 0.026329, 0.035467, 0.016286, 0.047113))     #interpreted as changes in snail carrying capacity 

plot(atra.df$atra, atra.df$growthrate, pch = 16, ylim = c(0,0.3),
     xlab = 'atrazine (ppb)', ylab = 'growth rate')
  segments(x0 = atra.df$atra, y0 = atra.df$growthrate + atra.df$st.err,
           x1 = atra.df$atra, y1 = atra.df$growthrate - atra.df$st.err)

bax.mod = lm(growthrate ~ logatra, weights = st.err^-1, data = atra.df)
  summary(bax.mod)
 
atra.df.no30 = subset(atra.df, atra != 30)
bax.mod.no30 = lm(growthrate ~ logatra, weights = st.err^-1, data = atra.df.no30)
  summary(bax.mod.no30)

  bax.phin.df = data.frame(atra = c(0:200),
                           logatra = log(c(0:200)+1),
                           Prediction = 0,
                           st.err = 0,
                           Prediction.no30 = 0,
                           st.err.no30 = 0)
  
  bax.phin.df[,3:4] <- predict(bax.mod, newdata = bax.phin.df, 
                               type = 'response', se.fit = T)[1:2]
  bax.phin.df[,5:6] <- predict(bax.mod.no30, newdata = bax.phin.df, 
                               type = 'response', se.fit = T)[1:2]
  
  lines(bax.phin.df$atra, bax.phin.df$Prediction, lty=2)
    lines(bax.phin.df$atra, bax.phin.df$Prediction + 1.96*bax.phin.df$st.err, lty=3)
    lines(bax.phin.df$atra, bax.phin.df$Prediction - 1.96*bax.phin.df$st.err, lty=3)  
  lines(bax.phin.df$atra, bax.phin.df$Prediction.no30, lty=2, col=2)
    lines(bax.phin.df$atra, bax.phin.df$Prediction.no30 + 1.96*bax.phin.df$st.err.no30, lty=3, col=2)
    lines(bax.phin.df$atra, bax.phin.df$Prediction.no30 - 1.96*bax.phin.df$st.err.no30, lty=3, col=2)  


phi_Nq_atr_baxrohr = function(He){
  if(He == 0){phiNq = 1} else{
    u = predict(bax.mod, newdata = data.frame(logatra = log(He+1)), type = 'response',
                se.fit = TRUE)[1:2]
    phiNq = rnorm(1, u$fit, u$se.fit) / predict(bax.mod, newdata = data.frame(logatra = 0), 
                                                type = 'response')
  }
  
  phiNq
  
}

phi_Nq_atr_baxrohr.no30 = function(He){
  if(He == 0){phiNq = 1} else{
    u = predict(bax.mod.no30, newdata = data.frame(logatra = log(He+1)), type = 'response',
                se.fit = TRUE)[1:2]
    phiNq = rnorm(1, u$fit, u$se.fit) / predict(bax.mod.no30, newdata = data.frame(logatra = 0), 
                                                type = 'response')
  }
  
  phiNq
  
}

plot(atra.df$atra, atra.df$growthrate / atra.df$growthrate[1], 
     pch = 16, ylim = c(0.8, 3.0),
     xlab = 'atrazine (ppb)', ylab = 'growth rate',
     main = 'atrazine:peak growth rate (~carrying capacity)')
#Add error bars
  segments(x0 = atra.df$atra, 
           y0 = (atra.df$growthrate + atra.df$st.err) / atra.df$growthrate[1],
           x1 = atra.df$atra, 
           y1 = (atra.df$growthrate - atra.df$st.err) / atra.df$growthrate[1])
#add model predictions  
    points(c(0:200), sapply(c(0:200), phi_Nq_atr_baxrohr), pch = 17, col=4, cex=0.6)
    points(c(0:200), sapply(c(0:200), phi_Nq_atr_baxrohr.no30), pch = 17, col=2, cex=0.6)

  legend('topleft', legend = c('30ppb included', '30ppb excluded'), pch = 17, col = c(4,2), cex = 0.7)
  

  keep.baxrohr = c('atra.df', 'phi_Nq_atr_baxrohr', 'bax.mod', 'phi_Nq_atr_baxrohr.no30', 'bax.mod.no30')