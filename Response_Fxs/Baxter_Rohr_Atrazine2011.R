require(ggplot2)
require(drc)

#Atrazine effect on phi_Nq (snail carrying capacity) from Baxter et al 2011 data with Rohr et al analysis####
atra.df<-data.frame('atra' = c(0,1,10,30,100),              #Raw atrazine concentration (ppb)
                    'logatra' = log(c(0,1,10,30,100)+1),    #Log atrazine concentration (ppb)
                    'phiNq' = c(1,1.2888,1.6535,0.9960,2.3215))  #Snail population response measured as peak snail growth rate and 
                                                            #interpreted as changes in snail carrying capacity

plot(atra.df$logatra, atra.df$phiNq, pch = 16, xlab = 'log+1 atrazine (ppb)', ylab = 'Relative peak growth rate')

atr_phiNq<-nls(phiNq ~ 1 + b * log(atra+1), data = atra.df, start = list(b = 0.25))
  summary(atr_phiNq)

atr.con = c(0:100)

phi_Nq_atr_baxrohr<-function(Fe){
  1 + (summary(atr_phiNq)$parameters[1] * log(Fe + 1))
}

lines(log(atr.con+1), phi_Nq_atr_baxrohr(atr.con), lty=2, col='red')