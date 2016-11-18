require(ggplot2)
require(drc)

#Atrazine effect on phi_Nq (snail carrying capacity) from Baxter et al 2011 data with Rohr et al analysis####
atra.df<-data.frame('atra' = c(0,1,10,30,100),              #Raw atrazine concentration (ppb)
                    'logatra' = log(c(0,1,10,30,100)+1),    #Log atrazine concentration (ppb)
                    'phiNq' = c(0,0.2888,0.6535,0,1.3215))  #Snail population response measured as peak snail growth rate and 
                                                            #interpreted as changes in snail carrying capacity

plot(atra.df$logatra, atra.df$phiNq, pch = 16)

atra_mod<-glm(phiNq ~ logatra+0, data=atra.df)

atra.slope = atra_mod$coefficients[1]

f_phi_Nq_at = function(He){ #Function to use in model
  phi_Nq = parameters['phi_N'] + parameters['phi_N'] * (atra.slope*log(He+1))
  return(phi_Nq)
}