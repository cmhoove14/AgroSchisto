#pool data from ibrahim 92 and hasheesh 11 to estimate reduction in fecundity of snail population by chlorpyrifos
source('Response_Fxs/Ibrahim92_fN_muN.R')
ib = c(4322/25, 3183/25, 1946/25) #total eggs/number snails from ibrahim
rels = c(1, ib[2]/ib[1], ib[3]/ib[1], 9.22/132.21) #last entry is eggs/snail in chlorP treatment / control from Hasheesh
concs = c(0, 125, 250, 720)

ch.fn<-nls(rels ~ exp(-b * concs), start = list(b = 0.003))

  ch.fN.meta = function(In){
    exp(-summary(ch.fn)$parameters[1] * In)
  }

plot(concs, rels, pch = 16, ylim = c(0,1), xlab = 'Chlorpyrifos concentration',
     ylab = 'relative fecundity')
  lines(c(0:750), f_Nq_chlor_ibrahim92(c(0:750)), lty = 2, col = 'red')
  lines(c(0:750), ch.fN.meta(c(0:750)), lty = 3, col = 'red')
  
#Especially in the realm of environmentally relevant concentrations, Ibrahim provides more info
  #and the additional data point from Hasheesh 2011 fits in with the relationship derived from Ibrahim data
  #so we'll use the ibrahim function 
  