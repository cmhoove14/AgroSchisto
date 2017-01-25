#Data extraction and model fitting to Ragab 2006 data
require(drc)

#Snail toxicity ##########
#Ammonium Nitrate
lc50.n.amm<-470
slp.n.amm<-1.3

mu_Nq_amm_rag06<-function(In){
  Ins = In/1000
  1 - 1/(1+exp(slp.n.amm*(log(Ins)-log(lc50.n.amm))))
} 

mu_Nq_amm_rag06(lc50.n.amm*1000)

plot(c(0:1000000), mu_Nq_amm_rag06(c(0:1000000)), lwd = 2, type = 'l', xlab = 'ammonium nitrate concentration (ppb)',
     ylab = 'mu_Nq', ylim = c(0,1), main = 'ammonium nitrate toxicity to snails, Ragab06')

#Pottasium Sulfate     
lc50.n.pot<-1900
slp.n.pot<-1.27

mu_Nq_pot_rag06<-function(In){
  Ins = In/1000
  1 - 1/(1+exp(slp.n.pot*(log(Ins)-log(lc50.n.pot))))
} 

mu_Nq_pot_rag06(lc50.n.pot*1000)

plot(c(0:1000000), mu_Nq_pot_rag06(c(0:1000000)), lwd = 2, type = 'l', xlab = 'potassium sulfate concentration (ppb)',
     ylab = 'mu_Nq', ylim = c(0,1), main = 'potassium sulfate toxicity to snails, Ragab06')

#Urea    
lc50.n.ure<-22000
slp.n.ure<-1.28

mu_Nq_ure_rag06<-function(In){
  Ins = In/1000
  1 - 1/(1+exp(slp.n.ure*(log(Ins)-log(lc50.n.ure))))
} 

mu_Nq_ure_rag06(lc50.n.ure*1000)

plot(c(0:10000000), mu_Nq_ure_rag06(c(0:10000000)), lwd = 2, type = 'l', xlab = 'urea concentration (ppb)',
     ylab = 'mu_Nq', ylim = c(0,1), main = 'urea toxicity to snails, Ragab06')
