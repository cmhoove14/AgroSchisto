require(drc)
source('Response_fxs/rohr08_piC.R')
source('Response_fxs/griggs08_piC.R')
source('Response_fxs/koprivnikar06_piC.R')

atr.rel.auc = c(1, as.numeric(rel.auc.rohr[2]),rel.auc.kop[c(2,3)], rel.auc.grg[c(2,3)])
atr.conc = c(0, 201, 20, 200, 15, 100)

plot(x=log(atr.conc+1), y=atr.rel.auc, pch = 16, xlab = 'log+1 Atrazine (ppb)', ylab = 'Relative daily auc', ylim = c(0,1))

plot(c(0,20,200), rel.auc.kop, pch = 16, ylim=c(0,1), ylab = 'relative auc', xlim = c(0,210),
     xlab = 'atrazine concentration (ppb)', main = 'Atrazine toxicity to E. trivolvis pooled')
  lines(atr.k.con, pi_C_atr_kop06(atr.k.con), lty=2)
 points(c(15,100), rel.auc.grg[c(2,3)], pch = 16, col='red')
  lines(atr.k.con, pi_C_atr_grg08(atr.k.con), lty=2, col='red')
 points(201, as.numeric(rel.auc.rohr[2]), pch = 16, col='blue')
  lines(atr.k.con, pi_C_atr_rohr08(atr.k.con), lty=2, col='blue')

atr.piC.fx<-nls(atr.rel.auc ~ exp(-b*atr.conc), start = list(b = 0.01))  
  summary(atr.piC.fx)
  
  pi_C_atr_all = function(He){
    exp(-summary(atr.piC.fx)$parameters[1]*(He)) 
  } 
  
    lines(atr.con, pi_C_atr_all(atr.con), lty=2, col='green')
    
  legend('bottomleft', legend = c('Koprivnikar06', 'Griggs08', 'Rohr08', 'Pooled'),
         pch = 16, col = c('black', 'red', 'blue', 'green'), cex=0.8)