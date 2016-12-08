require(drc)
source('Response_fxs/rohr08_piC.R')
source('Response_fxs/griggs08_piC.R')
source('Response_fxs/koprivnikar06_piC.R')

plot(c(0,20,200), rel.auc.kop, pch = 16, ylim=c(0,1), ylab = 'relative auc', xlim = c(0,210),
     xlab = 'atrazine conentration (ppb)')
  lines(atr.k.con, pi_C_atr_kop06(atr.k.con), lty=2)
  points(c(0,15,100), rel.auc.grg, pch = 16, col='red')
  lines(atr.k.con, pi_C_atr_grg08(atr.k.con), lty=2, col='red')
  points(201, as.numeric(rel.auc.rohr[2]), pch = 16, col='blue')
  
atr.rel.auc = c(1, as.numeric(rel.auc.rohr[2]),rel.auc.kop[c(2,3)], rel.auc.grg)
atr.conc = c(0, 201, 20, 200, 15, 100)

plot(x=log(atr.conc+1), y=atr.rel.auc, pch = 16, xlab = 'log+1 Atrazine (ppb)', ylab = 'Relative daily auc', ylim = c(0,1))

atr.auc = c(auc.grg.ctrl, auc.kop.ctrl, auc.rohr.ctrl,
            auc.grg.atr15, auc.grg.atr100,
            auc.kop.atr20, auc.kop.atr200, auc.rohr.atr201)

atr.conc2 = c(0,0,0,15,100,20,200,201)

plot(x=log(atr.conc2+1), y=atr.auc, pch = 16, xlab = 'log+1 Atrazine (ppb)', ylab = 'Daily auc')


atr.piC.fx<-glm(atr.rel.auc ~ log(atr.conc+1), family = binomial('logit'))