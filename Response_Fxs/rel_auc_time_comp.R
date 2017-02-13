source('Response_Fxs/tchounwou92_piC.R')


#See how integrating the survival curve at differnt times affects the relative auc estimates
auc.rels = data.frame(time = c(1:24),
                      mal0 = 0,
                      mal50 = 0,
                      mal100 = 0,
                      mal150 = 0,
                      mal200 = 0,
                      mal250 = 0)

for(i in 1:24){
  auc.rels[i,2] = integrate(f = tch.cerc.mal0.surv, lower=0, upper=i)[1]$value
  auc.rels[i,3] = integrate(f = tch.cerc.mal50.surv, lower=0, upper=i)[1]$value
  auc.rels[i,4] = integrate(f = tch.cerc.mal100.surv, lower=0, upper=i)[1]$value
  auc.rels[i,5] = integrate(f = tch.cerc.mal150.surv, lower=0, upper=i)[1]$value
  auc.rels[i,6] = integrate(f = tch.cerc.mal200.surv, lower=0, upper=i)[1]$value
  auc.rels[i,7] = integrate(f = tch.cerc.mal250.surv, lower=0, upper=i)[1]$value
}

plot(mal.con, auc.rels[24,2:7] / auc.rels[24,2], pch = 16, ylim = c(0,1),
     xlab = 'malathion conc (ppb)', ylab = 'rel auc')
  for(i in c(3,6,9,12,15,18,21)){
    points(mal.con, auc.rels[i,2:7] / auc.rels[i,2], pch = 17, col = i/3 + 1)
  }

legend('topright', title = 'auc time (hrs)', legend = c(3,6,9,12,15,18,21,24),
       pch = c(rep(17,7),16), col = c(2:8,1), cex = 0.7)