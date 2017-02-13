require(drc)

#Incorporate data from all three studies investigating atrazine/E. trivolvis cercarial die-off relationship
#This would probably work better with a random effects model or something of that sort, but will have to revisit later

source('Response_fxs/rohr08_piC.R')
source('Response_fxs/griggs08_piC.R')
source('Response_fxs/koprivnikar06_piC.R')

atr.rel.auc = c(1, as.numeric(rel.auc.rohr[2]), auc.kop[c(2,3)] / auc.kop[1], auc.grg[c(2,3)] / auc.grg[1])
atr.conc = c(0, 201, 20, 200, 15, 100)

atr.meta = drm(atr.rel.auc ~ atr.conc, fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                          fixed = c(NA, NA, NA, NA)))
  summary(atr.meta)
  plot(atr.meta)
  
  plot(x=atr.conc, y=atr.rel.auc, pch = 16, xlab = 'Atrazine (ppb)', ylab = 'Relative daily auc', ylim = c(0,1))
  
  piC.meta.df = data.frame(conc = c(0:200),
                           Prediction = 0,
                           Lower = 0,
                           Upper = 0)
  
  piC.meta.df[,2:4] <- predict(atr.meta, newdata = piC.meta.df, 
                               interval = 'confidence', level = 0.95)
  
  lines(piC.meta.df$conc, piC.meta.df$Prediction / piC.meta.df$Prediction[1], lty=2)
  lines(piC.meta.df$conc, piC.meta.df$Lower / piC.meta.df$Prediction[1], lty=3)
  lines(piC.meta.df$conc, piC.meta.df$Upper / piC.meta.df$Prediction[1], lty=3)
  
  title(main = expression(paste(pi[C], '(atrazine) meta')))