require(drc)
require(lme4)

#Incorporate data from all three studies investigating atrazine/E. trivolvis cercarial die-off relationship

source('Response_fxs/fin/rohr08_piC_beq_fin.R')
source('Response_fxs/fin/griggs08_piC_beq_fin.R')
source('Response_fxs/fin/koprivnikar06_piC_beq_fin.R')

#Griggs data
plot(x = cerc.g0$time_hrs[cerc.g0$conc == 0], y = cerc.g0$prop_surv[cerc.g0$conc == 0],
     xlab = 'time (hrs)', ylab = 'prop alive', pch = 16, cex = 1.2, xlim = c(0,25), ylim = c(0,1))
  lines(time, predict(grg.mod, 
                      data.frame(time_hrs=time, conc = 0)), lty = 2)
  for(i in unique(cerc.g0$conc)[c(2:3)]){
    points(cerc.g0$time_hrs[cerc.g0$conc == i], cerc.g0$prop_surv[cerc.g0$conc == i], 
           pch = 17, col = i+1)
    lines(time, predict(grg.mod, 
                        data.frame(time_hrs=time, conc = i)), lty = 2, col = i+1)
  } 
#Koprivnikar data  
  points(x = kop.c$time_hrs[kop.c$chem == 'control'], y = kop.c$surv[kop.c$chem == 'control'], col=2, pch = 16)
    lines(time, L.3.fx(time, lc50 = kopc.df[1,3], slp = kopc.df[1,5]), lty=2, col = 2)
          
  points(x = kop.c$time_hrs[kop.c$conc==20], y = kop.c$surv[kop.c$conc==20], pch = 17, col=2)
    lines(time, L.3.fx(time, lc50 = kopc.df[2,3], slp = kopc.df[2,5]), lty=2, col = 2)
    
  points(x = kop.c$time_hrs[kop.c$conc==200], y = kop.c$surv[kop.c$conc==200], pch = 15, col=2)
    lines(time, L.3.fx(time, lc50 = kopc.df[3,3], slp = kopc.df[3,5]), lty=2, col = 2)
#Rohr data  
  points(x = rohr$time_hrs[rohr$chem == 'control'], y = rohr$surv[rohr$chem == 'control']/100,
         pch = 16, col=4)
    lines(time, L.3.fx(time, coef(rohr.ctrl)[2], coef(rohr.ctrl)[1]), lty=2)
  
  points(x = rohr.atr$time_hrs, y = rohr.atr$surv/100,
         col = 4, pch = 15)
    lines(time, L.3.fx(time, coef(rohr.atrmod)[2], coef(rohr.atrmod)[1]), lty=2, col = 4)

  eb.meta = data.frame(atr = c(kopc.df$atr, grgc.df$atr, rohr.fin$conc[1:2]),
                       study = c(rep('kop', length(kopc.df$atr)), 
                                 rep('grg', length(grgc.df$atr)),
                                 rep('rohr', length(rohr.fin$conc[1:2]))),
                       e = c(kopc.df$e, grgc.df$e, rohr.fin$e[1:2]),
                       e.se = c(kopc.df$e.se, grgc.df$e.se, rohr.fin$e.se[1:2]),
                       b = c(kopc.df$b, grgc.df$b, rohr.fin$b[1:2]),
                       b.se = c(kopc.df$b.se, grgc.df$b.se, rohr.fin$b.se[1:2]))
  
#b (lc50) model and plots #################  
  plot(eb.meta$atr, eb.meta$b, pch = 17, ylim = c(0,20), ylab = 'tox params',
       xlab = 'atrazine (ppb)', xlim = c(0,500), col=2)
    points(eb.meta$atr, eb.meta$e, pch = 16)
    for(i in 1:length(eb.meta$atr)){
      segments(x0 = eb.meta$atr[i], y0 = eb.meta$e[i] + eb.meta$e.se[i],
               x1 = eb.meta$atr[i], y1 = eb.meta$e[i] - eb.meta$e.se[i])
      segments(x0 = eb.meta$atr[i], y0 = eb.meta$b[i] + eb.meta$b.se[i],
               x1 = eb.meta$atr[i], y1 = eb.meta$b[i] - eb.meta$b.se[i], col=2)
    }
    
    legend('topright', legend = c('slp', 'lc50'), pch = c(17, 16), col=c(2,1), cex = 0.7, bty = 'n')  
  
  #Model two-parameter log-logistic parameters as a function of atrazine concentration
    atrmeta.bme = lmer(b ~ atr + (1|study), weights = b.se^-1, data = eb.meta)
    atrmeta.eme = lmer(e ~ atr + (1|study), weights = e.se^-1, data = eb.meta)
    
    atrmeta.b = lm(b ~ atr, weights = b.se^-1, data = eb.meta)
    atrmeta.e = lm(e ~ atr, weights = e.se^-1, data = eb.meta)

  atrmeta.df = data.frame(atr = c(0:500),
                          e.pred = 0,
                          e.se = 0,
                          b.pred = 0,
                          b.se = 0,
                          eme.pred = 0,
                          eme.se = 0,
                          bme.pred = 0,
                          bme.se = 0)
  
  #Predictors using simple linear regression
    atrmeta.df[,2:3] <- predict(atrmeta.e, newdata = atrmeta.df, se.fit = TRUE)[1:2]
    atrmeta.df[,4:5] <- predict(atrmeta.b, newdata = atrmeta.df, se.fit = TRUE)[1:2]
  #Predictors using mixed effects model  
    atrmeta.df[,6] <- predict(atrmeta.eme, newdata = atrmeta.df, re.form = NA)
    atrmeta.df[,8] <- predict(atrmeta.bme, newdata = atrmeta.df, re.form = NA)
  
  lines(atrmeta.df$atr, atrmeta.df$b.pred, lty=2, col=2)
    lines(atrmeta.df$atr, atrmeta.df$b.pred + 1.96 * atrmeta.df$b.se, lty=3, col=2)
    lines(atrmeta.df$atr, atrmeta.df$b.pred - 1.96 * atrmeta.df$b.se, lty=3, col=2)
    
  lines(atrmeta.df$atr, atrmeta.df$e.pred, lty=2)
    lines(atrmeta.df$atr, atrmeta.df$e.pred + 1.96 * atrmeta.df$e.se, lty=3)
    lines(atrmeta.df$atr, atrmeta.df$e.pred - 1.96 * atrmeta.df$e.se, lty=3) 
  #Compare estimates from mixed effects model  
    lines(atrmeta.df$atr, atrmeta.df$bme.pred, lty=2, col=4, lwd = 2)
    lines(atrmeta.df$atr, atrmeta.df$eme.pred, lty=2, col=4, lwd = 2)
 
  legend(350, 20, legend = c('lm model slp', 'lm model lc50', 'ME model', '95% CIs'),
         lty = c(2,2,2,3), lwd = c(1,1,2,1), col = c(1,2,4,1), cex = 0.6, bty = 'n')  
  title('Pooled cercarial survival (atrazine) parameters')
  
#Create function to generate d-r function #####################
  meta.atr.piC.fx = function(He){
    e = as.numeric(predict(atrmeta.e, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(atrmeta.b, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    
    while(e.use <= 0)     e.use = rnorm(1, e[1], e[2])
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24)[1]$value

    return(auc)
  }  
  
#Final:generate relative cercariae-hrs function and plot sample output
  piC.meta_atr_unc = function(He){
      piC = meta.atr.piC.fx(He) / meta.atr.piC.fx(0)
    
    return(piC)
  }
  
  plot(x = c(0:200), y = sapply(c(0:200), piC.meta_atr_unc, simplify = T), pch = 17, col = 6, cex = 0.6)
  
  keep.meta.piC = c('piC.meta_atr_unc', 'meta.atr.piC.fx', 'atrmeta.e', 'atrmeta.b')  