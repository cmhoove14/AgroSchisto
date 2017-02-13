require(drc)
require(lme4)

#Incorporate data from all three studies investigating atrazine/E. trivolvis cercarial die-off relationship
#This would probably work better with a random effects model or something of that sort, but will revisit later

source('Response_fxs/rohr08_piC_beq.R')
source('Response_fxs/griggs08_piC_beq.R')
source('Response_fxs/koprivnikar06_piC_beq.R')

#Griggs data
plot(x = cerc.g0$time_hrs,y = cerc.g0$surv/100,
     xlab = 'time (hrs)', ylab = 'prop alive', pch = 16, xlim = c(0,25), ylim = c(0,1))
    lines(time, ll4(1,0,coef(grg.ctrl)[1], coef(grg.ctrl)[2], time), lty=2)
  points(cerc.g15$time_hrs, cerc.g15$surv/100, pch = 17)
    lines(time, ll4(1,0,coef(grg.15)[1], coef(grg.15)[2], time), lty=3)
  points(cerc.g100$time_hrs, cerc.g100$surv/100, pch = 15)
    lines(time, ll4(1,0,coef(grg.100)[1], coef(grg.100)[2], time), lty=4)
#Koprivnikar data  
  points(x = kop.cc$time_hrs, y = kop.cc$surv, col=2, pch = 16)
    lines(time, ll4(1,0,summary(kop.ctrl)$coefficients[1],
                    summary(kop.ctrl)$coefficients[2], time), lty=2, col=2)
  points(x = kop.20$time_hrs, y = kop.20$surv, pch = 17, col=2)
    lines(time, ll4(1,0,summary(kop.20mod)$coefficients[1],
                    summary(kop.20mod)$coefficients[2], time), lty=3, col=2)
  points(x = kop.200$time_hrs, y = kop.200$surv, pch = 15, col=2)
    lines(time, ll4(1,0,summary(kop.200mod)$coefficients[1],
                    summary(kop.200mod)$coefficients[2], time), lty=4, col=2)
#Rohr data  
  points(x = rohr$time_hrs[rohr$chem == 'control'], y = rohr$surv[rohr$chem == 'control']/100,
         pch = 16, col=4)
    lines(time, ll4(1,0,coef(rohr.ctrl)[1], coef(rohr.ctrl)[2], time), lty=2, col=4)
  points(x = rohr.atr$time_hrs, y = rohr.atr$surv/100,
         col = 4, pch = 15)
    lines(time, ll4(1,0,coef(rohr.atrmod)[1], coef(rohr.atrmod)[2], time), lty=2, col = 4)

  eb.meta = data.frame(atr = c(kopatr.df$atr, grgc.df$atr, rohr.fin$conc[1:2]),
                       study = c(rep('kop', length(kopatr.df$atr)), 
                                 rep('grg', length(grgc.df$atr)),
                                 rep('rohr', length(rohr.fin$conc[1:2]))),
                       e = c(kopatr.df$e, grgc.df$e, rohr.fin$e[1:2]),
                       e.se = c(kopatr.df$e.se, grgc.df$e.se, rohr.fin$e.se[1:2]),
                       b = c(kopatr.df$b, grgc.df$b, rohr.fin$b[1:2]),
                       b.se = c(kopatr.df$b.se, grgc.df$b.se, rohr.fin$b.se[1:2]))
  
#b (lc50) model and plots #################  
  plot(eb.meta$atr, eb.meta$b, pch = 17, ylim = c(0,20), ylab = 'tox params',
       xlab = 'atrazine (ppb)', xlim = c(0,210), col=2)
    points(eb.meta$atr, eb.meta$e, pch = 16)
    for(i in 1:length(eb.meta$atr)){
      segments(x0 = eb.meta$atr[i], y0 = eb.meta$e[i] + eb.meta$e.se[i],
               x1 = eb.meta$atr[i], y1 = eb.meta$e[i] - eb.meta$e.se[i])
      segments(x0 = eb.meta$atr[i], y0 = eb.meta$b[i] + eb.meta$b.se[i],
               x1 = eb.meta$atr[i], y1 = eb.meta$b[i] - eb.meta$b.se[i], col=2)
    }
    
    legend('topright', legend = c('slp', 'lc50'), pch = c(17, 16), col=c(2,1), cex = 0.7)  
  
  #Model two-parameter log-logistic parameters as a function of atrazine concentration
    atrmeta.bme = lmer(b ~ atr + (1|study), weights = b.se^-1, data = eb.meta)
    atrmeta.eme = lmer(e ~ atr + (1|study), weights = e.se^-1, data = eb.meta)
    
    atrmeta.b = lm(b ~ atr, weights = b.se^-1, data = eb.meta)
    atrmeta.e = lm(e ~ atr, weights = e.se^-1, data = eb.meta)

  atrmeta.df = data.frame(atr = c(0:210),
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
    #lines(atrmeta.df$atr, atrmeta.df$bme.pred, lty=2, col=4)
    #lines(atrmeta.df$atr, atrmeta.df$eme.pred, lty=2, col=4)
 
    
  title('Pooled cercarial survival (atrazine) parameters')
  
#Create function to generate d-r function #####################
  meta.atr.piC.fx = function(He){
    e = as.numeric(predict(atrmeta.e, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(atrmeta.b, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
    
    auc = integrate(f = function(t) {(1/(1+exp(b.use*(log(t / e.use)))))}, 
                    lower=0, upper=24)[1]$value
    auc
  }  
  
#Final:generate relative cercariae-hrs function  
  piC.meta_atr_unc = function(He){
    piC = meta.atr.piC.fx(He) / meta.atr.piC.fx(0)
    if(piC > 1) piC = 1
    else(return(piC))
  }
  
  keep.meta.piC = c('piC.meta_atr_unc', 'meta.atr.piC.fx', 'atrmeta.e', 'atrmeta.b')  