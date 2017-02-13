require(drc)

#Incorporate data from all three studies investigating atrazine/E. trivolvis cercarial die-off relationship
#This would probably work better with a random effects model or something of that sort, but will revisit later

source('Response_fxs/piC_atr_meta_beq.R')

   meta.atr.piC.plot = function(He, clr){
      e = as.numeric(predict(atrmeta.e, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
      b = as.numeric(predict(atrmeta.b, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
      
      e.use = rnorm(1, e[1], e[2])
      b.use = rnorm(1, b[1], b[2])
      
      lines(time, 1/(1+exp(b.use*(log(time / e.use)))), lty = 2, col=clr)
      
    }  
#Griggs data
plot(x = cerc.g0$time_hrs,y = cerc.g0$surv/100,
     xlab = 'time (hrs)', ylab = 'prop alive', pch = 16, xlim = c(0,25), ylim = c(0,1))
  points(cerc.g15$time_hrs, cerc.g15$surv/100, pch = 17, col=2)
  points(cerc.g100$time_hrs, cerc.g100$surv/100, pch = 15,col=4)
  
  set.seed(100)
  
  replicate(20, meta.atr.piC.plot(0,1))
  replicate(20, meta.atr.piC.plot(15,2))
  replicate(20, meta.atr.piC.plot(100,4))
  

 
#Koprivnikar data  
plot(x = kop.cc$time_hrs, y = kop.cc$surv, pch = 16,
     xlab = 'time (hrs)', ylab = 'prop alive', xlim = c(0,25), ylim = c(0,1))
  points(x = kop.20$time_hrs, y = kop.20$surv, pch = 17, col=2)
  points(x = kop.200$time_hrs, y = kop.200$surv, pch = 15, col=4)
  
  replicate(20, meta.atr.piC.plot(0,1))
  replicate(20, meta.atr.piC.plot(20,2))
  replicate(20, meta.atr.piC.plot(200,4))

#Rohr data  
plot(x = rohr$time_hrs[rohr$chem == 'control'], y = rohr$surv[rohr$chem == 'control']/100,
       pch = 16, xlab = 'time (hrs)', ylab = 'prop alive', xlim = c(0,25), ylim = c(0,1))
  points(x = rohr.atr$time_hrs, y = rohr.atr$surv/100,
         col = 2, pch = 15)
  replicate(20, meta.atr.piC.plot(0,1))
  replicate(20, meta.atr.piC.plot(201,2))

#These are just visual checks of how the modeled functions fit the original data
#They seem to fit well, though the Rohr study appears to have longer survival than predicted by the models
#Will look into more formal assessments of model fits when time permits