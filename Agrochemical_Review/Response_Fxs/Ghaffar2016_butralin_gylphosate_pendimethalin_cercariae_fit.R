#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Ghaffar 2016 LARVAL (S.mansoni) data
require(drc)

  but.vals = c(0, 556, 2417, 3906, 5560, 8703)/1000
  gly.vals = c(0, 1506, 3875, 9174, 15062, 26249)/1000
  pen.vals = c(0, 214.8, 535, 1299, 2148, 3762)/1000
  
L.3.fx = function(t, lc50 = lc50, slp = slp){
    1 / (1+exp(slp*log(t / lc50)))
}  
  
  time = seq(0,25,0.1)

#Cercarial (S. mansoni) toxicity ###############
cer<-read.csv('Agrochemical_Review/Response_Fxs/Data/ghaffar2016_cercariae.csv')
cer$conc = cer$conc/1000
  cer$alive = 100*cer$surv
  cer$dead = 100-cer$alive
  cer$total = 100
    
  cer.ctrl = subset(cer, chem == 'control')  
    
    gaf.ctrl.drc.cer = drm(alive/total ~ time_hrs, weights = total, data = cer.ctrl, type = 'binomial', 
                      fct = LL.3(names = c('b', 'd', 'e'),
                                fixed = c(NA, 1, NA)))

#butralin cercariae ###############
  cer.but = subset(cer, chem == 'butralin')
  
  gaf.but.drc.cer = drm(alive/total ~ time_hrs, conc, weights = total, data = cer.but, type = 'binomial', 
                    fct = LL.3(names = c('b', 'd', 'e'),
                               fixed = c(NA, 1, NA)))

#Get estimate of cercariae-hrs for each concentration    
gaf.but.aucs.cer = as.numeric()
  gaf.but.aucs.cer[1] = integrate(f = L.3.fx, lc50 = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], 
                                  lower=0, upper=24)[1]$value  
  
  for(j in 1:length(unique(cer.but$conc))){
    gaf.but.aucs.cer[j+1] = integrate(f = L.3.fx, lc50 = gaf.but.drc.cer$coefficients[j+5], slp = gaf.but.drc.cer$coefficients[j],
                                      lower=0, upper=24)[1]$value  
  }
  
#Compile LL.2 data for functional responses #############
but.cer = data.frame(butralin = but.vals,
                     logbutralin = log(but.vals+1),
                     e = c(summary(gaf.ctrl.drc.cer)$coefficients[2,1], summary(gaf.but.drc.cer)$coefficients[c(6:10),1]),
                     e.se = c(summary(gaf.ctrl.drc.cer)$coefficients[2,2], summary(gaf.but.drc.cer)$coefficients[c(6:10),2]),
                     b = c(summary(gaf.ctrl.drc.cer)$coefficients[1,1], summary(gaf.but.drc.cer)$coefficients[c(1:5),1]),
                     b.se = c(summary(gaf.ctrl.drc.cer)$coefficients[1,2], summary(gaf.but.drc.cer)$coefficients[c(1:5),2]))

#functions fit to L.2 parameters ############
but.cer.lm.e = lm(e ~ butralin, weights = e.se^-1, data = but.cer) #linear response of LC50
  but.cer.pred = function(but){
    predict(but.cer.lm.e, newdata = data.frame(butralin = but), 
            interval = 'confidence', level = 0.95)
  }
    
but.cer.lm.e2 = lm(e ~ logbutralin, weights = e.se^-1, data = but.cer) #log-linear response of LC50
  but.cer.pred2 = function(but){
    predict(but.cer.lm.e2, newdata = data.frame(logbutralin = log(but+1)), 
            interval = 'confidence', level = 0.95)
  }
  
  #AIC(but.cer.lm.e, but.cer.lm.e2)  #exponential is a better fit    

but.cer.lm.b = lm(b ~ butralin, weights = b.se^-1, data = but.cer)   
  but.cer.pred.b = function(but){
    predict(but.cer.lm.b, newdata = data.frame(butralin = but), interval = 'confidence', level = 0.95)
  }
  
auc.but.lin0 = function(He){
  e0 = as.numeric(predict(but.cer.lm.e, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(but.cer.lm.b, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
  auc0
}

piC.ghaf_butr.lin_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(but.cer.lm.e, newdata = data.frame(butralin = Heu), se.fit = TRUE)[1:2])
    b = as.numeric(predict(but.cer.lm.b, newdata = data.frame(butralin = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piC = auc/auc.but.lin0(0) #gaf.but.aucs.cer[1]

  return(piC)
}  #function to estimate AUC 

auc.but.exp0 = function(He){
  e0 = as.numeric(predict(but.cer.lm.e2, newdata = data.frame(logbutralin = log(He+1)), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(but.cer.lm.b, newdata = data.frame(butralin = He), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
  auc0
}

piC.ghaf_butr.exp_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(but.cer.lm.e2, newdata = data.frame(logbutralin = log(Heu+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(but.cer.lm.b, newdata = data.frame(butralin = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    #print(e.use)
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    #print(b.use)
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piC = auc/auc.but.exp0(0) #gaf.but.aucs.cer[1]
  
  return(piC)
}  #function to estimate AUC 

keep.gaf.but.cer = c('piC.ghaf_butr.lin_unc', 'but.cer.lm.e', 'but.cer.lm.b', 'L.3.fx', 'gaf.but.aucs.cer',
                     'piC.ghaf_butr.exp_unc', 'but.cer.lm.e2', 'auc.but.lin0', 'auc.but.exp0')

#glyphosate cercariae ###############
cer.gly = subset(cer, chem == 'glyphosate')
  
  gaf.gly.drc.cer = drm(alive/total ~ time_hrs, conc, weights = total, data = cer.gly, type = 'binomial', 
                    fct = LL.3(names = c('b', 'd', 'e'),
                              fixed = c(NA, 1, NA)))

#Get estimate of cercariae-hrs for each concentration    
gaf.gly.aucs.cer = as.numeric()
  gaf.gly.aucs.cer[1] = integrate(f = L.3.fx, lc50 = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], 
                                  lower=0, upper=24)[1]$value  
  
  for(j in 1:length(unique(cer.gly$conc))){
    gaf.gly.aucs.cer[j+1] = integrate(f = L.3.fx, lc50 = gaf.gly.drc.cer$coefficients[j+5], slp = gaf.gly.drc.cer$coefficients[j],
                                      lower=0, upper=24)[1]$value  
  }
  

#Compile L.2 parameters ######## 
gly.cer = data.frame(glyphosate = gly.vals,
                     logglyphosate = log(gly.vals+1),
                     e = c(summary(gaf.ctrl.drc.cer)$coefficients[2,1], summary(gaf.gly.drc.cer)$coefficients[c(6:10),1]),
                     e.se = c(summary(gaf.ctrl.drc.cer)$coefficients[2,2], summary(gaf.gly.drc.cer)$coefficients[c(6:10),2]),
                     b = c(summary(gaf.ctrl.drc.cer)$coefficients[1,1], summary(gaf.gly.drc.cer)$coefficients[c(1:5),1]),
                     b.se = c(summary(gaf.ctrl.drc.cer)$coefficients[1,2], summary(gaf.gly.drc.cer)$coefficients[c(1:5),2]))

#fit functions to LL.2 parameters ##########  
gly.cer.lm.e = lm(e ~ glyphosate, weights = e.se^-1, data = gly.cer) #linear response of LC50
  gly.cer.pred = function(gly){
    predict(gly.cer.lm.e, newdata = data.frame(glyphosate = gly), 
            interval = 'confidence', level = 0.95)
  }
  
gly.cer.lm.e2 = lm(e ~ logglyphosate, weights = e.se^-1, data = gly.cer) #log-linear response of LC50
  gly.cer.pred2 = function(gly){
    predict(gly.cer.lm.e2, newdata = data.frame(logglyphosate = log(gly+1)), 
            interval = 'confidence', level = 0.95)
  }
  
  # AIC(gly.cer.lm.e, gly.cer.lm.e2)  #Exponential is a better fit    

gly.cer.lm.b = lm(b ~ glyphosate, weights = b.se^-1, data = gly.cer)   
  gly.cer.pred.b = function(gly){
    predict(gly.cer.lm.b, newdata = data.frame(glyphosate = gly), interval = 'confidence', level = 0.95)
  }
  
auc.gly.lin0 = function(He){
    e0 = as.numeric(predict(gly.cer.lm.e, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(gly.cer.lm.b, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
    
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
    auc0
  }
  
  
piC.ghaf_gly.lin_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(gly.cer.lm.e, newdata = data.frame(glyphosate = Heu), se.fit = TRUE)[1:2])
    b = as.numeric(predict(gly.cer.lm.b, newdata = data.frame(glyphosate = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piC = auc/ auc.gly.lin0(0) #gaf.gly.aucs.cer[1]
  
  return(piC)
}  #function to estimate AUC 


auc.gly.exp0 = function(He){
  e0 = as.numeric(predict(gly.cer.lm.e2, newdata = data.frame(logglyphosate = log(He+1)), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(gly.cer.lm.b, newdata = data.frame(glyphosate = He), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
  auc0
}


piC.ghaf_gly.exp_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(gly.cer.lm.e2, newdata = data.frame(logglyphosate = log(Heu+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(gly.cer.lm.b, newdata = data.frame(glyphosate = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piC = auc/auc.gly.exp0(0) #gaf.gly.aucs.cer[1]
  
  return(piC)
}  #function to estimate AUC 

keep.gaf.gly.cer = c('piC.ghaf_gly.lin_unc', 'gly.cer.lm.e', 'gly.cer.lm.b', 'L.3.fx', 'gaf.gly.aucs.cer',
                     'piC.ghaf_gly.exp_unc', 'gly.cer.lm.e2', 'auc.gly.lin0', 'auc.gly.exp0')


#pendimethalin cercariae ###############
cer.pen = subset(cer, chem == 'pendimethalin')

  gaf.pen.drc.cer = drm(alive/total ~ time_hrs, conc, weights = total, data = cer.pen, type = 'binomial', 
                    fct = LL.3(names = c('b', 'd', 'e'),
                               fixed = c(NA, 1, NA)))
  
#Get estimate of cercariae-hrs for each concentration    
gaf.pen.aucs.cer = as.numeric()
gaf.pen.aucs.cer[1] = integrate(f = L.3.fx, lc50 = gaf.ctrl.drc.cer$coefficients[2], slp = gaf.ctrl.drc.cer$coefficients[1], 
                                lower=0, upper=24)[1]$value  

  for(j in 1:length(unique(cer.pen$conc))){
    gaf.pen.aucs.cer[j+1] = integrate(f = L.3.fx, lc50 = gaf.pen.drc.cer$coefficients[j+5], slp = gaf.pen.drc.cer$coefficients[j],
                                      lower=0, upper=24)[1]$value  
  }

#Compile L.2 parameters ######## 
pen.cer = data.frame(pendimethalin = pen.vals,
                     logpendimethalin = log(pen.vals+1),
                     e = c(summary(gaf.ctrl.drc.cer)$coefficients[2,1], summary(gaf.pen.drc.cer)$coefficients[c(6:10),1]),
                     e.se = c(summary(gaf.ctrl.drc.cer)$coefficients[2,2], summary(gaf.pen.drc.cer)$coefficients[c(6:10),2]),
                     b = c(summary(gaf.ctrl.drc.cer)$coefficients[1,1], summary(gaf.pen.drc.cer)$coefficients[c(1:5),1]),
                     b.se = c(summary(gaf.ctrl.drc.cer)$coefficients[1,2], summary(gaf.pen.drc.cer)$coefficients[c(1:5),2]))

#fit functions to LL.2 parameters ##########  
pen.cer.lm.e = lm(e ~ pendimethalin, weights = e.se^-1, data = pen.cer) #linear response of LC50
  pen.cer.pred = function(pen){
    predict(pen.cer.lm.e, newdata = data.frame(pendimethalin = pen), 
            interval = 'confidence', level = 0.95)
  }
  
pen.cer.lm.e2 = lm(e ~ logpendimethalin, weights = e.se^-1, data = pen.cer) #log-linear response of LC50
  pen.cer.pred2 = function(pen){
    predict(pen.cer.lm.e2, newdata = data.frame(logpendimethalin = log(pen+1)), 
            interval = 'confidence', level = 0.95)
  }
  
  #AIC(pen.cer.lm.e, pen.cer.lm.e2)  #Exponential is a better fit    

pen.cer.lm.b = lm(b ~ pendimethalin, weights = b.se^-1, data = pen.cer)   
  pen.cer.pred.b = function(pen){
    predict(pen.cer.lm.b, newdata = data.frame(pendimethalin = pen), interval = 'confidence', level = 0.95)
  }
  
auc.pen.lin0 = function(He){
    e0 = as.numeric(predict(pen.cer.lm.e, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(pen.cer.lm.b, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
    
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
    auc0
  }
  
piC.ghaf_pen.lin_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(pen.cer.lm.e, newdata = data.frame(pendimethalin = Heu), se.fit = TRUE)[1:2])
    b = as.numeric(predict(pen.cer.lm.b, newdata = data.frame(pendimethalin = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piC = auc/auc.pen.lin0(0) #gaf.pen.aucs.cer[1]
  
  return(piC)
}  #function to estimate AUC 
  
auc.pen.exp0 = function(He){
  e0 = as.numeric(predict(pen.cer.lm.e2, newdata = data.frame(logpendimethalin = log(He+1)), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(pen.cer.lm.b, newdata = data.frame(pendimethalin = He), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
  auc0
}


piC.ghaf_pen.exp_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(pen.cer.lm.e2, newdata = data.frame(logpendimethalin = log(Heu+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(pen.cer.lm.b, newdata = data.frame(pendimethalin = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                    rel.tol = 1e-5, stop.on.error = FALSE)[1]$value
    piC = auc/auc.pen.exp0(0) #gaf.pen.aucs.cer[1]
  
  return(piC)
}  #function to estimate AUC 

keep.gaf.pen.cer = c('piC.ghaf_pen.lin_unc', 'pen.cer.lm.e', 'pen.cer.lm.b', 'L.3.fx', 'gaf.pen.aucs.cer',
                     'piC.ghaf_pen.exp_unc', 'pen.cer.lm.e2', 'auc.pen.lin0', 'auc.pen.exp0')


#Final keep vector
keep.gaf.piC = c(keep.gaf.pen.cer, keep.gaf.gly.cer, keep.gaf.but.cer)  