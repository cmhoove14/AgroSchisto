#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
#Data extraction and model fitting to Tantawy 2002 data
require(drc)

L.3.fx = function(t, lc50 = lc50, slp = slp){
  1 / (1+exp(slp*log(t / lc50)))
}
time = c(0:25)


#miracidia toxicity ############
tant<-read.csv('Agrochemical_Review/Response_Fxs/Data/Tantawy2002.csv')

  tant$total = 100
  tant$conc = tant$conc/1000
    mir<-subset(tant, larv == 'miracidia')
    tant.mir.but.dat = subset(mir, chem == 'butachlor')
    tant.mir.fpb.dat = subset(mir, chem == 'fluazifop-p-butyl')
  
#butachlor toxicity to miracidia  ###############
tant.piM.but<-drm(surv/total ~ time_hrs, conc, weights = total,  
                  data = tant.mir.but.dat, type = 'binomial', 
                  fct = LL.3(names = c('b', 'd', 'e'),
                             fixed = c(NA, 1, NA)))

#Get estimate of cercariae-hrs for each concentration    
tant.but.piM.aucs = as.numeric()

for(j in 1:length(unique(tant.mir.but.dat$conc))){
  fx = function(t){
    predict(tant.piM.but, newdata = data.frame(time_hrs = t, conc = unique(tant.mir.but.dat$conc)[j]))
  }
  tant.but.piM.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
}

#Compile LL.2 data for functional responses #############
mir.but = data.frame(e = summary(tant.piM.but)$coefficients[c(6:10),1],
                     e.se = summary(tant.piM.but)$coefficients[c(6:10),2],
                     b = summary(tant.piM.but)$coefficients[c(1:5),1],
                     b.se = summary(tant.piM.but)$coefficients[c(1:5),2],
                     but = c(0,650,1500,4500,6500)/1000,
                     logbut = log(c(0,650,1500,4500,6500)/1000+1))

#fit models to LL.2 parameters across concentration ########
el.but.mir = lm(e ~ but, weights = e.se^-1, data = mir.but) #linear response of LC50
  el.pred.mir = function(but){
    predict(el.but.mir, newdata = data.frame(but = but), 
            interval = 'confidence', level = 0.95)
  }

el.but2.mir = lm(e ~ logbut, weights = e.se^-1, data = mir.but) #log-linear response of LC50
  el.pred2.mir = function(but){
    predict(el.but2.mir, newdata = data.frame(logbut = log(but+1)), 
            interval = 'confidence', level = 0.95)
  }

bl.but.mir = lm(b ~ but, weights = b.se^-1, data = mir.but)   
  bl.pred.mir = function(but){
    predict(bl.but.mir, newdata = data.frame(but = but), interval = 'confidence', level = 0.95)
  }

#develop functions to replicate data ###############  
piM.tant02_but.lin_auc0 = function(He){
      Heu = He/1000
      e0 = as.numeric(predict(el.but.mir, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
      b0 = as.numeric(predict(bl.but.mir, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
      
      e0.use = rnorm(1, e0[1], e0[2])
        while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
      b0.use = rnorm(1, b0[1], b0[2])
        while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
      
      auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value

    return(auc0)
  }  #function to estimate AUC when He = 0
  
piM.tant02_but.lin_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(el.but.mir, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bl.but.mir, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24)[1]$value
    piM = auc / piM.tant02_but.lin_auc0(0) #tant.but.piM.aucs[1]

  return(piM)
}  #function to estimate AUC 

piM.tant02_but.exp_auc0 = function(He){
  Heu = He/1000
  e0 = as.numeric(predict(el.but2.mir, newdata = data.frame(logbut = log(Heu+1)), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(bl.but.mir, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24)[1]$value
  
  return(auc0)
}  #function to estimate AUC when He = 0

piM.tant02_but.exp_unc = function(He){
  Heu = He/1000
  e = as.numeric(predict(el.but2.mir, newdata = data.frame(logbut = log(Heu+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bl.but.mir, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24)[1]$value
  piM = auc/piM.tant02_but.exp_auc0(0)  #tant.but.piM.aucs[1]
  
  return(piM)
}  #function to estimate AUC 

#keep vector for butachlor toxicity to miracidia
  keep.tant.but.piM = c('piM.tant02_but.lin_unc', 'piM.tant02_but.lin_auc0', 'piM.tant02_but.exp_unc', 'piM.tant02_but.exp_auc0', 'el.but.mir', 'el.but2.mir', 'bl.but.mir', 'L.3.fx')
    
#fluazifop-p-butyl toxicity to miracidia ###########
tant.piM.fpb<-drm(surv/total ~ time_hrs, conc, weights = total,  
                  data = tant.mir.fpb.dat, type = 'binomial', 
                  fct = LL.3(names = c('b', 'd', 'e'),
                             fixed = c(NA, 1, NA)))

#Get estimate of cercariae-hrs for each concentration    
  tant.fpb.piM.aucs = as.numeric()

  for(j in 1:length(unique(tant.mir.fpb.dat$conc))){
    fx = function(t){
      predict(tant.piM.fpb, newdata = data.frame(time_hrs = t, conc = unique(tant.mir.fpb.dat$conc)[j]))
    }
    tant.fpb.piM.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }

#Compile LL.2 data for functional responses #############
mir.fpb = data.frame(e = summary(tant.piM.fpb)$coefficients[c(6:10),1],
                     e.se = summary(tant.piM.fpb)$coefficients[c(6:10),2],
                     b = summary(tant.piM.fpb)$coefficients[c(1:5),1],
                     b.se = summary(tant.piM.fpb)$coefficients[c(1:5),2],
                     fpb = c(0, 1760,4500,9000,17600)/1000,
                     logfpb = log(c(0, 1760,4500,9000,17600)/1000+1))
  
#fit models to LL.2 parameters across concentration ########
el.fpb.mir = lm(e ~ fpb, weights = e.se^-1, data = mir.fpb) #linear response of LC50
  el.pred.mir.fpb = function(fpb){
    predict(el.fpb.mir, newdata = data.frame(fpb = fpb), 
            interval = 'confidence', level = 0.95)
  }

el.fpb2.mir = lm(e ~ logfpb, weights = e.se^-1, data = mir.fpb) #log-linear response of LC50
  el.pred2.mir.fpb = function(fpb){
    predict(el.fpb2.mir, newdata = data.frame(logfpb = log(fpb+1)), 
            interval = 'confidence', level = 0.95)
  }

bl.fpb.mir = lm(b ~ fpb, weights = b.se^-1, data = mir.fpb)   
  bl.pred.mir.fpb = function(fpb){
    predict(bl.fpb.mir, newdata = data.frame(fpb = fpb), interval = 'confidence', level = 0.95)
  }

  
#Functions to estimate die off ########  
piM.tant02_fpb.lin_auc0 = function(He){
    Heu = He/1000
    e0 = as.numeric(predict(el.fpb.mir, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(bl.fpb.mir, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e[1], e[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
    
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                    stop.on.error = FALSE)[1]$value

    return(auc0)
  }  #function to estimate AUC 
  
piM.tant02_fpb.lin_unc = function(He){
  Heu = He/1000
  e = as.numeric(predict(el.fpb.mir, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bl.fpb.mir, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                  stop.on.error = FALSE)[1]$value
  
  piM = auc/ piM.tant02_fpb.lin_auc0(0) #tant.fpb.piM.aucs[1]
  
  return(piM)
}  #function to estimate AUC 

piM.tant02_fpb.exp_auc0 = function(He){
  Heu = He/1000
  e0 = as.numeric(predict(el.fpb2.mir, newdata = data.frame(logfpb = log(Heu+1)), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(bl.fpb.mir, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e[1], e[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                  stop.on.error = FALSE)[1]$value
  return(auc0)
}  #function to estimate AUC 

piM.tant02_fpb.exp_unc = function(He){
  Heu = He/1000
  e = as.numeric(predict(el.fpb2.mir, newdata = data.frame(logfpb = log(Heu+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(bl.fpb.mir, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24)[1]$value
  piM = auc/piM.tant02_fpb.exp_auc0(0)  #tant.fpb.piM.aucs[1]

  return(piM)
}  #function to estimate AUC 

#keep vector for fluazifop-p-butyl toxicity to miracidia
  keep.tant.fpb.piM = c('piM.tant02_fpb.lin_unc', 'piM.tant02_fpb.lin_auc0', 'piM.tant02_fpb.exp_unc', 'piM.tant02_fpb.exp_auc0', 'el.fpb.mir', 'el.fpb2.mir', 'bl.fpb.mir', 'L.3.fx')
  
#keep vector ##########
  keep.tantawy.piM = c(keep.tant.but.piM, keep.tant.fpb.piM)    
  