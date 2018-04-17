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

#cercarial toxicity ############
tant<-read.csv('Agrochemical_Review/Response_Fxs/Data/Tantawy2002.csv')

  tant$total = 100
  tant$conc = tant$conc/1000
  cerc<-subset(tant, larv == 'cercariae')
    tant.cerc.but.dat = subset(cerc, chem == 'butachlor')
    tant.cerc.fpb.dat = subset(cerc, chem == 'fluazifop-p-butyl')

#butachlor toxicity to cercariae  ###############
tant.piC.but<-drm(surv/total ~ time_hrs, conc, weights = total,  data = tant.cerc.but.dat, type = 'binomial', 
                  fct = LL.3(names = c('b', 'd', 'e'),
                            fixed = c(NA, 1, NA)))

#Get estimate of cercariae-hrs for each concentration    
tant.but.pic.aucs = as.numeric()

  for(j in 1:length(unique(tant.cerc.but.dat$conc))){
    fx = function(t){
      predict(tant.piC.but, newdata = data.frame(time_hrs = t, conc = unique(tant.cerc.but.dat$conc)[j]))
    }
    tant.but.pic.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }

#compile butachlor data for function ###############
cerc.but = data.frame(e = summary(tant.piC.but)$coefficients[c(6:10),1],
                      e.se = summary(tant.piC.but)$coefficients[c(6:10),2],
                      b = summary(tant.piC.but)$coefficients[c(1:5),1],
                      b.se = summary(tant.piC.but)$coefficients[c(1:5),2],
                      but = unique(tant.cerc.but.dat$conc),
                      logbut = log(unique(tant.cerc.but.dat$conc)+1))
    
#fit models to L.2 parameters across concentration ########
  el.but = lm(e ~ but, weights = e.se^-1, data = cerc.but) #linear response of LC50
    el.pred = function(but){
      predict(el.but, newdata = data.frame(but = but), 
              interval = 'confidence', level = 0.95)
    }
    
  el.but2 = lm(e ~ logbut, weights = e.se^-1, data = cerc.but) #log-linear response of LC50
    el.pred2 = function(but){
      predict(el.but2, newdata = data.frame(logbut = log(but+1)), 
              interval = 'confidence', level = 0.95)
    }

  bl.but = lm(b ~ but, weights = b.se^-1, data = cerc.but)   
    bl.pred = function(but){
      predict(bl.but, newdata = data.frame(but = but), interval = 'confidence', level = 0.95)
    }
  
# Fit functions to predict change in AUC  #######       
piC.tant02_but.lin_auc0 = function(He){
  Heu = He/1000
  e0 = as.numeric(predict(el.but, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(bl.but, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                  stop.on.error = FALSE)[1]$value
  return(auc0)
}

piC.tant02_but.lin_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(el.but, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bl.but, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
        
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                    stop.on.error = FALSE)[1]$value
    piC = auc/piC.tant02_but.lin_auc0(0)  #tant.but.pic.aucs[1]
  
    return(piC)
} 

piC.tant02_but.exp_auc0 = function(He){
  Heu = He/1000
  e0 = as.numeric(predict(el.but2, newdata = data.frame(logbut = log(Heu+1)), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(bl.but, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                   stop.on.error = FALSE)[1]$value
  return(auc0)
}

piC.tant02_but.exp_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(el.but2, newdata = data.frame(logbut = log(Heu+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bl.but, newdata = data.frame(but = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                    stop.on.error = FALSE)[1]$value
    piC = auc/piC.tant02_but.exp_auc0(0)  #tant.but.pic.aucs[1]
  
  return(piC)
}  #Parameter estimate

keep.tant.but.pic = c('piC.tant02_but.lin_unc', 'piC.tant02_but.lin_auc0', 'piC.tant02_but.exp_unc', 'piC.tant02_but.exp_auc0', 'bl.but', 'el.but2', 'el.but', 'L.3.fx')      
    
#fluazifop-p-butyl toxicity to cercariae ###########
  tant.piC.fpb<-drm(surv/total ~ time_hrs, conc, weights = total,  data = tant.cerc.fpb.dat, type = 'binomial', fct = LL.3(names = c('b', 'd', 'e'),
             fixed = c(NA, 1, NA)))

#Get estimate of cercariae-hrs for each concentration    
  tant.fpb.pic.aucs = as.numeric()
  
  for(j in 1:length(unique(tant.cerc.fpb.dat$conc))){
    fx = function(t){
      predict(tant.piC.fpb, newdata = data.frame(time_hrs = t, conc = unique(tant.cerc.fpb.dat$conc)[j]))
    }
    tant.fpb.pic.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }
  
  
#compile fluazifop-p-butyl data for function ###############
  cerc.fpb = data.frame(e = summary(tant.piC.fpb)$coefficients[c(6:10),1],
                        e.se = summary(tant.piC.fpb)$coefficients[c(6:10),2],
                        b = summary(tant.piC.fpb)$coefficients[c(1:5),1],
                        b.se = summary(tant.piC.fpb)$coefficients[c(1:5),2],
                        fpb = unique(tant.cerc.fpb.dat$conc),
                        logfpb = log(unique(tant.cerc.fpb.dat$conc)+1))
  
#fit models to L.2 parameters across concentration ########
el.fpb = lm(e ~ fpb, weights = e.se^-1, data = cerc.fpb) #linear response of LC50
  el.pred.fpb = function(fpb){
    predict(el.fpb, newdata = data.frame(fpb = fpb), 
            interval = 'confidence', level = 0.95)
  }
    
el.fpb2 = lm(e ~ logfpb, weights = e.se^-1, data = cerc.fpb) #log-linear response of LC50
  el.pred.fpb2 = function(fpb){
    predict(el.fpb2, newdata = data.frame(logfpb = log(fpb+1)), 
            interval = 'confidence', level = 0.95)
  }
    
bl.fpb = lm(b ~ fpb, weights = b.se^-1, data = cerc.fpb)   
  bl.pred.fpb = function(fpb){
    predict(bl.fpb, newdata = data.frame(fpb = fpb), interval = 'confidence', level = 0.95)
  }
    
#Fit functions to model changes in AUC
piC.tant02_fpb.exp_auc0 = function(He){
    Heu = He/1000
    e0 = as.numeric(predict(el.fpb2, newdata = data.frame(logfpb = log(Heu+1)), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(bl.fpb, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], bo[2])
    
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value
    
  return(auc0)
}  #Parameter estimate with exponential function

piC.tant02_fpb.exp_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(el.fpb2, newdata = data.frame(logfpb = log(Heu+1)), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bl.fpb, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24)[1]$value
    piC = auc/piC.tant02_fpb.exp_auc0(0)  #tant.fpb.pic.aucs[1]

  return(piC)
}  #Parameter estimate with exponential function
  
piC.tant02_fpb.lin_auc0 = function(He){
  Heu = He/1000
  e0 = as.numeric(predict(el.fpb, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(bl.fpb, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], bo[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                   stop.on.error = FALSE)[1]$value
  
  return(auc0)
}  #Parameter estimate with linear function

piC.tant02_fpb.lin_unc = function(He){
    Heu = He/1000
    e = as.numeric(predict(el.fpb, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bl.fpb, newdata = data.frame(fpb = Heu), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                    stop.on.error = FALSE)[1]$value
    piC = auc/piC.tant02_fpb.lin_auc0(0)  #tant.fpb.pic.aucs[1]
  
   return(piC)
}  #Parameter estimate with linear function
    
#keep vector ########    
keep.tant.fpb.pic = c('piC.tant02_fpb.lin_unc', 'piC.tant02_fpb.lin_auc0', 'piC.tant02_fpb.exp_unc', 'piC.tant02_fpb.exp_auc0', 'bl.fpb',
                      'el.fpb2', 'el.fpb', 'L.3.fx') 
    
keep.tantawy.piC = c(keep.tant.but.pic, keep.tant.fpb.pic)    
    