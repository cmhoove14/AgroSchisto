#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(drc)

L.3.fx = function(t, lc50 = lc50, slp = slp){
  1 / (1+exp(slp*log(t / lc50)))
}
  time = seq(0,25,0.1)

#miracidial mortality (S. mansoni) from Tchounwou 1991 ####################################
mir2 = read.csv('Agrochemical_Review/Response_Fxs/Data/tchounwou91_urea_ammP_miracidial_mortality.csv')
  mir2$conc = mir2$conc/1000
  mir.amm = subset(mir2, chem == 'amm_sulph')
  mir.ure = subset(mir2, chem == 'urea')

#Ammonium Sulphate concentration ############
tch91.piM.amm<-drm(alive/total ~ time_hrs, conc, weights = total, data = mir.amm, type = 'binomial', 
                   fct = LL.3(names = c('b', 'd', 'e'),
                              fixed = c(NA, 1, NA)))

#Get estimate of miracidia-hrs for each concentration    
tch91.amm.pim.aucs = as.numeric()
  
  for(j in 1:length(unique(mir.amm$conc))){
    fx = function(t){
      predict(tch91.piM.amm, newdata = data.frame(time_hrs = t, conc = unique(mir.amm$conc)[j]))
    }
    tch91.amm.pim.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }

mir.auc.amm = data.frame(amm = unique(mir.amm$conc),
                         piM = c(tch91.amm.pim.aucs)/tch91.amm.pim.aucs[1])

#Create data frame with parameter values and ammonium sulphate concentrations #######################
miramm.df = data.frame(amm = unique(mir.amm$conc),
                     logamm = log(unique(mir.amm$conc)+1),
                     e = summary(tch91.piM.amm)$coefficients[c(7:12), 1],
                     e.se = summary(tch91.piM.amm)$coefficients[c(7:12), 2],
                     b = summary(tch91.piM.amm)$coefficients[c(1:6), 1],
                     b.se = summary(tch91.piM.amm)$coefficients[c(1:6), 2])

#linear fit
tch91.e.amm.lin = lm(e ~ amm, weights = e.se^-1, data = miramm.df) 
  tch91.e.amm.lin.pred = function(amm){
    predict(tch91.e.amm.lin, newdata = data.frame(amm = amm), interval = 'confidence', level = 0.95)
  }

#exponential fit
tch91.e.amm.exp = lm(e ~ logamm, weights = e.se^-1, data = miramm.df) 
  tch91.e.amm.exp.pred = function(amm){
      predict(tch91.e.amm.exp, newdata = data.frame(logamm = log(amm+1)), interval = 'confidence', level = 0.95)
    }

tch91.b.amm = lm(b ~ amm, weights = b.se^-1, data = miramm.df) 
  tch91.b.amm.pred = function(amm){
    predict(tch91.b.amm, newdata = data.frame(amm = amm), interval = 'confidence', level = 0.95)
  }

#Create function to generate d-r function #####################
piM.tch91_amm_exp_auc0 = function(In){
    Ins = In/1000
    e0 = as.numeric(predict(tch91.e.amm.exp, newdata = data.frame(logamm = log(Ins+1)), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(tch91.b.amm, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
    
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value

  return(auc0)
}  

piM.tch91_amm_unc = function(In){
    Ins = In/1000
  e = as.numeric(predict(tch91.e.amm.exp, newdata = data.frame(logamm = log(Ins+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch91.b.amm, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, e[1], e[2])
  
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                  stop.on.error = FALSE)[1]$value
  piM = auc/piM.tch91_amm_exp_auc0(0) #tch91.amm.pim.aucs[1]

  return(piM)
} 

piM.tch91_amm_lin_auc0 = function(In){
  Ins = In/1000
  e0 = as.numeric(predict(tch91.e.amm.lin, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(tch91.b.amm, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                   stop.on.error = FALSE)[1]$value
  
  return(auc0)
}  

piM.tch91_amm_lin = function(In){
    Ins = In/1000
    e = as.numeric(predict(tch91.e.amm.lin, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
    b = as.numeric(predict(tch91.b.amm, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, b[1], b[2])
    
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                    stop.on.error = FALSE)[1]$value
    piM = auc/piM.tch91_amm_lin_auc0(0) #tch91.amm.pim.aucs[1]
  
  return(piM)
}  
  
keep.tch91.amm = c('tch91.amm.pim.aucs', 'piM.tch91_amm_lin', 'piM.tch91_amm_exp_auc0', 'piM.tch91_amm_unc', 'piM.tch91_amm_lin_auc0', 'L.3.fx', 'miramm.df', 'tch91.e.amm.exp', 'tch91.e.amm.lin', 'tch91.b.amm', 'mir.auc.amm')

#Qualitative model validation ###############
#function to plot model predictions
predm.amm.plot = function(In, clr){
  Ins = In/1000
  e = as.numeric(predict(tch91.e.amm.lin, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch91.b.amm, newdata = data.frame(amm = Ins), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  lines(time, L.3.fx(time, slp = b.use, lc50 = e.use), lty=2, col = clr)
}   
  
#Urea concentration ############
tch91.piM.ure<-drm(alive/total ~ time_hrs, conc, weights = total, data = mir.ure, type = 'binomial', 
                   fct = LL.3(names = c('b', 'd', 'e'),
                              fixed = c(NA, 1, NA)))

#Get estimate of miracidia-hrs for each concentration    
tch91.ure.pim.aucs = as.numeric()

  for(j in 1:length(unique(mir.ure$conc))){
    fx = function(t){
      predict(tch91.piM.ure, newdata = data.frame(time_hrs = t, conc = unique(mir.ure$conc)[j]))
    }
    tch91.ure.pim.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }

mir.auc.ure = data.frame(ure = unique(mir.ure$conc),
                         piM = c(tch91.ure.pim.aucs)/tch91.ure.pim.aucs[1])

#Create data frame with parameter values and urea concentrations #######################
mirure.df = data.frame(ure = unique(mir.ure$conc),
                       logure = log(unique(mir.ure$conc)+1),
                       e = summary(tch91.piM.ure)$coefficients[c(7:12), 1],
                       e.se = summary(tch91.piM.ure)$coefficients[c(7:12), 2],
                       b = summary(tch91.piM.ure)$coefficients[c(1:6), 1],
                       b.se = summary(tch91.piM.ure)$coefficients[c(1:6), 2])

#linear fit
tch91.e.ure.lin = lm(e ~ ure, weights = e.se^-1, data = mirure.df) 
  tch91.e.ure.lin.pred = function(ure){
    predict(tch91.e.ure.lin, newdata = data.frame(ure = ure), interval = 'confidence', level = 0.95)
  }

#exponential fit
tch91.e.ure.exp = lm(e ~ logure, weights = e.se^-1, data = mirure.df) 
  tch91.e.ure.exp.pred = function(ure){
    predict(tch91.e.ure.exp, newdata = data.frame(logure = log(ure+1)), interval = 'confidence', level = 0.95)
  }

tch91.b.ure = lm(b ~ ure, weights = b.se^-1, data = mirure.df) 
  tch91.b.ure.pred = function(ure){
    predict(tch91.b.ure, newdata = data.frame(ure = ure), interval = 'confidence', level = 0.95)
  }

#Create function to generate d-r function #####################

piM.tch91_ure_lin_auc0 = function(In){
  Ins = In/1000
  e0 = as.numeric(predict(tch91.e.ure.lin, newdata = data.frame(ure = Ins), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(tch91.b.ure, newdata = data.frame(ure = Ins), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                   stop.on.error = FALSE)[1]$value
  
  return(auc0)
}  

piM.tch91_ure_unc = function(In){
  Ins = In/1000
  e = as.numeric(predict(tch91.e.ure.lin, newdata = data.frame(ure = Ins), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch91.b.ure, newdata = data.frame(ure = Ins), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                  stop.on.error = FALSE)[1]$value
  piM = auc/piM.tch91_ure_lin_auc0(0) #tch91.amm.pim.aucs[1]
  
  return(piM)
}  

keep.tch91.ure = c('tch91.ure.pim.aucs', 'tch91.e.ure.exp', 'piM.tch91_ure_unc', 'piM.tch91_ure_lin_auc0', 'L.3.fx', 'mirure.df', 'tch91.e.ure.lin', 'tch91.b.ure', 'mir.auc.ure')

#Qualitative model validation ###############
#function to plot model predictions
predm.ure.plot = function(In, clr){
  Ins = In/1000
  e = as.numeric(predict(tch91.e.ure.lin, newdata = data.frame(ure = Ins), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch91.b.ure, newdata = data.frame(ure = Ins), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  lines(time, L.3.fx(time, slp = b.use, lc50 = e.use), lty=2, col = clr)
}   

#Final keep vector #########
keep.tch91.Fe = c(keep.tch91.amm, keep.tch91.ure)