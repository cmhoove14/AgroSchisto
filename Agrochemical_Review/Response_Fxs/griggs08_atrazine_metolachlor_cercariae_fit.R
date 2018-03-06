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

cerc.g = read.csv('Agrochemical_Review/Response_Fxs/Data/Griggs2008.csv')
cerc.g$prop_surv = cerc.g$alive / cerc.g$total
  time = seq(0,25,0.1)
  
cerc.g0 = subset(cerc.g, chem != 'control')

#estimate time-dependent die-off function, and get auc from 0-24 hours ########
grg.mod = drm(prop_surv ~ time_hrs, conc, weights = total, data = cerc.g0, type = 'binomial', 
                 fct = LL.3(names = c('b', 'd', 'e'), fixed = c(NA, 1, NA)))

#Get estimate of cercariae-hrs for each concentration    
grg08_atr_aucs = as.numeric()
    
  for(j in 1:length(unique(cerc.g0$conc))){
    fx = function(t){
      predict(grg.mod, newdata = data.frame(time_hrs = t, conc = unique(cerc.g0$conc)[j]))
    }
    grg08_atr_aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }

#Create data frame with parameter values and atrazine concentrations #######################
grgc.df = data.frame(atr = c(0,15,100),
                     logatr = log(c(0,15,100)+1),
                     e = c(summary(grg.mod)$coefficients[c(4:6),1]),
                     e.se = c(summary(grg.mod)$coefficients[c(4:6),2]),
                     b = c(summary(grg.mod)$coefficients[c(1:3),1]),
                     b.se = c(summary(grg.mod)$coefficients[c(1:3),2]))
    
#parameters as function of atrazine  
  eg.mod = lm(e ~ atr, weights = e.se^-1, data = grgc.df) 
  eg.mod2 = lm(e ~ logatr, weights = e.se^-1, data = grgc.df)  
    AIC(eg.mod, eg.mod2) #exponential fits way better
  bg.mod = lm(b ~ atr, weights = b.se^-1, data = grgc.df) 
  
#Create function to generate d-r function with linear fit to lc50 parameter#####################
auc.grg.lin.atr0 = function(He){
  e0 = as.numeric(predict(eg.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(bg.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                   stop.on.error = FALSE)[1]$value
  auc0
}

piC.grg08_atr_unc = function(He){
    e = as.numeric(predict(eg.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    b = as.numeric(predict(bg.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
    
    e.use = rnorm(1, e[1], e[2])
      while(e.use < 0) e.use = rnorm(1, e[1], e[2])
    b.use = rnorm(1, b[1], b[2])
      while(b.use < 0) b.use = rnorm(1, e[1], e[2])
    auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, stop.on.error = FALSE)[1]$value
    piC = auc/auc.grg.lin.atr0(0) 
    
  return(piC)
  }  
  
#Create function to generate d-r function with exponential fit to lc50 parameter#####################
auc.grg.exp.atr0 = function(He){
  e0 = as.numeric(predict(eg.mod2, newdata = data.frame(logatr = log(He+1)), se.fit = TRUE)[1:2])
  b0 = as.numeric(predict(bg.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
  
  e0.use = rnorm(1, e0[1], e0[2])
    while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
  b0.use = rnorm(1, b0[1], b0[2])
    while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
  auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                   stop.on.error = FALSE)[1]$value
  auc0
}

piC.grg08_atr_unc2 = function(He){
      e = as.numeric(predict(eg.mod2, newdata = data.frame(logatr = log(He+1)), se.fit = TRUE)[1:2])
      b = as.numeric(predict(bg.mod, newdata = data.frame(atr = He), se.fit = TRUE)[1:2])
      
      e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
      b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
      auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24, 
                      stop.on.error = FALSE)[1]$value
      piC = auc/auc.grg.exp.atr0(0) 
      
    return(piC)
    
  }  
  
keep.grg08 = c('L.3.fx', 'grg08_atr_aucs', 'piC.grg08_atr_unc', 'auc.grg.lin.atr0', 'auc.grg.exp.atr0', 'piC.grg08_atr_unc2', 'grgc.df',
               'eg.mod', 'eg.mod2', 'bg.mod')
