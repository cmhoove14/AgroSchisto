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
mir = read.csv('Agrochemical_Review/Response_Fxs/Data/tchounwou08_miracidia.csv')
  mir$conc = mir$conc/1000
  mir.mal = subset(mir, chem == 'mal')

#Tchounwou Data plotted ############
tch91.piM.mal<-drm(alive/total ~ time_hrs, conc, weights = total, data = mir.mal, type = 'binomial', 
                   fct = LL.3(names = c('b', 'd', 'e'),
                              fixed = c(NA, 1, NA)))

#Get estimate of miracidia-hrs for each concentration    
tch91.mal.pim.aucs = as.numeric()
  
  for(j in 1:length(unique(mir.mal$conc))){
    fx = function(t){
      predict(tch91.piM.mal, newdata = data.frame(time_hrs = t, conc = unique(mir.mal$conc)[j]))
    }
    tch91.mal.pim.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
  }

mir.auc.mal = data.frame(mal = unique(mir.mal$conc),
                         piM = c(tch91.mal.pim.aucs)/tch91.mal.pim.aucs[1])    

#Create data frame with parameter values and malathion concentrations #######################
mirp.df = data.frame(mal = unique(mir.mal$conc),
                     logmal = log(unique(mir.mal$conc)+1),
                     e = summary(tch91.piM.mal)$coefficients[c(7:12), 1],
                     e.se = summary(tch91.piM.mal)$coefficients[c(7:12), 2],
                     b = summary(tch91.piM.mal)$coefficients[c(1:6), 1],
                     b.se = summary(tch91.piM.mal)$coefficients[c(1:6), 2])


tch91.em.mod = lm(e ~ mal, weights = e.se^-1, data = mirp.df) 
  tch91.em.fx = function(mal){
    predict(tch91.em.mod, newdata = data.frame(mal = mal), interval = 'confidence', level = 0.95)
  }

tch91.em.mod2 = lm(e ~ logmal, weights = e.se^-1, data = mirp.df) 
  tch91.em.fx2 = function(mal){
      predict(tch91.em.mod2, newdata = data.frame(logmal = log(mal+1)), interval = 'confidence', level = 0.95)
    }

tch91.bm.mod = lm(b ~ mal, weights = b.se^-1, data = mirp.df) 
  tch91.bm.fx = function(mal){
    predict(tch91.bm.mod, newdata = data.frame(mal = mal), interval = 'confidence', level = 0.95)
  }

#Create function to generate d-r function #####################
piM.tch91_mal_auc0 = function(In){
    Ins = In/1000
    e0 = as.numeric(predict(tch91.em.mod, newdata = data.frame(mal = Ins), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(tch91.bm.mod, newdata = data.frame(mal = Ins), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
    
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value
  
  return(auc0)
}  

piM.tch91_mal_unc = function(In){
  Ins = In/1000
  e = as.numeric(predict(tch91.em.mod, newdata = data.frame(mal = Ins), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch91.bm.mod, newdata = data.frame(mal = Ins), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                  stop.on.error = FALSE)[1]$value
  piM = auc/piM.tch91_mal_auc0(0)  #tch91.mal.pim.aucs[1]
  
  return(piM)
}  
  
keep.tch91.beq = c('tch91.mal.pim.aucs', 'piM.tch91_mal_unc', 'piM.tch91_mal_auc0', 'L.3.fx', 'mirp.df',
                   'tch91.em.mod', 'tch91.bm.mod', 'mir.auc.mal')

#Qualitative model validation ###############
#function to plot model predictions
predm.fx.plot = function(In, clr){
  Ins = In/1000
  e = as.numeric(predict(tch91.em.mod, newdata = data.frame(mal = Ins), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch91.bm.mod, newdata = data.frame(mal = Ins), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
  while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
  while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  lines(time, L.3.fx(time, lc50 = e.use, slp = b.use), lty=2, col = clr)
}
