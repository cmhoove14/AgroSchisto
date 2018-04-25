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
#Cercarial mortality (S. mansoni) from Tchounwou 1992 ##################
cerc = read.csv('Agrochemical_Review/Response_Fxs/Data/tchounwou92.csv')
cerc$conc = cerc$conc/1000
  cerc.mal = subset(cerc, chem == 'malathion')
  cerc.mal$dead = round(cerc.mal$dead)
  time = seq(0,25,0.1)

#Tchounwou Data plotted ############
  tch92.piC.mal<-drm(alive/total ~ time_hrs, conc, weights = total, data = cerc.mal, type = 'binomial', 
                     fct = LL.3(names = c('b', 'd', 'e'),
                                fixed = c(NA, 1, NA)))

#Get estimate of cercariae-hrs for each concentration    
tch92.mal.piC.aucs = as.numeric()

for(j in 1:length(unique(cerc.mal$conc))){
  fx = function(t){
    predict(tch92.piC.mal, newdata = data.frame(time_hrs = t, conc = unique(cerc.mal$conc)[j]))
  }
  tch92.mal.piC.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
}

cerc.auc.mal = data.frame(mal = unique(cerc.mal$conc),
                          piC = c(tch92.mal.piC.aucs)/tch92.mal.piC.aucs[1])

#Data frame of LL.2 parameters across maltahion concentrations ##################    
parms.df = data.frame(mal = unique(cerc.mal$conc),
                      logmal = log(unique(cerc.mal$conc)+1),
                      e = summary(tch92.piC.mal)$coefficients[c(7:12), 1],
                      e.se = summary(tch92.piC.mal)$coefficients[c(7:12), 2],
                      b = summary(tch92.piC.mal)$coefficients[c(1:6), 1],
                      b.se = summary(tch92.piC.mal)$coefficients[c(1:6), 2])

  tch92.mal.e.mod = lm(e ~ mal, weights = e.se^-1, data = parms.df) 
    fx.e.mod = function(In){
      predict(tch92.mal.e.mod, newdata = data.frame(mal = In), interval = 'confidence', level = 0.95)
    }

  tch92.mal.e.mod2 = lm(e ~ logmal, weights = e.se^-1, data = parms.df)
    fx.e.mod2 = function(In){
      predict(tch92.mal.e.mod2, newdata = data.frame(logmal = log(In+1)), interval = 'confidence', level = 0.95)
    }

  tch92.mal.b.mod = lm(b ~ mal, weights = b.se^-1, data = parms.df) 
    fx.b.mod = function(In){
      predict(tch92.mal.b.mod, newdata = data.frame(mal = In), interval = 'confidence', level = 0.95)
    }
    
#Function to estimate survival curve as function of malathion conc #####################################
piC.tch92_mal_auc0 = function(In){
    Ins = In/1000
    e0 = as.numeric(predict(tch92.mal.e.mod2, newdata = data.frame(logmal = log(Ins+1)), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(tch92.mal.b.mod, newdata = data.frame(mal = Ins), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value

  return(auc0)
}  


piC.tch92_mal_unc = function(In){
  Ins = In/1000
  e = as.numeric(predict(tch92.mal.e.mod2, newdata = data.frame(logmal = log(Ins+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch92.mal.b.mod, newdata = data.frame(mal = Ins), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                  stop.on.error = FALSE)[1]$value
  piC = auc/piC.tch92_mal_auc0(0)  #tch92.mal.piC.aucs[1]

  return(piC)
}  

keep.tch92.beq = c('tch92.mal.piC.aucs', 'piC.tch92_mal_unc', 'piC.tch92_mal_auc0', 'L.3.fx', 
                   'tch92.mal.e.mod2', 'tch92.mal.b.mod', 'cerc.auc.mal')  

#Qualitative model validation ###############

#function to plot model predictions
pred.fx.plot = function(In, clr){
  e = as.numeric(predict(tch92.mal.e.mod2, newdata = data.frame(logmal = log((In/1000)+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(tch92.mal.b.mod, newdata = data.frame(mal = (In/1000)), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  lines(time, L.3.fx(time, slp = b.use, lc50 = e.use), lty=2, col = clr)
}   