#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

# Miracidia and cercariae survival of Echinistoma paraensei exposed to roundup
# Roundup original has 41.0% active ingredient (gkyphosate) therefore we multiply concentrations by this proportion to get effective dose of glyphosate
source("Agrochemical_Review/Response_Fxs/tchounwou92_malathion_cercariae_fit.R")

rm(list = setdiff(ls(), c("parms.df", "tch92.piC.mal")))

require(drc)

L.3.fx = function(t, lc50 = lc50, slp = slp){
  1 / (1+exp(slp*log(t / lc50)))
}
  time = seq(0,25,0.1)

monte_dat <- read.csv("Agrochemical_Review/Response_Fxs/Data/monte2016_glyphosate_larvae.csv")

#Cercarial mortality from Monte 2016 ##################
cerc = subset(monte_dat, larvae == "cercariae" & conc < 900 & conc > 15) #Highest and lowest dose groups have no/instant mortality which makes drm upset
cerc$conc = cerc$conc*0.41 #Multiply by proportion active ingredient to get Glyphosate concentration
cerc$dead = round(cerc$dead) #round mortality to nearest whole number

#DRC model ############
  mont16.piC.gly<-drm(alive/total ~ time_hrs, conc, weights = total, data = cerc, type = 'binomial', 
                      fct = LL.3(names = c('b', 'd', 'e'),
                                 fixed = c(NA, 1, NA)))

#Get estimate of cercariae-hrs for each concentration    
mont16.gly.piC.aucs = as.numeric()

for(j in 1:length(unique(cerc$conc))){
  fx = function(t){
    predict(mont16.piC.gly, newdata = data.frame(time_hrs = t, conc = unique(cerc$conc)[j]))
  }
  mont16.gly.piC.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
}

cerc.auc.gly = data.frame(gly = unique(cerc$conc),
                          piC = c(mont16.gly.piC.aucs)/mont16.gly.piC.aucs[1])

#Data frame of LL.2 parameters across glyphosate concentrations ##################    
mont16_pars = data.frame(gly = unique(cerc$conc),
                         loggly = log(unique(cerc$conc)+1),
                         e = summary(mont16.piC.gly)$coefficients[c(6:10), 1],
                         e.se = summary(mont16.piC.gly)$coefficients[c(6:10), 2],
                         b = summary(mont16.piC.gly)$coefficients[c(1:5), 1],
                         b.se = summary(mont16.piC.gly)$coefficients[c(1:5), 2])

#Since no survival over observation period for this experiment, use control data from tchounwou study to get control dose group. Survival in the period of observation form this study is similar, tchounwou just observes for longer which allows us to fit a function to it
mont16_pars[nrow(mont16_pars)+1,] <- as.numeric(as.numeric(parms.df[1,]))

  mont16.gly.e.mod = lm(e ~ gly, weights = e.se^-1, data = mont16_pars) 
    fx.e.mod = function(He){
      predict(mont16.gly.e.mod, newdata = data.frame(gly = He), interval = 'confidence', level = 0.95)
    }

  mont16.gly.e.mod2 = lm(e ~ loggly, weights = e.se^-1, data = mont16_pars)
    fx.e.mod2 = function(He){
      predict(mont16.gly.e.mod2, newdata = data.frame(loggly = log(He+1)), interval = 'confidence', level = 0.95)
    }

  mont16.gly.b.mod = lm(b ~ gly, weights = b.se^-1, data = mont16_pars) 
    fx.b.mod = function(He){
      predict(mont16.gly.b.mod, newdata = data.frame(gly = He), interval = 'confidence', level = 0.95)
    }
    
#Function to estimate survival curve as function of glyphosate conc #####################################
piC.mont16_gly_auc0 = function(He){
    heu = He/1000
    e0 = as.numeric(predict(mont16.gly.e.mod2, newdata = data.frame(loggly = log(heu+1)), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(mont16.gly.b.mod, newdata = data.frame(gly = heu), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value

  return(auc0)
}  


piC.mont16_gly_unc = function(He){
  heu = He/1000
  e = as.numeric(predict(mont16.gly.e.mod2, newdata = data.frame(loggly = log(heu+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(mont16.gly.b.mod, newdata = data.frame(gly = heu), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                  stop.on.error = FALSE)[1]$value
  piC = auc/piC.mont16_gly_auc0(0)

  return(piC)
}  

keep.monte2016 <- c('mont16.gly.piC.aucs', 'piC.mont16_gly_unc', 'piC.mont16_gly_auc0', 'L.3.fx', 
                   'mont16.gly.e.mod2', 'mont16.gly.b.mod', 'cerc.auc.gly')  

#Qualitative model validation ###############
#function to plot model predictions
monte16_pred_ts = function(He, clr){
  e = as.numeric(predict(mont16.gly.e.mod2, newdata = data.frame(loggly = log((He/1000)+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(mont16.gly.b.mod, newdata = data.frame(gly = (He/1000)), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  lines(time, L.3.fx(time, slp = b.use, lc50 = e.use), lty=2, col = clr)
}   

#Miracidial mortality ######
  mir = subset(monte_dat, larvae == "miracidia") #Highest and lowest dose groups have no/instant mortality which makes drm upset
  mir$conc = mir$conc*0.41 #Multiply by proportion active ingredient to get Glyphosate concentration
  mir$dead = round(mir$dead) #round mortality to nearest whole number

#DRC model ############
  mont16.piM.gly<-drm(alive/total ~ time_hrs, conc, weights = total, data = mir, type = 'binomial', 
                      fct = LL.3(names = c('b', 'd', 'e'),
                                 fixed = c(NA, 1, NA)))

#Get estimate of cercariae-hrs for each concentration    
mont16.gly.piM.aucs = as.numeric()

for(j in 1:length(unique(mir$conc))){
  fx = function(t){
    predict(mont16.piM.gly, newdata = data.frame(time_hrs = t, conc = unique(mir$conc)[j]))
  }
  mont16.gly.piM.aucs[j] = integrate(f = fx, lower=0, upper=24)[1]$value  
}

mir.auc.gly = data.frame(gly = unique(mir$conc),
                         piM = c(mont16.gly.piM.aucs)/mont16.gly.piM.aucs[1])

#Data frame of LL.2 parameters across glyphosate concentrations ##################    
mont16_pars_mir = data.frame(gly = unique(mir$conc),
                             loggly = log(unique(mir$conc)+1),
                             e = summary(mont16.piM.gly)$coefficients[c(7:12), 1],
                             e.se = summary(mont16.piM.gly)$coefficients[c(7:12), 2],
                             b = summary(mont16.piM.gly)$coefficients[c(1:6), 1],
                             b.se = summary(mont16.piM.gly)$coefficients[c(1:6), 2])

  mont16.gly.e.mod.mir = lm(e ~ gly, weights = e.se^-1, data = mont16_pars_mir) 
    fx.e.mod.mir = function(He){
      predict(mont16.gly.e.mod.mir, newdata = data.frame(gly = He), interval = 'confidence', level = 0.95)
    }

  mont16.gly.e.mod2.mir = lm(e ~ loggly, weights = e.se^-1, data = mont16_pars_mir)
    fx.e.mod2.mir = function(He){
      predict(mont16.gly.e.mod2.mir, newdata = data.frame(loggly = log(He+1)), interval = 'confidence', level = 0.95)
    }

  mont16.gly.b.mod.mir = lm(b ~ gly, weights = b.se^-1, data = mont16_pars_mir) 
    fx.b.mod.mir = function(He){
      predict(mont16.gly.b.mod.mir, newdata = data.frame(gly = He), interval = 'confidence', level = 0.95)
    }
    
#Function to estimate survival curve as function of glyphosate conc #####################################
piM.mont16_gly_auc0 = function(He){
    heu = He/1000
    e0 = as.numeric(predict(mont16.gly.e.mod2.mir, newdata = data.frame(loggly = log(heu+1)), se.fit = TRUE)[1:2])
    b0 = as.numeric(predict(mont16.gly.b.mod.mir, newdata = data.frame(gly = heu), se.fit = TRUE)[1:2])
    
    e0.use = rnorm(1, e0[1], e0[2])
      while(e0.use < 0) e0.use = rnorm(1, e0[1], e0[2])
    b0.use = rnorm(1, b0[1], b0[2])
      while(b0.use < 0) b0.use = rnorm(1, b0[1], b0[2])
    auc0 = integrate(L.3.fx, lc50 = e0.use, slp = b0.use, lower=0, upper=24,
                     stop.on.error = FALSE)[1]$value

  return(auc0)
}  


piM.mont16_gly_unc = function(He){
  heu = He/1000
  e = as.numeric(predict(mont16.gly.e.mod2.mir, newdata = data.frame(loggly = log(heu+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(mont16.gly.b.mod.mir, newdata = data.frame(gly = heu), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  auc = integrate(L.3.fx, lc50 = e.use, slp = b.use, lower=0, upper=24,
                  stop.on.error = FALSE)[1]$value
  piC = auc/piM.mont16_gly_auc0(0)

  return(piC)
}  

keep.monte2016 = c(keep.monte2016, 'mont16.gly.piM.aucs', 'piM.mont16_gly_unc', 'piM.mont16_gly_auc0', 'L.3.fx', 
                   'mont16.gly.e.mod2.mir', 'mont16.gly.b.mod.mir', 'mir.auc.gly')  

#Qualitative model validation ###############

#function to plot model predictions
monte16_pred_ts_mir = function(He, clr){
  e = as.numeric(predict(mont16.gly.e.mod2.mir, newdata = data.frame(loggly = log((He/1000)+1)), se.fit = TRUE)[1:2])
  b = as.numeric(predict(mont16.gly.b.mod.mir, newdata = data.frame(gly = (He/1000)), se.fit = TRUE)[1:2])
  
  e.use = rnorm(1, e[1], e[2])
    while(e.use < 0) e.use = rnorm(1, e[1], e[2])
  b.use = rnorm(1, b[1], b[2])
    while(b.use < 0) b.use = rnorm(1, b[1], b[2])
  
  lines(time, L.3.fx(time, slp = b.use, lc50 = e.use), lty=2, col = clr)
}   