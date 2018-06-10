#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Ghaffar 2016 SNAIL (B. alexandrina) data
#Snail (30 B. alexandrina 6-8mm shell width) toxicity ##########
#Butralin ##############
but.dat = data.frame(lcs = c(0, 10, 25, 50, 90),
                     butralin = c(556, 2417, 3906, 5560, 8703))
  but.dat.no0 = subset(but.dat, lcs!=0)
  but.dat.no0$mort = but.dat.no0$lcs/100
  but.dat.no0$probit = qnorm(but.dat.no0$lcs/100)
  but.dat.no0$ppm = but.dat.no0$butralin/1000
  but.dat.no0$ppmlog10 = log10(but.dat.no0$butralin/1000)

#Create function based on LC values
  gaf_but_lin <- lm(probit~ppmlog10, data = but.dat.no0)
    gaf.b.lin <- coef(gaf_but_lin)[1]
    gaf.se.b.lin <- summary(gaf_but_lin)$coef[2,2]
      gaf.b.lin.up <- confint(gaf_but_lin)[1,1]
      gaf.b.lin.lo <- confint(gaf_but_lin)[1,2]
    gaf.m.lin <- coef(gaf_but_lin)[2]
  
muNq_but_Gaf16_lin <- function(He, b = gaf.b.lin){
  mun = pnorm(b + log10(He/1000)*gaf.m.lin)
  return(mun)
}  

muNq_butr_gaf16_uncertainty = function(He){
    heu = (He/1000)
    lc50 = rnorm(1, gaf.b.lin, gaf.se.b.lin)
    mun = pnorm(lc50 + log10(He/1000)*gaf.m.lin)
  
    return(mun)
}
  
#keep vector
keep.gaf.but = c('muNq_butr_gaf16_uncertainty',
                 'gaf.b.lin', 'gaf.se.b.lin', 'gaf.m.lin')    
        
#Glyphosate #############
gly.dat = data.frame(lcs = c(0, 0, 10, 25, 50, 90),
                     glyphosate = c(0, 1506, 3875, 9174, 15062, 26249))

  gly.dat.no0 = subset(gly.dat, lcs!=0)
  gly.dat.no0$probit = qnorm(gly.dat.no0$lcs/100)
  gly.dat.no0$ppm = gly.dat.no0$glyphosate/1000
  gly.dat.no0$ppmlog10 = log10(gly.dat.no0$glyphosate/1000)
  
#Create function based on LC values
  gaf_gly_lin <- lm(probit~ppmlog10, data = gly.dat.no0)
    gaf.b.lin.gly <- coef(gaf_gly_lin)[1]
    gaf.se.b.lin.gly <- summary(gaf_gly_lin)$coef[2,2]
      gaf.b.lin.up.gly <- confint(gaf_gly_lin)[1,1]
      gaf.b.lin.lo.gly <- confint(gaf_gly_lin)[1,2]
    gaf.m.lin.gly <- coef(gaf_gly_lin)[2]
  
muNq_gly_Gaf16_lin <- function(He, b = gaf.b.lin.gly){
  mun = pnorm(b + log10(He/1000)*gaf.m.lin.gly)
  return(mun)
}  

muNq_gly_gaf16_uncertainty = function(He){
  heu = (He/1000)
  lc50 = rnorm(1, gaf.b.lin.gly, gaf.se.b.lin.gly)
  mun = pnorm(lc50 + log10(He/1000)*gaf.m.lin.gly)
    
  return(mun)
}

keep.gaf.gly = c('muNq_gly_gaf16_uncertainty', 
                 'gaf.b.lin.gly', 'gaf.se.b.lin.gly', 'gaf.m.lin.gly')    
    
#Pendimethalin ##############
pen.dat = data.frame(lcs = c(0, 0, 10, 25, 50, 90),
                     pendimethalin = c(0, 214.8, 535, 1299, 2148, 3762))
    
pen.dat.no0 = subset(pen.dat, lcs!=0)
  pen.dat.no0$probit = qnorm(pen.dat.no0$lcs/100)
  pen.dat.no0$ppm = pen.dat.no0$pendimethalin/1000
  pen.dat.no0$ppmlog10 = log10(pen.dat.no0$pendimethalin/1000)
    
#Create function based on LC values
  gaf_pen_lin <- lm(probit~ppmlog10, data = pen.dat.no0)
    gaf.b.lin.pen <- coef(gaf_pen_lin)[1]
    gaf.se.b.lin.pen <- summary(gaf_pen_lin)$coef[2,2]
      gaf.b.lin.up.pen <- confint(gaf_pen_lin)[1,1]
      gaf.b.lin.lo.pen <- confint(gaf_pen_lin)[1,2]
    gaf.m.lin.pen <- coef(gaf_pen_lin)[2]
  
muNq_pen_Gaf16_lin <- function(He, b = gaf.b.lin.pen){
  mun = pnorm(b + log10(He/1000)*gaf.m.lin.pen)
  return(mun)
}  


muNq_pen_gaf16_uncertainty = function(He){
  heu = (He/1000)
  lc50 = rnorm(1, gaf.b.lin.pen, gaf.se.b.lin.pen)
  mun = pnorm(lc50 + log10(He/1000)*gaf.m.lin.pen)
    
  return(mun)
}

keep.gaf.pen = c('muNq_pen_gaf16_uncertainty',
                 'gaf.b.lin.pen', 'gaf.se.b.lin.pen', 'gaf.m.lin.pen')    
  

#Reproductive toxicity #########
  #Paper reports reproduction as R0: the summed product of live snails/week * eggs/snail/week produced
  #we just want reduction in eggs produced as the model already takes additional mortality 
  #into account; therefore we normalize reported R0s by relative survival between treatment groups
  #to get relative estimate of eggs/snail/week

#Longitudinal surival data from table 3
  gaf16.ad<-read.csv('Agrochemical_Review/Response_Fxs/Data/ghaffar2016_adult.csv')
    sn.wk = as.numeric()  #snails/week estimates for each treatment
  for(i in 1:length(unique(gaf16.ad$conc))){
    sn.wk[i] = sum(gaf16.ad$prop_surv[gaf16.ad$conc == unique(gaf16.ad$conc)[i]]) #proportion surviving over study period
  }  
  
  htch.wk = c(0.99, 0.92, 0.82, 0.24,
              0.89, 0.58, 0.46,
              0.9, 0.3, 0.05)   #vector of hatching proportions (hatches/egg) from table4

#vector of concentrations and reproduction measured as approximate hatchlings/snail/week
#hatchlings/snail/week = R0/live snails/week = eggs/snail/week * hatching proportion *live weeks = hatchlings/snail over study period
  gafrep = data.frame('but.conc'=c(0, 0.556, 2.41, 3.906),
                      'gly.conc'=c(0, 1.506, 3.87, 9.17),
                      'pen.conc'=c(0, 0.214, 0.535, 1.299),
                  #Estimates of hatchlings/snail as R0/live snails/week * weeks alive
                      'but.rep'=c((44.231/sn.wk[1])*htch.wk[1]*8, (7.05/sn.wk[2])*htch.wk[2]*8, 
                                  (4.52/sn.wk[3])*htch.wk[3]*6, (4.04/sn.wk[4])*htch.wk[4]*4),
                      'gly.rep'=c((44.231/sn.wk[1])*htch.wk[1]*8, (4.87/sn.wk[5])*htch.wk[5]*8, 
                                  (4.20/sn.wk[6])*htch.wk[6]*6, (4.23/sn.wk[7])*htch.wk[7]*4),
                      'pen.rep'=c((44.231/sn.wk[1])*htch.wk[1]*8, (5.18/sn.wk[8])*htch.wk[8]*8, 
                                  (4.75/sn.wk[9])*htch.wk[9]*5, (4.27/sn.wk[10])*htch.wk[10]*4))
  
  gafrep_ref <- gafrep[1,4]  #reference reproductive rate to scale to 0-1
  
#Butralin ##########
  but_gafrep <- gafrep[c(1:4),c(1,4)] 
  
  but_repro_mod= drm(but.rep ~ but.conc, data = but_gafrep, type = 'continuous',
                      fct = LL.3(names = c('b', 'd', 'e'),
                                 fixed = c(NA, max(but_gafrep$but.rep), NA)))

  but.r0.pred<-function(He){
    predict(but_repro_mod, data.frame(conc = He/1000), interval = 'confidence', level = 0.95)
  }  
    
  par.tricks.but.gaf = c(coef(but_repro_mod), 'd' = max(but_gafrep$but.rep))[c(1,3,2)]
    
  fNq.butr.fx.uncertainty<-function(He){
    fn <- rdrm(nosim = 1, fct = LL.3(), mpar = par.tricks.but.gaf, yerror = 'rnorm', xerror = He/1000,
               ypar = c(0, predict(but_repro_mod, data.frame(but.conc = He/1000), se.fit = T)[2]))$y / gafrep_ref
    
    if(fn < 0){
      return(0)
    } else {
      return(fn)
    }
  }
  
   keep.gaf.but = c(keep.gaf.but, 'fNq.butr.fx.uncertainty', 'par.tricks.but.gaf', 'gafrep_ref', "but_repro_mod")

#Glyphosate #########
  gly_gafrep <- gafrep[c(1:4),c(2,5)]

  gly_repro_mod= drm(gly.rep ~ gly.conc, data = gly_gafrep, type = 'continuous',
                      fct = LL.3(names = c('b', 'd', 'e'),
                                 fixed = c(NA, max(gly_gafrep$gly.rep), NA)))

  gly.r0.pred<-function(He){
    predict(gly_repro_mod, data.frame(conc = He/1000), interval = 'confidence', level = 0.95)
  }  
    
  par.tricks.gly.gaf = c(coef(gly_repro_mod), 'd' = max(gly_gafrep$gly.rep))[c(1,3,2)]
    
  fNq.gly.fx.uncertainty<-function(He){
    fn <- rdrm(nosim = 1, fct = LL.3(), mpar = par.tricks.gly.gaf, yerror = 'rnorm', xerror = He/1000,
               ypar = c(0, predict(gly_repro_mod, data.frame(gly.conc = He/1000), se.fit = T)[2]))$y / gafrep_ref
    
    if(fn < 0){
      return(0)
    } else {
      return(fn)
    }
  }
  
   keep.gaf.gly = c(keep.gaf.gly, 'fNq.gly.fx.uncertainty', 'par.tricks.gly.gaf', 'gafrep_ref', "gly_repro_mod")

#Pendimethalin #########
  pen_gafrep <- gafrep[c(1:4),c(3,6)]
  
  pen_repro_mod= drm(pen.rep ~ pen.conc, data = pen_gafrep, type = 'continuous',
                      fct = LL.3(names = c('b', 'd', 'e'),
                                 fixed = c(NA, max(pen_gafrep$pen.rep), NA)))

  pen.r0.pred<-function(He){
    predict(pen_repro_mod, data.frame(conc = He/1000), interval = 'confidence', level = 0.95)
  }  
    
  par.tricks.pen.gaf = c(coef(pen_repro_mod), 'd' = max(pen_gafrep$pen.rep))[c(1,3,2)]
    
  fNq.pen.fx.uncertainty<-function(He){
    fn <- rdrm(nosim = 1, fct = LL.3(), mpar = par.tricks.pen.gaf, yerror = 'rnorm', xerror = He/1000,
               ypar = c(0, predict(pen_repro_mod, data.frame(pen.conc = He/1000), se.fit = T)[2]))$y / gafrep_ref
    
    if(fn < 0){
      return(0)
    } else {
      return(fn)
    }
  }
  
   keep.gaf.pen = c(keep.gaf.pen, 'fNq.pen.fx.uncertainty', 'par.tricks.pen.gaf', 'gafrep_ref', "pen_repro_mod")

 
  