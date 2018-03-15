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
require(drc)
require(LW1949)
source("Agrochemical_Review/Models/litchfield_wilcoxon_get_b1_from_slope.R")

#Snail (30 B. alexandrina 6-8mm shell width) toxicity ##########
#Butralin ##############
but.dat = data.frame(lcs = c(0, 10, 25, 50, 90),
                     butralin = c(556, 2417, 3906, 5560, 8703))
  but.dat.no0 = subset(but.dat, lcs!=0)
  but.dat.no0$mort = but.dat.no0$lcs/100
  but.dat.no0$probit = qnorm(but.dat.no0$lcs/100)
  but.dat.no0$ppm = but.dat.no0$butralin/1000
  but.dat.no0$ppmlog10 = log10(but.dat.no0$butralin/1000)
  but.dat.no0$ppmlog = log(but.dat.no0$butralin/1000)

  plot(but.dat$butralin/1000, but.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,13),
       xlab = 'Butralin (ppm)', ylab = 'snail mortality rate')
  
  but.gaf.mod <- lm(mort ~ ppm, data = but.dat.no0)
    abline(but.gaf.mod, lty = 2)
    but.gaf.mod.lc84 = ((84/100) - coef(but.gaf.mod)[1])/coef(but.gaf.mod)[2]
    but.gaf.mod.lc50 = ((50/100) - coef(but.gaf.mod)[1])/coef(but.gaf.mod)[2]
    but.gaf.mod.lc16 = ((16/100) - coef(but.gaf.mod)[1])/coef(but.gaf.mod)[2]
  
  slp.manual = (but.gaf.mod.lc84/but.gaf.mod.lc50 + but.gaf.mod.lc50/but.gaf.mod.lc16)/2
  
  lc50.gaf.but.report = 5.56
  slp.gaf.but.report = 1.093 #reported slp of 1.093 is unreasonable, assume they didn't transform back from log scale
  b1.gaf.but = get_b1(slp.gaf.but.report)
  #get standard error from reported 95% CIs of lc50
    se.lc50.gaf.but = mean(c(log10(8.34/lc50.gaf.but.report), log10(lc50.gaf.but.report/3.7))) / 1.96
    
  mu_Nq_butr_gaf16_uncertainty = function(He){
    heu = (He/1000)
    lc50 = 10^(rnorm(1, log10(lc50.gaf.but.report), se.lc50.gaf.but))
    mun = pnorm(b1.gaf.but * log10(heu/lc50))

    return(mun)
  }
  
  plot(but.dat$butralin, but.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,13000),
       xlab = 'Butralin (ppb)', ylab = 'snail mortality rate', 
       main = 'D-R function based on reported values')
    segments(x0 = 3700, x1 = 8340, y0 = 0.5, y1 = 0.5, lty = 1, col = 1)
    segments(x0 = 6220, x1 = 12800, y0 = 0.9, y1 = 0.9, lty = 1, col = 1)
  
    points(seq(0,13000,50), sapply(seq(0,13000,50), mu_Nq_butr_gaf16_uncertainty, simplify = T), 
           pch = 5, col = 4, cex = 0.5)
    
#keep vector
keep.gaf.but = c('mu_Nq_butr_gaf16_uncertainty',
                 'lc50.gaf.but.report', 'se.lc50.gaf.but', 'slp.gaf.but.report')    
        
#Glyphosate #############
gly.dat = data.frame(lcs = c(0, 0, 10, 25, 50, 90),
                     glyphosate = c(0, 1506, 3875, 9174, 15062, 26249))

  gly.dat.no0 = subset(gly.dat, lcs!=0)
  gly.dat.no0$probit = qnorm(gly.dat.no0$lcs/100, mean = 5, sd = 1)
  gly.dat.no0$ppm = gly.dat.no0$glyphosate/1000
  gly.dat.no0$ppmlog10 = log10(gly.dat.no0$glyphosate/1000)
  
  lc50.gaf.gly.report = 15.062
  slp.gaf.gly.report = 0.335
  b1.gaf.gly = get_b1(exp(slp.gaf.gly.report))
  
  se.lc50.gaf.gly = log10(16.57/lc50.gaf.gly.report)  / 1.96
  
  mu_Nq_gly_gaf16_uncertainty = function(He){
      heu = He/1000
      lc50 = exp(rnorm(1, log(lc50.gaf.gly.report), se.lc50.gaf.gly))
      mun = pnorm(b1.gaf.gly * log(heu / lc50))

      return(mun)
    }

  plot(gly.dat$glyphosate, gly.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,33000),
       xlab = 'Glyphosate (ppb)', ylab = 'snail mortality rate',
       main = 'D-R function based on reported values')
    segments(x0 = 9130, x1 = 16570, y0 = 0.5, y1 = 0.5, lty = 1, col = 1)
    segments(x0 = 23870, x1 = 28900, y0 = 0.9, y1 = 0.9, lty = 1, col = 1) 
        
    points(seq(0,33000,100), sapply(seq(0,33000,100), mu_Nq_gly_gaf16_uncertainty, simplify = T), 
           pch = 5, col = 4, cex = 0.5)
    
keep.gaf.gly = c('mu_Nq_gly_gaf16_uncertainty', 
                 'lc50.gaf.gly.report', 'se.lc50.gaf.gly', 'slp.gaf.gly.report')    
    
#Pendimethalin ##############
pen.dat = data.frame(lcs = c(0, 0, 10, 25, 50, 90),
                     pendimethalin = c(0, 214.8, 535, 1299, 2148, 3762))
    
pen.dat.no0 = subset(pen.dat, lcs!=0)
  pen.dat.no0$probit = qnorm(pen.dat.no0$lcs/100, mean = 5, sd = 1)
  pen.dat.no0$ppm = pen.dat.no0$pendimethalin/1000
  pen.dat.no0$ppmlog10 = log10(pen.dat.no0$pendimethalin/1000)
    
  lc50.gaf.pen.report = 2.148
  slp.gaf.pen.report = 1.820 #non-transformed slope function here is the only one that makes sense
  b1.gaf.pen = get_b1(slp.gaf.pen.report)
  
  se.lc50.gaf.pen = mean(c(log(3.22/lc50.gaf.pen.report), log(lc50.gaf.pen.report/1.43))) / 1.96
  
  mu_Nq_pen_gaf16_uncertainty = function(He){
    heu = He/1000
    lc50 = exp(rnorm(1, log(lc50.gaf.pen.report), se.lc50.gaf.pen))
    mun = pnorm(slp.gaf.pen.report * log(heu / lc50))

    return(mun)
  }
  
  plot(pen.dat$pendimethalin, pen.dat$lcs/100, pch = 16, ylim = c(0,1), xlim = c(0,7000),
       xlab = 'pendimethalin (ppb)', ylab = 'snail mortality rate',
       main = 'D-R function based on reported values')
    segments(x0 = 1430, x1 = 3220, y0 = 0.5, y1 = 0.5, lty = 1, col = 1)
    segments(x0 = 2400, x1 = 6420, y0 = 0.9, y1 = 0.9, lty = 1, col = 1)
  
    points(seq(0,7500,25), sapply(seq(0,7500,25), mu_Nq_pen_gaf16_uncertainty, simplify = T), 
           pch = 5, col = 4, cex = 0.5)

keep.gaf.pen = c('mu_Nq_pen_gaf16_uncertainty',
                 'lc50.gaf.pen.report', 'se.lc50.gaf.pen', 'slp.gaf.pen.report')    
  
#Reductions in adult snail (B. alexandrina) reproduction ###########  
  #Paper reports reproduction as R0: the summed product of live snails/week * eggs/snail/week produced
  #we just want reduction in eggs produced as the model already takes additional mortality 
  #into account; therefore we normalize reported R0s by relative survival between treatment groups
  #to get relative estimate of eggs/snail/week
gaf16.ad<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/ghaffar2016_adult.csv')
  sn.wk = as.numeric()  #snails/week estimates for each treatment
for(i in 1:length(unique(gaf16.ad$conc))){
  sn.wk[i] = sum(gaf16.ad$prop_surv[gaf16.ad$conc == unique(gaf16.ad$conc)[i]] * 30) #proportion surviving * 30 snails per treatment
}  

htch.wk = c(0.99, 0.92, 0.82, 0.24,
            0.89, 0.58, 0.46,
            0.9, 0.3, 0.05)   #vector of hatching proportions (hatches/egg)

#vector of concentrations and reproduction measured as approximate hatchlings/snail/week
#hatchlings/snail/week = R0/live snails/week = eggs/snail/week * hatching proportion = hatchlings/snail/week
  gafrep = data.frame('but.conc'=but.dat$butralin[1:4]/1000,
                      'gly.conc'=gly.dat$glyphosate[1:4]/1000,
                      'pen.conc'=pen.dat$pendimethalin[1:4]/1000,
                      'but.rep'=c((44.231/sn.wk[1])*htch.wk[1], (7.05/sn.wk[2])*htch.wk[2], 
                                  (4.52/sn.wk[3])*htch.wk[3], (4.04/sn.wk[4])*htch.wk[4]),
                      'gly.rep'=c((44.231/sn.wk[1])*htch.wk[1], (4.87/sn.wk[5])*htch.wk[5], 
                                  (4.20/sn.wk[6])*htch.wk[6], (4.23/sn.wk[7])*htch.wk[7]),
                      'pen.rep'=c((44.231/sn.wk[1])*htch.wk[1], (5.18/sn.wk[8])*htch.wk[8], 
                                  (4.75/sn.wk[9])*htch.wk[9], (4.27/sn.wk[10])*htch.wk[10]))
  
    gafrep$but.rep = round(gafrep$but.rep, digits = 4)
    gafrep$gly.rep = round(gafrep$gly.rep, digits = 4)
    gafrep$pen.rep = round(gafrep$pen.rep, digits = 4)
    
#Reductions from butralin ##########
plot(gafrep$but.conc, gafrep$but.rep / gafrep$but.rep[1], pch = 16, ylim = c(0,1),
     xlab = 'Butralin (ppm)', ylab = 'relative reproduction rate')
    
  but.r0 = drm(but.rep ~ but.conc, data = gafrep, type = 'continuous',
               fct = LL.3(names = c("b", "d", "e"),
                          fixed = c(NA, gafrep$but.rep[1], NA))) 
    summary(but.r0)
    
  but.r0.pred = function(He){
    predict(but.r0, newdata = data.frame(but.conc = He), interval = 'confidence', level = 0.95)
  }
  
    lines(seq(0,4.500,.010), sapply(seq(0,4.500,.010), but.r0.pred, simplify = T)[1,] / gafrep$but.rep[1],
          lty = 2, col = 2) 
    lines(seq(0,4.500,.010), sapply(seq(0,4.500,.010), but.r0.pred, simplify = T)[2,] / gafrep$but.rep[1],
          lty = 3, col = 2) 
    lines(seq(0,4.500,.010), sapply(seq(0,4.500,.010), but.r0.pred, simplify = T)[3,] / gafrep$but.rep[1],
          lty = 3, col = 2) 
  
 fN.butr.fx.uncertainty = function(He){
  # if(He == 0) fn = 1 else{
     heu = He/1000
     init = predict(but.r0, newdata = data.frame(but.conc = heu), se.fit = T)
     fn = rnorm(1, init[1], init[2]) / gafrep$but.rep[1]
   #}
   
   #if(fn < 0) fn = 0
   #if(fn > 1) fn = 1
     
   return(fn)
 } #normalized to 1
 
 points(seq(0,4.500,.01), 
        sapply(seq(0,4500,10), fN.butr.fx.uncertainty, simplify = T),
        pch = 5, col = 4, cex = 0.5) 

 keep.gaf.but.fn = c('fN.butr.fx.uncertainty', 'but.r0', 'gafrep')
 
#Reductions from glyphosate ##########
 plot(gafrep$gly.conc, gafrep$gly.rep/gafrep$gly.rep[1], pch = 16, ylim = c(0,1),
      xlab = 'glyphosate (ppm)', ylab = 'relative reproduction rate')
 
   gly.r0 = drm(gly.rep ~ gly.conc, data = gafrep, type = 'continuous',
                fct = LL.4(names = c("b", "c", "d", "e"),
                          fixed = c(NA, 0.01, gafrep$gly.rep[1], NA)))  
   summary(gly.r0)
   
 gly.r0.pred = function(He){
   predict(gly.r0, newdata = data.frame(gly.conc = He), interval = 'confidence', level = 0.95)
 }
 
   lines(c(seq(0,1,0.01), seq(1.1,10,0.1)), sapply(c(seq(0,1,0.01), seq(1.1,10,0.1)), 
                                                 gly.r0.pred, simplify = T)[1,]/gafrep$but.rep[1],
         lty = 2, col = 2) 
   lines(c(seq(0,1,0.01), seq(1.1,10,0.1)), sapply(c(seq(0,1,0.01), seq(1.1,10,0.1)), 
                                                 gly.r0.pred, simplify = T)[2,]/gafrep$but.rep[1],
         lty = 3, col = 2) 
   lines(c(seq(0,1,0.01), seq(1.1,10,0.1)), sapply(c(seq(0,1,0.01), seq(1.1,10,0.1)), 
                                                 gly.r0.pred, simplify = T)[3,]/gafrep$but.rep[1],
         lty = 3, col = 2) 
 
   fN.gly.fx.uncertainty = function(He){
     #if(He == 0) fn = 1 else{
       heu = He/1000
       init = predict(gly.r0, newdata = data.frame(gly.conc = heu), se.fit = T)
       fn = rnorm(1, init[1], init[2]) / gafrep$but.rep[1]
     #}
     #if(fn < 0) fn = 0
     #if(fn > 1) fn = 1
     
     return(fn)
   } #normalized to 1
 
     points(seq(0, 10, 0.1), 
            sapply(seq(0, 10000,100), fN.gly.fx.uncertainty, simplify = T),
            pch = 5, col = 4, cex = 0.5) 
     
     keep.gaf.gly.fn = c('fN.gly.fx.uncertainty', 'gly.r0', 'gafrep')
     
#Reductions from pendimethalin ##########
plot(gafrep$pen.conc, gafrep$pen.rep/gafrep$but.rep[1], pch = 16, ylim = c(0,1),
     xlab = 'pendimethalin (ppb)', ylab = 'relative reproduction rate')
     
  pen.r0 = drm(pen.rep ~ pen.conc, data = gafrep, type = 'continuous',
               fct = LL.3(names = c("b", "d", "e"),
                         fixed = c(NA, gafrep$pen.rep[1], NA)))  
    summary(pen.r0)
     
  pen.r0.pred = function(He){
      predict(pen.r0, newdata = data.frame(pen.conc = He), interval = 'confidence', level = 0.95)
  }
     
     lines(seq(0,1.5,0.01), sapply(seq(0,1.5,0.01), pen.r0.pred, simplify = T)[1,]/gafrep$but.rep[1],
           lty = 2, col = 2) 
     lines(seq(0,1.5,0.01), sapply(seq(0,1.5,0.01), pen.r0.pred, simplify = T)[2,]/gafrep$but.rep[1],
           lty = 3, col = 2) 
     lines(seq(0,1.5,0.01), sapply(seq(0,1.5,0.01), pen.r0.pred, simplify = T)[3,]/gafrep$but.rep[1],
           lty = 3, col = 2) 
     
  fN.pen.fx.uncertainty = function(He){
    #if(He == 0) fn = 1 else{
      heu = He/1000
      init = predict(pen.r0, newdata = data.frame(pen.conc = heu), se.fit = T)
      fn = rnorm(1, init[1], init[2]) / gafrep$but.rep[1]
    #}
    #if(fn > 1) fn = 1
    #if(fn < 0) fn = 0
    
    return(fn)  
  } #normalized to 1
     
     points(seq(0,1.5,0.01), 
            sapply(seq(0,1500, 10), fN.pen.fx.uncertainty, simplify = T),
            pch = 5, col = 4, cex = 0.5) #plot with reference back to raw control value
  
    keep.gaf.pen.fn = c('fN.pen.fx.uncertainty', 'pen.r0', 'gafrep') 
    
keep.gaf.all = c(keep.gaf.but, keep.gaf.but.fn, keep.gaf.gly, keep.gaf.gly.fn, keep.gaf.pen, keep.gaf.pen.fn)    