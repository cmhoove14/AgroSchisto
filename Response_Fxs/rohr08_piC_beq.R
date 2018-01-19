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

rohr = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/cercarial Mortality/rohr2008.csv')
  time = c(0:25)

#Plot control survival curve, estimate time-dependent die-off function, and get auc from 0-24 hours ########
plot(x = rohr$time_hrs[rohr$chem == 'control'], y = rohr$surv[rohr$chem == 'control']/100,
     xlab = 'time (hrs)', ylab = 'prop dead', pch = 16)

rohr.ctrl = drm(alive/total ~ time_hrs, total, data = rohr, 
                type = 'binomial', fct = LL.4(names = c('b', 'c', 'd', 'e'),
                                                        fixed = c(NA, 0, 1, NA)))

lines(time, L.3.fx(time, coef(rohr.ctrl)[2], coef(rohr.ctrl)[1]), lty=2)

#Atrazine ########
rohr.atr = subset(rohr, chem == 'atrazine')
points(x = rohr.atr$time_hrs, y = rohr.atr$surv/100,
       col = 'gold', pch = 16)

  rohr.atrmod = drm(alive/total ~ time_hrs, total, data = rohr.atr, 
                    type = 'binomial', fct = LL.4(names = c('b', 'c', 'd', 'e'),
                                                 fixed = c(NA, 0, 1, NA)))
  
  lines(time, L.3.fx(time, coef(rohr.atrmod)[2], coef(rohr.atrmod)[1]), lty=2, col = 'gold')
  
#Malathion ########
  rohr.mal = subset(rohr, chem == 'malathion')
  points(x = rohr.mal$time_hrs, y = rohr.mal$surv/100,
         col = 2, pch = 16)
  
  rohr.malmod = drm(alive/total ~ time_hrs, total, data = rohr.mal, 
                    type = 'binomial', fct = LL.4(names = c('b', 'c', 'd', 'e'),
                                                 fixed = c(NA, 0, 1, NA)))
  
  lines(time, L.3.fx(time, coef(rohr.malmod)[2], coef(rohr.malmod)[1]), lty=2, col = 2)
  
#carbaryl ########
  rohr.car = subset(rohr, chem == 'carbaryl')
  points(x = rohr.car$time_hrs, y = rohr.car$surv/100,
         col = 'purple', pch = 16)
  
  rohr.carmod = drm(alive/total ~ time_hrs, total, data = rohr.car, 
                    type = 'binomial', fct = LL.4(names = c('b', 'c', 'd', 'e'),
                                                 fixed = c(NA, 0, 1, NA)))
  
  lines(time, L.3.fx(time, coef(rohr.carmod)[2], coef(rohr.carmod)[1]), lty=2, col = 'purple')
  
#glyphosate ########
rohr.gly = subset(rohr, chem == 'glyphosate')
  points(x = rohr.gly$time_hrs, y = rohr.gly$surv/100,
         col = 3, pch = 16)
  
  rohr.glymod = drm(alive/total ~ time_hrs, total, data = rohr.gly, 
                    type = 'binomial', fct = LL.4(names = c('b', 'c', 'd', 'e'),
                                                 fixed = c(NA, 0, 1, NA)))
  
  lines(time, L.3.fx(time, coef(rohr.glymod)[2], coef(rohr.glymod)[1]), lty=2, col = 3)
  
#Summaries ##########
  legend('bottomleft', legend = c('control', 'atrazine', 'malathion', 'glyphosate', 'carbaryl'),
         col = c('black', 'gold', 'red', 'green','purple'), pch = 16, cex=0.8, bty = 'n')
  title(main = 'Rohr 2008 cercarial mortality')
  
rohr.fin = data.frame(chem = c('control', 'atrazine', 'carbaryl', 'glyphosate', 'malathion'),
                      conc = c(0, 201, 33.5, 3700, 9.6),
                      e = c(coef(rohr.ctrl)[2], coef(rohr.atrmod)[2], coef(rohr.carmod)[2],
                            coef(rohr.glymod)[2], coef(rohr.malmod)[2]),
                      e.se = c(summary(rohr.ctrl)$coefficients[2,2], summary(rohr.atrmod)$coefficients[2,2],
                               summary(rohr.carmod)$coefficients[2,2], summary(rohr.glymod)$coefficients[2,2],
                               summary(rohr.malmod)$coefficients[2,2]),
                      b = c(coef(rohr.ctrl)[1], coef(rohr.atrmod)[1], coef(rohr.carmod)[1],
                            coef(rohr.glymod)[1], coef(rohr.malmod)[1]),
                      b.se = c(summary(rohr.ctrl)$coefficients[1,2], summary(rohr.atrmod)$coefficients[1,2],
                               summary(rohr.carmod)$coefficients[1,2], summary(rohr.glymod)$coefficients[1,2],
                               summary(rohr.malmod)$coefficients[1,2]))