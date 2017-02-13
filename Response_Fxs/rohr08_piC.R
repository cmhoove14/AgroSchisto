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

ll4 = function(hi,lo,slp,lc,x){
  lo + ((hi-lo)/(1+exp(slp*(log(x)-lc))))
}

rohr = read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/rohrarial Mortality/rohr2008.csv')
  time = c(0:25)

#Plot control survival curve, estimate time-dependent die-off function, and get auc from 0-24 hours ########
plot(x = rohr$time_hrs[rohr$chem == 'control'], y = rohr$surv[rohr$chem == 'control']/100,
     xlab = 'time (hrs)', ylab = 'prop dead', pch = 16)

rohr.ctrl = drm(alive/total ~ time_hrs, total, data = rohr, 
                type = 'binomial', fct = LL2.2())

lines(time, ll4(1,0,coef(rohr.ctrl)[1], coef(rohr.ctrl)[2], time), lty=2)

#Atrazine ########
rohr.atr = subset(rohr, chem == 'atrazine')
points(x = rohr.atr$time_hrs, y = rohr.atr$surv/100,
       col = 'gold', pch = 16)

  rohr.atrmod = drm(alive/total ~ time_hrs, total, data = rohr.atr, 
                    type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,coef(rohr.atrmod)[1], coef(rohr.atrmod)[2], time), lty=2, col ='gold')

#Malathion ########
  rohr.mal = subset(rohr, chem == 'malathion')
  points(x = rohr.mal$time_hrs, y = rohr.mal$surv/100,
         col = 'red', pch = 16)
  
  rohr.malmod = drm(alive/total ~ time_hrs, total, data = rohr.mal, 
                    type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,coef(rohr.malmod)[1], coef(rohr.malmod)[2], time), lty=2, col ='red')
  
#carbaryl ########
  rohr.car = subset(rohr, chem == 'carbaryl')
  points(x = rohr.car$time_hrs, y = rohr.car$surv/100,
         col = 'purple', pch = 16)
  
  rohr.carmod = drm(alive/total ~ time_hrs, total, data = rohr.car, 
                    type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,coef(rohr.carmod)[1], coef(rohr.carmod)[2], time), lty=2, col ='purple')
  
#glyphosate ########
rohr.gly = subset(rohr, chem == 'glyphosate')
  points(x = rohr.gly$time_hrs, y = rohr.gly$surv/100,
         col = 'green', pch = 16)
  
  rohr.glymod = drm(alive/total ~ time_hrs, total, data = rohr.gly, 
                    type = 'binomial', fct = LL2.2())
  
  lines(time, ll4(1,0,coef(rohr.glymod)[1], coef(rohr.glymod)[2], time), lty=2, col ='green')
  
#Summaries ##########
  legend('bottomleft', legend = c('control', 'atrazine', 'carbaryl', 'glyphosate', 'malathion'),
         col = c('black', 'gold', 'red', 'green','purple'), pch = 16, cex=0.8)
  title(main = 'Rohr 2008 cercarial mortality')
  
#Atrazine "model" (yes, I know there are only two data points)  ###########
#  atr.dat = data.frame('conc' = c(0,201),
#                       'auc24' = c(rel.auc.rohr[1], rel.auc.rohr[2]))
#  
#    atr.m<-nls(auc24 ~ exp(-b*conc), data = atr.dat,
#               start = list(b=0.0001))
#      summary(atr.m)
#    
#    atr.con = c(0:201)
#    
#   pi_C_atr_rohr08 = function(He){
#     exp(-summary(atr.m)$parameters[1]*(He)) 
#    } 
#    
#    plot(atr.dat$conc, atr.dat$auc24, pch = 16, xlab = 'Atrazine concentration', ylab = '24-hr AUC', ylim = c(0,1))
#      lines(atr.con, pi_C_atr_rohr08(atr.con), lty=2, col='red')
#      
#Carbaryl "model" (yes, I know there are only two data points)  ###########
#  carb.dat = data.frame('conc' = c(0,33.5),
#                       'auc24' = c(rel.auc.rohr[1], rel.auc.rohr[3]))
#  
#    carb.m<-nls(auc24 ~ exp(-b*conc), data = carb.dat,
#               start = list(b=0.01))
#      summary(carb.m)
#    
#    carb.con = c(0:40)
#    
#    pi_C_carb_rohr08 = function(In){
#      exp(-summary(carb.m)$parameters[1]*(In)) 
#    } 
#    
#    plot(carb.dat$conc, carb.dat$auc24, pch = 16, xlab = 'Carbaryl concentration', ylab = '24-hr AUC', ylim = c(0,1))
#      lines(carb.con, pi_C_carb_rohr08(carb.con), lty=2, col='red')
#    
#glyphosate "model" (yes, I know there are only two data points)  ###########
#  gly.dat = data.frame('conc' = c(0,3700),
#                       'auc24' = c(rel.auc.rohr[1], rel.auc.rohr[4]))
# 
#  gly.m<-nls(auc24 ~ exp(-b*conc), data = gly.dat,
#              start = list(b=0.001))
#    summary(gly.m)
#  
#  gly.con = c(0:4000)
#  
#  pi_C_gly_rohr08 = function(He){
#    exp(-summary(gly.m)$parameters[1]*(He)) 
#  } 
#  
#  plot(gly.dat$conc, gly.dat$auc24, pch = 16, xlab = 'glyphosate concentration', ylab = '24-hr AUC', ylim = c(0,1))
#    lines(gly.con, pi_C_gly_rohr08(gly.con), lty=2, col='red')
#    
#Malathion "model" (yes, I know there are only two data points)  ###########
#    
#  mal.dat = data.frame('conc' = c(0,9.6),
#                      'auc24' = c(rel.auc.rohr[1], rel.auc.rohr[5]))
#  
#  mal.m<-nls(auc24 ~ exp(-b*conc), data = mal.dat,
#             start = list(b=0.001))
#    summary(mal.m)
#  
#  mal.con = seq(0,10,0.1)
#  
#  pi_C_mal_rohr08 = function(In){
#    exp(-summary(mal.m)$parameters[1]*(In)) 
#  } 
#  
#  plot(mal.dat$conc, mal.dat$auc24, pch = 16, xlab = 'Malathion concentration', ylab = '24-hr AUC', ylim = c(0,1))
#    lines(mal.con, pi_C_mal_rohr08(mal.con), lty=2, col='red')