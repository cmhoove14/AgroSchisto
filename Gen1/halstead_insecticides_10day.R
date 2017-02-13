require(ggplot2)
require(drc)

#Insecticide info from Halstead et al chemosphere paper ################
#calculate decay rates (k) from hydrolysis half lives for each chemical in Halstead 2015 (from table S1 and pmep.cce.cornell.edu) 
mal.k = -log(0.5)/6.2         #from table s1 and in agreement of "less than 1 week in raw river water" from Cornell
chlor.k = -log(0.5)/25.5      #from table S1; within the range of reported half life from cornell
terb.k = -log(0.5)/6.5        #from table s1; in agreement with Cornell estimate of 5.5 days at pH of 7; degrades into formaldehyde
lamcy.k = -log(0.5)/1         #very fast; 0 in table S1; "Not expected to be prevalent in surface waters" according to cornell website
esfen.k = -log(0.5)/10        #cornell 4-15 days half life in water
perm.k = -log(0.5)/2          #cornell "half life of less than 2.5 days"

#median EECs from Halstead 2015
med.mal = 0.778
med.chlor = 5.810
med.terb = 1.435
med.lamcy = 0.649
med.esfen = 0.311
med.perm = 1.420

#suggested application intervals for each agrochemical from Halstead 2015
mal.days = c(1,6)             #two applications 5 days apart
chlor.days = c(1,11,21)       #three applications 10 days apart
terb.days = c(1)              #single application
lamcy.days = seq(1, 61, by=4) #16 applications 4 days apart
esfen.days = seq(1, 21, by=5) #5 applications 5 days apart
perm.days = seq(1, 36, by=5)  #8 applications 5 days apart

#Insecticide toxicity to crustaceans from Halstead et al; 4-day mortality endpoints ################
data<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/Cray.LC50.2012.csv')

#Malathion ********************************************************************************************
mal<-subset(data, Chem == 'Mal')
mal.sum = data.frame(mal.c = unique(mal$Conc),
                     mal.total = 5,
                     mal.d = 0,
                     mort = 0)

for(i in 1:length(mal.sum$mal.d)){
  mal.sum$mal.d[i] = sum(mal$Dead[mal$Conc == mal.sum$mal.c[i]])
}

mal.sum$mort = mal.sum$mal.d / mal.sum$mal.total

mal.mod<-drm(mal.d / mal.total ~ mal.c, weights = mal.total, data = mal.sum, type = 'binomial',  
             fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                        fixed = c(NA, 0, 1, NA)))
summary(mal.mod)
plot(mal.mod, ylim = c(0,1))

muPq_mal_Halstead<-function(In){
  predict(mal.mod, data.frame(mal.c = In))
}  

mal.hal.df = data.frame(mal.c = c(0:40000),
                        Prediction = 0,
                        Lower = 0,
                        Upper = 0)

mal.hal.df[,2:4] <- predict(mal.mod, newdata = mal.hal.df, 
                            interval = 'confidence', level = 0.95)

plot(mal.sum$mal.c, mal.sum$mort, pch = 16, xlab = 'Malathion Concentration (ppb)', ylab = 'prop dead', ylim = c(0,1))

lines(mal.hal.df$mal.c, mal.hal.df$Prediction, col = 2, lty=2)
lines(mal.hal.df$mal.c, mal.hal.df$Lower, col = 2, lty=3)
lines(mal.hal.df$mal.c, mal.hal.df$Upper, col = 2, lty=3)

muPq_mal_Halstead_uncertainty<-function(In){
  m1 = predict(mal.mod, newdata = data.frame(mal.c = In), se.fit = TRUE)[1]
  m2 = predict(mal.mod, newdata = data.frame(mal.c = In), se.fit = TRUE)[2]
  
  mup = rnorm(1,m1,m2)
  if(mup < 0){mup = 0}
  else if(mup > 1){mup = 1}
  
  mup
} 

#Chlorpyrifos   ********************************************************************************************
chlor<-subset(data, Chem == 'Chlor')
chlor.sum = data.frame(chlor.c = unique(chlor$Conc),
                       chlor.total = 5,
                       chlor.d = 0,
                       mort = 0)

for(i in 1:length(chlor.sum$chlor.d)){
  chlor.sum$chlor.d[i] = sum(chlor$Dead[chlor$Conc == chlor.sum$chlor.c[i]])
}

chlor.sum$mort = chlor.sum$chlor.d / chlor.sum$chlor.total

chlor.mod<-drm(chlor.d / chlor.total ~ chlor.c, weights = chlor.total, data = chlor.sum, type = 'binomial',  
               fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                          fixed = c(NA, 0, 1, NA)))
summary(chlor.mod)
plot(chlor.mod, ylim = c(0,1))

muPq_chlor_Halstead<-function(In){
  predict(chlor.mod, data.frame(chlor.c = In))
}  

chlor.hal.df = data.frame(chlor.c = seq(0,100,0.1),
                          Prediction = 0,
                          Lower = 0,
                          Upper = 0)

chlor.hal.df[,2:4] <- predict(chlor.mod, newdata = chlor.hal.df, 
                              interval = 'confidence', level = 0.95)

plot(chlor.sum$chlor.c, chlor.sum$mort, pch = 16, xlab = 'Chlorpyrifos Concentration (ppb)', 
     ylab = 'prop dead', ylim = c(0,1))

lines(chlor.hal.df$chlor.c, chlor.hal.df$Prediction, col = 2, lty=2)
lines(chlor.hal.df$chlor.c, chlor.hal.df$Lower, col = 2, lty=3)
lines(chlor.hal.df$chlor.c, chlor.hal.df$Upper, col = 2, lty=3)

muPq_chlor_Halstead_uncertainty<-function(In){
  m1 = predict(chlor.mod, newdata = data.frame(chlor.c = In), se.fit = TRUE)[1]
  m2 = predict(chlor.mod, newdata = data.frame(chlor.c = In), se.fit = TRUE)[2]
  
  mup = rnorm(1,m1,m2)
  if(mup < 0){mup = 0}
  else if(mup > 1){mup = 1}
  
  mup
} 

#Terbufos   ********************************************************************************************
terb<-subset(data, Chem == 'Terb')
terb.sum = data.frame(terb.c = unique(terb$Conc),
                      terb.total = 5,
                      terb.d = 0,
                      mort = 0)

for(i in 1:length(terb.sum$terb.d)){
  terb.sum$terb.d[i] = sum(terb$Dead[terb$Conc == terb.sum$terb.c[i]])
}

terb.sum$mort = terb.sum$terb.d / terb.sum$terb.total

terb.mod<-drm(terb.d / terb.total ~ terb.c, weights = terb.total, data = terb.sum, type = 'binomial', 
              fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                         fixed = c(NA, 0, 1, NA)))
summary(terb.mod)
plot(terb.mod, ylim = c(0,1))

muPq_terb_Halstead<-function(In){
  predict(terb.mod, data.frame(terb.c = In))
}  

terb.hal.df = data.frame(terb.c = seq(0,200,0.1),
                         Prediction = 0,
                         Lower = 0,
                         Upper = 0)

terb.hal.df[,2:4] <- predict(terb.mod, newdata = terb.hal.df, 
                             interval = 'confidence', level = 0.95)

plot(terb.sum$terb.c, terb.sum$mort, pch = 16, xlab = 'Terbufos Concentration (ppb)', 
     ylab = 'prop dead', ylim = c(0,1))

lines(terb.hal.df$terb.c, terb.hal.df$Prediction, col = 2, lty=2)
lines(terb.hal.df$terb.c, terb.hal.df$Lower, col = 2, lty=3)
lines(terb.hal.df$terb.c, terb.hal.df$Upper, col = 2, lty=3)

muPq_terb_Halstead_uncertainty<-function(In){
  m1 = predict(terb.mod, newdata = data.frame(terb.c = In), se.fit = TRUE)[1]
  m2 = predict(terb.mod, newdata = data.frame(terb.c = In), se.fit = TRUE)[2]
  
  mup = rnorm(1,m1,m2)
  if(mup < 0){mup = 0}
  else if(mup > 1){mup = 1}
  
  mup
} 

#Lambda-cyhalothrin   ********************************************************************************************
lamcy<-subset(data, Chem == 'Lambda')
lamcy.sum = data.frame(lamcy.c = unique(lamcy$Conc),
                       lamcy.total = 5,
                       lamcy.d = 0,
                       mort = 0)

for(i in 1:length(lamcy.sum$lamcy.d)){
  lamcy.sum$lamcy.d[i] = sum(lamcy$Dead[lamcy$Conc == lamcy.sum$lamcy.c[i]])
}

lamcy.sum$mort = lamcy.sum$lamcy.d / lamcy.sum$lamcy.total

lamcy.mod<-drm(lamcy.d / lamcy.total ~ lamcy.c, weights = lamcy.total, data = lamcy.sum, type = 'binomial', 
               fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                          fixed = c(NA, 0, 1, NA)))
summary(lamcy.mod)
plot(lamcy.mod, ylim = c(0,1))

muPq_lamcy_Halstead<-function(In){
  predict(lamcy.mod, data.frame(lamcy.c = In))
}  

lamcy.hal.df = data.frame(lamcy.c = seq(0,5,0.01),
                          Prediction = 0,
                          Lower = 0,
                          Upper = 0)

lamcy.hal.df[,2:4] <- predict(lamcy.mod, newdata = lamcy.hal.df, 
                              interval = 'confidence', level = 0.95)

plot(lamcy.sum$lamcy.c, lamcy.sum$mort, pch = 16, xlab = 'Lambda-Cyhalothrin Concentration (ppb)', 
     ylab = 'prop dead', ylim = c(0,1))

lines(lamcy.hal.df$lamcy.c, lamcy.hal.df$Prediction, col = 2, lty=2)
lines(lamcy.hal.df$lamcy.c, lamcy.hal.df$Lower, col = 2, lty=3)
lines(lamcy.hal.df$lamcy.c, lamcy.hal.df$Upper, col = 2, lty=3)

muPq_lamcy_Halstead_uncertainty<-function(In){
  m1 = predict(lamcy.mod, newdata = data.frame(lamcy.c = In), se.fit = TRUE)[1]
  m2 = predict(lamcy.mod, newdata = data.frame(lamcy.c = In), se.fit = TRUE)[2]
  
  mup = rnorm(1,m1,m2)
  if(mup < 0){mup = 0}
  else if(mup > 1){mup = 1}
  
  mup
} 

#esfenvalerate  ********************************************************************************************
esfen<-subset(data, Chem == 'Esfen')
esfen.sum = data.frame(esfen.c = unique(esfen$Conc),
                       esfen.total = 5,
                       esfen.d = 0,
                       mort = 0)

for(i in 1:length(esfen.sum$esfen.d)){
  esfen.sum$esfen.d[i] = sum(esfen$Dead[esfen$Conc == esfen.sum$esfen.c[i]])
}

esfen.sum$mort = esfen.sum$esfen.d / esfen.sum$esfen.total

esfen.mod<-drm(esfen.d / esfen.total ~ esfen.c, weights = esfen.total, data = esfen.sum, type = 'binomial', 
               fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                          fixed = c(NA, 0, 1, NA)))
summary(esfen.mod)
plot(esfen.mod, ylim = c(0,1))

muPq_esfen_Halstead<-function(In){
  predict(esfen.mod, data.frame(esfen.c = In))
}  

esfen.hal.df = data.frame(esfen.c = seq(0,50,0.1),
                          Prediction = 0,
                          Lower = 0,
                          Upper = 0)

esfen.hal.df[,2:4] <- predict(esfen.mod, newdata = esfen.hal.df, 
                              interval = 'confidence', level = 0.95)

plot(esfen.sum$esfen.c, esfen.sum$mort, pch = 16, xlab = 'Esfenvalerate Concentration (ppb)', 
     ylab = 'prop dead', ylim = c(0,1))

lines(esfen.hal.df$esfen.c, esfen.hal.df$Prediction, col = 2, lty=2)
lines(esfen.hal.df$esfen.c, esfen.hal.df$Lower, col = 2, lty=3)
lines(esfen.hal.df$esfen.c, esfen.hal.df$Upper, col = 2, lty=3)

muPq_esfen_Halstead_uncertainty<-function(In){
  m1 = predict(esfen.mod, newdata = data.frame(esfen.c = In), se.fit = TRUE)[1]
  m2 = predict(esfen.mod, newdata = data.frame(esfen.c = In), se.fit = TRUE)[2]
  
  mup = rnorm(1,m1,m2)
  if(mup < 0){mup = 0}
  else if(mup > 1){mup = 1}
  
  mup
} 

#Permethrin ********************************************************************************************
perm<-subset(data, Chem == 'Perm')
perm.sum = data.frame(perm.c = unique(perm$Conc),
                      perm.total = 5,
                      perm.d = 0,
                      mort = 0)

for(i in 1:length(perm.sum$perm.d)){
  perm.sum$perm.d[i] = sum(perm$Dead[perm$Conc == perm.sum$perm.c[i]])
}

perm.sum$mort = perm.sum$perm.d / perm.sum$perm.total

perm.mod<-drm(perm.d / perm.total ~ perm.c, weights = perm.total, data = perm.sum, type = 'binomial', 
              fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                         fixed = c(NA, 0, 1, NA)))
summary(perm.mod)
plot(perm.mod, ylim = c(0,1))

muPq_perm_Halstead<-function(In){
  predict(perm.mod, data.frame(perm.c = In))
}  

perm.hal.df = data.frame(perm.c = seq(0,10,0.01),
                         Prediction = 0,
                         Lower = 0,
                         Upper = 0)

perm.hal.df[,2:4] <- predict(perm.mod, newdata = perm.hal.df, 
                             interval = 'confidence', level = 0.95)

plot(perm.sum$perm.c, perm.sum$mort, pch = 16, xlab = 'Permethrin Concentration (ppb)', 
     ylab = 'prop dead', ylim = c(0,1))

lines(perm.hal.df$perm.c, perm.hal.df$Prediction, col = 2, lty=2)
lines(perm.hal.df$perm.c, perm.hal.df$Lower, col = 2, lty=3)
lines(perm.hal.df$perm.c, perm.hal.df$Upper, col = 2, lty=3)

muPq_perm_Halstead_uncertainty<-function(In){
  m1 = predict(perm.mod, newdata = data.frame(perm.c = In), se.fit = TRUE)[1]
  m2 = predict(perm.mod, newdata = data.frame(perm.c = In), se.fit = TRUE)[2]
  
  mup = rnorm(1,m1,m2)
  if(mup < 0){mup = 0}
  else if(mup > 1){mup = 1}
  
  mup
}



#*$*$*$*$*$*$*$*$*$*******ALL BELOW NEEDS TO BE UPDATED*******$*$*$*$*$*$*$*$*$*$*$* #####################
#Insecticide toxicity to crustaceans from Halstead et al; 10-day mortality endpoints 
data2<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/Cray.10d.LC50.2012.act2.csv')
#Malathion
mal10<-subset(data2, Chem == 'Mal')
mal10.c<-unique(mal10$Conc)
mal10.d<-as.numeric()
for(i in 1:length(mal10.c)){
  mal10.d[i] = sum(mal10$Dead[mal10$Conc == mal10.c[i]])/5
}
mal.mod10day<-drm(mal10.d ~ mal10.c, type = 'binomial', fct = LL2.2())

muPq_mal_Halstead10day<-function(In){
  1/(1+exp(mal.mod10day$coefficients[1]*(log(In)-mal.mod10day$coefficients[2])))
}  

dose.m10 = seq(min(mal10.c), max(mal10.c), 1)
plot(mal10.c, mal10.d, pch = 16, xlab = 'Concentration', ylab = 'mortality')
lines(dose.m10, muPq_mal_Halstead10day(dose.m10), lty=2, col='red')


#Chlorpyrifos   
chlor10<-subset(data2, Chem == 'Chlor')
chlor10.c<-unique(chlor10$Conc)
chlor10.d<-as.numeric()
for(i in 1:length(chlor10.c)){
  chlor10.d[i] = sum(chlor10$Dead[chlor10$Conc == chlor10.c[i]])/5
}
chlor.mod10day<-drm(chlor10.d ~ chlor10.c, type = 'binomial', fct = LL2.2())

muPq_ch_Halstead10day<-function(In){
  1/(1+exp(chlor.mod10day$coefficients[1]*(log(In)-chlor.mod10day$coefficients[2])))
}

dose.ch10 = seq(min(chlor10.c), max(chlor10.c), 1)
plot(chlor10.c, chlor10.d, pch = 16, xlab = 'Concentration', ylab = 'mortality')
lines(dose.ch10, muPq_ch_Halstead10day(dose.ch10), lty=2, col='red')

#Terbufos   
terb10<-subset(data2, Chem == 'Terb')
terb10.c<-unique(terb10$Conc)
terb10.d<-as.numeric()
for(i in 1:length(terb10.c)){
  terb10.d[i] = sum(terb10$Dead[terb10$Conc == terb10.c[i]])/5
}
terb.mod10day<-drm(terb10.d ~ terb10.c, type = 'binomial', fct = LL2.2())

muPq_terb_Halstead10day<-function(In){
  1/(1+exp(terb.mod10day$coefficients[1]*(log(In)-terb.mod10day$coefficients[2])))
}

dose.tb10 = seq(min(terb10.c), max(terb10.c), 1)
plot(terb10.c, terb10.d, pch = 16, xlab = 'Concentration', ylab = 'mortality')
lines(dose.tb10, muPq_terb_Halstead10day(dose.tb10), lty=2, col='red')


#Lambda-cyhalothrin   
lamcy10<-subset(data2, Chem == 'Lambda')
lamcy10.c<-unique(lamcy10$Conc)
lamcy10.d<-as.numeric()
for(i in 1:length(lamcy10.c)){
  lamcy10.d[i] = sum(lamcy10$Dead[lamcy10$Conc == lamcy10.c[i]])/5
}
lamcy.mod10day<-drm(lamcy10.d ~ lamcy10.c, type = 'binomial', fct = LL2.2())

muPq_lamcy_Halstead10day<-function(In){
  1/(1+exp(lamcy.mod10day$coefficients[1]*(log(In)-lamcy.mod10day$coefficients[2])))
}

dose.lamcy10= seq(min(lamcy10.c), max(lamcy10.c), 0.01)
plot(lamcy10.c, lamcy10.d, pch = 16, xlab = 'Concentration', ylab = 'mortality')
lines(dose.lamcy10, muPq_lamcy_Halstead10day(dose.lamcy10), lty=2, col='red')

#esfenvalerate  
esfen10<-subset(data2, Chem == 'Esfen')
esfen10.c<-unique(esfen10$Conc)
esfen10.d<-as.numeric()
for(i in 1:length(esfen10.c)){
  esfen10.d[i] = sum(esfen10$Dead[esfen10$Conc == esfen10.c[i]])/5
}
esfen.mod10day<-drm(esfen10.d ~ esfen10.c, type = 'binomial', fct = LL2.2())

muPq_esfen_Halstead10day<-function(In){
  1/(1+exp(esfen.mod10day$coefficients[1]*(log(In)-esfen.mod10day$coefficients[2])))
}

dose.esfen10= seq(min(esfen10.c), max(esfen10.c), 0.01)
plot(esfen10.c, esfen10.d, pch = 16, xlab = 'Concentration', ylab = 'mortality')
lines(dose.esfen10, muPq_esfen_Halstead10day(dose.esfen10), lty=2, col='red')

#Permethrin
perm10<-subset(data2, Chem == 'Perm')
perm10.c<-unique(perm10$Conc)
perm10.d<-as.numeric()
for(i in 1:length(perm10.c)){
  perm10.d[i] = sum(perm10$Dead[perm10$Conc == perm10.c[i]])/5
}
perm.mod10day<-drm(perm10.d ~ perm10.c, type = 'binomial', fct = LL2.2())

muPq_perm_Halstead10day<-function(In){
  1/(1+exp(perm.mod10day$coefficients[1]*(log(In)-perm.mod10day$coefficients[2])))
}

dose.perm10= seq(min(perm10.c), max(perm10.c), 0.01)
plot(perm10.c, perm10.d, pch = 16, xlab = 'Concentration', ylab = 'mortality')
lines(dose.perm10, muPq_perm_Halstead10day(dose.perm10), lty=2, col='red')

#Use 4- & 10- day LC50 values from Halstead et al to estimate 1 day LC50 according to Sanchez-Bayo 2009 ########
lcs<-data.frame(chem = c(rep('esfen',2),
                         rep('lamcy',2),
                         rep('perm',2),
                         rep('chlor',2),
                         rep('mal',2),
                         rep('terb',2)),
                days = c(rep(c(4,10), 6)),
                hours = c(rep(c(96,240), 6)),
                lc50 = c(0.22, 0.22,
                         0.21, 0.21,
                         0.58, 0.58,
                         29.3, 6.26, 
                         48936, 32935, 
                         8.89, 8.89),
                lc50_hi = c(0.26, 0.26,
                            0.28, 0.28,
                            1.39, 1.39,
                            32.9, 7.66,
                            159822, 85157,
                            12.8, 12.8),
                lc50_lo = c(0.14, 0.14,
                            0.14, 0.14,
                            0.54, 0.54,
                            19.2, 5.12, 
                            14984, 12738,
                            7.77, 7.77))    

ggplot(data = lcs, aes(x = log(lc50), y = log(days), col = chem, group = chem)) +  
  theme_bw() +
  geom_point() +
  ylim(0, log(10.5)) +
  geom_hline(yintercept = log(1), linetype = 2)

#Get slope of log(time) over log(LC50) for Chlorpyrifos
ch.mb = lm(log(lcs$days[lcs$chem == 'chlor']) ~ log(lcs$lc50[lcs$chem == 'chlor']))

plot(x = log(lcs$lc50[lcs$chem == 'chlor']), y = log(lcs$days[lcs$chem == 'chlor']),
     xlim = c(0,7), ylim = c(0,3), ylab = 'log(days)', xlab = 'log(lc50)',
     pch = 16)
abline(a = ch.mb$coefficients[1], b = ch.mb$coefficients[2], lty = 2, col = 'red')
segments(x0 = (log(lcs$lc50_lo[lcs$chem == 'chlor'])[1]), x1 = (log(lcs$lc50_hi[lcs$chem == 'chlor'])[1]),
         y0 = (log(lcs$days[lcs$chem == 'chlor'])[1]), y1 = (log(lcs$days[lcs$chem == 'chlor'])[1]))
segments(x0 = (log(lcs$lc50_lo[lcs$chem == 'chlor'])[2]), x1 = (log(lcs$lc50_hi[lcs$chem == 'chlor'])[2]),
         y0 = (log(lcs$days[lcs$chem == 'chlor'])[2]), y1 = (log(lcs$days[lcs$chem == 'chlor'])[2]))

ch.1day = exp((0 - ch.mb$coefficients[1]) / ch.mb$coefficients[2])