#Data extraction and model fitting to Bakry 2012 data
require(drc)

#Snail (B. alexandrina) toxicity; create data frame, clean, etc. ##########
muN.bak = data.frame(atr = c(.330, 1.250, 4.750),
                     atr.se = 0,
                     gly = c(.840, 3.150, 12.600),
                     gly.se = 0,
                     mort = c(0.1, 0.5, 0.9))

  muN.bak$surv = 1 - muN.bak$mort

#Standard errors imputed assuming proportional error to concentration
  
  se.atr = (log10(1.88) - log10(1.25)) / 1.96 #st. err of lc50 in ppm
  se.gly = (log10(4.82) - log10(3.15)) / 1.96 #error bars are assymetrical in this data even after log transformation
                                              #so not entirely sure how to handle that...

  muN.bak$atr.se = (muN.bak$atr / 1.25) * se.atr #st. err proportional to concentration
  muN.bak$gly.se = (muN.bak$gly / 3.15) * se.gly #st. err proportional to concentration

#visualize ##############
  plot(muN.bak$atr, muN.bak$mort, ylim = c(0,1), xlim = c(0,13),
       pch = 16, col = 'gold', xlab = 'Herbicide (ppm)', ylab = 'mortality')
    points(muN.bak$gly, muN.bak$mort, pch = 16, col = 3)
    for(i in 1:length(muN.bak$atr)){
      segments(y0 = muN.bak$mort[i], x0 = muN.bak$atr[i] + muN.bak$atr.se[i],
               y1 = muN.bak$mort[i], x1 = muN.bak$atr[i] - muN.bak$atr.se[i], col='gold')
      segments(y0 = muN.bak$mort[i], x0 = muN.bak$gly[i] + muN.bak$gly.se[i],
               y1 = muN.bak$mort[i], x1 = muN.bak$gly[i] - muN.bak$gly.se[i], col=3)
    }
    
#fit functions using drm assuming cohort size of 50 snails for each LC outcome ################
  #50 snails used in longitudinal LC10 exposure studies, so seems like a reasonable assusmption  
  bak12.atr.drm = drm(mort ~ atr, weights = rep(50,3), data = muN.bak, type = 'binomial',
                      fct = LL.2())

  bak12.gly.drm = drm(mort ~ gly, weights = rep(50,3), data = muN.bak, type = 'binomial',
                      fct = LL.2())

#visualize fits ##########
  bak.test.df = data.frame(atr = seq(0,13,0.1),
                           gly = seq(0,13,0.1),
                           atr.muN = 0,
                           atr.se = 0,
                           gly.muN = 0,
                           gly.se = 0)
  
#plot functions to visualize fit to data
  bak.test.df[,3:4] = predict(bak12.atr.drm, newdata = bak.test.df, se.fit = TRUE)
  bak.test.df[,5:6] = predict(bak12.gly.drm, newdata = bak.test.df, se.fit = TRUE)
  
    lines(bak.test.df$atr, bak.test.df$atr.muN, col='gold', lty=2)
      lines(bak.test.df$atr, bak.test.df$atr.muN + 1.96*bak.test.df$atr.se, col='gold', lty=3)
      lines(bak.test.df$atr, bak.test.df$atr.muN - 1.96*bak.test.df$atr.se, col='gold', lty=3)
      
    lines(bak.test.df$gly, bak.test.df$gly.muN, col=3, lty=2)
      lines(bak.test.df$gly, bak.test.df$gly.muN + 1.96*bak.test.df$gly.se, col=3, lty=3)
      lines(bak.test.df$gly, bak.test.df$gly.muN - 1.96*bak.test.df$gly.se, col=3, lty=3)
      
    legend('bottomright', legend = c('atrazine', 'glyphosate', '95% CI'), col = c('gold', 3, 1), 
           cex = 0.7, lty=c(2,2,3))  
#derive functions for each snail mortality response ####################    
  muNq_atr_bak12<-function(He){
    He.use = He/1000
    rdrm(1, LL.2(), coef(bak12.atr.drm), He.use, yerror = 'rbinom', ypar = 100)$y / 100 #estimate deaths / live 
  }
    #for(i in seq(0,13000,10)){
    #  points(i/1000, muNq_atr_bak12(i), pch = 17, cex=0.5, col=2)
    #}

  muNq_gly_bak12<-function(He){
    He.use = He/1000
    rdrm(1, LL.2(), coef(bak12.gly.drm), He.use, yerror = 'rbinom', ypar = 100)$y / 100 #estimate deaths / live 
  }
    #for(i in seq(0,13000,10)){
    #  points(i/1000, muNq_gly_bak12(i), pch = 17, cex=0.5, col=4)
    #}
  
  keep.bak12 = c('muN.bak', 'muNq_atr_bak12', 'muNq_gly_bak12', 'bak12.atr.drm', 'bak12.gly.drm')

#The paper also provides info on longitudinal survival of snail cohorts exposed to LC10 of each herbicide ################
#So let's compare that data with the expected long. survival from the model at the same concentration
bakry12<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/bakry2012.csv')

#Plot data from paper  
plot(bakry12$time[bakry12$chem == 'control'], bakry12$surv[bakry12$chem == 'control'], ylim = c(0,50), pch = 16, cex = 1.2,
     xlab = 'time (days)', ylab = 'n alive')
  points(bakry12$time[bakry12$chem == 'atrazine'], bakry12$surv[bakry12$chem == 'atrazine'], col = 'goldenrod', pch = 16, cex = 1.2)
  points(bakry12$time[bakry12$chem == 'glyphosate'], bakry12$surv[bakry12$chem == 'glyphosate'], col = 'forestgreen', pch = 16, cex = 1.2)
  title('long. snail survival exposed to LC10s')
  legend('bottomleft', legend = c('control', 'atrazine', 'glyphosate'), col = c(1,'goldenrod','forestgreen'), pch = 16, cex=0.7, 
         title = 'observed')

#Generate model predictions    
#Background mortality seems to be pretty constant (constant slope in control group), so let's use control group to check that value
m = (bakry12$surv[bakry12$chem == 'control'][7] - bakry12$surv[bakry12$chem == 'control'][1]) /
    (bakry12$time[bakry12$chem == 'control'][7] - bakry12$time[bakry12$chem == 'control'][1])

m.df = data.frame('days' = c(0:42),
                  'control' = 0,
                  'atr' = 0,
                  'gly' = 0)

m.df[1,c(2:4)] = 50

for(i in 1:42){
  m.df[i+1, 2] = m.df[i,2] + m  #subtract deaths per day (from control slope derived above)
  m.df[i+1, 3] = m.df[i,3] - m.df[i,3] * muNq_atr_bak12(330)  #for agrochemical toxicity, additional mortality as per capita deaths
  m.df[i+1, 4] = m.df[i,4] - m.df[i,4] * muNq_gly_bak12(840) #for agrochemical toxicity, additional mortality as per capita deaths
}

lines(m.df$days, m.df$control, lty = 2, lwd=2)
lines(m.df$days, m.df$atr, lty = 2, col='goldenrod', lwd=2)
lines(m.df$days, m.df$gly, lty = 2, col='forestgreen', lwd=2)
legend('topright', legend = c('control', 'atrazine', 'glyphosate'), col = c(1,'goldenrod','forestgreen'),  lty=2, cex=0.7, 
       title = 'modeled')

#Now lets look at reproduction over time ###########
  plot(bakry12$time[bakry12$chem == 'control'], bakry12$hatch[bakry12$chem == 'control'], type='l', lwd=2, ylim = c(0,450),
       xlab = 'time (days)', ylab = 'total hatches')    
    lines(bakry12$time[bakry12$chem == 'atrazine'], bakry12$hatch[bakry12$chem == 'atrazine'], col='goldenrod', lwd=2)
    lines(bakry12$time[bakry12$chem == 'glyphosate'], bakry12$hatch[bakry12$chem == 'glyphosate'], col='forestgreen', lwd=2)
    legend('topright', legend = c('control', 'atrazine', 'glyphosate'), col = c(1,'goldenrod','forestgreen'), lwd=2, cex=0.7)

#Relative decreases in eggs/snail for 4 week study period
  fn.atr = sum(bakry12$hatch[bakry12$chem == 'atrazine']) / sum(bakry12$hatch[bakry12$chem == 'control'])
  fn.gly = sum(bakry12$hatch[bakry12$chem == 'glyphosate']) / sum(bakry12$hatch[bakry12$chem == 'control'])

plot(1,1, ylab = 'relative hatchlings', xlab = 'herbicide concentration', ylim = c(0,1), xlim = c(0,1000), pch = 16)
  points(330, fn.atr, col='goldenrod', pch=16)
  points(840, fn.gly, col='forestgreen', pch=16)

legend('topright', legend = c('control', 'atrazine', 'glyphosate'), col = c(1,'goldenrod','forestgreen'), pch=16, cex=0.7)

#Fit function to atrazine data points  
#nls estimate for a negative exponential would not converge (only two data points to fit to);
#Estimate of 0.0078 as the coefficient comes from fit in excel, appears to fit pretty well
  lines(c(0:1000), exp(-0.0078*c(0:1000)), lty = 2, col = 'goldenrod')
  
  f_Nq_atr_bak12_exp<-function(In){
    exp(-0.0078 * In)
  }  
  
  #Try a linear fit as well
  lines(c(0:1000), 1-((fn.atr-1)/(-330))*c(0:1000), lty = 3, col = 'goldenrod')
  
  f_Nq_atr_bak12_lin<-function(In){
    1 - 0.0028*In
  } 
  
  #Fit function to deltamethrin data points  
  #nls estimate for a negative exponential would not converge (only two data points to fit to);
  #Estimate of 0.0028 as the coefficient comes from fit in excel, appears to fit pretty well
  lines(c(0:1000), exp(-0.0028*c(0:1000)), lty = 2, col = 'forestgreen')
  
  f_Nq_gly_bak12_exp<-function(In){
    exp(-0.0028 * In)
  }  
  
  #Try a linear fit as well
  lines(c(0:1000), 1-((fn.gly-1)/(-840))*c(0:1000), lty = 3, col = 'forestgreen')
  
  f_Nq_gly_bak12_lin<-function(In){
    1 - 0.0011*In
  } 