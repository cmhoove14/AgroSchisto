#Data extraction and model fitting to Bakry 2012 data
require(drc)

#Snail (B. alexandrina) toxicity ##########
#Atrazine
lc50.n.atr<-1.25
slp.n.atr<-2.48

mu_Nq_atr_bak12<-function(In){
  Ins = In/1000
  1 - 1/(1+exp(slp.n.atr*(log(Ins)-log(lc50.n.atr))))
} 

mu_Nq_atr_bak12(lc50.n.atr*1000)

plot(c(0:10000), mu_Nq_atr_bak12(c(0:10000)), lwd = 2, type = 'l', xlab = 'atrazine concentration (ppb)',
     ylab = 'mu_Nq', ylim = c(0,1), main = 'atrazine toxicity to snails, bakry2012')

#Deltamethrin     
lc50.n.gly<-3.15
slp.n.gly<-2.16

mu_Nq_gly_bak12<-function(In){
  Ins = In/1000
  1 - 1/(1+exp(slp.n.gly*(log(Ins)-log(lc50.n.gly))))
} 

mu_Nq_gly_bak12(lc50.n.gly*1000)

plot(c(0:10000), mu_Nq_gly_bak12(c(0:10000)), lwd = 2, type = 'l', xlab = 'glyphosate concentration (ppb)',
     ylab = 'mu_Nq', ylim = c(0,1), main = 'glyphosate toxicity to snails, bakry2012')

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
  m.df[i+1, 3] = m.df[i,3] - m.df[i,3] * mu_Nq_atr_bak12(330)  #for agrochemical toxicity, additional mortality as per capita deaths
  m.df[i+1, 4] = m.df[i,4] - m.df[i,4] * mu_Nq_gly_bak12(840) #for agrochemical toxicity, additional mortality as per capita deaths
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