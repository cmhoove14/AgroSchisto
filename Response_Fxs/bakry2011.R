#Data extraction and model fitting to Bakry 2011 data
require(drc)

#Snail (H. duryi) toxicity ##########
#Malathion
lc50.n.mal<-1.76
slp.n.mal<-2.68

mu_Nq_mal_bak11<-function(In){
  Ins = In/1000
  1 - 1/(1+exp(slp.n.mal*(log(Ins)-log(lc50.n.mal))))
} 

mu_Nq_mal_bak11(lc50.n.mal*1000)

plot(c(0:10000), mu_Nq_mal_bak11(c(0:10000)), lwd = 2, type = 'l', xlab = 'malathion concentration (ppb)',
     ylab = 'mu_Nq', ylim = c(0,1), main = 'malathion toxicity to snails, bakry2011')

#Deltamethrin     
lc50.n.del<-4.82
slp.n.del<-2.74

mu_Nq_del_bak11<-function(In){
  Ins = In/1000
  1 - 1/(1+exp(slp.n.del*(log(Ins)-log(lc50.n.del))))
} 

mu_Nq_del_bak11(lc50.n.del*1000)

plot(c(0:10000), mu_Nq_del_bak11(c(0:10000)), lwd = 2, type = 'l', xlab = 'deltamethrin concentration (ppb)',
     ylab = 'mu_Nq', ylim = c(0,1), main = 'deltamethrin toxicity to snails, tantawy2002')
 
#The paper also provides info on longitudinal survival of snail cohorts exposed to LC10 of each insecticide ################
  #So let's compare that data with the expected long. survival from the model at the same concentration
  bakry11<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/bakry2011.csv')

#Plot data from paper  
  plot(bakry11$time[bakry11$chem == 'control'], bakry11$surv[bakry11$chem == 'control'], ylim = c(0,50), pch = 16, cex = 1.2,
       xlab = 'time (days)', ylab = 'n alive')
    points(bakry11$time[bakry11$chem == 'malathion'], bakry11$surv[bakry11$chem == 'malathion'], col = 2, pch = 16, cex = 1.2)
    points(bakry11$time[bakry11$chem == 'deltamethrin'], bakry11$surv[bakry11$chem == 'deltamethrin'], col = 6, pch = 16, cex = 1.2)
    title('long. snail survival exposed to LC10s')
    legend('bottomleft', legend = c('control', 'deltamethrin', 'malathion'), col = c(1,2,6), pch = 16, cex=0.7, 
           title = 'observed')
    
#Generate model predictions    
  #Background mortality seems to be pretty constant (constant slope in control group), so let's use control group to check that value
    m = (bakry11$surv[bakry11$chem == 'control'][7] - bakry11$surv[bakry11$chem == 'control'][1]) /
        (bakry11$time[bakry11$chem == 'control'][7] - bakry11$time[bakry11$chem == 'control'][1])

    m.df = data.frame('days' = c(0:28),
                      'control' = 0,
                      'mal' = 0,
                      'delt' = 0)
    
    m.df[1,c(2:4)] = 50
    
    for(i in 1:28){
      m.df[i+1, 2] = m.df[i,2] + m  #subtract deaths per day (from control slope derived above)
      m.df[i+1, 3] = m.df[i,3] + m - m.df[i,3] * mu_Nq_mal_bak11(480)  #for agrochemical toxicity, additional mortality as per capita deaths
      m.df[i+1, 4] = m.df[i,4] + m - m.df[i,4] * mu_Nq_del_bak11(1210) #for agrochemical toxicity, additional mortality as per capita deaths
    }
    
    lines(m.df$days, m.df$control, lty = 2, lwd=2)
    lines(m.df$days, m.df$mal, lty = 2, col=2, lwd=2)
    lines(m.df$days, m.df$delt, lty = 2, col=6, lwd=2)
      legend('left', legend = c('control', 'deltamethrin', 'malathion'), col = c(1,2,6),  lty=2, cex=0.7, 
             title = 'modeled')
      
#Now lets look at reproduction over time ###########
  plot(bakry11$time[bakry11$chem == 'control'], bakry11$eggs_snail[bakry11$chem == 'control'], type='l', lwd=2, ylim = c(0,12),
       xlab = 'time (days)', ylab = 'mean eggs/snail')    
    lines(bakry11$time[bakry11$chem == 'malathion'], bakry11$eggs_snail[bakry11$chem == 'malathion'], col=2, lwd=2)
    lines(bakry11$time[bakry11$chem == 'deltamethrin'], bakry11$eggs_snail[bakry11$chem == 'deltamethrin'], col=6, lwd=2)
    legend('topright', legend = c('control', 'deltamethrin', 'malathion'), col = c(1,2,6), lwd=2, cex=0.7)
    
#Relative decreases in eggs/snail for 4 week study period
  fn.mal = 3.9/54.4
  fn.delt = 14.33/54.4
  
  plot(1,1, ylab = 'relative eggs/snail', xlab = 'insecticide concentration', ylim = c(0,1), xlim = c(0,1300), pch = 16)
    points(480, fn.mal, col=2, pch=16)
    points(1210, fn.delt, col=6, pch=16)
    
    legend('topright', legend = c('control', 'malathion', 'deltamethrin'), col = c(1,2,6), pch=16, cex=0.7)
    
  #Fit function to malathion data points  
    #nls estimate for a negative exponential would not converge (only two data points to fit to);
    #Estimate of 0.0055 as the coefficient comes from fit in excel, appears to fit pretty well
    lines(c(0:1500), exp(-0.0055*c(0:1500)), lty = 2, col = 2)
    
    f_Nq_mal_bak11_exp<-function(In){
      exp(-0.0055 * In)
    }  
    
    #Try a linear fit as well
    lines(c(0:1500), 1-((fn.mal-1)/(-480))*c(0:1500), lty = 3, col = 2)
    
    f_Nq_mal_bak11_lin<-function(In){
      1 - 0.0019*In
    } 
    
  #Fit function to deltamethrin data points  
    #nls estimate for a negative exponential would not converge (only two data points to fit to);
    #Estimate of 0.0011 as the coefficient comes from fit in excel, appears to fit pretty well
    lines(c(0:1500), exp(-0.0011*c(0:1500)), lty = 2, col = 6)
    
    f_Nq_delt_bak11_exp<-function(In){
      exp(-0.0011 * In)
    }  
    
    #Try a linear fit as well
    lines(c(0:1500), 1-((fn.delt-1)/(-1210))*c(0:1500), lty = 3, col = 6)
    
    f_Nq_mal_bak11_lin<-function(In){
      1 - 0.000609*In
    } 