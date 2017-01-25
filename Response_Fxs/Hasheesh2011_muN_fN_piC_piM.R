#Hasheesh and Mohamed 2011 data and analysis assessing toxicity of Chlorpyrifos and Profenofos ###############
# Direct toxicity to snails affecting snail mortality rate (Table 1); daily mortality rate ########

#Chlorpyrifos
  lc50.ch.n = 1.32
  slp.ch.n = 2.5
  
  mu_Nq_ch_has11<-function(In){
    Ins = In/1000
    1 - 1/(1+exp(slp.ch.n*(log(Ins)-log(lc50.ch.n))))
  } 
  
  mu_Nq_ch_has11(64)
  
    plot(c(0:2820), mu_Nq_ch_has11(c(0:2820)), lwd = 2, type = 'l', xlab = 'chlorpyrifos concentration (ppb)',
         ylab = 'mu_Nq', ylim = c(0,1), main = 'chlorpyrifos toxicity to snails, hasheesh2011')
    
#Profenofos    
  lc50.pr.n = 2.5
  slp.pr.n = 1.6
  mu_Nq_pr_has11<-function(In){
    Ins = In/1000
    1 - 1/(1+exp(slp.pr.n*(log(Ins)-log(lc50.pr.n))))
  } 
  
  mu_Nq_pr_has11(100)
  
    plot(c(0:3720), mu_Nq_pr_has11(c(0:3720)), lwd = 2, type = 'l', xlab = 'profenofos concentration (ppb)',
         ylab = 'mu_Nq', ylim = c(0,1), main = 'profenofos toxicity to snails, hasheesh2011')

#Compare models above to observed data from longitudinal exposure to LC25 concnentrations
fn<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail reproduction/Hasheesh2011_snail_mort_repro_weekly.csv')
 
  plot(fn$time[fn$chem == 'control'], fn$surv[fn$chem == 'control'], pch = 16, cex=1.2,
       xlab = 'time (days)', ylab = 'number alive', ylim = c(0,50),
       main = 'Snail survival over time Hasheesh 2011')
    points(fn$time[fn$chem == 'chlorpyrifos'], fn$surv[fn$chem == 'chlorpyrifos'], col = 'red', pch = 16, cex=1.2)
    points(fn$time[fn$chem == 'profenofos'], fn$surv[fn$chem == 'profenofos'], col = 'orange', pch = 16, cex=1.2)
    legend('topright', legend = c('control', 'chlorpyrifos', 'profenofos'),
           pch=16 , col = c(1,2, 'orange'), cex = 0.8)
    
#Generate model predictions    
  #Background mortality seems to be pretty constant (constant slope in control group), so let's use control group to check that value
    m = (fn$surv[fn$chem == 'control'][11] - fn$surv[fn$chem == 'control'][1]) /
        (fn$time[fn$chem == 'control'][11] - fn$time[fn$chem == 'control'][1])
    
    m.df = data.frame('days' = c(0:70),
                      'control' = 0,
                      'chlor' = 0,
                      'prof' = 0)
    
    m.df[1,c(2:4)] = 50
    
    for(i in 1:70){
      m.df[i+1, 2] = m.df[i,2] + m  #subtract deaths per day (from control slope derived above)
      m.df[i+1, 3] = m.df[i,3] + m - m.df[i,3] * mu_Nq_ch_has11(720)  #for agrochemical toxicity, additional mortality as per capita deaths
      m.df[i+1, 4] = m.df[i,4] + m - m.df[i,4] * mu_Nq_pr_has11(1400) #for agrochemical toxicity, additional mortality as per capita deaths
    }
    
    lines(m.df$days, m.df$control, lty = 2, lwd=2)
    lines(m.df$days, m.df$chlor, lty = 2, col=2, lwd=2)
    lines(m.df$days, m.df$prof, lty = 2, col='orange', lwd=2)
    
    
# Direct toxicity to snails affecting reproduction (Table 2) ###########
  plot(fn$time[fn$chem == 'control'], fn$repro[fn$chem == 'control'], type = 'l', lwd=2,
        xlab = 'time (weeks)', ylab = 'eggs/snail', ylim = c(0,max(fn$repro)+2),
        main = 'Snail reproduction over time Hasheesh 2011')
    lines(fn$time[fn$chem == 'chlorpyrifos'], fn$repro[fn$chem == 'chlorpyrifos'], col = 'red', lwd=2)
    lines(fn$time[fn$chem == 'profenofos'], fn$repro[fn$chem == 'profenofos'], col = 'orange', lwd=2) 
    legend('topleft', legend = c('control', 'chlorpyrifos', 'profenofos'),
           lwd = 2, col = c(1,2, 'orange'), cex = 0.8)
#Derive function as percent reduction in mean eggs/snail over 10 week experiment 
  #this is slightly different from Ibrahim which was relative juveniles produced; i.e. incorporated hatchability
  #chlorpyrifos
    ch.red = sum(fn$repro[fn$chem == 'chlorpyrifos']) / sum(fn$repro[fn$chem == 'control'])
    
    ch.rel = c(1, ch.red)
    ch.conc = c(0, 720)
    
    plot(ch.conc, ch.rel, ylim = c(0,1), xlab = 'chlorpyrifos concentration (ppb)', 
         ylab = 'relative decrase in fecundity', pch = 16)
    
    #nls estimate for a negative exponential would not converge (only two data points to fit to);
    #Estimate of 0.004 as the coefficient comes from fit in excel, appears to fit pretty well
    lines(c(0:750), exp(-0.0037*c(0:750)), lty = 2, col = 2)
       
      f_Nq_chlor_hash11_exp<-function(In){
        exp(-0.0037 * In)
      }  
      
    #Try a linear fit as well
    lines(c(0:750), 1-((ch.red-1)/(-720))*c(0:750), lty = 3, col = 2)
      
      f_Nq_chlor_hash11_lin<-function(In){
        1 - 0.0013*In
      } 
   
  #profenofos
  pr.red = sum(fn$repro[fn$chem == 'profenofos']) / sum(fn$repro[fn$chem == 'control'])
  
  pr.rel = c(1, pr.red)
  pr.conc = c(0, 1400)
  
  plot(pr.conc, pr.rel, ylim = c(0,1), xlab = 'profenofos concentration (ppb)', 
       ylab = 'relative decrase in fecundity', pch = 16)
  
  #nls estimate for a negative exponential would not converge (only two data points to fit to);
  #Estimate of 0.001 as the coefficient comes from fit in excel, appears to fit pretty well
  lines(c(0:1500), exp(-0.0012*c(0:1500)), lty = 2, col = 2)
  
    f_Nq_prof_hash11_exp<-function(In){
      exp(-0.0012 * In)
    }  
  
  #Try a linear fit as well
  lines(c(0:1500), 1-((pr.red-1)/(-1400))*c(0:1500), lty = 3, col = 2)
 
    f_Nq_prof_hash11_lin<-function(In){
      1 - 0.00057*In
    } 
  
    
#Toxicity to miracidia and cercariae from table 5 ###############
#ChlorP miracidia
  lc50m.chlor = 0.78
  slpm.chlor = 1.86
  
    pi_Mq_ch_has11<-function(In){
      Ins = In/1000
      1 - 1/(1+exp(slpm.chlor*(log(Ins)-log(lc50m.chlor))))
    }  
    
    pi_Mq_ch_has11(780)
    
    in.test = seq(0,10000,1)
    
    plot(in.test, pi_Mq_ch_has11(in.test), type = 'l', lwd = 2, xlim = c(0,10000),
        xlab = 'chlorpyrifos concentration (ppb)', ylab = 'pi_M', 
        main = 'Hasheesh 2011 miracidia chlorpyrifos data')

#Profenofos miracidia
  lc50m.prof = 1.5
  slpm.prof = 1.64
  
    pi_Mq_pr_has11<-function(In){
      Ins = In/1000
      1 - 1/(1+exp(slpm.prof*(log(Ins)-log(lc50m.prof))))
    }  
    
    pi_Mq_pr_has11(1500)
    
    plot(in.test, pi_Mq_pr_has11(in.test), type = 'l', lwd = 2, xlim = c(0,10000), ylim = c(0,1),
         xlab = 'profenofos concentration (ppb)', ylab = 'pi_M', 
         main = 'Hasheesh 2011 miracidia profenofos data')
    
#ChlorP cercariae
  lc50c.chlor = 0.96
  slpc.chlor = 2.3
  
    pi_Cq_ch_has11<-function(In){
      Ins = In/1000
     1 - 1/(1+exp(slpc.chlor*(log(Ins)-log(lc50c.chlor))))
    }  
    
    pi_Cq_ch_has11(960)
    
    plot(in.test, pi_Cq_ch_has11(in.test), type = 'l', lwd = 2, xlim = c(0,10000),
         xlab = 'chlorpyrifos concentration (ppb)', ylab = 'pi_C', 
         main = 'Hasheesh 2011 cercariae chlorpyrifos data')
    
#profenofos cercariae
lc50c.prof = 1.85
slpc.prof = 1.84

  pi_Cq_pr_has11<-function(In){
    Ins = In/1000
    1 - 1/(1+exp(slpc.prof*(log(Ins)-log(lc50c.prof))))
  }  
  
  pi_Cq_pr_has11(1850)
  
  plot(in.test, pi_Cq_pr_has11(in.test), type = 'l', lwd = 2, xlim = c(0,10000),
       xlab = 'profenofos concentration (ppb)', ylab = 'pi_C', 
       main = 'Hasheesh 2011 cercariae profenofos data')
  