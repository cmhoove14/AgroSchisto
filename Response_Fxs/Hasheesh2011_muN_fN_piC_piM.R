#Hasheesh and Mohamed 2011 data and analysis assessing toxicity of Chlorpyrifos and Profenofos ###############
# Direct toxicity to snails (Bu. truncatus); daily mortality rate ########

#Chlorpyrifos #########
mun.ch = data.frame(conc = c(.72, 1.32, 2.82),
                     mort = c(.25, .50 , .90),
                     surv = 0)
mun.ch$surv = 1 - mun.ch$mort
 
  se = (log10(1.98) - log10(0.88)) / 1.96 #st. err of lc50 in ppm

  mun.ch$se = (mun.ch$conc / 1.32) * se #st. err proportional to concentration
 
plot(mun.ch$conc, mun.ch$mort, pch = 16, cex = 1.2, ylim = c(0,1), xlim = c(0,3.5),
     xlab = 'Chlorpyrifos (ppm)', ylab = 'prop dead')
  for(i in 1:length(unique(mun.ch$conc))){
    segments(x0 = mun.ch$conc[i] + mun.ch$se[i], y0 = mun.ch$mort[i],
             x1 = mun.ch$conc[i] - mun.ch$se[i], y1 = mun.ch$mort[i])
  }
  
  has11ch.mod = drm(mort ~ conc, data = mun.ch, weights = se^-1, type = 'binomial',
                  fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                             fixed = c(NA, 0, 1, NA)))
  
  muNq_ch_Hash11<-function(In){
    Ins = In/1000
    predict(has11ch.mod, data.frame(conc = Ins))
  }  

  lines(c(0:4000)/1000, sapply(c(0:4000), muNq_ch_Hash11, simplify = T),
        lty = 2, col = 2)    
  
  muNq_ch_Hash11_uncertainty<-function(In){
    Ins = In/1000
    rdrm(1, LL.2(), coef(has11ch.mod), Ins, yerror = 'rbinom', ypar = 100)$y / 100 #estimate deaths / live 
  }
    points(seq(0, 4000, 10)/1000, sapply(seq(0, 4000, 10), muNq_ch_Hash11_uncertainty, simplify = T),
           pch = 5, cex = 0.6, col = 4)
    title(main = expression(paste('Hasheesh 2011 - Chlorpyrifos direct toxicity to ',
                                  italic('Bu. truncatus', sep = ''))))
      legend('bottomright', pch = c(16, 5), legend = c('Obs. points', 'Est. points'), col = c(1,4),
             cex = 0.7)
  
#Profenofos #########
mun.prof = data.frame(conc = c(1.4, 2.5, 3.72),
                    mort = c(.25, .50 , .90),
                    surv = 0)
  mun.prof$surv = 1 - mun.prof$mort
    
  se.pr = (log10(3.33) - log10(1.88)) / 1.96 #st. err of lc50 in ppm
    
  mun.prof$se = (mun.prof$conc / 2.5) * se.pr #st. err proportional to concentration
    
  plot(mun.prof$conc, mun.prof$mort, pch = 16, cex = 1.2, ylim = c(0,1), xlim = c(0,5),
       xlab = 'Profenofos (ppm)', ylab = 'prop dead')
    for(i in 1:length(unique(mun.prof$conc))){
      segments(x0 = mun.prof$conc[i] + mun.prof$se[i], y0 = mun.prof$mort[i],
               x1 = mun.prof$conc[i] - mun.prof$se[i], y1 = mun.prof$mort[i])
    }
    
    has11pr.mod = drm(mort ~ conc, data = mun.prof, weights = se^-1, type = 'binomial',
                    fct = LL.4(names = c("Slope","Lower Limit","Upper Limit", "ED50"),
                               fixed = c(NA, 0, 1, NA)))
    
    muNq_prof_Hash11<-function(In){
      Ins = In/1000
      predict(has11pr.mod, data.frame(conc = Ins))
    }  
    
    lines(c(0:5000)/1000, sapply(c(0:5000), muNq_prof_Hash11, simplify = T),
          lty = 2, col = 2)   
    
    muNq_prof_Hash11_uncertainty<-function(In){
      Ins = In/1000
      rdrm(1, LL.2(), coef(has11pr.mod), Ins, yerror = 'rbinom', ypar = 100)$y / 100 #estimate deaths / live 
    }
    
    points(seq(0, 5000, 10)/1000, sapply(seq(0, 5000, 10), muNq_prof_Hash11_uncertainty, simplify = T),
           pch = 5, cex = 0.5, col = 4)
      title(main = expression(paste('Hasheesh 2011 - Profenofos direct toxicity to ',
                                    italic('Bu. truncatus', sep = ''))))
      legend('bottomright', pch = c(16, 5), legend = c('Obs. points', 'Est. points'), col = c(1,4),
             cex = 0.7)
    
#Compare models above to observed data from longitudinal exposure to LC25 concnentrations #####
fn<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Snail reproduction/Hasheesh2011_snail_mort_repro_weekly.csv')
 
  plot(fn$time[fn$chem == 'control'], fn$surv[fn$chem == 'control'], pch = 16, cex=1.2,
       xlab = 'time (days)', ylab = 'number alive', ylim = c(0,50),
       main = 'Snail survival over time Hasheesh 2011')
    points(fn$time[fn$chem == 'chlorpyrifos'], fn$surv[fn$chem == 'chlorpyrifos'], col = 'red', pch = 16, cex=1.2)
    points(fn$time[fn$chem == 'profenofos'], fn$surv[fn$chem == 'profenofos'], col = 'orange', pch = 16, cex=1.2)
    legend('topright', legend = c('control', 'chlorpyrifos', 'profenofos'),
           pch=16 , col = c(1,2, 'orange'), cex = 0.8)
    
#Generate model predictions    
  #Background mortality seems to be pretty constant (constant slope in control group), 
    #so let's use control group to check that value
    m = (fn$surv[fn$chem == 'control'][11] - fn$surv[fn$chem == 'control'][1]) /
        (fn$time[fn$chem == 'control'][11] - fn$time[fn$chem == 'control'][1])
    
    m.df = data.frame('days' = c(0:70),
                      'control' = 0,
                      'chlor' = 0,
                      'prof' = 0)
    
    m.df[1,c(2:4)] = 50
    
    for(i in 1:70){
      m.df[i+1, 2] = m.df[i,2] + m  #subtract deaths per day (from control slope derived above)
      m.df[i+1, 3] = m.df[i,3] - m.df[i,3] * muNq_ch_Hash11(720)  #for agrochemical toxicity, additional mortality as per capita deaths
      m.df[i+1, 4] = m.df[i,4] - m.df[i,4] * muNq_prof_Hash11(1400) #for agrochemical toxicity, additional mortality as per capita deaths
    }
    
    lines(m.df$days, m.df$control, lty = 2, lwd=2)
    lines(m.df$days, m.df$chlor, lty = 2, col=2, lwd=2)
    lines(m.df$days, m.df$prof, lty = 2, col='orange', lwd=2)
    
#investigate other effects; but most have only have two data points, 
  #and unclear what D-R model was used to generate, so no D-R function possible #######    
    
# Direct toxicity to snails affecting reproduction (Table 2)
  plot(fn$time[fn$chem == 'control'], fn$repro[fn$chem == 'control'], type = 'l', lwd=2,
        xlab = 'time (weeks)', ylab = 'eggs/snail', ylim = c(0,max(fn$repro)+2),
        main = 'Snail reproduction over time Hasheesh 2011')
    lines(fn$time[fn$chem == 'chlorpyrifos'], fn$repro[fn$chem == 'chlorpyrifos'], col = 'red', lwd=2)
    lines(fn$time[fn$chem == 'profenofos'], fn$repro[fn$chem == 'profenofos'], col = 'orange', lwd=2) 
    legend('topleft', legend = c('control', 'chlorpyrifos', 'profenofos'),
           lwd = 2, col = c(1,2, 'orange'), cex = 0.8)
#Derive function as percent reduction in mean eggs/snail over 10 week experiment 
#this is slightly different from Ibrahim which was relative juveniles produced; 
    #i.e. incorporated hatchability
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
  