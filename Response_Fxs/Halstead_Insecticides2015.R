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
  #Malathion
  mal<-subset(data, Chem == 'Mal')
    mal.c<-unique(mal$Conc)
    mal.d<-as.numeric()
    for(i in 1:length(mal.c)){
      mal.d[i] = sum(mal$Dead[mal$Conc == mal.c[i]])/5
    }
      mal.mod<-drm(mal.d ~ mal.c, type = 'binomial', fct = LL2.2())
      
      muPq_mal_Halstead<-function(In){
        1/(1+exp(mal.mod$coefficients[1]*(log(In)-mal.mod$coefficients[2])))
      }  
      
      dose.m = seq(min(mal.c), max(mal.c), 1)
      plot(mal.c, mal.d, pch = 16, xlab = 'Concentration', ylab = 'mortality')
        lines(dose.m, muPq_mal_Halstead(dose.m), lty=2, col='red')
        

  #Chlorpyrifos   
  chlor<-subset(data, Chem == 'Chlor')
    chlor.c<-unique(chlor$Conc)
    chlor.d<-as.numeric()
    for(i in 1:length(chlor.c)){
      chlor.d[i] = sum(chlor$Dead[chlor$Conc == chlor.c[i]])/5
    }
      chlor.mod<-drm(chlor.d ~ chlor.c, type = 'binomial', fct = LL2.2())
  
      muPq_ch_Halstead<-function(In){
        1/(1+exp(chlor.mod$coefficients[1]*(log(In)-chlor.mod$coefficients[2])))
      }
      
      dose.ch = seq(min(chlor.c), max(chlor.c), 1)
      plot(chlor.c, chlor.d, pch = 16, xlab = 'Concentration', ylab = 'mortality')
        lines(dose.ch, muPq_ch_Halstead(dose.ch), lty=2, col='red')
      
  #Terbufos   
  terb<-subset(data, Chem == 'Terb')
    terb.c<-unique(terb$Conc)
    terb.d<-as.numeric()
    for(i in 1:length(terb.c)){
      terb.d[i] = sum(terb$Dead[terb$Conc == terb.c[i]])/5
    }
      terb.mod<-drm(terb.d ~ terb.c, type = 'binomial', fct = LL2.2())
  
      muPq_terb_Halstead<-function(In){
        1/(1+exp(terb.mod$coefficients[1]*(log(In)-terb.mod$coefficients[2])))
      }
      
      dose.tb = seq(min(terb.c), max(terb.c), 1)
      plot(terb.c, terb.d, pch = 16, xlab = 'Concentration', ylab = 'mortality')
        lines(dose.tb, muPq_terb_Halstead(dose.tb), lty=2, col='red')
      
  
  #Lambda-cyhalothrin   
  lamcy<-subset(data, Chem == 'Lambda')
    lamcy.c<-unique(lamcy$Conc)
    lamcy.d<-as.numeric()
    for(i in 1:length(lamcy.c)){
      lamcy.d[i] = sum(lamcy$Dead[lamcy$Conc == lamcy.c[i]])/5
    }
      lamcy.mod<-drm(lamcy.d ~ lamcy.c, type = 'binomial', fct = LL2.2())
  
      muPq_lamcy_Halstead<-function(In){
        1/(1+exp(lamcy.mod$coefficients[1]*(log(In)-lamcy.mod$coefficients[2])))
      }
      
      dose.lamcy= seq(min(lamcy.c), max(lamcy.c), 0.01)
      plot(lamcy.c, lamcy.d, pch = 16, xlab = 'Concentration', ylab = 'mortality')
        lines(dose.lamcy, muPq_lamcy_Halstead(dose.lamcy), lty=2, col='red')
      
  #esfenvalerate  
  esfen<-subset(data, Chem == 'Esfen')
    esfen.c<-unique(esfen$Conc)
    esfen.d<-as.numeric()
    for(i in 1:length(esfen.c)){
      esfen.d[i] = sum(esfen$Dead[esfen$Conc == esfen.c[i]])/5
    }
      esfen.mod<-drm(esfen.d ~ esfen.c, type = 'binomial', fct = LL2.2())
  
      muPq_esfen_Halstead<-function(In){
        1/(1+exp(esfen.mod$coefficients[1]*(log(In)-esfen.mod$coefficients[2])))
      }
      
      dose.esfen= seq(min(esfen.c), max(esfen.c), 0.01)
      plot(esfen.c, esfen.d, pch = 16, xlab = 'Concentration', ylab = 'mortality')
        lines(dose.esfen, muPq_esfen_Halstead(dose.esfen), lty=2, col='red')
      
  #Permethrin
  perm<-subset(data, Chem == 'Perm')
    perm.c<-unique(perm$Conc)
    perm.d<-as.numeric()
    for(i in 1:length(perm.c)){
      perm.d[i] = sum(perm$Dead[perm$Conc == perm.c[i]])/5
    }
      perm.mod<-drm(perm.d ~ perm.c, type = 'binomial', fct = LL2.2())
      
      muPq_perm_Halstead<-function(In){
        1/(1+exp(perm.mod$coefficients[1]*(log(In)-perm.mod$coefficients[2])))
      }
      
      dose.perm= seq(min(perm.c), max(perm.c), 0.01)
      plot(perm.c, perm.d, pch = 16, xlab = 'Concentration', ylab = 'mortality')
        lines(dose.perm, muPq_perm_Halstead(dose.perm), lty=2, col='red')
      
#Insecticide toxicity to crustaceans from Halstead et al; 10-day mortality endpoints ################
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