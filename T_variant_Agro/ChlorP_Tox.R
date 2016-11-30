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

#Toxicity to Biomphalaria snails from Ibrahim 1992 ###################
  snail.repro = data.frame(dose = c(0,125,250,500),
                           f_red = c(1, 3005/4225, 1749/4225, 1313/4225))
  
  
  plot(snail.repro$dose, snail.repro$f_red, ylim = c(0,1), pch = 16, 
       ylab = 'reduction in snail recruitment rate', xlab = 'ChlorP ppb', main='Ibrahim92 snail recruitment')
  
  mod1<-nls(f_red ~ exp(-b*(dose+1e-6)), data=snail.repro, start = list(b=0.01))
    summary(mod1)
    
  dose = c(0:500)
  
  f_Nq_chlor_ibrahim92 = function(In){
    exp(-0.0028479*(In))
  }
  
  lines(dose, f_Nq_chlor_ibrahim92(dose), lty=2, col='red')
  text(75,0.1, labels = expression(paste('f'[N],'(q) = ', 'e'^'-bq', '     b=0.00285', sep='')))
  
  snail.mort = data.frame(dose = c(0,125,250,500),
                          mort = c(16.7, 23.3, 26.7, 100))
  
  plot(snail.mort$dose, (snail.mort$mort - snail.mort$mort[1])/100, ylim = c(0,1), pch = 16, 
       ylab = 'increased snail mortality', xlab = 'ChlorP ppb', main = 'Ibrahim 92 snail mortality')
  
  ibr_muNq<-drm((snail.mort$mort - snail.mort$mort[1])/100 ~ snail.mort$dose,
                  data = snail.mort, type = 'binomial', fct = LL.2())
  
  muNq_chlor_ibrahim92<-function(In){
    1/(1+exp(ibr_muNq$coefficients[1]*(log(In)-log(ibr_muNq$coefficients[2]))))
  }  
  
  lines(dose, muNq_chlor_ibrahim92(dose), lty=2, col='red')
  text(150,0.5, labels = expression(paste(mu[N],'(q) = ', 'f(SF,LC'[50],',q)', '     SF=-3.87, LC'[50],'=358.5', sep='')), cex=0.7)
  
#Toxicity to prawns from Halstead 2015 ############################
  pred.mort<-data.frame('dose'=c(0,0.64,3.2,6.4,32,64),
                        'death'=c(0,0,0,0,0.8,1))
  
  plot(pred.mort$dose, pred.mort$death, pch = 16, 
       ylab = 'predator mortality', xlab = 'chlorP ppb', main = 'Halstead 15 pred mortality')
  
  hal.chlor<-drm(pred.mort$death ~ pred.mort$dose, data = pred.mort, type = 'binomial', fct = LL.2())
  
  muPq_ch_Halstead4day<-function(In){
    1/(1+exp(hal.chlor$coefficients[1]*(log(In)-log(hal.chlor$coefficients[2]))))
  }
  
  lines(dose, muPq_ch_Halstead4day(dose), lty=2, col='red')
  text(50,0.5, labels = expression(paste(mu[N],'(q) = ', 'f(SF,LC'[50],',q)', '     SF=-11.26, LC'[50],'=28.29', sep='')), cex=0.7)
  
  
  muPq_ch_Halstead<-function(In){
    (1/(1+exp(hal.chlor$coefficients[1]*(log(In)-log(hal.chlor$coefficients[2])))))/4
  }
  
#Toxicity to prawns from Satapornvanit 2009 #############
  dose = seq(0,100,0.1)
  sat.mort<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_mort.csv")
  
  ch.mod.mupq<-drm(sat.mort$mort[sat.mort$chem == 'chlorpyrifos']/100 ~ sat.mort$conc[sat.mort$chem == 'chlorpyrifos'],
                   data = sat.mort, type = 'binomial', fct = LL.2())
  
  muPq_ch_sat09<-function(In){
    1/(1+exp(ch.mod.mupq$coefficients[1]*(log(In)-log(ch.mod.mupq$coefficients[2]))))
  }  
  
  plot(sat.mort$conc[sat.mort$chem == 'chlorpyrifos'], sat.mort$mort[sat.mort$chem == 'chlorpyrifos']/100, 
       pch = 16, ylim = c(0,1),
       ylab = 'predator mortality', xlab = 'chlorP ppb', main = 'Satapornvanit 09 pred mortality')
    
  lines(dose/100, muPq_ch_sat09(dose/100), lty=2, col='red')
  
sap.fr<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/AgroData/Data/Predator Mortality/satapornvanit2009_m.rosenbergii_feed_rate.csv")
  fr.ch = subset(sap.fr, chem == 'chlorpyrifos' & feed_rate > 0)  

  ch.mod.fr<-nls(per_control ~ exp(-b*(conc+1e-6)), data=fr.ch, start = list(b=0.01))
    summary(ch.mod.fr)
  
  plot(fr.ch$conc, fr.ch$per_control, pch = 16, ylim = c(0,1), 
       xlab = 'Chlorpyrfios (ppb)', ylab = 'Reduction in predator feeding rate')
  
  alpha_q_ch_sat09<-function(In){
    exp(-2.3*In)
  }
  
  lines(dose, alpha_q_ch_sat09(dose), lty=2, col='red')
  text(0.15,0.1, labels = expression(paste(alpha,'(q) = ', 'e'^'-bq', '     b=2.30', sep='')))
 
#Toxicity to miracidia and cercariae from Hasheesh 2011 ###############
  lc50m = 0.78
  lc50c = 0.96
  lc50m.s = 1.86
  lc50c.s = 2.3
  
  pi_Mq_ch_has11<-function(In){
    Ins = In/1000
    1/(1+exp(1.86*(log(Ins)-log(0.78))))
  }  
  
  pi_Mq_ch_has11(780)
    
  in.test = seq(0,10000,1)
  
  plot(in.test, pi_Mq_ch_has11(in.test), type = 'l', xlim = c(0,100))
  
  pi_Cq_ch_has11<-function(In){
    Ins = In/1000
    1/(1+exp(2.3*(log(Ins)-log(0.96))))
  }  
  
  pi_Cq_ch_has11(960)
  
  in.test = seq(0,10000,1)
  
  plot(in.test, pi_Cq_ch_has11(in.test), type = 'l', xlim = c(0,100))