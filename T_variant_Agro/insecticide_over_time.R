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
#Model to observe expected insecticide concentrations over time in hypothetical application scenario
decay<-data.frame(time = c(1:365),
                  mal = 0,
                  chlor = 0,
                  terb = 0,
                  lamcy = 0,
                  esfen = 0,
                  perm = 0)

#calculate decay rates (k) from hydrolysis half lives for each chemical in Halstead 2015 (from table S1 and pmep.cce.cornell.edu)
  mal.k = -log(0.5)/6.2         #from table s1 and in agreement of "less than 1 week in raw river water" from Cornell
  chlor.k = -log(0.5)/25.5      #from table S1; within the range of reported half life from cornell
  terb.k = -log(0.5)/6.5        #from table s1; in agreement with Cornell estimate of 5.5 days at pH of 7; degrades into formaldehyde
  lamcy.k = -log(0.5)/1         #very fast; 0 in table S1; "Not expected to be prevalent in surface waters" according to cornell website
  esfen.k = -log(0.5)/10        #cornell 4-15 days half life in water
  perm.k = -log(0.5)/2          #cornell "half life of less than 2.5 days"

#at time=0, add agrochemical at concentration =to median EEC from Halstead 2015
  med.mal = 0.778
  med.chlor = 5.810
  med.terb = 1.435
  med.lamcy = 0.649
  med.esfen = 0.311
  med.perm = 1.420
  
  decay$mal[1] = med.mal
  decay$chlor[1] = med.chlor
  decay$terb[1] = med.terb
  decay$lamcy[1] = med.lamcy
  decay$esfen[1] = med.esfen
  decay$perm[1] = med.perm
  
#fill chemical concentration over the year based on application reccomendations in Halstead table S1
  for(i in 2:nrow(decay)){
    if(i == 5){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 6){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k + med.mal
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k 
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k + med.esfen
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k + med.perm
    }
    else if(i == 9){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 11){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k + med.chlor
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k + med.esfen
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k + med.perm
    }
    else if(i == 13){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 16){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k + med.esfen
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k + med.perm
    }
    else if(i == 17){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 21){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k + med.chlor
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k + med.esfen
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k + med.perm
    }
    else if(i == 25){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 26){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k + med.perm
    }
    else if(i == 29){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 31){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k + med.perm
    }
    else if(i == 33){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 36){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k + med.perm
    }
    else if(i == 37){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 41){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 45){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 49){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 53){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 57){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else if(i == 61){
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k + med.lamcy
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
    else{
      decay[i,2] = decay[(i-1),2] - decay[(i-1),2]*mal.k
      decay[i,3] = decay[(i-1),3] - decay[(i-1),3]*chlor.k
      decay[i,4] = decay[(i-1),4] - decay[(i-1),4]*terb.k
      decay[i,5] = decay[(i-1),5] - decay[(i-1),5]*lamcy.k
      decay[i,6] = decay[(i-1),6] - decay[(i-1),6]*esfen.k
      decay[i,7] = decay[(i-1),7] - decay[(i-1),7]*perm.k
    }
  }

plot(x=decay$time,y=decay$chlor, xlab = 'time', ylab = 'insecticide concentration (ppb)',
     type = 'l', lwd=2, col='red') 
  lines(x=decay$time, y=decay$mal, lwd=2, col = 'blue')
  lines(x=decay$time, y=decay$perm, lwd=2, col = 'orange')
  lines(x=decay$time, y=decay$terb, lwd=2, col = 'purple')
  lines(x=decay$time, y=decay$lamcy, lwd=2, col = 'green')
  lines(x=decay$time, y=decay$esfen, lwd=2, col = 'black')
  legend('topright', legend = c('Chlorpyrifos', 'Malathion', 'Terbufos',
                                'L_Cyhalothrin', 'Esfenvalerate', 'Permethrin'),
         lwd = 2, col = c('red', 'blue', 'purple', 'green', 'black', 'orange'), cex = 0.6)
#REad in modeled toxicity to crustaceans #################
data<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/Cray.10d.LC50.2012.act2.csv')
  mal<-subset(data, Chem == 'Mal')
  mal.c<-unique(mal$Conc)
  mal.d<-as.numeric()
  for(i in 1:length(mal.c)){
    mal.d[i] = sum(mal$Dead[mal$Conc == mal.c[i]])/5
  }
  
  mal_mod<-drm(Dead ~ Conc, data = mal, fct = LL2.2())
  
  chlor<-subset(data, Chem == 'Chlor')
  chlor.c<-unique(chlor$Conc)
  chlor.d<-as.numeric()
  for(i in 1:length(chlor.c)){
    chlor.d[i] = sum(chlor$Dead[chlor$Conc == chlor.c[i]])/5
  }
  
  chlor_mod<-drm(Dead ~ Conc, data = chlor, fct = LL2.2())
  
  terb<-subset(data, Chem == 'Terb')
  terb.c<-unique(terb$Conc)
  terb.d<-as.numeric()
  for(i in 1:length(terb.c)){
    terb.d[i] = sum(terb$Dead[terb$Conc == terb.c[i]])/5
  }
  
  terb_mod<-drm(Dead ~ Conc, data = terb, fct = LL2.2())
  
  lamcy<-subset(data, Chem == 'Lambda')
  lamcy.c<-unique(lamcy$Conc)
  lamcy.d<-as.numeric()
  for(i in 1:length(lamcy.c)){
    lamcy.d[i] = sum(lamcy$Dead[lamcy$Conc == lamcy.c[i]])/5
  }
  
  lamcy_mod<-drm(Dead ~ Conc, data = lamcy, fct = LL2.2())
  
  esfen<-subset(data, Chem == 'Esfen')
  esfen.c<-unique(esfen$Conc)
  esfen.d<-as.numeric()
  for(i in 1:length(esfen.c)){
    esfen.d[i] = sum(esfen$Dead[esfen$Conc == esfen.c[i]])/5
  }
  
  esfen_mod<-drm(Dead ~ Conc, data = esfen, fct = LL2.2())
  
  perm<-subset(data, Chem == 'Perm')
  perm.c<-unique(perm$Conc)
  perm.d<-as.numeric()
  for(i in 1:length(perm.c)){
    perm.d[i] = sum(perm$Dead[perm$Conc == perm.c[i]])/5
  }
  
  perm_mod<-drm(Dead ~ Conc, data = perm, fct = LL2.2())

#Fill decay data frame with insecticide toxicity ###############
  m2<-data.frame(Conc = decay$mal)
  decay$mal_muPq = predict(mal_mod, m2)/10
  
  c2<-data.frame(Conc = decay$chlor)
  decay$chlor_muPq = predict(chlor_mod, c2)/10
  
  t2<-data.frame(Conc = decay$terb)
  decay$terb_muPq = predict(terb_mod, t2)/10
  
  l2<-data.frame(Conc = decay$lamcy)
  decay$lamcy_muPq = predict(lamcy_mod, l2)/10
  
  e2<-data.frame(Conc = decay$esfen)
  decay$esfen_muPq = predict(esfen_mod, e2)/10
  
  p2<-data.frame(Conc = decay$perm)
  decay$perm_muPq = predict(perm_mod, p2)/10
  
#Fill decay data frame with R0 estimates #############
  for(i in 1:nrow(decay)){
    decay[i,14] = get_Ro_beta_lamda(muPq = decay[i,8],
                                    beta = beta.use,
                                    lamda = lamda.use)[3]
    decay[i,15] = get_Ro_beta_lamda(muPq = decay[i,9],
                                    beta = beta.use,
                                    lamda = lamda.use)[3]
    decay[i,16] = get_Ro_beta_lamda(muPq = decay[i,10],
                                    beta = beta.use,
                                    lamda = lamda.use)[3]
    decay[i,17] = get_Ro_beta_lamda(muPq = decay[i,11],
                                    beta = beta.use,
                                    lamda = lamda.use)[3]
    decay[i,18] = get_Ro_beta_lamda(muPq = decay[i,12],
                                    beta = beta.use,
                                    lamda = lamda.use)[3]
    decay[i,19] = get_Ro_beta_lamda(muPq = decay[i,13],
                                    beta = beta.use,
                                    lamda = lamda.use)[3]
  }
  
  colnames(decay)[14:19] = c('mal_r0', 'chlor_r0', 'terb_r0',
                             'lamcy_r0', 'esfen_r0', 'perm_r0')
  
plot(x=decay$time,y=decay$chlor_r0, xlab = 'time', ylab = 'R0 estimates', ylim = c(0,4),
       type = 'l', lwd=2, col='red', main = 'organophosphates') 
  lines(x=decay$time, y=decay$mal_r0+0.025, lwd=2, col = 'blue')
  lines(x=decay$time, y=decay$terb_r0, lwd=2, col = 'purple')
  legend('topright', legend = c('Chlorpyrifos', 'Malathion', 'Terbufos'),
         lwd = 2, col = c('red', 'blue', 'purple'), cex = 0.6)
  
plot(x=decay$time, y=decay$perm_r0, lwd=2, col = 'orange', ylim = c(0,4), type = 'l',
      xlab = 'time', ylab = 'R0 estimates', main = 'pyrethroids')
  lines(x=decay$time, y=decay$lamcy_r0, lwd=2, col = 'green')
  lines(x=decay$time, y=decay$esfen_r0, lwd=2, col = 'black')
  legend('topright', legend = c('L_Cyhalothrin', 'Esfenvalerate', 'Permethrin'),
         lwd = 2, col = c('green', 'black', 'orange'), cex = 0.6)
  
r0.mal = mean(decay$mal_r0)
r0.chlor = mean(decay$chlor_r0)
r0.terb = mean(decay$terb_r0)
r0.lamcy = mean(decay$lamcy_r0)
r0.esfen = mean(decay$esfen_r0)
r0.perm = mean(decay$perm_r0)

mal.days = length(decay$time[decay$mal_r0 >= 1])
chlor.days = length(decay$time[decay$chlor_r0 >= 1])
terb.days = length(decay$time[decay$terb_r0 >= 1])
lamcy.days = length(decay$time[decay$lamcy_r0 >= 1])
esfen.days = length(decay$time[decay$esfen_r0 >= 1])
perm.days = length(decay$time[decay$perm_r0 >= 1])

barplot(height = c(chlor.days, mal.days, terb.days, lamcy.days, esfen.days, perm.days),
        names.arg = c('Chlorpyrifos', 'Malathion', 'Terbufos',
                       'L_Cyhalothrin', 'Esfenvalerate', 'Permethrin'),
        col = c('red', 'blue', 'purple', 'green', 'black', 'orange'), horiz = TRUE, xlab = 'Days with R0>1', xlim = c(0,50))