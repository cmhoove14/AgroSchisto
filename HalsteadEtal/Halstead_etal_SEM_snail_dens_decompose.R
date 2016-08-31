#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

###Load Libraries
library(lavaan)
library(AICcmodavg)
library(semPlot)
library(ggplot2)
library(MASS)
library(lattice)
st.er <- function(x) {
  sd(x)/sqrt(length(x))
}

setwd('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal')
ind <- read.csv("SEM_data.csv") # individually-transformed variables
dat<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/R_use.csv")
derp<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/snail_repro.csv")


#Create latent variables for algae from agal data ######################
alg5.2 <- '# Create latent variables for algae
Phyt.1 =~ Ph.F0.1 + Ph.F0.2 + Ph.F0.4 + Ph.F0.8 + Ph.QY.1 + Ph.QY.2 + Ph.QY.4 + Ph.QY.8
PerF0 =~ Pe.F0.1 + Pe.F0.2 + Pe.F0.4
PerQY =~ Pe.QY.1 + Pe.QY.2 + Pe.QY.4

Ph.F0.1 ~~ Ph.QY.1
Ph.F0.4 ~~ Ph.QY.4
Ph.F0.8 ~~ Ph.QY.8

# Regressions
Phyt.1 ~ At + Fe
PerF0 ~ At + Fe + AF
PerQY ~ At + Fe'
alg5.2.fit <- sem(alg5.2, data=ind)
summary(alg5.2.fit, rsq=T, standardized=T, modindices=F)
# Extract latent variables
algae <- predict(alg5.2.fit)
ind <- cbind(ind, algae)
 

# Full structural equation model ########################
mod2.2 <- '
Snails =~ TBgE1 + TBtE1 + TPhE1 + TBgH4.8 + TBtH1 + TPhH1 + BgL + BtL + PhyL
Algae <~ 1*Phyt.1 + PerF0 + PerQY

PerF0 ~~ PerQY
TBtE1 ~~ TPhE1
TBtH1 ~~ PhyL
TBgH4.8 ~~ TPhH1
TBgE1 ~~ TPhE1
TBgE1 ~~ TBtE1

Pred ~ Ch
Phyt.1 ~ At + Fe
PerF0 ~ At + Fe + AF
PerQY ~ At + Fe
Snails ~ Pred + Algae'
mod2.2.fit <- sem(mod2.2, data=ind)
summary(mod2.2.fit, rsq=T, standardized=T, modindices=F)
# Extract latent variables
latent <- predict(mod2.2.fit)
latent <- as.data.frame(latent)
ind <- cbind(ind,latent)

# Plots of linear relationships between predictor and latent variables #######################
# Predators
PredSurv <- -ind$Pred
ind <- cbind(ind,PredSurv)
ind$Pred.2 <- ind$Pred+0.08067405

#Algae  response to treatments ##########################
lm.alg<-lm(ind$Alg2.2~ind$At+ind$Fe+ind$AF)
  summary(lm.alg)
  confint(lm.alg)
  
  alg.df<-data.frame('mean'=c(mean(ind$Alg2.2[ind$At==0 & ind$Fe==0 & ind$Ch==0]), 
                              mean(ind$Alg2.2[ind$At==0 & ind$Fe==0 & ind$Ch==1]),
                              mean(ind$Alg2.2[ind$At==0 & ind$Fe==1 & ind$Ch==0]),
                              mean(ind$Alg2.2[ind$At==1 & ind$Fe==0 & ind$Ch==0]),
                              mean(ind$Alg2.2[ind$At==1 & ind$Fe==1 & ind$Ch==0])),
                     'st.err'=c(st.er(ind$Alg2.2[ind$At==0 & ind$Fe==0 & ind$Ch==0]), 
                                st.er(ind$Alg2.2[ind$At==0 & ind$Fe==0 & ind$Ch==1]),
                                st.er(ind$Alg2.2[ind$At==0 & ind$Fe==1 & ind$Ch==0]),
                                st.er(ind$Alg2.2[ind$At==1 & ind$Fe==0 & ind$Ch==0]),
                                st.er(ind$Alg2.2[ind$At==1 & ind$Fe==1 & ind$Ch==0])),
                     'Treatment'=c('Control', 'ChlorP Only', 'Fert Only', 
                                   'Atra Only', 'Fert & Atra Only'))
    alg.df$Treatment<-factor(alg.df$Treatment, levels=c('Control', 'ChlorP Only', 'Fert Only', 
                                                        'Atra Only', 'Fert & Atra Only'))
  
  ggplot(alg.df, aes(x=Treatment, y=mean, fill=Treatment)) +
    theme_bw() +
    ylab('mean algal production') +
    theme(axis.title=element_text(size=18), axis.text=element_text(size=15)) +
    scale_fill_manual(values=c('Gray', 'red',  'green', 'gold2','orange')) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Agro treatments effects on Algal production")
  
  #Does ChlorP have a negative effect on ChlorP?
    t.test(ind$Alg2.2[ind$At==0 & ind$Ch==0 & ind$Fe==0], 
           ind$Alg2.2[ind$At==0 & ind$Ch==1 & ind$Fe==0]) #Nope
  
#Proportional increases over mean algal prodcution latent variable in control  
  
  mean.control<-mean(ind$Alg2.2[ind$At==0 & ind$Fe==0 & ind$Ch==0])
  at.only<-ind$Alg2.2[ind$At==1 & ind$Fe==0 & ind$Ch==0]
  fe.only<-ind$Alg2.2[ind$At==0 & ind$Fe==1 & ind$Ch==0]
  at.fe<-ind$Alg2.2[ind$At==1 & ind$Fe==1 & ind$Ch==0]
  
  at.prop<-rep(0,10)
    for(i in 1:length(ind$Alg2.2[ind$At==1 & ind$Fe==0 & ind$Ch==0])){
      at.prop[i]<-at.only[i]/mean.control
  }
  
  fe.prop<-rep(0,10)
    for(i in 1:length(ind$Alg2.2[ind$At==0 & ind$Fe==1 & ind$Ch==0])){
      fe.prop[i]<-fe.only[i]/mean.control
  }
  
  at.fe.prop<-rep(0,5)
    for(i in 1:length(ind$Alg2.2[ind$At==1 & ind$Fe==1 & ind$Ch==0])){
      at.fe.prop[i]<-at.fe[i]/mean.control
  }
    
    alg.prop<-data.frame('mean'=c(mean(at.prop), mean(fe.prop), mean(at.fe.prop)),
                        'st.err'=c(st.er(at.prop), st.er(fe.prop), st.er(at.fe.prop)),
                     'Treatment'=c('Atra Only','Fert Only', 'Fert & Atra Only'))
  alg.prop$Treatment<-factor(alg.prop$Treatment, levels=c('Fert Only','Atra Only', 
                                                      'Fert & Atra Only'))
  
  ggplot(alg.prop, aes(x=Treatment, y=mean, fill=Treatment)) +
    theme_bw() +
    ylab('Proportion increase over mean algal production in control tanks') +
    theme(axis.title=element_text(size=18), axis.text=element_text(size=15)) +
    scale_fill_manual(values=c('green', 'gold2','darkolivegreen')) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Agro treatments effects on Algal production")
  
  
# Merge with other data sets to investigate effects of algal dynamics/pred density on raw variables ################
  #Merge with all results
    ind1<-merge(ind, dat, by.x='Tank', by.y='tank')
      for(i in 1:length(ind1$Tank)){
        ind1[i,8]<-paste(ind1[i,2],'_',ind1[i,3],'_',ind1[i,4])
      }
      colnames(ind1)[8]<-'Treats'
      
      ind1$Treats<-as.character(ind1$Treats)
      
      ind1$Treats[ind1$Treats=="0 _ 0 _ 0"]<- "Control"
      ind1$Treats[ind1$Treats=="1 _ 0 _ 0"]<- "Atra"
      ind1$Treats[ind1$Treats=="0 _ 1 _ 0"]<- "ChlorP"
      ind1$Treats[ind1$Treats=="0 _ 0 _ 1"]<- "Fert"
      ind1$Treats[ind1$Treats=="1 _ 1 _ 0"]<- "Atra_Chlor"
      ind1$Treats[ind1$Treats=="1 _ 0 _ 1"]<- "Atra_Fert"
      ind1$Treats[ind1$Treats=="0 _ 1 _ 1"]<- "Chlor_Fert"
      ind1$Treats[ind1$Treats=="1 _ 1 _ 1"]<- "All_Three"
      
      ind1$Treats<-as.factor(ind1$Treats)
      
      ind1$Treats<- factor(ind1$Treats, levels=c("Control","Atra","ChlorP",
                                               "Fert","Atra_Chlor",
                                               "Atra_Fert",
                                               "Chlor_Fert","All_Three"))
      
      ind1$At<-as.factor(ind1$At)
      ind1$Ch<-as.factor(ind1$Ch)
      ind1$Fe<-as.factor(ind1$Fe)
      
  #Merge with longitudinal reproduction results
    ind2<-merge(ind, derp, by.x='Tank', by.y='tank')


  
# Live b. truncatus at end as a function of predators ######################
  
  #fit relationship to predator mortality
  bt.liv.pred<- lm(ind1$bt_liv_fin~ind1$Pred.2)
    summary(bt.liv.pred)
    confint(bt.liv.pred)
    
  #Visualize
  plot(x=ind1$Pred.2, y=ind1$bt_liv_fin, bty="l", xlab="Predator mortality", ylab="B. truncatus alive at end",
        pch=19, bg="black", cex.axis=1.3, cex.lab=1.55)
    points(x=ind1$Pred.2[ind1$Ch==1], y=ind1$bt_liv_fin[ind1$Ch==1], col='red', pch=15, cex=1.5)
    points(x=ind1$Pred.2[ind1$At==1], y=ind1$bt_liv_fin[ind1$At==1], col='gold2', pch=17, cex=1.2)
    points(x=ind1$Pred.2[ind1$Fe==1], y=ind1$bt_liv_fin[ind1$Fe==1], col='green', pch=16)
    lines(sort(ind1$Pred.2), fitted(bt.liv.pred)[order(ind1$Pred.2)])
    legend("topleft", legend=c('None (Control)', 'Atrazine','Chlorpyrifos', 'Fertilizer'),
           pch=c(20,17,15,16), col=c('black', 'gold2','red','green'), title="AgroCs Present")
   
  
  
  
#Regression models simultaneously considering predation and bottom up effects ############
    
  #~_~_~_~B. truncatus alive at the end~_~_~_~#
  
  #Predation effects and algae effects  
  lm.bt.end1<-lm(ind1$bt_liv_fin~ind1$Pred.2+ind1$Alg2.2)
    summary(lm.bt.end1)
    confint((lm.bt.end1))
    
  predict.bt.end1<-predict(lm(ind1$bt_liv_fin~ind1$Pred.2+ind1$Alg2.2), se.fit=T)
    bt.end1.data<-data.frame('At'=ind1$At,
                             'Ch'=ind1$Ch,
                             'Fe'=ind1$Fe,
                             'Treat'=rep(0,60),
                             'bt_end'=ind1$bt_liv_fin,
                             'Alg_Prod'=ind1$Alg2.2,
                             'pred_mort'=ind1$Pred.2,
                             'prediction'=predict.bt.end1$fit,
                             'st.err_predict'=predict.bt.end1$se.fit)
    for(i in 1:60){
      bt.end1.data$Treat=paste(bt.end1.data$At,'_',bt.end1.data$Ch,'_',bt.end1.data$Fe)
    }
    
      bt.end1.data$Treat[bt.end1.data$Treat=="0 _ 0 _ 0"]<- "Control"
      bt.end1.data$Treat[bt.end1.data$Treat=="1 _ 0 _ 0"]<- "Atra"
      bt.end1.data$Treat[bt.end1.data$Treat=="0 _ 1 _ 0"]<- "ChlorP"
      bt.end1.data$Treat[bt.end1.data$Treat=="0 _ 0 _ 1"]<- "Fert"
      bt.end1.data$Treat[bt.end1.data$Treat=="1 _ 1 _ 0"]<- "Atra_Chlor"
      bt.end1.data$Treat[bt.end1.data$Treat=="1 _ 0 _ 1"]<- "Atra_Fert"
      bt.end1.data$Treat[bt.end1.data$Treat=="0 _ 1 _ 1"]<- "Chlor_Fert"
      bt.end1.data$Treat[bt.end1.data$Treat=="1 _ 1 _ 1"]<- "All_Three"
        
      bt.end1.data$Treat<-as.factor(bt.end1.data$Treat)
        
      bt.end1.data$Treat<- factor(bt.end1.data$Treat, levels=c("Control","Atra","ChlorP",
                                                   "Fert","Atra_Chlor",
                                                   "Atra_Fert",
                                                   "Chlor_Fert","All_Three"))
      
      ggplot(bt.end1.data, aes(x=Alg_Prod, y=prediction, color=Treat)) +
        theme_bw()+
        theme(axis.title=element_text(size=20),
              axis.text=element_text(size=15))+
        ylab('Predicted final B. truncatus numbers')+
        xlab('Algal production')+
        scale_color_manual(values=cbPalette)+
        geom_point(size=3.5)+
        geom_errorbar(aes(ymin=prediction-st.err_predict,
                          ymax=prediction+st.err_predict))
    
    plot(x=bt.end1.data$Alg_Prod, y=bt.end1.data$prediction, bty='l', pch=19,
         xlab="Algal Production", ylab="Predicted B. truncatus end")
    
    plot(x=bt.end1.data$pred_mort, y=bt.end1.data$prediction, bty='l', pch=19,
         xlab="Pred Mortality", ylab="Predicted B. truncatus end")
    
    #bt.end1.data$Alg_Prod<-bt.end1.data$Alg_Prod-min(bt.end1.data$Alg_Prod)
    
    wireframe(bt_liv_fin~Alg2.2+Pred.2, data=ind1,
              xlab='Algal production', ylab='Predator mortality',
              drape=TRUE, colorkey=TRUE)
  
  #Visualize lm.bt.end1 fit
    plot(x=ind1$Alg2.2[ind1$Ch==1], y=ind1$bt_liv_fin[ind1$Ch==1], bty="l", 
         xlab="algal production", ylab="B. truncatus alive at end",
         pch=15, cex=1.4, col='red', bg="black", cex.axis=1.3, cex.lab=1.55) 
    points(x=ind1$Alg2.2[ind1$Ch==1 & ind1$At==1], y=ind1$bt_liv_fin[ind1$Ch==1 & ind1$At==1], 
           col='gold2', pch=17, cex=1.2)
    points(x=ind1$Alg2.2[ind1$Ch==1 & ind1$Fe==1], y=ind1$bt_liv_fin[ind1$Ch==1 & ind1$Fe==1], 
           col='green', pch=16)
    abline(a=lm.bt.end1$coefficients[1]+mean(ind1$Pred.2)*lm.bt.end1$coefficients[2], 
           b=lm.bt.end1$coefficients[3], lty=2, lwd=2)#mean effect of predation included
    legend("topleft", legend=c('Chlorpyrifos', 'Atrazine','Fertilizer'),
           pch=c(15,17,16), col=c('red','gold2','green'), title="AgroCs Present")
    
  #Skip algal effects, look at Atrazine/Fert controlling for predation  
  lm.bt.end2<-lm(ind1$bt_liv_fin~ind1$Pred.2+ind1$At*ind1$Fe)  
    summary(lm.bt.end2) 
  
  #Visualize lm.bt.end2 fit
    plot(x=ind1$Alg2.2, y=ind1$bt_liv_fin, bty="l", 
          xlab="algal production within", ylab="B. truncatus alive at end",
          pch=15, cex=1.4, col='red', bg="black", cex.axis=1.3, cex.lab=1.55) 
      points(x=ind1$Alg2.2[ind1$Ch==1 & ind1$At==1], y=ind1$bt_liv_fin[ind1$Ch==1 & ind1$At==1], 
             col='gold2', pch=17, cex=1.2)
      points(x=ind1$Alg2.2[ind1$Ch==1 & ind1$Fe==1], y=ind1$bt_liv_fin[ind1$Ch==1 & ind1$Fe==1], 
             col='green', pch=16)
      abline(a=0, b=lm.bt.end2$coefficients[2], lty=2, lwd=2)
      legend("topleft", legend=c('Chlorpyrifos', 'Atrazine','Fertilizer'),
             pch=c(15,17,16), col=c('red','gold2','green'), title="AgroCs Present")
    
  #Effects of Atra/Fert within ChlorP studies   
  lm.bt.end3<-lm(ind1$bt_liv_fin[ind1$Ch==1]~ind1$At[ind1$Ch==1]*ind1$Fe[ind1$Ch==1])
    summary(lm.bt.end3) #Only intercept is significant (aka snail numbers when chlorpyrifos is present)
    
  #Visualize lm.bt.end3 fit  
    plot(x=ind1$Alg2.2[ind1$Ch==1], y=ind1$bt_liv_fin[ind1$Ch==1], bty="l", 
       xlab="algal production within ChlorP treatments", ylab="B. truncatus alive at end",
       pch=15, cex=1.4, col='red', bg="black", cex.axis=1.3, cex.lab=1.55) 
    points(x=ind1$Alg2.2[ind1$Ch==1 & ind1$At==1], y=ind1$bt_liv_fin[ind1$Ch==1 & ind1$At==1], 
           col='gold2', pch=17, cex=1.2)
    points(x=ind1$Alg2.2[ind1$Ch==1 & ind1$Fe==1], y=ind1$bt_liv_fin[ind1$Ch==1 & ind1$Fe==1], 
           col='green', pch=16)
    abline(a=lm.bt.end3$coefficients[1], b=lm.bt.end3$coefficients[2], lty=2, lwd=2)
    
  #~_~_~_~B. truncatus sampled hatchlings~_~_~_~#
  
  plot(x=ind1$Pred.2, y=ind1$bt_hatch, bty="l", 
       xlab="Predator mortality", ylab="B. truncatus hatchlings sampled",
       pch=19, bg="black", cex.axis=1.3, cex.lab=1.55)
    points(x=ind1$Pred.2[ind1$Ch==1], y=ind1$bt_hatch[ind1$Ch==1], col='red', pch=15, cex=1.5)
    points(x=ind1$Pred.2[ind1$At==1], y=ind1$bt_hatch[ind1$At==1], col='gold2', pch=17, cex=1.2)
    points(x=ind1$Pred.2[ind1$Fe==1], y=ind1$bt_hatch[ind1$Fe==1], col='green', pch=16)
    legend("topleft", legend=c('None (Control)', 'Atrazine','Chlorpyrifos', 'Fertilizer'),
           pch=c(20,17,15,16), col=c('black', 'gold2','red','green'), title="AgroCs Present") 
  
  #Predator effects only, accounting for adult (i.e. reproducing) population size    
  lm.bt.hat1<-lm(ind1$bt_hatch~ind1$Pred.2+ind1$bt_liv_fin)
    summary(lm.bt.hat1)
  #Add regression line to plot; intercept insignificant, coded instead as mean effect of Bt_liv_fin (slightly significant)  
    abline(a=mean(ind1$bt_liv_fin)*lm.bt.hat1$coefficients[3], 
           b=lm.bt.hat1$coefficients[2], lty=2, lwd=2)
    
  #Algal effects only, accounting for adult (i.e. reproducing) population size    
  lm.bt.hat1.2<-lm(ind1$bt_hatch~ind1$Alg2.2+ind1$bt_liv_fin)
    summary(lm.bt.hat1.2) #NOPE
    
  #Predation effects along with algal effects, accounting for reproducing population size  
  lm.bt.hat2<-lm(ind1$bt_hatch~ind1$Pred.2+ind1$Alg2.2+ind1$total_bt_fin)
    summary(lm.bt.hat2) 
    confint(lm.bt.hat2)
    
  plot(ind1$Alg2.2, ind1$bt_hatch, bty='l', 
       xlab='Algal production', ylab='B. truncatus hatchlings sampled',
       pch=19, bg="black", cex.axis=1.3, cex.lab=1.55)
    points(x=ind1$Alg2.2[ind1$Ch==1], y=ind1$bt_hatch[ind1$Ch==1], col='red', pch=15, cex=1.5)
    points(x=ind1$Alg2.2[ind1$At==1], y=ind1$bt_hatch[ind1$At==1], col='gold2', pch=17, cex=1.2)
    points(x=ind1$Alg2.2[ind1$Fe==1], y=ind1$bt_hatch[ind1$Fe==1], col='green', pch=16)
    legend("topleft", legend=c('None (Control)', 'Atrazine','Chlorpyrifos', 'Fertilizer'),
           pch=c(20,17,15,16), col=c('black', 'gold2','red','green'), title="AgroCs Present") 
    abline(a=lm.bt.hat2$coefficients[1]+mean(ind1$Pred.2)*lm.bt.hat2$coefficients[2],
           b=lm.bt.hat2$coefficients[3], lty=2, lwd=2)
    
  #See if we can skip algal production to assess "direct" effects of Atra/Fert on hatchlings   
  lm.bt.hat3<-lm(ind1$bt_hatch~ind1$Pred.2+ind1$At*ind1$Fe +ind1$bt_liv_fin)  
    summary(lm.bt.hat3) 
    
  #~_~_~_~B. truncatus sampled hatchlings~_~_~_~#
    
  plot(x=ind1$Pred.2, y=ind1$bt_eggs, bty="l", 
         xlab="Predator mortality", ylab="B. truncatus eggs sampled",
         pch=19, bg="black", cex.axis=1.3, cex.lab=1.55)
    points(x=ind1$Pred.2[ind1$Ch==1], y=ind1$bt_eggs[ind1$Ch==1], col='red', pch=15, cex=1.5)
    points(x=ind1$Pred.2[ind1$At==1], y=ind1$bt_eggs[ind1$At==1], col='gold2', pch=17, cex=1.2)
    points(x=ind1$Pred.2[ind1$Fe==1], y=ind1$bt_eggs[ind1$Fe==1], col='green', pch=16)
    legend("topleft", legend=c('None (Control)', 'Atrazine','Chlorpyrifos', 'Fertilizer'),
           pch=c(20,17,15,16), col=c('black', 'gold2','red','green'), title="AgroCs Present")
    
  #Predator effects only, accounting for adult (i.e. reproducing) population size    
    lm.bt.egg1<-lm(ind1$bt_eggs~ind1$Pred.2+ind1$bt_liv_fin)
      summary(lm.bt.egg1)
  #Add regression line to plot; intercept insignificant, coded instead as mean effect of Bt_liv_fin (slightly significant)  
    abline(a=lm.bt.egg1$coefficients[1], 
           b=lm.bt.hat1$coefficients[2], lty=2, lwd=2)
    
  #Algal effects only, accounting for adult (i.e. reproducing) population size    
    lm.bt.egg1.2<-lm(ind1$bt_eggs~ind1$Alg2.2+ind1$bt_liv_fin)
      summary(lm.bt.egg1.2) #Intercept significant but so what
    
  #Predation effects along with algal effects, accounting for reproducing population size  
    lm.bt.egg2<-lm(ind1$bt_eggs~ind1$Pred.2+ind1$Alg2.2+ind1$total_bt_fin)
      summary(lm.bt.egg2) 
      confint(lm.bt.egg2)
      
    plot(ind1$Alg2.2, ind1$bt_eggs, bty='l', 
           xlab='Algal production', ylab='B. truncatus eggs sampled',
           pch=19, bg="black", cex.axis=1.3, cex.lab=1.55)
      points(x=ind1$Alg2.2[ind1$Ch==1], y=ind1$bt_eggs[ind1$Ch==1], col='red', pch=15, cex=1.5)
      points(x=ind1$Alg2.2[ind1$At==1], y=ind1$bt_eggs[ind1$At==1], col='gold2', pch=17, cex=1.2)
      points(x=ind1$Alg2.2[ind1$Fe==1], y=ind1$bt_eggs[ind1$Fe==1], col='green', pch=16)
      legend("topleft", legend=c('None (Control)', 'Atrazine','Chlorpyrifos', 'Fertilizer'),
             pch=c(20,17,15,16), col=c('black', 'gold2','red','green'), title="AgroCs Present") 
      abline(a=lm.bt.egg2$coefficients[1]+mean(ind1$Pred.2)*lm.bt.egg2$coefficients[2],
             b=lm.bt.egg2$coefficients[3], lty=2, lwd=2)
      
  #See if we can skip algal production to assess "direct" effects of Atra/Fert on hatchlings   
    lm.bt.egg3<-lm(ind1$bt_eggs~ind1$Pred.2+ind1$At*ind1$Fe +ind1$bt_liv_fin)  
      summary(lm.bt.egg3)  
      
    #Add this line to plot
      abline(a=mean(ind1$Pred.2)*lm.bt.egg3$coefficients[2],
             b=lm.bt.egg3$coefficients[4], lty=2, lwd=2, col='darkgreen')
      

#B. truncatus hatchlings sampled during 12 weeks as a function of predators##################
  #fit relationship to predator mortality
    bt.hatch.pred<- lm(ind1$bt_hatch~ind1$Pred.2)
      summary(bt.hatch.pred)
    
  #Visualize
    plot(x=ind1$Pred.2, y=ind1$bt_hatch, bty="l", xlab="Predator mortality", ylab="B. truncatus hatchlings",
         pch=19, bg="black", cex.axis=1.3, cex.lab=1.55)
      points(x=ind1$Pred.2[ind1$Ch==1], y=ind1$bt_hatch[ind1$Ch==1], pch=1, col='red')
      lines(sort(ind1$Pred.2), fitted(bt.hatch.pred)[order(ind1$Pred.2)])
      legend("topleft", legend=c('ChlorP Present'), pch=1, col='red')
      
    #Get residuals of pred/ b. truncatus realtionship
      ind1$bt.hatch.pre.resid<-resid(lm(ind1$bt_hatch~ind1$PredSurv))
    
    #Fit function to algae/ residuals of predator relationship
      bt.hatch.resid.fit<-lm(ind1$bt.hatch.pre.resid~ind1$Alg2.2+I(ind1$Alg2.2^2))
    
    #Visualize
      plot(ind1$Alg2.2, ind1$bt.hatch.pre.resid, pch=19, bty="l", xlab="Algal production",
           ylab="Residual B. truncatus hatchlings", bg="black", cex.axis=1.3, cex.lab=1.55)
        points(x=ind1$Alg2.2[ind1$Ch==1], y=ind1$bt.hatch.pre.resid[ind1$Ch==1], col='red', pch=15, cex=1.5)
        points(x=ind1$Alg2.2[ind1$At==1], y=ind1$bt.hatch.pre.resid[ind1$At==1], col='gold2', pch=17, cex=1.5)
        points(x=ind1$Alg2.2[ind1$Fe==1], y=ind1$bt.hatch.pre.resid[ind1$Fe==1], col='green', pch=16)
        legend("topleft", legend=c('None (Control)', 'Atrazine','Chlorpyrifos', 'Fertilizer'),
               pch=c(20,17,15,16), col=c('black', 'gold2','red','green'), title="AgroCs Present")
        lines(sort(ind$Alg2.2), fitted(bt.hatch.resid.fit)[order(ind$Alg2.2)])
    
    #Get rid of negative residuals to calculate means and standard errors    
    ind1$bt.hatch.pre.resid<-ind1$bt.hatch.pre.resid-min(ind1$bt.hatch.pre.resid)
    
    #Anova to see if there are significant differences in resid snail density between treatments   
    anv.bt.hatch<-aov(bt.hatch.pre.resid~Treats, data=ind1)
    summary(anv.bt.hatch)
    boxplot(bt.hatch.pre.resid~Treats, data=ind1, ylab="Residual B. truncatus hatchlings", 
            xlab="Treatments", cex.lab=1.5, cex.axis=1.2)
    
    #calc means/st. dev for treatment presence/absence
    atra.bt.hatch.1<-mean(ind1$bt.hatch.pre.resid[ind1$At==1])
      atra.bt.hatch.1.se<-st.er(ind1$bt.hatch.pre.resid[ind1$At==1])
    atra.bt.hatch.0<-mean(ind1$bt.hatch.pre.resid[ind1$At==0])
      atra.bt.hatch.0.se<-st.er(ind1$bt.hatch.pre.resid[ind1$At==0])
    fert.bt.hatch.1<-mean(ind1$bt.hatch.pre.resid[ind1$Fe==1])
      fert.bt.hatch.1.se<-st.er(ind1$bt.hatch.pre.resid[ind1$Fe==1])
    fert.bt.hatch.0<-mean(ind1$bt.hatch.pre.resid[ind1$Fe==0])
      fert.bt.hatch.0.se<-st.er(ind1$bt.hatch.pre.resid[ind1$Fe==0])
    atra.fert.bt.hatch.1<-mean(ind1$bt.hatch.pre.resid[ind1$Fe==1 & ind1$At==1])
      atra.fert.bt.hatch.1.se<-st.er(ind1$bt.hatch.pre.resid[ind1$Fe==1 & ind1$At==1])
    
    bar.bt.hatch.dat<-data.frame('Treatment'=c('Fert_Absent', 'Atra_Absent', 'Fert_Present', 'Atra_Present', 'Both_Present'),
                               'mean'=c(atra.bt.hatch.0, fert.bt.hatch.0, atra.bt.hatch.1, fert.bt.hatch.1, atra.fert.bt.hatch.1),
                               'st.err'=c(atra.bt.hatch.0.se, fert.bt.hatch.0.se, atra.bt.hatch.1.se, 
                                          fert.bt.hatch.1.se, atra.fert.bt.hatch.1.se)) 
    bar.bt.hatch.dat$Treatment<-factor(bar.bt.hatch.dat$Treatment, levels=c('Fert_Absent', 'Atra_Absent', 'Fert_Present', 
                                                                        'Atra_Present', 'Both_Present'))
    
    bar.cols<-c('gray75', 'gray75', 'green2', 'gold2',  'darkolivegreen')
    
    ggplot(bar.bt.hatch.dat, aes(x=Treatment, y=mean, fill=Treatment)) +
      theme_bw()+
      theme(axis.title.x=element_text(size=20),
            axis.title.y=element_text(size=20),
            axis.text=element_text(size=15)) +
      ylab('Residual B. truncatus hatchling counts +/- St.err') +
      scale_fill_manual(values=bar.cols) +
      geom_bar(position=position_dodge(), stat="identity", width=.7) +
      geom_errorbar(aes(ymin=mean-st.err,
                        ymax=mean+st.err),
                    width=.2, position=position_dodge(.7))
    
    t.test(y=ind1$bt.hatch.pre.resid[ind1$At==0], x=ind1$bt.hatch.pre.resid[ind1$At==1], alternative = 'greater')#p=0.0963    
    t.test(y=ind1$bt.hatch.pre.resid[ind1$Fe==0], x=ind1$bt.hatch.pre.resid[ind1$Fe==1], alternative = 'greater')#p=0.0286
    t.test(y=ind1$bt.hatch.pre.resid[ind1$Fe==1 & ind1$At==1], 
           x=ind1$bt.hatch.pre.resid[ind1$Fe==1], alternative = 'greater')#p=0.6348

#B. truncatus egg masses as a function of treatments/ predation ########################
  #fit relationship to predator mortality
    bt.eggs.pred<- lm(ind1$bt_eggs~ind1$Pred.2)
      summary(bt.eggs.pred)
    
  #Visualize
    plot(x=ind1$Pred.2, y=ind1$bt_eggs, bty="l", xlab="Predator mortality", ylab="B. truncatus egg masses",
         pch=19, bg="black", cex.axis=1.3, cex.lab=1.55)
    points(x=ind1$Pred.2[ind1$Ch==1], y=ind1$bt_eggs[ind1$Ch==1], pch=1, col='red')
    lines(sort(ind1$Pred.2), fitted(bt.eggs.pred)[order(ind1$Pred.2)])
    legend("topleft", legend=c('ChlorP Present'), pch=1, col='red')
    
  #Get residuals of pred/ b. truncatus realtionship
    ind1$bt.eggs.pre.resid<-resid(lm(ind1$bt_eggs~ind1$PredSurv))
    
  #Fit function to algae/ residuals of predator relationship
    bt.eggs.resid.fit<-lm(ind1$bt.eggs.pre.resid~ind1$Alg2.2+I(ind1$Alg2.2^2))
    
  #Visualize
    plot(ind1$Alg2.2, ind1$bt.eggs.pre.resid, pch=19, bty="l", xlab="Algal production",
         ylab="Residual B. truncatus egg masses", bg="black", cex.axis=1.3, cex.lab=1.55)
      points(x=ind1$Alg2.2[ind1$Ch==1], y=ind1$bt.eggs.pre.resid[ind1$Ch==1], col='red', pch=15, cex=1.5)
      points(x=ind1$Alg2.2[ind1$At==1], y=ind1$bt.eggs.pre.resid[ind1$At==1], col='gold2', pch=17, cex=1.5)
      points(x=ind1$Alg2.2[ind1$Fe==1], y=ind1$bt.eggs.pre.resid[ind1$Fe==1], col='green', pch=16)
      legend("topleft", legend=c('None (Control)', 'Atrazine','Chlorpyrifos', 'Fertilizer'),
             pch=c(20,17,15,16), col=c('black', 'gold2','red','green'), title="AgroCs Present")
      lines(sort(ind$Alg2.2), fitted(bt.eggs.resid.fit)[order(ind$Alg2.2)])
      
    #Get rid of negative residuals to calculate means and standard errors    
    ind1$bt.eggs.pre.resid<-ind1$bt.eggs.pre.resid-min(ind1$bt.eggs.pre.resid)
    
    #Anova to see if there are significant differences in resid snail density between treatments   
    anv.bt.eggs<-aov(bt.eggs.pre.resid~Treats, data=ind1)
    summary(anv.bt.eggs)
    boxplot(bt.eggs.pre.resid~Treats, data=ind1, ylab="Residual B. truncatus eggslings", 
            xlab="Treatments", cex.lab=1.5, cex.axis=1.2)
    
    #calc means/st. dev for treatment presence/absence
    atra.bt.eggs.1<-mean(ind1$bt.eggs.pre.resid[ind1$At==1])
      atra.bt.eggs.1.se<-st.er(ind1$bt.eggs.pre.resid[ind1$At==1])
    atra.bt.eggs.0<-mean(ind1$bt.eggs.pre.resid[ind1$At==0])
      atra.bt.eggs.0.se<-st.er(ind1$bt.eggs.pre.resid[ind1$At==0])
    fert.bt.eggs.1<-mean(ind1$bt.eggs.pre.resid[ind1$Fe==1])
      fert.bt.eggs.1.se<-st.er(ind1$bt.eggs.pre.resid[ind1$Fe==1])
    fert.bt.eggs.0<-mean(ind1$bt.eggs.pre.resid[ind1$Fe==0])
      fert.bt.eggs.0.se<-st.er(ind1$bt.eggs.pre.resid[ind1$Fe==0])
    atra.fert.bt.eggs.1<-mean(ind1$bt.eggs.pre.resid[ind1$Fe==1 & ind1$At==1])
      atra.fert.bt.eggs.1.se<-st.er(ind1$bt.eggs.pre.resid[ind1$Fe==1 & ind1$At==1])
    
    bar.bt.eggs.dat<-data.frame('Treatment'=c('Fert_Absent', 'Atra_Absent', 'Fert_Present', 'Atra_Present', 'Both_Present'),
                                 'mean'=c(atra.bt.eggs.0, fert.bt.eggs.0, atra.bt.eggs.1, fert.bt.eggs.1, atra.fert.bt.eggs.1),
                                 'st.err'=c(atra.bt.eggs.0.se, fert.bt.eggs.0.se, atra.bt.eggs.1.se, 
                                            fert.bt.eggs.1.se, atra.fert.bt.eggs.1.se)) 
    bar.bt.eggs.dat$Treatment<-factor(bar.bt.eggs.dat$Treatment, levels=c('Fert_Absent', 'Atra_Absent', 'Fert_Present', 
                                                                            'Atra_Present', 'Both_Present'))
    
    bar.cols<-c('gray75', 'gray75', 'green2', 'gold2',  'darkolivegreen')
    
    ggplot(bar.bt.eggs.dat, aes(x=Treatment, y=mean, fill=Treatment)) +
      theme_bw()+
      theme(axis.title.x=element_text(size=20),
            axis.title.y=element_text(size=20),
            axis.text=element_text(size=15)) +
      ylab('Residual B. truncatus egg masses +/- St.err') +
      scale_fill_manual(values=bar.cols) +
      geom_bar(position=position_dodge(), stat="identity", width=.7) +
      geom_errorbar(aes(ymin=mean-st.err,
                        ymax=mean+st.err),
                    width=.2, position=position_dodge(.7))
    
    t.test(y=ind1$bt.eggs.pre.resid[ind1$At==0], x=ind1$bt.eggs.pre.resid[ind1$At==1], alternative = 'greater')#p=0.0963    
    t.test(y=ind1$bt.eggs.pre.resid[ind1$Fe==0], x=ind1$bt.eggs.pre.resid[ind1$Fe==1], alternative = 'greater')#p=0.0286
    t.test(y=ind1$bt.eggs.pre.resid[ind1$Fe==1 & ind1$At==1], 
           x=ind1$bt.eggs.pre.resid[ind1$Fe==1], alternative = 'greater')#p=0.6348
    
#Carry on with same process of analysis as in Neal's code ############################    
  
  #Get residuals of pred/ b. truncatus relationship
  ind1$bt.liv.pre.resid<-resid(bt.liv.pred)#Chris changed to refer to lm from above instead of new linear model
  
  #Fit function to algae/ residuals of predator relationship
  bt.liv.resid.fit<-lm(ind1$bt.liv.pre.resid~ind1$Alg2.2+I(ind1$Alg2.2^2))

  #Visualize
  plot(ind1$Alg2.2, ind1$bt.liv.pre.resid, pch=19, bty="l", xlab="Algal production",
       ylab="Residual B. truncatus finals", bg="black", cex.axis=1.3, cex.lab=1.55)
  points(x=ind1$Alg2.2[ind1$Ch==1], y=ind1$bt.liv.pre.resid[ind1$Ch==1], col='red', pch=15, cex=1.5)
  points(x=ind1$Alg2.2[ind1$At==1], y=ind1$bt.liv.pre.resid[ind1$At==1], col='gold2', pch=17, cex=1.5)
  points(x=ind1$Alg2.2[ind1$Fe==1], y=ind1$bt.liv.pre.resid[ind1$Fe==1], col='green', pch=16)
  legend("topleft", legend=c('None (Control)', 'Atrazine','Chlorpyrifos', 'Fertilizer'),
         pch=c(20,17,15,16), col=c('black', 'gold2','red','green'), title="AgroCs Present")
    lines(sort(ind$Alg2.2), fitted(bt.liv.resid.fit)[order(ind$Alg2.2)])
  
  #Get rid of negative residuals to calculate means and standard errors    
  ind1$bt.liv.pre.resid<-ind1$bt.liv.pre.resid-min(ind1$bt.liv.pre.resid)
    
  #Anova to see if there are significant differences in resid snail density between treatments   
    anv.bt.liv<-aov(bt.liv.pre.resid~Treats, data=ind1)
      summary(anv.bt.liv)
    boxplot(bt.liv.pre.resid~Treats, data=ind1, ylab="Residual Final B. truncatus", 
            xlab="Treatments", cex.lab=1.5, cex.axis=1.2)
    
    #calc means/st. dev for treatment presence/absence
    atra.bt.liv.1<-mean(ind1$bt.liv.pre.resid[ind1$At==1])
      atra.bt.liv.1.se<-st.er(ind1$bt.liv.pre.resid[ind1$At==1])
    atra.bt.liv.0<-mean(ind1$bt.liv.pre.resid[ind1$At==0])
      atra.bt.liv.0.se<-st.er(ind1$bt.liv.pre.resid[ind1$At==0])
    fert.bt.liv.1<-mean(ind1$bt.liv.pre.resid[ind1$Fe==1])
      fert.bt.liv.1.se<-st.er(ind1$bt.liv.pre.resid[ind1$Fe==1])
    fert.bt.liv.0<-mean(ind1$bt.liv.pre.resid[ind1$Fe==0])
      fert.bt.liv.0.se<-st.er(ind1$bt.liv.pre.resid[ind1$Fe==0])
    atra.fert.bt.liv.1<-mean(ind1$bt.liv.pre.resid[ind1$Fe==1 & ind1$At==1])
      atra.fert.bt.liv.1.se<-st.er(ind1$bt.liv.pre.resid[ind1$Fe==1 & ind1$At==1])
    
    bar.bt.liv.dat<-data.frame('Treatment'=c('Fert_Absent', 'Atra_Absent', 'Fert_Present', 'Atra_Present', 'Both_Present'),
                        'mean'=c(atra.bt.liv.0, fert.bt.liv.0, atra.bt.liv.1, fert.bt.liv.1, atra.fert.bt.liv.1),
                        'st.err'=c(atra.bt.liv.0.se, fert.bt.liv.0.se, atra.bt.liv.1.se, 
                                   fert.bt.liv.1.se, atra.fert.bt.liv.1.se)) 
    bar.bt.liv.dat$Treatment<-factor(bar.bt.liv.dat$Treatment, levels=c('Fert_Absent', 'Atra_Absent', 'Fert_Present', 
                                                          'Atra_Present', 'Both_Present'))
    
    bar.cols<-c('gray75', 'gray75', 'green2', 'gold2',  'darkolivegreen')
    
    ggplot(bar.bt.liv.dat, aes(x=Treatment, y=mean, fill=Treatment)) +
      theme_bw()+
      theme(axis.title.x=element_text(size=20),
            axis.title.y=element_text(size=20),
            axis.text=element_text(size=15)) +
      ylab('Residual B. truncatus final counts +/- St.err') +
      scale_fill_manual(values=bar.cols) +
      geom_bar(position=position_dodge(), stat="identity", width=.7) +
      geom_errorbar(aes(ymin=mean-st.err,
                        ymax=mean+st.err),
                    width=.2, position=position_dodge(.7))
    