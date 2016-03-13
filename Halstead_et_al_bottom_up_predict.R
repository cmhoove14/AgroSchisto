#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

###Load Libraries and data #############
library(lavaan)
library(AICcmodavg)
library(semPlot)
library(ggplot2)
library(MASS)
library(lattice)
library(fitdistrplus)
library(logspline)
st.er <- function(x) {
  sd(x)/sqrt(length(x))
}
cols<-c('purple', 'green3','gold2', 'orange')

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
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
# Predators
PredSurv <- -ind$Pred
ind <- cbind(ind,PredSurv)
ind$Pred.2 <- ind$Pred+0.08067405

# Merge SEModel with mesocosm results ###################
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

#Get summary statistics for algal production and predator mortality variables ################
  #Beta distribution function found at http://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
      
  pm.mean<-mean(ind1$Pred.2[ind1$Ch==1])
  pm.sd<-sd(ind1$Pred.2[ind1$Ch==1])
  
  descdist(ind1$Pred.2[ind1$Ch==1], discrete = FALSE)#Let's try a beta distribution
  estBetaParams(mu=pm.mean-1, var=pm.sd^2)
      
      hist(ind1$Pred.2[ind1$Ch==1]-1, breaks=10, 
           ylim=c(0,12), xlim = c(0,0.3))
        lines(density(rnorm(5000, pm.mean-1, pm.sd)), col='blue')
        lines(density(rbeta(5000, 11.671, 77.412)), col='red')
        legend('topleft', legend=c('normal', 'beta'), lty=1, 
               col=c('red', 'blue'), title='Distribution')
    #Beta distribution doesn't appear to do too much better than normal
        
  ap0.mean<-mean(ind1$Alg2.2[ind1$At==0 & ind1$Fe==0])
  ap0.sd<-sd(ind1$Alg2.2[ind1$At==0 & ind1$Fe==0])
  
  descdist(ind1$Alg2.2[ind1$At==0 & ind1$Fe==0], discrete = FALSE)#Close to normal but we'll try beta as well
  
    hist(ind1$Alg2.2[ind1$At==0 & ind1$Fe==0], breaks=10, 
         xlim=c(-0.2,0.2), ylim=c(0,10))
      lines(density(rnorm(5000, ap0.mean, ap0.sd)))
  #Normal distribution looks good
      
  ap1.mean<-mean(ind1$Alg2.2[ind1$At==1 & ind1$Fe==0])
  ap1.sd<-sd(ind1$Alg2.2[ind1$At==1 & ind1$Fe==0])
  
  descdist(ind1$Alg2.2[ind1$At==1 & ind1$Fe==0], discrete = FALSE)#Let's try a beta distribution
  estBetaParams(mu=ap1.mean, var=ap1.sd^2)
  
    hist(ind1$Alg2.2[ind1$At==1 & ind1$Fe==0], breaks=10,
         ylim=c(0,12), xlim=c(0,0.4))
      lines(density(rnorm(5000, ap1.mean, ap1.sd)), col='blue')
      lines(density(rbeta(5000, 17.5, 91.61)), col='red')
      legend('topleft', legend=c('normal', 'beta'), lty=1, 
             col=c('red', 'blue'), title='Distribution')
    #Beta and normal appear pretty similar    
      
  ap2.mean<-mean(ind1$Alg2.2[ind1$At==0 & ind1$Fe==1])
  ap2.sd<-sd(ind1$Alg2.2[ind1$At==0 & ind1$Fe==1])
  
  descdist(ind1$Alg2.2[ind1$At==0 & ind1$Fe==1], discrete = FALSE)#Let's try a beta distribution
  estBetaParams(mu=ap2.mean, var=ap2.sd^2)
  
    hist(ind1$Alg2.2[ind1$At==0 & ind1$Fe==1], breaks=10,
         xlim=c(0,0.3), ylim=c(0,20))
      lines(density(rnorm(5000, ap2.mean, ap2.sd)), col='blue')
      lines(density(rbeta(5000, 20.29, 176.55)), col='red')
      legend('topleft', legend=c('normal', 'beta'), lty=1, 
             col=c('red', 'blue'), title='Distribution')
  #Beta and normal appear pretty similar    
  
  ap3.mean<-mean(ind1$Alg2.2[ind1$At==1 & ind1$Fe==1])
  ap3.sd<-sd(ind1$Alg2.2[ind1$At==1 & ind1$Fe==1])
  
  descdist(ind1$Alg2.2[ind1$At==1 & ind1$Fe==1], discrete = FALSE)#Close to uniform/normal
  
    hist(ind1$Alg2.2[ind1$At==1 & ind1$Fe==1], breaks=10,
         xlim=c(0.15, 0.35), ylim=c(0,16))
      lines(density(rnorm(5000, ap3.mean, ap3.sd)), col='blue')
      lines(density(runif(5000,
                          min=min(ind1$Alg2.2[ind1$At==1 & ind1$Fe==1]),
                          max=max(ind1$Alg2.2[ind1$At==1 & ind1$Fe==1]))), col='red')
      legend('topleft', legend=c('normal', 'uniform'), lty=1, 
             col=c('red', 'blue'), title='Distribution')
    #uniform distribution might be better approximation, but 
    #conservative estimate would probably be normal distribution
  
  #Try a 3d plot of bt_fin response to alg production and pred mortality
    cloud(bt_liv_fin ~ Alg2.2+Pred.2, data=ind1)  
      
#Now lets generate the regression equation to be used in predicting snail densities ##########################
 #(exploration of regression equations done in 'Halstead et al SEM snail dens decompose' script)    

  lm.bt.end1<-lm(bt_liv_fin~Pred.2+Alg2.2, data=ind1) #End B. truncatus numbers; prediction function
    summary(lm.bt.end1)
    confint((lm.bt.end1))
    
  
#Now generate some data, informed by mesocosm data, to determine relative effect of algal production ########
 #on final B. truncatus counts
  bt.fin<-data.frame('Pred.2'=rnorm(10000, pm.mean, pm.sd), #Pred mortality for Chlorpyrifos studies to eliminate predation effect on snail numbers
                     'Fe'=as.factor(c(rep(0,5000), rep(1, 5000))),
                     'At'=as.factor(c(rep(0,2500), rep(1, 2500),rep(0,2500), rep(1, 2500))),
                     'Treatment'=c(rep('No_Atra_Fert', 2500),
                                   rep('Atra', 2500),
                                   rep('Fert', 2500),
                                   rep('Atra+Fert', 2500)),
                     'Alg2.2'=c(rnorm(2500, ap0.mean, ap0.sd),
                                rnorm(2500, ap1.mean, ap1.sd),
                                rnorm(2500, ap2.mean, ap2.sd),
                                rnorm(2500, ap3.mean, ap3.sd)),
                     'Bt_liv_fin'=rep(0,10000))    
    bt.fin$Treatment<-factor(bt.fin$Treatment, levels=c('No_Atra_Fert',
                                                     'Atra',
                                                     'Fert',
                                                     'Atra+Fert'))
    
  #Generate resulting algal production density plots  
    pred2.dens<-density(bt.fin$Pred.2)
    alg0.dens<-density(bt.fin$Alg2.2[bt.fin$Fe==0 & bt.fin$At==0])
    alg1.dens<-density(bt.fin$Alg2.2[bt.fin$Fe==1 & bt.fin$At==0])
    alg2.dens<-density(bt.fin$Alg2.2[bt.fin$Fe==0 & bt.fin$At==1])
    alg3.dens<-density(bt.fin$Alg2.2[bt.fin$Fe==1 & bt.fin$At==1])
    
  #Plot pred mortality distribution  
  ggplot(bt.fin, aes(x=Pred.2))+
    theme_bw()+
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15))+
    xlim(1.0, 1.3) +
    ylim(0,12)+
    xlab('Predator mortality')+
    geom_density(aes(fill=Treatment, border='black'), alpha=0.5)+
    scale_fill_manual(values=rep('red',4))
  
  #Plot algal distributions  
  ggplot(bt.fin, aes(x=Alg2.2)) + 
    theme_bw()+
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15))+
    scale_fill_manual(values=c('purple', 'green2', 'gold2', 'orange'))  +
    geom_density(aes(group=Treatment, colour=Treatment, fill=Treatment, 
                     border='black'), alpha=0.5)
    
  #Use regression parameters to predict snail densities
    bt.fin.sim<-predict(lm.bt.end1, bt.fin, se.fit=T)

    bt.fin$Bt_liv_fin<-bt.fin.sim$fit
    
    bt.fin$Bt_liv_fin.se<-bt.fin.sim$se.fit
    
    bt.fin$Treatment<-paste(bt.fin$At, '_', bt.fin$Fe)
  
  ggplot(bt.fin, aes(x=Alg2.2, y=Bt_liv_fin, colour=Treatment)) +
    theme_bw()+
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15))+
    scale_color_manual(values=cols) +
    ylab("Mean +/- SEM predicted B. truncatus") +
    geom_point(size=4)+
    geom_errorbar(aes(ymin=Bt_liv_fin-Bt_liv_fin.se,
                      ymax=Bt_liv_fin+Bt_liv_fin.se), width=.01)
  
  ggplot(bt.fin, aes(x=Bt_liv_fin))+ 
    theme_bw()+
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15))+
    scale_fill_manual(values=c('purple', 'green2', 'gold2', 'orange'))  +
    geom_density(aes(group=Treatment, colour=Treatment, fill=Treatment, 
                     border='black'), alpha=0.5)
  
#Get relative increases in carrying capacity for each treatment########
  phi.n.ref<-rep(0,10000)
  for(i in 1:length(phi.n.ref)){
    phi.n.ref[i] = 
      (predict(lm.bt.end1, newdata=data.frame
               ('Alg2.2'=rnorm(1, ap0.mean, ap0.sd), 
               'Pred.2'=rnorm(1, mean=pm.mean, sd=pm.sd)))) / 
      (predict(lm.bt.end1, newdata=data.frame
               ('Alg2.2'=rnorm(1, ap0.mean, ap0.sd),
               'Pred.2'=rnorm(1, mean=pm.mean, sd=pm.sd))))
  }
  mean(phi.n.ref)
  sd(phi.n.ref)
  
  phi.n.atra<-rep(0,10000)
    for(i in 1:length(phi.n.atra)){
    phi.n.atra[i] = 
      (predict(lm.bt.end1, newdata=data.frame
              ('Alg2.2'=rnorm(1, ap1.mean, ap1.sd), 
              'Pred.2'=rnorm(1, mean=pm.mean, sd=pm.sd)))) / 
      (predict(lm.bt.end1, newdata=data.frame
              ('Alg2.2'=rnorm(1, ap0.mean, ap0.sd),
              'Pred.2'=rnorm(1, mean=pm.mean, sd=pm.sd))))
    } #Predict proportion increase in carrying capacity for atrazine treatments
        #by randomly sampling from algal distributions and predator mortality
        #distribution
  phi.n.fert<-rep(0,10000)
    for(i in 1:length(phi.n.fert)){
      phi.n.fert[i] = 
        (predict(lm.bt.end1, newdata=data.frame
                 ('Alg2.2'=rnorm(1, ap2.mean, ap2.sd), 
                 'Pred.2'=rnorm(1, mean=pm.mean, sd=pm.sd)))) / 
        (predict(lm.bt.end1, newdata=data.frame
                 ('Alg2.2'=rnorm(1, ap0.mean, ap0.sd),
                 'Pred.2'=rnorm(1, mean=pm.mean, sd=pm.sd))))
    } #Predict proportion increase in carrying capacity for fertilizer treatments
    #by randomly sampling from algal distributions and predator mortality
    #distribution
    
  phi.n.atfe<-rep(0,10000)
    for(i in 1:length(phi.n.atfe)){
      phi.n.atfe[i] = 
        (predict(lm.bt.end1, newdata=data.frame
                 ('Alg2.2'=rnorm(1, ap3.mean, ap3.sd), 
                 'Pred.2'=rnorm(1, mean=pm.mean, sd=pm.sd)))) / 
        (predict(lm.bt.end1, newdata=data.frame
                 ('Alg2.2'=rnorm(1, ap0.mean, ap0.sd),
                 'Pred.2'=rnorm(1, mean=pm.mean, sd=pm.sd))))
    } #Predict proportion increase in carrying capacity for fertilizer treatments
    #by randomly sampling from algal distributions and predator mortality
    #distribution
    
  phi.n.plot<-data.frame('mean'=c(mean(phi.n.fert),
                                  mean(phi.n.atra),
                                  mean(phi.n.atfe)),
                         'st.d'=c(sd(phi.n.fert),
                                  sd(phi.n.atra),
                                  sd(phi.n.atfe)),
                         'Treatment'=c('Fertilizer',
                                              'Atrazine',
                                              'Both'))
  phi.n.plot$Treatment<-factor(phi.n.plot$Treatment, levels=c('Fertilizer',
                                                              'Atrazine',
                                                              'Both'))
  
  ggplot(phi.n.plot, aes(x=Treatment, y=mean, fill=Treatment))+
    theme_bw()+
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15))+
    scale_fill_manual(values=c('green2', 'gold2', 'orange')) +
    ylab("Scalar of Carrying capacity")+
    geom_bar(position=position_dodge(), stat="identity", width = .7) +
    geom_errorbar(aes(ymin=mean-st.d,
                      ymax=mean+st.d),
                  width=.2, position=position_dodge(.7))
    
    
#Aggregate and merge for bar plot ###################     
  bt.fin.agg1<-aggregate(bt.fin, by=list(bt.fin[,7]), FUN = mean)

  bt.fin.agg2<-aggregate(bt.fin, by=list(bt.fin[,7]), FUN = st.er)
  
  bt.fin.agg<-merge(bt.fin.agg1, bt.fin.agg2, by='Group.1')
  
  bt.fin.agg<-bt.fin.agg[,c(1,6,7,13)]
  
  colnames(bt.fin.agg)<-c("Treatment", "Mean", "st.err", "st.err2")#st.err is mean of model standard errors in each treatment
    bt.fin.agg$Treatment<-as.character(bt.fin.agg$Treatment)
      
      bt.fin.agg$Treatment[bt.fin.agg$Treatment=='0 _ 0']<-"Pred-free base"
      bt.fin.agg$Treatment[bt.fin.agg$Treatment=='1 _ 0']<-"Pred-free atra"
      bt.fin.agg$Treatment[bt.fin.agg$Treatment=='0 _ 1']<-"Pred-free fert"
      bt.fin.agg$Treatment[bt.fin.agg$Treatment=='1 _ 1']<-"Pred-free atra+fert"
      
    bt.fin.agg$Treatment<-factor(bt.fin.agg$Treatment, levels=c("Pred-free base", "Pred-free fert",
                                                         "Pred-free atra","Pred-free atra+fert")) 

  ggplot(bt.fin.agg, aes(x=Treatment, y=Mean, fill=Treatment)) +
    theme_bw()+
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15))+
    scale_fill_manual(values=cols) +
    ylab("Mean +/- SEM predicted B. truncatus")+
    geom_bar(position=position_dodge(), stat="identity", width = .7) +
    geom_errorbar(aes(ymin=Mean-st.err,
                      ymax=Mean+st.err),
                  width=.2, position=position_dodge(.7))
#Compare to observed data and estimate parameters########################
  bt.obs<-data.frame('Treatment'=
                       c('ChlorP','ChlorP+Fert',
                                   'ChlorP+Atra', 'ChlorP+Atra+Fert'),
                     'Mean_Bt_fin'=
                       c(mean(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 &dat$fert==0]),
                         mean(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 &dat$fert==1]),
                         mean(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 &dat$fert==0]),
                         mean(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 &dat$fert==1])),
                     'St.err_Bt_fin'=
                       c(st.er(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 &dat$fert==0]),
                         st.er(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 &dat$fert==1]),
                         st.er(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 &dat$fert==0]),
                         st.er(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 &dat$fert==1])))
  
  bt.obs$Treatment<-factor(bt.obs$Treatment, levels = c('ChlorP',
                                                        'ChlorP+Fert',
                                                        'ChlorP+Atra',
                                                        'ChlorP+Atra+Fert'))
  
#Are there significant differences between treatments when ChlorP is present?
  alg.sub<-subset(ind1, Ch==1)
    alg.sub$Treats<-factor(alg.sub$Treats, levels=c("ChlorP",
                                                    "Chlor_Fert",
                                                    "Atra_Chlor",
                                                    "All_Three"))
  
  bt.anova<-aov(bt_liv_fin ~ Treats, data=alg.sub)
    summary(bt.anova)
    plot(bt_liv_fin ~ Treats, data=alg.sub)
    
#Plot observed bt_fin counts between treatments    
  ggplot(bt.obs, aes(x=Treatment, y=Mean_Bt_fin, fill=Treatment))+
    theme_bw()+
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15))+
    scale_fill_manual(values=cols) +
    ylab("Mean +/- SEM observed B. truncatus")+
    geom_bar(position=position_dodge(), stat="identity", width = .7) +
    geom_errorbar(aes(ymin=Mean_Bt_fin-St.err_Bt_fin,
                      ymax=Mean_Bt_fin+St.err_Bt_fin),
                  width=.2, position=position_dodge(.7))
  
#Chlorpyrifos only tanks  
  obs1.mean<-mean(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==0])
  obs1.sd<-sd(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==0])
    
  #Check distribution
    hist(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==0], 
           breaks=10, xlim=c(0,550), ylim=c(0,5))
    
    descdist(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==0])
    #Let's try a Weibull? Gamma? distribution
    
  #Log transformed distribution
    hist(log(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==0]), 
         breaks=5, ylim=c(0,5))
        
    descdist(log(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==0]))
    #Let's try a Weibull? Gamma? distribution
      
      obs1.weib<-fitdistr(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==0],
                        'weibull')
    
#Chlorpyrifos + Fertilizer tanks
  obs2.mean<-mean(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==1])
  obs2.sd<-sd(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==1])
    
  #Plot distribution
    hist(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==1], 
           breaks=10, xlim=c(0,550), ylim=c(0,5))
    
    descdist(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==1])
    #Let's try an exponential distribution
    
  #Plot log transformed    
    hist(log(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==1]), 
         breaks=10, ylim=c(0,5))
    
    descdist(log(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==1]))
    #Let's try an exponential distribution
    
#Chlorpyrifos + atrazine tanks  
  obs3.mean<-mean(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==0])
  obs3.sd<-sd(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==0])
  
  #Plot histogram
    hist(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==0], 
           breaks=10, xlim=c(0,550), ylim=c(0,5))
      
    descdist(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==0])
    #Let's try an exponential distribution
  
  #Log transform
    hist(log(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==0]), 
         breaks=10, ylim=c(0,5))
    
    descdist(log(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==0]))
    #Let's try an exponential distribution
    
#Chlorpyrifos +Atrazine +Fertilizer treatments      
  obs4.mean<-mean(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==1])
  obs4.sd<-sd(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==1])
  
  #Plot histogram
    hist(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==1], 
           breaks=10, xlim=c(0,550), ylim=c(0,5))
      
    descdist(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==1])
    #Let's try an exponential distribution
    
  #Plot log-transformed histogram
    hist(log(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==1]), 
         breaks=10, ylim=c(0,5))
    
    descdist(log(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==1]))
    #Let's try an exponential distribution
    
#3/11/2016: Decided to estimate scalar of carrying capacity from mean values in observed data###################
  bt.obs$Scalar<-c((obs1.mean/obs1.mean),
                   (obs2.mean/obs1.mean),
                   (obs3.mean/obs1.mean),
                   (obs4.mean/obs1.mean))
  
    blerb<-predict(lm.bt.end1,data.frame('Pred.2'=c(rep(pm.mean, 4)),
                                       'Alg2.2'=c(ap0.mean, ap1.mean, ap2.mean, ap3.mean)),
                 se.fit=T) #Estimates with st.err from regression
    
  bt.obs$Scalar_max<-c(((obs1.mean+blerb$se.fit[[1]])/obs1.mean),
                       ((obs2.mean+blerb$se.fit[[2]])/obs1.mean),
                       ((obs3.mean+blerb$se.fit[[3]])/obs1.mean),
                       ((obs4.mean+blerb$se.fit[[4]])/obs1.mean)) 
  
  
  bt.obs$Scalar_min<-c(((obs1.mean-blerb$se.fit[[1]])/obs1.mean),
                       ((obs2.mean-blerb$se.fit[[2]])/obs1.mean),
                       ((obs3.mean-blerb$se.fit[[3]])/obs1.mean),
                       ((obs4.mean-blerb$se.fit[[4]])/obs1.mean)) 
  
  ggplot(bt.obs, aes(x=Treatment, y=Scalar, fill=Treatment))+
    theme_bw()+
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15))+
    scale_fill_manual(values=cols) +
    ylab("Mean +/- SEM observed B. truncatus")+
    geom_bar(position=position_dodge(), stat="identity", width = .7) +
    geom_errorbar(aes(ymin=Scalar_min,
                      ymax=Scalar_max),
                  width=.2, position=position_dodge(.7))
  
  
#Estimate scalar of carrying capacity based on distribution of observed data ###############
    phi.n.obs1<-rep(0,100000) #Fertilizer treatments
    for(i in 1:length(phi.n.obs1)){
      phi.n.obs1[i] = 
        (rnorm(1, obs2.mean, obs2.sd) #Random sample from observed Bt_fin in fertilizer+chlorP treatments
          /
        rnorm(1, obs1.mean, obs1.sd))#Random sample from observed Bt_fin in chlorP only treatments
      }
    mean(phi.n.obs1)
    sd(phi.n.obs1)
    
    phi.n.obs2<-rep(0,100000) #Atrazine treatments
    for(i in 1:length(phi.n.obs2)){
      phi.n.obs2[i] = 
        (rnorm(1, obs3.mean, obs3.sd) #Random sample from observed Bt_fin in atrazine+chlorP treatments
      /
        rnorm(1, obs1.mean, obs1.sd)) #Random sample from observed Bt_fin in chlorP only treatments
    }
    mean(phi.n.obs2)
    sd(phi.n.obs2)
    
    phi.n.obs3<-rep(0,100000) #Atrazine+Fertilizer treatments
    for(i in 1:length(phi.n.obs3)){
      phi.n.obs3[i] = 
        (rnorm(1, obs4.mean, obs4.sd) #Random sample from observed Bt_fin in atrazine+chlorP treatments
      /
        rnorm(1, obs1.mean, obs1.sd)) #Random sample from observed Bt_fin in chlorP only treatments
    }
    mean(phi.n.obs3)
    sd(phi.n.obs3)
    
    phi.obs.plot<-data.frame('mean'=c(mean(phi.n.obs1),
                                    mean(phi.n.obs2),
                                    mean(phi.n.obs3)),
                           'st.d'=c(sd(phi.n.obs1),
                                    sd(phi.n.obs2),
                                    sd(phi.n.obs3)),
                           'Treatment'=c('Fertilizer',
                                         'Atrazine',
                                         'Both'))
    phi.obs.plot$Treatment<-factor(phi.obs.plot$Treatment, levels=c('Fertilizer',
                                                                'Atrazine',
                                                                'Both'))
    
    ggplot(phi.obs.plot, aes(x=Treatment, y=mean, fill=Treatment))+
      theme_bw()+
      theme(axis.title=element_text(size=20),
            axis.text=element_text(size=15))+
      scale_fill_manual(values=c('green2', 'gold2', 'orange')) +
      ylab("Scalar of Carrying capacity")+
      geom_bar(position=position_dodge(), stat="identity", width = .7) +
      geom_errorbar(aes(ymin=mean-st.d,
                        ymax=mean+st.d),
                    width=.2, position=position_dodge(.7))  
    
#Compare to regression prediction as opposed to generated data prediction ##########
  bt.reg<-data.frame('Tank'=ind1$Tank,
                     'Atra'=ind1$At,
                     'Chlor'=ind1$Ch,
                     'Fert'=ind1$Fe,
                     'Treatment'=ind1$Treats,
                     'Alg2.2'=ind1$Alg2.2,
                     'Pred.2'=ind1$Pred.2,
                     'Bt_liv_fin_obs'=ind1$bt_liv_fin)
  bt.reg.fit<-predict(lm.bt.end1, bt.reg, se.fit=T)
  bt.reg$fit<-bt.reg.fit$fit
  bt.reg$fit.se<-bt.reg.fit$se.fit
  
  bt.reg.plot<-subset(bt.reg, Treatment=='All_Three' | Treatment=='ChlorP' |
                        Treatment=='Chlor_Fert' | Treatment=='Atra_Chlor')
  
  bt.reg.plot1<-aggregate(bt.reg.plot, by=list(bt.reg.plot$Treatment), FUN=mean)
  
  bt.reg.plot2<-aggregate(bt.reg.plot, by=list(bt.reg.plot$Treatment), FUN=st.er)
  bt.reg.plot<-merge(bt.reg.plot1, bt.reg.plot2, by='Group.1')
  bt.reg.plot$Group.1<-as.factor(c('Pred-free Atra+fert', 'Pred-free atra',
                                 'Pred-free fert', 'pred-free only'))
  bt.reg.plot$Group.1<-factor(bt.reg.plot$Group.1, levels=c('pred-free only',
                                                            'Pred-free fert',
                                                            'Pred-free atra',
                                                            'Pred-free Atra+fert'))
  colnames(bt.reg.plot)[1]<-'Treatment'
  
ggplot(bt.reg.plot, aes(x=Treatment, y=fit.x, fill=Treatment))+
  theme_bw()+
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=15))+
  scale_fill_manual(values=cols) +
  ylab("Mean +/- SEM B. truncatus")+
  geom_bar(position=position_dodge(), stat="identity", width = .7) +
  geom_errorbar(aes(ymin=fit.x-fit.se.x,
                    ymax=fit.x+fit.se.x),
                width=.2, position=position_dodge(.7))