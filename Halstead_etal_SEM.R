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
source("lavaan.modavg.R")
library(semPlot)
library(ggplot2)
st.er <- function(x) {
  sd(x)/sqrt(length(x))
}


setwd('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal')
ind <- read.csv("SEM_data.csv") # individually-transformed variables

#Create latent variables for algae from agal data
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
#write.csv(ind, file="indSchistoFinal2010.csv")

# Full structural equation model
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

#Chris' addition ##################
  plot(x=ind$Algae, y=ind$Snails)
    points(x=ind$Algae[ind$Ch==1], y=ind$Snails[ind$Ch==1], col='red', pch=15, cex=2)
    points(x=ind$Algae[ind$At==1], y=ind$Snails[ind$At==1], col='yellow', pch=17, cex=1.5)
    points(x=ind$Algae[ind$Fe==1], y=ind$Snails[ind$Fe==1], col='green', pch=16)
#End Chris' addition ####################

# Plots of linear relationships between predictor and latent variables
# Predators
PredSurv <- -ind$Pred
ind <- cbind(ind,PredSurv)
Pred.2 <- ind$Pred+0.08067405

#Chris' addition ##################
ind$Pred.2<-ind$Pred+0.08067405
#End Chris' addition ####################

# calculate the mean of each group using tapply, which returns a matrix, and save it as an object
mean.preds <- tapply(ind$Pred+0.08067405, list(ind$Ch), mean)
# calculate the standard deviation of each group, and save it as an object
sd.preds <- tapply(ind$Pred+0.08067405, list(ind$Ch), sd)
# calculate the sample size in each group, and save it as an object
n.preds <- tapply(ind$Pred, list(ind$Ch), length)
sem.preds <- sd.preds/(n.preds^0.5)

mids <- barplot(mean.preds,
                xlab="Chlorpyrifos",
                ylab="Predator mortality",
                ylim=c(-0.1, 1.2),
                col="white", cex.axis=1.5,
                names.arg=c("Absent","Present"),
                cex.names=1.55, cex.lab=1.55)
arrows(mids, mean.preds-sem.preds, mids, mean.preds+sem.preds, code=3, angle=90, length=0.1)
text(mids, 0.04, paste("n=", n.preds), cex=1.2)

# Snails as a function of predators
snpr.fit<- lm(ind$Snail2.2~Pred.2)
summary(snpr.fit) # R^2=0.9374
plot(Pred.2, ind$Snail2.2, bty="l", xlab="Predator mortality", ylab="Snail density",
     pch=19, bg="black", cex.axis=1.5, cex.lab=1.55)
lines(sort(Pred.2), fitted(snpr.fit)[order(Pred.2)])

#Chris' addition ##################
plot(ind$Pred.2, ind$Snail2.2, bty="l", xlab="Predator mortality", ylab="Snail density",
     pch=19, bg="black", cex.axis=1.5, cex.lab=1.55)
  points(x=ind$Pred.2[ind$Ch==1], y=ind$Snail2.2[ind$Ch==1], col='red', pch=15, cex=2)
  points(x=ind$Pred.2[ind$At==1], y=ind$Snail2.2[ind$At==1], col='yellow', pch=17, cex=1.5)
  points(x=ind$Pred.2[ind$Fe==1], y=ind$Snail2.2[ind$Fe==1], col='green', pch=16)
#End Chris' addition ####################

# Snails as a function of algae
# Use residuals of relationship with predators
Sn2.2resid <- resid(lm(ind$Snail2.2~ind$PredSurv))

#Chris' addition ##################
  plot(x=ind$PredSurv, y=ind$Snail2.2, pch=19)
ind$Sn2.2resid<-resid(lm(ind$Snail2.2~ind$PredSurv))
#End Chris' addition ####################
  
resid.fit <- lm(Sn2.2resid~ind$Alg2.2+I(ind$Alg2.2^2))
summary(resid.fit) # R^2=0.467
plot(ind$Alg2.2, ind$Sn2.2resid, pch=19, bty="l", xlab="Algal production",
     ylab="Residual snail density", bg="black", cex.axis=1.5, cex.lab=1.55)
lines(sort(ind$Alg2.2), fitted(resid.fit)[order(ind$Alg2.2)])

#Chris' addition ##################
plot(ind$Alg2.2, ind$Sn2.2resid, pch=19, bty="l", xlab="Algal production",
     ylab="Residual snail density", bg="black", cex.axis=1.5, cex.lab=1.55)
  points(x=ind$Alg2.2[ind$Ch==1], y=ind$Sn2.2resid[ind$Ch==1], col='red', pch=15, cex=2)
  points(x=ind$Alg2.2[ind$At==1], y=ind$Sn2.2resid[ind$At==1], col='yellow', pch=17, cex=1.5)
  points(x=ind$Alg2.2[ind$Fe==1], y=ind$Sn2.2resid[ind$Fe==1], col='green', pch=16)
    #Line to overall fit
      lines(sort(ind$Alg2.2), fitted(resid.fit)[order(ind$Alg2.2)], lwd=2)
    #Line to just chlorpyrifos fit  
      resid.fit.Ch <- lm(ind$Sn2.2resid[ind$Ch==1]~ind$Alg2.2[ind$Ch==1]+I(ind$Alg2.2[ind$Ch==1]^2))
      lines(sort(ind$Alg2.2[ind$Ch==1]), fitted(resid.fit.Ch)[order(ind$Alg2.2[ind$Ch==1])], col='red', lwd=2)
    #Line to just atrazine fit
      resid.fit.At <- lm(ind$Sn2.2resid[ind$At==1]~ind$Alg2.2[ind$At==1]+I(ind$Alg2.2[ind$At==1]^2))
      lines(sort(ind$Alg2.2[ind$At==1]), fitted(resid.fit.At)[order(ind$Alg2.2[ind$At==1])], col='yellow', lwd=2)
    #Line to just Fertilizer fit
      resid.fit.Fe <- lm(ind$Sn2.2resid[ind$Fe==1]~ind$Alg2.2[ind$Fe==1]+I(ind$Alg2.2[ind$Fe==1]^2))
      lines(sort(ind$Alg2.2[ind$Fe==1]), fitted(resid.fit.Fe)[order(ind$Alg2.2[ind$Fe==1])], col='green', lwd=2)
    legend("topleft", legend=c('None (Control)', 'Atrazine','Chlorpyrifos', 'Fertilizer'),
           pch=c(20,17,15,16), col=c('black', 'yellow','red','green'), title="AgroCs Present")

#Add unique treatment variable    
for(i in 1:length(ind$Tank)){
  ind[i,8]<-paste(ind[i,2],'_',ind[i,3],'_',ind[i,4])
}
    colnames(ind)[8]<-'Treats'
    
    ind$Treats<-as.character(ind$Treats)
    
    ind$Treats[ind$Treats=="0 _ 0 _ 0"]<- "Control"
    ind$Treats[ind$Treats=="1 _ 0 _ 0"]<- "Atra"
    ind$Treats[ind$Treats=="0 _ 1 _ 0"]<- "ChlorP"
    ind$Treats[ind$Treats=="0 _ 0 _ 1"]<- "Fert"
    ind$Treats[ind$Treats=="1 _ 1 _ 0"]<- "Atra_Chlor"
    ind$Treats[ind$Treats=="1 _ 0 _ 1"]<- "Atra_Fert"
    ind$Treats[ind$Treats=="0 _ 1 _ 1"]<- "Chlor_Fert"
    ind$Treats[ind$Treats=="1 _ 1 _ 1"]<- "All_Three"
    
    ind$Treats<-as.factor(ind$Treats)
    
    ind$Treats<- factor(ind$Treats, levels=c("Control","Atra","ChlorP",
                                                                         "Fert","Atra_Chlor",
                                                                         "Atra_Fert",
                                                                         "Chlor_Fert","All_Three"))
    
    
#Get rid of negative residuals to calculate means and standard errors    
  ind$Sn2.2resid<-ind$Sn2.2resid-min(ind$Sn2.2resid)
#Anova to see if there are significant differences in resid snail density between treatments   
  anv<-aov(Sn2.2resid~Treats, data=ind)
    summary(anv)
  boxplot(Sn2.2resid~Treats, data=ind, ylab="Residual Snail Density", xlab="Treatments", cex.lab=1.5, cex.axis=1.2)
#calc means/st. dev for treatment presence/absence
  atra.1<-mean(ind$Sn2.2resid[ind$At==1])
    atra.1.se<-st.er(ind$Sn2.2resid[ind$At==1])
  atra.0<-mean(ind$Sn2.2resid[ind$At==0])
    atra.0.se<-st.er(ind$Sn2.2resid[ind$At==0])
  fert.1<-mean(ind$Sn2.2resid[ind$Fe==1])
    fert.1.se<-st.er(ind$Sn2.2resid[ind$Fe==1])
  fert.0<-mean(ind$Sn2.2resid[ind$Fe==0])
    fert.0.se<-st.er(ind$Sn2.2resid[ind$Fe==0])
  atra.fert.1<-mean(ind$Sn2.2resid[ind$Fe==1 & ind$At==1])
    atra.fert.1.se<-st.er(ind$Sn2.2resid[ind$Fe==1 & ind$At==1])
    
t.test(y=ind$Sn2.2resid[ind$At==0], x=ind$Sn2.2resid[ind$At==1], alternative = 'greater')    
t.test(y=ind$Sn2.2resid[ind$Fe==0], x=ind$Sn2.2resid[ind$Fe==1], alternative = 'greater')   
t.test(y=ind$Sn2.2resid[ind$At==0], x=ind$Sn2.2resid[ind$Fe==1])
t.test(y=ind$Sn2.2resid[ind$At==0], x=ind$Sn2.2resid[ind$Fe==0])    
t.test(y=ind$Sn2.2resid[ind$Fe==1 & ind$At==1], x=ind$Sn2.2resid[ind$At==1], alternative = 'greater')    

    
bar.dat<-data.frame('Treatment'=c('Fert_Absent', 'Atra_Absent', 'Fert_Present', 'Atra_Present', 'Both_Present'),
                      'mean'=c(atra.0,fert.0,atra.1,fert.1,atra.fert.1),
                      'st.err'=c(atra.0.se,fert.0.se,atra.1.se,fert.1.se, atra.fert.1.se)) 
  bar.dat$Treatment<-factor(bar.dat$Treatment, levels=c('Fert_Absent', 'Atra_Absent', 'Fert_Present', 
                                                        'Atra_Present', 'Both_Present'))

  bar.cols<-c('gray75', 'gray75', 'green2', 'gold2',  'darkolivegreen')
  
ggplot(bar.dat, aes(x=Treatment, y=mean, fill=Treatment)) +
  theme_bw()+
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text=element_text(size=15)) +
  ylab('Residual Snail Density Means +/- St.err') +
  scale_fill_manual(values=bar.cols) +
  geom_bar(position=position_dodge(), stat="identity", width=.7) +
  geom_errorbar(aes(ymin=mean-st.err,
                    ymax=mean+st.err),
                width=.2, position=position_dodge(.7))

    
#End Chris' addition ####################

# Snails as a function of algae
plot(ind$Alg2.2, ind$Snail2.2, xlab="Algae", ylab="Snail Density", pch=19, bg="black")