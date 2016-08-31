#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

###Load Libraries & data
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
snail<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/snail_repro_totals.csv")
  snail$treats=paste(snail$atra, snail$chlor, snail$fert, sep='_')
  snail$treat<-as.character(snail$treat)
    snail$treat[snail$treat=="S"]<-"Control"
    snail$treat[snail$treat=="W"]<-"Control"
  snail$treat<-factor(snail$treat, levels = c("A 1x", "A 2x", "F 1x", "F 2x", "C 1x", "C 2x",
                                              "C + A", "A + F", "F + C", "A + C + F", "Control"))  
  snail.treats<-unique(snail$treat)
  
cbPallete=c("yellow", "gold3", "green", "darkgreen", "red", "darkred",
            "orange", "purple", "blue", "black", "grey")
#Reshape to long format and make subset data frames#########################
snail2<-reshape(snail, varying = colnames(snail[-c(1:6)]),
                v.names="measure",
                timevar="variable",
                times=colnames(snail[-c(1:6)]),
                new.row.names=1:4020,
                direction="long")
snail2<-snail2[,-c(2:4,9)] #Get rid of needless variables
  
eggs<-subset(snail2, variable=='tot_egg1' | 
                     variable=='tot_egg2' | 
                     variable=='tot_egg3' | 
                     variable=='tot_egg4' |
                     variable=='tot_egg8')  
  #Replace variable names with time descriptors
    eggs$variable[eggs$variable=='tot_egg1']<-'Week1'
    eggs$variable[eggs$variable=='tot_egg2']<-'Week2'
    eggs$variable[eggs$variable=='tot_egg3']<-'Week3'
    eggs$variable[eggs$variable=='tot_egg4']<-'Week4'
    eggs$variable[eggs$variable=='tot_egg8']<-'Week8'
  eggs$tank<-as.factor(eggs$tank)  
    
bt.eggs<-subset(snail2, variable=='TBtegg1' | 
                        variable=='TBtegg2' | 
                        variable=='TBtegg3' | 
                        variable=='TBtegg4' |
                        variable=='TBtegg8')  
  #Replace variable names with time descriptors
    bt.eggs$variable[bt.eggs$variable=='TBtegg1']<-'Week1'
    bt.eggs$variable[bt.eggs$variable=='TBtegg2']<-'Week2'
    bt.eggs$variable[bt.eggs$variable=='TBtegg3']<-'Week3'
    bt.eggs$variable[bt.eggs$variable=='TBtegg4']<-'Week4'
    bt.eggs$variable[bt.eggs$variable=='TBtegg8']<-'Week8'    
    

hatchs<-subset(snail2, variable=='tot_hatch1' | 
                       variable=='tot_hatch2' | 
                       variable=='tot_hatch3' | 
                       variable=='tot_hatch4' |
                       variable=='tot_hatch8')  
  #Replace variable names with time descriptors
    hatchs$variable[hatchs$variable=='tot_hatch1']<-'Week1'
    hatchs$variable[hatchs$variable=='tot_hatch2']<-'Week2'
    hatchs$variable[hatchs$variable=='tot_hatch3']<-'Week3'
    hatchs$variable[hatchs$variable=='tot_hatch4']<-'Week4'
    hatchs$variable[hatchs$variable=='tot_hatch8']<-'Week8'

adults<-subset(snail2, variable=='tot_adult1' | 
                       variable=='tot_adult2' | 
                       variable=='tot_adult3' | 
                       variable=='tot_adult4' |
                       variable=='tot_adult8')  
  #Replace variable names with time descriptors
    adults$variable[adults$variable=='tot_adult1']<-'Week1'
    adults$variable[adults$variable=='tot_adult2']<-'Week2'
    adults$variable[adults$variable=='tot_adult3']<-'Week3'
    adults$variable[adults$variable=='tot_adult4']<-'Week4'
    adults$variable[adults$variable=='tot_adult8']<-'Week8'

#Aggregate data frames ###############      
 #Egg production over time
  eggs.ag1<-aggregate.data.frame(eggs, by=list(eggs[,2], eggs[,4]), FUN=mean) 
    eggs.ag1<-eggs.ag1[,-c(3:6)] #remove unneeded variables 
      colnames(eggs.ag1)<-c("Treatment", "variable", "mean")
  
  eggs.ag2<-aggregate.data.frame(eggs, by=list(eggs[,2], eggs[,4]), FUN=st.er) #calculate st.error of treatment groups
    eggs.ag2<-eggs.ag2[,-c(3:6)] #remove unneeded variables 
      colnames(eggs.ag2)<-c("Treatment", "variable", "st.err")  
  
  eggs.ag<-merge(eggs.ag1, eggs.ag2, by=c("Treatment", "variable"))
  
# B. truncatus egg production over time
  bteggs.ag1<-aggregate.data.frame(bt.eggs, by=list(bt.eggs[,2], bt.eggs[,4]), FUN=mean) 
    bteggs.ag1<-bteggs.ag1[,-c(3:6)] #remove unneeded variables 
      colnames(bteggs.ag1)<-c("Treatment", "variable", "mean")
  
  bteggs.ag2<-aggregate.data.frame(bt.eggs, by=list(bt.eggs[,2], bt.eggs[,4]), FUN=st.er) 
    bteggs.ag2<-bteggs.ag2[,-c(3:6)] #remove unneeded variables 
      colnames(bteggs.ag2)<-c("Treatment", "variable", "st.err")
  
  bteggs.ag<-merge(bteggs.ag1, bteggs.ag2, by=c("Treatment", "variable"))
#Hatchling production over time  
  hatchs.ag1<-aggregate.data.frame(hatchs, by=list(hatchs[,2], hatchs[,4]), FUN=mean) 
    hatchs.ag1<-hatchs.ag1[,-c(3:6)] #remove unneeded variables 
      colnames(hatchs.ag1)<-c("Treatment", "variable", "mean")
  
  hatchs.ag2<-aggregate.data.frame(hatchs, by=list(hatchs[,2], hatchs[,4]), FUN=st.er) #calculate st.error of treatment groups
    hatchs.ag2<-hatchs.ag2[,-c(3:6)] #remove unneeded variables 
      colnames(hatchs.ag2)<-c("Treatment", "variable", "st.err")  
  
    hatchs.ag<-merge(hatchs.ag1, hatchs.ag2, by=c("Treatment", "variable"))

  
#Plot to examine differences between treatment groups ###################
 #Eggs over time    
  ggplot(eggs.ag, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
    scale_color_manual(values=cbPallete) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.25))
    
 #Eggs over time in individual tanks
  ggplot(eggs, aes(x=variable, y=measure, group=tank, color=treat)) +
    theme_bw()+
    scale_color_manual(values=cbPallete) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5)
  
 #Eggs over time in atrazine containing and control tanks
  eggs.at<-subset(eggs, treat=="A 1x" | 
                        treat=="A 2x" |
                        treat=="C + A" |
                        treat=="A + F" |
                        treat=="A + C + F" |
                        treat=="Control")
  
  ggplot(eggs.at, aes(x=variable, y=measure, group=tank, color=treat)) +
    theme_bw()+
    scale_color_manual(values=c("yellow", "gold2", "orange", "purple", "black", "grey")) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5)
  
 #Eggs over time in chlorP and control tanks
  eggs.ch<-subset(eggs, treat=="C 1x" | 
                    treat=="C 2x" |
                    treat=="C + A" |
                    treat=="F + C" |
                    treat=="A + C + F" |
                    treat=="Control")
  
  ggplot(eggs.ch, aes(x=variable, y=measure, group=tank, color=treat)) +
    theme_bw()+
    scale_color_manual(values=c("red", "darkred", "orange", "blue", "black", "grey")) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5)
  
  #Just control and two atrazine doses
    eggs.ag.at<-subset(eggs.ag, Treatment=="A 1x" | 
                                Treatment=="A 2x" |
                                Treatment=="Control")
    ggplot(eggs.ag.at, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
      theme_bw()+
      ylab("Mean +/- st.error sampled egg masses")+
      xlab("Time")+
      scale_color_manual(values=c("yellow", #At 1x
                                  "gold3", #At 2x
                                  "grey" #Controls
                                  )) +
      geom_line(position=position_dodge(.25), size=1) +
      geom_point(position=position_dodge(.25), size=3.5) +
      geom_errorbar(aes(ymin=mean-st.err,
                        ymax=mean+st.err),
                    width=.2, position=position_dodge(.25))
    
  #Just chlorP-free treatment groups
    eggs.ag.ch<-subset(eggs.ag, Treatment=="A 1x" |
                                Treatment=="A 2x" |
                                Treatment=="F 1x" |
                                Treatment=="F 2x" |
                                Treatment=="A + F" |
                                Treatment=="Control")
                         
    ggplot(eggs.ag.ch, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
      theme_bw()+
      scale_color_manual(values=c("yellow", #At 1x
                                  "gold3", #At 2x
                                  "green", #Fe 1x
                                  "darkgreen", #Fe 2x
                                  "purple", #At:Fe
                                  "grey" #Controls
                                  )) +
      geom_line(position=position_dodge(.25), size=1) +
      geom_point(position=position_dodge(.25), size=3.5) +
      geom_errorbar(aes(ymin=mean-st.err,
                        ymax=mean+st.err),
                    width=.2, position=position_dodge(.25))
    
 #Just bt eggs over time in atrazine treatments
    bteggs.ag.at<-subset(bteggs.ag, Treatment=="A 1x" | 
                         Treatment=="A 2x" |
                         Treatment=="Control")
    ggplot(bteggs.ag.at, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
      theme_bw()+
      scale_color_manual(values=c("yellow", #At 1x
                                  "gold3", #At 2x
                                  "grey" #Controls
      )) +
      geom_line(position=position_dodge(.25), size=1) +
      geom_point(position=position_dodge(.25), size=3.5) +
      geom_errorbar(aes(ymin=mean-st.err,
                        ymax=mean+st.err),
                    width=.2, position=position_dodge(.25))
    
 #Hatchlings over time 
  ggplot(hatchs.ag, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
    scale_color_manual(values=cbPallete) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.25))
  
 #Hatchlings over time in individual tanks
  ggplot(hatchs, aes(x=variable, y=measure, group=tank, color=treat)) +
    theme_bw()+
    scale_color_manual(values=cbPallete) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5)
  
  
  
  #Linear models of time as predictor to peak growth (in terms of hatchlings produced)
  #before or at 7 weeks post first measurement #########################
    days<-c(0,7,14,21)
      #ChlorP only
        c.ch<-log(c(0.1,0.6,5.4,9.6))
        c.ch.m<-glm(c.ch ~ days)
          summary(c.ch.m) #slope = 0.227 +/- 0.03427
        plot(x=days, y=c.ch) 
          abline(a=c.ch.m$coefficients[1], b=c.ch.m$coefficients[2])
        
      #Ch:At
        c.chat<-log(c(0.1,1.2,14.8,26.2))
        c.chat.m<-glm(c.chat ~ days)
          summary(c.chat.m) #slope = 0.27453 +/- 0.04758
        plot(x=days, y=c.chat) 
          abline(a=c.chat.m$coefficients[1], b=c.chat.m$coefficients[2])
          
      #Ch:Fe
        c.chfe<-log(c(0.1,5,38.2,95))
        c.chfe.m<-glm(c.chfe ~ days)
          summary(c.chfe.m) #slope = 0.32290 +/- 0.06821
        plot(x=days, y=c.chfe) 
          abline(a=c.chfe.m$coefficients[1], b=c.chfe.m$coefficients[2])    
  
      #At:Ch:Fe
        c.tre<-log(c(0.1,3.6,26.4))
        c.tre.m<-glm(c.tre ~ days[c(1:3)])
          summary(c.tre.m) #slope = 0.39828 +/- 0.06562
        plot(x=days[c(1:3)], y=c.tre) 
          abline(a=c.tre.m$coefficients[1], b=c.tre.m$coefficients[2]) 
    
#What if we account for total number of snails in growth rates? Doesn't seem right....###############
  mean.fin<-data.frame("Treatment"=rep(0,12), 
                       "All_Snails"=rep(0,12))
  for(i in 1:length(snail.treats)){
      mean.fin[i,1]=as.character(snail.treats[i])
      mean.fin[i,2]=mean(snail$fin_total[snail$treat==snail.treats[i]])
    }
  
  eggs$per_cap=0
  for(i in 1:length(snail.treats)){
    eggs$per_cap[eggs$treat==snail.treats[i]]=(eggs$measure[eggs$treat==snail.treats[i]] /
                                                 mean.fin$All_Snails[mean.fin$Treatment==snail.treats[i]])
  }
  
  eggs.ag1.2<-aggregate.data.frame(eggs, by=list(eggs[,2], eggs[,4]), FUN=mean) 
    eggs.ag1.2<-eggs.ag1.2[,-c(3:7)] #remove unneeded variables 
      colnames(eggs.ag1.2)<-c("Treatment", "variable", "mean")
  
  eggs.ag2.2<-aggregate.data.frame(eggs, by=list(eggs[,2], eggs[,4]), FUN=st.er) #calculate st.error of treatment groups
    eggs.ag2.2<-eggs.ag2.2[,-c(3:7)] #remove unneeded variables 
      colnames(eggs.ag2.2)<-c("Treatment", "variable", "st.err")  
  
  eggs.ag.2<-merge(eggs.ag1.2, eggs.ag2.2, by=c("Treatment", "variable"))
  
  ggplot(eggs.ag.2, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
    scale_color_manual(values=cbPallete) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.25))
  