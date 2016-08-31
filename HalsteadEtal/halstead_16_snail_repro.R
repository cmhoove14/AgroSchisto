#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Investigating results of Halstead et al 2016 mesocosm experiment snail reproduction data
#Data file "snail_repro.csv" contains results of experiment investigated here
#email choover@berkeley.edu to obtain

require(ggplot2)
require(reshape2)

st.er <- function(x) {
  sd(x)/sqrt(length(x))
} #Function to calculate standard error of the mean
cbPalette <- c("#999999", "yellow", "red", "green", "orange", "darkgreen", "pink", "blue")


derp<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/snail_repro.csv")

for(i in 1:nrow(derp)){
  derp[i,5]=paste(derp[i,2], derp[i,3], derp[i,4], sep="_")
} #Add unique treatment code to data frame
treats<-unique(derp[,5])
colnames(derp)[5]<-"Treatment"

varbs<-colnames(derp)

#Reshape to long format #########################
derp2<-reshape(derp, varying = varbs[c(6:65)],
              v.names="measure",
              timevar="variable",
              times=varbs[c(6:65)],
              new.row.names=1:3600,
              direction="long")

derp2<-derp2[,-c(1:4,8)] #Get rid of needless variables


#Function to check out distribution of certain variable measures within a treatment ##############
  dist<-function(data=derp2, varb, treat){
    hist(data$measure[data$variable == varb & data$Treatment==treat], breaks=10, main=varb)
  }
  dist(varb='TBthatch1', treat='1_1_1')
  dist(varb='TBthatch2', treat='1_1_1')
  dist(varb='TBthatch3', treat='1_1_1')
  dist(varb='TBthatch4', treat='1_1_1')
  dist(varb='TBthatch8', treat='1_1_1')
    
#Before merging, check out snail hatchling numbers between groups to inform scalar of snail repro #######
  derp3<-subset(derp, Treatment == "0_1_0" | Treatment == "1_1_0"
                | Treatment == "1_1_1" | Treatment == "0_1_1")

  derp3<-reshape(derp3, varying = varbs[c(6:65)],
               v.names="measure",
               timevar="variable",
               times=varbs[c(6:65)],
               new.row.names=1:1500,
               direction="long")
    derp3<-derp3[,-c(2:4,8)]
      chlor_tanks<-unique(derp3$tank)
    derp3$tank<-rep(seq(from=1, to=25, by=1), times=60)
    
    snail.scale<-data.frame('tank'=chlor_tanks, 
                            "treat"= derp3[c(1:25),2], 
                            "bt_hatch", 
                            "bt_eggs", 
                            "bg_hatch", 
                            "bg_eggs")

#merge data set to prepare for plotting ##########################
aggdata1<-aggregate.data.frame(derp2, by=list(derp2[,1], derp2[,2]), FUN=mean) #calculate means of treatment groups
aggdata1<-aggdata1[,-c(3,4)] #remove unneeded variables 
colnames(aggdata1)<-c("Treatment", "variable", "mean")

aggdata2<-aggregate.data.frame(derp2, by=list(derp2[,1], derp2[,2]), FUN=st.er) #calculate st.error of treatment groups
aggdata2<-aggdata2[,-c(3,4)] #remove unneeded variables 
colnames(aggdata2)<-c("Treatment", "variable", "st.err")  

aggdata<-merge(aggdata1, aggdata2, by=c("Treatment", "variable"))

  aggdata$Treatment<-as.character(aggdata$Treatment)
  
  aggdata$Treatment[aggdata$Treatment=="0_0_0"]<- "Control"
  aggdata$Treatment[aggdata$Treatment=="1_0_0"]<- "Atrazine"
  aggdata$Treatment[aggdata$Treatment=="0_1_0"]<- "ChlorP"
  aggdata$Treatment[aggdata$Treatment=="0_0_1"]<- "Fertilizer"
  aggdata$Treatment[aggdata$Treatment=="1_1_0"]<- "Atrazine_ChlorP"
  aggdata$Treatment[aggdata$Treatment=="1_0_1"]<- "Atrazine_Fertilizer"
  aggdata$Treatment[aggdata$Treatment=="0_1_1"]<- "ChlorP_Fertilizer"
  aggdata$Treatment[aggdata$Treatment=="1_1_1"]<- "All_Three"
  
  aggdata$Treatment<-as.factor(aggdata$Treatment)
#Plot to visualize differences between treatment groups #######################
#Biomphalaria glabrata egg masses #################
bg.eggs<-subset(aggdata, 
                  variable=="TBgegg1" | variable=="TBgegg2" | variable=="TBgegg3" | 
                  variable=="TBgegg4" | variable=="TBgegg8") 
  
  bg.eggs$variable[bg.eggs$variable=="TBgegg1"]<-"Week 1"
  bg.eggs$variable[bg.eggs$variable=="TBgegg2"]<-"Week 2"
  bg.eggs$variable[bg.eggs$variable=="TBgegg3"]<-"Week 3"
  bg.eggs$variable[bg.eggs$variable=="TBgegg4"]<-"Week 4"
  bg.eggs$variable[bg.eggs$variable=="TBgegg8"]<-"Week 8"
  colnames(bg.eggs)[2]<-"Time"
  
  bg.eggs$Treatment<- factor(bg.eggs$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                 "Fertilizer","Atrazine_ChlorP",
                                                                 "Atrazine_Fertilizer",
                                                                 "ChlorP_Fertilizer","All_Three"))

  ggplot(bg.eggs, aes(x=Time, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
    scale_color_manual(values=cbPalette) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.25)) +
    ggtitle("B. glabrata egg masses sampled over time")  
#Biomphalaria glabrata hatchlings #################
bg.hatch<-subset(aggdata, 
                  variable=="TBghatch1" | variable=="TBghatch2" | variable=="TBghatch3" | 
                    variable=="TBghatch4" | variable=="TBghatch8") 
  
  bg.hatch$variable[bg.hatch$variable=="TBghatch1"]<-"Week 1"
  bg.hatch$variable[bg.hatch$variable=="TBghatch2"]<-"Week 2"
  bg.hatch$variable[bg.hatch$variable=="TBghatch3"]<-"Week 3"
  bg.hatch$variable[bg.hatch$variable=="TBghatch4"]<-"Week 4"
  bg.hatch$variable[bg.hatch$variable=="TBghatch8"]<-"Week 8"
  colnames(bg.hatch)[2]<-"Time"
  
  bg.hatch$Treatment<- factor(bg.hatch$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                         "Fertilizer","Atrazine_ChlorP",
                                                         "Atrazine_Fertilizer",
                                                         "ChlorP_Fertilizer","All_Three"))
  
  ggplot(bg.hatch, aes(x=Time, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
    scale_color_manual(values=cbPalette) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.25)) +
    ggtitle("B. glabrata hatchlings sampled over time")  
#Biomphalaria glabrata adults #################
bg.adult<-subset(aggdata, 
                   variable=="TBgadult1" | variable=="TBgadult2" | variable=="TBgadult3" | 
                     variable=="TBgadult4" | variable=="TBgadult8") 
  
  bg.adult$variable[bg.adult$variable=="TBgadult1"]<-"Week 1"
  bg.adult$variable[bg.adult$variable=="TBgadult2"]<-"Week 2"
  bg.adult$variable[bg.adult$variable=="TBgadult3"]<-"Week 3"
  bg.adult$variable[bg.adult$variable=="TBgadult4"]<-"Week 4"
  bg.adult$variable[bg.adult$variable=="TBgadult8"]<-"Week 8"
  colnames(bg.adult)[2]<-"Time"
  
  bg.adult$Treatment<- factor(bg.adult$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                           "Fertilizer","Atrazine_ChlorP",
                                                           "Atrazine_Fertilizer",
                                                           "ChlorP_Fertilizer","All_Three"))
  
  ggplot(bg.adult, aes(x=Time, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
    scale_color_manual(values=cbPalette) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.25)) +
    ggtitle("Longitudinal sampling of B. glabrata adults")  
#Bulinus truncatus egg masses #################
bt.eggs<-subset(aggdata, 
                  variable=="TBtegg1" | variable=="TBtegg2" | variable=="TBtegg3" | 
                    variable=="TBtegg4" | variable=="TBtegg8") 
  
  bt.eggs$variable[bt.eggs$variable=="TBtegg1"]<-"Week 1"
  bt.eggs$variable[bt.eggs$variable=="TBtegg2"]<-"Week 2"
  bt.eggs$variable[bt.eggs$variable=="TBtegg3"]<-"Week 3"
  bt.eggs$variable[bt.eggs$variable=="TBtegg4"]<-"Week 4"
  bt.eggs$variable[bt.eggs$variable=="TBtegg8"]<-"Week 8"
  colnames(bt.eggs)[2]<-"Time"
  
  bt.eggs$Treatment<- factor(bt.eggs$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                         "Fertilizer","Atrazine_ChlorP",
                                                         "Atrazine_Fertilizer",
                                                         "ChlorP_Fertilizer","All_Three"))
  
  ggplot(bt.eggs, aes(x=Time, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=15))+
    scale_color_manual(values=cbPalette) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.25)) +
    ggtitle("B. truncatus egg masses sampled over time")
#Bulinus truncatus hatchlings #################
bt.hatch<-subset(aggdata, 
                   variable=="TBthatch1" | variable=="TBthatch2" | variable=="TBthatch3" | 
                     variable=="TBthatch4" | variable=="TBthatch8") 
  
  bt.hatch$variable[bt.hatch$variable=="TBthatch1"]<-"Week 1"
  bt.hatch$variable[bt.hatch$variable=="TBthatch2"]<-"Week 2"
  bt.hatch$variable[bt.hatch$variable=="TBthatch3"]<-"Week 3"
  bt.hatch$variable[bt.hatch$variable=="TBthatch4"]<-"Week 4"
  bt.hatch$variable[bt.hatch$variable=="TBthatch8"]<-"Week 8"
  colnames(bt.hatch)[2]<-"Time"
  
  bt.hatch$Treatment<- factor(bt.hatch$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                           "Fertilizer","Atrazine_ChlorP",
                                                           "Atrazine_Fertilizer",
                                                           "ChlorP_Fertilizer","All_Three"))
  
  ggplot(bt.hatch, aes(x=Time, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=15))+
    scale_color_manual(values=cbPalette) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.25)) +
    ggtitle("B. truncatus hatchlings sampled over time")  
  
  bt.hatch$log_mean= log(bt.hatch$mean+1)
  
  lm.ch<-lm(bt.hatch$log_mean[bt.hatch$Treatment=="ChlorP"][c(1:4)] ~ c(0,7,14,21))
    summary(lm.ch) #slope = 0.10395
    
  lm.ch.at<-lm(bt.hatch$log_mean[bt.hatch$Treatment=="Atrazine_ChlorP"][c(1:4)] ~ c(0,7,14,21))
    summary(lm.ch.at) #slope = 0.13433
    
  lm.ch.fe<-lm(bt.hatch$log_mean[bt.hatch$Treatment=="ChlorP_Fertilizer"][c(1:4)] ~ c(0,7,14,21))
    summary(lm.ch.fe) #slope = 0.14115
    
  lm.ch.atfe<-lm(bt.hatch$log_mean[bt.hatch$Treatment=="All_Three"][c(1:3)] ~ c(0,7,14))
    summary(lm.ch.atfe) #slope = 0.1895
    
#Bulinus truncatus adults #################
bt.adult<-subset(aggdata, 
                   variable=="TBtadult1" | variable=="TBtadult2" | variable=="TBtadult3" | 
                     variable=="TBtadult4" | variable=="TBtadult8") 
  
  bt.adult$variable[bt.adult$variable=="TBtadult1"]<-"Week 1"
  bt.adult$variable[bt.adult$variable=="TBtadult2"]<-"Week 2"
  bt.adult$variable[bt.adult$variable=="TBtadult3"]<-"Week 3"
  bt.adult$variable[bt.adult$variable=="TBtadult4"]<-"Week 4"
  bt.adult$variable[bt.adult$variable=="TBtadult8"]<-"Week 8"
  colnames(bt.adult)[2]<-"Time"
  
  
  bt.adult$Treatment<- factor(bt.adult$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                           "Fertilizer","Atrazine_ChlorP",
                                                           "Atrazine_Fertilizer",
                                                           "ChlorP_Fertilizer","All_Three"))
  
  ggplot(bt.adult, aes(x=Time, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
    scale_color_manual(values=cbPalette) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.25)) +
    ggtitle("Longitudinal sampling of B. truncatus adults")  
  
  
#Investigate/plot reproduction only in tanks containing chlorpyrifos ################
cbPalette2 <- c("red", "purple","orange", "pink", "blue")
  
d8a<-subset(derp, Treatment == "0_1_0" |
              Treatment == "0_1_1" |
              Treatment == "1_1_0" |
              Treatment == "1_1_1")
  
  d8a$Treatment<-as.character(d8a$Treatment)
  
  d8a$Treatment[d8a$Treatment=="0_1_0"]<- "ChlorP"
  d8a$Treatment[d8a$Treatment=="1_1_0"]<- "Atrazine_ChlorP"
  d8a$Treatment[d8a$Treatment=="0_1_1"]<- "ChlorP_Fertilizer"
  d8a$Treatment[d8a$Treatment=="1_1_1"]<- "All_Three"
  d8a$Treatment[d8a$tank==11]<- "ChlorP2x"
  d8a$Treatment[d8a$tank==21]<- "ChlorP2x"
  d8a$Treatment[d8a$tank==36]<- "ChlorP2x"
  d8a$Treatment[d8a$tank==47]<- "ChlorP2x"
  d8a$Treatment[d8a$tank==55]<- "ChlorP2x"
  
  d8a$Treatment<-as.factor(d8a$Treatment)
#Are there direct effects of chlorpyrifos on snail reproduction? ########### 
  
  #B. truncatus hatchling reproduction
    bt.rp<-d8a[,c(1,50,9,18,27,36,45)]
      colnames(bt.rp)[c(3:7)]<-c("Week1","Week2","Week3","Week4","Week8")
      
    bt.rp<-reshape(bt.rp, varying = colnames(bt.rp[3:7]),
                     v.names="Hatchlings",
                     timevar="Time",
                     times=colnames(bt.rp[3:7]),
                     new.row.names=1:300,
                     direction="long")
    
    bt.rp$Treatment<- factor(bt.rp$Treatment, levels=c("ChlorP", "ChlorP2x", "Atrazine_ChlorP",
                                                       "ChlorP_Fertilizer","All_Three"))  
   
      
    
#Tanks 36, 39, and 55 had predators survive ChlorP treatments. Did this affect hatchlings sampled?  ####
 bt.rp$pred.surv<-'no'
      bt.rp$pred.surv[bt.rp$tank==36]<-'yes'
      bt.rp$pred.surv[bt.rp$tank==39]<-'yes'
      bt.rp$pred.surv[bt.rp$tank==55]<-'yes'
      
    ggplot(bt.rp, aes(x=Time, y=Hatchlings, group=tank, color=pred.surv))+
      theme_bw()+
      geom_point(size=1.5)+
      geom_line()
#Aggregate data and plot #####################    
    bt.rp.agg1<-aggregate.data.frame(bt.rp, by=list(bt.rp[,2], bt.rp[,3]), FUN=mean) #calculate means of treatment groups
    bt.rp.agg1<-bt.rp.agg1[,-c(3:5,7)] #remove unneeded variables 
    colnames(bt.rp.agg1)<-c("Treatment", "Time", "mean")
    
    bt.rp.agg2<-aggregate.data.frame(bt.rp, by=list(bt.rp[,2], bt.rp[,3]), FUN=st.er) #calculate st.error of treatment groups
    bt.rp.agg2<-bt.rp.agg2[,-c(3:5,7)] #remove unneeded variables 
    colnames(bt.rp.agg2)<-c("Treatment", "Time", "st.err")  
    
    bt.rp.agg<-merge(bt.rp.agg1, bt.rp.agg2, by=c("Treatment", "Time"))
    
    ggplot(bt.rp.agg, aes(x=Time, y=mean, group=Treatment, color=Treatment)) +
        theme_bw()+
        scale_color_manual(values=cbPalette2) +
        geom_line(position=position_dodge(.25), size=1) +
        geom_point(position=position_dodge(.25), size=3.5) +
        geom_errorbar(aes(ymin=mean-st.err,
                        ymax=mean+st.err),
                    width=.2, position=position_dodge(.25)) +
        ggtitle("B. truncatus hatchlings sampled over time")
    
  #B. glabrata hatchling reproduction
    bg.rp<-d8a[,c(1,50,6,15,24,33,42)]
    colnames(bg.rp)[c(3:7)]<-c("Week1","Week2","Week3","Week4","Week8")
    
    bg.rp<-reshape(bg.rp, varying = colnames(bg.rp[3:7]),
                   v.names="Hatchlings",
                   timevar="Time",
                   times=colnames(bg.rp[3:7]),
                   new.row.names=1:300,
                   direction="long")
    
    bg.rp$Treatment<- factor(bg.rp$Treatment, levels=c("ChlorP", "ChlorP2x", "Atrazine_ChlorP",
                                                       "ChlorP_Fertilizer","All_Three"))  
    
    bg.rp.agg1<-aggregate.data.frame(bg.rp, by=list(bg.rp[,2], bg.rp[,3]), FUN=mean) #calculate means of treatment groups
    bg.rp.agg1<-bg.rp.agg1[,-c(3:5,7)] #remove unneeded variables 
    colnames(bg.rp.agg1)<-c("Treatment", "Time", "mean")
    
    bg.rp.agg2<-aggregate.data.frame(bg.rp, by=list(bg.rp[,2], bg.rp[,3]), FUN=st.er) #calculate st.error of treatment groups
    bg.rp.agg2<-bg.rp.agg2[,-c(3:5,7)] #remove unneeded variables 
    colnames(bg.rp.agg2)<-c("Treatment", "Time", "st.err")  
    
    bg.rp.agg<-merge(bg.rp.agg1, bg.rp.agg2, by=c("Treatment", "Time"))
    
    ggplot(bg.rp.agg, aes(x=Time, y=mean, group=Treatment, color=Treatment)) +
      theme_bw()+
      scale_color_manual(values=cbPalette2) +
      geom_line(position=position_dodge(.25), size=1) +
      geom_point(position=position_dodge(.25), size=3.5) +
      geom_errorbar(aes(ymin=mean-st.err,
                        ymax=mean+st.err),
                    width=.2, position=position_dodge(.25)) +
      ggtitle("B. glabrata hatchlings sampled over time")
#Does not appear to be any significant effect of 2x chlorpyrifos dose on generation of hatchlings, implying
# no effect of chlorpyrifos on reproduction as was found in Ibrahim 1992 where reproduction declined with increasing dose 
#Plot individual tank results before aggregating ###############################
  #B. truncatus egg prpoduction
    bt.egg<-derp[,c(1,50,8,17,26,35,44)]
    colnames(bt.egg)[c(3:7)]<-c("Week1","Week2","Week3","Week4","Week8")
    
    bt.egg<-reshape(bt.egg, varying = colnames(bt.egg[3:7]),
                   v.names="Eggs",
                   timevar="Time",
                   times=colnames(bt.egg[3:7]),
                   new.row.names=1:300,
                   direction="long")
    
    bt.egg$Treatment[bt.egg$Treatment=="0_0_0"]<- "Control"
    bt.egg$Treatment[bt.egg$Treatment=="1_0_0"]<- "Atrazine"
    bt.egg$Treatment[bt.egg$Treatment=="0_1_0"]<- "ChlorP"
    bt.egg$Treatment[bt.egg$Treatment=="0_0_1"]<- "Fertilizer"
    bt.egg$Treatment[bt.egg$Treatment=="1_1_0"]<- "Atrazine_ChlorP"
    bt.egg$Treatment[bt.egg$Treatment=="1_0_1"]<- "Atrazine_Fertilizer"
    bt.egg$Treatment[bt.egg$Treatment=="0_1_1"]<- "ChlorP_Fertilizer"
    bt.egg$Treatment[bt.egg$Treatment=="1_1_1"]<- "All_Three"
    
    bt.egg$Treatment<- factor(bt.egg$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                         "Fertilizer","Atrazine_ChlorP",
                                                         "Atrazine_Fertilizer",
                                                         "ChlorP_Fertilizer","All_Three"))
    
    ggplot(bt.egg, aes(x=Time, y=Eggs, group=tank, colour=Treatment))+
      theme_bw()+
      scale_color_manual(values=cbPalette) +
      geom_line(position=position_dodge(.25), size=1) +
      geom_point(position=position_dodge(.25), size=3.5) +
      ggtitle("B. truncatus egg samplers")  
    
  control<-subset(bt.egg, Treatment=='Control')  
    ggplot(control, aes(x=Time, y=Eggs, group=tank, color=tank))+
      theme_bw()+
      geom_line(position=position_dodge(.25), size=1) +
      geom_point(position=position_dodge(.25), size=3.5)
  atraz<-subset(bt.egg, Treatment=='Atrazine')
    ggplot(atraz, aes(x=Time, y=Eggs, group=tank, color=tank))+
      theme_bw()+
      geom_line(position=position_dodge(.25), size=1) +
      geom_point(position=position_dodge(.25), size=3.5)
  chlorpyr<-subset(bt.egg, Treatment=='ChlorP')
    ggplot(chlorpyr, aes(x=Time, y=Eggs, group=tank, color=tank))+
      theme_bw()+
      geom_line(position=position_dodge(.25), size=1) +
      geom_point(position=position_dodge(.25), size=3.5)
  fertiliz<-subset(bt.egg, Treatment=='Fertilizer')
    ggplot(fertiliz, aes(x=Time, y=Eggs, group=tank, color=tank))+
      theme_bw()+
      geom_line(position=position_dodge(.25), size=1) +
      geom_point(position=position_dodge(.25), size=3.5)
  atra_chlorP<-subset(bt.egg, Treatment=='Atrazine_ChlorP')
    ggplot(atra_chlorP, aes(x=Time, y=Eggs, group=tank, color=tank))+
      theme_bw()+
      geom_line(position=position_dodge(.25), size=1) +
      geom_point(position=position_dodge(.25), size=3.5)
    