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

derp<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/snail_repro.csv")

for(i in 1:nrow(derp)){
  derp[i,50]=paste(derp[i,2], derp[i,3], derp[i,4], sep="_")
} #Add unique treatment code to data frame
treats<-unique(derp[,50])
colnames(derp)[50]<-"Treatment"

varbs<-colnames(derp)

#Reshape to long format and merge data set to prepare for plotting #########################
derp2<-reshape(derp, varying = varbs[c(5:49)],
              v.names="measure",
              timevar="variable",
              times=varbs[c(5:49)],
              new.row.names=1:2700,
              direction="long")

derp2<-derp2[,-c(1:4,8)] #Get rid of needless variables

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
  cbPalette <- c("#999999", "yellow", "red", "green", "orange", "darkgreen", "pink", "blue")
  
#Biomphalaria glabrata egg masses #################
bg.eggs<-subset(aggdata, 
                  variable=="TBgegg1" | variable=="TBgegg2" | variable=="TBgegg3" | 
                  variable=="TBgegg4" | variable=="TBgegg8") 
  
  bg.eggs$variable[bg.eggs$variable=="TBgegg1"]<-"Week 1"
  bg.eggs$variable[bg.eggs$variable=="TBgegg2"]<-"Week 2"
  bg.eggs$variable[bg.eggs$variable=="TBgegg3"]<-"Week 3"
  bg.eggs$variable[bg.eggs$variable=="TBgegg4"]<-"Week 4"
  bg.eggs$variable[bg.eggs$variable=="TBgegg8"]<-"Week 8"
  
  bg.eggs$Treatment<- factor(bg.eggs$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                 "Fertilizer","Atrazine_ChlorP",
                                                                 "Atrazine_Fertilizer",
                                                                 "ChlorP_Fertilizer","All_Three"))

  ggplot(bg.eggs, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
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
  
  bg.hatch$Treatment<- factor(bg.hatch$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                         "Fertilizer","Atrazine_ChlorP",
                                                         "Atrazine_Fertilizer",
                                                         "ChlorP_Fertilizer","All_Three"))
  
  ggplot(bg.hatch, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
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
  
  bg.adult$Treatment<- factor(bg.adult$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                           "Fertilizer","Atrazine_ChlorP",
                                                           "Atrazine_Fertilizer",
                                                           "ChlorP_Fertilizer","All_Three"))
  
  ggplot(bg.adult, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
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
  
  bt.eggs$Treatment<- factor(bt.eggs$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                         "Fertilizer","Atrazine_ChlorP",
                                                         "Atrazine_Fertilizer",
                                                         "ChlorP_Fertilizer","All_Three"))
  
  ggplot(bt.eggs, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
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
  
  bt.hatch$Treatment<- factor(bt.hatch$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                           "Fertilizer","Atrazine_ChlorP",
                                                           "Atrazine_Fertilizer",
                                                           "ChlorP_Fertilizer","All_Three"))
  
  ggplot(bt.hatch, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
    scale_color_manual(values=cbPalette) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.25)) +
    ggtitle("B. truncatus hatchlings sampled over time")  
#Bulinus truncatus adults #################
bt.adult<-subset(aggdata, 
                   variable=="TBtadult1" | variable=="TBtadult2" | variable=="TBtadult3" | 
                     variable=="TBtadult4" | variable=="TBtadult8") 
  
  bt.adult$variable[bt.adult$variable=="TBtadult1"]<-"Week 1"
  bt.adult$variable[bt.adult$variable=="TBtadult2"]<-"Week 2"
  bt.adult$variable[bt.adult$variable=="TBtadult3"]<-"Week 3"
  bt.adult$variable[bt.adult$variable=="TBtadult4"]<-"Week 4"
  bt.adult$variable[bt.adult$variable=="TBtadult8"]<-"Week 8"
  
  bt.adult$Treatment<- factor(bt.adult$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                           "Fertilizer","Atrazine_ChlorP",
                                                           "Atrazine_Fertilizer",
                                                           "ChlorP_Fertilizer","All_Three"))
  
  ggplot(bt.adult, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
    scale_color_manual(values=cbPalette) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.25)) +
    ggtitle("Longitudinal sampling of B. truncatus adults")  