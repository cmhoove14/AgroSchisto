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
colnames(derp)[50]<-"atra_chlor_fert"

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
colnames(aggdata1)<-c("atra_chlor_fert", "variable", "mean")

aggdata2<-aggregate.data.frame(derp2, by=list(derp2[,1], derp2[,2]), FUN=st.er) #calculate st.error of treatment groups
aggdata2<-aggdata2[,-c(3,4)] #remove unneeded variables 
colnames(aggdata2)<-c("atra_chlor_fert", "variable", "st.err")  

aggdata<-merge(aggdata1, aggdata2, by=c("atra_chlor_fert", "variable"))

#Plot to visualize differences between treatment groups #######################
cbPalette <- c("#999999", "green", "red", "darkgreen", "yellow", "pink", "orange", "blue")

#Biomphalaria glabrata egg masses #################
bg.eggs<-subset(aggdata, 
                  variable=="TBgegg1" | variable=="TBgegg2" | variable=="TBgegg3" | 
                  variable=="TBgegg4" | variable=="TBgegg8") 

  ggplot(bg.eggs, aes(x=variable, y=mean, fill=atra_chlor_fert)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width = .7) +
    geom_errorbar(aes(ymin=mean-st.err,
                    ymax=mean+st.err),
                width=.2, position=position_dodge(.7)) +
    ggtitle("Longitudinal B. glabrata egg masses")
#Biomphalaria glabrata hatchlings #################
bg.hatch<-subset(aggdata, 
                  variable=="TBghatch1" | variable=="TBghatch2" | variable=="TBghatch3" | 
                    variable=="TBghatch4" | variable=="TBghatch8") 
  
  ggplot(bg.hatch, aes(x=variable, y=mean, fill=atra_chlor_fert)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width = .7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Longitudinal B. glabrata hatchlings")  
#Biomphalaria glabrata adults #################
bg.adult<-subset(aggdata, 
                   variable=="TBgadult1" | variable=="TBgadult2" | variable=="TBgadult3" | 
                     variable=="TBgadult4" | variable=="TBgadult8") 
  
  ggplot(bg.adult, aes(x=variable, y=mean, fill=atra_chlor_fert)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width = .7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Longitudinal B. glabrata adults")  
#Bulinus truncatus egg masses #################
bt.eggs<-subset(aggdata, 
                  variable=="TBtegg1" | variable=="TBtegg2" | variable=="TBtegg3" | 
                    variable=="TBtegg4" | variable=="TBtegg8") 
  
  ggplot(bt.eggs, aes(x=variable, y=mean, fill=atra_chlor_fert)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width = .7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Longitudinal B. truncatus egg masses")
#Bulinus truncatus hatchlings #################
bt.hatch<-subset(aggdata, 
                   variable=="TBthatch1" | variable=="TBthatch2" | variable=="TBthatch3" | 
                     variable=="TBthatch4" | variable=="TBthatch8") 
  
  ggplot(bt.hatch, aes(x=variable, y=mean, fill=atra_chlor_fert)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width = .7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Longitudinal B. truncatus hatchlings")  
#Bulinus truncatus adults #################
bt.adult<-subset(aggdata, 
                   variable=="TBtadult1" | variable=="TBtadult2" | variable=="TBtadult3" | 
                     variable=="TBtadult4" | variable=="TBtadult8") 
  
  ggplot(bt.adult, aes(x=variable, y=mean, fill=atra_chlor_fert)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width = .7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Longitudinal B. truncatus adults")  