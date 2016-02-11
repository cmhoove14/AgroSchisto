#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Investigating results of Halstead et al 2016 mesocosm experiment
#Data file "R_use.csv" contains results of experiment investigated here
#email choover@berkeley.edu to obtain

require(ggplot2)
require(reshape2)

st.er <- function(x) {
  sd(x)/sqrt(length(x))
} #Function to calculate standard error of the mean

dat<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/R_use.csv")

#Add some variables ##########################
for(i in 1:nrow(dat)){
  dat[i,43]=paste(dat[i,2], dat[i,3], dat[i,4], sep="_")
} #Add unique treatment code to data frame
  treats<-unique(dat[,43])
  
for(i in 1:nrow(dat)){
    dat[i,44]=100*(dat[i,30]/dat[i,29])
} #Calculate final prevalence of B. glabrata
  colnames(dat)[44]<-"bg_prev"
  dat$bg_prev[is.na(dat$bg_prev)]<-0
  
for(i in 1:nrow(dat)){
    dat[i,45]=100*(dat[i,33]/dat[i,32])
} #Calculate final prevalence of B. truncatus
  colnames(dat)[45]<-"bt_prev"
  dat$bt_prev[is.na(dat$bt_prev)]<-0
    
varbs<-colnames(dat) #retrieve variable names

#Reshape to long format and merge data set to prepare for plotting #########################
dat2<-reshape(dat, varying = varbs[c(13:42,44,45)],
              v.names="measure",
              timevar="variable",
              times=varbs[c(13:42,44,45)],
              new.row.names=1:1920,
              direction="long")

dat2<-dat2[,-c(1:12,16)] #Get rid of needless variables

aggdata1<-aggregate.data.frame(dat2, by=list(dat2[,1], dat2[,2]), FUN=mean) #calculate means of treatment groups
  aggdata1<-aggdata1[,-c(3,4)] #remove unneeded variables 
  colnames(aggdata1)<-c("Treatment", "variable", "mean")
  
aggdata2<-aggregate.data.frame(dat2, by=list(dat2[,1], dat2[,2]), FUN=st.er) #calculate st.error of treatment groups
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

  

#Obtain mean and st. error of all variables within each treatment group in wide format ###########
aggdata3<-aggregate.data.frame(dat, by=list(dat[,43]), FUN=mean) #calculate means of treatment groups
  aggdata3<-aggdata3[,-c(2:13,44)] #remove unneeded variables 
  agg3.v<-paste(varbs[13:42], "mean", sep="_")
  colnames(aggdata3)[c(2:31)]<-agg3.v #change columns names to include mean
  
aggdata4<-aggregate.data.frame(dat, by=list(dat[,43]), FUN=st.er) #calculate st. errors of treatment groups
  aggdata4<-aggdata4[,-c(2:13,44)] #remove unneeded variables 
  agg4.v<-paste(varbs[13:42], "st.er", sep="_")
  colnames(aggdata4)[c(2:31)]<-agg4.v #change columns names to include mean

aggdata5<-merge(aggdata3, aggdata4, by.x="Group.1",by.y="Group.1") 

#Plot to visualize differences between treatment groups #######################
cbPalette <- c("#999999", "yellow", "red", "green", "orange", "darkgreen", "pink", "blue")

#Predators alive at end #################
preds.fin<-subset(aggdata, 
                  variable=="p.all_fin" | variable=="b.flu_fin")

  preds.fin$variable[preds.fin$variable=="p.all_fin"]<-"P. alleni"
  preds.fin$variable[preds.fin$variable=="b.flu_fin"]<-"B. flumineum"

  preds.fin$Treatment<- factor(preds.fin$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                           "Fertilizer","Atrazine_ChlorP",
                                                           "Atrazine_Fertilizer",
                                                           "ChlorP_Fertilizer","All_Three"))

  ggplot(preds.fin, aes(x=variable, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width = .7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Predators alive at end")
  
#Predators alive at 24 hours ########################
preds.24<-subset(aggdata, 
                 variable=="p.all_24" | variable=="b.flu_24")

  preds.24$variable[preds.24$variable=="p.all_24"]<-"P. alleni"
  preds.24$variable[preds.24$variable=="b.flu_24"]<-"B. flumineum"

  preds.24$Treatment<- factor(preds.24$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                             "Fertilizer","Atrazine_ChlorP",
                                                             "Atrazine_Fertilizer",
                                                             "ChlorP_Fertilizer","All_Three"))
  
  ggplot(preds.24, aes(x=variable, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                    ymax=mean+st.err),
                width=.2, position=position_dodge(.7)) +
    ggtitle("Predators alive at 24 hrs")
#Total snail reproduction ######################
snail.repro<-subset(aggdata, 
                 variable=="bg_eggs" | variable=="bg_hatch" | variable=="bt_eggs" | variable=="bt_hatch")

  snail.repro$variable[snail.repro$variable=="bg_eggs"]<-"B. glabrata eggs"
  snail.repro$variable[snail.repro$variable=="bg_hatch"]<-"B. glabrata hatchlings"
  snail.repro$variable[snail.repro$variable=="bt_eggs"]<-"B. truncatus eggs"
  snail.repro$variable[snail.repro$variable=="bt_hatch"]<-"B. truncatus hatchlings"
  
  snail.repro$Treatment<- factor(snail.repro$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                           "Fertilizer","Atrazine_ChlorP",
                                                           "Atrazine_Fertilizer",
                                                           "ChlorP_Fertilizer","All_Three"))
  
  ggplot(snail.repro, aes(x=variable, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                    ymax=mean+st.err),
                width=.2, position=position_dodge(.7)) +
    ggtitle("Snail reproduction at end")
  
#End snails alive between two species ########################
snail.live_fin<-subset(aggdata, 
                  variable=="bg_liv_fin" |  variable=="bt_liv_fin")

  snail.live_fin$variable[snail.live_fin$variable=="bg_liv_fin"]<-"B. glabrata"
  snail.live_fin$variable[snail.live_fin$variable=="bt_liv_fin"]<-"B. truncatus"

  snail.live_fin$Treatment<- factor(snail.live_fin$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                 "Fertilizer","Atrazine_ChlorP",
                                                                 "Atrazine_Fertilizer",
                                                                 "ChlorP_Fertilizer","All_Three"))
  
  
  ggplot(snail.live_fin, aes(x=variable, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                    ymax=mean+st.err),
                width=.2, position=position_dodge(.7)) +
    ggtitle("Snails alive at end")
#End snails total between two species ########################
snail.tot_fin<-subset(aggdata, 
                         variable=="total_bg_fin" |  variable=="total_bt_fin")
  
  snail.tot_fin$variable[snail.tot_fin$variable=="total_bg_fin"]<-"B. glabrata"
  snail.tot_fin$variable[snail.tot_fin$variable=="total_bt_fin"]<-"B. truncatus"
  
  snail.tot_fin$Treatment<- factor(snail.tot_fin$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                       "Fertilizer","Atrazine_ChlorP",
                                                                       "Atrazine_Fertilizer",
                                                                       "ChlorP_Fertilizer","All_Three"))
  
  ggplot(snail.tot_fin, aes(x=variable, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Total snails (alive/dead) at end")  
#End snails infected between two species ###################### 
snail.inf_fin<-subset(aggdata, 
                         variable=="bg_inf_fin" |  variable=="bt_inf_fin")
  
  snail.inf_fin$variable[snail.inf_fin$variable=="bg_inf_fin"]<-"B. glabrata"
  snail.inf_fin$variable[snail.inf_fin$variable=="bt_inf_fin"]<-"B. truncatus"
  
  snail.inf_fin$Treatment<- factor(snail.inf_fin$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                     "Fertilizer","Atrazine_ChlorP",
                                                                     "Atrazine_Fertilizer",
                                                                     "ChlorP_Fertilizer","All_Three"))
  
  ggplot(snail.inf_fin, aes(x=variable, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Snails infected at end")
#End snail infection prevalence (infections/100 snails) ####################
snail.inf_prev<-subset(aggdata, 
                        variable=="bg_prev" |  variable=="bt_prev")
 
  snail.inf_prev$variable[snail.inf_prev$variable=="bg_prev"]<-"B. glabrata"
  snail.inf_prev$variable[snail.inf_prev$variable=="bt_prev"]<-"B. truncatus"
  
  snail.inf_prev$Treatment<- factor(snail.inf_prev$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                     "Fertilizer","Atrazine_ChlorP",
                                                                     "Atrazine_Fertilizer",
                                                                     "ChlorP_Fertilizer","All_Three"))
  
  ggplot(snail.inf_prev, aes(x=variable, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Snail infection prev (inf/100 snails)")
#Periphyton levels measured across time ########################
peri_time<-subset(aggdata, 
                      variable=="peri0" |  variable=="peri1" |  variable=="peri2"
                  |  variable=="peri4" |  variable=="peri8")
  
  peri_time$variable[peri_time$variable=="peri0"]<-"Week 0"
  peri_time$variable[peri_time$variable=="peri1"]<-"Week 1"
  peri_time$variable[peri_time$variable=="peri2"]<-"Week 2"
  peri_time$variable[peri_time$variable=="peri4"]<-"Week 4"
  peri_time$variable[peri_time$variable=="peri8"]<-"Week 8"
  
  peri_time$Treatment<- factor(peri_time$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                       "Fertilizer","Atrazine_ChlorP",
                                                                       "Atrazine_Fertilizer",
                                                                       "ChlorP_Fertilizer","All_Three"))
  
ggplot(peri_time, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
  theme_bw()+
  scale_color_manual(values=cbPalette) +
  geom_line(position=position_dodge(.25), size=1) +
  geom_point(position=position_dodge(.25), size=3.5) +
  geom_errorbar(aes(ymin=mean-st.err,
                    ymax=mean+st.err),
                width=.2, position=position_dodge(.25)) +
  ggtitle("Periphyton chlorophyl a over time")  
#Phytoplankton levels measured across time ########################
phyto_time<-subset(aggdata, 
                  variable=="phyto0" |  variable=="phyto1" |  variable=="phyto2"
                  |  variable=="phyto4" |  variable=="phyto8")

phyto_time$variable[phyto_time$variable=="phyto0"]<-"Week 0"
phyto_time$variable[phyto_time$variable=="phyto1"]<-"Week 1"
phyto_time$variable[phyto_time$variable=="phyto2"]<-"Week 2"
phyto_time$variable[phyto_time$variable=="phyto4"]<-"Week 4"
phyto_time$variable[phyto_time$variable=="phyto8"]<-"Week 8"

phyto_time$Treatment<- factor(phyto_time$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                           "Fertilizer","Atrazine_ChlorP",
                                                           "Atrazine_Fertilizer",
                                                           "ChlorP_Fertilizer","All_Three"))

ggplot(phyto_time, aes(x=variable, y=mean, group=Treatment, color=Treatment)) +
  theme_bw()+
  scale_color_manual(values=cbPalette) +
  geom_line(position=position_dodge(.25), size=1) +
  geom_point(position=position_dodge(.25), size=3.5) +
  geom_errorbar(aes(ymin=mean-st.err,
                    ymax=mean+st.err),
                width=.2, position=position_dodge(.25)) +
  ggtitle("Phytoplankton chlorophyl a over time")  

#Does number of snails predict number of infected snails? ########################
plot(x=dat$bg_liv_fin, y=dat$bg_inf_fin, cex=0.3, 
     xlab = "Number snails", ylab ="Infected snails",
     main="Biomphalaria glabrata")

mod_bg<-lm(dat$bg_inf_fin ~ dat$bg_liv_fin)
summary(mod_bg) #Linear regression for all data points

mod_bg2<-lm(dat$bg_inf_fin[dat$bg_liv_fin >0] ~ dat$bg_liv_fin[dat$bg_liv_fin >0])
summary(mod_bg2) #Linear regression for replicates without snail pop crash

  abline(a=mod_bg$coefficients[1], b=mod_bg$coefficients[2], col="red")
  points(x=dat$bg_liv_fin[dat$atra==1], y=dat$bg_inf_fin[dat$atra==1],
       pch=17, cex=1.5, col="darkgray")
  points(x=dat$bg_liv_fin[dat$chlor==1], y=dat$bg_inf_fin[dat$chlor==1],
       pch=16, col="red")
  points(x=dat$bg_liv_fin[dat$fert==1], y=dat$bg_inf_fin[dat$fert==1],
         pch=4, cex=1.1, col="green")
  legend("bottomright",legend=c("atrazine present", "chlorpyrifos present", "fertilizer present"),
         pch=c(17,16,4), col=c("darkgray", "red", "green"))

plot(x=dat$bt_liv_fin, y=dat$bt_inf_fin, cex=0.3,
     xlab = "Number snails", ylab ="Infected snails")
mod_bt<-lm(dat$bt_inf_fin ~ dat$bt_liv_fin)
summary(mod_bt)
  abline(a=mod_bt$coefficients[1], b=mod_bt$coefficients[2], col="red")
  points(x=dat$bt_liv_fin[dat$atra==1], y=dat$bt_inf_fin[dat$atra==1],
       pch=17, cex=1.5, col="darkgray")
  points(x=dat$bt_liv_fin[dat$chlor==1], y=dat$bt_inf_fin[dat$chlor==1],
       pch=16, col="red")
  points(x=dat$bt_liv_fin[dat$fert==1], y=dat$bt_inf_fin[dat$fert==1],
       pch=4, cex=1.1, col="green")
  legend("bottomright",legend=c("atrazine present", "chlorpyrifos present", "fertilizer present"),
       pch=c(17,16,4), col=c("darkgray", "red", "green"))

#Create subset data frames for each treatment ############################
control<-subset(dat, atra==0 & chlor==0 & fert==0)

atra<-subset(dat, atra==1 & chlor==0 & fert==0)

chlor<-subset(dat, atra==0 & chlor==1 & fert==0)

fert<-subset(dat, atra==0 & chlor==0 & fert==1)

atra_fert<-subset(dat, atra==1 & chlor==0 & fert==1)

atra_chlor<-subset(dat, atra==1 & chlor==1 & fert==0)

chlor_fert<-subset(dat, atra==0 & chlor==1 & fert==1)

acf<-subset(dat, atra==1 & chlor==1 & fert==1)

treatments<-list(c("control", "atrazine_only", "chlorpyrifos_only", "fertilizer_only",
                   "atrazine_chlorpyrifos", "atrazine_fertilizer", "chlorpyrifos_fertilizer",
                   "all_three"), c("mean", "st.error"))