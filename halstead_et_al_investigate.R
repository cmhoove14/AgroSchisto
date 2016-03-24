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
require(reshape)

st.er <- function(x) {
  sd(x)/sqrt(length(x))
} #Function to calculate standard error of the mean
cbPalette <- c("#999999", "yellow", "red", "green", "orange", "darkgreen", "pink", "blue")

dat<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/R_use.csv")

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
  
for(i in 1:nrow(dat)){
    dat[i,46]=dat[i,31]-dat[i,29]
} #Calculate number of dead B. glabrata
  colnames(dat)[46]<-"bg_dead"

for(i in 1:nrow(dat)){
    dat[i,47]=dat[i,34]-dat[i,32]
} #Calculate number of dead B. truncatus
  colnames(dat)[47]<-"bt_dead"
  
for(i in 1:nrow(dat)){
    dat[i,48]=(log((dat[i,29]+1)/27)/84)
} #Calculate intrinsic reproduction rate (r) of B. glabrata based on 12 week increase in population
  colnames(dat)[48]<-"bg_r"
  dat$bg_r[dat$bg_r <= 0]<-0 #relabel negative population growth as 0
  
for(i in 1:nrow(dat)){
    dat[i,49]=(log((dat[i,32]+1)/11)/84)
} #Calculate intrinsic reproduction rate (r) of B. glabrata based on 12 week increase in population
  colnames(dat)[49]<-"bt_r"
  dat$bt_r[dat$bt_r <= 0]<-0#relabel negative population growth as 0

for(i in 1:nrow(dat)){
  dat[i,50]=mean(dat[i,18],dat[i,19],dat[i,20],dat[i,21],dat[i,22])
}  #Calculate average periphyton chlorophyl-a over time period
  colnames(dat)[50]<-"peri_ave"
  
for(i in 1:nrow(dat)){
    dat[i,51]=dat[i,25]/dat[i,32]
}  #Calculate estimate of eggs/snail produced through whole experiment to account for higher number of snails in Chlor treatments
  colnames(dat)[51]<-"bt_egg_perCap"
    dat$bt_egg_perCap[dat$bt_egg_perCap==Inf]<-NA
    dat$bt_egg_perCap[is.na(dat$bt_egg_perCap)==TRUE]<-0
  
varbs<-colnames(dat) #retrieve variable names

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

#Check out frequency of predators surviving insecticide treatments #################
chlor.all<-subset(dat, chlor==1)
  chlor.all$tank[chlor.all$all.pred_fin >=1] #Tanks 36, 39, & 55 had preds despite ChlorP presence

#24 hour mortality of prawns in treatment groups
  ggplot(dat, aes(x=V43, y=p.all_24, fill=V43, group=V43))+
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width = .7)

#Do periphyton levels predict final snail density (absent predation influence)? ########################
  chlorP<-subset(dat, chlor==1)
    plot(x=chlorP$peri_ave, y=chlorP$bt_liv_fin, 
         pch=16, col="blue",
         xlab="Mean Periphyton chlorophyl-a", ylab="B. truncatus alive at end") 
    points(x=chlorP$peri_ave, y=chlorP$bg_liv_fin,
           pch=16, col="red")

#Reshape to long format and merge data set to prepare for plotting #########################
dat2<-reshape(dat, varying = varbs[c(13:42,44:49,51)],
              v.names="measure",
              timevar="variable",
              times=varbs[c(13:42,44:49,51)],
              new.row.names=1:2220,
              direction="long")

dat2<-dat2[,-c(1:12,14,17)] #Get rid of needless variables

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

#Predators alive at end #################
preds.fin<-subset(aggdata, 
                  variable=="p.all_fin" | variable=="b.flu_fin")

  preds.fin$variable[preds.fin$variable=="p.all_fin"]<-"P. alleni"
  preds.fin$variable[preds.fin$variable=="b.flu_fin"]<-"B. flumineum"
  colnames(preds.fin)[2]<-"Species"
  
  preds.fin$Treatment<- factor(preds.fin$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                           "Fertilizer","Atrazine_ChlorP",
                                                           "Atrazine_Fertilizer",
                                                           "ChlorP_Fertilizer","All_Three"))

  ggplot(preds.fin, aes(x=Species, y=mean, fill=Treatment)) +
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
  colnames(preds.24)[2]<-"Species"
  
  preds.24$Treatment<- factor(preds.24$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                             "Fertilizer","Atrazine_ChlorP",
                                                             "Atrazine_Fertilizer",
                                                             "ChlorP_Fertilizer","All_Three"))
  
  ggplot(preds.24, aes(x=Species, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                    ymax=mean+st.err),
                width=.2, position=position_dodge(.7)) +
    ggtitle("Predators alive at 24 hrs")
#What is the observed daily prawn mortality rate in chlorP-free tanks over the full 12 weeks?
  p.0<-3*length(dat$tank[dat$chlor==0])
  p.12<-sum(dat$p.all_fin[dat$chlor==0])
    #daily mortality rate assuming constant death throughout 12 weeks =ln(Nt/N0)/-t with t=12 weeks *7 days=84 days
      p.r<-log(p.12/p.0)/-84 #=0.006862177
  
#Check out proportion of P. alleni surviving after 24 hours in ChlorP presence/absence ####################
  
  prawn.tox<-data.frame('chlorP'=dat$chlor,
                        'chlorP2'=dat$treat,
                        'prawn.0'=rep(3,length(dat$chlor)),
                        'prawn.24'=dat$p.all_24)
  
  prawn.tox$chlorP[prawn.tox$chlorP2=="C 2x"]<-2
  
  prawn.tox<-prawn.tox[,-2]
  
  prawn.tox2<-untable(prawn.tox[,c(1,3)], num=prawn.tox[,2])
  
  prawn.tox2$dose[prawn.tox2$chlorP==2]<-128
  prawn.tox2$dose[prawn.tox2$chlorP==1]<-64
  prawn.tox2$dose[prawn.tox2$chlorP==0]<-0
  
  prawn.tox2$outcome<-rep(0, length(prawn.tox2$chlorP))
  
  prawn.tox2$outcome[prawn.tox2$prawn.24==0]<-rep(c(1,1,1), length(prawn.tox2$prawn.24[prawn.tox2$prawn.24==0])/3)
  prawn.tox2$outcome[prawn.tox2$prawn.24==1]<-rep(c(1,1,0), length(prawn.tox2$prawn.24[prawn.tox2$prawn.24==1])/3)
  prawn.tox2$outcome[prawn.tox2$prawn.24==2]<-rep(c(1,0,0), length(prawn.tox2$prawn.24[prawn.tox2$prawn.24==2])/3)
  prawn.tox2$outcome[prawn.tox2$prawn.24==3]<-rep(c(0,0,0), length(prawn.tox2$prawn.24[prawn.tox2$prawn.24==3])/3)

#~*~*~*~DAILY MORTALITY RATES FOR PRAWN SPECIES~*~*~*~*~
  sum(prawn.tox2$outcome[prawn.tox2$chlorP==2])/length(prawn.tox2$outcome[prawn.tox2$chlorP==2]) #60% mortality when chlorP present
  -log(1-0.6) #Observed daily mortality rate = 0.9162907 when ChlorP @128 ug/L is present
  
  sum(prawn.tox2$outcome[prawn.tox2$chlorP==1])/length(prawn.tox2$outcome[prawn.tox2$chlorP==1]) #73.3% mortality when chlorP present
    -log(1-0.766667) #Observed daily mortality rate = 1.455289 when ChlorP @64 ug/L is present
  
  sum(prawn.tox2$outcome[prawn.tox2$chlorP!=0])/length(prawn.tox2$outcome[prawn.tox2$chlorP!=0]) #73.3% mortality when chlorP present
  -log(1-0.7333333) #Observed daily mortality rate = 1.321756 when ChlorP is present
  
  sum(prawn.tox2$outcome[prawn.tox2$chlorP==0])/length(prawn.tox2$outcome[prawn.tox2$chlorP==0]) #3.8% mortality when chlorP absent
    -log(1-0.03809524) #Observed daily mortality rate = 0.03883984 when ChlorP is absent

#Probit model with three tested doses
  pr<-glm(outcome ~ dose, family=binomial(link="probit"),data=prawn.tox2)
    summary(pr)
    
  tox.predict<-data.frame('dose'=seq(1,150,1),
                          'mort'=rep(0,150),
                          'st.er'=rep(0,150))
  
  tox.predict$mort<-predict(pr, tox.predict, type = "response", se.fit=TRUE)$fit
  tox.predict$st.er<-predict(pr, tox.predict, type = "response", se.fit=TRUE)$se.fit
  
  plot(tox.predict$dose, tox.predict$mort, type='l', lwd=1.5,
       xlab='Dose', ylab = '% mortality')
    lines(tox.predict$dose, tox.predict$mort+tox.predict$st.er, col='red', lwd=0.8, lty=2)
    lines(tox.predict$dose, tox.predict$mort-tox.predict$st.er, col='red', lwd=0.8, lty=2)
    
#Total snail reproduction ######################
snail.repro<-subset(aggdata, 
                 variable=="bg_eggs" | variable=="bg_hatch" | variable=="bt_eggs" | variable=="bt_hatch")

  snail.repro$variable[snail.repro$variable=="bg_eggs"]<-"B. glabrata eggs"
  snail.repro$variable[snail.repro$variable=="bg_hatch"]<-"B. glabrata hatchlings"
  snail.repro$variable[snail.repro$variable=="bt_eggs"]<-"B. truncatus eggs"
  snail.repro$variable[snail.repro$variable=="bt_hatch"]<-"B. truncatus hatchlings"
  colnames(snail.repro)[2]<-"Species"
  
  snail.repro$Treatment<- factor(snail.repro$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                           "Fertilizer","Atrazine_ChlorP",
                                                           "Atrazine_Fertilizer",
                                                           "ChlorP_Fertilizer","All_Three"))
  
  ggplot(snail.repro, aes(x=Species, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                    ymax=mean+st.err),
                width=.2, position=position_dodge(.7)) +
    ggtitle("Snail reproduction at end")
  
#"Global" estimate of reproduction rate (r) based on final number of live snails ##########
rep_r<-subset(aggdata, 
                         variable=="bt_r" |  variable=="bg_r")
  
  rep_r$variable[rep_r$variable=="bg_r"]<-"B. glabrata"
  rep_r$variable[rep_r$variable=="bt_r"]<-"B. truncatus"
  colnames(rep_r)[2]<-"Species"
  
  rep_r$Treatment<- factor(rep_r$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                       "Fertilizer","Atrazine_ChlorP",
                                                                       "Atrazine_Fertilizer",
                                                                       "ChlorP_Fertilizer","All_Three"))
  
  
  ggplot(rep_r, aes(x=Species, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Approximate daily reproduction rate over 12 week period")
#End snails alive between two species ########################
snail.live_fin<-subset(aggdata, 
                  variable=="bg_liv_fin" |  variable=="bt_liv_fin")

  snail.live_fin$variable[snail.live_fin$variable=="bg_liv_fin"]<-"B. glabrata"
  snail.live_fin$variable[snail.live_fin$variable=="bt_liv_fin"]<-"B. truncatus"
  colnames(snail.live_fin)[2]<-"Species"
  
  snail.live_fin$Treatment<- factor(snail.live_fin$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                 "Fertilizer","Atrazine_ChlorP",
                                                                 "Atrazine_Fertilizer",
                                                                 "ChlorP_Fertilizer","All_Three"))
  
  
  ggplot(snail.live_fin, aes(x=Species, y=mean, fill=Treatment)) +
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
  colnames(snail.tot_fin)[2]<-"Species"
  
  snail.tot_fin$Treatment<- factor(snail.tot_fin$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                       "Fertilizer","Atrazine_ChlorP",
                                                                       "Atrazine_Fertilizer",
                                                                       "ChlorP_Fertilizer","All_Three"))
  
  ggplot(snail.tot_fin, aes(x=Species, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Total snails (alive+dead) at end")  
#End snails dead between two species ########################
snail.dead_fin<-subset(aggdata, 
                        variable=="bg_dead" |  variable=="bt_dead")
  
  snail.dead_fin$variable[snail.dead_fin$variable=="bg_dead"]<-"B. glabrata"
  snail.dead_fin$variable[snail.dead_fin$variable=="bt_dead"]<-"B. truncatus"
  colnames(snail.dead_fin)[2]<-"Species"
  
  snail.dead_fin$Treatment<- factor(snail.dead_fin$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                     "Fertilizer","Atrazine_ChlorP",
                                                                     "Atrazine_Fertilizer",
                                                                     "ChlorP_Fertilizer","All_Three"))
  
  ggplot(snail.dead_fin, aes(x=Species, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Dead snails at end")  
#End snails infected between two species ###################### 
snail.inf_fin<-subset(aggdata, 
                         variable=="bg_inf_fin" |  variable=="bt_inf_fin")
  
  snail.inf_fin$variable[snail.inf_fin$variable=="bg_inf_fin"]<-"B. glabrata"
  snail.inf_fin$variable[snail.inf_fin$variable=="bt_inf_fin"]<-"B. truncatus"
  colnames(snail.inf_fin)[2]<-"Species"
  
  snail.inf_fin$Treatment<- factor(snail.inf_fin$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                     "Fertilizer","Atrazine_ChlorP",
                                                                     "Atrazine_Fertilizer",
                                                                     "ChlorP_Fertilizer","All_Three"))
  
  ggplot(snail.inf_fin, aes(x=Species, y=mean, fill=Treatment)) +
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
  colnames(snail.inf_prev)[2]<-"Species"
  
  snail.inf_prev$Treatment<- factor(snail.inf_prev$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                     "Fertilizer","Atrazine_ChlorP",
                                                                     "Atrazine_Fertilizer",
                                                                     "ChlorP_Fertilizer","All_Three"))
  
  ggplot(snail.inf_prev, aes(x=Species, y=mean, fill=Treatment)) +
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
  colnames(peri_time)[2]<-"Time"
  
  peri_time$Treatment<- factor(peri_time$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                                       "Fertilizer","Atrazine_ChlorP",
                                                                       "Atrazine_Fertilizer",
                                                                       "ChlorP_Fertilizer","All_Three"))
  
ggplot(peri_time, aes(x=Time, y=mean, group=Treatment, color=Treatment)) +
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
  colnames(phyto_time)[2]<-"Time"
  
  phyto_time$Treatment<- factor(phyto_time$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                             "Fertilizer","Atrazine_ChlorP",
                                                             "Atrazine_Fertilizer",
                                                             "ChlorP_Fertilizer","All_Three"))
  
  ggplot(phyto_time, aes(x=Time, y=mean, group=Treatment, color=Treatment)) +
    theme_bw()+
    scale_color_manual(values=cbPalette) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.25)) +
    ggtitle("Phytoplankton chlorophyl a over time")  
#B. truncatus per capita eggs ##################
  bt_eggs_per_cap<-subset(aggdata,variable=='bt_egg_perCap')
  
  bt_eggs_per_cap$Treatment<- factor(bt_eggs_per_cap$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                               "Fertilizer","Atrazine_ChlorP",
                                                               "Atrazine_Fertilizer",
                                                               "ChlorP_Fertilizer","All_Three"))
  
  ggplot(bt_eggs_per_cap, aes(x=variable, y=mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Bt eggs per cap at end")

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