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

for(i in 1:nrow(dat)){
  dat[i,43]=paste(dat[i,2], dat[i,3], dat[i,4], sep="_")
} #Add unique treatment code to data frame
  treats<-unique(dat[,43])
varbs<-colnames(dat) #retrieve variable names

#Reshape to long format and merge data set to prepare for plotting #########################
dat2<-reshape(dat, varying = varbs[13:42],
              v.names="measure",
              timevar="variable",
              times=varbs[13:42],
              new.row.names=1:1800,
              direction="long")

dat2<-dat2[,-c(1:12,16)] #Get rid of needless variables

aggdata1<-aggregate.data.frame(dat2, by=list(dat2[,1], dat2[,2]), FUN=mean) #calculate means of treatment groups
  aggdata1<-aggdata1[,-c(3,4)] #remove unneeded variables 
  colnames(aggdata1)<-c("atra_chlor_fert", "variable", "mean")
  
aggdata2<-aggregate.data.frame(dat2, by=list(dat2[,1], dat2[,2]), FUN=st.er) #calculate st.error of treatment groups
  aggdata2<-aggdata2[,-c(3,4)] #remove unneeded variables 
  colnames(aggdata2)<-c("atra_chlor_fert", "variable", "st.err")  
  
aggdata<-merge(aggdata1, aggdata2, by=c("atra_chlor_fert", "variable"))
  

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
cbPalette <- c("#999999", "green", "red", "darkgreen", "yellow", "pink", "orange", "blue")

#Predators alive at end #################
preds.fin<-subset(aggdata, 
                  variable=="p.all_fin" | variable=="b.flu_fin" | variable=="all.pred_fin") 

  ggplot(preds.fin, aes(x=variable, y=mean, fill=atra_chlor_fert)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width = .7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Predators alive at end")
  
#Predators alive at 24 hours ########################
preds.24<-subset(aggdata, 
                 variable=="p.all_24" | variable=="b.flu_24" | variable=="all.pred_24")

  ggplot(preds.24, aes(x=variable, y=mean, fill=atra_chlor_fert)) +
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

  ggplot(snail.repro, aes(x=variable, y=mean, fill=atra_chlor_fert)) +
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

  ggplot(snail.live_fin, aes(x=variable, y=mean, fill=atra_chlor_fert)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                    ymax=mean+st.err),
                width=.2, position=position_dodge(.7)) +
    ggtitle("Snails alive at end")
  
#End snails infected between two species ###################### 
snail.inf_fin<-subset(aggdata, 
                         variable=="bg_inf_fin" |  variable=="bt_inf_fin")
  
  ggplot(snail.inf_fin, aes(x=variable, y=mean, fill=atra_chlor_fert)) +
    theme_bw()+
    scale_fill_manual(values=cbPalette) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=mean-st.err,
                      ymax=mean+st.err),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Snails infected at end")

  
  
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

###################
pred_surv<-matrix(ncol=2, nrow=8, dimnames=treatments)
  pred_surv[1,1]=mean(control$all.pred_fin)
  pred_surv[1,2]=st.er(control$all.pred_fin)
  
  pred_surv[2,1]=mean(atra$all.pred_fin)
  pred_surv[2,2]=st.er(atra$all.pred_fin)
  
  pred_surv[3,1]=mean(chlor$all.pred_fin)
  pred_surv[3,2]=st.er(chlor$all.pred_fin)
  
  pred_surv[4,1]=mean(fert$all.pred_fin)
  pred_surv[4,2]=st.er(fert$all.pred_fin)
  
  pred_surv[5,1]=mean(atra_chlor$all.pred_fin)
  pred_surv[5,2]=st.er(atra_chlor$all.pred_fin)
  
  pred_surv[6,1]=mean(atra_fert$all.pred_fin)
  pred_surv[6,2]=st.er(atra_fert$all.pred_fin)
  
  pred_surv[7,1]=mean(chlor_fert$all.pred_fin)
  pred_surv[7,2]=st.er(chlor_fert$all.pred_fin)
  
  pred_surv[8,1]=mean(acf$all.pred_fin)
  pred_surv[8,2]=st.er(acf$all.pred_fin)
  
tot_snails<-matrix(ncol=2, nrow=8, dimnames=treatments)  
  tot_snails[1,1]=mean(control$total_all_fin)
  tot_snails[1,2]=st.er(control$total_all_fin)
  
  tot_snails[2,1]=mean(atra$total_all_fin)
  tot_snails[2,2]=st.er(atra$total_all_fin)
  
  tot_snails[3,1]=mean(chlor$total_all_fin)
  tot_snails[3,2]=st.er(chlor$total_all_fin)
  
  tot_snails[4,1]=mean(fert$total_all_fin)
  tot_snails[4,2]=st.er(fert$total_all_fin)
  
  tot_snails[5,1]=mean(atra_chlor$total_all_fin)
  tot_snails[5,2]=st.er(atra_chlor$total_all_fin)
  
  tot_snails[6,1]=mean(atra_fert$total_all_fin)
  tot_snails[6,2]=st.er(atra_fert$total_all_fin)
  
  tot_snails[7,1]=mean(chlor_fert$total_all_fin)
  tot_snails[7,2]=st.er(chlor_fert$total_all_fin)
  
  tot_snails[8,1]=mean(acf$total_all_fin)
  tot_snails[8,2]=st.er(acf$total_all_fin)
