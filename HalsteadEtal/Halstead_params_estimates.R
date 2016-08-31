#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Load libraries and small functions #############
require(ggplot2)
require(reshape)
st.er <- function(x) {
  sd(x)/sqrt(length(x))
} #Function to calculate standard error of the mean
cbPalette <- c("#999999", "yellow", "red", "green", "orange", "darkgreen", "pink", "blue")

#Table of Contents of parameter estimates ########################
  # 1) Estimate snail population dynamics parameters
    # Reproduction rate (agroC sensitive and constant options)
    # Carrying capacity (agroC sensitive and constant options)
    # Mortality rate
    # Infected mortality rate
    #Proportion of exposed that reproduce
    #

#Daily B. truncatus reproduction rate ###########
dat<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/R_use.csv")

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
bt.obs$Scalar<-c((bt.obs[1,2]/bt.obs[1,2]),
                 (bt.obs[2,2]/bt.obs[1,2]),
                 (bt.obs[3,2]/bt.obs[1,2]),
                 (bt.obs[4,2]/bt.obs[1,2]))


bt.obs$Scalar_max<-c(((bt.obs[1,2]+bt.obs[1,3])/bt.obs[1,2]),
                     ((bt.obs[2,2]+bt.obs[2,3])/bt.obs[1,2]),
                     ((bt.obs[3,2]+bt.obs[3,3])/bt.obs[1,2]),
                     ((bt.obs[4,2]+bt.obs[4,3])/bt.obs[1,2])) 


bt.obs$Scalar_min<-c(((bt.obs[1,2]-bt.obs[1,3])/bt.obs[1,2]),
                     ((bt.obs[2,2]-bt.obs[2,3])/bt.obs[1,2]),
                     ((bt.obs[3,2]-bt.obs[3,3])/bt.obs[1,2]),
                     ((bt.obs[4,2]-bt.obs[4,3])/bt.obs[1,2])) 

bt.obs$k<-bt.obs$Mean_Bt_fin*40

ggplot(bt.obs, aes(x=Treatment, y=Scalar, fill=Treatment))+
  theme_bw()+
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=15))+
  scale_fill_manual(values=c("red", "pink", "orange", "blue")) +
  ylab("Mean +/- SEM observed B. truncatus")+
  geom_bar(position=position_dodge(), stat="identity", width = .7) +
  geom_errorbar(aes(ymin=Scalar_min,
                    ymax=Scalar_max),
                width=.2, position=position_dodge(.7))

#VECTOR TO BE USED IN MODEL
phi_Nqs<-bt.obs$Scalar

#Carrying capacity of snail population ##############
fin_num.c<-data.frame('Treatment' = dat$treat[dat$chlor==1], 
                      'Bt_fin' = dat$bt_liv_fin[dat$chlor==1],
                      'Bg_fin' = dat$bg_liv_fin[dat$chlor==1],
                      'Hs_fin' = dat$hs_liv_fin[dat$chlor==1],
                      'All_fin' = 0)

fin_num.c$Treatment<-factor(fin_num.c$Treatment, levels = c('C 1x', 'C 2x',
                                                            'C + A', 'F + C',
                                                            'A + C + F'))

  for(i in 1:25){
    fin_num.c[i,5] = sum(fin_num.c[i,2], fin_num.c[i,3], fin_num.c[i,4])
  }

fin_num.c.agg<-aggregate.data.frame(fin_num.c, by = list(fin_num.c[,1]), FUN = mean)