#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

peri<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/R_use.csv")
  for(i in 1:nrow(peri)){
    peri[i,43]=paste(peri[i,2], peri[i,3], peri[i,4], sep="_")
    } #Add unique treatment code to data frame
    treats<-unique(peri[,43])
  peri<-peri[,-c(2:17,23:42)]  

snail<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/snail_repro.csv")
  
 peri.rep<-merge(peri, snail, by='tank')
 peri.rep<-peri.rep[,-c(8:10)]
 colnames(peri.rep)[7]<-"Treatment"
   
 
 peri.rep$Treatment<-as.character(peri.rep$Treatment)
 
   peri.rep$Treatment[peri.rep$Treatment=="0_0_0"]<- "Control"
   peri.rep$Treatment[peri.rep$Treatment=="1_0_0"]<- "Atrazine"
   peri.rep$Treatment[peri.rep$Treatment=="0_1_0"]<- "ChlorP"
   peri.rep$Treatment[peri.rep$Treatment=="0_0_1"]<- "Fertilizer"
   peri.rep$Treatment[peri.rep$Treatment=="1_1_0"]<- "Atrazine_ChlorP"
   peri.rep$Treatment[peri.rep$Treatment=="1_0_1"]<- "Atrazine_Fertilizer"
   peri.rep$Treatment[peri.rep$Treatment=="0_1_1"]<- "ChlorP_Fertilizer"
   peri.rep$Treatment[peri.rep$Treatment=="1_1_1"]<- "All_Three"
   
    peri.rep$Treatment<-as.factor(peri.rep$Treatment)
    
    peri.rep$Treatment<- factor(peri.rep$Treatment, levels=c("Control","Atrazine","ChlorP",
                                                           "Fertilizer","Atrazine_ChlorP",
                                                           "Atrazine_Fertilizer",
                                                           "ChlorP_Fertilizer","All_Three"))
    
  cbPalette <- c("#999999", "yellow", "red", "green", "orange", "darkgreen", "pink", "blue")
#Baseline periphyton chlorA levels with week1 egg counts
 ggplot(peri.rep, aes(x=peri0, y=TBtegg1, color=Treatment))+
   theme_bw()+
   scale_color_manual(values=cbPalette)+
   geom_point(size=4)
#Week1 periphyton chlorA levels with week2 egg counts
 ggplot(peri.rep, aes(x=peri1, y=TBtegg2, color=Treatment))+
   theme_bw()+
   scale_color_manual(values=cbPalette)+
   geom_point(size=4)
#Week2 periphyton chlorA levels with week3 egg counts
 ggplot(peri.rep, aes(x=peri2, y=TBtegg3, color=Treatment))+
   theme_bw()+
   scale_color_manual(values=cbPalette)+
   geom_point(size=4)
 
peri.rep