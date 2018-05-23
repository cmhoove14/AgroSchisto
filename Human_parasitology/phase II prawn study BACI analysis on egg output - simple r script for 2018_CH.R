
#Phase II prawn study interim analysis
#date last updated: May 2018
#author: S. Sokolow

library(glmmADMB)
library(lme4)
library(ggplot2)

# read data
SH.df=read.csv("~/RemaisWork/Schisto/Stanford/Human_Parasitology_data/SH egg output stacked.csv", header=T, na.strings=c("NA","Null")) 

SM.df=read.csv('SM worm burden stacked removed extraneous.csv', header=T, na.strings=c("NA","Null")) 

#change cluster to factor:

SM.df=transform(SM.df, Cluster=as.factor(Cluster))
SH.df=transform(SH.df, Cluster=as.factor(Cluster))

# remove na's
SH.df=na.omit(SH.df)
SM.df=na.omit(SM.df)

# plot effects overall, by cluster, by village nested in cluster:

#all effects:
plot1= ggplot(data=SH.df,  aes(Label, log10(Data+1), colour=treatment))+stat_summary(fun.data = "mean_cl_boot")+stat_summary(fun.data = "mean_cl_boot", geom = "line", aes(group = treatment), linetype = "dotted") +theme_bw()+theme(legend.justification=c(1,1), legend.position="top")+theme(axis.title.x = element_blank(), panel.grid = element_blank())+
  ylab("log egg counts")
#all effects by cluster:
plot2= ggplot(data=SH.df,  aes(Label, log10(Data+1), colour=treatment))+stat_summary(fun.data = "mean_cl_boot")+stat_summary(fun.data = "mean_cl_boot", geom = "line", aes(group = treatment), linetype = "dotted") +theme_bw()+theme(legend.justification=c(1,1), legend.position="top")+theme(axis.title.x = element_blank(), panel.grid = element_blank())+
  ylab("log egg counts")+facet_grid(facets=~Cluster)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
#all effects by treatment and village nested in cluster:
plot3= ggplot(data=SH.df,  aes(Label, log10(Data+1), colour=treatment))+stat_summary(fun.data = "mean_cl_boot",aes(group = Village))+stat_summary(fun.data = "mean_cl_boot", geom = "line", aes(group = Village), linetype = "dotted") +theme_bw()+theme(legend.justification=c(1,1), legend.position="top")+theme(axis.title.x = element_blank(), panel.grid = element_blank())+
  ylab("log egg counts")+facet_grid(facets=~Cluster)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

#repeat above replacing SH for SM with plot 1b = same as plot 1 but for SM

# print plots for SH and SM on same screen
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(plot1, vp = vplayout(1, 1))
print(plot1b, vp = vplayout(1, 2))



#neg binomial BACI analysis on egg output w/ individual and village random effects
nbinom.SM<-glmmadmb(Data~(1|ID)+(1|Village)+Label*treatment+sex, data = SM.df, family="nbinom")
nbinom.SH<-glmmadmb(Data~(1|ID)+(1|Village)+Label*treatment+sex, data = SH.df, family="nbinom")

#neg binomial BACI analysis on egg output w/ individual and cluster random effects
nbinom2.SM<-glmmadmb(Data~(1|ID)+(1|Cluster)+Label*treatment+sex, data = SM.df, family="nbinom")
nbinom2.SH<-glmmadmb(Data~(1|ID)+(1|Cluster)+Label*treatment+sex, data = SH.df, family="nbinom")

# neg binomial BACI analysis with all 3 (cluster/village/ID) random effects does not converge

