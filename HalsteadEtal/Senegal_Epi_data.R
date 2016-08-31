
require(ggplot2)

#Data and plot for two villages with 4 measured burdens ###############################
lampsar<-as.data.frame(matrix(ncol=4, nrow=8))
  colnames(lampsar)<-c("Site","Date","eggs_10mL","st.err")

lampsar$Site<-c(rep(c("Lampsar1"), times=4), rep(c("Lampsar2"), times=4))

lampsar$Date<-c(rep(c("2011-7-1","2012-7-1","2013-2-1","2013-9-1"), times=2))
  lampsar$Date<-as.Date(lampsar$Date)

lampsar$eggs_10mL<-c(32,1.91,82.5,11.4,6.5,1.55,161.9,17.8)

lampsar$st.err<-c(8.9,1.19,23.8,3,2.9,.86,68,3.23)

cols<-c("red", "blue")

ggplot(lampsar, aes(x=Date, y=eggs_10mL, group=Site, color=Site)) +
  theme_bw()+
  scale_color_grey() +
  geom_line(position=position_dodge(.25), size=1) +
  geom_point(position=position_dodge(.25), size=3.5) +
  geom_errorbar(aes(ymin=eggs_10mL-st.err,
                    ymax=eggs_10mL+st.err),
                width=.2, position=position_dodge(.25)) +
  ggtitle("Longitudinal mean egg burden")  

#Plot including PZQ treatment assuming 100% coverage and efficacy ############################
lampsar2<-as.data.frame(matrix(ncol=4, nrow=14))
  colnames(lampsar2)<-c("Site","Date","eggs_10mL","st.err")

lampsar2$Site<-c(rep(c("Lampsar1"), times=7), rep(c("Lampsar2"), times=7))

lampsar2$Date<-c(rep(c("2011-7-1","2012-2-1","2012-7-1","2012-8-1","2013-2-1","2013-3-1","2013-9-1"), times=2))
  lampsar2$Date<-as.Date(lampsar2$Date)

lampsar2$eggs_10mL<-c(32,0,1.91,0,82.5,0,11.4,6.5,0,1.55,0,161.9,0,17.8)

lampsar2$st.err<-c(8.9,2,1.19,2,23.8,2,3,2.9,2,.86,2,68,2,3.23)

cols<-c("red", "blue")

ggplot(lampsar2, aes(x=Date, y=eggs_10mL, group=Site, color=Site)) +
  theme_bw()+
  scale_color_grey() +
  geom_line(position=position_dodge(.25), size=1) +
  geom_point(position=position_dodge(.25), size=3.5) +
  geom_errorbar(aes(ymin=eggs_10mL-st.err,
                    ymax=eggs_10mL+st.err),
                width=.2, position=position_dodge(.25)) +
  ggtitle("Longitudinal mean egg burden") 

#Data and plot for 5 villages with two measured burdens #########################################
burd<-as.data.frame(matrix(ncol=3,nrow=10))
  colnames(burd)<-c("Site", "Date", "Mean_Worm_Burden")
  
burd$Site<-c(rep(c("Lampsar1"), times=2), rep(c("Lampsar2"), times=2), rep(c("Mbarigo1"), times=2),
             rep(c("Mbarigo2"), times=2), rep(c("Ndiol Maure"), times=2))  

burd$Date<-c(rep(c("2011-7-1", "2012-7-1")))
  burd$Date<-as.Date(burd$Date)

burd$Mean_Worm_Burden<-c(25,1.12,6.5,1.5,17.9,6.2,35,12.1,29.2,8.5)

cols2<-c("red", "blue", "darkgreen", "orange", "purple")

ggplot(burd, aes(x=Date, y=Mean_Worm_Burden, group=Site, color=Site)) +
  theme_bw()+
  scale_color_grey() +
  geom_line(position=position_dodge(.25), size=1) +
  geom_point(position=position_dodge(.25), size=3.5) +
  ggtitle("Longitudinal mean worm burden")  

#Plot with PZQ treatment assuming 100% coverage and efficacy (i.e. W drops to 0) ################################
burd2<-as.data.frame(matrix(ncol=3,nrow=15))
colnames(burd2)<-c("Site", "Date", "Mean_Worm_Burden")

burd2$Site<-c(rep(c("Lampsar1"), times=3), rep(c("Lampsar2"), times=3), rep(c("Mbarigo1"), times=3),
             rep(c("Mbarigo2"), times=3), rep(c("Ndiol Maure"), times=3))  

burd2$Date<-c(rep(c("2011-7-1","2012-2-1", "2012-7-1")))
burd2$Date<-as.Date(burd2$Date)

burd2$Mean_Worm_Burden<-c(25,0, 1.12,6.5,0,1.5,17.9,0,6.2,35,0,12.1,29.2,0,8.5)

cols2<-c("red", "blue", "darkgreen", "orange", "purple")

ggplot(burd2, aes(x=Date, y=Mean_Worm_Burden, group=Site, color=Site)) +
  theme_bw()+
  #scale_color_manual(cols2) +
  geom_line(position=position_dodge(.25), size=1) +
  geom_point(position=position_dodge(.25), size=3.5) +
  ggtitle("Longitudinal mean worm burden") 

#All Epi data put together, with "Rainy" seasons included ############################

inf<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/human_infection_R.csv")
  inf$date<-as.Date(inf$date, format="%m/%d/%Y")
  inf$Site<-as.factor(inf$Site)
  inf.0<-subset(inf, wrm_burd ==0)
  inf.pos<-subset(inf, wrm_burd !=0)
  
plot(inf.0$wrm_burd ~ inf.0$date, xlim=c(as.Date("2011-06-01"),as.Date("2013-9-1")), ylim=c(0,175),
     xlab="Date", ylab="Mean Worm Burden (W)", main="Known epi data",
     col="black", pch=16, cex=2)
  rect(xleft=as.Date("2011-06-01"), xright=as.Date("2011-10-31"), #2011 rainy season
       ybottom=0, ytop=200, col="blue", density=20)  
  rect(xleft=as.Date("2012-06-01"), xright=as.Date("2012-10-31"), #2012 rainy season
       ybottom=0, ytop=200, col="blue", density=20)  
  rect(xleft=as.Date("2013-06-01"), xright=as.Date("2013-10-31"), #2013 rainy season
       ybottom=0, ytop=200, col="blue", density=20)
  points(inf.pos$date[inf.pos$Site=="Lampsar2"], inf.pos$wrm_burd[inf.pos$Site=="Lampsar2"],
         pch=17, col="red", cex=2) #Add Lampsar2 measured worm burdens
  points(inf.pos$date[inf.pos$Site=="Mbarigo1"], inf.pos$wrm_burd[inf.pos$Site=="Mbarigo1"],
         pch=17, col="blue", cex=2) #Add Mbarigo1 measured worm burdens
  points(inf.pos$date[inf.pos$Site=="Mbarigo2"], inf.pos$wrm_burd[inf.pos$Site=="Mbarigo2"],
         pch=17, col="darkgreen", cex=2) #Add Mbarigo2 measured worm burdens
  points(inf.pos$date[inf.pos$Site=="Ndiol_Maure"], inf.pos$wrm_burd[inf.pos$Site=="Ndiol_Maure"],
         pch=17, col="orange", cex=2) #Add Ndiol_Maure measured worm burdens
  segments(x0=as.Date("2012-2-7"), y0=0,
           x1=as.Date("2012-6-22"), y1=8,
           lty=2, lwd=2, col="red") #Transmission from Feb '12 to July '12
  segments(x0=as.Date("2012-7-5"), y0=0,
           x1=as.Date("2013-2-7"), y1=161,
           lty=2, lwd=2, col="red") #Lampsar2 transmission from July '12 to Feb '13
  segments(x0=as.Date("2013-2-7"), y0=0,
           x1=as.Date("2013-8-25"), y1=17.5,
           lty=2, lwd=2, col="red") #Lampsar2 transmission from July '12 to Feb '13
  legend(x=as.Date("2011-11-21"), y=130, legend=c("PZQ", "Lampsar2", "Mbarigo1", "Mbarigo2", "Ndiol_Maure"),
         pch=c(16,17,17,17,17), col=c("black", "red", "blue", "darkgreen", "orange"))
  legend(x=as.Date("2012-11-20"), y=80, legend=c("transmission"), lty=2, lwd=2, col="red")
  