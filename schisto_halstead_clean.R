#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Cleaned version of schisto halstead modeling
  #Step 1: estimate agrochemical sensitive parameters from mesocosm data
    #Step 1.1: Scalar of snail carrying capacity (bottom-up effect)
    #Step 1.2: Scalar of snail reproduction rate (bottom-up effect)
    #Step 1.3: Daily prawn mortality rate (top-down effect)
    #Step 1.4: Excess daily prawn mortality rate from chlorP exposure (top-down effect)
  #Step 2: R0 runs with mesocosm data
    #Step 2.1: 
  #Step 3: R0 runs with additional outside data on chlorP and atrazine dose-response
    #Step 3.1: Predator mortality from Halstead et al 2015
    #Step 3.2: Enhanced growth rate from Atrazine; dose-response from Rohr analysis of Baxter data

#Load libraries and small functions #############
require(deSolve)
require(drc)
require(ggplot2)
require(reshape)
require(lattice)
require(ggthemes)
require(gridExtra)
st.er <- function(x) {
  sd(x)/sqrt(length(x))
} #Function to calculate standard error of the mean
cbPalette <- c("#999999", "yellow", "red", "green", "orange", "darkgreen", "pink", "blue")

#Step 1.1 Estimate change in snail carrying capacity by relative increases in final snail numbers ###########
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
#Step 1.2 Estimate change in snail reproduction rate by peak growth rates in tanks ###############
derp<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/snail_repro.csv")
    
  for(i in 1:nrow(derp)){ 
      derp[i,5]=paste(derp[i,2], derp[i,3], derp[i,4], sep="_")
  } #Add unique treatment code to data frame
    treats<-unique(derp[,5])
    colnames(derp)[5]<-"Treatment"
    
    varbs<-colnames(derp)
    
  bt.chlor<-subset(derp, chlor == 1, 
                   select = c(tank,Treatment, TBthatch1, TBthatch2, TBthatch3, TBthatch4, TBthatch8))  
  
  bt.chlor<-reshape(bt.chlor, varying = c("TBthatch1", "TBthatch2", "TBthatch3", "TBthatch4", "TBthatch8"),
                    v.names="measure",
                    timevar="variable",
                    times= c("TBthatch1", "TBthatch2", "TBthatch3", "TBthatch4", "TBthatch8"),
                    new.row.names=1:3600,
                    direction="long")
  
  bt.chlor$Treatment[bt.chlor$Treatment=="0_1_0"]<- "ChlorP"
  bt.chlor$Treatment[bt.chlor$Treatment=="1_1_0"]<- "Atrazine_ChlorP"
  bt.chlor$Treatment[bt.chlor$Treatment=="0_1_1"]<- "ChlorP_Fertilizer"
  bt.chlor$Treatment[bt.chlor$Treatment=="1_1_1"]<- "All_Three"
  
  bt.chlor$Treatment<-factor(bt.chlor$Treatment, levels=c("ChlorP", "ChlorP_Fertilizer", "Atrazine_ChlorP",
                                                          "All_Three"))
  
  bt.chlor$variable[bt.chlor$variable=="TBthatch1"]<-"Week 1"
  bt.chlor$variable[bt.chlor$variable=="TBthatch2"]<-"Week 2"
  bt.chlor$variable[bt.chlor$variable=="TBthatch3"]<-"Week 3"
  bt.chlor$variable[bt.chlor$variable=="TBthatch4"]<-"Week 4"
  bt.chlor$variable[bt.chlor$variable=="TBthatch8"]<-"Week 8"
  colnames(bt.chlor)[3]<-"Time"
  
  ggplot(bt.chlor, aes(x=Time, y=measure, group=tank, color=Treatment)) +
    theme_bw()+
    scale_color_manual(values=c("red", "pink", "orange", "blue")) +
    geom_line(position=position_dodge(.25), size=1) +
    geom_point(position=position_dodge(.25), size=3.5) +
    ggtitle("B. truncatus hatchlings sampled over time")  
  
  #Create data frame to fill with weakly growth rates
  
  bt.growth<-subset(derp, chlor == 1, 
                    select = c(tank,Treatment, TBthatch1, TBthatch2, TBthatch3, TBthatch4))  
  #Calculate peak mean daily per capita growth rate from week 1 - 2; 1 - 3; 1 - 4 
    #growth rate =(#births/population)/time
      #Assume population = number of B. truncatus added at beginning of mesocosm and that hatchlings do not mature
      #to sexual maturity by 4 weeks
  for(i in 1:25){
    bt.growth[i,7] = (bt.growth[i,4]/11) / 7
    bt.growth[i,8] = (bt.growth[i,5]/11) / 14
    bt.growth[i,9] = (bt.growth[i,6]/11) / 21
  }
  
  for(i in 1:25){
    bt.growth[i,10] = max(c(bt.growth[i,7], bt.growth[i,8], bt.growth[i,9]))
  }
  #Aggregate to mean for each treatment and plot
  ag1<-aggregate.data.frame(bt.growth, by=list(bt.growth[,2]), FUN=mean) #calculate means of treatment groups
  ag1<-ag1[,-c(2:7)] #remove unneeded variables 
  colnames(ag1)<-c("Treatment", "Week2.mean", "Week3.mean", "Week4.mean", "Peak.mean")
  
  ag2<-aggregate.data.frame(bt.growth, by=list(bt.growth[,2]), FUN=st.er) #calculate st.error of treatment groups
  ag2<-ag2[,-c(2:7)] #remove unneeded variables 
  colnames(ag2)<-c("Treatment", "Week2.ster", "Week3.ster", "Week4.ster", "Peak.ster")  
  
  growth.ag<-merge(ag1, ag2, by=c("Treatment"))
    growth.ag$Treatment[growth.ag$Treatment=="0_1_0"]<- "ChlorP"
    growth.ag$Treatment[growth.ag$Treatment=="1_1_0"]<- "Atrazine_ChlorP"
    growth.ag$Treatment[growth.ag$Treatment=="0_1_1"]<- "ChlorP_Fertilizer"
    growth.ag$Treatment[growth.ag$Treatment=="1_1_1"]<- "All_Three"
    
    growth.ag$Treatment<-factor(growth.ag$Treatment, levels=c("ChlorP", "ChlorP_Fertilizer", "Atrazine_ChlorP",
                                                            "All_Three"))
    
  ggplot(growth.ag, aes(x=Treatment, y=Peak.mean, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=c("red", "pink", "orange", "blue")) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=Peak.mean-Peak.ster,
                      ymax=Peak.mean+Peak.ster),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Peak growth rate (from hatchling samplers)")
  
  growth.ag$Scalar<-c((growth.ag[1,5]/growth.ag[1,5]),
                   (growth.ag[2,5]/growth.ag[1,5]),
                   (growth.ag[3,5]/growth.ag[1,5]),
                   (growth.ag[4,5]/growth.ag[1,5]))
  
  
  growth.ag$Scalar_max<-c(((growth.ag[1,5]+growth.ag[1,9])/growth.ag[1,5]),
                       ((growth.ag[2,5]+growth.ag[2,9])/growth.ag[1,5]),
                       ((growth.ag[3,5]+growth.ag[3,9])/growth.ag[1,5]),
                       ((growth.ag[4,5]+growth.ag[4,9])/growth.ag[1,5])) 
  
  
  growth.ag$Scalar_min<-c(((growth.ag[1,5]-growth.ag[1,9])/growth.ag[1,5]),
                       ((growth.ag[2,5]-growth.ag[2,9])/growth.ag[1,5]),
                       ((growth.ag[3,5]-growth.ag[3,9])/growth.ag[1,5]),
                       ((growth.ag[4,5]-growth.ag[4,9])/growth.ag[1,5])) 
  
  ggplot(growth.ag, aes(x=Treatment, y=Scalar, fill=Treatment)) +
    theme_bw()+
    scale_fill_manual(values=c("red", "pink", "orange", "blue")) +
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=Scalar_min,
                      ymax=Scalar_max),
                  width=.2, position=position_dodge(.7)) +
    ggtitle("Peak growth rate SCALAR (from hatchling samplers)")
  
  #VECTOR TO BE USED IN MODEL
    f_Nqs<-growth.ag$Scalar
    
#Step 1.2.2 Estimate change in snail reproduction rate from beginning and final counts in treatment groups ################
  bt.fin.ch = sum(dat$bt_liv_fin[dat$chlor==1]) #Final number B. truncatus in ChlorP tanks
  bt.start.ch = length(dat$tank[dat$chlor==1])*27  #Starting number of B. truncatus
  
  bt.fin.ch0 = sum(dat$bt_liv_fin[dat$chlor==0]) #Final number B. truncatus in ChlorP tanks
  bt.start.ch0 = length(dat$tank[dat$chlor==0])*27  #Starting number of B. truncatus

  
  ch.fn<-((sum(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==0]) - #Final number B. truncatus
         length(dat$tank[dat$chlor==1 & dat$atra==0 & dat$fert==0])*27) / #Starting number of B. truncatus
         (length(dat$tank[dat$chlor==1 & dat$atra==0 & dat$fert==0])*27)) / #births / starting pop gives per capita birthrate
         84 # / 84 days gives per daily per capita birthrate over study period
    
  ch.fe.fn<-((sum(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 & dat$fert==1]) - #Final number B. truncatus
              length(dat$tank[dat$chlor==1 & dat$atra==0 & dat$fert==1])*27) / #Starting number of B. truncatus
             (length(dat$tank[dat$chlor==1 & dat$atra==0 & dat$fert==1])*27)) / #births / starting pop gives per capita birthrate
               84 # / 84 days gives per daily per capita birthrate over study period  
  
  ch.at.fn<-((sum(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==0]) - #Final number B. truncatus
              length(dat$tank[dat$chlor==1 & dat$atra==1 & dat$fert==0])*27) / #Starting number of B. truncatus
             (length(dat$tank[dat$chlor==1 & dat$atra==1 & dat$fert==0])*27)) / #births / starting pop gives per capita birthrate
                84 # / 84 days gives per daily per capita birthrate over study period 
  
  ch.at.fe.fn<-((sum(dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 & dat$fert==1]) - #Final number B. truncatus
                length(dat$tank[dat$chlor==1 & dat$atra==1 & dat$fert==1])*27) / #Starting number of B. truncatus
               (length(dat$tank[dat$chlor==1 & dat$atra==1 & dat$fert==1])*27)) / #births / starting pop gives per capita birthrate
                84 # / 84 days gives per daily per capita birthrate over study period 
  
  fNqs.ch<-c(ch.fn / ch.fn,
             ch.fe.fn / ch.fn,
             ch.at.fn / ch.fn,
             ch.at.fe.fn / ch.fn)
  
  plot.fnchs<-data.frame('Treatment' = factor(c('Ch', 'Ch:Fe', 'Ch:At', 'At:Ch:Fe'), 
                                             levels = c('Ch', 'Ch:Fe', 'Ch:At', 'At:Ch:Fe')),
                         'f_N_Scalar' = fNqs.ch)
  
    ggplot(plot.fnchs, aes(x = Treatment, y = f_N_Scalar, fill = Treatment)) +
      theme_bw()+
      scale_fill_manual(values=c("red", "pink", "orange", "blue")) +
      geom_bar(position=position_dodge(), stat="identity", width=.7) +
      ggtitle("f_Nq scalar from 84 day total reproduction")
      
#Step 1.3 Estimate daily prawn mortality in chlorP-free tanks ###################
  p.0<-3*length(dat$tank[dat$chlor==0]) #Total starting number of prawns in chloP-free tanks (3 in each)
  p.1<-sum(dat$p.all_24[dat$chlor==0]) #Total surviving prawns in chlorP-free tanks after day 1
  p.84<-sum(dat$p.all_fin[dat$chlor==0]) #Total surviving prawns in chlorP-free tanks at end of experiment (12 weeks)
  muP<-(p.0-p.1)/p.0 #Per capita daily death rate by 1-day endpoint = number of deaths/ starting population
  
#Step 1.4 Estimate daily prawn mortality in chlorP tanks ################### 
  p.0.q<-3*length(dat$tank[dat$chlor==1]) #Total starting number of prawns in chloP tanks (3 in each)
  p.1.q<-sum(dat$p.all_24[dat$chlor==1]) #Total surviving prawns in chlorP tanks after day 1
  p.84.q<-sum(dat$p.all_fin[dat$chlor==1]) #Total surviving prawns in chlorP tanks at end of experiment (12 weeks)
  muPq<-(p.0.q-p.1.q)/p.0.q #Per capita death rate = number of deaths / starting population
#Step 2 Read in parameter values and R0 function ############
  parameters=c( #Excluding beta, Phi_Nq, f_Nq, and muPq which will be read into the R0 function
    ##standard snail parameters 
      f_N=0.16, # recruitment rate (from sokolow et al)
      phi_N=(1-1/(80*0.16))/10000, # carrying capacity from sokolow et al
      z=0.5, #Proportion of exposed snails that reproduce from sokolow et al
      mu_N=1/80, #Mortality rate from Sokolow et al
      sigma=1/50, #Transition rate from exposed to infected from sokolow et al
      mu_I=1/20, #additional snail death due to infection from sokolow et al

    #prawn parameters
      alpha=0.003, #attack rate
      Th=0.1,#~Prawn predation limit
      f_P=0.234/2, #prawn birth rate from Cervantes-Santiago Aquaculture 2010 paper (/2 for 1:1 female-male ratio)
      phi_P=1/(40*3),  #prawn carrying capacity
      mu_P= muP, #observed daily 24 hr prawn mortality rate in chlorP-free tanks of mesocosm#0.00137, previous rate
      
    #Adult Worm, Miracidia and Circariae Parameters
      lamda=0.00004, #probability of snail shedding a cercariae that infects a human host and survives to reproduction
      mu_W=1/(3*365), # death rate of adult worms
      m=0.5, #miracidial shedding rate per reproductive female divided by miracidial mortality; from sokolow et al
      
    #Human parameters
      H=300, #number of humans
      mu_H=1/(60*365) #Assumes 60 year lifespan
    )
  
  get_Ro<-function(muPq = 0, phi_Nq = 1, beta, f_Nq = 1) #variable parameters to be manipulated
  { 
    f_N<-parameters["f_N"]
    phi_N<-parameters["phi_N"]
    z<-parameters["z"]
    mu_N<-parameters["mu_N"]
    sigma<-parameters["sigma"]
    mu_I<-parameters["mu_I"]
    alpha<-parameters["alpha"]
    Th<-parameters["Th"]
    f_P<-parameters["f_P"]
    phi_P<-parameters["phi_P"]
    mu_P<-parameters["mu_P"]
    lamda<-parameters["lamda"]
    mu_W<-parameters["mu_W"]
    m<-parameters["m"]
    H<-parameters["H"]
    mu_H<-parameters["mu_H"]
    
    P_eq<-(1-((muPq+mu_P)/f_P))/phi_P #Equilibrium estimate of P given prawn predator parameters
    if(P_eq<0){
      P_eq=0
    }
    #Equilibrium estimate of N given snail parameters
      #Shorthand values to use in N_eq expression
      a= -(alpha*Th*f_N*f_Nq*phi_N*phi_Nq^-1)
      b= f_N*f_Nq*alpha*Th - mu_N*alpha*Th - f_N*f_Nq*phi_N*phi_Nq^-1
      c= f_N*f_Nq - mu_N - alpha*P_eq
    
      if((b^2-4*a*c)<0){ #If prawn population sufficient to eliminate snails, N_eq=0
        N_eq=0
      } else {
        N_eq <- (-b - sqrt(b^2-4*a*c)) / (2*a) #Function to solve quadratic expression for N_eq
      }
    
    
    pred<-(alpha*P_eq)/(1+(alpha*N_eq*Th))#death rate of snails due to predators given equilibrium estimates of P and N
    
    T1<-0.5*beta*m*H*N_eq
    T2<-lamda*sigma
    T3<- (mu_W+mu_H)*(mu_N+pred+sigma)*(mu_N+pred+mu_I)
    
    Ro_est <- sqrt((T1*T2)/T3)
    
    print(N_eq)
    print(P_eq)
    Ro_est 
    
  }  
  
p.dead = parameters["f_P"] - parameters["mu_P"] #muPq value at which predator population is eliminated 
  
#Step 2.1 Get R0 estimates across agrochemical treatment groups #####################
  #Beta estimates from model fitting to epi data
    beta_0= 2.8372e-05
    beta_up=4e-05 #NEED TO UPDATE FROM ARATHI< THIS IS JUST A PLACEHOLDER
    beta_lo=1e-05

  #village R0        
    R0_vil = get_Ro(beta = beta_0, muPq = p.dead) 
    R0_up = get_Ro(beta = beta_up, muPq = p.dead) 
    R0_lo = get_Ro(beta = beta_lo, muPq = p.dead) 
  
  #village R0 with prawns
    R0_vil.p<-get_Ro(beta=beta_0)
    R0_up.p<-get_Ro(beta=beta_up)
    R0_lo.p<-get_Ro(beta=beta_lo)
  
  
  #atrazine only R0
    R0_atra0 = get_Ro(phi_Nq = phi_Nqs[3], 
                      #f_Nq = f_Nqs[3], 
                      #f_Nq = fNqs.ch[3],
                      beta = beta_0)
    R0_atra_up = get_Ro(phi_Nq = phi_Nqs[3], 
                        #f_Nq = f_Nqs[3], 
                        #f_Nq = fNqs.ch[3],
                        beta = beta_up)
    R0_atra_lo = get_Ro(phi_Nq = phi_Nqs[3], 
                        #f_Nq = f_Nqs[3], 
                        #f_Nq = fNqs.ch[3],
                        beta = beta_lo)
  
  #chlorpyrifos only R0
    R0_chlor0 = get_Ro(muPq = muPq, 
                       beta = beta_0)
    R0_chlor_up = get_Ro(muPq = muPq, 
                         beta = beta_up)
    R0_chlor_lo = get_Ro(muPq = muPq, 
                         beta = beta_lo)
  
  #fertilizer only R0
    R0_fert0 = get_Ro(phi_Nq = phi_Nqs[2], 
                      #f_Nq = f_Nqs[2], 
                      #f_Nq = fNqs.ch[2],
                      beta = beta_0)
    R0_fert_up = get_Ro(phi_Nq = phi_Nqs[2], 
                        #f_Nq = f_Nqs[2], 
                        #f_Nq = fNqs.ch[2],
                        beta = beta_up)
    R0_fert_lo = get_Ro(phi_Nq = phi_Nqs[2], 
                        #f_Nq = f_Nqs[2],
                        #f_Nq = fNqs.ch[2],
                        beta = beta_lo)  
  
  #atrazine+chlorpyrifos R0
    R0_atch0 = get_Ro(muPq = muPq, 
                      phi_Nq = phi_Nqs[3], 
                      #f_Nq = f_Nqs[3], 
                      #f_Nq = fNqs.ch[3],
                      beta = beta_0)
    R0_atch_up = get_Ro(muPq = muPq, 
                        phi_Nq = phi_Nqs[3], 
                        #f_Nq = f_Nqs[3], 
                        #f_Nq = fNqs.ch[3],
                        beta = beta_up)
    R0_atch_lo = get_Ro(muPq = muPq, 
                        phi_Nq = phi_Nqs[3], 
                        #f_Nq = f_Nqs[3], 
                        #f_Nq = fNqs.ch[3],
                        beta = beta_lo) 
  
  #atrazine+fertilizer R0
    R0_atfe0 = get_Ro(phi_Nq = phi_Nqs[4], 
                      #f_Nq = f_Nqs[4], 
                      #f_Nq = fNqs.ch[4],
                      beta = beta_0)
    R0_atfe_up = get_Ro(phi_Nq = phi_Nqs[4], 
                        #f_Nq = f_Nqs[4], 
                        #f_Nq = fNqs.ch[4],
                        beta = beta_up)
    R0_atfe_lo = get_Ro(phi_Nq = phi_Nqs[4], 
                        #f_Nq = f_Nqs[4], 
                        #f_Nq = fNqs.ch[4],
                        beta = beta_lo) 
  
  #chlorpyrifos+fertilizer R0
    R0_chfe0 = get_Ro(muPq = muPq, 
                      phi_Nq = phi_Nqs[2], 
                      #f_Nq = f_Nqs[2],
                      #f_Nq = fNqs.ch[2],
                      beta = beta_0)
    R0_chfe_up = get_Ro(muPq = muPq, 
                        phi_Nq = phi_Nqs[2], 
                        #f_Nq = f_Nqs[2], 
                        #f_Nq = fNqs.ch[2],
                        beta = beta_up)
    R0_chfe_lo = get_Ro(muPq = muPq, 
                        phi_Nq = phi_Nqs[2], 
                        #f_Nq = f_Nqs[2], 
                        #f_Nq = fNqs.ch[2],
                        beta = beta_lo)  
  
  #atrazine+chlorpyrifos+fertilizer R0
    R0_tre0 = get_Ro(muPq = muPq, 
                     phi_Nq = phi_Nqs[4], 
                     #f_Nq = f_Nqs[4], 
                     #f_Nq = fNqs.ch[4],
                     beta = beta_0)
    R0_tre_up = get_Ro(muPq = muPq, 
                       phi_Nq = phi_Nqs[4], 
                       #f_Nq = f_Nqs[4], 
                       #f_Nq = fNqs.ch[4],
                       beta = beta_up)
    R0_tre_lo = get_Ro(muPq = muPq, 
                       phi_Nq = phi_Nqs[4], 
                       #f_Nq = f_Nqs[4], 
                       #f_Nq = fNqs.ch[4],
                       beta = beta_lo) 
  
  r0s.3<-data.frame("Treatment"=c('At', 'Ch', 'Fe',  
                                  'At:Ch', 'At:Fe', 'Ch:Fe', 'At:Ch:Fe'),
                    'r0_0'=c(R0_atra0, R0_chlor0, R0_fert0, R0_atch0, R0_atfe0, R0_chfe0, R0_tre0),
                    'r0_up'=c(R0_atra_up, R0_chlor_up, R0_fert_up, R0_atch_up, R0_atfe_up, R0_chfe_up, R0_tre_up),
                    'r0_lo'=c(R0_atra_lo,  R0_chlor_lo, R0_fert_lo, R0_atch_lo, R0_atfe_lo, R0_chfe_lo, R0_tre_lo))  
  
  r0s.3$Treatment<-factor(r0s.3$Treatment, levels = c('Fe', 'At', 'At:Fe', #Bottom up effects only
                                                      'Ch', #Top-down effects only
                                                      'Ch:Fe', 'At:Ch',   'At:Ch:Fe')) #Both top-down and bottom-up effects
#Step 2.2 Plot R0 estimates across agrochemical treatment groups ##################
  gg1<-ggplot(r0s.3, aes(x=Treatment, y=r0_0))+
    #Theme formatting
      theme_bw()+
      theme(axis.title=element_text(size=14),
            axis.text=element_text(size=10))+
      scale_y_continuous(breaks=c(0, 1.0, 2.0, 3.0, 4.0, 5.0), limits=c(0,5.5))+
      xlab("")+
      ylab(expression('R'[0]))+
    #Village R0 lines
      geom_hline(aes(yintercept=R0_vil), colour='grey30', size=1, linetype=3)+
      #geom_hline(aes(yintercept=R0_up), colour='red', linetype=2, alpha=0.5, size=1)+
      #geom_hline(aes(yintercept=R0_lo), colour='red', linetype=2, alpha=0.5, size=1)+
      #geom_hline(aes(yintercept=R0_vil.p), colour='grey30', size=1, linetype=3)+
      #geom_hline(aes(yintercept=R0_up.p), colour='green2', linetype=2, alpha=0.5, size=1)+
      #geom_hline(aes(yintercept=R0_lo.p), colour='green2', linetype=2, alpha=0.5, size=1)+
    #Adding data from data frame
      geom_point(size=3) +
      geom_errorbar(aes(ymin=r0_lo,
                        ymax=r0_up),
                    width=.1, position=position_dodge(.7))+
      #Add labels to R0 lines
      geom_label(x=1.475, y=3.2, label="R ",  size=5, colour = 'grey30', fill='white', label.size = NA)+
      geom_text(x=1.525, y=3.09, label="0,f",  size=3, colour = 'grey30')+
      #geom_label(x=5.475, y=0.53, label='R ',  size=5, colour = 'grey30', fill='white', label.size = NA)+
      #geom_text(x=5.525, y=0.508, label="0,f",  size=3, colour = 'grey30')+
      #geom_text(x=5.585, y=0.50, label="p",  size=3, colour = 'grey30')+
    #Add treatment labels
      geom_segment(x=0.7, xend=3.3, y=5.25, yend=5.25, colour='grey40', lineend='square')+
      geom_text(x=2, y=5.45, label='bottom-up effects', size=5, colour='grey40')+
      geom_segment(x=3.7, xend=4.3, y=5.25, yend=5.25, colour='grey40', lineend='square')+
      geom_text(x=4, y=5.65, label='top-down', size=5, colour='grey40')+
      geom_text(x=4, y=5.45, label='effects', size=5, colour='grey40')+
      geom_segment(x=4.7, xend=7.3, y=5.25, yend=5.25, colour='grey40', lineend='square')+
      geom_text(x=6, y=5.45, label='bottom-up & top-down effects', size=5, colour='grey40')+
    #Add plot label
      geom_text(label='A', x=0.57, y=5.5, size=10)
  
  

#Step 3 incorporate agrochemical dose-response from outside sources ################
#Step 3.1.1 ChlorP dose-response from Halstead 2015 Chemosphere paper ################
#Observed 4-day mortality endpoints ###########
  ecotox_mod4<-data.frame('chem'=rep("Chlorpyrifos", 30),
                           'dose'=c(rep(0,5), rep(0.64,5), rep(3.2,5), 
                                    rep(6.4,5), rep(32,5), rep(64,5)),
                           'dead'=c(rep(0,5), rep(0,5), rep(0,5), 
                                    rep(0,6), rep(1,9)),
                           "per_dead"=c(rep(0,5), rep(0,5), rep(0,5), 
                                        rep(0,5), rep(4/5,5), rep(1,5)))
  
  plot(x = unique(log(ecotox_mod4$dose+1)), y = c(0,0,0,0,0.8,1),
       xlab = 'concentration (log+1)', ylab = "%mortality (#dead/5)", pch = 16)
  
  ecotox4<-glm(dead ~ dose, family=binomial(link="probit"),data=ecotox_mod4)
  summary(ecotox4)
  
  #Extrapolate response to constant gradient of Chlorpyrifos concentration
  p.ecotox4<-data.frame(dose=seq(from=0, to=150, by=0.1))
  p.ecotox4[, c('mortality', 'st.er')]<-predict(ecotox4, p.ecotox4, 
                                                 type = "response", se.fit=TRUE)
  
  plot(x=c(0,0.32,0.64,3.2,6.4,32,64), pch = 16,
       y=c(0,0,   0   ,0  ,0  ,4/5,5/5)/4, xlab = "ChlorP Concentration", ylab = "%Mortality")
    lines(x=p.ecotox4$dose, y=p.ecotox4$mortality/4, col='red')
    lines(x=p.ecotox4$dose, y=p.ecotox4$mortality/4+p.ecotox4$st.er/4, lty=2, col='red', cex=0.8)
    lines(x=p.ecotox4$dose, y=p.ecotox4$mortality/4-p.ecotox4$st.er/4, lty=2, col='red', cex=0.8)
    abline(a = 0.079, b = 0)
  
  for (i in 1:nrow(p.ecotox4)){
    p.ecotox4$R0[i] <- get_Ro(muPq = (p.ecotox4$mortality[i]/4), phi_Nq = 1, beta = beta_0)
  }
  
  for (i in 1:nrow(p.ecotox4)){
    p.ecotox4$R0_lo[i] <- get_Ro(muPq = ((p.ecotox4$mortality[i]/4 - p.ecotox4$st.er[i])/4), 
                                  phi_Nq = 1, beta = beta_0)
  }
  
  for (i in 1:nrow(p.ecotox4)){
    p.ecotox4$R0_up[i] <- get_Ro(muPq = ((p.ecotox4$mortality[i]/4 + p.ecotox4$st.er[i])/4), 
                                  phi_Nq = 1, beta = beta_0)
  }
  
  
#Observed 10-day mortality endpoints #########
  ecotox_10<-data.frame('chem'=rep("Chlorpyrifos", 30),
                          'dose'=c(rep(0,5), rep(0.64,5), rep(3.2,5), 
                                   rep(6.4,5), rep(32,5), rep(64,5)),
                          'dead'=c(rep(0,5), rep(0,5), rep(0,5), 
                                       rep(0,2), rep(1,13)),
                          "per_dead"=c(rep(0,5), rep(0,5), rep(0,5), 
                                       rep(3/5,5), rep(1,10)))
 


  ecotox10<-glm(dead ~ dose, family=binomial(link="probit"),data=ecotox_10)
    summary(ecotox10)
  
  #Extrapolate response to constant gradient of Chlorpyrifos concentration
  p.ecotox10<-data.frame('dose'=seq(from=0, to=150, by=0.01),
                         'log_dose'=log(seq(from=0, to=150, by=0.01)+1))
    p.ecotox10[, c('mortality', 'st.er')]<-predict(ecotox10, p.ecotox10, 
                                                type = "response", se.fit=TRUE)
    
    p.ecotox10$st.er = log(p.ecotox10$st.er+1)
  
  plot(x=log(c(0,0.32,0.64,3.2,6.4,32,64)+1), 
       y=c(0,0,   0   ,0  ,3/5,5/5,5/5)/10, xlab = "ChlorP Concentration", ylab = "%Mortality")
    lines(x=p.ecotox10$log_dose, y=p.ecotox10$mortality/10, col='red')
    lines(x=p.ecotox10$log_dose, y=p.ecotox10$mortality+p.ecotox10$st.er/10, lty=2, col='red', cex=0.8)
    lines(x=p.ecotox10$log_dose, y=p.ecotox10$mortality-p.ecotox10$st.er/10, lty=2, col='red', cex=0.8)
    parameters["mu_P"]=0 #let's just make all mortality due to mortality observed in ecotox study
    
  for (i in 1:nrow(p.ecotox10)){
    p.ecotox10$R0[i] <- get_Ro(muPq = (p.ecotox10$mortality[i]/10), phi_Nq = 1, beta = beta_0)
  }
  
  for (i in 1:nrow(p.ecotox10)){
    p.ecotox10$R0_lo[i] <- get_Ro(muPq = ((p.ecotox10$mortality[i]/10 + p.ecotox10$st.er[i])/10), 
                                 phi_Nq = 1, beta = beta_0)
  }
  
  for (i in 1:nrow(p.ecotox10)){
    p.ecotox10$R0_up[i] <- get_Ro(muPq = ((p.ecotox10$mortality[i]/10 - p.ecotox10$st.er[i])/10), 
                                 phi_Nq = 1, beta = beta_0)
  }
#10-day data with enhanced sample size and 99% mortality in last group ################  
  ecotox_mod10<-data.frame('dose'=c(rep(0,100), rep(0.64,100), rep(3.2,100), rep(6.4,100), rep(32,100), rep(64,100)),
                          'response'=c(rep(0,97), 1,1,1, rep(0,97), 1,1,1, rep(0,97), 1,1,1, 
                                       rep(0,40), rep(1,160),rep(1,100)))
  
  ecotox10_mod<-glm(response ~ dose, family=binomial(link="probit"),data=ecotox_mod10)
    summary(ecotox10_mod)
  
  #Extrapolate response to constant gradient of Chlorpyrifos concentration
  p.ecotox_mod10<-data.frame(dose=seq(from=0, to=150, by=1))
  p.ecotox_mod10[, c('mortality', 'st.er')]<-predict(ecotox10_mod, p.ecotox_mod10, 
                                               type = "response", se.fit=TRUE)
  
  plot(x=c(0,0.32,0.64,3.2,6.4,32,64), pch = 16,
       y=c(0,0,   0   ,0  ,60/100  ,99/100,99/100)/10, xlab = "ChlorP Concentration", ylab = "%Mortality")
    lines(x=p.ecotox_mod10$dose, y=p.ecotox_mod10$mortality/10, col='red')
    lines(x=p.ecotox_mod10$dose, y=p.ecotox_mod10$mortality/10+p.ecotox_mod10$st.er/10, lty=2, col='red', cex=0.8)
    lines(x=p.ecotox_mod10$dose, y=p.ecotox_mod10$mortality/10-p.ecotox_mod10$st.er/10, lty=2, col='red', cex=0.8)
    abline(a = 0.079, b = 0)
  plot(p.ecotox_mod10$dose, p.ecotox_mod10$mortality/10, type = 'l', xlab='concentration', ylab = 'daily mortality rate')
    abline(a = parameters['f_P'] - parameters['mu_P'], b = 0, lty=2, col='red')
  
    for (i in 1:nrow(p.ecotox_mod10)){
      p.ecotox_mod10$R0[i] <- get_Ro(muPq = (p.ecotox_mod10$mortality[i]/10), phi_Nq = 1, beta = beta_0)
    }
    
    for (i in 1:nrow(p.ecotox_mod10)){
      p.ecotox_mod10$R0_lo[i] <- get_Ro(muPq = ((p.ecotox_mod10$mortality[i]/10 + p.ecotox_mod10$st.er[i])/10), 
                                    phi_Nq = 1, beta = beta_0)
    }
    
    for (i in 1:nrow(p.ecotox_mod10)){
      p.ecotox_mod10$R0_up[i] <- get_Ro(muPq = ((p.ecotox_mod10$mortality[i]/10 - p.ecotox_mod10$st.er[i])/10), 
                                    phi_Nq = 1, beta = beta_0)
    }
  
  
#Step 3.1.2 Plot R0 response to modeled ChlorP mortality across concentrations ##############
  
  gg2<-ggplot(p.ecotox10, aes(x=dose, y=R0))+
    theme_bw()+
    theme(axis.text=element_text(size=10), axis.title=element_text(size=14))+
    #scale_y_continuous(breaks=c(0, 0.44, 0.5,1.0,1.07,1.5), limits=c(0,2))+
    scale_x_continuous(breaks=c(0,20,40,60,64), limits=c(0,70))+
    ylab(expression('R'[0]))+
    xlab(expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep = '')))+
    geom_line()+
    geom_line(aes(y=R0_lo), linetype=2)+
    geom_line(aes(y=R0_up), linetype=2)+
    geom_text(label='B', x=0, y=3, size=10)#+
  #geom_label(x=64, y=1.15, label="concentration", colour='grey40', size=4, fill='white',label.size = NA)+
  #geom_label(x=64, y=1.225, label="tested", colour='grey40', size=4, fill='white', label.size = NA)+
  #geom_segment(x=64, xend=64, y=-Inf, yend=1.08, linetype=3, colour='grey40')
  
#Step 3.2 Enhanced growth rate from Atrazine ##############
    
  r0.atra.chlor<-data.frame("Atra" = rep(c(0,1,10,100), 6),
                            "dose" = rep(c(0,0.64,3.2,6.4,32,64), each=4),
                            "f_Nq" = rep(c(1,1.2888,1.6535,2.3215), 6),#Growth rate scalars from fit to peak snails
                            "muPq" = rep(0,24),
                            "rate" = rep(0,24),
                            "R0" = rep(0,24))
  
  r0.atra.chlor$rate<-(predict(ecotox10, r0.atra.chlor, 
                                type = "response", se.fit=TRUE)$fit)/10 #fill mortality rate data from model
    
    
  for(i in 1:nrow(r0.atra.chlor)){
      r0.atra.chlor[i,6] = get_Ro(muPq = r0.atra.chlor[i,5],
                                  beta = beta_0,
                                  f_N = r0.atra.chlor[i,3])
    }
    
    r0.atra.chlor$Atra<-factor(r0.atra.chlor$Atra, levels=c(0, 1, 10, 100))
    r0.atra.chlor$dose<-factor(r0.atra.chlor$dose, levels=c(0, 0.64, 3.2, 6.4, 32, 64))
    
#Step 3.2.1 plot the heat map #################
    ggplot(r0.atra.chlor, aes(x=dose, y=Atra, fill=R0))+
      theme_bw()+
      geom_tile(color='white', size=0.1)+
      scale_fill_continuous(low='green', high='red')+
      coord_equal()+
      labs(x=expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep = '')), 
           #(x=expression(paste('Predator mortality rate (', mu[P][,][q], ')', sep = '')), axis label for predator mortality rate
           y=expression(paste('Atrazine concentration (', mu, 'g/L)', sep = '')))+
      theme(axis.ticks=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=14),
            legend.title=element_text(size=15), legend.text=element_text(size=12))+
      geom_text(label='C', x=0.75, y=7.25, size=10, alpha=.50)  