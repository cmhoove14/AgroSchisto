#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Cleaned version of schisto halstead modeling WITH HOLLING"S TYPE III FUNCIONTAL RESPONSE
#Step 1: estimate agrochemical sensitive parameters from mesocosm data
#Step 1.1: Scalar of snail carrying capacity (bottom-up effect)
#Step 1.2: Daily prawn mortality rate (top-down effect)
#Step 1.3: Daily prawn mortality rate from chlorP exposure (top-down effect)
#Step 2: R0 runs with mesocosm data
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
require(fitdistrplus)
require(rootSolve)

  st.er <- function(x) {
    sd(x)/sqrt(length(x))
  } #Function to calculate standard error of the mean
  
#Step 1.1 Estimate change in snail carrying capacity by relative increases in final snail numbers ###########
  dat<-read.csv("C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/R_use.csv")
  #Bootstrapping to obtain distribution of bottom-up effects 
    Fes<-dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 &dat$fert==1]
    Ats<-dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 &dat$fert==0]
    Fe.Ats<-dat$bt_liv_fin[dat$chlor==1 & dat$atra==1 &dat$fert==1]
    refs<-dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 &dat$fert==0]
    ref.mean<-mean(dat$bt_liv_fin[dat$chlor==1 & dat$atra==0 &dat$fert==0])
  
#Fertilizer bottom up effects  
  Fe.samples <- lapply(1:5000, function(i)
    (sample(Fes, replace = T)/ref.mean))
      Fe.mean<-sapply(Fe.samples, mean)
    #Check distribution 
        descdist(Fe.mean) #normal
      plot(density(Fe.mean))
      lines(density(rnorm(10000, 
                          mean = fitdistr(Fe.mean, 'normal')$estimate[1],
                          sd = fitdistr(Fe.mean, 'normal')$estimate[2])),
            col='red')
    #Extract distribution parameters
      Fe.mu = fitdistr(Fe.mean, 'normal')$estimate[1]
      Fe.sd = fitdistr(Fe.mean, 'normal')$estimate[2]
      
    #Generate parameter distribution to sample from
      Fe.use<-rnorm(5000, mean = Fe.mu, sd = Fe.sd)
      Fe.use<-Fe.use
        plot(density(Fe.use))
        
#Atrazine bottom up effects  
  At.samples <- lapply(1:5000, function(i)
    (sample(Ats, replace = T)/ref.mean))
      At.mean<-sapply(At.samples, mean)
    #Check distribution
        descdist(At.mean) #Normal
      plot(density(At.mean))
      lines(density(rnorm(10000, 
                          mean = fitdistr(At.mean, 'normal')$estimate[1],
                          sd = fitdistr(At.mean, 'normal')$estimate[2])),
            col='red')
    #Extract distribution parameters
      At.mu = fitdistr(At.mean, 'normal')$estimate[1]
      At.sd = fitdistr(At.mean, 'normal')$estimate[2]
      
    #Generate parameter distribution to sample from
      At.use<-rnorm(5000, mean = At.mu, sd = At.sd)
        plot(density(At.use))
      
#Atrazine and fertilizer combined bottom up effects  
  FeAt.samples <- lapply(1:5000, function(i)
    (sample(Fe.Ats, replace = T)/ref.mean))
      FeAt.mean<-sapply(FeAt.samples, mean)
    #Check distribution
        descdist(FeAt.mean) #Normal
      plot(density(FeAt.mean))
      lines(density(rnorm(10000, 
                          mean = fitdistr(FeAt.mean, 'normal')$estimate[1],
                          sd = fitdistr(FeAt.mean, 'normal')$estimate[2])),
            col='red')
    #Extract distribution parameters
      FeAt.mu = fitdistr(FeAt.mean, 'normal')$estimate[1]
      FeAt.sd = fitdistr(FeAt.mean, 'normal')$estimate[2]
      
    #Generate parameter distribution to sample from
      FeAt.use<-rnorm(5000, mean = FeAt.mu, sd = FeAt.sd)
      plot(density(FeAt.use))
      
#Bottom-up effects distributions for each treatment
  plot(density(Fe.use), col = 'green', xlab = 'Carrying capacity scalar', lwd=2, ylim = c(0,1.5), 
       main = 'Bottom-up parameter value distributions')
    lines(density(At.use), col = 'gold2', lwd=2)
    lines(density(FeAt.use), col = 'purple', lwd=2)
    abline(v=1, lty=2, lwd=2, col='red')
    legend('topright', legend = c('Fe', 'At', 'Fe:At'), lwd = 2, col = c('green', 'gold2', 'purple'))
      
#Step 1.2 Estimate daily prawn mortality in chlorP-free tanks ###################
  #Bootstrapping approach for distribution in MC approach
  p.24<-dat$p.all_24[dat$chlor==0] 
    p.24<-3-p.24 #convert to deaths
  #Bootstrap
    p.24.samples <- lapply(1:5000, function(i)
      (sample(p.24, replace = T)/3))
        muP.mean<-sapply(p.24.samples, mean)
        muP = mean(muP.mean)
  #Check distribution      
        descdist(muP.mean) #beta distribution since mortality rate restricted from 0 to 1
          plot(density(muP.mean))
          
          fitdist(muP.mean, distr = 'beta', start = list(shape1=1, shape2=1), method = 'mme')
          
          lines(density(rbeta(10000, 
                              shape1 = fitdist(muP.mean, distr = 'beta', 
                                               start = list(shape1=1, shape2=1), method = 'mme')$estimate[1],
                              shape2 = fitdist(muP.mean, distr = 'beta', 
                                               start = list(shape1=1, shape2=1), method = 'mme')$estimate[2])),
                col='red')
  #Extract distribution parameters
          muP.1 = fitdist(muP.mean, distr = 'beta', 
                            start = list(shape1=1, shape2=1), method = 'mme')$estimate[1]
          muP.2 = fitdist(muP.mean, distr = 'beta', 
                          start = list(shape1=1, shape2=1), method = 'mme')$estimate[2]
          
  #Generate parameter distribution to sample from
          muP.use<-rbeta(10000, 
                         shape1 = muP.1, 
                         shape2 = muP.2)
          plot(density(muP.use))
          
          
          
    #Point estimate of daily predator mortality    
      p.0<-3*length(dat$tank[dat$chlor==0]) #Total starting number of prawns in chloP-free tanks (3 in each)
      p.1<-sum(dat$p.all_24[dat$chlor==0]) #Total surviving prawns in chlorP-free tanks after day 1
          muP<-(p.0-p.1)/p.0 #Per capita daily death rate by 1-day endpoint = number of deaths/ starting population
          
          
#Step 1.3 Estimate daily prawn mortality in chlorP tanks ################### 
  #Bootstrapping approach for distribution in MC approach
  p.24.ch<-dat$p.all_24[dat$chlor==1] 
    p.24.ch<-3-p.24.ch #convert to deaths
  #Bootstrap
    p.24.ch.samples <- lapply(1:5000, function(i)
      (sample(p.24.ch, replace = T)/3 - muP))
          muPq.mean<-sapply(p.24.ch.samples, mean)
          muPq = mean(muPq.mean)
  #Check distribution      
          descdist(muPq.mean) #beta distribution since mortality rate restricted from 0 to 1
            plot(density(muPq.mean))
            
          fitdist(muPq.mean, distr = 'beta', start = list(shape1=1, shape2=1), method = 'mme')
            
          lines(density(rbeta(10000, 
                              shape1 = fitdist(muPq.mean, distr = 'beta', 
                                               start = list(shape1=1, shape2=1), method = 'mme')$estimate[1],
                              shape2 = fitdist(muPq.mean, distr = 'beta', 
                                               start = list(shape1=1, shape2=1), method = 'mme')$estimate[2])),
                col='red')
  #Extract distribution parameters
            muPq.1 = fitdist(muPq.mean, distr = 'beta', 
                             start = list(shape1=1, shape2=1), method = 'mme')$estimate[1]
            muPq.2 = fitdist(muPq.mean, distr = 'beta', 
                             start = list(shape1=1, shape2=1), method = 'mme')$estimate[2]
            
  #Generate parameter distribution to sample from
          muPq.use<-rbeta(10000, 
                          shape1 = muPq.1,
                          shape2 = muPq.2)
            plot(density(muPq.use))
            
#Plot distribution of predator mortalities to use
  plot(density(muP.use), xlim = c(0,1), ylim = c(0,20), lwd = 2, main = 'Predator mortality distributions',
       xlab = 'daily predator mortality rate') 
    lines(density(muPq.use), col = 'red', lwd=2)
    legend('topright', legend = c('Absent', 'Present'), title = 'ChlorP', lwd=2, col = c('black', 'red'))

#Step 2 Read in parameter values, transmission parameters, and R0 function ############
 parameters=c(#Excluding beta, Phi_Nq, f_Nq, and muPq which will be read into the R0 function
              ##standard snail parameters 
              f_N=0.10, # recruitment rate (from sokolow et al)
              phi_N=10000, # carrying capacity from sokolow et al
              z=0.5, #Proportion of exposed snails that reproduce from sokolow et al
              mu_N=1/60, #Mortality rate from Sokolow et al
              sigma=1/40, #Transition rate from exposed to infected from sokolow et al
              mu_I=1/10 - 1/60, #additional snail death due to infection from sokolow et al
              
              #predator parameters
              alpha=0.003, #attack rate
              Th=0.067,#~crayfish predation limit
              nn=2,
              f_P=0.234/2, #crayfish birth rate from Cervantes-Santiago Aquaculture 2010 paper (/2 for 1:1 female-male ratio)
              phi_P=120,  #crayfish carrying capacity
              mu_P= muP, #observed daily 24 hr crayfish mortality rate in chlorP-free tanks of mesocosm#0.00137, previous rate
              
              #Adult Worm, Miracidia and Circariae Parameters
              mu_W=1/(3.3*365), # death rate of adult worms
              m=0.5, #miracidial shedding rate per reproductive female divided by miracidial mortality; from sokolow et al
              
              #Human parameters
              H=300, #number of humans
              mu_H=1/(60*365) #Assumes 60 year lifespan
            )
#R0 function ###########
    
    get_Ro_beta_lamda<-function(muPq = 0, phi_Nq = 1, beta, lamda, f_Nq = 1) #variable parameters to be sampled
    { 
      f_N<-parameters["f_N"]
      phi_N<-parameters["phi_N"]
      z<-parameters["z"]
      mu_N<-parameters["mu_N"]
      sigma<-parameters["sigma"]
      mu_I<-parameters["mu_I"]
      alpha<-parameters["alpha"]
      nn<-parameters["nn"]
      Th<-parameters["Th"]
      f_P<-parameters["f_P"]
      phi_P<-parameters["phi_P"]
      mu_P<-parameters["mu_P"]
      mu_W<-parameters["mu_W"]
      m<-parameters["m"]
      H<-parameters["H"]
      mu_H<-parameters["mu_H"]
      
      P_eq<-(1-((muPq+mu_P)/f_P))*phi_P #Equilibrium estimate of P given prawn predator parameters
      
      if(P_eq<0){P_eq=0}
      
      #Equilibrium estimate of N given snail parameters
      N_eq = max(uniroot.all(f = function(y){(f_N*f_Nq)*(1-y/(phi_N*phi_Nq)) - 
                                             mu_N - 
                                             (P_eq*alpha*(y/200)^(nn-1))/(1+alpha*Th*(y/200)^nn)}, 
                            c(0, as.numeric(phi_N*phi_Nq))))
      
      if(N_eq < 0){
        N_eq = 0
      }
      
      #print(N_eq)
      
      pred<-(alpha*P_eq*(N_eq/200)^(nn-1))/(1+(alpha*(N_eq/200)^nn*Th))#death rate of snails due to predators given equilibrium estimates of P and N
      
      Ro_est <- sqrt((0.5*beta*H*N_eq*lamda*sigma)/((mu_W+mu_H)*(mu_N+(pred/5)+sigma)*(mu_N+(pred/10)+mu_I)))
      
      return(c(N_eq,P_eq,Ro_est ))
      
    } #End R0 function  
  
  p.dead = parameters["f_P"] - parameters["mu_P"] #muPq value at which predator population is eliminated 
            
#Transmission parameters from max likelihood estimation################# 
  #Obtained using optim function fit to three follow-up epi points; see code 'fit_fin.R'
  fin<-read.csv('Gen1/shortlist_transParams.csv') 
            
  beta.use<- fin$beta[fin$negLL == min(fin$negLL)] #Best fit value of beta
  
  lamda.lo<-fin$lamda1[fin$negLL == min(fin$negLL)] #Best fit value of lo lamda
  lamda.hi<-fin$lamda2[fin$negLL == min(fin$negLL)] #Best fit value of hi lamda
  
  lamda.use<-fin$lamda.twa[fin$negLL == min(fin$negLL)] #Best fit value of twa lamda (use in R0 sims)
  
  #R0 corresponding to best fit parameters
    get_Ro_beta_lamda(muPq = p.dead, beta = beta.use, lamda = lamda.use)
   
  fin$prob <- fin$likelihood / sum(fin$likelihood)  
    plot(density(fin$prob))
       
  beta.025 = min(fin$beta)
  beta.975 = max(fin$beta)
  
  lamda1.025 = min(fin$lamda1)
  lamda1.975 = max(fin$lamda1)
  
  lamda2.025 = min(fin$lamda2)
  lamda2.975 = max(fin$lamda2)
    
  hist(fin$R0)
    r0.025 = min(fin$R0)
    r0.975 = max(fin$R0)
      
#R0 monte carlo simulations with uncertainty from agro parameters and model fitting #######################
  parameters["mu_P"]<-0 #predator mortality corresponds to distributions above in simulations

#MC for control: regular predator mortality, no bottom up effects, transmission parameters      
  control<-as.numeric()
  set.seed(043017)
  
    for(i in 1:1000){
      trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]

      control[i] = get_Ro_beta_lamda(muPq = sample(muP.use, 1),
                                     phi_Nq = 1,
                                     beta = trans$beta,
                                     lamda = trans$lamda.twa)[3] 
      
    }
  
        hist(control)
        abline(v=median(control), col = 'red', lty=2)
        median(control)
        
#MC for atrazine only: regular predator mortality, bottom up effects from atrazine, transmission parameters              
  At.MC<-as.numeric()
  set.seed(043017)
  
        for(i in 1:1000){
          trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]

          At.MC[i] = get_Ro_beta_lamda(muPq = sample(muP.use, 1),
                                         phi_Nq = sample(At.use, 1),
                                         beta = trans$beta,
                                         lamda = trans$lamda.twa)[3] 
          
        }
        
        hist(At.MC)      
        abline(v=median(At.MC), col = 'red', lty=2)
        median(At.MC)
      
#MC for fertilizer only: regular predator mortality, bottom up effects from fertilizer, transmission parameters              
  Fe.MC<-as.numeric()
  set.seed(043017)
  
        for(i in 1:1000){
          trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]

          Fe.MC[i] = get_Ro_beta_lamda(muPq = sample(muP.use, 1),
                                       phi_Nq = sample(Fe.use, 1),
                                       beta = trans$beta,
                                       lamda = trans$lamda.twa)[3] 
          
        }
        
        hist(Fe.MC) 
        abline(v=median(Fe.MC), col = 'red', lty=2)
        median(Fe.MC)
        
#MC for fertilizer/atrazine: regular predator mortality, bottom up effects from atrazine and fertilizer, transmission parameters              
  FeAt.MC<-as.numeric()
  set.seed(043017)
  
        for(i in 1:1000){
          trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]

          FeAt.MC[i] = get_Ro_beta_lamda(muPq = sample(muP.use, 1),
                                       phi_Nq = sample(FeAt.use, 1),
                                       beta = trans$beta,
                                       lamda = trans$lamda.twa)[3] 
          
        }
        
        hist(FeAt.MC) 
        abline(v=median(FeAt.MC), col = 'red', lty=2)
        median(FeAt.MC)
        
#MC for chlorpyrifos only: chlorP predator mortality, no bottom-up effects, transmission parameters              
  Ch.MC<-as.numeric()
  set.seed(043017)
  
        for(i in 1:1000){
          trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]

          Ch.MC[i] = get_Ro_beta_lamda(muPq = sample(muPq.use, 1),
                                         phi_Nq = 1,
                                         beta = trans$beta,
                                         lamda = trans$lamda.twa)[3] 
          
        }
        
        hist(Ch.MC) 
        abline(v=median(Ch.MC), col = 'red', lty=2)
        median(Ch.MC)
        
#MC for atrazine/chlorpyrifos: chlorP predator mortality, atrazine bottom-up effects, transmission parameters              
  At.Ch.MC<-as.numeric()
  set.seed(043017)
  
        for(i in 1:1000){
          trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]

          At.Ch.MC[i] = get_Ro_beta_lamda(muPq = sample(muPq.use, 1),
                                       phi_Nq = sample(At.use, 1),
                                       beta = trans$beta,
                                       lamda = trans$lamda.twa)[3] 
          
        }
        
        hist(At.Ch.MC)
        abline(v=median(At.Ch.MC), col = 'red', lty=2)
        median(At.Ch.MC)
     
#MC for fertilizer/chlorpyrifos: chlorP predator mortality, fertilizer bottom-up effects, transmission parameters              
  Fe.Ch.MC<-as.numeric()
  set.seed(043017)
  
        for(i in 1:1000){
          trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]

          Fe.Ch.MC[i] = get_Ro_beta_lamda(muPq = sample(muPq.use, 1),
                                          phi_Nq = sample(Fe.use, 1),
                                          beta = trans$beta,
                                          lamda = trans$lamda.twa)[3] 
          
        }
        
        hist(Fe.Ch.MC)
        abline(v=median(Fe.Ch.MC), col = 'red', lty=2)
        median(Fe.Ch.MC)
        
#MC for triplets: chlorP predator mortality, fertilizer/atrazine bottom-up effects, transmission parameters              
  Tre.MC<-as.numeric()
  set.seed(043017)
  
        for(i in 1:1000){
          trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]
          
          Tre.MC[i] = get_Ro_beta_lamda(muPq = sample(muPq.use, 1),
                                          phi_Nq = sample(FeAt.use, 1),
                                          beta = trans$beta,
                                          lamda = trans$lamda.twa)[3] 
          
        }
        
        hist(Tre.MC)
        abline(v=median(Tre.MC), col = 'red', lty=2)
        median(Tre.MC)
        
      
#Boxplot for panel A: estimated R0 in treatment groups  ################         
     MC.df<-data.frame('R0' = c(control,
                                At.MC, Fe.MC, FeAt.MC, 
                                Ch.MC, 
                                At.Ch.MC, Fe.Ch.MC, Tre.MC),
                       'treat' = c(rep('Control', times=length(control)),
                                   rep('At', times=length(At.MC)), 
                                   rep('Fe', times=length(Fe.MC)), 
                                   rep('At:Fe', times=length(FeAt.MC)), 
                                   rep('Ch', times=length(Ch.MC)), 
                                   rep('At:Ch', times=length(At.Ch.MC)), 
                                   rep('Ch:Fe', times=length(Fe.Ch.MC)), 
                                   rep('At:Ch:Fe', times=length(Tre.MC))))
     
     MC.df$treat<-factor(MC.df$treat, levels = c('Control','Fe', 'At', 'At:Fe', #Bottom up effects only
                                                 'Ch', #Top-down effects only
                                                 'Ch:Fe', 'At:Ch', 'At:Ch:Fe')) #Both top-down and bottom-up effects

     MC.df<-MC.df[complete.cases(MC.df),]
     
  gg1<-ggplot(MC.df, aes(x=treat, y=R0, group = treat))+
       #Theme formatting
         theme_bw()+
         theme(axis.title=element_text(size=14),
               axis.text=element_text(size=10))+
         scale_y_continuous(breaks = c(0,1,3,5,7), limits=c(0,max(MC.df$R0)+1.5))+
         xlab("")+
         ylab(expression('R'[0]))+
       #Adding data from data frame
         geom_boxplot(width = 0.25, outlier.colour = 'grey60', outlier.shape = 1)+
       #Add treatment labels
         geom_segment(x=1.7, xend=4.3, y=max(MC.df$R0)+0.5, yend=max(MC.df$R0)+0.5, 
                      colour='grey40', lineend='square')+
         geom_text(x=3, y=max(MC.df$R0)+0.85, 
                   label='bottom-up effects', size=5, colour='grey40')+
         geom_segment(x=4.7, xend=5.3, y=max(MC.df$R0)+0.5, yend=max(MC.df$R0)+0.5, 
                      colour='grey40', lineend='square')+
         geom_text(x=5, y=max(MC.df$R0)+1.3, 
                   label='top-down', size=5, colour='grey40')+
         geom_text(x=5, y=max(MC.df$R0)+0.85, 
                   label='effects', size=5, colour='grey40')+
         geom_segment(x=5.7, xend=8.3, y=max(MC.df$R0)+0.5, yend=max(MC.df$R0)+0.5, 
                      colour='grey40', lineend='square')+
         geom_text(x=7, y=max(MC.df$R0)+0.85, 
                   label='bottom-up & top-down effects', size=5, colour='grey40')+
       #Add plot label
         geom_text(label='A', x=0.75, y=max(MC.df$R0)+1.4, size=6)

#Boxplot export in .eps format for nature pub ##############    
  gg1export<-ggplot(MC.df, aes(x=treat, y=R0, group = treat))+
    #Theme formatting
    theme_bw()+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=10))+
    scale_y_continuous(breaks = c(0,1,3,5,7), limits=c(0,max(MC.df$R0)+1.5))+
    xlab("")+
    ylab(expression('R'[0]))+
    #Adding data from data frame
    geom_boxplot(width = 0.25, outlier.colour = 'grey60', outlier.shape = 1)+
    #Add treatment labels
    geom_segment(x=1.7, xend=4.3, y=max(MC.df$R0)+0.5, yend=max(MC.df$R0)+0.5, 
                 colour='grey40', lineend='square')+
    #geom_text(x=3, y=max(MC.df$R0)+0.85, 
    #          label='bottom-up effects', size=5, colour='grey40')+
    geom_segment(x=4.7, xend=5.3, y=max(MC.df$R0)+0.5, yend=max(MC.df$R0)+0.5, 
                 colour='grey40', lineend='square')+
    #geom_text(x=5, y=max(MC.df$R0)+1.3, 
    #          label='top-down', size=5, colour='grey40')+
    #geom_text(x=5, y=max(MC.df$R0)+0.85, 
    #          label='effects', size=5, colour='grey40')+
    geom_segment(x=5.7, xend=8.3, y=max(MC.df$R0)+0.5, yend=max(MC.df$R0)+0.5, 
                 colour='grey40', lineend='square')#+
    #geom_text(x=7, y=max(MC.df$R0)+0.85, 
    #          label='bottom-up & top-down effects', size=5, colour='grey40')+
    #Add plot label
    #geom_text(label='A', x=0.75, y=max(MC.df$R0)+1.4, size=6)
  
  postscript('Panel_A_Boxplot.eps', width = 7.5, height = 4, horizontal = FALSE)
    gg1export
    dev.off()
    
#R0 monte carlo simulations agrochemical uncertainty only #######################
  for(i in 1:length(MC.df)){
    MC.df$Uncertainty = 'Total'
  }

#MC for control with agrochemical uncertainty      
  control.ag<-as.numeric()
    for(i in 1:1000){
      control.ag[i] = get_Ro_beta_lamda(muPq = sample(muP.use, 1),
                                       phi_Nq = 1,
                                       beta = beta.use,
                                       lamda = lamda.use)[3] 
      
    }
    
    hist(control.ag)
    abline(v=median(control.ag), col = 'red', lty=2)
    median(control.ag)
    
#MC for atrazine only: regular predator mortality, bottom up effects from atrazine, transmission parameters              
  At.MC.ag<-as.numeric()
  for(i in 1:1000){
    At.MC.ag[i] = get_Ro_beta_lamda(muPq = sample(muP.use, 1),
                                 phi_Nq = sample(At.use, 1),
                                 beta = beta.use,
                                 lamda = lamda.use)[3] 
    
  }
  
  hist(At.MC.ag)
  abline(v=median(At.MC.ag), col = 'red', lty=2)
  median(At.MC.ag)
  
#MC for fertilizer only: regular predator mortality, bottom up effects from fertilizer, transmission parameters              
  Fe.MC.ag<-as.numeric()
  for(i in 1:1000){
    Fe.MC.ag[i] = get_Ro_beta_lamda(muPq = sample(muP.use, 1),
                                 phi_Nq = sample(Fe.use, 1),
                                 beta = beta.use,
                                 lamda = lamda.use)[3] 
    
  }
  
  hist(Fe.MC.ag) 
  abline(v=median(Fe.MC.ag), col = 'red', lty=2)
  median(Fe.MC.ag)
  
#MC for fertilizer/atrazine: regular predator mortality, bottom up effects from atrazine and fertilizer, transmission parameters              
  FeAt.MC.ag<-as.numeric()
  for(i in 1:1000){
    FeAt.MC.ag[i] = get_Ro_beta_lamda(muPq = sample(muP.use, 1),
                                   phi_Nq = sample(FeAt.use, 1),
                                   beta = beta.use,
                                   lamda = lamda.use)[3] 
    
  }
  
  hist(FeAt.MC.ag) 
  abline(v=median(FeAt.MC.ag), col = 'red', lty=2)
  median(FeAt.MC.ag)
  
#MC for chlorpyrifos only: chlorP predator mortality, no bottom-up effects, transmission parameters              
  Ch.MC.ag<-as.numeric()
  for(i in 1:1000){
    Ch.MC.ag[i] = get_Ro_beta_lamda(muPq = sample(muPq.use, 1),
                                 phi_Nq = 1,
                                 beta = beta.use,
                                 lamda = lamda.use)[3] 
    
  }
  
  hist(Ch.MC.ag)       
  abline(v=median(Ch.MC.ag), col = 'red', lty=2)
  median(Ch.MC.ag)
  
#MC for atrazine/chlorpyrifos: chlorP predator mortality, atrazine bottom-up effects, transmission parameters              
  At.Ch.MC.ag<-as.numeric()
  for(i in 1:1000){
    At.Ch.MC.ag[i] = get_Ro_beta_lamda(muPq = sample(muPq.use, 1),
                                    phi_Nq = sample(At.use, 1),
                                    beta = beta.use,
                                    lamda = lamda.use)[3] 
    
  }
  
  hist(At.Ch.MC.ag)
  abline(v=median(At.Ch.MC.ag), col = 'red', lty=2)
  median(At.Ch.MC.ag)
  
#MC for fertilizer/chlorpyrifos: chlorP predator mortality, fertilizer bottom-up effects, transmission parameters              
  Fe.Ch.MC.ag<-as.numeric()
  for(i in 1:1000){
    Fe.Ch.MC.ag[i] = get_Ro_beta_lamda(muPq = sample(muPq.use, 1),
                                    phi_Nq = sample(Fe.use, 1),
                                    beta = beta.use,
                                    lamda = lamda.use)[3] 
    
  }
  
  hist(Fe.Ch.MC.ag)
  abline(v=median(Fe.Ch.MC.ag), col = 'red', lty=2)
  median(Fe.Ch.MC.ag)
  
#MC for triplets: chlorP predator mortality, fertilizer/atrazine bottom-up effects, transmission parameters              
  Tre.MC.ag<-as.numeric()
  for(i in 1:1000){
    trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]
    
    Tre.MC.ag[i] = get_Ro_beta_lamda(muPq = sample(muPq.use, 1),
                                  phi_Nq = sample(FeAt.use, 1),
                                  beta = beta.use,
                                  lamda = lamda.use)[3] 
    
  }
  
  hist(Tre.MC.ag)
  abline(v=median(Tre.MC.ag), col = 'red', lty=2)
  median(Tre.MC.ag)
  
#Compile in data frame to plot
  MC.df.ag<-data.frame('R0' = c(control.ag,
                             At.MC.ag, Fe.MC.ag, FeAt.MC.ag, 
                             Ch.MC.ag, 
                             At.Ch.MC.ag, Fe.Ch.MC.ag, Tre.MC.ag),
                    'treat' = c(rep('Control', times=length(control.ag)),
                                rep('At', times=length(At.MC.ag)), 
                                rep('Fe', times=length(Fe.MC.ag)), 
                                rep('At:Fe', times=length(FeAt.MC.ag)), 
                                rep('Ch', times=length(Ch.MC.ag)), 
                                rep('At:Ch', times=length(At.Ch.MC.ag)), 
                                rep('Ch:Fe', times=length(Fe.Ch.MC.ag)), 
                                rep('At:Ch:Fe', times=length(Tre.MC.ag))))
  
  MC.df.ag$treat<-factor(MC.df.ag$treat, levels = c('Control','Fe', 'At', 'At:Fe', #Bottom up effects only
                                              'Ch', #Top-down effects only
                                              'Ch:Fe', 'At:Ch', 'At:Ch:Fe')) #Both top-down and bottom-up effects
  
  MC.df.ag<-MC.df.ag[complete.cases(MC.df.ag),]
#R0 monte carlo simulations transmission uncertainty only #######################
  for(i in 1:length(MC.df.ag)){
    MC.df.ag$Uncertainty = 'Agrochemical parameterization'
  }

#Control    
  control.trans<-as.numeric()
  for(i in 1:1000){
    trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]
    
    control.trans[i] = get_Ro_beta_lamda(muPq = muP,
                                   phi_Nq = 1,
                                   beta = trans$beta,
                                   lamda = trans$lamda.twa)[3] 
    
  }
  
  hist(control.trans)
  abline(v=median(control.trans), col = 'red', lty=2)
  median(control.trans)
  
#MC for atrazine only: regular predator mortality, bottom up effects from atrazine, transmission parameters              
  At.MC.trans<-as.numeric()
  for(i in 1:1000){
    trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]
    
    At.MC.trans[i] = get_Ro_beta_lamda(muPq = muP,
                                 phi_Nq = At.mu,
                                 beta = trans$beta,
                                 lamda = trans$lamda.twa)[3] 
    
  }
  
  hist(At.MC.trans)      
  abline(v=median(At.MC.trans), col = 'red', lty=2)
  median(At.MC.trans)
  
#MC for fertilizer only: regular predator mortality, bottom up effects from fertilizer, transmission parameters              
  Fe.MC.trans<-as.numeric()
  for(i in 1:1000){
    trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]
    
    Fe.MC.trans[i] = get_Ro_beta_lamda(muPq = muP,
                                 phi_Nq = Fe.mu,
                                 beta = trans$beta,
                                 lamda = trans$lamda.twa)[3] 
    
  }
  
  hist(Fe.MC.trans) 
  abline(v=median(Fe.MC.trans), col = 'red', lty=2)
  median(Fe.MC.trans)
  
#MC for fertilizer/atrazine: regular predator mortality, bottom up effects from atrazine and fertilizer, transmission parameters              
  FeAt.MC.trans<-as.numeric()
  for(i in 1:1000){
    trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]
    
    FeAt.MC.trans[i] = get_Ro_beta_lamda(muPq = muP,
                                   phi_Nq = FeAt.mu,
                                   beta = trans$beta,
                                   lamda = trans$lamda.twa)[3] 
    
  }
  
  hist(FeAt.MC.trans) 
  abline(v=median(FeAt.MC.trans), col = 'red', lty=2)
  median(FeAt.MC.trans)
  
#MC for chlorpyrifos only: chlorP predator mortality, no bottom-up effects, transmission parameters              
  Ch.MC.trans<-as.numeric()
  for(i in 1:1000){
    trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]
    
    Ch.MC.trans[i] = get_Ro_beta_lamda(muPq = muPq,
                                 phi_Nq = 1,
                                 beta = trans$beta,
                                 lamda = trans$lamda.twa)[3] 
    
  }
  
  hist(Ch.MC.trans) 
  abline(v=median(Ch.MC.trans), col = 'red', lty=2)
  median(Ch.MC.trans)
  
#MC for atrazine/chlorpyrifos: chlorP predator mortality, atrazine bottom-up effects, transmission parameters              
  At.Ch.MC.trans<-as.numeric()
  for(i in 1:1000){
    trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]
    
    At.Ch.MC.trans[i] = get_Ro_beta_lamda(muPq = muPq,
                                    phi_Nq = At.mu,
                                    beta = trans$beta,
                                    lamda = trans$lamda.twa)[3] 
    
  }
  
  hist(At.Ch.MC.trans)
  abline(v=median(At.Ch.MC.trans), col = 'red', lty=2)
  median(At.Ch.MC.trans)
  
#MC for fertilizer/chlorpyrifos: chlorP predator mortality, fertilizer bottom-up effects, transmission parameters              
  Fe.Ch.MC.trans<-as.numeric()
  for(i in 1:1000){
    trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]
    
    Fe.Ch.MC.trans[i] = get_Ro_beta_lamda(muPq = muPq,
                                    phi_Nq = Fe.mu,
                                    beta = trans$beta,
                                    lamda = trans$lamda.twa)[3] 
    
  }
  
  hist(Fe.Ch.MC.trans)
  abline(v=median(Fe.Ch.MC.trans), col = 'red', lty=2)
  median(Fe.Ch.MC.trans)
  
#MC for triplets: chlorP predator mortality, fertilizer/atrazine bottom-up effects, transmission parameters              
  Tre.MC.trans<-as.numeric()
  for(i in 1:1000){
    trans = fin[sample(nrow(fin), size = 1, prob = fin$prob),]
    
    Tre.MC.trans[i] = get_Ro_beta_lamda(muPq = muPq,
                                  phi_Nq = FeAt.mu,
                                  beta = trans$beta,
                                  lamda = trans$lamda.twa)[3] 
    
  }
  
  hist(Tre.MC.trans)
  abline(v=median(Tre.MC.trans), col = 'red', lty=2)
  median(Tre.MC.trans)
  
#Compile in data frame to plot
  MC.df.trans<-data.frame('R0' = c(control.trans,
                                At.MC.trans, Fe.MC.trans, FeAt.MC.trans, 
                                Ch.MC.trans, 
                                At.Ch.MC.trans, Fe.Ch.MC.trans, Tre.MC.trans),
                       'treat' = c(rep('Control', times=length(control.trans)),
                                   rep('At', times=length(At.MC.trans)), 
                                   rep('Fe', times=length(Fe.MC.trans)), 
                                   rep('At:Fe', times=length(FeAt.MC.trans)), 
                                   rep('Ch', times=length(Ch.MC.trans)), 
                                   rep('At:Ch', times=length(At.Ch.MC.trans)), 
                                   rep('Ch:Fe', times=length(Fe.Ch.MC.trans)), 
                                   rep('At:Ch:Fe', times=length(Tre.MC.trans))))
  
  MC.df.trans$treat<-factor(MC.df.trans$treat, levels = c('Control','Fe', 'At', 'At:Fe', #Bottom up effects only
                                                    'Ch', #Top-down effects only
                                                    'Ch:Fe', 'At:Ch', 'At:Ch:Fe')) #Both top-down and bottom-up effects
  
  MC.df.trans<-MC.df.trans[complete.cases(MC.df.trans),]
  
  for(i in 1:length(MC.df.trans)){
    MC.df.trans$Uncertainty = 'Model fitting'
  }
  
#Plot boxplot (alternative panel A) with agrochemical, transmission, and total uncertainty separated ###########
  
  MC.df.all<-rbind(MC.df, MC.df.ag, MC.df.trans)
  MC.df.all$Uncertainty<-factor(MC.df.all$Uncertainty, levels = c('Total',
                                                                  'Model fitting',
                                                                  'Agrochemical parameterization'))
  
gg1.1<-ggplot(MC.df.all, aes(x=treat, y=R0))+
    #Theme formatting
    theme_bw()+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=10))+
    scale_y_continuous(breaks = c(0,1,3,5,7), limits=c(0,max(MC.df$R0)+1.5))+
    xlab("")+
    ylab(expression('R'[0]))+
    #Adding data from data frame
    geom_boxplot(aes(fill = Uncertainty), width = 0.5, outlier.colour = 'grey60', outlier.shape = 1)+
    #Add treatment labels
    geom_segment(x=1.7, xend=4.3, y=max(MC.df$R0)+0.5, yend=max(MC.df$R0)+0.5, 
                 colour='grey40', lineend='square')+
    geom_text(x=3, y=max(MC.df$R0)+0.85, 
              label='bottom-up effects', size=5, colour='grey40')+
    geom_segment(x=4.7, xend=5.3, y=max(MC.df$R0)+0.5, yend=max(MC.df$R0)+0.5, 
                 colour='grey40', lineend='square')+
    geom_text(x=5, y=max(MC.df$R0)+1.3, 
              label='top-down', size=5, colour='grey40')+
    geom_text(x=5, y=max(MC.df$R0)+0.85, 
              label='effects', size=5, colour='grey40')+
    geom_segment(x=5.7, xend=8.3, y=max(MC.df$R0)+0.5, yend=max(MC.df$R0)+0.5, 
                 colour='grey40', lineend='square')+
    geom_text(x=7, y=max(MC.df$R0)+0.85, 
              label='bottom-up & top-down effects', size=5, colour='grey40')+
    #Add plot label
    geom_text(label='A', x=0.75, y=max(MC.df$R0)+1.4, size=6)
  
#Step 3: Bring in chlorP d/r data from chemosphere 2015 paper #####################
  #10-day data with 99% mortality in last group
  ecotox_mod10<-data.frame('dose'=c(rep(0,100), rep(0.64,100), rep(3.2,100), 
                                    rep(6.4,100), rep(32,100), rep(64,100)),
                           'response'=c(rep(0,100), rep(0,100), rep(0,100),
                                        rep(0,40), rep(1,159), 0, 0, rep(1,99)))
     
     ecotox10_mod<-glm(response ~ dose, family=binomial(link="probit"),data=ecotox_mod10)
      summary(ecotox10_mod)

     #Extrapolate response to constant gradient of Chlorpyrifos concentration
     p.ecotox_mod10<-data.frame(dose=seq(from=0, to=150, by=0.1))
      p.ecotox_mod10[, c('mortality', 'st.er')]<-predict(ecotox10_mod, p.ecotox_mod10, 
                                                        type = "response", se.fit=TRUE)[-3]
      
      p.ecotox_mod10$logchlor<-log(p.ecotox_mod10$dose+1)  
    #normalize to 0 mortality at 0 concentration
      p.ecotox_mod10$mortality = p.ecotox_mod10$mortality - p.ecotox_mod10[1,2] 
     
     plot(x=c(0,0.32,0.64,3.2,6.4,32,64), pch = 16,
          y=c(0,0,   0   ,0  ,60/100  ,99/100,99/100)/10, xlab = "ChlorP Concentration", ylab = "%Mortality")
       lines(x=p.ecotox_mod10$dose, y=p.ecotox_mod10$mortality/10, col='red')
       lines(x=p.ecotox_mod10$dose, 
             y=p.ecotox_mod10$mortality/10+1.96*p.ecotox_mod10$st.er/10, 
             lty=2, col='red', cex=0.8)
       lines(x=p.ecotox_mod10$dose, 
             y=p.ecotox_mod10$mortality/10-1.96*p.ecotox_mod10$st.er/10, 
             lty=2, col='red', cex=0.8)
        abline(a = parameters['f_P'] - parameters['mu_P'], b = 0, lty=2, col='blue')
     
     parameters['mu_P']=muP #Treat as additional mortality to observed daily rate
     
  #Calculate R0 with ChlorP only  
     set.seed(043017)
     
     for (i in 1:nrow(p.ecotox_mod10)){
       p.ecotox_mod10$R0[i] <- get_Ro_beta_lamda(muPq = (p.ecotox_mod10$mortality[i]/10), phi_Nq = 1, 
                                                 beta = beta.use, lamda = lamda.use)[3]
     }
     
     set.seed(043017)
     
     for (i in 1:nrow(p.ecotox_mod10)){
       p.ecotox_mod10$R0_lo[i] <- get_Ro_beta_lamda(muPq = (p.ecotox_mod10[i,2]/10 - 1.96*p.ecotox_mod10[i,3]/10), 
                                                    phi_Nq = 1, beta = beta.use, lamda = lamda.use)[3]
     }
     
     set.seed(043017)
     
     for (i in 1:nrow(p.ecotox_mod10)){
       p.ecotox_mod10$R0_up[i] <- get_Ro_beta_lamda(muPq = (p.ecotox_mod10[i,2]/10 + 1.96*p.ecotox_mod10[i,3]/10), 
                                                    phi_Nq = 1, beta = beta.use, lamda = lamda.use)[3]
     }
     
  #Calculate equilibrium P estimates  
     set.seed(043017)
     
     for (i in 1:nrow(p.ecotox_mod10)){
       p.ecotox_mod10$P_eq[i] <- get_Ro_beta_lamda(muPq = (p.ecotox_mod10[i,2]/10),phi_Nq = 1, 
                                                   beta = beta.use, lamda = lamda.use)[2]
     }
      p.ecotox_mod10$P_eq = p.ecotox_mod10$P_eq/200 #convert numbers to densities
     
      set.seed(043017)
      
     for (i in 1:nrow(p.ecotox_mod10)){
       p.ecotox_mod10$P_eq_up[i] <- get_Ro_beta_lamda(muPq = (p.ecotox_mod10[i,2]/10 + 1.96*p.ecotox_mod10[i,3]/10), 
                                                      phi_Nq = 1, beta = beta.use, lamda = lamda.use)[2]
     }
      p.ecotox_mod10$P_eq_up = p.ecotox_mod10$P_eq_up/200 #convert numbers to densities
      
      set.seed(043017)
      
     for (i in 1:nrow(p.ecotox_mod10)){
       p.ecotox_mod10$P_eq_lo[i] <- get_Ro_beta_lamda(muPq = (p.ecotox_mod10[i,2]/10 - 1.96*p.ecotox_mod10[i,3]/10), 
                                                      phi_Nq = 1, beta = beta.use, lamda = lamda.use)[2]
     }
      p.ecotox_mod10$P_eq_lo = p.ecotox_mod10$P_eq_lo/200 #convert numbers to densities
      
     
  #Calculate equilibrium N estimates 
      set.seed(043017)
      
     for (i in 1:nrow(p.ecotox_mod10)){
       p.ecotox_mod10$N_eq[i] <- get_Ro_beta_lamda(muPq = (p.ecotox_mod10[i,2]/10), phi_Nq = 1, 
                                                   beta = beta.use, lamda = lamda.use)[1]
     }
      p.ecotox_mod10$N_eq = p.ecotox_mod10$N_eq/200 #convert numbers to densities
      
      set.seed(043017)
      
     for (i in 1:nrow(p.ecotox_mod10)){
       p.ecotox_mod10$N_eq_up[i] <- get_Ro_beta_lamda(muPq = (p.ecotox_mod10[i,2]/10 + 1.96*p.ecotox_mod10[i,3]/10), 
                                                      phi_Nq = 1, beta = beta.use, lamda = lamda.use)[1]
     }
      p.ecotox_mod10$N_eq_up = p.ecotox_mod10$N_eq_up/200 #convert numbers to densities
      
      set.seed(043017)
      
     for (i in 1:nrow(p.ecotox_mod10)){
       p.ecotox_mod10$N_eq_lo[i] <- get_Ro_beta_lamda(muPq = (p.ecotox_mod10[i,2]/10 - 1.96*p.ecotox_mod10[i,3]/10), 
                                                      phi_Nq = 1, beta = beta.use, lamda = lamda.use)[1]
     }
      p.ecotox_mod10$N_eq_lo = p.ecotox_mod10$N_eq_lo/200 #convert numbers to densities
     
#Step 3.1 Plot R0, P_eq, and N_eq response to modeled ChlorP P-mortality across concentrations ##############
#R0 response to chlorP concentration
gg2.0<-ggplot(p.ecotox_mod10, aes(x=dose, y=R0))+
       theme_bw()+
       theme(axis.text=element_text(size=10), axis.title=element_text(size=14))+
       scale_x_continuous(breaks=c(0,15,30,35), labels = c('0','15','30','64'), 
                           limits=c(0,35))+
       geom_vline(aes(xintercept = 35), colour = 'grey90') + 
       scale_y_continuous(breaks=c(0,1,
                                   round(max(p.ecotox_mod10$R0), digits = 2)), 
                          limits=c(0,max(p.ecotox_mod10$R0+0.2)))+
       ylab(expression('R'[0]))+
       xlab(expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep = '')))+
       geom_line(size=1)+
       geom_line(aes(y=R0_lo), linetype=2, alpha = 0.5)+
       geom_line(aes(y=R0_up), linetype=2, alpha = 0.5)+
       geom_segment(y = 3.6, yend = round(max(p.ecotox_mod10$R0), digits = 2), x = 31.5, xend = 35, 
                    col = 'white', size=2)+
       geom_segment(y = 3.6, yend = round(max(p.ecotox_mod10$R0), digits = 2), x = 31.5, xend = 35, 
                    col = 'black', linetype = 3, size=1)+
       geom_segment(y = 3.6, yend = round(max(p.ecotox_mod10$R0), digits = 2), x = 34, xend = 35, 
                     col = 'black', size=1)+
       geom_text(label='B', x=0, y=3.6, size=6)

  #Panel B export in .eps format for Nature pub    
      gg2export<-ggplot(p.ecotox_mod10, aes(x=dose, y=R0))+
        theme_bw()+
        theme(axis.text=element_text(size=10), axis.title=element_text(size=14))+
        scale_x_continuous(breaks=c(0,15,30,35), labels = c('0','15','30','64'), 
                           limits=c(0,35))+
        geom_vline(aes(xintercept = 35), colour = 'grey90') + 
        scale_y_continuous(breaks=c(0,1,
                                    round(max(p.ecotox_mod10$R0), digits = 2)), 
                           limits=c(0,max(p.ecotox_mod10$R0+0.2)))+
        ylab(expression('R'[0]))+
        xlab(expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep = '')))+
        geom_line(aes(y=R0_lo), linetype=2, col = 'grey70')+
        geom_line(aes(y=R0_up), linetype=2, col = 'grey70')+
        geom_line(size=1)+
        geom_segment(y = 3.6, yend = round(max(p.ecotox_mod10$R0), digits = 2), x = 31.5, xend = 35, 
                     col = 'white', size=2)+
        geom_segment(y = 3.6, yend = round(max(p.ecotox_mod10$R0), digits = 2), x = 31.5, xend = 35, 
                     col = 'black', linetype = 3, size=1)+
        geom_segment(y = 3.6, yend = round(max(p.ecotox_mod10$R0), digits = 2), x = 34, xend = 35, 
                     col = 'black', size=1)#+
        #geom_text(label='B', x=0, y=3.75, size=6)
      
      postscript('Panel_B_ChlorP_R0.eps', width = 6, height = 3, horizontal = FALSE)
      gg2export
      dev.off()     
      
#Predator population response to chlorP concentration
gg2.1<-ggplot(p.ecotox_mod10, aes(x=dose, y=P_eq))+
       theme_bw()+
       theme(axis.text=element_text(size=8.5), axis.title=element_text(size=14))+
       scale_x_continuous(breaks=c(0, 10,
                                   min(p.ecotox_mod10$dose[p.ecotox_mod10$R0 > 1]), #min concentration required for disease elimination
                                   round(min(p.ecotox_mod10$dose[p.ecotox_mod10$P_eq == 0])), #min concentration for predator extirpation
                                   64), 
                          limits=c(0,round(min(p.ecotox_mod10$dose[p.ecotox_mod10$P_eq_lo == 0]))+0.25))+
       scale_y_continuous(breaks=c(0,0.2,0.4), limits=c(0,0.425))+
       xlab('')+
       ylab(expression(paste('Predator density (P/m'^'2',')', sep = '')))+
       geom_line(size=1, colour = 'red')+
       geom_line(aes(y=P_eq_lo), linetype=2, colour = 'red', alpha = 0.5)+
       geom_line(aes(y=P_eq_up), linetype=2, colour = 'red', alpha = 0.5)+
       geom_text(label='A', x=-0.5, y=0.43, size=6)

#Predator pop by chlorP in .eps format for Nature pub
  gg2.1export<-ggplot(p.ecotox_mod10, aes(x=dose, y=P_eq))+
    theme_bw()+
    theme(axis.text=element_text(size=8.5), axis.title=element_text(size=14))+
    scale_x_continuous(breaks=c(0, 10,
                                min(p.ecotox_mod10$dose[p.ecotox_mod10$R0 > 1]), #min concentration required for disease elimination
                                round(min(p.ecotox_mod10$dose[p.ecotox_mod10$P_eq == 0])), #min concentration for predator extirpation
                                64), 
                       limits=c(0,round(min(p.ecotox_mod10$dose[p.ecotox_mod10$P_eq_lo == 0]))+0.25))+
    scale_y_continuous(breaks=c(0,0.2,0.4), limits=c(0,0.425))+
    xlab('')+
    ylab(expression(paste('Predator density (P/m'^'2',')', sep = '')))+
    geom_line(size=1, colour = 'red')+
    geom_line(aes(y=P_eq_lo), linetype=2, colour = 'red')+
    geom_line(aes(y=P_eq_up), linetype=2, colour = 'red')#+
    #geom_text(label='A', x=-0.5, y=0.43, size=6)

  postscript('Supp_FigA_Pred_Dens_ChlorP.eps', width = 6, height = 3, horizontal = FALSE)
    gg2.1export
    dev.off()     
  
#Snail population response to chlorP concentration (through predator population)
gg2.2<-ggplot(p.ecotox_mod10, aes(x=dose, y=N_eq))+
         theme_bw()+
         theme(axis.text=element_text(size=8.5), axis.title=element_text(size=14))+
         scale_x_continuous(breaks=c(0, 10,
                                     min(p.ecotox_mod10$dose[p.ecotox_mod10$R0 > 1]), #min concentration required for disease elimination
                                     round(min(p.ecotox_mod10$dose[p.ecotox_mod10$P_eq == 0])), #min concentration for predator extirpation
                                   64), 
                            limits=c(0,round(min(p.ecotox_mod10$dose[p.ecotox_mod10$P_eq_lo == 0]))+0.25))+
         scale_y_continuous(breaks=c(0,20,40), limits=c(0,45))+
         xlab(expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep = '')))+
         ylab(expression(paste('Snail density (N/m'^'2',')', sep = '')))+
         geom_line(size=1, colour = 'blue')+
         geom_line(aes(y=N_eq_lo), linetype=2, colour = 'blue', alpha = 0.5)+
         geom_line(aes(y=N_eq_up), linetype=2, colour = 'blue', alpha = 0.5)+
         geom_text(label='B', x=-0.5, y=44.5, size=6)

#Snail pop by chlorP in .eps format for Nature pub
gg2.2export<-ggplot(p.ecotox_mod10, aes(x=dose, y=N_eq))+
  theme_bw()+
  theme(axis.text=element_text(size=8.5), axis.title=element_text(size=14))+
  scale_x_continuous(breaks=c(0, 10,
                              min(p.ecotox_mod10$dose[p.ecotox_mod10$R0 > 1]), #min concentration required for disease elimination
                              round(min(p.ecotox_mod10$dose[p.ecotox_mod10$P_eq == 0])), #min concentration for predator extirpation
                              64), 
                     limits=c(0,round(min(p.ecotox_mod10$dose[p.ecotox_mod10$P_eq_lo == 0]))+0.25))+
  scale_y_continuous(breaks=c(0,20,40), limits=c(0,45))+
  xlab(expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep = '')))+
  ylab(expression(paste('Snail density (N/m'^'2',')', sep = '')))+
  geom_line(size=1, colour = 'blue')+
  geom_line(aes(y=N_eq_lo), linetype=2, colour = 'blue')+
  geom_line(aes(y=N_eq_up), linetype=2, colour = 'blue')#+
  #geom_text(label='B', x=-0.5, y=44.5, size=6)

  postscript('Supp_FigB_Snail_Dens_ChlorP.eps', width = 6, height = 3, horizontal = FALSE)
  gg2.2export
  dev.off()     

     
#Step 3.2 Enhanced snail pop dynamics from Atrazine from Baxter data ##############
  atra.df<-data.frame('atra' = c(0,1,10,100),
                      'logatra' = log(c(0,1,10,100)+1),
                      'phiNq' = c(1,1.2888,1.6535,2.3215))

    plot(atra.df$logatra, atra.df$phiNq, pch = 16)
    
    atra_mod<-glm(phiNq ~ logatra, data=atra.df)
      summary(atra_mod)
      confint(atra_mod)
    
    abline(a = atra_mod$coefficients[1], b = atra_mod$coefficients[2], lty=2, col='red')
    
    #Extrapolate response to constant gradient of Atrazine concentration
    atra.predict<-data.frame('atra' = seq(from=0, to=100, by=0.1),
                             'logatra' = log(seq(from=0, to=100, by=0.1)+1))
    atra.predict[, c('phiNq', 'st.er')]<-predict(atra_mod, atra.predict, 
                                                       type = "response", se.fit=TRUE)[-3]
    
    atra.predict$phiNq<-atra.predict$phiNq + (1 - atra.predict[1,3]) #Normalize to 1
    
    set.seed(043017)
    
    for (i in 1:nrow(atra.predict)){
      atra.predict$R0[i] <- get_Ro_beta_lamda(muPq = p.dead, 
                                              phi_Nq = atra.predict[i,3], 
                                              beta = beta.use, lamda = lamda.use)[3]
    }
    
    set.seed(043017)
    
    for (i in 1:nrow(atra.predict)){
      atra.predict$R0_lo[i] <- get_Ro_beta_lamda(muPq = p.dead, 
                                                 phi_Nq = atra.predict[i,3] - 1.96*atra.predict[i,4], 
                                                 beta = beta.use, lamda = lamda.use)[3]
    }
    
    set.seed(043017)
    
    for (i in 1:nrow(atra.predict)){
      atra.predict$R0_up[i] <- get_Ro_beta_lamda(muPq = p.dead, 
                                                 phi_Nq = atra.predict[i,3] + 1.96*atra.predict[i,4], 
                                                 beta = beta.use, lamda = lamda.use)[3]
    }
    
gg2.5<-ggplot(atra.predict, aes(x=atra, y=R0))+
        theme_bw()+
        theme(axis.text=element_text(size=10), axis.title=element_text(size=14))+
        scale_x_continuous(limits=c(0,100))+
        scale_y_continuous(breaks=c(3,4,5,6), 
                           limits=c(3,6))+
        ylab(expression('R'[0]))+
        xlab(expression(paste('Atrazine concentration (', mu, 'g/L)', sep = '')))+
        geom_line(size=1)+
        geom_line(aes(y=R0_lo), linetype=2, alpha = 0.5)+
        geom_line(aes(y=R0_up), linetype=2, alpha = 0.5)+
        #geom_vline(xintercept = round(log(100+1), digits = 1), lty=4, col = 'red')+
        geom_text(label='C', x=0, y=6, size=6)
#Panel C with log atrazine scale ###############
gg2.5.1<-ggplot(atra.predict, aes(x=logatra, y=R0))+
          theme_bw()+
          theme(axis.text=element_text(size=10), axis.title=element_text(size=14))+
          scale_x_continuous(breaks=log(c(0,1,10,100)+1),
                             limits = log(c(0,100)+1),
                             labels = c('0', '1', '10', '100'))+
          scale_y_continuous(breaks=c(3,4,5,6), 
                             limits=c(3,6))+
          ylab(expression('R'[0]))+
          xlab(expression(paste('log+1 Atrazine concentration (', mu, 'g/L)', sep = '')))+
          geom_line(size=1)+
          geom_line(aes(y=R0_lo), linetype=2, alpha = 0.5)+
          geom_line(aes(y=R0_up), linetype=2, alpha = 0.5)+
          #geom_vline(xintercept = round(log(100+1), digits = 1), lty=4, col = 'red')+
          geom_text(label='C', x=0, y=5.9, size=6)
    
#Panel C export in .eps for nature pub 
gg2.5export<-ggplot(atra.predict, aes(x=logatra, y=R0))+
  theme_bw()+
  theme(axis.text=element_text(size=10), axis.title=element_text(size=14))+
  scale_x_continuous(breaks=log(c(0,1,10,100)+1),
                     limits = log(c(0,100)+1),
                     labels = c('0', '1', '10', '100'))+
  scale_y_continuous(breaks=c(3,4,5,6), 
                     limits=c(3,6))+
  ylab(expression('R'[0]))+
  xlab(expression(paste('log+1 Atrazine concentration (', mu, 'g/L)', sep = '')))+
  geom_line(aes(y=R0_lo), linetype=2, colour = 'grey60')+
  geom_line(aes(y=R0_up), linetype=2, colour = 'grey60')+
  geom_line(size=1)#+
  #geom_vline(xintercept = round(log(100+1), digits = 1), lty=4, col = 'red')+
  #geom_text(label='C', x=0, y=6, size=6)

  postscript('Panel_C_Atr_R0.eps', width = 6, height = 3, horizontal = FALSE)
  gg2.5export
  dev.off()     

#Step 3.3 combine chlorpyrifos and atrazine d/r data in heat map ##############
  r0.atra.chlor<-data.frame("Atra" = rep(seq(from=0, to=100, by=0.1), times = 641),
                            "logatra" = rep(log(seq(from=0, to=100, by=0.1)+1), times = 641),
                            "dose" = rep(seq(0,64,0.1), each = 1001),#Chlorpyrifos concentration
                            "logdose" = rep(log(seq(0,64,0.1)+1), each = 1001),
                            "phi_Nq" = 0,
                            "rate" = 0,
                            "R0" =0,
                            'P_eq'=0,
                            'N_eq' = 0)

     
     r0.atra.chlor$rate<-(predict(ecotox10_mod, r0.atra.chlor, 
                                  type = "response", se.fit=TRUE)$fit)/10 #fill mortality rate data from model
     
     r0.atra.chlor$rate = r0.atra.chlor$rate - r0.atra.chlor[1,6] #normalize to 0; no excess death with no concentration
     
     r0.atra.chlor$phi_Nq<-(predict(atra_mod, r0.atra.chlor, 
                                  type = "response", se.fit=TRUE)$fit) #fill bottom-up effect data from model
     
     r0.atra.chlor$phi_Nq<-r0.atra.chlor$phi_Nq + (1 - r0.atra.chlor[1,5]) #Normalize to 1
     
     set.seed(043017)
     
      r0.atra.chlor[ ,c(9,8,7)] = t(mapply(FUN = get_Ro_beta_lamda, muPq = r0.atra.chlor[,6], phi_Nq = r0.atra.chlor[,5],
                                       MoreArgs = list(beta = beta.use, lamda = lamda.use), SIMPLIFY = T))
     
     plot(x = r0.atra.chlor$logatra[r0.atra.chlor$dose == 64],
          y = r0.atra.chlor$R0[r0.atra.chlor$dose == 64],
          xlab = 'atrazine concentration', ylab = 'R0',
          type = 'l', lwd = 2)
     
     
#Step 3.3.1 plot the heat map #################
  gg3<-ggplot(r0.atra.chlor, aes(x=Atra, y=dose, fill=R0))+
       theme_bw()+
       #scale_fill_brewer(type = 'div', palette = 'RdYlGn', direction = -1)+
       scale_fill_distiller(palette = "Spectral")+
       scale_x_continuous(breaks = c(0,25,50,75,100), limits = c(0,100))+
       scale_y_continuous(breaks = c(0,20,40,60), limits = c(0,64))+
       #scale_y_continuous(breaks = c(0,10,20), limits = c(0,25))+
       geom_tile(size=0.01)+
       coord_equal()+
       geom_vline(xintercept = c(0,25,50,75,100), col = 'lightgrey', alpha = 0.5)+
       geom_hline(yintercept = c(0,20,40,60), col = 'lightgrey', alpha = 0.25)+
       labs(y=expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep = '')), 
            #(x=expression(paste('Predator mortality rate (', mu[P][,][q], ')', sep = '')), axis label for predator mortality rate
            x=expression(paste('Atrazine concentration (', mu, 'g/L)', sep = '')))+
       theme(axis.ticks=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=14),
             legend.title=element_text(size=15), legend.text=element_text(size=12))+
       geom_text(label='D', x=-3, y=65, size=6, alpha=.50)    
     
#Plot in .eps format for nature pub
     gg3export<-ggplot(r0.atra.chlor, aes(x=Atra, y=dose, fill=R0))+
       theme_bw()+
       #scale_fill_brewer(type = 'div', palette = 'RdYlGn', direction = -1)+
       scale_fill_distiller(palette = "Spectral")+
       scale_x_continuous(breaks = c(0,25,50,75,100), limits = c(0,100))+
       scale_y_continuous(breaks = c(0,20,40,60), limits = c(0,64))+
       #scale_y_continuous(breaks = c(0,10,20), limits = c(0,25))+
       geom_tile(size=0.01)+
       coord_equal()+
       geom_vline(xintercept = c(0,25,50,75,100), col = 'grey10')+
       geom_hline(yintercept = c(0,20,40,60), col = 'grey10')+
       labs(y=expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep = '')), 
            #(x=expression(paste('Predator mortality rate (', mu[P][,][q], ')', sep = '')), axis label for predator mortality rate
            x=expression(paste('Atrazine concentration (', mu, 'g/L)', sep = '')))+
       theme(axis.ticks=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=14),
             legend.title=element_text(size=15), legend.text=element_text(size=12))#+
       #geom_text(label='D', x=-3, y=65, size=6, alpha=.50)    
     
     postscript('Panel_D_HeatMap.eps', width = 6, height = 5, horizontal = FALSE)
     gg3export
     dev.off()     
     
#Plot heat map with log atrazine scale ##############
     r0.atra.chlor<-data.frame("Atra" = 0,
                               "logatra" = rep(seq(from=0, to=log(101), by=0.01), times = 641),
                               "dose" = rep(seq(0,64,0.1), each = 462),
                               "logdose" = rep(log(seq(0,64,0.1)+1), each = 462),
                               "phi_Nq" = 0,
                               "rate" = 0,
                               "R0" =0)
     
     r0.atra.chlor$Atra<-exp(r0.atra.chlor$logatra)
     
     r0.atra.chlor$rate<-(predict(ecotox10_mod, r0.atra.chlor, 
                                  type = "response", se.fit=TRUE)$fit)/10 #fill mortality rate data from model
     
     r0.atra.chlor$rate = r0.atra.chlor$rate - r0.atra.chlor[1,6] #normalize to 0; no excess death with no concentration
     
     r0.atra.chlor$phi_Nq<-(predict(atra_mod, r0.atra.chlor, 
                                    type = "response", se.fit=TRUE)$fit) #fill bottom-up effect data from model
     
     r0.atra.chlor$phi_Nq<-r0.atra.chlor$phi_Nq + (1 - r0.atra.chlor[1,5]) #Normalize to 1
     
     #r0.atra.chlor[ ,c(9,8,7)] = t(mapply(FUN = get_Ro_beta_lamda, muPq = r0.atra.chlor[,6], phi_Nq = r0.atra.chlor[,5],
     #                                     MoreArgs = list(beta = beta.use, lamda = lamda.use), SIMPLIFY = T))
     
     for(i in 1:nrow(r0.atra.chlor)){
       r0.atra.chlor[i,7] = get_Ro_beta_lamda(muPq = r0.atra.chlor[i,6],
                                              beta = beta.use,
                                              lamda = lamda.use,
                                              phi_Nq = r0.atra.chlor[i,5])[3]
       if(i %% 1000==0) print(i)
     }
     
     gg3.1<-ggplot(r0.atra.chlor, aes(x=logatra, y=dose, fill=R0))+
       theme_bw()+
       #scale_fill_brewer(type = 'div', palette = 'RdYlGn', direction = -1)+
       scale_fill_distiller(palette = "Spectral")+
       scale_x_continuous(breaks = log(c(0,1,10,100)+1), 
                          limits = log(c(0,100)+1),
                          labels = c('0','1','10','100'))+
       scale_y_continuous(breaks = c(0,20,40,60), limits = c(0,64))+
       geom_raster(interpolate = TRUE)+
       coord_equal(ratio = 1/20)+
       #geom_vline(xintercept = c(unique(r0.atra.chlor$logatra)), col = 'black', alpha = 0.5)+
       geom_hline(yintercept = c(0,20,40,60), col = 'lightgrey', alpha = 0.25)+
       labs(y=expression(paste('Chlorpyrifos concentration (', mu, 'g/L)', sep = '')), 
            x=expression(paste('log+1 Atrazine concentration (', mu, 'g/L)', sep = '')))+
       theme(axis.ticks=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=14),
             legend.title=element_text(size=15), legend.text=element_text(size=12))+
       geom_text(label='D', x=-3, y=65, size=6, alpha=.50)     
#Plot heat map with both axes log transformed ################
  gg3.2<-ggplot(r0.atra.chlor, aes(x=logatra, y=logdose, fill=R0))+
    theme_bw()+
    #scale_fill_brewer(type = 'div', palette = 'RdYlGn', direction = -1)+
    scale_fill_distiller(palette = "Spectral")+
    scale_x_continuous(breaks = log(c(0,1,10,100)+1), 
                       limits = log(c(0,100)+1),
                       labels = c('0','1','10','100'))+
    scale_y_continuous(breaks = log(c(0,1,10,60)+1), 
                       limits = log(c(0,65)+1),
                       labels = c('0','1','10','60'))+
    geom_raster(interpolate = TRUE)+
    coord_equal()+
    #geom_vline(xintercept = c(unique(r0.atra.chlor$logatra)), col = 'black', alpha = 0.5)+
    geom_hline(yintercept = c(0,20,40,60), col = 'lightgrey', alpha = 0.25)+
    labs(y=expression(paste('log+1 Chlorpyrifos concentration (', mu, 'g/L)', sep = '')), 
         x=expression(paste('log+1 Atrazine concentration (', mu, 'g/L)', sep = '')))+
    theme(axis.ticks=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=14),
          legend.title=element_text(size=15), legend.text=element_text(size=12))+
    geom_text(label='D', x=-3, y=65, size=6, alpha=.50)   
     
#Step 4 put panels together for final figure ###########################
  grid.arrange(gg1, gg2.0, gg2.1, gg2.2, gg2.5, gg3, ncol=7, nrow=6, 
               layout_matrix=rbind(c(1,1,1,1,1,1,1), 
                                   c(1,1,1,1,1,1,1),
                                   c(2,2,2,3,3,3,3),
                                   c(2,2,2,4,4,4,4),
                                   c(5,5,5,6,6,6,6),
                                   c(5,5,5,6,6,6,6)))
     
#Move pred and snail pop trajectories to supplementary figure ###########################
  grid.arrange(gg1, gg2.0, gg2.5, gg3, ncol=2, nrow=5, 
               layout_matrix=rbind(c(1,1),
                                   c(1,1),
                                   c(2,3),
                                   c(4,4),
                                   c(4,4)))
     
  grid.arrange(gg2.1, gg2.2, ncol=1, nrow=2, 
               layout_matrix=rbind(c(1),
                                   c(2)))
#Plots with log atrazine scale in panels C&D ###################
  pdf('Fig3.pdf', height = 11, width = 8)
  grid.arrange(gg1, gg2.0, gg2.5.1, gg3.1, ncol=2, nrow=5, 
               layout_matrix=rbind(c(1,1),
                                   c(1,1),
                                   c(2,3),
                                   c(4,4),
                                   c(4,4)))
  dev.off()