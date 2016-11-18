#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Model to compare influence of different insecticides on R0

#Get data, provided by Neal; load packages####################
data10<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/Cray.10d.LC50.2012.act2.csv')
data4<-read.csv('C:/Users/chris_hoover/Documents/RemaisWork/Schisto/Data/Halstead_etal/Cray.LC50.2012.csv')
require(drc)

data10ag<-aggregate(data10$Dead, by = list(Chem = data10$Chem,
                                           Conc = data10$Conc), FUN = sum)

data4ag<-aggregate(data4$Dead, by = list(Chem = data4$Chem,
                                         Conc = data4$Conc), FUN = sum)

data4ag$total = 5
  colnames(data4ag)[3] = 'dead'
data10ag$total = 5
  colnames(data10ag)[3] = 'dead'

mod4d <- drm(dead/total~Conc, Chem, weights=total, data=data4ag, fct=LL.2(), 
              type="binomial")

fx<-function(df1, ins){
  dfname = paste(ins, 'df', sep = '')
  dfname = subset(df1, Chem == ins)
  dfname
}


#Malathion ###########################
mal<-subset(data10, Chem == 'Mal')
  mal.c<-unique(mal$Conc)
  mal.d<-as.numeric()
    for(i in 1:length(mal.c)){
    mal.d[i] = sum(mal$Dead[mal$Conc == mal.c[i]])/5
    }
  
mal_mod<-drm(Dead ~ Conc, data = mal, fct = LL2.2())

  summary(mal_mod)  
  plot(mal_mod)

#Extrapolate response to constant gradient of Malathion concentration
mal.df<-data.frame(Conc=seq(from=0, to=max(mal.c), by=1))
mal.df[, c('mortality', 'st.er')]<-predict(mal_mod, mal.df, 
                                             type = "response", se.fit=TRUE)

plot(x=mal.c, y = mal.d/10, pch = 16,
     xlab = "Concentration (ppb)", ylab = "%Mortality",
     main = 'Malathion', ylim = c(0,0.1))
  lines(x=mal.df$Conc, y=mal.df$mortality/10, col='red')
  lines(x=mal.df$Conc, 
        y=mal.df$mortality/10+1.96*mal.df$st.er/10, 
        lty=2, col='red', cex=0.8)
  lines(x=mal.df$Conc, 
        y=mal.df$mortality/10-1.96*mal.df$st.er/10, 
        lty=2, col='red', cex=0.8)  
#Chlorpyrifos #####################  
chlor<-subset(data10, Chem == 'Chlor')
  chlor.c<-unique(chlor$Conc)
  chlor.d<-as.numeric()
  for(i in 1:length(chlor.c)){
    chlor.d[i] = sum(chlor$Dead[chlor$Conc == chlor.c[i]])/5
  }
  
chlor_mod<-drm(Dead ~ Conc, data = chlor, fct = LL2.2())

  summary(chlor_mod)  
  plot(chlor_mod)
  
  b = coef(chlor_mod)[2]
  a = coef(chlor_mod)[1]
  x=c(0:64)
  
  y = ((a/b)*(x/b)^(a-1))*(1+(x/b)^a)^-2
  
#Extrapolate response to constant gradient of Chlorpyrifos concentration
chlor.df<-data.frame(Conc=seq(from=0, to=max(chlor.c), by=1))
  chlor.df[, c('mortality', 'st.er')]<-predict(chlor_mod, chlor.df, 
                                             type = "response", se.fit=TRUE)
  
plot(x=chlor.c, y=chlor.d/10, pch = 16,
     xlab = "Concentration (ppb)", ylab = "%Mortality",
     main = 'Chlorpyrifos')
  lines(x=chlor.df$Conc, y=chlor.df$mortality/10, col='red')
  lines(x=chlor.df$Conc, 
        y=chlor.df$mortality/10+1.96*chlor.df$st.er/10, 
        lty=2, col='red', cex=0.8)
  lines(x=chlor.df$Conc, 
        y=chlor.df$mortality/10-1.96*chlor.df$st.er/10, 
        lty=2, col='red', cex=0.8)
  
#Terbufos #####################  
terb<-subset(data10, Chem == 'Terb')
  terb.c<-unique(terb$Conc)
  terb.d<-as.numeric()
  for(i in 1:length(terb.c)){
    terb.d[i] = sum(terb$Dead[terb$Conc == terb.c[i]])/5
  }
  
terb_mod<-drm(Dead ~ Conc, data = terb, fct = LL2.2())
  
  summary(terb_mod)  
  plot(terb_mod)
  
  #Extrapolate response to constant gradient of Terbufos concentration
terb.df<-data.frame(Conc=seq(from=0, to=max(terb.c), by=1))
  terb.df[, c('mortality', 'st.er')]<-predict(terb_mod, terb.df, 
                                               type = "response", se.fit=TRUE)
  
plot(x=terb.c, y=terb.d/10, pch = 16,
     xlab = "Concentration (ppb)", ylab = "%Mortality",
     main = 'Terbufos')
  lines(x=terb.df$Conc, y=terb.df$mortality/10, col='red')
  lines(x=terb.df$Conc, 
        y=terb.df$mortality/10+1.96*terb.df$st.er/10, 
        lty=2, col='red', cex=0.8)
  lines(x=terb.df$Conc, 
        y=terb.df$mortality/10-1.96*terb.df$st.er/10, 
        lty=2, col='red', cex=0.8)  
#Lambda-cyhalothrin #####################  
lamcy<-subset(data10, Chem == 'Lambda')
  lamcy.c<-unique(lamcy$Conc)
  lamcy.d<-as.numeric()
  for(i in 1:length(lamcy.c)){
    lamcy.d[i] = sum(lamcy$Dead[lamcy$Conc == lamcy.c[i]])/5
  }
  
lamcy_mod<-drm(Dead ~ Conc, data = lamcy, fct = LL2.2())
  
  summary(lamcy_mod)  
  plot(lamcy_mod)
  
  #Extrapolate response to constant gradient of Lambda-cyhalothrin concentration
lamcy.df<-data.frame(Conc=seq(from=0, to=max(lamcy.c), by=0.01))
  lamcy.df[, c('mortality', 'st.er')]<-predict(lamcy_mod, lamcy.df, 
                                              type = "response", se.fit=TRUE)
  
plot(x=lamcy.c, y=lamcy.d/10, pch = 16,
     xlab = "Concentration (ppb)", ylab = "%Mortality",
     main = 'lambda-cyhalothrin')
  lines(x=lamcy.df$Conc, y=lamcy.df$mortality/10, col='red')
  lines(x=lamcy.df$Conc, 
        y=lamcy.df$mortality/10+1.96*lamcy.df$st.er/10, 
        lty=2, col='red', cex=0.8)
  lines(x=lamcy.df$Conc, 
        y=lamcy.df$mortality/10-1.96*lamcy.df$st.er/10, 
        lty=2, col='red', cex=0.8)  
#esfenvalerate #####################  
esfen<-subset(data10, Chem == 'Esfen')
  esfen.c<-unique(esfen$Conc)
  esfen.d<-as.numeric()
  for(i in 1:length(esfen.c)){
    esfen.d[i] = sum(esfen$Dead[esfen$Conc == esfen.c[i]])/5
  }
  
esfen_mod<-drm(Dead ~ Conc, data = esfen, fct = LL2.2())
  
  summary(esfen_mod)  
  plot(esfen_mod)
  
  #Extrapolate response to constant gradient of Esfenvalerate concentration
esfen.df<-data.frame(Conc=seq(from=0, to=max(esfen.c), by=0.1))
  esfen.df[, c('mortality', 'st.er')]<-predict(esfen_mod, esfen.df, 
                                              type = "response", se.fit=TRUE)
  
plot(x=esfen.c, y=esfen.d/10, pch = 16,
     xlab = "Concentration (ppb)", ylab = "%Mortality",
     main = 'Esfenvalerate')
  lines(x=esfen.df$Conc, y=esfen.df$mortality/10, col='red')
  lines(x=esfen.df$Conc, 
        y=esfen.df$mortality/10+1.96*esfen.df$st.er/10, 
        lty=2, col='red', cex=0.8)
  lines(x=esfen.df$Conc, 
        y=esfen.df$mortality/10-1.96*esfen.df$st.er/10, 
        lty=2, col='red', cex=0.8)  
  
#Permethrin #####################  
perm<-subset(data10, Chem == 'Perm')
  perm.c<-unique(perm$Conc)
  perm.d<-as.numeric()
  for(i in 1:length(perm.c)){
    perm.d[i] = sum(perm$Dead[perm$Conc == perm.c[i]])/5
  }
  
perm_mod<-drm(Dead ~ Conc, data = perm, fct = LL2.2())
  
  summary(perm_mod)  
  plot(perm_mod)
  
#Extrapolate response to constant gradient of Permethrin concentration
perm.df<-data.frame(Conc=seq(from=0, to=max(perm.c), by=0.01))
  perm.df[, c('mortality', 'st.er')]<-predict(perm_mod, perm.df, 
                                              type = "response", se.fit=TRUE)
  
plot(x=perm.c, y=perm.d/10, pch = 16,
     xlab = "Concentration (ppb)", ylab = "%Mortality",
     main = 'Permethrin')
  lines(x=perm.df$Conc, y=perm.df$mortality/10, col='red')
  lines(x=perm.df$Conc, 
        y=perm.df$mortality/10+1.96*perm.df$st.er/10, 
        lty=2, col='red', cex=0.8)
  lines(x=perm.df$Conc, 
        y=perm.df$mortality/10-1.96*perm.df$st.er/10, 
        lty=2, col='red', cex=0.8)  
  
  
  
#Get parameters, R0 function ####################
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
    mu_P= 0.0380952, #observed daily 24 hr crayfish mortality rate in chlorP-free tanks of mesocosm#0.00137, previous rate
    
    #Adult Worm, Miracidia and Circariae Parameters
    mu_W=1/(3.3*365), # death rate of adult worms
    m=0.5, #miracidial shedding rate per reproductive female divided by miracidial mortality; from sokolow et al
    
    #Human parameters
    H=300, #number of humans
    mu_H=1/(60*365) #Assumes 60 year lifespan
  )
p.dead = parameters["f_P"] - parameters["mu_P"] #muPq value at which predator population is eliminated 
  
fin<-read.csv('shortlist_transParams.csv') 
  beta.use<- fin$beta[fin$negLL == min(fin$negLL)] #Best fit value of beta
  lamda.use<-fin$lamda.twa[fin$negLL == min(fin$negLL)] #Best fit value of twa lamda (use in R0 sims)

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
    
    print(N_eq)
    
    pred<-(alpha*P_eq*(N_eq/200)^(nn-1))/(1+(alpha*(N_eq/200)^nn*Th))#death rate of snails due to predators given equilibrium estimates of P and N
    
    Ro_est <- sqrt((0.5*beta*H*N_eq*lamda*sigma)/((mu_W+mu_H)*(mu_N+(pred/5)+sigma)*(mu_N+(pred/10)+mu_I)))
    
    return(c(N_eq,P_eq,Ro_est ))
    
  } #End R0 function 

#Get atrazine d-r data ##############
atra.df<-data.frame('atra' = c(0,1,10,100),
                    'logatra' = log(c(0,1,10,100)+1),
                    'phiNq' = c(1,1.2888,1.6535,2.3215))

atra_mod<-glm(phiNq ~ logatra, data=atra.df)

#Extrapolate response to constant gradient of Atrazine concentration
atra.predict<-data.frame('atra' = seq(from=0, to=100, by=0.1),
                         'logatra' = log(seq(from=0, to=100, by=0.1)+1))
atra.predict[, c('phiNq', 'st.er')]<-predict(atra_mod, atra.predict, 
                                             type = "response", se.fit=TRUE)[-3]

#plot heat maps with organophosphate insecticides coupled with atrazine ###################
#Malathion #############
  mal.df<-data.frame(Conc=rep(seq(from=0, to=40, by=0.1), times=93))
  mal.df[, c('mortality', 'st.er')]<-predict(mal_mod, mal.df, 
                                             type = "response", se.fit=TRUE)
  mal.df$logatra = rep(seq(from=0, to=log(101), by=0.05), each = length(seq(from=0, to=40, by=0.1)))
  mal.df$Atra = exp(mal.df$logatra)
  mal.df$phi_Nq<-(predict(atra_mod, mal.df, 
                                 type = "response", se.fit=TRUE)$fit) #fill bottom-up effect data from model
  
  for(i in 1:nrow(mal.df)){
    mal.df[i,7] = get_Ro_beta_lamda(muPq = mal.df[i,2]/10,
                                           beta = beta.use,
                                           lamda = lamda.use,
                                           phi_Nq = mal.df[i,6])[3]
  }
  colnames(mal.df)[7]<-'R0'
  mal.df$chem = 'Malathion'
  
#Chlorpyrifos #################
  chlor.df<-data.frame(Conc=rep(seq(from=0, to=40, by=0.1), times=93))
  chlor.df[, c('mortality', 'st.er')]<-predict(chlor_mod, chlor.df, 
                                             type = "response", se.fit=TRUE)
  chlor.df$logatra = rep(seq(from=0, to=log(101), by=0.05), each = length(seq(from=0, to=40, by=0.1)))
  chlor.df$Atra = exp(chlor.df$logatra)
  chlor.df$phi_Nq<-(predict(atra_mod, chlor.df, 
                          type = "response", se.fit=TRUE)$fit) #fill bottom-up effect data from model
  
  for(i in 1:nrow(chlor.df)){
    chlor.df[i,7] = get_Ro_beta_lamda(muPq = chlor.df[i,2]/10,
                                    beta = beta.use,
                                    lamda = lamda.use,
                                    phi_Nq = chlor.df[i,6])[3]
  }
  colnames(chlor.df)[7]<-'R0'
  chlor.df$chem = 'Chlorpyrifos'
  
#Terbufos ###################
  terb.df<-data.frame(Conc=rep(seq(from=0, to=40, by=0.1), times=93))
  terb.df[, c('mortality', 'st.er')]<-predict(terb_mod, terb.df, 
                                               type = "response", se.fit=TRUE)
  terb.df$logatra = rep(seq(from=0, to=log(101), by=0.05), each = length(seq(from=0, to=40, by=0.1)))
  terb.df$Atra = exp(terb.df$logatra)
  terb.df$phi_Nq<-(predict(atra_mod, terb.df, 
                            type = "response", se.fit=TRUE)$fit) #fill bottom-up effect data from model
  
  for(i in 1:nrow(terb.df)){
    terb.df[i,7] = get_Ro_beta_lamda(muPq = terb.df[i,2]/10,
                                      beta = beta.use,
                                      lamda = lamda.use,
                                      phi_Nq = terb.df[i,6])[3]
  }
  colnames(terb.df)[7]<-'R0'
  terb.df$chem = 'Terbufos'  

#put together and plot #################  
org.all<-rbind(mal.df, chlor.df, terb.df) 
  org.all$chem<-factor(org.all$chem, levels = c('Chlorpyrifos', 'Malathion', 'Terbufos'))
  
peak.org<-data.frame(z = c(18.4,28,36.6), 
                        chem = c('Malathion', 'Chlorpyrifos', 'Terbufos'))
  peak.org$chem<-factor(peak.org$chem, levels = c('Chlorpyrifos', 'Malathion', 'Terbufos'))
  
median.org<-data.frame(z = c(0.778,5.81,1.435), 
                          chem = c('Malathion', 'Chlorpyrifos', 'Terbufos'))
  median.org$chem<-factor(median.org$chem, levels = c('Chlorpyrifos', 'Malathion', 'Terbufos'))
  
heat.org<-ggplot(org.all, aes(x=logatra, y=Conc, fill=R0))+
            theme_bw()+
            scale_fill_distiller(palette = "Spectral")+
            geom_raster(interpolate = TRUE)+
            coord_equal(1/10)+
            facet_wrap(~chem, ncol=3)+
            geom_hline(aes(yintercept = z), peak.org, linetype = 2)+
            geom_hline(aes(yintercept = z), median.org, linetype = 2, col = 'lightgrey')+
            scale_x_continuous(breaks = log(c(0,1,10,100)+1), 
                               limits = log(c(0,100)+1),
                               labels = c('0','1','10','100'))+
            scale_y_continuous(breaks = c(0,10,20,30,40), limits = c(0,40))+
            labs(y=expression(paste('insecticide concentration (', mu, 'g/L)', sep = '')), 
                 x=expression(paste('log+1 Atrazine concentration (', mu, 'g/L)', sep = '')))+
            theme(axis.ticks=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=14),
                  legend.title=element_text(size=15), legend.text=element_text(size=12))
  
#plot heat maps with pyrethroid insecticides coupled with atrazine ###################  
#Lambda-cyhalothrin ###############
  lamcy.df<-data.frame(Conc=rep(seq(from=0, to=10, by=0.05), times=93))
  lamcy.df[, c('mortality', 'st.er')]<-predict(lamcy_mod, lamcy.df, 
                                              type = "response", se.fit=TRUE)
  lamcy.df$logatra = rep(seq(from=0, to=log(101), by=0.05), each = length(seq(from=0, to=10, by=0.05)))
  lamcy.df$Atra = exp(lamcy.df$logatra)
  lamcy.df$phi_Nq<-(predict(atra_mod, lamcy.df, 
                           type = "response", se.fit=TRUE)$fit) #fill bottom-up effect data from model
  
  for(i in 1:nrow(lamcy.df)){
    lamcy.df[i,7] = get_Ro_beta_lamda(muPq = lamcy.df[i,2]/10,
                                     beta = beta.use,
                                     lamda = lamda.use,
                                     phi_Nq = lamcy.df[i,6])[3]
  }
  colnames(lamcy.df)[7]<-'R0'
  lamcy.df$chem = 'L_Cyhalothrin'  
  
#Esfenvalerate #####################
  esfen.df<-data.frame(Conc=rep(seq(from=0, to=10, by=0.05), times=93))
  esfen.df[, c('mortality', 'st.er')]<-predict(esfen_mod, esfen.df, 
                                               type = "response", se.fit=TRUE)
  esfen.df$logatra = rep(seq(from=0, to=log(101), by=0.05), each = length(seq(from=0, to=10, by=0.05)))
  esfen.df$Atra = exp(esfen.df$logatra)
  esfen.df$phi_Nq<-(predict(atra_mod, esfen.df, 
                            type = "response", se.fit=TRUE)$fit) #fill bottom-up effect data from model
  
  for(i in 1:nrow(esfen.df)){
    esfen.df[i,7] = get_Ro_beta_lamda(muPq = esfen.df[i,2]/10,
                                      beta = beta.use,
                                      lamda = lamda.use,
                                      phi_Nq = esfen.df[i,6])[3]
  }
  colnames(esfen.df)[7]<-'R0'
  esfen.df$chem = 'Esfenvalerate' 
  
#Permethrin ###################
  perm.df<-data.frame(Conc=rep(seq(from=0, to=10, by=0.05), times=93))
  perm.df[, c('mortality', 'st.er')]<-predict(perm_mod, perm.df, 
                                               type = "response", se.fit=TRUE)
  perm.df$logatra = rep(seq(from=0, to=log(101), by=0.05), each = length(seq(from=0, to=10, by=0.05)))
  perm.df$Atra = exp(perm.df$logatra)
  perm.df$phi_Nq<-(predict(atra_mod, perm.df, 
                            type = "response", se.fit=TRUE)$fit) #fill bottom-up effect data from model
  
  for(i in 1:nrow(perm.df)){
    perm.df[i,7] = get_Ro_beta_lamda(muPq = perm.df[i,2]/10,
                                      beta = beta.use,
                                      lamda = lamda.use,
                                      phi_Nq = perm.df[i,6])[3]
  }
  colnames(perm.df)[7]<-'R0'
  perm.df$chem = 'Permethrin' 
  
#put together and plot #################  
py.all<-rbind(lamcy.df, esfen.df, perm.df)
  py.all$chem<-factor(py.all$chem, levels = c('L_Cyhalothrin', 'Esfenvalerate', 'Permethrin'))


peak.py<-data.frame(z = c(1.77,1.03,5.98), 
                      chem = c('L_Cyhalothrin', 'Esfenvalerate', 'Permethrin'))
  peak.py$chem<-factor(peak.py$chem, levels = c('L_Cyhalothrin', 'Esfenvalerate', 'Permethrin'))

median.py<-data.frame(z = c(0.649,0.311,1.42), 
                        chem = c('L_Cyhalothrin', 'Esfenvalerate', 'Permethrin'))
  median.py$chem<-factor(median.py$chem, levels = c('L_Cyhalothrin', 'Esfenvalerate', 'Permethrin'))

heat.py<-ggplot(py.all, aes(x=logatra, y=Conc, fill=R0))+
          theme_bw()+
          scale_fill_distiller(palette = "Spectral")+
          geom_raster(interpolate = TRUE)+
          coord_equal(1/3)+
          facet_wrap(~chem, ncol=3)+
          geom_hline(aes(yintercept = z), peak.py, linetype = 2)+
          geom_hline(aes(yintercept = z), median.py, linetype = 2, col = 'lightgrey')+
          scale_x_continuous(breaks = log(c(0,1,10,100)+1), 
                             limits = log(c(0,100)+1),
                             labels = c('0','1','10','100'))+
          scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0,10))+
          labs(y=expression(paste('insecticide concentration (', mu, 'g/L)', sep = '')), 
               x=expression(paste('log+1 Atrazine concentration (', mu, 'g/L)', sep = '')))+
          theme(axis.ticks=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=14),
                legend.title=element_text(size=15), legend.text=element_text(size=12))
