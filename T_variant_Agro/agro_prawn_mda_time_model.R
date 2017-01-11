#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(deSolve)
source('T_variant_Agro/agroC_data.R')

#Model structure to dynamically vary agrochemical application, prawn introduction, and MDA ####################
agroc_prawn_mda<-function(chem, chem.days, chem.k, chem.med, 
                          prawn.add, prawn.harvest, prawn.n, prawn.l, 
                          mda.add, mda.eff, mda.cov,
                          n.years){
  
  parameters2['k_Q1'] = chem.k
  parameters2['cov'] = mda.cov
  time = seq(1, 365*n.years, by=1)
  
  nstart = c(S = eqbm.mda$S,
             E = eqbm.mda$E,
             I = eqbm.mda$I,
             Wt = eqbm.mda$Wt,
             Wu = eqbm.mda$Wu,
             P = eqbm.mda$P,
             L = eqbm.mda$L,
             Q1 = eqbm.mda$Q1)
  
  events = c(chem.days, mda.add, prawn.add, prawn.add, prawn.add+prawn.harvest)
  events.all<-matrix(ncol = length(events), nrow = n.years)
  
  for(i in 1:length(events)){
    for(j in 1:n.years){
      events.all[j,i] = events[i] + 365*(j-1)
    }
  }
  
  events.fill = c(t(events.all))
  harvests<-events.all[,length(events)]
  
  events.df<-data.frame(var=rep(c(rep('Q1', length(chem.days)), 'Wt', 'P', 'L', 'P'), times = n.years),
                        time = events.fill,
                        value = rep(c(rep(chem.med, length(chem.days)), (1-mda.eff), prawn.n, prawn.l, 0), times = n.years),
                        method = rep(c(rep("add", length(chem.days)), "mult", "rep", "rep", "rep"), times = n.years))
  
  output.all = as.data.frame(ode(nstart, time, agro_time_wt_wu, parameters2,
                                 events = list(data = events.df)))
  output.all$N = output.all$S + output.all$E + output.all$I
  output.all$W = output.all$Wt*mda.cov + output.all$Wu*(1-mda.cov)
  output.all$P.B = output.all$P * (as.numeric(parameters['a.p'])*(output.all$L/10)^as.numeric(parameters['b.p']))/10
  harvest.size = sum(output.all$P.B[output.all$time%in%harvests])
  
  par(mfrow = c(2,2), mar = c(2.5,4,1.5,2.0))
  plot(x = output.all$time, y = output.all$N, lwd=2, xlab = 'time', ylab = 'snail dynamics', 
       type = 'l', ylim = c(0, max(output.all$N)))
  lines(output.all$time, output.all$S, lwd=2, col = 'blue')
  lines(output.all$time, output.all$E, lwd=2, col = 'orange')
  lines(output.all$time, output.all$I, lwd=2, col = 'red')
  legend('topleft', legend = c('N', 'S', 'E', 'I'), lty = 1, col = c('black', 'blue', 'orange', 'red'), cex = 0.6)
  
  
  plot(output.all$time, output.all$P, lwd=2, type = 'l', col = 'green', 
       xlab = 'time', ylab = 'prawn dynamics', ylim = c(0, 2100),
       main = paste('total harvest = ', round(harvest.size/1000, digits=1), 'kg', sep = ''))
  lines(output.all$time[output.all$time >= prawn.add], output.all$L[output.all$time >= prawn.add], col = 'lightblue', lty=2)
  lines(output.all$time, output.all$P.B, lwd=2, col = 'darkgreen')
  title(main = paste('total harvest = ', round(harvest.size/1000, digits=1), 'kg', sep = ''), cex=0.1)
  legend('topright', legend = c('P', 'P.B', 'L'), lty = 1, col = c('green', 'darkgreen', 'lightblue'), cex = 0.6)
  
  plot(output.all$time, output.all$W, lwd=2, col = 'purple', type = 'l', 
       xlab = 'time', ylab = 'mean worm burden', ylim = c(0,120),
       main = paste('average worm burden/person = ', round(mean(output.all$W[output.all$time>mda.add]), digits=1), sep=''))  
    lines(output.all$time, output.all$Wt, lwd=2, col = 'purple', lty=2)
    lines(output.all$time, output.all$Wu, lwd=2, col = 'purple', lty=3)
    legend('topright', legend = c('W', 'Wt', 'Wu'), lty = c(1,2,3), col = 'purple', cex = 0.6)
    legend('top', legend = c(paste('coverage = ', mda.cov, sep=''),
                             paste('efficacy = ', mda.eff), sep=''), cex = 0.6)
  
  plot(output.all$time, output.all$Q1,lwd=2, col = 'gold2', xlab = 'time', ylab = 'chem concentration (ppb)', type = 'l',
       main = chem)
  
  return(c(mean(output.all$W[output.all$time>mda.add]),harvest.size))
  
}

#Simulations with prawn, agrochemical, and MDA introduction ###################
#Run with no MDA, prawns, or agrochemical for baseline
agroc_prawn_mda(chem = 'Chlorpyrifos',
                chem.days = 0,
                chem.k = chlor.k,
                chem.med = 0,
                prawn.add = 1,
                prawn.harvest = harvest.time,
                prawn.n = 0,
                prawn.l = 25,
                mda.add = 1,
                mda.eff = 0,
                mda.cov = 0.8,
                n.years = 3)
#Run with annual MDA
agroc_prawn_mda(chem = 'Chlorpyrifos',
                chem.days = 0,
                chem.k = chlor.k,
                chem.med = 0,
                prawn.add = 1,
                prawn.harvest = harvest.time,
                prawn.n = 0,
                prawn.l = 25,
                mda.add = 1,
                mda.eff = 0.95,
                mda.cov = 0.8,
                n.years = 3)

#Run with annual MDA, annual prawn introduction
agroc_prawn_mda(chem = 'Chlorpyrifos',
                chem.days = 0,
                chem.k = chlor.k,
                chem.med = 0,
                prawn.add = 1,
                prawn.harvest = harvest.time,
                prawn.n = 2*area,
                prawn.l = 25,
                mda.add = 1,
                mda.eff = 0.95,
                mda.cov = 0.8,
                n.years = 3)

#Run with annual MDA, annual prawn introduction, annual ChlorP introducation
agroc_prawn_mda(chem = 'Chlorpyrifos',
                chem.days = chlor.days,
                chem.k = chlor.k,
                chem.med = med.chlor,
                prawn.add = 1,
                prawn.harvest = harvest.time,
                prawn.n = 2*area,
                prawn.l = 25,
                mda.add = 1,
                mda.eff = 0.95,
                mda.cov = 0.8,
                n.years = 3)

#Run with annual MDA, annual prawn introduction DELAYED 200 days, annual ChlorP introducation
agroc_prawn_mda(chem = 'Chlorpyrifos',
                chem.days = chlor.days,
                chem.k = chlor.k,
                chem.med = med.chlor,
                prawn.add = 200,
                prawn.harvest = harvest.time,
                prawn.n = 2*area,
                prawn.l = 25,
                mda.add = 1,
                mda.eff = 0.95,
                mda.cov = 0.8,
                n.years = 3)

#Run with annual MDA, annual prawn introduction, annual Terbufos introducation MUST CHANGE TOXICITY MODEL
agroc_prawn_mda(chem = 'Terbufos',
                chem.days = terb.days,
                chem.k = terb.k,
                chem.med = med.terb,
                prawn.add = 1,
                prawn.harvest = harvest.time,
                prawn.n = 2*area,
                prawn.l = 25,
                mda.add = 1,
                mda.eff = 0.95,
                mda.cov = 0.8,
                n.years = 3)

#Take out plotting functions and return values of interest for plots investigating different outcomes ##################
agroc_prawn_mda_noplot<-function(chem, chem.days, chem.k, chem.med, 
                          prawn.add, prawn.harvest, prawn.n, prawn.l, 
                          mda.add, mda.eff, mda.cov,
                          n.years){
  
  parameters2['k_Q1'] = chem.k
  parameters2['cov'] = mda.cov
  time = seq(1, 365*n.years, by=1)
  
  nstart = c(S = eqbm.mda$S,
             E = eqbm.mda$E,
             I = eqbm.mda$I,
             Wt = eqbm.mda$Wt,
             Wu = eqbm.mda$Wu,
             P = eqbm.mda$P,
             L = eqbm.mda$L,
             Q1 = eqbm.mda$Q1)
  
  events = c(chem.days, mda.add, prawn.add, prawn.add, prawn.add+prawn.harvest)
  events.all<-matrix(ncol = length(events), nrow = n.years)
  
  for(i in 1:length(events)){
    for(j in 1:n.years){
      events.all[j,i] = events[i] + 365*(j-1)
    }
  }
  
  events.fill = c(t(events.all))
  harvests<-events.all[,length(events)]
  
  events.df<-data.frame(var=rep(c(rep('Q1', length(chem.days)), 'Wt', 'P', 'L', 'P'), times = n.years),
                        time = events.fill,
                        value = rep(c(rep(chem.med, length(chem.days)), (1-mda.eff), prawn.n, prawn.l, 0), times = n.years),
                        method = rep(c(rep("add", length(chem.days)), "mult", "rep", "rep", "rep"), times = n.years))
  
  output.all = as.data.frame(ode(nstart, time, agro_time_wt_wu, parameters2,
                                 events = list(data = events.df)))
  output.all$N = output.all$S + output.all$E + output.all$I
  output.all$W = output.all$Wt*mda.cov + output.all$Wu*(1-mda.cov)
  output.all$P.B = output.all$P * (as.numeric(parameters['a.p'])*(output.all$L/10)^as.numeric(parameters['b.p']))/10
  harvest.size = sum(output.all$P.B[output.all$time%in%harvests])
  
  
  return(c(mean(output.all$W[output.all$time>mda.add]),harvest.size))
  
}

#When should prawns be added given annual MDA and agro introduction? ####################
prawn.add.chlor<-data.frame(add = c(1:365),
                      worm = 0,
                      harvest = 0)

for(i in 1:nrow(prawn.add.chlor)){
  prawn.add.chlor[i,c(2,3)] = agroc_prawn_mda_noplot(chem = 'Chlorpyrifos',
                                             chem.days = chlor.days,
                                             chem.k = chlor.k,
                                             chem.med = med.chlor,
                                             prawn.add = prawn.add.chlor[i,1],
                                             prawn.harvest = harvest.time,
                                             prawn.n = 2*area,
                                             prawn.l = 25,
                                             mda.add = 1,
                                             mda.eff = 0.95,
                                             mda.cov = 0.8,
                                             n.years = 3)[c(1,2)]
  
}

plot(x = prawn.add.chlor$add-1, y = prawn.add.chlor$harvest/1000, type = 'l', lwd = 2,
     xlab = 'days since insecticide', xlim = c(0, 365 - harvest.time),
     ylab = 'outputs', col = 'red', ylim = c(0,max(c(prawn.add.chlor$harvest/1000, prawn.add.chlor$worm))+1))
  lines(x = prawn.add.chlor$add-1, y = prawn.add.chlor$worm, col = "red", lty=2, lwd=2)
  legend('bottom', lty = c(1,2), legend = c('harvest size (kg)', 'average worm burden'), cex = 0.75)
  
  min(prawn.add.chlor$add[prawn.add.chlor$harvest >= max(round(prawn.add.chlor$harvest))-1])
  p.add.best.chlor<-prawn.add.chlor$add[prawn.add.chlor$worm == min(prawn.add.chlor$worm)]

prawn.add.terb<-data.frame(add = c(1:365),
                            worm = 0,
                            harvest = 0)

for(i in 1:nrow(prawn.add.terb)){ #MUST CHANGE TOXICITY MODEL IN DIFF EQUATIONS
  prawn.add.terb[i,c(2,3)] = agroc_prawn_mda_noplot(chem = 'Terbufos',
                                                     chem.days = terb.days,
                                                     chem.k = terb.k,
                                                     chem.med = med.terb,
                                                     prawn.add = prawn.add.terb[i,1],
                                                     prawn.harvest = harvest.time,
                                                     prawn.n = 2*area,
                                                     prawn.l = 25,
                                                     mda.add = 1,
                                                     mda.eff = 0.95,
                                                     mda.cov = 0.8,
                                                     n.years = 3)[c(1,2)]
  
}

lines(x = prawn.add.terb$add-1, y = prawn.add.terb$worm, col = "purple", lty=2, lwd=2)
lines(x = prawn.add.terb$add-1, y = prawn.add.terb$harvest/1000+0.1, col = "purple", lty=1, lwd=2)


prawn.add.esfen<-data.frame(add = c(1:365),
                           worm = 0,
                           harvest = 0)

for(i in 1:nrow(prawn.add.esfen)){ #MUST CHANGE TOXICITY MODEL IN DIFF EQUATIONS
  prawn.add.esfen[i,c(2,3)] = agroc_prawn_mda_noplot(chem = 'Terbufos',
                                                    chem.days = esfen.days,
                                                    chem.k = esfen.k,
                                                    chem.med = med.esfen,
                                                    prawn.add = prawn.add.esfen[i,1],
                                                    prawn.harvest = harvest.time,
                                                    prawn.n = 2*area,
                                                    prawn.l = 25,
                                                    mda.add = 1,
                                                    mda.eff = 0.95,
                                                    mda.cov = 0.8,
                                                    n.years = 3)[c(1,2)]
  
}

lines(x = prawn.add.esfen$add-1, y = prawn.add.esfen$worm, col = "black", lty=2, lwd=2)
lines(x = prawn.add.esfen$add-1, y = prawn.add.esfen$harvest/1000+0.2, col = "black", lty=1, lwd=2)
legend('bottomright', lwd=2, col = c('red', 'purple', 'black'), legend = c('Chlorpyrifos', 'Terbufos', 'Esfenvalerate'), cex = 0.6)

#When should MDA administration occur given fixed agrochemical and prawn introduction? ################
when.mda<-data.frame(mda = c(1:365),
                     out = 0)

for(i in 1:nrow(when.mda)){
  when.mda[i,2] = agroc_prawn_mda_noplot(chem = 'Chlorpyrifos',
                                         chem.days = chlor.days,
                                         chem.k = chlor.k,
                                         chem.med = med.chlor,
                                         prawn.add = p.add.best.chlor,
                                         prawn.harvest = harvest.time,
                                         prawn.n = 2*area,
                                         prawn.l = 25,
                                         mda.add = when.mda[i,1],
                                         mda.eff = 0.95,
                                         mda.cov = 0.8,
                                         n.years = 3)[1]
  
}

par(mfrow=c(1,1))
plot(x = when.mda$mda, y = when.mda$out, xlab = 'day of MDA administration', ylab = 'harvest/worm burden ratio',
     type='l', lwd=2)

op.mda.add<-when.mda$mda[when.mda$out == min(when.mda$out)]


agroc_prawn_mda(chem = 'Chlorpyrifos',
                       chem.days = chlor.days,
                       chem.k = chlor.k,
                       chem.med = med.chlor,
                       prawn.add = p.add.best.chlor,
                       prawn.harvest = harvest.time,
                       prawn.n = 2*area,
                       prawn.l = 25,
                       mda.add = op.mda.add,
                       mda.eff = 0.95,
                       mda.cov = 0.8,
                       n.years = 3)


#optimization function to find optimal timing of MDA and prawn introduction, holding chlorpyrifos introduction constant######
#alter function to have only two timing inputs
agroc_prawn_mda_optim<-function(prawn.add=1, mda.add=1){
  
  chem = 'Chlorpyrifos'
  chem.days = chlor.days
  chem.k = chlor.k
  chem.med = med.chlor
  
  prawn.harvest = harvest.time
  prawn.n = 2*area
  prawn.l = 25
  
  mda.eff = 0.95
  mda.cov = 0.8
  n.years = 3
  
  parameters2['k_Q1'] = chem.k
  parameters2['cov'] = mda.cov
  time = seq(1, 365*n.years, by=1)
  
  nstart = c(S = eqbm.mda$S,
             E = eqbm.mda$E,
             I = eqbm.mda$I,
             Wt = eqbm.mda$Wt,
             Wu = eqbm.mda$Wu,
             P = eqbm.mda$P,
             L = eqbm.mda$L,
             Q1 = eqbm.mda$Q1)
  
  events = c(chem.days, mda.add, prawn.add, prawn.add, prawn.add+prawn.harvest)
  events.all<-matrix(ncol = length(events), nrow = n.years)
  
  for(i in 1:length(events)){
    for(j in 1:n.years){
      events.all[j,i] = events[i] + 365*(j-1)
    }
  }
  
  events.fill = c(t(events.all))
  harvests<-events.all[,length(events)]
  
  events.df<-data.frame(var=rep(c(rep('Q1', length(chem.days)), 'Wt', 'P', 'L', 'P'), times = n.years),
                        time = events.fill,
                        value = rep(c(rep(chem.med, length(chem.days)), (1-mda.eff), prawn.n, prawn.l, 0), times = n.years),
                        method = rep(c(rep("add", length(chem.days)), "mult", "rep", "rep", "rep"), times = n.years))
  
  output.all = as.data.frame(ode(nstart, time, agro_time_wt_wu, parameters2,
                                 events = list(data = events.df)))
  output.all$N = output.all$S + output.all$E + output.all$I
  output.all$W = output.all$Wt*mda.cov + output.all$Wu*(1-mda.cov)
  output.all$P.B = output.all$P * (as.numeric(parameters['a.p'])*(output.all$L/10)^as.numeric(parameters['b.p']))/10
  harvest.size = sum(output.all$P.B[output.all$time%in%harvests])
  
  return(mean(output.all$W[output.all$time>mda.add])/harvest.size)
  
}

agroc_prawn_mda_optim(prawn.add = 200, mda.add = 5)

when.mda<-data.frame(mda = c(1:365),
                     out = 0)

for(i in 1:nrow(when.mda)){
  when.mda[i,2] = agroc_prawn_mda_optim(prawn.add = 200, mda.add = when.mda[i,1])

}

par(mfrow=c(1,1))
plot(x = when.mda$mda, y = when.mda$out, xlab = 'day of MDA administration', ylab = 'worm burden',
     type='l', lwd=2, ylim = c(15,20), col = 'purple')

op.mda.add<-when.mda$mda[when.mda$out == min(when.mda$out)]

agroc_prawn_mda(chem = 'Chlorpyrifos',
                chem.days = chlor.days,
                chem.k = chlor.k,
                chem.med = med.chlor,
                prawn.add = 200,
                prawn.harvest = harvest.time,
                prawn.n = 2*area,
                prawn.l = 25,
                mda.add = op.mda.add,
                mda.eff = 0.95,
                mda.cov = 0.8,
                n.years = 3)