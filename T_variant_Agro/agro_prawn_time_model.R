#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Schisto model with time-variant agrochemical concentration
  #To do list
    #
    #
    #How to include transport from agrochemical application site to water contact site? ###################

require(deSolve)
source('T_variant_Agro/agroC_data.R')

#Model structure to dynamically vary agrochemical application and prawn introduction ####################
agroc_prawn<-function(chem, chem.days, chem.k, chem.med, 
                      prawn.add, prawn.harvest, prawn.n, prawn.l, 
                      n.years){
 parameters['k_Q1'] = chem.k
 time = seq(1, 365*n.years, by=1)
 
 nstart = c(S = eqbm$S,
            E = eqbm$E,
            I = eqbm$I,
            W = eqbm$W,
            P = eqbm$P,
            L = eqbm$L,
            Q1 = eqbm$Q1)
 
 events = c(chem.days, prawn.add, prawn.add, prawn.add+prawn.harvest)
 events.all<-matrix(ncol = length(events), nrow = n.years)
 
 for(i in 1:length(events)){
   for(j in 1:n.years){
     events.all[j,i] = events[i] + 365*(j-1)
   }
 }
 
 events.fill = c(t(events.all))
   harvests<-events.all[,length(events)]

  events.df<-data.frame(var=rep(c(rep('Q1', length(chem.days)),'P', 'L', 'P'), times = n.years),
                     time = events.fill,
                     value = rep(c(rep(chem.med, length(chem.days)), prawn.n, prawn.l, 0), times = n.years),
                     method = rep(c(rep("add", length(chem.days)), "rep", "rep", "rep"), times = n.years))
  
  output.all = as.data.frame(ode(nstart, time, agro_time, parameters,
                                 events = list(data = events.df)))
    output.all$N = output.all$S + output.all$E + output.all$I
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
  lines(output.all$time[output.all$time >= 200], output.all$L[output.all$time >= 200], col = 'lightblue', lty=2)
  lines(output.all$time, output.all$P.B, lwd=2, col = 'darkgreen')
  legend('topright', legend = c('P', 'P.B', 'L'), lty = 1, col = c('green', 'darkgreen', 'lightblue'), cex = 0.6)
  
  plot(output.all$time, output.all$W, lwd=2, col = 'purple', type = 'l', 
       xlab = 'time', ylab = 'mean worm burden', ylim = c(0,100),
       main = paste('average worm burden/person = ', round(mean(output.all$W), digits=1), sep=''))  
  
  plot(output.all$time, output.all$Q1,lwd=2, col = 'gold2', xlab = 'time', ylab = 'chem concentration (ppb)', type = 'l',
       main = chem)
  
}

#Simulations with different agrochemicals, prawn addition timing ###################
#Run with no agrochemical introducation
agroc_prawn(chem = 'Chlorpyrifos',
            chem.days = 0,
            chem.k = chlor.k,
            chem.med = 0,
            prawn.add = 1,
            prawn.harvest = harvest.time,
            prawn.n = 2*area,
            prawn.l = 25,
            n.years = 3)
#Run with chlorpyrifos
agroc_prawn(chem = 'Chlorpyrifos',
            chem.days = chlor.days,
            chem.k = chlor.k,
            chem.med = med.chlor,
            prawn.add = 1,
            prawn.harvest = harvest.time,
            prawn.n = 2*area,
            prawn.l = 25,
            n.years = 3)
#Delay harvest to avoid agrochemical
agroc_prawn(chem = 'Chlorpyrifos',
            chem.days = chlor.days,
            chem.k = chlor.k,
            chem.med = med.chlor,
            prawn.add = 200,
            prawn.harvest = harvest.time,
            prawn.n = 2*area,
            prawn.l = 25,
            n.years = 3)

#Switch to terbufos, less toxic agrochemical (according to Halstead results) CHANGE TOXICITY MODEL
agroc_prawn(chem = 'Terbufos',
            chem.days = terb.days,
            chem.k = terb.k,
            chem.med = med.terb,
            prawn.add = 1,
            prawn.harvest = harvest.time,
            prawn.n = 2*area,
            prawn.l = 25,
            n.years = 3)

#Switch to esfenvalerate, more toxic agrochemical (according to Halstead results) CHANGE TOXICITY MODEL
agroc_prawn(chem = 'Esfenvalerate',
            chem.days = esfen.days,
            chem.k = esfen.k,
            chem.med = med.esfen,
            prawn.add = 1,
            prawn.harvest = harvest.time,
            prawn.n = 2*area,
            prawn.l = 25,
            n.years = 3)