#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
require(rootSolve)
require(deSolve)
source('Review_models/r0_of_q.R')

#R0(q) code to investigate influence of chlorpyrifos

#Chlorpyrifos functional responses ##############
  source('T_variant_Agro/ChlorP_Tox.R')


parameters['phi_P'] = 0.1*area  
  
r0.conc = matrix(ncol = 101, nrow = 12)

for(i in 1:100){
  r0.conc[1,i] = r0.q(conc = (i-1))[3]                       #Baseline r0 with no agrochemical effects included
  r0.conc[2,i] = r0.q(conc = (i-1),
                       f.f_Nq = f_Nq_chlor_ibrahim92)[3]      #Toxic effects on snail reproduction
  r0.conc[3,i] = r0.q(conc = (i-1),
                       f.mu_Pq = muPq_ch_sat09)[3]            #Mortality effects on predator population from Sat'09
  r0.conc[4,i] = r0.q(conc = (i-1),
                       f.alpha_q = alpha_q_ch_sat09)[3]       #Feeding reduction effects on predator population
  r0.conc[5,i] = r0.q(conc = (i-1),
                       f.mu_Nq = muNq_chlor_ibrahim92)[3]     #Mortality effects on snail population
  r0.conc[6,i] = r0.q(conc = (i-1),
                       f.mu_Nq = muNq_chlor_ibrahim92,
                       f.f_Nq = f_Nq_chlor_ibrahim92)[3]      #Both repro and mortality effects in snails
  r0.conc[7,i] = r0.q(conc = (i-1),
                       f.mu_Pq = muPq_ch_sat09,
                       f.alpha_q = alpha_q_ch_sat09)[3]       #Both attack rate and mortality effects in preds
  r0.conc[8,i] = r0.q(conc = (i-1),
                       f.mu_Pq = muPq_ch_Halstead)[3]         #Mortality effects on predator population from Halstead'15
  r0.conc[9,i] = r0.q(conc = (i-1),
                       f.pi_Mq = pi_Mq_ch_has11)[3]           #Mortality effects on miracidia from Hasheesh11
  r0.conc[10,i] = r0.q(conc = (i-1),
                       f.pi_Cq = pi_Cq_ch_has11)[3]           #Mortality effects on cercariae from Hasheesh11
  r0.conc[11,i] = r0.q(conc = (i-1),
                        f.pi_Mq = pi_Mq_ch_has11,
                        f.pi_Cq = pi_Cq_ch_has11)[3]          #Mortality effects on miracidia & cercariae from Hasheesh11
  r0.conc[12,i] = r0.q(conc = (i-1),
                        f.f_Nq = f_Nq_chlor_ibrahim92,
                        f.mu_Pq = muPq_ch_Halstead,
                        f.mu_Nq = muNq_chlor_ibrahim92,
                        f.pi_Mq = pi_Mq_ch_has11,
                        f.pi_Cq = pi_Cq_ch_has11)[3]          #Combined effects
}

  plot(x=c(0:100), y=r0.conc[1,], type = 'l', lwd = 3, ylim = c(0,6), col = 'grey40',
       ylab = expression('R'[0]), xlab = 'Chlorpyrifos Concentration (ppb)',
       main = expression(paste('R'[0],'(q) simulations with Chlorpyrifos', sep='')))
    legend('bottomleft', legend = expression('baseline R'[0]), lwd=3, col = 'grey40', cex=0.75)
 #Snail effects
  lines(x=c(0:100), y=r0.conc[2,], lwd=2, col='lightblue')
  lines(x=c(0:100), y=r0.conc[5,], lwd=2, col='navy')
  lines(x=c(0:100), y=r0.conc[6,], lwd=2, col='blue', lty=2)
    legend('bottomleft', legend = c(expression('baseline R'[0]),
                                    expression(paste('f'[N],'(q)', sep='')),
                                    expression(paste(mu[N],'(q)', sep='')),
                                    'combined'), 
           lwd=c(3,2,2,2), lty = c(1,1,1,2), 
           col = c('grey40', 'lightblue', 'navy', 'blue'), cex=0.75)
  
 #Predator effects 
  lines(x=c(0:100), y=r0.conc[3,], lwd=2, col='olivedrab1')
  lines(x=c(0:100), y=r0.conc[4,], lwd=2, col='olivedrab4')
  lines(x=c(0:100), y=r0.conc[7,], lwd=2, col='palegreen', lty=2)
  lines(x=c(0:100), y=r0.conc[8,], lwd=2, col='green')
    legend('bottomleft', legend = c(expression('baseline R'[0]),
                                    expression(paste(mu[P],'(q) Halstead', sep='')),
                                    expression(paste(alpha,'(q) Sat.', sep='')),
                                    expression(paste(mu[P],'(q) Sat.', sep='')),
                                    'combined'), 
           lwd=c(3,2,2,2,2), lty = c(1,1,1,1,2), 
           col = c('grey40', 'green', 'olivedrab4', 'olivedrab1', 'palegreen'), cex=0.75)
 #Larval effects
  lines(x=c(0:100), y=r0.conc[9,], lwd=2, col='indianred4')
  lines(x=c(0:100), y=r0.conc[10,], lwd=2, col='lightcoral')
  lines(x=c(0:100), y=r0.conc[11,], lwd=2, col='red', lty=2)
    legend('bottomleft', legend = c(expression('baseline R'[0]),
                                    expression(paste(pi[M],'(q)', sep='')),
                                    expression(paste(pi[C],'(q)', sep='')),
                                    'combined'), 
           lwd=c(3,2,2,2), lty = c(1,1,1,2), 
           col = c('grey40', 'indianred4', 'lightcoral', 'red'), cex=0.75)
 #All effects combined
  lines(x=c(0:100), y=r0.conc[12,], lwd=3, col='black', lty=3)
    legend('bottomleft', legend = c(expression('baseline R'[0]),
                                    expression(paste('R'[0],'(q)', sep=''))), 
           lwd=c(3,3), lty = c(1,3), col = c('grey40', 'black'), cex=0.75)
    