#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Function to assess influence of simultaneous chemotoxic/direct effects and trophic/indirect effects ####
load('Review_models/r0_of_Atrazine_ws.RData')
  require(rootSolve)

#3d plot of direct effect (muN) and bottom up effect (phi_N) influence on R0 #################
muN.range = seq(min(conc.atr.means.par[, 5]), max(conc.atr.means.par[, 5]), 
                length.out = 101)
phiN.range = seq(min(conc.atr.means.par[, 6]), max(conc.atr.means.par[, 6]), 
                 length.out = 101)

m.munphin = matrix(nrow = 101,
                   ncol = 101)
for(i in 1:length(muN.range)){
  for(j in 1:length(phiN.range)){
    m.munphin[i,j] = r0.fix(mu_Nqx = muN.range[i],
                            phi_Nqx = phiN.range[j])[3]
  }
}

persp(x = muN.range + parameters['mu_N'], xlim = range(muN.range) + parameters['mu_N'], 
      y = phiN.range * parameters['phi_N']/parameters['A'], ylim = range(phiN.range) * parameters['phi_N']/parameters['A'],
      z = m.munphin - r0.He(0)[3], ticktype = 'detailed', nticks = 4, zlim = c(-3,3),
      xlab = 'Snail mortality rate', ylab = 'Snail carrying capacity', zlab = 'Change in R0',
      phi = 30, theta = 30, shade = 0.4, col = 'lightblue')

#3d plot of direct effect (muN) and top down effect (mu_P) influence on R0 #################
rm(list = ls())
load('Review_models/r0_of_malathion2_ws.RData')

piC.range = seq(min(conc.mal.means.par[, 1]), max(conc.mal.means.par[, 1]), 
                length.out = 101)
muP.range = seq(min(conc.mal.means.par[, 5]), max(conc.mal.means.par[, 5]), 
                length.out = 101)

m.picmup = matrix(nrow = 101,
                  ncol = 101)

for(i in 1:length(piC.range)){
  for(j in 1:length(muP.range)){
    m.picmup[i,j] = r0.fix(pi_Cqx = piC.range[i],
                           mu_Pqx = muP.range[j])[3]
  }
}

persp(x = piC.range, xlim = range(piC.range), 
      y = (muP.range + parameters['mu_P']), ylim = (range(muP.range) + parameters['mu_P']),
      z = m.picmup - r0.In(0)[3], ticktype = 'detailed', nticks = 4, zlim = c(-3,3),
      xlab = 'Cercarial mortality estimate', ylab = 'Predator Mortality estimate', zlab = 'Change in R0',
      phi = 30, theta = 300, shade = 0.4, col = 'lightblue')