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
require(ggplot2)

source('Review_models/r0_of_q.R')

#Assess sensitivity of model results to agrochemically influenced parameters
r0.sens = data.frame(var = c('Snail birth rate', 'Snail carrying capacity',
                             'Snail mortality rate', 'Predator mortality rate',
                             'Predator attack rate', 'Cercarial shedding rate',
                             'Miracidial survival', 'Cercarial survival',
                             'Schistosome egg viability'),
                     high = NA,
                     point = rep(r0.fix()[3], 9),
                     low = NA,
                     rel.hi = NA,
                     rel.pt = rep(0, 9),
                     rel.lo = NA)

r0.sens$high = c(r0.fix(f_Nqx = 2)[3],
                 r0.fix(phi_Nqx = 2)[3],
                 r0.fix(mu_Nqx = parameters['mu_N'])[3],
                 r0.fix(mu_Pqx = parameters['mu_P'])[3],
                 r0.fix(alpha_qx = 2)[3],
                 r0.fix(theta_qx = 2)[3],
                 r0.fix(pi_Mqx = 2)[3],
                 r0.fix(pi_Cqx = 2)[3],
                 r0.fix(v_qx = 2)[3])

r0.sens$low = c(r0.fix(f_Nqx = 0.5)[3],
                 r0.fix(phi_Nqx = 0.5)[3],
                 r0.fix(mu_Nqx = parameters['mu_N'] * -0.5)[3],
                 r0.fix(mu_Pqx = parameters['mu_P'] * -0.5)[3],
                 r0.fix(alpha_qx = 0.5)[3],
                 r0.fix(theta_qx = 0.5)[3],
                 r0.fix(pi_Mqx = 0.5)[3],
                 r0.fix(pi_Cqx = 0.5)[3],
                 r0.fix(v_qx = 0.5)[3])

r0.sens$rel.lo = c(r0.fix(f_Nqx = 0.5)[3] / r0.sens$point[1] - 1,
                    r0.fix(phi_Nqx = 0.5)[3] / r0.sens$point[1] - 1,
                    r0.fix(mu_Nqx = parameters['mu_N'] * -0.5)[3] / r0.sens$point[1] - 1,
                    r0.fix(mu_Pqx = parameters['mu_P'] * -0.5)[3] / r0.sens$point[1] - 1,
                    r0.fix(alpha_qx = 0.5)[3] / r0.sens$point[1] - 1,
                    r0.fix(theta_qx = 0.5)[3] / r0.sens$point[1] - 1,
                    r0.fix(pi_Mqx = 0.5)[3] / r0.sens$point[1] - 1,
                    r0.fix(pi_Cqx = 0.5)[3] / r0.sens$point[1] - 1,
                    r0.fix(v_qx = 0.5)[3] / r0.sens$point[1] - 1)

r0.sens$rel.hi = c(r0.fix(f_Nqx = 2)[3] / r0.sens$point[1] - 1,
                   r0.fix(phi_Nqx = 2)[3] / r0.sens$point[1] - 1,
                   r0.fix(mu_Nqx = parameters['mu_N'])[3] / r0.sens$point[1] - 1,
                   r0.fix(mu_Pqx = parameters['mu_P'])[3] / r0.sens$point[1] - 1,
                   r0.fix(alpha_qx = 2)[3] / r0.sens$point[1] - 1,
                   r0.fix(theta_qx = 2)[3] / r0.sens$point[1] - 1,
                   r0.fix(pi_Mqx = 2)[3] / r0.sens$point[1] - 1,
                   r0.fix(pi_Cqx = 2)[3] / r0.sens$point[1] - 1,
                   r0.fix(v_qx = 2)[3] / r0.sens$point[1] - 1)

#reorder data frame for compatibility with axis labels
r0.sens$var = factor(r0.sens$var, levels = as.character(r0.sens$var[c(order(r0.sens$var, decreasing = T))])) 

par.labs = c(expression(paste('Snail mortality rate (', italic(mu[N]),')')), #axis labels
             expression(paste('Snail carrying capacity (', italic(Phi[N]),')')),
             expression(paste('Snail recruitment rate (', italic(f[N]),')')),
             expression(paste('Schistosome egg viability (', italic(v),')')),
             expression(paste('Predator mortality rate (', italic(mu[P]),')')),
             expression(paste('Predator attack rate (', italic(alpha),')')),
             expression(paste('Miracidial survival (', italic(pi[M]),')')),
             expression(paste('Cercarial survival (', italic(pi[C]),')')),
             expression(paste('Cercarial shedding rate of infected snails (', italic(theta),')')))

sens.bar = ggplot(r0.sens, aes(x = var, y = rel.hi)) +
            theme_bw()+
            scale_y_continuous(limits = c(-1,1), breaks = seq(-1, 1, 0.25))+
            scale_x_discrete(labels = par.labs) +
            coord_flip()+
            geom_bar(fill = 'blue', stat = 'identity', width = 0.25)+
            geom_bar(data = r0.sens, aes(x = var, y = rel.lo), 
                     fill = 'red', stat = 'identity', width = 0.25)+
            #geom_hline(xintercept = 0, lty = 2)+
            labs(title = expression(paste('Sensitivity of R'[0], ' to Agrochemically Influenced Parameters')), 
                 y = expression(paste(Delta, 'R'[0], ' (%)'), sep = ''), 
                 x = 'Parameter')