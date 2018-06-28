#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############
source('Agrochemical_Review/Models/mod_q.R')

ac.sim = as.data.frame(ode(ac.start, ac.time, agrochem_mod, ac.pars))
  ac.sim$N = ac.sim$S.A + ac.sim$E.A + ac.sim$I.A
  ac.sim$Wf = ac.sim$W * sapply(ac.sim$W, mateprob, k=ac.pars['k']) * 0.5
  
  ac.eqbm = ac.sim[dim(ac.sim)[1],]
    
ac.sim %>% ggplot(aes(x = time)) +
  theme_bw() +
  geom_line(aes(y = N), size = 1.25, col = 'black') +
  geom_line(aes(y = S.A), size = 1.25, col = 'green') +
  geom_line(aes(y = E.A), size = 1.25, col = 'orange') +
  geom_line(aes(y = I.A), size = 1.25, col = 'red') 

ggplot(ac.sim, aes(x = time)) +
  theme_bw() +
  geom_line(aes(y = P.A), size = 1.25, col = 'blue') 

ggplot(ac.sim, aes(x = time)) +
  theme_bw() +
  geom_line(aes(y = W), size = 1.25, col = 'purple') 
