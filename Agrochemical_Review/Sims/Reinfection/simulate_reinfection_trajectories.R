source('Agrochemical_Review/Models/models.R')
require(tidyverse)

ac.sim = as.data.frame(ode(ac.start, ac.time, agrochem_mod, ac.pars))
  ac.sim$N = ac.sim$S.A + ac.sim$E.A + ac.sim$I.A
  ac.sim$W = (1-cvrg)*ac.sim$Wu + cvrg*ac.sim$Wt
  ac.sim$Wf = ac.sim$W * sapply(ac.sim$W, mateprob, k=ac.pars['k']) * 0.5
  
  ac.eqbm = ac.sim[dim(ac.sim)[1],]

base.start = c(S = ac.eqbm$S.A,
               E = ac.eqbm$E.A,
               I = ac.eqbm$I.A,
               Wu = ac.eqbm$Wu,
               Wt = ac.eqbm$Wt,
               P = ac.eqbm$P.A)

base.time = c(0:(365*30))

mda.eff = 0.92

mda.event = data.frame(var = "Wt",
                       time = 30,
                       value = 1-mda.eff,
                       method = "mult")
  
base.sim = as.data.frame(ode(base.start, base.time, agrochem_mod, ac.pars,
                             events = list(data = mda.event)))
  base.sim$N = base.sim$S + base.sim$E + base.sim$I
  base.sim$W = (1-cvrg)*base.sim$Wu + cvrg*base.sim$Wt
  base.sim$Wf = base.sim$W * sapply(base.sim$W, mateprob, k=ac.pars['k']) * 0.5

base.sim %>% ggplot(aes(x = time)) +
  theme_bw() +
  geom_line(aes(y = Wt), size = 1.25, lty = 2, col = 'purple') +
  geom_line(aes(y = Wu), size = 1.25, lty = 3, col = 'purple') +
  geom_line(aes(y = W), size = 1.25, col = 'purple')
  
  