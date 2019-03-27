source("Agrochemical_Review/Models/initial_parameters.R")
source("Agrochemical_Review/Models/model_helper_functions.R")
source("Agrochemical_Review/Models/r0_functions.R")

#Do some tests #############  
parameters = init_pars

#Influence of variable predator densities on r0 #########
pred.dens = seq(0,1,0.01)
  r0.pd = as.numeric()
  n.pd = as.numeric()

for(i in 1:length(pred.dens)){
  parameters['K_P'] = pred.dens[i]*area
  r0.pd[i] = r0.Ag(parameters = parameters)[3]
  n.pd[i] = r0.Ag(parameters = parameters)[1]
}
 
  plot(pred.dens, r0.pd, type = 'l', lwd = 2, xlab = 'pred CC', ylab = 'R0')
  plot(pred.dens, n.pd, type = 'l', lwd = 2, col = 3, xlab = 'pred CC', ylab = 'N*')
  
  parameters['K_P'] = 0.125*area
  
#influence of various mortality rates on r0 ###############
mup.dens = seq(0.01,1,0.01) - parameters['mu_P']
r0.mup = as.numeric()

for(i in 1:length(mup.dens)){
  parameters['mu_P'] = mup.dens[i]+0.038
  r0.mup[i] = r0.Ag(parameters = parameters)[3]
}

plot(mup.dens, r0.mup, type = 'l', lwd = 2, xlab = 'pred mortality rate', ylab = 'R0')

parameters['mu_P'] = 0.038 #Reset mortality rate to original

#influence of various snail carrying capacities on r0 ###############
Kn.dens = seq(10, 200, 10)*area
r0.Kn = as.numeric()

for(i in 1:length(Kn.dens)){
  parameters['K_N'] = Kn.dens[i]
  r0.Kn[i] = r0.Ag(parameters = parameters)[3]
}

plot(Kn.dens, r0.Kn, type = 'l', lwd = 2, xlab = 'snail carrying capacity', ylab = 'R0')

parameters['K_N'] = 50*area #Reset carrying capacity to original 

#influence of various snail reproductive rates on r0 ###############
fn.dens = seq(parameters['f_N'], 0, length.out = 50)
  r0.fn = as.numeric()

  for(i in 1:length(fn.dens)){
    parameters['f_N'] = fn.dens[i]
    r0.fn[i] = r0.Ag(parameters = parameters)[3]
  }

plot(fn.dens, r0.fn, type = 'l', lwd = 2, xlab = 'snail reproductive rate', ylab = 'R0')

parameters['f_N'] = 0.1 #Reset reproductive rate to original

#influence of various snail mortality rates on r0 ###############
mun.dens = seq(0.0001, 1, length.out = 100)
r0.mun = as.numeric()

for(i in 1:length(mun.dens)){
  parameters['mu_N'] = mun.dens[i]
  r0.mun[i] = r0.Ag(parameters = parameters)[3]
}

plot(mun.dens, r0.mun, type = 'l', lwd = 2, xlab = 'snail mortality rate', ylab = 'R0')

parameters['mu_N'] = 0.017 #Reset reproductive rate to original


#influence of cercarial hours exposure on r0 ###############
pic.dens = seq(parameters['pi_C']/10, parameters['pi_C'], length.out = 100)
r0.pic = as.numeric()

for(i in 1:length(pic.dens)){
  parameters['pi_C'] = pic.dens[i]
  r0.pic[i] = r0.Ag(parameters = parameters)[3]
}

plot(pic.dens, r0.pic, type = 'l', lwd = 2, xlab = 'cercariae hours exposure per day', ylab = 'R0')

parameters['mu_N'] = 14.21 #Reset cercariae hours exposure to original
