#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

##################################################################################################

#50 simulations with 0<kap<2 and 1.2e-4<lam<3.7e-4 #############
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/eps_surface_lamRange_kRange.Rdata")
load("Elimination_Feasibility/Organize/Models/Outputs_Refs/model_fit_profile_likelihood_parameters.Rdata")

  sim.range = 50
  lam.range = seq(min(fin_pars95ci$lamda_twa), max(fin_pars95ci$lamda_twa), length.out = sim.range)  #transmission intensity
  kap.range = seq(0, 2, length.out = sim.range)     #Pos. density dependence

windows(width = 13, height = 13.2)

persp(y = lam.range, ylim = range(lam.range), x = kap.range, xlim = range(kap.range),
      z = eps.fill, zlim = c(-0.05, 0.05), ticktype = 'detailed', nticks = 4, 
      xlab = 'Pos. Density Dependence',
      ylab = 'Transmission Intensity',
      zlab = 'Elimination Feasibility Estimator',
      phi = 18, theta = 45, shade = 0.4, cex.lab = 1.1)

#Save as tiff ###########
tiff("Elimination_Feasibility/plots/PLoS_Figs_PostReview/Fig5.tiff", 
     height = 13.2, width = 13.2, units = 'cm', compression = "lzw", res = 300) 

persp(y = lam.range, ylim = range(lam.range), x = kap.range, xlim = range(kap.range),
      z = eps.fill, zlim = c(-0.05, 0.05), ticktype = 'detailed', nticks = 4, 
      xlab = 'Pos. Density Dependence',
      ylab = 'Transmission Intensity',
      zlab = 'Elimination Feasibility Estimator',
      phi = 18, theta = 45, shade = 0.4, cex.lab = 0.75, cex.axis = 0.5)

dev.off()

rm(list = ls())
