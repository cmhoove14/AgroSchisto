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
#Packages and source other scripts for model and parameter fitting values #####
require(reshape2)

epss<-read.csv("Elimination_Feasibility/eps_estimates_100bestfit_kap0-3.csv")
  kaps = unique(epss$kap.vec)
  trans = sort(unique(epss$r0.vec))
    trans.lo = sort(unique(epss$r0.vec[epss$lam.vec < 0.00011]))

#get initial idea of relationsip between eps, r0, and kappa
  plot(epss$beta.vec, epss$lam.vec, pch = 16, cex = 0.7)
    plot(epss$lam.vec[epss$lam.vec < 0.00011], epss$r0.vec[epss$lam.vec < 0.00011], 
         pch = 16, cex = 0.7, ylim = c(0,6), xlim = c(4e-5,2e-4))
      points(epss$lam.vec[epss$lam.vec > 0.00011], epss$r0.vec[epss$lam.vec > 0.00011], 
             pch = 16, cex = 0.7, col=2)
  
  plot(epss$r0.vec, epss$eps.vec, pch=16, cex=0.7)
  plot(epss$kap.vec, epss$eps.vec, pch=16, cex=0.7)
  
  plot(epss$r0.vec[epss$kap.vec==0], epss$eps.vec[epss$kap.vec==3], pch=16, cex=0.7)
  plot(epss$lam.vec[epss$kap.vec==0], epss$eps.vec[epss$kap.vec==3], pch=16, cex=0.7)
  plot(epss$beta.vec[epss$kap.vec==0], epss$eps.vec[epss$kap.vec==3], pch=16, cex=0.7)
  
#Reshape data frame and surface plot ######
  epss.rs = reshape(epss,
                    timevar = 'kap.vec',
                    idvar = c('beta.vec', 'lam.vec', 'r0.vec'),
                    direction = 'wide')
  
  
  zm = as.matrix(epss.rs[,c(4:103)])
    persp(x = trans, xlim = range(trans), 
          y = kaps, ylim = range(kaps),
          z = zm, ticktype = 'detailed', nticks = 4, 
          ylab = 'Pos. Density Dependence',
          xlab = 'Transmission Intensity',
          zlab = 'Elimination Feasibility Estimator',
          phi = 30, theta = 0, shade = 0.4, col = 'lightblue')
    
  zm.lamlo = as.matrix(subset(epss.rs, lam.vec < 0.00011)[,c(4:103)])
    persp(x = trans.lo, xlim = range(trans.lo), 
          y = kaps, ylim = range(kaps),
          z = zm.lamlo, ticktype = 'detailed', nticks = 4, 
          ylab = 'Pos. Density Dependence',
          xlab = 'Transmission Intensity',
          zlab = 'Elimination Feasibility Estimator',
          phi = 30, theta = 0, shade = 0.4, col = 'lightblue')
  
  
     