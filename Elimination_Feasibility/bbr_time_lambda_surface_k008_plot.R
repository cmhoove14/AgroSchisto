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

bbrm = read.csv('Elimination_Feasibility/bbr_surface_time_lambda_k008.csv')
  lam.range = seq(1e-4, 3e-4, length.out = 50)  #Range of transmission intensities to test
  mda.years = c(1:20) # annual MDA for 20 years

windows()    
  persp(y = lam.range, ylim = range(lam.range), x = mda.years[-20], xlim = range(mda.years[-20]),
        z = as.matrix(bbrm), zlim = c(-0.5, 1.0), ticktype = 'detailed', nticks = 4, 
        xlab = 'Time (yrs)',
        ylab = 'Transmission Intensity',
        zlab = 'BBR',
        phi = 10, theta = 360-45, shade = 0.4)