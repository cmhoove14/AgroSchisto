#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Gustafson et al 2016 atrazine influence on Physa Acuta cercarial production 
gust_dat <- data.frame(atr = c(0,3,30), 
                       logatr = log(c(0,3,30)+1),
                       cerc = c(792.8, 536.3, 236.9),
                       sd = c(108.2, 124.9, 153.1))

gust_theta_mod <- lm(cerc ~ logatr, weights = 1/sd, data = gust_dat)
gust_int <- coef(gust_theta_mod)[1]

gust_theta_pred <- predict(gust_theta_mod, newdata = data.frame(logatr = log(c(0:30)+1)), interval = "confidence", level = 0.95)

gust_theta_uncertainty <- function(He){
  heu <- log(He+1)
  theta.init <- predict(gust_theta_mod, newdata = data.frame(logatr = heu), se.fit = TRUE)
  theta <- rnorm(1, theta.init[[1]], theta.init[[2]]) / gust_int
  
  if(theta < 0){
    return(0)
  } else {

   return(theta)
   
  }
  
}

keep.gust2016 <- c("gust_theta_uncertainty", "gust_theta_mod", "gust_int")
