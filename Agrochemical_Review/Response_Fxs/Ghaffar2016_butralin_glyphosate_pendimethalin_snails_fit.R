#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Ghaffar 2016 SNAIL (B. alexandrina) data
#Snail (30 B. alexandrina 6-8mm shell width) toxicity ##########
#Butralin ##############
but.dat = data.frame(lcs = c(0, 10, 25, 50, 90),
                     butralin = c(556, 2417, 3906, 5560, 8703))
  but.dat.no0 = subset(but.dat, lcs!=0)
  but.dat.no0$mort = but.dat.no0$lcs/100
  but.dat.no0$probit = qnorm(but.dat.no0$lcs/100)
  but.dat.no0$ppm = but.dat.no0$butralin/1000
  but.dat.no0$ppmlog10 = log10(but.dat.no0$butralin/1000)

#Create function based on LC values
  gaf_but_lin <- lm(probit~ppmlog10, data = but.dat.no0)
    gaf.b.lin <- coef(gaf_but_lin)[1]
    gaf.se.b.lin <- summary(gaf_but_lin)$coef[2,2]
      gaf.b.lin.up <- confint(gaf_but_lin)[1,1]
      gaf.b.lin.lo <- confint(gaf_but_lin)[1,2]
    gaf.m.lin <- coef(gaf_but_lin)[2]
  
muNq_but_Gaf16_lin <- function(He, b = gaf.b.lin){
  mun = pnorm(b + log10(He/1000)*gaf.m.lin)
  return(mun)
}  

mu_Nq_butr_gaf16_uncertainty = function(He){
    heu = (He/1000)
    lc50 = rnorm(1, gaf.b.lin, gaf.se.b.lin)
    mun = pnorm(lc50 + log10(He/1000)*gaf.m.lin)
  
    return(mun)
}
  
#keep vector
keep.gaf.but = c('mu_Nq_butr_gaf16_uncertainty',
                 'gaf.b.lin', 'gaf.se.b.lin', 'gaf.m.lin')    
        
#Glyphosate #############
gly.dat = data.frame(lcs = c(0, 0, 10, 25, 50, 90),
                     glyphosate = c(0, 1506, 3875, 9174, 15062, 26249))

  gly.dat.no0 = subset(gly.dat, lcs!=0)
  gly.dat.no0$probit = qnorm(gly.dat.no0$lcs/100)
  gly.dat.no0$ppm = gly.dat.no0$glyphosate/1000
  gly.dat.no0$ppmlog10 = log10(gly.dat.no0$glyphosate/1000)
  
#Create function based on LC values
  gaf_gly_lin <- lm(probit~ppmlog10, data = gly.dat.no0)
    gaf.b.lin.gly <- coef(gaf_gly_lin)[1]
    gaf.se.b.lin.gly <- summary(gaf_gly_lin)$coef[2,2]
      gaf.b.lin.up.gly <- confint(gaf_gly_lin)[1,1]
      gaf.b.lin.lo.gly <- confint(gaf_gly_lin)[1,2]
    gaf.m.lin.gly <- coef(gaf_gly_lin)[2]
  
muNq_gly_Gaf16_lin <- function(He, b = gaf.b.lin.gly){
  mun = pnorm(b + log10(He/1000)*gaf.m.lin.gly)
  return(mun)
}  

mu_Nq_gly_gaf16_uncertainty = function(He){
  heu = (He/1000)
  lc50 = rnorm(1, gaf.b.lin.gly, gaf.se.b.lin.gly)
  mun = pnorm(lc50 + log10(He/1000)*gaf.m.lin.gly)
    
  return(mun)
}

keep.gaf.gly = c('mu_Nq_gly_gaf16_uncertainty', 
                 'gaf.b.lin.gly', 'gaf.se.b.lin.gly', 'gaf.m.lin.gly')    
    
#Pendimethalin ##############
pen.dat = data.frame(lcs = c(0, 0, 10, 25, 50, 90),
                     pendimethalin = c(0, 214.8, 535, 1299, 2148, 3762))
    
pen.dat.no0 = subset(pen.dat, lcs!=0)
  pen.dat.no0$probit = qnorm(pen.dat.no0$lcs/100)
  pen.dat.no0$ppm = pen.dat.no0$pendimethalin/1000
  pen.dat.no0$ppmlog10 = log10(pen.dat.no0$pendimethalin/1000)
    
#Create function based on LC values
  gaf_pen_lin <- lm(probit~ppmlog10, data = pen.dat.no0)
    gaf.b.lin.pen <- coef(gaf_pen_lin)[1]
    gaf.se.b.lin.pen <- summary(gaf_pen_lin)$coef[2,2]
      gaf.b.lin.up.pen <- confint(gaf_pen_lin)[1,1]
      gaf.b.lin.lo.pen <- confint(gaf_pen_lin)[1,2]
    gaf.m.lin.pen <- coef(gaf_pen_lin)[2]
  
muNq_pen_Gaf16_lin <- function(He, b = gaf.b.lin.pen){
  mun = pnorm(b + log10(He/1000)*gaf.m.lin.pen)
  return(mun)
}  


mu_Nq_pen_gaf16_uncertainty = function(He){
  heu = (He/1000)
  lc50 = rnorm(1, gaf.b.lin.pen, gaf.se.b.lin.pen)
  mun = pnorm(lc50 + log10(He/1000)*gaf.m.lin.pen)
    
  return(mun)
}

keep.gaf.pen = c('mu_Nq_pen_gaf16_uncertainty',
                 'gaf.b.lin.pen', 'gaf.se.b.lin.pen', 'gaf.m.lin.pen')    
  