#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Data extraction and model fitting to Bajet et al 2012 uinsecticide toxicity to Macrobrachium lar##########
#Paper reports LC50, LC10, NOEC, and LOEC values. To obtain d-r function, slope is estimated
#by linear regression fit to the LC50 and LC10 values and incorporated with reported LC50 and its uncertainty
#All data points are then used in plots for qualitative validation of the resulting function

get_m <- function(lc50, lc10){
  (qnorm(0.5) - qnorm(0.1))/(log10(lc50) - log10(lc10))
}

#Lambda cyhalothrin ##############
noec.lamcy.baj <- 0.005
lc10.lamcy.baj <- 0.009
lc50.lamcy.baj <- 0.012
lc50.lamcy.baj.up <- 0.015
lc50.lamcy.baj.lo <- 0.009

lc50.lamcy.baj.se <- mean(c(log10(lc50.lamcy.baj.up/lc50.lamcy.baj), 
                            log10(lc50.lamcy.baj/lc50.lamcy.baj.lo))) / 1.96
  
baj.lamcy.m.lin <- get_m(lc50.lamcy.baj, lc10.lamcy.baj)
  
muPq_lamcy_Bajet12_uncertainty = function(In){
    lc50 = 10^rnorm(1, log10(lc50.lamcy.baj), lc50.lamcy.baj.se)
    mun = pnorm(baj.lamcy.m.lin * log10(In/lc50))
  
    return(mun)
}

plot(c(noec.lamcy.baj, lc10.lamcy.baj, lc50.lamcy.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.lamcy.baj*2), ylim = c(0,1), pch = 16,
     xlab = "Lambda-Cyhalothrin (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.lamcy.baj.lo, x1 = lc50.lamcy.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.lamcy.baj*2,lc50.lamcy.baj/20), 
         sapply(seq(0,lc50.lamcy.baj*2,lc50.lamcy.baj/20), muPq_lamcy_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

#keep vector
keep.Bajet.lamcy = c('muPq_lamcy_Bajet12_uncertainty', 'baj.lamcy.m.lin', 'lc50.lamcy.baj', 'lc50.lamcy.baj.se')    

#Deltamethrin ##############
noec.delt.baj <- 0.01
lc10.delt.baj <- 0.034
lc50.delt.baj <- 0.26
lc50.delt.baj.up <- 0.41
lc50.delt.baj.lo <- 0.16

lc50.delt.baj.se <- mean(c(log10(lc50.delt.baj.up/lc50.delt.baj), 
                            log10(lc50.delt.baj/lc50.delt.baj.lo))) / 1.96
  
baj.delt.m.lin <- get_m(lc50.delt.baj, lc10.delt.baj)
  
muPq_deltamethrin_Bajet12_uncertainty = function(In){
    lc50 = 10^rnorm(1, log10(lc50.delt.baj), lc50.delt.baj.se)
    mun = pnorm(baj.delt.m.lin * log10(In/lc50))
  
    return(mun)
}

plot(c(noec.delt.baj, lc10.delt.baj, lc50.delt.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.delt.baj*2), ylim = c(0,1), pch = 16,
     xlab = "Deltamethrin (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.delt.baj.lo, x1 = lc50.delt.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.delt.baj*2,lc50.delt.baj/20), 
         sapply(seq(0,lc50.delt.baj*2,lc50.delt.baj/20), muPq_deltamethrin_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

  
#keep vector
keep.Bajet.delt = c('muPq_deltamethrin_Bajet12_uncertainty', 'baj.delt.m.lin', 'lc50.delt.baj', 'lc50.delt.baj.se')    

#Cypermethrin ##############
noec.cyper.baj <- 0.2
lc10.cyper.baj <- 0.42
lc50.cyper.baj <- 1.33
lc50.cyper.baj.up <- 1.89
lc50.cyper.baj.lo <- 0.92

lc50.cyper.baj.se <- mean(c(log10(lc50.cyper.baj.up/lc50.cyper.baj), 
                            log10(lc50.cyper.baj/lc50.cyper.baj.lo))) / 1.96
  
baj.cyper.m.lin <- get_m(lc50.cyper.baj, lc10.cyper.baj)
  
muPq_cypermethrin_Bajet12_uncertainty = function(In){
    lc50 = 10^rnorm(1, log10(lc50.cyper.baj), lc50.cyper.baj.se)
    mun = pnorm(baj.cyper.m.lin * log10(In/lc50))
  
    return(mun)
}

plot(c(noec.cyper.baj, lc10.cyper.baj, lc50.cyper.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.cyper.baj*2), ylim = c(0,1), pch = 16,
     xlab = "cypermethrin (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.cyper.baj.lo, x1 = lc50.cyper.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.cyper.baj*2,lc50.cyper.baj/20), 
         sapply(seq(0,lc50.cyper.baj*2,lc50.cyper.baj/20), muPq_cypermethrin_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

#keep vector
keep.Bajet.cyper = c('muPq_cypermethrin_Bajet12_uncertainty', 'baj.cyper.m.lin', 'lc50.cyper.baj', 'lc50.cyper.baj.se')    

#Chlorpyrifos ##############
noec.chlor.baj <- 0.05
lc10.chlor.baj <- 0.11
lc50.chlor.baj <- 1.71
lc50.chlor.baj.up <- 2.16
lc50.chlor.baj.lo <- 1.36

lc50.chlor.baj.se <- mean(c(log10(lc50.chlor.baj.up/lc50.chlor.baj), 
                            log10(lc50.chlor.baj/lc50.chlor.baj.lo))) / 1.96
  
baj.chlor.m.lin <- get_m(lc50.chlor.baj, lc10.chlor.baj)
  
muPq_chlorpyrifos_Bajet12_uncertainty = function(In){
    lc50 = 10^rnorm(1, log10(lc50.chlor.baj), lc50.chlor.baj.se)
    mun = pnorm(baj.chlor.m.lin * log10(In/lc50))
  
    return(mun)
}

plot(c(noec.chlor.baj, lc10.chlor.baj, lc50.chlor.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.chlor.baj*2), ylim = c(0,1), pch = 16,
     xlab = "chlorpyrifos (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.chlor.baj.lo, x1 = lc50.chlor.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.chlor.baj*2,lc50.chlor.baj/20), 
         sapply(seq(0,lc50.chlor.baj*2,lc50.chlor.baj/20), muPq_chlorpyrifos_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

#keep vector
keep.Bajet.chlor = c('muPq_chlorpyrifos_Bajet12_uncertainty', 'baj.chlor.m.lin', 'lc50.chlor.baj', 'lc50.chlor.baj.se')    

#Profenofos ##############
noec.prof.baj <- 1
lc10.prof.baj <- 1.29
lc50.prof.baj <- 3.77
lc50.prof.baj.up <- 4.8
lc50.prof.baj.lo <- 2.8

lc50.prof.baj.se <- mean(c(log10(lc50.prof.baj.up/lc50.prof.baj), 
                            log10(lc50.prof.baj/lc50.prof.baj.lo))) / 1.96
  
baj.prof.m.lin <- get_m(lc50.prof.baj, lc10.prof.baj)
  
muPq_profenofos_Bajet12_uncertainty = function(In){
    lc50 = 10^rnorm(1, log10(lc50.prof.baj), lc50.prof.baj.se)
    mun = pnorm(baj.prof.m.lin * log10(In/lc50))
  
    return(mun)
}

plot(c(noec.prof.baj, lc10.prof.baj, lc50.prof.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.prof.baj*2), ylim = c(0,1), pch = 16,
     xlab = "profenofos (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.prof.baj.lo, x1 = lc50.prof.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.prof.baj*2,lc50.prof.baj/20), 
         sapply(seq(0,lc50.prof.baj*2,lc50.prof.baj/20), muPq_profenofos_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

#keep vector
keep.Bajet.prof = c('muPq_profenofos_Bajet12_uncertainty', 'baj.prof.m.lin', 'lc50.prof.baj', 'lc50.prof.baj.se')    

#Malathion ##############
noec.mal.baj <- 100
lc10.mal.baj <- 614
lc50.mal.baj <- 1484
lc50.mal.baj.up <- 1983
lc50.mal.baj.lo <- 1103

lc50.mal.baj.se <- mean(c(log10(lc50.mal.baj.up/lc50.mal.baj), 
                        log10(lc50.mal.baj/lc50.mal.baj.lo))) / 1.96
  
baj.mal.m.lin <- get_m(lc50.mal.baj, lc10.mal.baj)
  
muPq_malathion_Bajet12_uncertainty = function(In){
    lc50 = 10^rnorm(1, log10(lc50.mal.baj), lc50.mal.baj.se)
    mun = pnorm(baj.mal.m.lin * log10(In/lc50))
  
    return(mun)
}

plot(c(noec.mal.baj, lc10.mal.baj, lc50.mal.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.mal.baj*2), ylim = c(0,1), pch = 16,
     xlab = "malathion (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.mal.baj.lo, x1 = lc50.mal.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.mal.baj*2,lc50.mal.baj/20), 
         sapply(seq(0,lc50.mal.baj*2,lc50.mal.baj/20), muPq_malathion_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

#keep vector
keep.Bajet.mal = c('muPq_malathion_Bajet12_uncertainty', 'baj.mal.m.lin', 'lc50.mal.baj', 'lc50.mal.baj.se')    

#Carbaryl ##############
noec.carb.baj <- 5
lc10.carb.baj <- 18.2
lc50.carb.baj <- 48.6
lc50.carb.baj.up <- 77.1
lc50.carb.baj.lo <- 30.7

lc50.carb.baj.se <- mean(c(log10(lc50.carb.baj.up/lc50.carb.baj), 
                        log10(lc50.carb.baj/lc50.carb.baj.lo))) / 1.96
  
baj.carb.m.lin <- get_m(lc50.carb.baj, lc10.carb.baj)
  
muPq_carbaryl_Bajet12_uncertainty = function(In){
    lc50 = 10^rnorm(1, log10(lc50.carb.baj), lc50.carb.baj.se)
    mun = pnorm(baj.carb.m.lin * log10(In/lc50))
  
    return(mun)
}

plot(c(noec.carb.baj, lc10.carb.baj, lc50.carb.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.carb.baj*2), ylim = c(0,1), pch = 16,
     xlab = "carbaryl (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.carb.baj.lo, x1 = lc50.carb.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.carb.baj*2,lc50.carb.baj/20), 
         sapply(seq(0,lc50.carb.baj*2,lc50.carb.baj/20), muPq_carbaryl_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

#keep vector
keep.Bajet.carb = c('muPq_carbaryl_Bajet12_uncertainty', 'baj.carb.m.lin', 'lc50.carb.baj', 'lc50.carb.baj.se')    

#2,4-D ##############
# No toxicity observed, therefore function returns 0
muPq_24d_Bajet12_uncertainty <- function(He){
  return(0)
}

#Butachlor ##############
noec.but.baj <- 6000
lc10.but.baj <- 4439
lc50.but.baj <- 22390
lc50.but.baj.up <- 25560
lc50.but.baj.lo <- 19610

lc50.but.baj.se <- mean(c(log10(lc50.but.baj.up/lc50.but.baj), 
                        log10(lc50.but.baj/lc50.but.baj.lo))) / 1.96
  
baj.but.m.lin <- get_m(lc50.but.baj, lc10.but.baj)
  
muPq_butachlor_Bajet12_uncertainty = function(He){
    lc50 = 10^rnorm(1, log10(lc50.but.baj), lc50.but.baj.se)
    mun = pnorm(baj.but.m.lin * log10(He/lc50))
  
    return(mun)
}

plot(c(noec.but.baj, lc10.but.baj, lc50.but.baj), c(0,10,50)/100, 
     xlim = c(0, lc50.but.baj*2), ylim = c(0,1), pch = 16,
     xlab = "butachlor (ppb)", ylab = "mortality",
     main = "Bajet et al 2012 Macrobrachium lar")
  segments(x0 = lc50.but.baj.lo, x1 = lc50.but.baj.up, y0 = 0.5, y1 = 0.5)
  
  set.seed(43093)
  
  points(seq(0,lc50.but.baj*2,lc50.but.baj/20), 
         sapply(seq(0,lc50.but.baj*2,lc50.but.baj/20), muPq_butachlor_Bajet12_uncertainty),
         pch = 5, col = 4, cex = 0.5)
  legend("bottomright", bty="n", pch = c(16, 5), col = c(1,4), 
         legend = c("Reported", "Sampled"), cex = 0.75)

#keep vector
keep.Bajet.but = c('muPq_butachlor_Bajet12_uncertainty', 'baj.but.m.lin', 'lc50.but.baj', 'lc50.but.baj.se')    

#Full keep vector
keep.Bajet.all <- c(keep.Bajet.lamcy, keep.Bajet.delt, keep.Bajet.cyper, keep.Bajet.chlor,
                    keep.Bajet.prof, keep.Bajet.mal, keep.Bajet.carb, "muPq_24d_Bajet12_uncertainty",
                    keep.Bajet.but)

