#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#herbicide toxicity to Bi. alexandrina from Omran and Salama 2013 ###############

#Load packages, clean/supplement data, fit initial model #############
require(drc)
require(LW1949)

ons = read.csv('~/RemaisWork/Schisto/Data/AgroData/Data/Snail Mortality/Omran&Salama2013.csv')
  ons$mort = round(pnorm(ons$probit, mean = 5), digits = 1)
  ons$total = 10
  ons$dead = ons$total * ons$mort

#atrazine mortality ##########      
  atr.dat = dataprep(dose = ons$ppm[ons$chem == 'atrazine'],
                     ntot = ons$total[ons$chem == 'atrazine'],
                     nfx = ons$dead[ons$chem == 'atrazine'])
  
  atr.lw = fitLWauto(atr.dat)
  fatr = LWestimate(atr.lw, atr.dat)
  
    lc50.ons.atr = 40.1425 #fatr$LWest[1]
    slp.ons.atr = 2.637754 #fatr$LWest[12]
    se.lc50.ons.atr = 0.1459167 #log10(fatr$LWest[3]/lc50.ons.atr) / qnorm(0.975)
    
  ons.munq.atr = function(He){
    #if(In == 0) mun = 0 else{
    heu = He/1000
    lc50 = 10^(rnorm(1, log10(lc50.ons.atr), se.lc50.ons.atr))
    mun = pnorm((slp.ons.atr) * log10(heu/lc50)) 
    #}
      
    return(mun)
  }
  
  plot(ons.atr$conc, ons.atr$mort, pch = 16, ylim = c(0,1), #xlim = c(0,1000),
       xlab = 'atrazine (ppb)', ylab = c('mortality')) 
    points(seq(0,500000,500), sapply(seq(0,500000,500), ons.munq.atr),
           pch = 5, col = 4, cex = 0.5)
    
  plotDE(atr.dat)
  predLines(fatr)
  points(seq(0,500,1), sapply(seq(0,500000,1000), ons.munq.atr)*100,
         pch = 5, col = 4, cex = 0.5)

keep.ons.atr = c('ons.munq.atr', 'lc50.ons.atr', 'slp.ons.atr', 'se.lc50.ons.atr')

#Glyphosate mortality #########
  gly.dat = dataprep(dose = ons$ppm[ons$chem == 'glyphosate'],
                     ntot = ons$total[ons$chem == 'glyphosate'],
                     nfx = ons$dead[ons$chem == 'glyphosate'])
  
  gly.lw = fitLWauto(gly.dat)
  fgly = LWestimate(gly.lw, gly.dat)
  
  lc50.ons.gly = 118.709 #fgly$LWest[1]
  slp.ons.gly = 3.847287 # fgly$LWest[12]
  se.lc50.ons.gly = 0.1564323 #log10(fgly$LWest[3]/lc50.ons.gly) / qnorm(0.975)
  
  ons.munq.gly = function(He){
    #if(In == 0) mun = 0 else{
    heu = He/1000
    lc50 = 10^(rnorm(1, log10(lc50.ons.gly), se.lc50.ons.gly))
    mun = pnorm((slp.ons.gly) * log10(heu/lc50)) 
    #}
    
    return(mun)
  }
  
  plot(ons.gly$conc, ons.gly$mort, pch = 16, ylim = c(0,1), #xlim = c(0,1000),
       xlab = 'glyphosate (ppb)', ylab = c('mortality')) 
    points(seq(0,500000,500), sapply(seq(0,500000,500), ons.munq.gly),
           pch = 5, col = 4, cex = 0.5)
  
  plotDE(gly.dat)
  predLines(fgly)
  points(seq(0,500,1), sapply(seq(0,500000,1000), ons.munq.gly)*100,
         pch = 5, col = 4, cex = 0.5)
  
  keep.ons.gly = c('ons.munq.gly', 'lc50.ons.gly', 'slp.ons.gly', 'se.lc50.ons.gly')

keep.ons.all = c(keep.ons.atr, keep.ons.gly)   