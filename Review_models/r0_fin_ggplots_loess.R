#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

load('Review_models/Savio/Atrazine/r0_atrazine_savio.RData')
load('Review_models/Savio/Malathion/r0_malathion_savio.RData')
load('Review_models/Savio/Chlorpyrifos/r0_chlorpyrifos_savio.RData')

  keepers = c('conc.atr','conc.atr.means.r0', 'conc.atr.sds.r0','r0.atr.fix.n0',  
              'conc.ch', 'conc.ch.means.r0', 'conc.ch.sds.r0',
              'conc.mal', 'conc.mal.means.r0', 'conc.mal.sds.r0','r0.mal.fix.n0', 
              'parameters','r0.He', 'r0.In','nil1', 'nil0')
  
  rm(list = setdiff(ls(), keepers))

require(rootSolve)
require(ggplot2)
  rec.atr = 102*2   #From Halstead et al Nature Comm 2017 table S6
  rec.mal = 28*2    #From Halstead et al Chemosphere 2015 table S1
  rec.ch = 64*2   #From Halstead et al Nature Comm 2017 table S6

#Make gg compatible df for atrazine results #############
atrm = matrix(ncol = 5, nrow = length(conc.atr)*4)
  colnames(atrm) = c('conc','Modeled','val','UL', 'LL')

  atrm[,1] = rep(conc.atr, 4)
  atrm[,2] = c(rep('Cercarial mortality', length(conc.atr)),     #From Rohr 2008 study
               rep('Snail mortality', length(conc.atr)),         #from Omran & Salama study
               rep('Snail carrying capacity', length(conc.atr)), #from Rohr reanalysis of Baxter study
               rep('all', length(conc.atr)))                     #combined of the three above
#Fill cercarial mortality From Rohr 2008 study
  atrm[c(1:length(conc.atr)),3] = lowess(conc.atr, conc.atr.means.r0[,4] - r0.He(0)[3], f = 0.01)$y
  atrm[c(1:length(conc.atr)),4] = lowess(conc.atr, 
                            conc.atr.means.r0[,4] - r0.He(0)[3] + conc.atr.sds.r0[,4], f = 0.01)$y
  atrm[c(1:length(conc.atr)),5] = lowess(conc.atr, 
                            conc.atr.means.r0[,4] - r0.He(0)[3] - conc.atr.sds.r0[,4], f = 0.01)$y
#Fill snail mortality from Omran & Salama study
  atrm[c((length(conc.atr)+1):(2*length(conc.atr))),3] = lowess(conc.atr, conc.atr.means.r0[,7] - r0.He(0)[3], f = 0.01)$y
  atrm[c((length(conc.atr)+1):(2*length(conc.atr))),4] = lowess(conc.atr, 
                               conc.atr.means.r0[,7] - r0.He(0)[3] + conc.atr.sds.r0[,7], f = 0.01)$y
  atrm[c((length(conc.atr)+1):(2*length(conc.atr))),5] = lowess(conc.atr, 
                               conc.atr.means.r0[,7] - r0.He(0)[3] - conc.atr.sds.r0[,7], f = 0.01)$y
#Fill snail carrying capacity from Rohr reanalysis of Baxter study  
  atrm[c((2*length(conc.atr)+1):(3*length(conc.atr))),3] = lowess(conc.atr, conc.atr.means.r0[,6] - r0.He(0)[3], f = 0.01)$y
  atrm[c((2*length(conc.atr)+1):(3*length(conc.atr))),4] = lowess(conc.atr, 
                                conc.atr.means.r0[,6] - r0.He(0)[3] + conc.atr.sds.r0[,6], f = 0.01)$y
  atrm[c((2*length(conc.atr)+1):(3*length(conc.atr))),5] = lowess(conc.atr, 
                                conc.atr.means.r0[,6] - r0.He(0)[3] - conc.atr.sds.r0[,6], f = 0.01)$y
#Fill combo of three  
  atrm[c((3*length(conc.atr)+1):(4*length(conc.atr))),3] = lowess(conc.atr, conc.atr.means.r0[,17] - r0.He(0)[3], f = 0.01)$y
  atrm[c((3*length(conc.atr)+1):(4*length(conc.atr))),4] = lowess(conc.atr, 
                                conc.atr.means.r0[,17] - r0.He(0)[3] + conc.atr.sds.r0[,17], f = 0.01)$y
  atrm[c((3*length(conc.atr)+1):(4*length(conc.atr))),5] = lowess(conc.atr, 
                                conc.atr.means.r0[,17] - r0.He(0)[3] - conc.atr.sds.r0[,17], f = 0.01)$y
  
  atr.gg<-as.data.frame(atrm, stringsAsFactors = FALSE)
    atr.gg$conc<-as.numeric(atr.gg$conc)
    atr.gg$val<-as.numeric(atr.gg$val)
    atr.gg$UL<-as.numeric(atr.gg$UL)
    atr.gg$LL<-as.numeric(atr.gg$LL)
    
 atr.ggpars = subset(atr.gg, Modeled != 'all')  
  atr.ggpars$Chemical = 'Atrazine'
 atr.ggall = subset(atr.gg, Modeled == 'all')
  atr.ggall$Chemical = 'Atrazine'
  atr.ggall$rec = atr.ggall$conc / rec.atr
 
#Edit atrazine observed points df for ggplot ##############      
  r0.atr.gg = subset(r0.atr.fix.n0, study != 'Koprivnikar06' & study != 'Griggs08' & study != 'Bakry12')
    r0.atr.gg$par = as.character(r0.atr.gg$par)
    r0.atr.gg$par[r0.atr.gg$par == 'phi_N'] = 'Snail carrying capacity'
    r0.atr.gg$par[r0.atr.gg$par == 'muN'] = 'Snail mortality'
    r0.atr.gg$par[r0.atr.gg$par == 'piC'] = 'Cercarial mortality'
    
    colnames(r0.atr.gg)[2] = 'Observed'

#Plot atrazine component effects############        
gatr<-ggplot(atr.ggpars) +
      theme_bw() +
      theme(legend.position = 'bottom', legend.direction = 'vertical') +
      xlim(0,2000) +
      ylim(-4,4) +
      labs(title = expression(paste('Component effects of Atrazine on ', 'R'[0])),
           x = 'Atrazine (ppb)',
           y = expression(paste(Delta, R[0]))) +
      geom_hline(yintercept = 0, lty=2) + 
      geom_ribbon(aes(x = conc, ymin = LL, ymax = UL, fill = Modeled), alpha = 0.4) +
      geom_line(aes(x = conc, y = val, col = Modeled), size = 1)
gatr      
  #add observed points
    gatr + geom_point(data = r0.atr.gg, aes(x = atr, y = r0 - r0.He(0)[3], shape = Observed), size = 2)

#Make gg compatible df for malathion results #############     
malm = matrix(ncol = 5, nrow = length(conc.mal)*6)
colnames(malm) = c('conc','Modeled','val','UL', 'LL')

malm[,1] = rep(conc.mal, 6)
malm[,2] = c(rep('Cercarial mortality', length(conc.mal)),   #from Tchounwou92
             rep('Miracidial mortality', length(conc.mal)),  #from Tchounwou91
             rep('Predator mortality', length(conc.mal)),    #from Halstead 2015
             rep('Snail mortality', length(conc.mal)),       #From tchounwou 91(uses bulinus species)
             rep('Snail reproduction', length(conc.mal)),    #From tchounwou 91(uses bulinus species)
             rep('all', length(conc.mal)))                   #all of the above combined
#fill cercarial mortality from Tchounwou92
  malm[c(1:length(conc.mal)),3] = lowess(conc.mal, conc.mal.means.r0[,1] - r0.In(0)[3], f = 0.01)$y
  malm[c(1:length(conc.mal)),4] = lowess(conc.mal, 
                            conc.mal.means.r0[,1] - r0.In(0)[3] + conc.mal.sds.r0[,1], f = 0.01)$y
  malm[c(1:length(conc.mal)),5] = lowess(conc.mal, 
                            conc.mal.means.r0[,1] - r0.In(0)[3] - conc.mal.sds.r0[,1], f = 0.01)$y
#fill miracidial mortality from Tchounwou91
  malm[c((length(conc.mal)+1):(length(conc.mal)*2)),3] = lowess(conc.mal, conc.mal.means.r0[,2] - r0.In(0)[3], f = 0.01)$y
  malm[c((length(conc.mal)+1):(length(conc.mal)*2)),4] = lowess(conc.mal, 
                               conc.mal.means.r0[,2] - r0.In(0)[3] + conc.mal.sds.r0[,2], f = 0.01)$y
  malm[c((length(conc.mal)+1):(length(conc.mal)*2)),5] = lowess(conc.mal, 
                               conc.mal.means.r0[,2] - r0.In(0)[3] - conc.mal.sds.r0[,2], f = 0.01)$y
#fill predator mortality from Halstead15  
  malm[c((2*length(conc.mal)+1):(length(conc.mal)*3)),3] = lowess(conc.mal, conc.mal.means.r0[,5] - r0.In(0)[3], f = 0.01)$y
  malm[c((2*length(conc.mal)+1):(length(conc.mal)*3)),4] = lowess(conc.mal, 
                               conc.mal.means.r0[,5] - r0.In(0)[3] + conc.mal.sds.r0[,5], f = 0.01)$y
  malm[c((2*length(conc.mal)+1):(length(conc.mal)*3)),5] = lowess(conc.mal, 
                               conc.mal.means.r0[,5] - r0.In(0)[3] - conc.mal.sds.r0[,5], f = 0.01)$y
#fill snail mortality from Tchounwou91  
  malm[c((3*length(conc.mal)+1):(length(conc.mal)*4)),3] = lowess(conc.mal, conc.mal.means.r0[,6] - r0.In(0)[3], f = 0.01)$y
  malm[c((3*length(conc.mal)+1):(length(conc.mal)*4)),4] = lowess(conc.mal, 
                                conc.mal.means.r0[,6] - r0.In(0)[3] + conc.mal.sds.r0[,6], f = 0.01)$y
  malm[c((3*length(conc.mal)+1):(length(conc.mal)*4)),5] = lowess(conc.mal, 
                                conc.mal.means.r0[,6] - r0.In(0)[3] - conc.mal.sds.r0[,6], f = 0.01)$y
#fill snail reproduction reduction from Tchounwou92  
  malm[c((4*length(conc.mal)+1):(length(conc.mal)*5)),3] = lowess(conc.mal, conc.mal.means.r0[,7] - r0.In(0)[3], f = 0.01)$y
  malm[c((4*length(conc.mal)+1):(length(conc.mal)*5)),4] = lowess(conc.mal, 
                                conc.mal.means.r0[,7] - r0.In(0)[3] + conc.mal.sds.r0[,7], f = 0.01)$y
  malm[c((4*length(conc.mal)+1):(length(conc.mal)*5)),5] = lowess(conc.mal, 
                                conc.mal.means.r0[,7] - r0.In(0)[3] - conc.mal.sds.r0[,7], f = 0.01)$y
#fill combo of all  
  malm[c((5*length(conc.mal)+1):(length(conc.mal)*6)),3] = lowess(conc.mal, conc.mal.means.r0[,13] - r0.In(0)[3], f = 0.01)$y
  malm[c((5*length(conc.mal)+1):(length(conc.mal)*6)),4] = lowess(conc.mal, 
                                conc.mal.means.r0[,13] - r0.In(0)[3] + conc.mal.sds.r0[,13], f = 0.01)$y
  malm[c((5*length(conc.mal)+1):(length(conc.mal)*6)),5] = lowess(conc.mal, 
                                conc.mal.means.r0[,13] - r0.In(0)[3] - conc.mal.sds.r0[,13], f = 0.01)$y

mal.gg<-as.data.frame(malm, stringsAsFactors = FALSE)
  mal.gg$conc<-as.numeric(mal.gg$conc)
  mal.gg$val<-as.numeric(mal.gg$val)
  mal.gg$UL<-as.numeric(mal.gg$UL)
  mal.gg$LL<-as.numeric(mal.gg$LL)

  mal.ggpars = subset(mal.gg, Modeled != 'all' & conc <=150000) 
    mal.ggpars$Chemical = 'Malathion'
  mal.ggall = subset(mal.gg, Modeled == 'all' & conc <=150000)
    mal.ggall$Chemical = 'Malathion'
    mal.ggall$rec = mal.ggall$conc / rec.mal

#Edit malathion observed points df for ggplot ##############      
r0.mal.gg = subset(r0.mal.fix.n0, study != 'Bakry11')
  r0.mal.gg$par = as.character(r0.mal.gg$par)
  r0.mal.gg$par[r0.mal.gg$par == 'piM'] = 'Miracidial mortality'
  r0.mal.gg$par[r0.mal.gg$par == 'muN'] = 'Snail mortality'
  r0.mal.gg$par[r0.mal.gg$par == 'muP'] = 'Predator mortality'
  r0.mal.gg$par[r0.mal.gg$par == 'piC'] = 'Cercarial mortality'
  r0.mal.gg$par[r0.mal.gg$par == 'fN'] = 'Snail reproduction'
  
  colnames(r0.mal.gg)[2] = 'Observed'
  
#Plot malathion component effects############        
  gmal<-ggplot(mal.ggpars) +
    theme_bw() +
    theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank()) +
    xlim(0,150000) +
    ylim(-4,4) +
    labs(title = expression(paste('Component effects of Malathion on ', 'R'[0])),
         x = 'Malathion (ppb)',
         y = expression(paste(Delta, R[0]))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc, ymin = LL, ymax = UL, fill = Modeled), alpha = 0.4) +
    geom_line(aes(x = conc, y = val, col = Modeled), size = 1)

gmal    
  #add observed points
  gmal + geom_point(data = r0.mal.gg, aes(x = mal, y = r0 - r0.In(0)[3], shape = Observed), size = 2)
  
  
#Make gg compatible df for chlorpyrifos results #############     
chm = matrix(ncol = 5, nrow = length(conc.ch)*8)
colnames(chm) = c('conc','Modeled','val','UL', 'LL')

chm[,1] = rep(conc.ch, 8)
chm[,2] = c(rep('Predator mortality (Halstead)', length(conc.ch)),   #from Halstead 2015
            rep('Predator mortality (Satapornvanit)', length(conc.ch)),   #from Satapornvanit
            rep('Predator attack rate', length(conc.ch)), #from Satapornvanit 
            rep('Cercarial mortality', length(conc.ch)),  #from Hasheesh 2011
            rep('Miracidial mortality', length(conc.ch)), #from Hasheesh 2011
            rep('Snail mortality', length(conc.ch)),      #from Hasheesh 2011
            rep('Snail reproduction', length(conc.ch)),   #From Ibrahim 1992
            rep('all', length(conc.ch)))                  #all of the above combined

#fill predator mortality from Halstead 2015 *******************************************************************************
  chm[c(1:length(conc.ch)),3] = lowess(conc.ch, conc.ch.means.r0[,1] - r0.In(0)[3], f = 0.01)$y
  chm[c(1:length(conc.ch)),4] = lowess(conc.ch, 
                                         conc.ch.means.r0[,1] - r0.In(0)[3] + conc.ch.sds.r0[,1], f = 0.01)$y
  chm[c(1:length(conc.ch)),5] = lowess(conc.ch, 
                                         conc.ch.means.r0[,1] - r0.In(0)[3] - conc.ch.sds.r0[,1], f = 0.01)$y
#fill predator mortality from Satapornvanit ******************************************************************************
  chm[c((length(conc.ch)+1):(length(conc.ch)*2)),3] = lowess(conc.ch, conc.ch.means.r0[,2] - r0.In(0)[3], f = 0.01)$y
  chm[c((length(conc.ch)+1):(length(conc.ch)*2)),4] = lowess(conc.ch, 
                                                                conc.ch.means.r0[,2] - r0.In(0)[3] + conc.ch.sds.r0[,2], f = 0.01)$y
  chm[c((length(conc.ch)+1):(length(conc.ch)*2)),5] = lowess(conc.ch, 
                                                                conc.ch.means.r0[,2] - r0.In(0)[3] - conc.ch.sds.r0[,2], f = 0.01)$y
#fill predator attack rate from Satapornvanit ******************************************************************************
  chm[c((2*length(conc.ch)+1):(length(conc.ch)*3)),3] = lowess(conc.ch, conc.ch.means.r0[,3] - r0.In(0)[3], f = 0.01)$y
  chm[c((2*length(conc.ch)+1):(length(conc.ch)*3)),4] = lowess(conc.ch, 
                                                             conc.ch.means.r0[,3] - r0.In(0)[3] + conc.ch.sds.r0[,3], f = 0.01)$y
  chm[c((2*length(conc.ch)+1):(length(conc.ch)*3)),5] = lowess(conc.ch, 
                                                             conc.ch.means.r0[,3] - r0.In(0)[3] - conc.ch.sds.r0[,3], f = 0.01)$y
#fill cercarial mortality from Hasheesh 2011 ******************************************************************************
  chm[c((3*length(conc.ch)+1):(length(conc.ch)*4)),3] = lowess(conc.ch, conc.ch.means.r0[,4] - r0.In(0)[3], f = 0.01)$y
  chm[c((3*length(conc.ch)+1):(length(conc.ch)*4)),4] = lowess(conc.ch, 
                                                               conc.ch.means.r0[,4] - r0.In(0)[3] + conc.ch.sds.r0[,4], f = 0.01)$y
  chm[c((3*length(conc.ch)+1):(length(conc.ch)*4)),5] = lowess(conc.ch, 
                                                               conc.ch.means.r0[,4] - r0.In(0)[3] - conc.ch.sds.r0[,4], f = 0.01)$y
#fill miracidial mortality from Hasheesh 2011 ******************************************************************************
  chm[c((4*length(conc.ch)+1):(length(conc.ch)*5)),3] = lowess(conc.ch, conc.ch.means.r0[,5] - r0.In(0)[3], f = 0.01)$y
  chm[c((4*length(conc.ch)+1):(length(conc.ch)*5)),4] = lowess(conc.ch, 
                                                                  conc.ch.means.r0[,5] - r0.In(0)[3] + conc.ch.sds.r0[,5], f = 0.01)$y
  chm[c((4*length(conc.ch)+1):(length(conc.ch)*5)),5] = lowess(conc.ch, 
                                                                  conc.ch.means.r0[,5] - r0.In(0)[3] - conc.ch.sds.r0[,5], f = 0.01)$y
#fill snail mortality from Hasheesh 2011  ******************************************************************************
  chm[c((5*length(conc.ch)+1):(length(conc.ch)*6)),3] = lowess(conc.ch, conc.ch.means.r0[,6] - r0.In(0)[3], f = 0.01)$y
  chm[c((5*length(conc.ch)+1):(length(conc.ch)*6)),4] = lowess(conc.ch, 
                                                                  conc.ch.means.r0[,6] - r0.In(0)[3] + conc.ch.sds.r0[,6], f = 0.01)$y
  chm[c((5*length(conc.ch)+1):(length(conc.ch)*6)),5] = lowess(conc.ch, 
                                                                  conc.ch.means.r0[,6] - r0.In(0)[3] - conc.ch.sds.r0[,6], f = 0.01)$y
#fill snail reproduction reduction from Ibrahim 1992  ******************************************************************************
  chm[c((6*length(conc.ch)+1):(length(conc.ch)*7)),3] = lowess(conc.ch, conc.ch.means.r0[,9] - r0.In(0)[3], f = 0.01)$y
  chm[c((6*length(conc.ch)+1):(length(conc.ch)*7)),4] = lowess(conc.ch, 
                                                                  conc.ch.means.r0[,9] - r0.In(0)[3] + conc.ch.sds.r0[,9], f = 0.01)$y
  chm[c((6*length(conc.ch)+1):(length(conc.ch)*7)),5] = lowess(conc.ch, 
                                                                  conc.ch.means.r0[,9] - r0.In(0)[3] - conc.ch.sds.r0[,9], f = 0.01)$y
#fill combo of all  ******************************************************************************  
  chm[c((7*length(conc.ch)+1):(length(conc.ch)*8)),3] = lowess(conc.ch, conc.ch.means.r0[,19] - r0.In(0)[3], f = 0.01)$y
  chm[c((7*length(conc.ch)+1):(length(conc.ch)*8)),4] = lowess(conc.ch, 
                                                                  conc.ch.means.r0[,19] - r0.In(0)[3] + conc.ch.sds.r0[,19], f = 0.01)$y
  chm[c((7*length(conc.ch)+1):(length(conc.ch)*8)),5] = lowess(conc.ch, 
                                                                  conc.ch.means.r0[,19] - r0.In(0)[3] - conc.ch.sds.r0[,19], f = 0.01)$y

ch.gg<-as.data.frame(chm, stringsAsFactors = FALSE)
  ch.gg$conc<-as.numeric(ch.gg$conc)
  ch.gg$val<-as.numeric(ch.gg$val)
  ch.gg$UL<-as.numeric(ch.gg$UL)
  ch.gg$LL<-as.numeric(ch.gg$LL)
  
ch.ggpars = subset(ch.gg, Modeled != 'all' & conc <=150000) 
  ch.ggpars$Chemical = 'Chlorpyrifos'
ch.ggall = subset(ch.gg, Modeled == 'all' & conc <=150000)
  ch.ggall$Chemical = 'Chlorpyrifos'
  ch.ggall$rec = ch.ggall$conc / rec.ch

#Edit chlorpyrifos observed points df for ggplot ##############   
  #Data frame doesn't exist yet, will generate if need be
#r0.ch.gg = subset(r0.ch.fix.n0, study != 'Bakry11')
#r0.ch.gg$par = as.character(r0.ch.gg$par)
#r0.ch.gg$par[r0.ch.gg$par == 'piM'] = 'Miracidial mortality'
#r0.ch.gg$par[r0.ch.gg$par == 'muN'] = 'Snail mortality'
#r0.ch.gg$par[r0.ch.gg$par == 'muP'] = 'Predator mortality'
#r0.ch.gg$par[r0.ch.gg$par == 'piC'] = 'Cercarial mortality'
#r0.ch.gg$par[r0.ch.gg$par == 'fN'] = 'Snail reproduction'

#colnames(r0.ch.gg)[2] = 'Observed'

#Plot chlorpyrifos component effects############        
gch<-ggplot(ch.ggpars) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank()) +
  xlim(0,4000) +
  ylim(-5,5) +
  labs(title = expression(paste('Component effects of Chlorpyrifos on ', 'R'[0])),
       x = 'Chlorpyrifos (ppb)',
       y = expression(paste(Delta, R[0]))) +
  geom_hline(yintercept = 0, lty=2) + 
  geom_ribbon(aes(x = conc, ymin = LL, ymax = UL, fill = Modeled), alpha = 0.4) +
  geom_line(aes(x = conc, y = val, col = Modeled), size = 1)

gch    
#add observed points
#gch + geom_point(data = r0.ch.gg, aes(x = ch, y = r0 - r0.In(0)[3], shape = Observed), size = 2)


#Plot atrazine, malathion, and chlorpyrifos component effects as facets ################
  comp = rbind(atr.ggpars, mal.ggpars, ch.ggpars)
  gcomp = ggplot(data = comp) +
    facet_grid(. ~ Chemical, scales = 'free_x') +
    theme_bw() +
    theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ylim(-3.3,3.3) +
    labs(x = 'Agrochemical conc (ppb)',
         y = expression(paste(Delta, R[0]))) +
    ggtitle(label = expression(paste('Component agrochemical effects on ', 'R'[0]))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc, ymin = LL, ymax = UL, fill = Modeled), alpha = 0.4) +
    geom_line(aes(x = conc, y = val, col = Modeled), size = 1)
  
  gcomp
    
#Plot atrazine and malathion combined effects together ###########
  comb = rbind(atr.ggall, mal.ggall, ch.ggall)
  gall<-ggplot(comb) +
    theme_bw() +
    theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = c(0,1,2), labels = c('0', 'REC','2REC'), limits = c(0,2.1)) +
    ylim(-4,4) +
    labs(title = expression(paste('Combined agrochemical effects on ', 'R'[0])),
         x = 'Normalized Agrochemical Concentration',
         y = expression(paste(Delta, R[0]))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = rec, ymin = LL, ymax = UL, fill = Chemical), alpha = 0.4) +
    geom_line(aes(x = rec, y = val, col = Chemical), size = 1)
  
  gall