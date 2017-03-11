#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

load('Review_models/r0_of_Atrazine_ws.RData')
load('Review_models/r0_of_malathion2_ws2.RData')

  keepers = c('conc.atr', 'conc.mal', 'conc.atr.means.r0', 'conc.atr.sds.r0', 'r0.He', 'r0.In',
              'conc.mal.means.r0', 'conc.mal.sds.r0', 'r0.atr.fix.n0', 'r0.mal.fix.n0', 'parameters',
              'nil1', 'nil0')
  
  rm(list = setdiff(ls(), keepers))

require(rootSolve)
require(ggplot2)
  rec.atr = 500
  rec.mal = 100000

#Make gg compatible df for atrazine results #############
atrm = matrix(ncol = 5, nrow = length(conc.atr)*4)
  colnames(atrm) = c('conc','Modeled','val','UL', 'LL')

  atrm[,1] = rep(conc.atr, 4)
  atrm[,2] = c(rep('Cercarial mortality', length(conc.atr)),
               rep('Snail mortality', length(conc.atr)),
               rep('Snail carrying capacity', length(conc.atr)),
               rep('all', length(conc.atr)))
  atrm[c(1:501),3] = lowess(conc.atr, conc.atr.means.r0[,4] - r0.He(0)[3], f = 0.01)$y
  atrm[c(1:501),4] = lowess(conc.atr, 
                            conc.atr.means.r0[,4] - r0.He(0)[3] + conc.atr.sds.r0[,4], f = 0.01)$y
  atrm[c(1:501),5] = lowess(conc.atr, 
                            conc.atr.means.r0[,4] - r0.He(0)[3] - conc.atr.sds.r0[,4], f = 0.01)$y
  atrm[c(502:1002),3] = lowess(conc.atr, conc.atr.means.r0[,5] - r0.He(0)[3], f = 0.01)$y
  atrm[c(502:1002),4] = lowess(conc.atr, 
                               conc.atr.means.r0[,5] - r0.He(0)[3] + conc.atr.sds.r0[,5], f = 0.01)$y
  atrm[c(502:1002),5] = lowess(conc.atr, 
                               conc.atr.means.r0[,5] - r0.He(0)[3] - conc.atr.sds.r0[,5], f = 0.01)$y
  atrm[c(1003:1503),3] = lowess(conc.atr, conc.atr.means.r0[,6] - r0.He(0)[3], f = 0.01)$y
  atrm[c(1003:1503),4] = lowess(conc.atr, 
                                conc.atr.means.r0[,6] - r0.He(0)[3] + conc.atr.sds.r0[,6], f = 0.01)$y
  atrm[c(1003:1503),5] = lowess(conc.atr, 
                                conc.atr.means.r0[,6] - r0.He(0)[3] - conc.atr.sds.r0[,6], f = 0.01)$y
  atrm[c(1504:2004),3] = lowess(conc.atr, conc.atr.means.r0[,13] - r0.He(0)[3], f = 0.01)$y
  atrm[c(1504:2004),4] = lowess(conc.atr, 
                                conc.atr.means.r0[,13] - r0.He(0)[3] + conc.atr.sds.r0[,13], f = 0.01)$y
  atrm[c(1504:2004),5] = lowess(conc.atr, 
                                conc.atr.means.r0[,13] - r0.He(0)[3] - conc.atr.sds.r0[,13], f = 0.01)$y
  
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
  r0.atr.gg = subset(r0.atr.fix.n0, study != 'Koprivnikar06' & study != 'Griggs08')
    r0.atr.gg$par = as.character(r0.atr.gg$par)
    r0.atr.gg$par[r0.atr.gg$par == 'phi_N'] = 'Snail carrying capacity'
    r0.atr.gg$par[r0.atr.gg$par == 'muN'] = 'Snail mortality'
    r0.atr.gg$par[r0.atr.gg$par == 'piC'] = 'Cercarial mortality'
    
    colnames(r0.atr.gg)[2] = 'Observed'

#Plot atrazine component effects############        
gatr<-ggplot(atr.ggpars) +
      theme_bw() +
      theme(legend.position = 'bottom', legend.direction = 'vertical') +
      xlim(0,500) +
      ylim(-3,3) +
      labs(title = expression(paste('Component effects of Atrazine on ', 'R'[0])),
           x = 'Atrazine (ppb)',
           y = expression(paste(Delta, R[0]))) +
      geom_hline(yintercept = 0, lty=2) + 
      geom_ribbon(aes(x = conc, ymin = LL, ymax = UL, fill = Modeled), alpha = 0.4) +
      geom_line(aes(x = conc, y = val, col = Modeled), size = 1)
      
  #add observed points
    gatr + geom_point(data = r0.atr.gg, aes(x = atr, y = r0 - r0.He(0)[3], shape = Observed), size = 2)

#Make gg compatible df for malathion results #############     
malm = matrix(ncol = 5, nrow = length(conc.mal)*6)
colnames(malm) = c('conc','Modeled','val','UL', 'LL')

malm[,1] = rep(conc.mal, 6)
malm[,2] = c(rep('Cercarial mortality', length(conc.mal)),
             rep('Miracidial mortality', length(conc.mal)),
             rep('Predator mortality', length(conc.mal)),
             rep('Snail mortality', length(conc.mal)),
             rep('Snail reproduction', length(conc.mal)),
             rep('all', length(conc.mal)))

malm[c(1:length(conc.mal)),3] = lowess(conc.mal, conc.mal.means.r0[,1] - r0.In(0)[3], f = 0.01)$y
malm[c(1:length(conc.mal)),4] = lowess(conc.mal, 
                          conc.mal.means.r0[,1] - r0.In(0)[3] + conc.mal.sds.r0[,1], f = 0.01)$y
malm[c(1:length(conc.mal)),5] = lowess(conc.mal, 
                          conc.mal.means.r0[,1] - r0.In(0)[3] - conc.mal.sds.r0[,1], f = 0.01)$y

malm[c(1993:(length(conc.mal)*2)),3] = lowess(conc.mal, conc.mal.means.r0[,2] - r0.In(0)[3], f = 0.01)$y
malm[c(1993:(length(conc.mal)*2)),4] = lowess(conc.mal, 
                             conc.mal.means.r0[,2] - r0.In(0)[3] + conc.mal.sds.r0[,2], f = 0.01)$y
malm[c(1993:(length(conc.mal)*2)),5] = lowess(conc.mal, 
                             conc.mal.means.r0[,2] - r0.In(0)[3] - conc.mal.sds.r0[,2], f = 0.01)$y

malm[c(3985:(length(conc.mal)*3)),3] = lowess(conc.mal, conc.mal.means.r0[,5] - r0.In(0)[3], f = 0.01)$y
malm[c(3985:(length(conc.mal)*3)),4] = lowess(conc.mal, 
                             conc.mal.means.r0[,5] - r0.In(0)[3] + conc.mal.sds.r0[,5], f = 0.01)$y
malm[c(3985:(length(conc.mal)*3)),5] = lowess(conc.mal, 
                             conc.mal.means.r0[,5] - r0.In(0)[3] - conc.mal.sds.r0[,5], f = 0.01)$y

malm[c(5977:(length(conc.mal)*4)),3] = lowess(conc.mal, conc.mal.means.r0[,6] - r0.In(0)[3], f = 0.01)$y
malm[c(5977:(length(conc.mal)*4)),4] = lowess(conc.mal, 
                              conc.mal.means.r0[,6] - r0.In(0)[3] + conc.mal.sds.r0[,6], f = 0.01)$y
malm[c(5977:(length(conc.mal)*4)),5] = lowess(conc.mal, 
                              conc.mal.means.r0[,6] - r0.In(0)[3] - conc.mal.sds.r0[,6], f = 0.01)$y

malm[c(7969:(length(conc.mal)*5)),3] = lowess(conc.mal, conc.mal.means.r0[,7] - r0.In(0)[3], f = 0.01)$y
malm[c(7969:(length(conc.mal)*5)),4] = lowess(conc.mal, 
                              conc.mal.means.r0[,7] - r0.In(0)[3] + conc.mal.sds.r0[,7], f = 0.01)$y
malm[c(7969:(length(conc.mal)*5)),5] = lowess(conc.mal, 
                              conc.mal.means.r0[,7] - r0.In(0)[3] - conc.mal.sds.r0[,7], f = 0.01)$y

malm[c(9961:(length(conc.mal)*6)),3] = lowess(conc.mal, conc.mal.means.r0[,13] - r0.In(0)[3], f = 0.01)$y
malm[c(9961:(length(conc.mal)*6)),4] = lowess(conc.mal, 
                              conc.mal.means.r0[,13] - r0.In(0)[3] + conc.mal.sds.r0[,13], f = 0.01)$y
malm[c(9961:(length(conc.mal)*6)),5] = lowess(conc.mal, 
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
    theme(legend.position = 'bottom', legend.direction = 'vertical') +
    xlim(0,150000) +
    ylim(-3,3.1) +
    labs(title = expression(paste('Component effects of Malathion on ', 'R'[0])),
         x = 'Malathion (ppb)',
         y = expression(paste(Delta, R[0]))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc, ymin = LL, ymax = UL, fill = Modeled), alpha = 0.4) +
    geom_line(aes(x = conc, y = val, col = Modeled), size = 1)
  
  #add observed points
  gmal + geom_point(data = r0.mal.gg, aes(x = mal, y = r0 - r0.In(0)[3], shape = Observed), size = 2)
  
  
#Plot atrazine and malathion component effects as facets ################
  comp = rbind(atr.ggpars, mal.ggpars)
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
    
#Plot atrazine and malathion combined effects together ###########
  comb = rbind(atr.ggall, mal.ggall)
  gall<-ggplot(comb) +
    theme_bw() +
    theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = c(0,1), labels = c('0', 'REC'), limits = c(0,1)) +
    ylim(-4,4) +
    labs(title = expression(paste('Combined agrochemical effects on ', 'R'[0])),
         x = 'Normalized Agrochemical Concentration',
         y = expression(paste(Delta, R[0]))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = rec, ymin = LL, ymax = UL, fill = Chemical), alpha = 0.4) +
    geom_line(aes(x = rec, y = val, col = Chemical), size = 1)