#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

#Plot individual chem trajectories first ########
require(ggplot2)
#Atrazine ################
load('Review_models/Savio/Atrazine/atr5000sims_r0med_2017-12-12.Rdata')
load('Review_models/Savio/Atrazine/atr5000sims_r0.75_2017-12-12.Rdata')
load('Review_models/Savio/Atrazine/atr5000sims_r0.25_2017-12-12.Rdata')
load('Review_models/Savio/Atrazine/atr5000sims_pars_2017-12-12.Rdata')
 
#Net effect 
  atr.df = data.frame(conc = seq(0,204, length.out = 500),
                      conc.eec = seq(0,204, length.out = 500)/102,
                      chem = 'Atrazine',
                      botup = 'Present',
                      med = atr.med.r0[,7],
                      iqr.25 = atr.med.r0.25[,7],
                      iqr.75 = atr.med.r0.75[,7])
  
#summary of net effect at EEC
  atr.df[which(round(atr.df$conc.eec, digits = 2) == 1.00),][2,]
  
  gatr<-ggplot(atr.df) +
          theme_bw() +
          #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
          xlim(0,204) +
          ylim(-50,100) +
          labs(title = expression(paste('Net effect of Atrazine on ', 'R'[0])),
               x = 'Atrazine (ppb)',
               y = expression(paste(Delta, R[0], ' (%)'))) +
          geom_hline(yintercept = 0, lty=2) + 
          geom_ribbon(aes(x = conc, ymin = iqr.25, ymax = iqr.75), alpha = 0.4) +
          geom_line(aes(x = conc, y = med), size = 1)
  gatr     

#Individual parameter effects    
  par.atr.df = data.frame(conc = seq(0,204, length.out = 500),
                          conc.eec = seq(0,204, length.out = 500)/102,
                          chem = 'Atrazine',
                          botup = 'Present')
  
  par.atr.df = cbind(par.atr.df, atr.med.r0[,c(1:6)], atr.med.r0.25[,c(1:6)], atr.med.r0.75[,c(1:6)])
  
    colnames(par.atr.df)[c(5:22)] = c(paste0('r0.',c('piCgrg', 'piCkop', 'piCror', 'munbak', 'phinbax', 'munons') ),
                                      paste0('r025.' ,c('piCgrg', 'piCkop', 'piCror', 'munbak', 'phinbax', 'munons')),
                                      paste0('r075.' ,c('piCgrg', 'piCkop', 'piCror', 'munbak', 'phinbax', 'munons')))
    
  par.atr.df = reshape(data = par.atr.df, idvar = c('conc', 'conc.eec', 'chem'), drop = 'botup', varying = 5:22, direction = 'long',
                       timevar = 'study', sep = '.', new.row.names = c(1:3000))  
  par.atr.df$keep = 0
  par.atr.df$keep[which(par.atr.df$study %in% c('munons', 'piCror', 'phinbax'))] = 1
    
  gatr.par<-ggplot(par.atr.df) +
              theme_bw() +
              #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
              xlim(0,204) +
              ylim(-100,100) +
              labs(title = expression(paste('Parameter effect of Atrazine on ', 'R'[0])),
                   x = 'Atrazine (ppb)',
                   y = expression(paste(Delta, R[0], ' (%)'))) +
              geom_hline(yintercept = 0, lty=2) + 
              geom_ribbon(aes(x = conc, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
              geom_line(aes(x = conc, y = r0, col = study), size = 1)
    
  gatr.par
  
#restrict to studies used to produce net effect  
  par.atr.df$keep = 0
  par.atr.df$keep[which(par.atr.df$study %in% c('munons', 'piCror', 'phinbax'))] = 1
  par.atr.keeps = subset(par.atr.df, keep == 1)
  
  gatr.par.keeps<-ggplot(par.atr.keeps) +
    theme_bw() +
    theme(legend.position = 'top', legend.direction = 'horizontal') +
    xlim(0,2) +
    ylim(-100,100) +
    labs(x = 'Atrazine (EEC)',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc.eec, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc.eec, y = r0, col = study), size = 1)
  
  gatr.par.keeps
  
#Butralin ################
load('Review_models/Savio/Butralin/but5000sims_r0med_2017-12-12.Rdata')
load('Review_models/Savio/Butralin/but5000sims_r0.75_2017-12-12.Rdata')
load('Review_models/Savio/Butralin/but5000sims_r0.25_2017-12-12.Rdata')
load('Review_models/Savio/Butralin/but5000sims_pars_2017-12-12.Rdata')

  but.df = data.frame(conc = rep(seq(0, 33.8, length.out = 500), 2),
                      conc.eec = rep(seq(0, 33.8, length.out = 500), 2)/16.9,
                      chem = 'Butralin',
                      botup = c(rep('Absent', 500), rep('Present', 500)),
                      med = c(but.med.r0[,8], but.med.r0[,9]),
                      iqr.25 = c(but.med.r0.25[,8], but.med.r0.25[,9]),
                      iqr.75 = c(but.med.r0.75[,8], but.med.r0.75[,9]))
  
#summary of net effect at EEC
  but.df[which(round(but.df$conc.eec, digits = 2) == 1.00),][4,]
  
  
  gbut<-ggplot(but.df) +
          theme_bw() +
          #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
          xlim(0,34) +
          ylim(-100,100) +
          labs(title = expression(paste('Net effect of Butralin on ', 'R'[0])),
               x = 'Butralin (ppb)',
               y = expression(paste(Delta, R[0], ' (%)'))) +
          geom_hline(yintercept = 0, lty=2) + 
          geom_ribbon(aes(x = conc, ymin = iqr.25, ymax = iqr.75, fill = botup), alpha = 0.4) +
          geom_line(aes(x = conc, y = med, col = botup), size = 1)
  gbut   
  
#Individual parameter effects    
  par.but.df = data.frame(conc = seq(0, 33.8, length.out = 500),
                          conc.eec = seq(0, 33.8, length.out = 500)/16.9,
                          chem = 'Butralin')
  
  par.but.df = cbind(par.but.df, but.med.r0[,c(1:5)], but.med.r0.25[,c(1:5)], but.med.r0.75[,c(1:5)])
  
  colnames(par.but.df)[c(4:18)] = c(paste0('r0.',c('piMgaf', 'piCgaf', 'mungaf', 'fnbak', 'phinbax')),
                                    paste0('r025.' ,c('piMgaf', 'piCgaf', 'mungaf', 'fnbak', 'phinbax')),
                                    paste0('r075.' ,c('piMgaf', 'piCgaf', 'mungaf', 'fnbak', 'phinbax')))
  
  par.but.df = reshape(data = par.but.df, idvar = c('conc', 'conc.eec', 'chem'), varying = 4:18, direction = 'long',
                       timevar = 'study', sep = '.', new.row.names = c(1:2500))  
  
  gbut.par<-ggplot(par.but.df) +
    theme_bw() +
    #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
    xlim(0,33.8) +
    ylim(-100,100) +
    labs(title = expression(paste('Parameter effect of Butralin on ', 'R'[0])),
         x = 'Butralin (ppb)',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc, y = r0, col = study), size = 1)
  
  gbut.par
  
#restrict to studies used to produce net effect  
  par.but.df$keep = 0
  par.but.df$keep[which(par.but.df$study %in% c('piCgaf', 'piMgaf', 'phinbax', 'mungaf'))] = 1
  par.but.keeps = subset(par.but.df, keep == 1)
  
  gbut.par.keeps<-ggplot(par.but.keeps) +
    theme_bw() +
    theme(legend.position = 'top', legend.direction = 'horizontal') +
    xlim(0,2) +
    ylim(-100,100) +
    labs(x = 'Butralin (EEC)',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc.eec, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc.eec, y = r0, col = study), size = 1)
  
  gbut.par.keeps
  
#Glyphosate ################
  load('Review_models/Savio/Glyphosate/gly5000sims_r0med_2017-12-12.Rdata')
  load('Review_models/Savio/Glyphosate/gly5000sims_r0.75_2017-12-12.Rdata')
  load('Review_models/Savio/Glyphosate/gly5000sims_r0.25_2017-12-12.Rdata')
  load('Review_models/Savio/Glyphosate/gly5000sims_pars_2017-12-12.Rdata')
  
  gly.df = data.frame(conc = rep(seq(0, 7400, length.out = 500), 2),
                      conc.eec = rep(seq(0, 7400, length.out = 500), 2)/3700,
                      chem = 'Glyphosate',
                      botup = c(rep('Absent', 500), rep('Present', 500)),
                      med = c(gly.med.r0[,11], gly.med.r0[,12]),
                      iqr.25 = c(gly.med.r0.25[,11], gly.med.r0.25[,12]),
                      iqr.75 = c(gly.med.r0.75[,11], gly.med.r0.75[,12]))
  
#summary of net effect at EEC
  gly.df[which(round(gly.df$conc.eec, digits = 2) == 1.00),][4,]
  
  ggly<-ggplot(gly.df) +
          theme_bw() +
          #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
          xlim(0,7400) +
          ylim(-100,100) +
          labs(title = expression(paste('Net effect of Glyphosate on ', 'R'[0])),
               x = 'Glyphosate (ppb)',
               y = expression(paste(Delta, R[0], ' (%)'))) +
          geom_hline(yintercept = 0, lty=2) + 
          geom_ribbon(aes(x = conc, ymin = iqr.25, ymax = iqr.75, fill = botup), alpha = 0.4) +
          geom_line(aes(x = conc, y = med, col = botup), size = 1)
  ggly    

#Individual parameter effects    
par.gly.df = data.frame(conc = seq(0, 7400, length.out = 500),
                        conc.eec = seq(0, 7400, length.out = 500)/3700,
                        chem = 'Glyphosate')

par.gly.df = cbind(par.gly.df, gly.med.r0[,c(1:8)], gly.med.r0.25[,c(1:8)], gly.med.r0.75[,c(1:8)])

colnames(par.gly.df)[c(4:27)] = c(paste0('r0.',c('piMgaf', 'piCgaf', 'mungaf', 'fNgaf','munbak', 'fnbak','munons', 'phinbax')),
                                  paste0('r025.' ,c('piMgaf', 'piCgaf', 'mungaf', 'fNgaf','munbak', 'fnbak','munons','phinbax')),
                                  paste0('r075.' ,c('piMgaf', 'piCgaf', 'mungaf', 'fNgaf','munbak', 'fnbak','munons','phinbax')))

par.gly.df = reshape(data = par.gly.df, idvar = c('conc', 'conc.eec', 'chem'), varying = 4:27, direction = 'long',
                     timevar = 'study', sep = '.', new.row.names = c(1:4000))  

  ggly.par<-ggplot(par.gly.df) +
    theme_bw() +
    #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
    xlim(0,7400) +
    ylim(-100,100) +
    labs(title = expression(paste('Parameter effect of Glyphosate on ', 'R'[0])),
         x = 'Glyphosate (ppb)',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc, y = r0, col = study), size = 1)
  
  ggly.par
  
#restrict to studies used to produce net effect  
  par.gly.df$keep = 0
  par.gly.df$keep[which(par.gly.df$study %in% c('piCgaf', 'piMgaf', 'phinbax', 'mungaf'))] = 1
  par.gly.keeps = subset(par.gly.df, keep == 1)
  
  ggly.par.keeps<-ggplot(par.gly.keeps) +
    theme_bw() +
    theme(legend.position = 'top', legend.direction = 'horizontal') +
    xlim(0,2) +
    ylim(-100,100) +
    labs(x = 'Glyphosate (EEC)',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc.eec, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc.eec, y = r0, col = study), size = 1)
  
  ggly.par.keeps
  
#Chlorpyrifos ################
  load('Review_models/Savio/Chlorpyrifos/ch5000sims_r0med_2017-12-12.Rdata')
  load('Review_models/Savio/Chlorpyrifos/ch5000sims_r0.75_2017-12-12.Rdata')
  load('Review_models/Savio/Chlorpyrifos/ch5000sims_r0.25_2017-12-12.Rdata')
  #load('Review_models/Savio/Chlorpyrifos/ch5000sims_pars_2017-12-12.Rdata')

  ch.df = data.frame(conc = seq(0,128, length.out = 500),
                     conc.eec = seq(0,128, length.out = 500)/64,
                     chem = 'Chlorpyrifos',
                     med = chlor.med.r0[,10],
                     iqr.25 = chlor.med.r0.25[,10],
                     iqr.75 = chlor.med.r0.75[,10])
  
#summary of net effect at EEC
  ch.df[which(round(ch.df$conc.eec, digits = 2) == 1.00),][2,]
  
  gch<-ggplot(ch.df) +
          theme_bw() +
          #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
          xlim(0,128) +
          ylim(-100,100) +
          labs(title = expression(paste('Net effect of Chlorpyrifos on ', 'R'[0])),
               x = 'Chlorpyrifos (ppb)',
               y = expression(paste(Delta, R[0], ' (%)'))) +
          geom_hline(yintercept = 0, lty=2) + 
          geom_ribbon(aes(x = conc, ymin = iqr.25, ymax = iqr.75), alpha = 0.4) +
          geom_line(aes(x = conc, y = med), size = 1)
  gch      

#Individual parameter effects    
par.ch.df = data.frame(conc = seq(0, 128, length.out = 500),
                        conc.eec = seq(0, 128, length.out = 500)/64,
                        chem = 'Chlorpyrifos')

par.ch.df = cbind(par.ch.df, chlor.med.r0[,c(1:9)], chlor.med.r0.25[,c(1:9)], chlor.med.r0.75[,c(1:9)])

colnames(par.ch.df)[c(4:30)] = c(paste0('r0.',c('muPhal', 'muPsat', 'psisat', 'piChash','piMhash', 'muNhash','muNibr', 'fNhash', 'fNibr')),
                                  paste0('r025.',c('muPhal', 'muPsat', 'psisat', 'piChash','piMhash', 'muNhash','muNibr', 'fNhash', 'fNibr')),
                                  paste0('r075.',c('muPhal', 'muPsat', 'psisat', 'piChash','piMhash', 'muNhash','muNibr', 'fNhash', 'fNibr')))

par.ch.df = reshape(data = par.ch.df, idvar = c('conc', 'conc.eec', 'chem'), varying = 4:30, direction = 'long',
                     timevar = 'study', sep = '.', new.row.names = c(1:4500))  

gch.par<-ggplot(par.ch.df) +
  theme_bw() +
  #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  xlim(0,128) +
  ylim(-100,100) +
  labs(title = expression(paste('Parameter effect of Chlorpyrifos on ', 'R'[0])),
       x = 'Chlorpyrifos (ppb)',
       y = expression(paste(Delta, R[0], ' (%)'))) +
  geom_hline(yintercept = 0, lty=2) + 
  geom_ribbon(aes(x = conc, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
  geom_line(aes(x = conc, y = r0, col = study), size = 1)

gch.par

#restrict to studies used to produce net effect  
  par.ch.df$keep = 0
  par.ch.df$keep[which(par.ch.df$study %in% c('piChash', 'piMhash', 'muPhal', 'fNibr', 'muNhash'))] = 1
  par.ch.keeps = subset(par.ch.df, keep == 1)
  
  gch.par.keeps<-ggplot(par.ch.keeps) +
    theme_bw() +
    theme(legend.position = 'top', legend.direction = 'horizontal') +
    xlim(0,2) +
    ylim(-100,100) +
    labs(x = 'Chlorpyrifos (EEC)',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc.eec, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc.eec, y = r0, col = study), size = 1)
  
  gch.par.keeps

#Malathion ################
  load('Review_models/Savio/Malathion/mal5000sims_r0med_2017-12-12.Rdata')
  load('Review_models/Savio/Malathion/mal5000sims_r0.75_2017-12-12.Rdata')
  load('Review_models/Savio/Malathion/mal5000sims_r0.25_2017-12-12.Rdata')
  load('Review_models/Savio/Malathion/mal5000sims_pars_2017-12-12.Rdata')

  mal.df = data.frame(conc = seq(0, 36.8, length.out = 500),
                      conc.eec = seq(0, 36.8, length.out = 500)/18.4,
                      chem = 'Malathion',
                      med = mal.med.r0[,8],
                      iqr.25 = mal.med.r0.25[,8],
                      iqr.75 = mal.med.r0.75[,8])
  
#summary of net effect at EEC
  mal.df[which(round(mal.df$conc.eec, digits = 2) == 1.00),][2,]
  
  gmal<-ggplot(mal.df) +
          theme_bw() +
          #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
          xlim(0,37) +
          ylim(-100,100) +
          labs(title = expression(paste('Net effect of Malathion on ', 'R'[0])),
               x = 'Malathion (ppb)',
               y = expression(paste(Delta, R[0], ' (%)'))) +
          geom_hline(yintercept = 0, lty=2) + 
          geom_ribbon(aes(x = conc, ymin = iqr.25, ymax = iqr.75), alpha = 0.4) +
          geom_line(aes(x = conc, y = med), size = 1)
  gmal 
 
#Individual parameter effects    
par.mal.df = data.frame(conc = seq(0, 36.8, length.out = 500),
                       conc.eec = seq(0, 36.8, length.out = 500)/18.4,
                       chem = 'Malathion')

par.mal.df = cbind(par.mal.df, mal.med.r0[,c(1:7)], mal.med.r0.25[,c(1:7)], mal.med.r0.75[,c(1:7)])

colnames(par.mal.df)[c(4:24)] = c(paste0('r0.',c('piCtch', 'piMtch', 'fNbak', 'munbak', 'muPhal','muNtch', 'fNtch')),
                                 paste0('r025.',c('piCtch', 'piMtch', 'fNbak', 'munbak', 'muPhal','muNtch', 'fNtch')),
                                 paste0('r075.',c('piCtch', 'piMtch', 'fNbak', 'munbak', 'muPhal','muNtch', 'fNtch')))

par.mal.df = reshape(data = par.mal.df, idvar = c('conc', 'conc.eec', 'chem'), varying = 4:24, direction = 'long',
                    timevar = 'study', sep = '.', new.row.names = c(1:3500))  

  gmal.par<-ggplot(par.mal.df) +
    theme_bw() +
    #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
    xlim(0,36.8) +
    ylim(-100,100) +
    labs(title = expression(paste('Parameter effect of Malathion on ', 'R'[0])),
         x = 'Malathion (ppb)',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc, y = r0, col = study), size = 1)
  
  gmal.par

#restrict to studies used to produce net effect  
  par.mal.df$keep = 0
  par.mal.df$keep[which(par.mal.df$study %in% c('piCtch', 'piMtch', 'muPhal', 'fNtch', 'muNtch'))] = 1
  par.mal.keeps = subset(par.mal.df, keep == 1)
  
  gmal.par.keeps<-ggplot(par.mal.keeps) +
    theme_bw() +
    theme(legend.position = 'top', legend.direction = 'horizontal') +
    xlim(0,2) +
    ylim(-100,100) +
    labs(x = 'Malathion (EEC)',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc.eec, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc.eec, y = r0, col = study), size = 1)
  
  gmal.par.keeps
  
#Profenofos ################
  load('Review_models/Savio/Profenofos/pr5000sims_r0med_2017-12-12.Rdata')
  load('Review_models/Savio/Profenofos/pr5000sims_r0.75_2017-12-12.Rdata')
  load('Review_models/Savio/Profenofos/pr5000sims_r0.25_2017-12-12.Rdata')
  load('Review_models/Savio/Profenofos/pr5000sims_pars_2017-12-12.Rdata')
  
  prof.df = data.frame(conc = seq(0, 6.02, length.out = 500),
                      conc.eec = seq(0, 6.02, length.out = 500)/3.01,
                      chem = 'Profenofos',
                      med = prof.med.r0[,5],
                      iqr.25 = prof.med.r0.25[,5],
                      iqr.75 = prof.med.r0.75[,5])
  
#summary of net effect at EEC
  prof.df[which(round(prof.df$conc.eec, digits = 2) == 1.00),][2,]
  
  gprof<-ggplot(prof.df) +
          theme_bw() +
          #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
          xlim(0,6.1) +
          ylim(-100,150) +
          labs(title = expression(paste('Net effect of Profenofos on ', 'R'[0])),
               x = 'Profenofos (ppb)',
               y = expression(paste(Delta, R[0], ' (%)'))) +
          geom_hline(yintercept = 0, lty=2) + 
          geom_ribbon(aes(x = conc, ymin = iqr.25, ymax = iqr.75), alpha = 0.4) +
          geom_line(aes(x = conc, y = med), size = 1)
  gprof      
  
#Individual parameter effects    
par.prof.df = data.frame(conc = seq(0, 6.02, length.out = 500),
                        conc.eec = seq(0, 6.02, length.out = 500)/3.01,
                        chem = 'Profenofos')

par.prof.df = cbind(par.prof.df, prof.med.r0[,c(1:4)], prof.med.r0.25[,c(1:4)], prof.med.r0.75[,c(1:4)])

colnames(par.prof.df)[c(4:15)] = c(paste0('r0.',c('muPsat','piChash', 'piMhash', 'muNhash')),
                                  paste0('r025.',c('muPsat','piChash', 'piMhash', 'muNhash')),
                                  paste0('r075.',c('muPsat','piChash', 'piMhash', 'muNhash')))

par.prof.df = reshape(data = par.prof.df, idvar = c('conc', 'conc.eec', 'chem'), varying = 4:15, direction = 'long',
                     timevar = 'study', sep = '.', new.row.names = c(1:2000))  

gprof.par<-ggplot(par.prof.df) +
  theme_bw() +
  #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  xlim(0,6.02) +
  ylim(-100,100) +
  labs(title = expression(paste('Parameter effect of Profenofos on ', 'R'[0])),
       x = 'Profenofos (ppb)',
       y = expression(paste(Delta, R[0], ' (%)'))) +
  geom_hline(yintercept = 0, lty=2) + 
  geom_ribbon(aes(x = conc, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
  geom_line(aes(x = conc, y = r0, col = study), size = 1)

gprof.par

#restrict to studies used to produce net effect  
  par.prof.df$keep = 0
  par.prof.df$keep[which(par.prof.df$study %in% c('piChash', 'piMhash', 'muPsat', 'muNhash'))] = 1
  par.prof.keeps = subset(par.prof.df, keep == 1)
  
  gprof.par.keeps<-ggplot(par.prof.keeps) +
    theme_bw() +
    theme(legend.position = 'top', legend.direction = 'horizontal') +
    xlim(0,2) +
    ylim(-100,100) +
    labs(x = 'Profenofos (EEC)',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc.eec, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc.eec, y = r0, col = study), size = 1)
  
  gprof.par.keeps

#Plot trajectories of herbicides and insecticides together############
#Herbicides ###############
#Net effects #########

hrbdf = rbind(atr.df, but.df, gly.df)  
  hrbdf = subset(hrbdf, botup == 'Present')
  
  ghrb<-ggplot(hrbdf) +
          theme_bw() +
          theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank()) +
          ylim(-100,150) +
          labs(title = expression(paste('Net effect of herbicides on ', 'R'[0])),
               x = 'Herbicide concentration',
               y = expression(paste(Delta, R[0], ' (%)'))) +
          scale_x_continuous(breaks = c(0,1,2), labels = c('0', 'EEC', '2EEC')) +
          geom_hline(yintercept = 0, lty=2) + 
          geom_ribbon(aes(x = conc.eec, ymin = iqr.25, ymax = iqr.75, fill = chem), alpha = 0.4) +
          geom_line(aes(x = conc.eec, y = med, col = chem), size = 1)
  #tiff('~/RemaisWork/Schisto/Agro_Review/Figures/r0_ConcRange/herbicides_net_response.tiff', 
  #     width = 2550, height = 1650, res = 300)
  ghrb   
  
  #dev.off()
  
#Parameter effects ###########
hrbdf.par = rbind(par.atr.df, par.but.df, par.gly.df)  

  ghrb.par<-ggplot(hrbdf.par) +
    theme_bw() +
    facet_grid(. ~ chem, scales = 'free_x') +
    theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank()) +
    ylim(-100,100) +
    labs(title = expression(paste('Parameter effects of herbicides on ', 'R'[0])),
         x = 'Herbicide concentration',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    scale_x_continuous(breaks = c(0,1,2), labels = c('0', 'EEC', '2EEC')) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc.eec, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc.eec, y = r0, col = study), size = 1)

  windows(width = 30, height = 15)
  
  ghrb.par  
  
  
#Kept parameter effects along with net effects ###########
  #first make master dataframe
  hrbdf.master = rbind(par.atr.keeps, par.but.keeps, par.gly.keeps)[,-8] 
    hrbdf.master$effect = 'Direct'
    hrbdf.master$effect[which(hrbdf.master$study == "phinbax")] = "Bottom Up"
    
  hrbdf.edit = hrbdf
  colnames(hrbdf.edit)[4:7] = colnames(hrbdf.master)[4:7] #make colnames compatible with master df
  hrbdf.edit$study = 'Combined'
  hrbdf.edit$effect = 'Net'
  
  hrbdf.master = rbind(hrbdf.master, hrbdf.edit)
  
#plot as matrix with chemicals along columns and component/net effects along rows  
  ghrb.master<-ggplot(hrbdf.master) +
    theme_bw() +
    facet_grid(effect ~ chem, switch = 'y') +
    theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(),
          strip.background = element_blank(), strip.text = element_text(size = 15),
          axis.title = element_text(size = 12)) +
    ylim(-100,100) +
    labs(x = 'Herbicide concentration (EEC)',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    #ggtitle('Component and net effects of herbicides') +
    scale_x_continuous(breaks = c(0.1,1,2), labels = c('0.1', '1', '2'), trans = "log") +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc.eec, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc.eec, y = r0, col = study), size = 1) +
    scale_color_manual(labels = c("Combined",
                                  expression(paste('Abdel-Ghaffar et al. 2016 (', mu[N],')', sep = '')),
                                  expression(paste('Omran & Salama 2013 (', mu[N], ')', sep = '')),
                                  expression(paste('Baxter et al. 2011 (', Phi[N], ')', sep = '')),
                                  expression(paste('Abdel-Ghaffar et al. 2016 (', pi[C], ')', sep = '')),
                                  expression(paste('Rohr et al. 2008 (', pi[C], ')', sep = '')),
                                  expression(paste('Abdel-Ghaffar et al. 2016 (', pi[M], ')', sep = ''))),
                       values = c("Combined" = 'black', 
                                  "mungaf" = 'red', 
                                  "munons" = 'grey', 
                                  "phinbax" = 'darkgreen', 
                                  "piCgaf" = 'violet', 
                                  "piCror" = 'purple', 
                                  "piMgaf" = 'gold')) +
    scale_fill_manual(labels = c("Combined",
                                  expression(paste('Abdel-Ghaffar et al. 2016 (', mu[N],')', sep = '')),
                                  expression(paste('Omran & Salama 2013 (', mu[N], ')', sep = '')),
                                  expression(paste('Baxter et al. 2011 (', Phi[N], ')', sep = '')),
                                  expression(paste('Abdel-Ghaffar et al. 2016 (', pi[C], ')', sep = '')),
                                  expression(paste('Rohr et al. 2008 (', pi[C], ')', sep = '')),
                                  expression(paste('Abdel-Ghaffar et al. 2016 (', pi[M], ')', sep = ''))),
                       values = c("Combined" = 'black', 
                                  "mungaf" = 'red', 
                                  "munons" = 'grey', 
                                  "phinbax" = 'darkgreen', 
                                  "piCgaf" = 'violet', 
                                  "piCror" = 'purple', 
                                  "piMgaf" = 'gold'))
  
  windows(width = 25, height = 25)
  
  ghrb.master 
  
#Insecticides ###############
#Net effects ##########  
insdf = rbind(ch.df, mal.df, prof.df)  

  gins<-ggplot(insdf) +
          theme_bw() +
          theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank()) +
          ylim(-100,150) +
          labs(title = expression(paste('Net effect of insecticides on ', 'R'[0])),
               x = 'Insecticide concentration',
               y = expression(paste(Delta, R[0], ' (%)'))) +
          scale_x_continuous(breaks = c(0,1,2), labels = c('0', 'EEC', '2EEC')) +
          geom_hline(yintercept = 0, lty=2) + 
          geom_ribbon(aes(x = conc.eec, ymin = iqr.25, ymax = iqr.75, fill = chem), alpha = 0.4) +
          geom_line(aes(x = conc.eec, y = med, col = chem), size = 1)
  
  #tiff('~/RemaisWork/Schisto/Agro_Review/Figures/r0_ConcRange/insecticides_net_response.tiff', 
  #     width = 2550, height = 1650, res = 300)
     
  gins
  
  #dev.off()
  
#Parameter effects #############
  insdf.par = rbind(par.ch.df, par.mal.df, par.prof.df)  
  
  gins.par<-ggplot(insdf.par) +
    theme_bw() +
    facet_grid(. ~ chem, scales = 'free_x') +
    theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank()) +
    ylim(-100,100) +
    labs(title = expression(paste('Parameter effects of insecticides on ', 'R'[0])),
         x = 'Insecticide concentration',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    scale_x_continuous(breaks = c(0,1,2), labels = c('0', 'EEC', '2EEC')) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc.eec, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc.eec, y = r0, col = study), size = 1)
  
  windows(width = 30, height = 15)
  
  gins.par 
  
#Kept parameter effects along with net effects ###########
  #first make master dataframe
  insdf.master = rbind(par.ch.keeps, par.mal.keeps, par.prof.keeps)[,-8] 
  insdf.master$effect = 'Direct'
    insdf.master$effect[grep("muP", insdf.master$study)] = "Top Down"
  
  insdf.edit = insdf
  insdf.edit$study = 'Combined'
  insdf.edit = insdf.edit[c(1:3,7,4:6)]
  colnames(insdf.edit)[4:7] = colnames(insdf.master)[4:7] #make colnames compatible with master df
  insdf.edit$effect = 'Net'
  
  insdf.master = rbind(insdf.master, insdf.edit)
  insdf.master$effect = factor(insdf.master$effect, levels = c("Top Down", "Direct", "Net"))
  
#plot as matrix with chemicals along columns and component/net effects along rows  
  gins.master<-ggplot(insdf.master) +
    theme_bw() +
    facet_grid(effect ~ chem, switch = 'y') +
    theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(),
          strip.background = element_blank(), strip.text = element_text(size = 15),
          axis.title = element_text(size = 12)) +
    ylim(-100,100) +
    labs(x = 'Insecticide concentration (EEC)',
         y = expression(paste(Delta, R[0], ' (%)'))) +
    #ggtitle('Component and net effects of insecticides') +
    scale_x_continuous(breaks = c(0.1,1,2), labels = c('0.1', '1', '2'), trans = "log") +
    coord_trans(y = "log")
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc.eec, ymin = r025, ymax = r075, fill = study), alpha = 0.4) +
    geom_line(aes(x = conc.eec, y = r0, col = study), size = 1) +
    scale_color_manual(labels = c("Combined",
                                  expression(paste('Ibrahim et al. 1992 (', f[N],')', sep = '')),
                                  expression(paste('Tchounwou et al. 1991 (', f[N], ')', sep = '')),
                                  expression(paste('Hasheesh & Mohamed 2011 (', mu[N], ')', sep = '')),
                                  expression(paste('Tchounwou et al. 1991A (', mu[N], ')', sep = '')),
                                  expression(paste('Halstead et al. 2015 (', mu[P], ')', sep = '')),
                                  expression(paste('Satapornvanit et al. 2009 (', mu[P], ')', sep = '')),
                                  expression(paste('Hasheesh & Mohamed 2011 (', pi[C], ')', sep = '')),
                                  expression(paste('Tchounwou et al 1992 (', pi[C], ')', sep = '')),
                                  expression(paste('Hasheesh & Mohamed 2011 (', pi[M], ')', sep = '')),
                                  expression(paste('Tchounwou et al 1991B (', pi[M], ')', sep = ''))),
                       values = c("Combined" = 'black', 
                                  "fNibr" = 'orange', 
                                  "fNtch" = 'grey', 
                                  "muNhash" = 'grey', 
                                  "muNtch" = 'grey', 
                                  "muPhal" = 'red', 
                                  "muPsat" = 'darkred',
                                  "piChash" = 'purple',
                                  "piCtch" = 'violet',
                                  "piMhash" = 'grey',
                                  "piMtch" = 'grey')) +
    scale_fill_manual(labels = c("Combined",
                                 expression(paste('Ibrahim et al. 1992 (', f[N],')', sep = '')),
                                 expression(paste('Tchounwou et al. 1991 (', f[N], ')', sep = '')),
                                 expression(paste('Hasheesh & Mohamed 2011 (', mu[N], ')', sep = '')),
                                 expression(paste('Tchounwou et al. 1991A (', mu[N], ')', sep = '')),
                                 expression(paste('Halstead et al. 2015 (', mu[P], ')', sep = '')),
                                 expression(paste('Satapornvanit et al. 2009 (', mu[P], ')', sep = '')),
                                 expression(paste('Hasheesh & Mohamed 2011 (', pi[C], ')', sep = '')),
                                 expression(paste('Tchounwou et al 1992 (', pi[C], ')', sep = '')),
                                 expression(paste('Hasheesh & Mohamed 2011 (', pi[M], ')', sep = '')),
                                 expression(paste('Tchounwou et al 1991B (', pi[M], ')', sep = ''))),
                      values = c("Combined" = 'black', 
                                 "fNibr" = 'orange', 
                                 "fNtch" = 'grey', 
                                 "muNhash" = 'grey', 
                                 "muNtch" = 'grey', 
                                 "muPhal" = 'red', 
                                 "muPsat" = 'darkred',
                                 "piChash" = 'purple',
                                 "piCtch" = 'violet',
                                 "piMhash" = 'grey',
                                 "piMtch" = 'grey'))
  
  windows(width = 25, height = 25)
  
  gins.master 
  
  
#Combined herbicides and insecticides plot ############
hrbdf = subset(hrbdf, botup == 'Present')  
  hrbdf = hrbdf[,-4]
  hrbdf$Class = 'Herbicides'
  insdf$Class = 'Insecticides'
 
compdf = rbind(hrbdf, insdf) 

gcomp = ggplot(data = compdf) +
  facet_grid(. ~ Class, scales = 'free_x') +
  theme_bw() +
  theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ylim(-100,100) +
  labs(x = 'Agrochemical conc (EEC)',
       y = expression(paste(Delta, R[0], '%'))) +
  ggtitle(label = expression(paste('Net agrochemical effects on ', 'R'[0]))) +
  geom_hline(yintercept = 0, lty=2) + 
  geom_ribbon(aes(x = conc.eec, ymin = iqr.25, ymax = iqr.75, fill = chem), alpha = 0.4) +
  geom_line(aes(x = conc.eec, y = med, col = chem), size = 1)

windows(width = 30, height = 15)

gcomp

tiff('~/RemaisWork/Schisto/Agro_Review/Figures/r0_ConcRange/ins&herb_net_response.tiff', 
     width = 3300, height = 1650, res = 300)

gcomp

dev.off()
