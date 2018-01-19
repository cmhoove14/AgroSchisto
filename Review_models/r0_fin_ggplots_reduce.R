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
load('Review_models/Savio/Atrazine/atr5000sims_r0med_2017-11-28.Rdata')
load('Review_models/Savio/Atrazine/atr5000sims_r0.75_2017-11-28.Rdata')
load('Review_models/Savio/Atrazine/atr5000sims_r0.25_2017-11-28.Rdata')
load('Review_models/Savio/Atrazine/atr5000sims_pars_2017-11-28.Rdata')
 
#Net effect 
  atr.df = data.frame(conc = seq(0,204, length.out = 500),
                      conc.eec = seq(0,204, length.out = 500)/102,
                      chem = 'Atrazine',
                      botup = 'Present',
                      med = atr.med.r0[,7],
                      iqr.25 = atr.med.r0.25[,7],
                      iqr.75 = atr.med.r0.75[,7])
  
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
#Butralin ################
load('Review_models/Savio/Butralin/but5000sims_r0med_2017-11-29.Rdata')
load('Review_models/Savio/Butralin/but5000sims_r0.75_2017-11-29.Rdata')
load('Review_models/Savio/Butralin/but5000sims_r0.25_2017-11-29.Rdata')
load('Review_models/Savio/Butralin/but5000sims_pars_2017-11-29.Rdata')

  but.df = data.frame(conc = rep(seq(0, 33.8, length.out = 500), 2),
                      conc.eec = rep(seq(0, 33.8, length.out = 500), 2)/16.9,
                      chem = 'Butralin',
                      botup = c(rep('Absent', 500), rep('Present', 500)),
                      med = c(but.med.r0[,8], but.med.r0[,9]),
                      iqr.25 = c(but.med.r0.25[,8], but.med.r0.25[,9]),
                      iqr.75 = c(but.med.r0.75[,8], but.med.r0.75[,9]))
  
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
  
#Glyphosate ################
  load('Review_models/Savio/Glyphosate/gly5000sims_r0med_2017-11-29.Rdata')
  load('Review_models/Savio/Glyphosate/gly5000sims_r0.75_2017-11-29.Rdata')
  load('Review_models/Savio/Glyphosate/gly5000sims_r0.25_2017-11-29.Rdata')
  load('Review_models/Savio/Glyphosate/gly5000sims_pars_2017-11-29.Rdata')
  
  gly.df = data.frame(conc = rep(seq(0, 7400, length.out = 500), 2),
                      conc.eec = rep(seq(0, 7400, length.out = 500), 2)/3700,
                      chem = 'Glyphosate',
                      botup = c(rep('Absent', 500), rep('Present', 500)),
                      med = c(gly.med.r0[,11], gly.med.r0[,12]),
                      iqr.25 = c(gly.med.r0.25[,11], gly.med.r0.25[,12]),
                      iqr.75 = c(gly.med.r0.75[,11], gly.med.r0.75[,12]))
  
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

#Chlorpyrifos ################
  load('Review_models/Savio/Chlorpyrifos/ch5000sims_r0med_2017-11-29.Rdata')
  load('Review_models/Savio/Chlorpyrifos/ch5000sims_r0.75_2017-11-29.Rdata')
  load('Review_models/Savio/Chlorpyrifos/ch5000sims_r0.25_2017-11-29.Rdata')
  load('Review_models/Savio/Chlorpyrifos/ch5000sims_pars_2017-11-29.Rdata')

  ch.df = data.frame(conc = seq(0,128, length.out = 500),
                     conc.eec = seq(0,128, length.out = 500)/64,
                     chem = 'Chlorpyrifos',
                     med = chlor.med.r0[,10],
                     iqr.25 = chlor.med.r0.25[,10],
                     iqr.75 = chlor.med.r0.75[,10])

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

#get rid of munHash to see results better
par.ch.sub = subset(par.ch.df, study != 'muNhash')

gch.par.sub<-ggplot(par.ch.sub) +
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

gch.par.sub

#Malathion ################
  load('Review_models/Savio/Malathion/mal5000sims_r0med_2017-11-28.Rdata')
  load('Review_models/Savio/Malathion/mal5000sims_r0.75_2017-11-28.Rdata')
  load('Review_models/Savio/Malathion/mal5000sims_r0.25_2017-11-28.Rdata')
  load('Review_models/Savio/Malathion/mal5000sims_pars_2017-11-28.Rdata')

  mal.df = data.frame(conc = seq(0, 36.8, length.out = 500),
                      conc.eec = seq(0, 36.8, length.out = 500)/18.4,
                      chem = 'Malathion',
                      med = mal.med.r0[,8],
                      iqr.25 = mal.med.r0.25[,8],
                      iqr.75 = mal.med.r0.75[,8])
  
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
  
#Profenofos ################
  load('Review_models/Savio/Profenofos/pr5000sims_r0med_2017-11-28.Rdata')
  load('Review_models/Savio/Profenofos/pr5000sims_r0.75_2017-11-28.Rdata')
  load('Review_models/Savio/Profenofos/pr5000sims_r0.25_2017-11-28.Rdata')
  load('Review_models/Savio/Profenofos/pr5000sims_pars_2017-11-28.Rdata')
  
  prof.df = data.frame(conc = seq(0, 6.02, length.out = 500),
                      conc.eec = seq(0, 6.02, length.out = 500)/3.01,
                      chem = 'Profenofos',
                      med = prof.med.r0[,5],
                      iqr.25 = prof.med.r0.25[,5],
                      iqr.75 = prof.med.r0.75[,5])
  
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
  
#Plot trajectories of herbicides and insecticides together############
#Herbicides ###############
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
  ghrb      
  
#Insecticides ###############
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
  gins   
  
  
  
  
load('Review_models/Savio/Chlorpyrifos/r0_Ch_6-30-17.RData')
load('Review_models/Savio/Glyphosate/r0_Gly_6-30-17.RData')

  keepers = c('conc.atr','conc.atr.means.r0', 'conc.atr.sds.r0','r0.atr.fix.n0',  
              'conc.ch', 'conc.ch.means.r0', 'conc.ch.sds.r0',
              'conc.mal', 'conc.mal.means.r0', 'conc.mal.sds.r0','r0.mal.fix.n0', 
              'conc.gly', 'conc.gly.means.r0', 'conc.gly.sds.r0',
              'parameters','r0.He', 'r0.In','nil1', 'nil0')
  
  rm(list = setdiff(ls(), keepers))

require(rootSolve)
require(ggplot2)
  rec.atr = 102   #From Halstead et al Nature Comm 2017 table S6
  rec.mal = 583    #From Halstead et al Chemosphere 2015 table S1
  rec.ch = 64   #From Halstead et al Nature Comm 2017 table S6
  rec.gly = 3700

  
  
#Make gg compatible df for atrazine results #############
atrm = matrix(ncol = 5, nrow = length(conc.atr)*4)
  colnames(atrm) = c('conc','Modeled','val','UL', 'LL')

  atrm[,1] = rep(conc.atr, 4)
  atrm[,2] = c(rep('Cercarial mortality', length(conc.atr)),     #From Rohr 2008 study
               rep('Snail mortality', length(conc.atr)),         #from Omran & Salama study
               rep('Snail carrying capacity', length(conc.atr)), #from Rohr reanalysis of Baxter study
               rep('all', length(conc.atr)))                     #combined of the three above
#Fill cercarial mortality From Rohr 2008 study
  atrm[c(1:length(conc.atr)),3] = conc.atr.means.r0[,4] - r0.He(0)[3]
  atrm[c(1:length(conc.atr)),4] = conc.atr.means.r0[,4] - r0.He(0)[3] + conc.atr.sds.r0[,4]
  atrm[c(1:length(conc.atr)),5] = conc.atr.means.r0[,4] - r0.He(0)[3] - conc.atr.sds.r0[,4]
#Fill snail mortality from Omran & Salama study
  atrm[c((length(conc.atr)+1):(2*length(conc.atr))),3] = conc.atr.means.r0[,7] - r0.He(0)[3]
  atrm[c((length(conc.atr)+1):(2*length(conc.atr))),4] = conc.atr.means.r0[,7] - r0.He(0)[3] + conc.atr.sds.r0[,7]
  atrm[c((length(conc.atr)+1):(2*length(conc.atr))),5] = conc.atr.means.r0[,7] - r0.He(0)[3] - conc.atr.sds.r0[,7]
#Fill snail carrying capacity from Rohr reanalysis of Baxter study  
  atrm[c((2*length(conc.atr)+1):(3*length(conc.atr))),3] = conc.atr.means.r0[,6] - r0.He(0)[3]
  atrm[c((2*length(conc.atr)+1):(3*length(conc.atr))),4] = conc.atr.means.r0[,6] - r0.He(0)[3] + conc.atr.sds.r0[,6]
  atrm[c((2*length(conc.atr)+1):(3*length(conc.atr))),5] = conc.atr.means.r0[,6] - r0.He(0)[3] - conc.atr.sds.r0[,6]
#Fill combo of three  
  atrm[c((3*length(conc.atr)+1):(4*length(conc.atr))),3] = conc.atr.means.r0[,17] - r0.He(0)[3]
  atrm[c((3*length(conc.atr)+1):(4*length(conc.atr))),4] = conc.atr.means.r0[,17] - r0.He(0)[3] + conc.atr.sds.r0[,17]
  atrm[c((3*length(conc.atr)+1):(4*length(conc.atr))),5] = conc.atr.means.r0[,17] - r0.He(0)[3] - conc.atr.sds.r0[,17]
  
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
      theme(legend.position = 'bottom', legend.direction = 'horizontal') +
      xlim(0,1000) +
      ylim(-4,5) +
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
  malm[c(1:length(conc.mal)),3] = conc.mal.means.r0[,1] - r0.In(0)[3]
  malm[c(1:length(conc.mal)),4] = conc.mal.means.r0[,1] - r0.In(0)[3] + conc.mal.sds.r0[,1]
  malm[c(1:length(conc.mal)),5] = conc.mal.means.r0[,1] - r0.In(0)[3] - conc.mal.sds.r0[,1]
#fill miracidial mortality from Tchounwou91
  malm[c((length(conc.mal)+1):(length(conc.mal)*2)),3] = conc.mal.means.r0[,2] - r0.In(0)[3]
  malm[c((length(conc.mal)+1):(length(conc.mal)*2)),4] = conc.mal.means.r0[,2] - r0.In(0)[3] + conc.mal.sds.r0[,2]
  malm[c((length(conc.mal)+1):(length(conc.mal)*2)),5] = conc.mal.means.r0[,2] - r0.In(0)[3] - conc.mal.sds.r0[,2]
#fill predator mortality from Halstead15  
  malm[c((2*length(conc.mal)+1):(length(conc.mal)*3)),3] = conc.mal.means.r0[,5] - r0.In(0)[3]
  malm[c((2*length(conc.mal)+1):(length(conc.mal)*3)),4] = conc.mal.means.r0[,5] - r0.In(0)[3] + conc.mal.sds.r0[,5]
  malm[c((2*length(conc.mal)+1):(length(conc.mal)*3)),5] = conc.mal.means.r0[,5] - r0.In(0)[3] - conc.mal.sds.r0[,5]
#fill snail mortality from Tchounwou91  
  malm[c((3*length(conc.mal)+1):(length(conc.mal)*4)),3] = conc.mal.means.r0[,6] - r0.In(0)[3]
  malm[c((3*length(conc.mal)+1):(length(conc.mal)*4)),4] = conc.mal.means.r0[,6] - r0.In(0)[3] + conc.mal.sds.r0[,6]
  malm[c((3*length(conc.mal)+1):(length(conc.mal)*4)),5] = conc.mal.means.r0[,6] - r0.In(0)[3] - conc.mal.sds.r0[,6]
#fill snail reproduction reduction from Tchounwou92  
  malm[c((4*length(conc.mal)+1):(length(conc.mal)*5)),3] = conc.mal.means.r0[,7] - r0.In(0)[3]
  malm[c((4*length(conc.mal)+1):(length(conc.mal)*5)),4] = conc.mal.means.r0[,7] - r0.In(0)[3] + conc.mal.sds.r0[,7]
  malm[c((4*length(conc.mal)+1):(length(conc.mal)*5)),5] = conc.mal.means.r0[,7] - r0.In(0)[3] - conc.mal.sds.r0[,7]
#fill combo of all  
  malm[c((5*length(conc.mal)+1):(length(conc.mal)*6)),3] = conc.mal.means.r0[,13] - r0.In(0)[3]
  malm[c((5*length(conc.mal)+1):(length(conc.mal)*6)),4] = conc.mal.means.r0[,13] - r0.In(0)[3] + conc.mal.sds.r0[,13]
  malm[c((5*length(conc.mal)+1):(length(conc.mal)*6)),5] = conc.mal.means.r0[,13] - r0.In(0)[3] - conc.mal.sds.r0[,13]

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
    xlim(0,50000) +
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
chm = matrix(ncol = 5, nrow = length(conc.ch)*10)
colnames(chm) = c('conc','Modeled','val','UL', 'LL')

chm[,1] = rep(conc.ch, 10)
chm[,2] = c(rep('Predator mortality (Halstead)', length(conc.ch)),   #from Halstead 2015
            rep('Predator mortality (Satapornvanit)', length(conc.ch)),   #from Satapornvanit
            rep('Predator attack rate', length(conc.ch)), #from Satapornvanit 
            rep('Cercarial mortality', length(conc.ch)),  #from Hasheesh 2011
            rep('Miracidial mortality', length(conc.ch)), #from Hasheesh 2011
            rep('Snail mortality (Hasheesh)', length(conc.ch)),      #from Hasheesh 2011
            rep('Snail mortality (Ibrahim)', length(conc.ch)),      #from Hasheesh 2011
            rep('Snail reproduction (Hasheesh)', length(conc.ch)),      #from Hasheesh 2011
            rep('Snail reproduction (Ibrahim)', length(conc.ch)),   #From Ibrahim 1992
            rep('all', length(conc.ch)))                  #all of the above combined

#fill predator mortality from Halstead 2015 *******************************************************************************
  chm[c(1:length(conc.ch)),3] = conc.ch.means.r0[,1] - r0.In(0)[3]
  chm[c(1:length(conc.ch)),4] = conc.ch.means.r0[,1] - r0.In(0)[3] + conc.ch.sds.r0[,1]
  chm[c(1:length(conc.ch)),5] = conc.ch.means.r0[,1] - r0.In(0)[3] - conc.ch.sds.r0[,1]
#fill predator mortality from Satapornvanit ******************************************************************************
  chm[c((length(conc.ch)+1):(length(conc.ch)*2)),3] = conc.ch.means.r0[,2] - r0.In(0)[3]
  chm[c((length(conc.ch)+1):(length(conc.ch)*2)),4] = conc.ch.means.r0[,2] - r0.In(0)[3] + conc.ch.sds.r0[,2]
  chm[c((length(conc.ch)+1):(length(conc.ch)*2)),5] = conc.ch.means.r0[,2] - r0.In(0)[3] - conc.ch.sds.r0[,2]
#fill predator attack rate from Satapornvanit ******************************************************************************
  chm[c((2*length(conc.ch)+1):(length(conc.ch)*3)),3] = conc.ch.means.r0[,3] - r0.In(0)[3]
  chm[c((2*length(conc.ch)+1):(length(conc.ch)*3)),4] = conc.ch.means.r0[,3] - r0.In(0)[3] + conc.ch.sds.r0[,3]
  chm[c((2*length(conc.ch)+1):(length(conc.ch)*3)),5] = conc.ch.means.r0[,3] - r0.In(0)[3] - conc.ch.sds.r0[,3]
#fill cercarial mortality from Hasheesh 2011 ******************************************************************************
  chm[c((3*length(conc.ch)+1):(length(conc.ch)*4)),3] = conc.ch.means.r0[,4] - r0.In(0)[3]
  chm[c((3*length(conc.ch)+1):(length(conc.ch)*4)),4] = conc.ch.means.r0[,4] - r0.In(0)[3] + conc.ch.sds.r0[,4]
  chm[c((3*length(conc.ch)+1):(length(conc.ch)*4)),5] = conc.ch.means.r0[,4] - r0.In(0)[3] - conc.ch.sds.r0[,4]
#fill miracidial mortality from Hasheesh 2011 ******************************************************************************
  chm[c((4*length(conc.ch)+1):(length(conc.ch)*5)),3] = conc.ch.means.r0[,5] - r0.In(0)[3]
  chm[c((4*length(conc.ch)+1):(length(conc.ch)*5)),4] = conc.ch.means.r0[,5] - r0.In(0)[3] + conc.ch.sds.r0[,5]
  chm[c((4*length(conc.ch)+1):(length(conc.ch)*5)),5] = conc.ch.means.r0[,5] - r0.In(0)[3] - conc.ch.sds.r0[,5]
#fill snail mortality from Hasheesh 2011  ******************************************************************************
  chm[c((5*length(conc.ch)+1):(length(conc.ch)*6)),3] = conc.ch.means.r0[,6] - r0.In(0)[3]
  chm[c((5*length(conc.ch)+1):(length(conc.ch)*6)),4] = conc.ch.means.r0[,6] - r0.In(0)[3] + conc.ch.sds.r0[,6]
  chm[c((5*length(conc.ch)+1):(length(conc.ch)*6)),5] = conc.ch.means.r0[,6] - r0.In(0)[3] - conc.ch.sds.r0[,6]
#fill snail mortality from Ibrahim 1992  ******************************************************************************
  chm[c((6*length(conc.ch)+1):(length(conc.ch)*7)),3] = conc.ch.means.r0[,7] - r0.In(0)[3]
  chm[c((6*length(conc.ch)+1):(length(conc.ch)*7)),4] = conc.ch.means.r0[,7] - r0.In(0)[3] + conc.ch.sds.r0[,7]
  chm[c((6*length(conc.ch)+1):(length(conc.ch)*7)),5] = conc.ch.means.r0[,7] - r0.In(0)[3] - conc.ch.sds.r0[,7]
#fill snail reproduction reduction from Hasheesh 2011  ******************************************************************************
  chm[c((7*length(conc.ch)+1):(length(conc.ch)*8)),3] = conc.ch.means.r0[,8] - r0.In(0)[3]
  chm[c((7*length(conc.ch)+1):(length(conc.ch)*8)),4] = conc.ch.means.r0[,8] - r0.In(0)[3] + conc.ch.sds.r0[,8]
  chm[c((7*length(conc.ch)+1):(length(conc.ch)*8)),5] = conc.ch.means.r0[,8] - r0.In(0)[3] - conc.ch.sds.r0[,8]
#fill snail reproduction reduction from Ibrahim 1992  ******************************************************************************
  chm[c((8*length(conc.ch)+1):(length(conc.ch)*9)),3] = conc.ch.means.r0[,9] - r0.In(0)[3]
  chm[c((8*length(conc.ch)+1):(length(conc.ch)*9)),4] = conc.ch.means.r0[,9] - r0.In(0)[3] + conc.ch.sds.r0[,9]
  chm[c((8*length(conc.ch)+1):(length(conc.ch)*9)),5] = conc.ch.means.r0[,9] - r0.In(0)[3] - conc.ch.sds.r0[,9]
#fill combo of all  ******************************************************************************  
  chm[c((9*length(conc.ch)+1):(length(conc.ch)*10)),3] = conc.ch.means.r0[,19] - r0.In(0)[3]
  chm[c((9*length(conc.ch)+1):(length(conc.ch)*10)),4] = conc.ch.means.r0[,19] - r0.In(0)[3] + conc.ch.sds.r0[,19]
  chm[c((9*length(conc.ch)+1):(length(conc.ch)*10)),5] = conc.ch.means.r0[,19] - r0.In(0)[3] - conc.ch.sds.r0[,19]

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
  xlim(0,200) +
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


#Make gg compatible df for chlorpyrifos results #############     
glym = matrix(ncol = 5, nrow = length(conc.gly)*8)
colnames(glym) = c('conc','Modeled','val','UL', 'LL')

glym[,1] = rep(conc.gly, 8)
glym[,2] = c(rep('Miracidial mortality', length(conc.gly)),   #from Halstead 2015
            rep('Cercarial mortality', length(conc.gly)),   #from Satapornvanit
            rep('Snail mortality (Abdel-Ghaffar)', length(conc.gly)), #from Satapornvanit 
            rep('Snail reproduction (Abdel-Ghaffar)', length(conc.gly)),  #from Hasheesh 2011
            rep('Snail mortality (Bakry)', length(conc.gly)), #from Hasheesh 2011
            rep('Snail reproduction (Bakry)', length(conc.gly)),      #from Hasheesh 2011
            rep('Snail mortality (Omran & Salama)', length(conc.gly)),   #From Ibrahim 1992
            rep('all', length(conc.gly)))                  #combined

#fill miracidial mortality from Abdel Ghaffar  *******************************************************************************
glym[c(1:length(conc.gly)),3] = conc.gly.means.r0[,1] - r0.In(0)[3]
glym[c(1:length(conc.gly)),4] = conc.gly.means.r0[,1] - r0.In(0)[3] + conc.gly.sds.r0[,1]
glym[c(1:length(conc.gly)),5] = conc.gly.means.r0[,1] - r0.In(0)[3] - conc.gly.sds.r0[,1]
#fill cercarial mortality from Abdel Ghaffar  ******************************************************************************
glym[c((length(conc.gly)+1):(length(conc.gly)*2)),3] = conc.gly.means.r0[,2] - r0.In(0)[3]
glym[c((length(conc.gly)+1):(length(conc.gly)*2)),4] = conc.gly.means.r0[,2] - r0.In(0)[3] + conc.gly.sds.r0[,2]
glym[c((length(conc.gly)+1):(length(conc.gly)*2)),5] = conc.gly.means.r0[,2] - r0.In(0)[3] - conc.gly.sds.r0[,2]
#fill snail mortality from Abdel Ghaffar ******************************************************************************
glym[c((2*length(conc.gly)+1):(length(conc.gly)*3)),3] = conc.gly.means.r0[,3] - r0.In(0)[3]
glym[c((2*length(conc.gly)+1):(length(conc.gly)*3)),4] = conc.gly.means.r0[,3] - r0.In(0)[3] + conc.gly.sds.r0[,3]
glym[c((2*length(conc.gly)+1):(length(conc.gly)*3)),5] = conc.gly.means.r0[,3] - r0.In(0)[3] - conc.gly.sds.r0[,3]
#fill snail reproduction reduction from Abdel Ghaffar ******************************************************************************
glym[c((3*length(conc.gly)+1):(length(conc.gly)*4)),3] = conc.gly.means.r0[,4] - r0.In(0)[3]
glym[c((3*length(conc.gly)+1):(length(conc.gly)*4)),4] = conc.gly.means.r0[,4] - r0.In(0)[3] + conc.gly.sds.r0[,4]
glym[c((3*length(conc.gly)+1):(length(conc.gly)*4)),5] = conc.gly.means.r0[,4] - r0.In(0)[3] - conc.gly.sds.r0[,4]
#fill snail mortality from Bakry ******************************************************************************
glym[c((4*length(conc.gly)+1):(length(conc.gly)*5)),3] = conc.gly.means.r0[,5] - r0.In(0)[3]
glym[c((4*length(conc.gly)+1):(length(conc.gly)*5)),4] = conc.gly.means.r0[,5] - r0.In(0)[3] + conc.gly.sds.r0[,5]
glym[c((4*length(conc.gly)+1):(length(conc.gly)*5)),5] = conc.gly.means.r0[,5] - r0.In(0)[3] - conc.gly.sds.r0[,5]
#fill snail reproduction reduction from Bakry  ******************************************************************************
glym[c((5*length(conc.gly)+1):(length(conc.gly)*6)),3] = conc.gly.means.r0[,6] - r0.In(0)[3]
glym[c((5*length(conc.gly)+1):(length(conc.gly)*6)),4] = conc.gly.means.r0[,6] - r0.In(0)[3] + conc.gly.sds.r0[,6]
glym[c((5*length(conc.gly)+1):(length(conc.gly)*6)),5] = conc.gly.means.r0[,6] - r0.In(0)[3] - conc.gly.sds.r0[,6]
#fill snail mortality from Omran & Salama  ******************************************************************************
glym[c((6*length(conc.gly)+1):(length(conc.gly)*7)),3] = conc.gly.means.r0[,7] - r0.In(0)[3]
glym[c((6*length(conc.gly)+1):(length(conc.gly)*7)),4] = conc.gly.means.r0[,7] - r0.In(0)[3] + conc.gly.sds.r0[,7]
glym[c((6*length(conc.gly)+1):(length(conc.gly)*7)),5] = conc.gly.means.r0[,7] - r0.In(0)[3] - conc.gly.sds.r0[,7]
#fill combo of all  ******************************************************************************  
glym[c((7*length(conc.gly)+1):(length(conc.gly)*8)),3] = conc.gly.means.r0[,8] - r0.In(0)[3]
glym[c((7*length(conc.gly)+1):(length(conc.gly)*8)),4] = conc.gly.means.r0[,8] - r0.In(0)[3] + conc.gly.sds.r0[,8]
glym[c((7*length(conc.gly)+1):(length(conc.gly)*8)),5] = conc.gly.means.r0[,8] - r0.In(0)[3] - conc.gly.sds.r0[,8]

gly.gg<-as.data.frame(glym, stringsAsFactors = FALSE)
gly.gg$conc<-as.numeric(gly.gg$conc)
gly.gg$val<-as.numeric(gly.gg$val)
gly.gg$UL<-as.numeric(gly.gg$UL)
gly.gg$LL<-as.numeric(gly.gg$LL)

gly.ggpars = subset(gly.gg, Modeled != 'all') 
gly.ggpars$Chemical = 'Glyphosate'
gly.ggall = subset(gly.gg, Modeled == 'all')
gly.ggall$Chemical = 'Glyphosate'
gly.ggall$rec = gly.ggall$conc / rec.gly

#Plot glyphosate component effects############        
ggly<-ggplot(gly.ggpars) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank()) +
  xlim(0,4000) +
  ylim(-5,5) +
  labs(title = expression(paste('Component effects of Glyphosate on ', 'R'[0])),
       x = 'Glyphosate (ppb)',
       y = expression(paste(Delta, R[0]))) +
  geom_hline(yintercept = 0, lty=2) + 
  geom_ribbon(aes(x = conc, ymin = LL, ymax = UL, fill = Modeled), alpha = 0.4) +
  geom_line(aes(x = conc, y = val, col = Modeled), size = 1)

ggly    

#Plot atrazine, malathion, chlorpyrifos, and glyphosate component effects as facets ################
comp = rbind(atr.ggpars, mal.ggpars, ch.ggpars, gly.ggpars)
  gcomp = ggplot(data = comp) +
    facet_grid(. ~ Chemical, scales = 'free_x') +
    theme_bw() +
    theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ylim(-5,5) +
    labs(x = 'Agrochemical conc (ppb)',
         y = expression(paste(Delta, R[0]))) +
    ggtitle(label = expression(paste('Component agrochemical effects on ', 'R'[0]))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = conc, ymin = LL, ymax = UL, fill = Modeled), alpha = 0.4) +
    geom_line(aes(x = conc, y = val, col = Modeled), size = 1)
  
  windows(width = 30, height = 15)
  
  gcomp
    
#Plot atrazine, malathion, chlorpyrifos, and glyphosate combined effects together ###########
comb = rbind(atr.ggall, mal.ggall, ch.ggall, gly.ggall)
  gall<-ggplot(comb) +
    theme_bw() +
    theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = c(0,1,2), labels = c('0', 'EEC','2x EEC'), limits = c(0,2.0)) +
    ylim(-3,3) +
    labs(title = expression(paste('Combined agrochemical effects on ', 'R'[0])),
         x = 'Normalized Agrochemical Concentration',
         y = expression(paste(Delta, R[0]))) +
    geom_hline(yintercept = 0, lty=2) + 
    geom_ribbon(aes(x = rec, ymin = LL, ymax = UL, fill = Chemical), alpha = 0.4) +
    geom_line(aes(x = rec, y = val, col = Chemical), size = 1)
  
  windows(width = 30, height = 15)
  gall