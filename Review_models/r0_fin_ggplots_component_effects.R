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
load('Review_models/Savio/Butralin/but5000sims_r0med_2017-12-12.Rdata')
load('Review_models/Savio/Butralin/but5000sims_r0.75_2017-12-12.Rdata')
load('Review_models/Savio/Butralin/but5000sims_r0.25_2017-12-12.Rdata')
load('Review_models/Savio/Butralin/but5000sims_pars_2017-12-12.Rdata')


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
  load('Review_models/Savio/Glyphosate/gly5000sims_r0med_2017-12-12.Rdata')
  load('Review_models/Savio/Glyphosate/gly5000sims_r0.75_2017-12-12.Rdata')
  load('Review_models/Savio/Glyphosate/gly5000sims_r0.25_2017-12-12.Rdata')
  load('Review_models/Savio/Glyphosate/gly5000sims_pars_2017-12-12.Rdata')
  
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
  load('Review_models/Savio/Chlorpyrifos/ch5000sims_r0med_2017-12-12.Rdata')
  load('Review_models/Savio/Chlorpyrifos/ch5000sims_r0.75_2017-12-12.Rdata')
  load('Review_models/Savio/Chlorpyrifos/ch5000sims_r0.25_2017-12-12.Rdata')
  #load('Review_models/Savio/Chlorpyrifos/ch5000sims_pars_2017-12-12.Rdata')

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
  load('Review_models/Savio/Malathion/mal5000sims_r0med_2017-12-12.Rdata')
  load('Review_models/Savio/Malathion/mal5000sims_r0.75_2017-12-12.Rdata')
  load('Review_models/Savio/Malathion/mal5000sims_r0.25_2017-12-12.Rdata')
  load('Review_models/Savio/Malathion/mal5000sims_pars_2017-12-12.Rdata')

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

#Profenofos ################
  load('Review_models/Savio/Profenofos/pr5000sims_r0med_2017-12-12.Rdata')
  load('Review_models/Savio/Profenofos/pr5000sims_r0.75_2017-12-12.Rdata')
  load('Review_models/Savio/Profenofos/pr5000sims_r0.25_2017-12-12.Rdata')
  load('Review_models/Savio/Profenofos/pr5000sims_pars_2017-12-12.Rdata')
  
#Individual parameter effects    
par.prof.df = data.frame(conc = seq(0, 6.02, length.out = 500),
                        conc.eec = seq(0, 6.02, length.out = 500)/3.01,
                        chem = 'Profenofos')

par.prof.df = cbind(par.prof.df, prof.med.r0[,c(1:4)], prof.med.r0.25[,c(1:4)], prof.med.r0.75[,c(1:4)])

colnames(par.prof.df)[c(4:15)] = c(paste0('r0.',c('muPsat','piChash', 'piMhash', 'munhash')),
                                  paste0('r025.',c('muPsat','piChash', 'piMhash', 'munhash')),
                                  paste0('r075.',c('muPsat','piChash', 'piMhash', 'munhash')))

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

#Plot trajectories of herbicides and insecticides together############
#Herbicides ###############
#Parameter effects
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
  
#Insecticides ###############
#Parameter effects
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
  
#Combined herbicides and insecticides plot ############
  hrbdf.par$Class = 'Herbicides'
  insdf.par$Class = 'Insecticides'
 
compdf = rbind(hrbdf.par, insdf.par) 

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
  geom_ribbon(aes(x = conc.eec, ymin = r025, ymax = r075, fill = chem), alpha = 0.4) +
  geom_line(aes(x = conc.eec, y = r0, col = chem), size = 1)

windows(width = 30, height = 15)

gcomp

tiff('~/RemaisWork/Schisto/Agro_Review/Figures/r0_ConcRange/ins&herb_net_response.tiff', 
     width = 3300, height = 1650, res = 300)

gcomp

dev.off()
