#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

require(forestplot)
source('Review_models/fin/r0_of_q_fin.R')

date = '2018-02-19'  #updated functions across the board
date000 = '2017-12-11'  #updated functions across the board
date00 = '2017-11-27'  #runs as below, but with no implied 3rd data points for any d-r functions
date0 = '2017-11-22'  #Runs with 5000 sims and median/IQR values & Glyphosate EEC = 3700
date1 = '2017-11-11' #Runs with 5000 sims and median/IQR values
date2 = '2017-11-10' #Runs with 5000 sims
date3 = '2017-11-09' #Runs with 1000 sims

#Load the eec data frames and reorder for inclusion in forestplot ##########
#EEC data frames
load(paste('Review_models/r0_EECs/eec.p1.df', date, '.RData', sep = ''))
  eec.p1.df = eec.p1.df[-1,] #Get rid of simulation with atrazine d-r function that excludes 30ppb datapoint
  eec.p1.df = eec.p1.df[c(1,5,4,2,6,3),]
  eec.p1.df$relr0.025 = (eec.p1.df$r0 - 1.96*eec.p1.df$r0.sd) / r0.He()[3] * 100 - 100
  eec.p1.df$relr0.975 = (eec.p1.df$r0 + 1.96*eec.p1.df$r0.sd) / r0.He()[3] * 100 - 100
  
load(paste('Review_models/r0_EECs/eec.p2.df', date,'.RData', sep = ''))
  ind = sapply(eec.p2.df, is.factor)
  eec.p2.df[ind] = lapply(eec.p2.df[ind], as.character)
  eec.p2.df = eec.p2.df[order(eec.p2.df[,1], eec.p2.df[,4], eec.p2.df[,2]),]
  
  eec.p2.df$relr0.025 = (eec.p2.df$r0 - 1.96*eec.p2.df$r0.sd) / r0.He()[3] * 100 - 100
  eec.p2.df$relr0.975 = (eec.p2.df$r0 + 1.96*eec.p2.df$r0.sd) / r0.He()[3] * 100 - 100
  
  eec.p2.df = rbind(subset(eec.p2.df, study != 'Ragab & Shoukry 2006'), #reorder to put fertilizer estimates at the bottom
                    subset(eec.p2.df, study == 'Ragab & Shoukry 2006'))
  
load(paste('Review_models/r0_EECs/eec.p3.df', date,'.RData', sep = ''))
  ind = sapply(eec.p3.df, is.factor)
  eec.p3.df[ind] = lapply(eec.p3.df[ind], as.character)
  eec.p3.df = eec.p3.df[order(eec.p3.df[,1], eec.p3.df[,4], eec.p3.df[,2]),]

  eec.p3.df$relr0.025 = (eec.p3.df$r0 - 1.96*eec.p3.df$r0.sd) / r0.He()[3] * 100 - 100
  eec.p3.df$relr0.975 = (eec.p3.df$r0 + 1.96*eec.p3.df$r0.sd) / r0.He()[3] * 100 - 100
  
  eec.p3.df = subset(eec.p3.df, chem != 'Zinc')
  
load(paste('Review_models/r0_EECs/eec.p4.df', date,'.RData', sep = ''))
  ind = sapply(eec.p4.df, is.factor)
  eec.p4.df[ind] = lapply(eec.p4.df[ind], as.character)
  eec.p4.df = eec.p4.df[order(eec.p4.df[,1], eec.p4.df[,4], eec.p4.df[,2]),]

  eec.p4.df = subset(eec.p4.df, study != 'Meta')
  
  eec.p4.df$relr0.025 = (eec.p4.df$r0 - 1.96*eec.p4.df$r0.sd) / r0.He()[3] * 100 - 100
  eec.p4.df$relr0.975 = (eec.p4.df$r0 + 1.96*eec.p4.df$r0.sd) / r0.He()[3] * 100 - 100
  
  eec.p4.df = rbind(subset(eec.p4.df, study != 'Tchounwou et al 1991b'), #reorder to put fertilizer estimates at the bottom
                    subset(eec.p4.df, study == 'Tchounwou et al 1991b')[c(3,4,1,2),])
  
  eec.all = rbind(eec.p1.df, eec.p2.df, eec.p3.df, eec.p4.df)

#10% EEC data frames  
load(paste('Review_models/r0_EECs/eec0.1.p1.df', date,'.RData', sep = ''))
  eec0.1.p1.df = eec0.1.p1.df[-1,] #Get rid of simulation with atrazine d-r function that excludes 30ppb datapoint
  eec0.1.p1.df = eec0.1.p1.df[c(1,5,4,2,3,6),]
  
  eec0.1.p1.df$relr0.025 = (eec0.1.p1.df$r0 - 1.96*eec0.1.p1.df$r0.sd) / r0.He()[3] * 100 - 100
  eec0.1.p1.df$relr0.975 = (eec0.1.p1.df$r0 + 1.96*eec0.1.p1.df$r0.sd) / r0.He()[3] * 100 - 100
  
  eec0.1.p1.df[c(2:6),c(7:23)] = NA  #replace estimates with NAs for studies tested at one dose (no d-r fucntion)
  
  
load(paste('Review_models/r0_EECs/eec0.1.p2.df', date,'.RData', sep = ''))
  ind = sapply(eec0.1.p2.df, is.factor)
  eec0.1.p2.df[ind] = lapply(eec0.1.p2.df[ind], as.character)
  eec0.1.p2.df = eec0.1.p2.df[order(eec0.1.p2.df[,1], eec0.1.p2.df[,4], eec0.1.p2.df[,2]),]
  eec0.1.p2.df$relr0.025 = (eec0.1.p2.df$r0 - 1.96*eec0.1.p2.df$r0.sd) / r0.He()[3] * 100 - 100
  eec0.1.p2.df$relr0.975 = (eec0.1.p2.df$r0 + 1.96*eec0.1.p2.df$r0.sd) / r0.He()[3] * 100 - 100
  
    eec0.1.p2.df$relr0.lo[eec0.1.p2.df$relr0.lo < -100] = -100
  
    eec0.1.p2.df = rbind(subset(eec0.1.p2.df, study != 'Ragab & Shoukry 2006'), #reorder to put fertilizer estimates at the bottom
                         subset(eec0.1.p2.df, study == 'Ragab & Shoukry 2006'))
    
    eec0.1.p2.df[c(1, 7, 11, 15, 19, 25),c(7:23)] = NA   #replace estimates with NAs for studies tested at one dose (no d-r fucntion)
    
    
load(paste('Review_models/r0_EECs/eec0.1.p3.df', date,'.RData', sep = ''))
  ind = sapply(eec0.1.p3.df, is.factor)
  eec0.1.p3.df[ind] = lapply(eec0.1.p3.df[ind], as.character)
  eec0.1.p3.df = eec0.1.p3.df[order(eec0.1.p3.df[,1], eec0.1.p3.df[,4], eec0.1.p3.df[,2]),]
  
  eec0.1.p3.df$relr0.025 = (eec0.1.p3.df$r0 - 1.96*eec0.1.p3.df$r0.sd) / r0.He()[3] * 100 - 100
  eec0.1.p3.df$relr0.975 = (eec0.1.p3.df$r0 + 1.96*eec0.1.p3.df$r0.sd) / r0.He()[3] * 100 - 100
  
  eec0.1.p3.df = subset(eec0.1.p3.df, chem != 'Zinc')
  
load(paste('Review_models/r0_EECs/eec0.1.p4.df', date,'.RData', sep = ''))
  ind = sapply(eec0.1.p4.df, is.factor)
  eec0.1.p4.df[ind] = lapply(eec0.1.p4.df[ind], as.character)
  eec0.1.p4.df = eec0.1.p4.df[order(eec0.1.p4.df[,1], eec0.1.p4.df[,4], eec0.1.p4.df[,2]),]

  eec0.1.p4.df = subset(eec0.1.p4.df, study != 'Meta')
  
  eec0.1.p4.df$relr0.025 = (eec0.1.p4.df$r0 - 1.96*eec0.1.p4.df$r0.sd) / r0.He()[3] * 100 - 100
  eec0.1.p4.df$relr0.975 = (eec0.1.p4.df$r0 + 1.96*eec0.1.p4.df$r0.sd) / r0.He()[3] * 100 - 100
  
  eec0.1.p4.df = rbind(subset(eec0.1.p4.df, study != 'Tchounwou et al 1991b'), #reorder to put fertilizer estimates at the bottom
                       subset(eec0.1.p4.df, study == 'Tchounwou et al 1991b')[c(3,4,1,2),])
  
dev.off()  

#Text list for forest plot ########
tabtext =list(list('Study', #Studies ##############
                   
                   'Pathway 1: Bottom-up ecological effects',
                   expression(bold('    Baxter et al, 2011'^'2')),
                   expression('    Halstead et al 2017'^'3'),
                   expression('    Rohr et al, 2008'^'3'),
                   expression('    Johnson et al, 2007'^'3'),
                   expression('    Halstead et al 2017'^'3'),
                   expression('    Johnson et al, 2007'^'3'),
                   
                   ' ',
                   
                   'Pathway 2: Direct effects on snails',
                   expression('    Bakry et al 2012'^'3'),  #only provide two data points
                   '    Bakry et al 2012', 
                   expression(bold('    Omran & Salama 2013')), 
                   expression(bold('    Abdel-Ghaffar et al 2016')),
                   expression(bold('    Abdel-Ghaffar et al 2016')),
                   '    Tantawy 2002',
                   expression('    Hasheesh & Mohamed 2011'^'3'),
                   expression(bold('    Ibrahim et al 1992')),
                   expression(bold('    Hasheesh & Mohamed 2011')),
                   '    Ibrahim et al 1992',
                   expression('    Bakry et al 2011'^'3'),
                   '    Bakry et al 2011',
                   '    Tantawy 2002',
                   expression(bold('    Abdel-Ghaffar et al 2016')),
                   expression('    Bakry et al 2012'^'3'), #only provide two data points
                   expression(bold('    Abdel-Ghaffar et al 2016')),
                   '    Bakry et al 2012', 
                   '    Omran & Salama 2013', 
                   expression('    Bakry et al 2011'^'3'),
                   expression(bold('    Tchounwou et al 1991')), 
                   '    Bakry et al 2011', 
                   expression(bold('    Tchounwou et al 1991')), 
                   '    Abdel-Ghaffar et al 2016',
                   '    Abdel-Ghaffar et al 2016',
                   expression('    Hasheesh & Mohamed 2011'^'3'),
                   expression(bold('    Hasheesh & Mohamed 2011')), 
                   '    Ragab & Shoukry 2006',
                   '    Ragab & Shoukry 2006',
                   '    Ragab & Shoukry 2006',
                   ' ',
                   
                   'Pathway 3: Top-down ecological effects',
                   expression(bold('    Halstead et al 2015')), 
                   '    Satapornvanit et al 2009', 
                   '    Satapornvanit et al 2009', 
                   '    Satapornvanit et al 2009', 
                   '    Halstead et al 2015',
                   '    Halstead et al 2015',
                   expression(bold('    Halstead et al 2015')),
                   '    Halstead et al 2015',
                   expression(bold('    Satapornvanit et al 2009')),
                   '    Halstead et al 2015',
                   #'    Satapornvanit et al 2009', #zinc not considered an agrochemical
                   #'    Satapornvanit et al 2009', #zinc not considered an agrochemical
                   ' ',
                   
                   'Pathway 4: Direct effects on schistosomes',
                   '    Griggs et al 2008', 
                   '    Koprivnikar et al 2006', 
                   #'    Meta',
                   expression(bold('    Rohr et al 2008')), 
                   '    Tantawy 2002', 
                   '    Tantawy 2002', 
                   expression(bold('    Abdel-Ghaffar et al 2016')),
                   expression(bold('    Abdel-Ghaffar et al 2016')),
                   expression(bold('    Hasheesh & Mohamed 2011')),
                   expression(bold('    Hasheesh & Mohamed 2011')),
                   '    Tantawy 2002',
                   '    Tantawy 2002', 
                   expression(bold('    Abdel-Ghaffar et al 2016')),
                   expression(bold('    Abdel-Ghaffar et al 2016')),
                   expression(bold('    Tchounwou et al 1992')),
                   expression(bold('    Tchounwou et al 1991a')),
                   '    Abdel-Ghaffar et al 2016',
                   '    Abdel-Ghaffar et al 2016',
                   expression(bold('    Hasheesh & Mohamed 2011')),
                   expression(bold('    Hasheesh & Mohamed 2011')),
                   '    Tchounwou et al 1991b',
                   '    Tchounwou et al 1991b',
                   '    Tchounwou et al 1991b',
                   '    Tchounwou et al 1991b'),
              
              list(expression('Agrochemical'^'1'), #Chemicals ##############
                   
                   ' ',
                   'Atrazine',
                   'Atrazine',
                   'Atrazine',
                   expression(paste('NH'[4],'NO'[3],' & H'[3],'PO'[4])),
                   expression(paste('NaNO'[3], ' & NaH'[2], 'PO'[4])),
                   expression(paste('NH'[4],'NO'[3],' & H'[3],'PO'[4])),
                   ' ',
                   
                   ' ',
                   'Atrazine', 
                   'Atrazine', 
                   'Atrazine', 
                   'Butralin', 
                   'Butralin', 
                   'Butralin', 
                   'Chlorpyrifos', 
                   'Chlorpyrifos', 
                   'Chlorpyrifos', 
                   'Chlorpyrifos',
                   'Deltamethrin', 
                   'Deltamethrin', 
                   'Fluazifop-p-butyl',
                   'Glyphosate', 
                   'Glyphosate',
                   'Glyphosate', 
                   'Glyphosate', 
                   'Glyphosate', 
                   'Malathion', 
                   'Malathion', 
                   'Malathion', 
                   'Malathion',
                   'Pendimethalin', 
                   'Pendimethalin', 
                   'Profenofos',
                   'Profenofos', 
                   expression(paste('CO(NH'[2],')'[2])), 
                   expression(paste('NH'[4],'NO'[3])),
                   expression(paste('K'[2],'SO'[4])), 
                   ' ',
                   
                   ' ',
                   'Chlorpyrifos', 
                   'Chlorpyrifos', 
                   'Chlorpyrifos', 
                   'Dimethoate', 
                   'Esfenvalerate', 
                   expression(paste(lambda,'-Cyhalothrin', sep = '')),
                   'Malathion', 
                   'Permethrin', 
                   'Profenofos', 
                   'Terbufos', 
                   #'Zinc', 
                   #'Zinc',
                   ' ',
                   
                   ' ',
                   'Atrazine',
                   'Atrazine',
                   #'Atrazine',
                   'Atrazine',
                   'Butachlor',
                   'Butachlor', 
                   'Butralin', 
                   'Butralin', 
                   'Chlorpyrifos', 
                   'Chlorpyrifos', 
                   'Fluazifop-p-butyl', 
                   'Fluazifop-p-butyl', 
                   'Glyphosate', 
                   'Glyphosate', 
                   'Malathion', 
                   'Malathion', 
                   'Pendimethalin',
                   'Pendimethalin', 
                   'Profenofos', 
                   'Profenofos', 
                   expression(paste('CO(NH'[2],')'[2])), 
                   expression(paste('CO(NH'[2],')'[2])), 
                   expression(paste('(NH'[4],')'[2],'SO'[4])), 
                   expression(paste('(NH'[4],')'[2],'SO'[4]))), 
                   
              list('Species', #Species ##############
                   ' ', #Pathway 1
                   expression(italic('Physella spp')),
                   expression(italic('Bulinus truncatus')), 
                   expression(italic('Planorbella trivolvis')),
                   expression(italic('Planorbella trivolvis')),
                   expression(italic('Bulinus truncatus')), 
                   expression(italic('Planorbella trivolvis')),
                   ' ', ' ',  #Pathway 2
                   
                   expression(italic('Biomphalaria alexandrina')),
                   expression(italic('Biomphalaria alexandrina')), 
                   expression(italic('Biomphalaria alexandrina')),
                   expression(italic('Biomphalaria alexandrina')), 
                   expression(italic('Biomphalaria alexandrina')),
                   expression(italic('Biomphalaria alexandrina')), 
                   expression(italic('Bulinus truncatus')), 
                   expression(italic('Biomphalaria alexandrina')),
                   expression(italic('Bulinus truncatus')),
                   expression(italic('Biomphalaria alexandrina')), 
                   expression(italic('Helisoma duryi')),
                   expression(italic('Helisoma duryi')),
                   expression(italic('Biomphalaria alexandrina')),
                   expression(italic('Biomphalaria alexandrina')), 
                   expression(italic('Biomphalaria alexandrina')),
                   expression(italic('Biomphalaria alexandrina')), 
                   expression(italic('Biomphalaria alexandrina')),
                   expression(italic('Biomphalaria alexandrina')),
                   expression(italic('Helisoma duryi')),
                   expression(italic('Bulinus havenensis')), 
                   expression(italic('Helisoma duryi')),
                   expression(italic('Bulinus havenensis')),
                   expression(italic('Biomphalaria alexandrina')), 
                   expression(italic('Biomphalaria alexandrina')),
                   expression(italic('Bulinus truncatus')), 
                   expression(italic('Bulinus truncatus')),
                   expression(italic('Biomphalaria alexandrina')), 
                   expression(italic('Biomphalaria alexandrina')), 
                   expression(italic('Biomphalaria alexandrina')), 
                   ' ', ' ',   #Pathway 3
                
                   expression(italic('Procambarus clarkii')),
                   expression(italic('Macrobrachium rosenbergii')),
                   expression(italic('Macrobrachium rosenbergii')),
                   expression(italic('Macrobrachium rosenbergii')),
                   expression(italic('Procambarus clarkii')),
                   expression(italic('Procambarus clarkii')),
                   expression(italic('Procambarus clarkii')),
                   expression(italic('Procambarus clarkii')),
                   expression(italic('Macrobrachium rosenbergii')),
                   expression(italic('Procambarus clarkii')),
                   #expression(italic('Macrobrachium rosenbergii')),
                   #expression(italic('Macrobrachium rosenbergii')),
                   ' ', ' ',   #Pathway 4
                   
                   expression(italic('Echinistoma trivolvis')),
                   expression(italic('Echinistoma trivolvis')),
                   #expression(italic('Echinistoma trivolvis')),
                   expression(italic('Echinistoma trivolvis')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma haemotobium')),
                   expression(italic('Schistosoma haemotobium')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma haemotobium')),
                   expression(italic('Schistosoma haemotobium')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni')),
                   expression(italic('Schistosoma mansoni'))),
               
              list('Parameter', #Parameters ##############
                   ' ', #Pathway 1
                   expression(paste(Phi[N])),
                   expression(paste(Phi[N])),
                   expression(paste(Phi[N])),
                   expression(paste(Phi[N])),
                   expression(paste(Phi[N])),
                   expression(paste(theta)),
                   
                   ' ', ' ', #Pathway 2
                   expression(paste('f'[N])),
                   expression(paste(mu[N])),
                   expression(paste(mu[N])),
                   expression(paste('f'[N])),
                   expression(paste(mu[N])), 
                   expression(paste(mu[N])),
                   expression(paste('f'[N])),
                   expression(paste('f'[N])),
                   expression(paste(mu[N])),
                   expression(paste(mu[N])),
                   expression(paste('f'[N])),
                   expression(paste(mu[N])), 
                   expression(paste(mu[N])),
                   expression(paste('f'[N])),
                   expression(paste('f'[N])),
                   expression(paste(mu[N])),
                   expression(paste(mu[N])),
                   expression(paste(mu[N])),
                   expression(paste('f'[N])),
                   expression(paste('f'[N])),
                   expression(paste(mu[N])),
                   expression(paste(mu[N])),
                   expression(paste('f'[N])),
                   expression(paste(mu[N])),
                   expression(paste('f'[N])),
                   expression(paste(mu[N])),
                   expression(paste(mu[N])),
                   expression(paste(mu[N])), 
                   expression(paste(mu[N])),
                   
                   ' ', ' ', #Pathway 3
                   expression(paste(mu[P])),
                   expression(paste(mu[P])),
                   expression(paste(psi[P])),
                   expression(paste(mu[P])),
                   expression(paste(mu[P])),
                   expression(paste(mu[P])),
                   expression(paste(mu[P])),
                   expression(paste(mu[P])),
                   expression(paste(mu[P])),
                   expression(paste(mu[P])),
                   #expression(paste(mu[P])),
                   #expression(paste(psi[P])),
                   
                   ' ', ' ', #Pathway 4
                   expression(paste(pi[C])),
                   expression(paste(pi[C])),
                   expression(paste(pi[C])), 
                   #expression(paste(pi[C])),
                   expression(paste(pi[C])),
                   expression(paste(pi[M])), 
                   expression(paste(pi[C])),
                   expression(paste(pi[M])),
                   expression(paste(pi[C])), 
                   expression(paste(pi[M])),
                   expression(paste(pi[C])),
                   expression(paste(pi[M])),
                   expression(paste(pi[C])),
                   expression(paste(pi[M])),
                   expression(paste(pi[C])),
                   expression(paste(pi[M])),
                   expression(paste(pi[C])),
                   expression(paste(pi[M])),
                   expression(paste(pi[C])),
                   expression(paste(pi[M])),
                   expression(paste(pi[M])), 
                   expression(italic('v')), 
                   expression(paste(pi[M])),
                   expression(italic('v')))
               
)

#Forest plots #########
#With relative change in median r0 and IQR at both EEC and 10% EEC values #########
tiff(paste('~/RemaisWork/Schisto/Agro_Review/Figures/EEC_forest/Pathways_forest_medIQR', date, '.tiff', sep = ''),
     width = 4100*0.9, height = 5893*0.9, res = 300)
  forestplot(labeltext = tabtext, 
             legend_args = fpLegend(pos = list(x = 0.25, y = 0.99)),
             legend = c('EEC', '10% EEC'),
             fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
             line.margin = 0.1,
             mean = cbind(c(NA,NA, eec.p1.df$relr0.med,   
                            NA,NA, eec.p2.df$relr0.med,   
                            NA,NA, eec.p3.df$relr0.med,   
                            NA,NA, eec.p4.df$relr0.med), 
                          c(NA,NA, eec0.1.p1.df$relr0.med,   
                            NA,NA, eec0.1.p2.df$relr0.med,   
                            NA,NA, eec0.1.p3.df$relr0.med,   
                            NA,NA, eec0.1.p4.df$relr0.med)),
             lower =cbind(c(NA,NA, eec.p1.df$relr0.25,
                            NA,NA, eec.p2.df$relr0.25,
                            NA,NA, eec.p3.df$relr0.25,
                            NA,NA, eec.p4.df$relr0.25),
                          c(NA,NA, eec0.1.p1.df$relr0.25,
                            NA,NA, eec0.1.p2.df$relr0.25,
                            NA,NA, eec0.1.p3.df$relr0.25,
                            NA,NA, eec0.1.p4.df$relr0.25)),
             upper =cbind(c(NA,NA, eec.p1.df$relr0.75,
                            NA,NA, eec.p2.df$relr0.75,
                            NA,NA, eec.p3.df$relr0.75,
                            NA,NA, eec.p4.df$relr0.75),
                          c(NA,NA, eec0.1.p1.df$relr0.75,
                            NA,NA, eec0.1.p2.df$relr0.75,
                            NA,NA, eec0.1.p3.df$relr0.75,
                            NA,NA, eec0.1.p4.df$relr0.75)),
             new_page = TRUE,
             is.summary=c(TRUE,TRUE,
                          rep(FALSE,nrow(eec.p1.df)+1),
                          TRUE,
                          rep(FALSE,nrow(eec.p2.df)+1),
                          TRUE,
                          rep(FALSE,nrow(eec.p3.df)+1),
                          TRUE,
                          rep(FALSE,nrow(eec.p4.df)+1)),
             hrzl_lines = list('2'= gpar(lwd = 1, col="black", columns = 1:4),
                               '3'= gpar(col = 'grey'),
                               '11'= gpar(col = 'grey'),
                               '42'= gpar(col = 'grey'),
                               '54'= gpar(col = 'grey')),
             txt_gp = fpTxtGp(xlab = gpar(cex = 1.2),
                              ticks = gpar(cex = 1.1)),
             vertices = TRUE,
             boxsize = 0.3,
             col = fpColors(box = c('black', 'blue'), lines = 'black'),
             clip=c(-Inf,Inf),
             xlab = expression(paste(Delta, R['0'], ' (%)')),
             xticks = c(-100, -50,0,50,100))
  dev.off()
  
#With raw change in r0 at both EEC and 10% EEC values #########
  tiff(paste('~/RemaisWork/Schisto/Agro_Review/Figures/EEC_forest/Pathways_forest_rawR0', date, '.tiff', sep = ''),
       width = 4100*0.9, height = 5893*0.9, res = 300)
  forestplot(labeltext = tabtext, 
             legend_args = fpLegend(pos = list(x = 0.25, y = 0.99)),
             legend = c('EEC', '10% EEC'),
             fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
             line.margin = 0.1,
             mean = cbind(c(NA,NA, eec.p1.df$r0.med,   
                            NA,NA, eec.p2.df$r0.med,   
                            NA,NA, eec.p3.df$r0.med,   
                            NA,NA, eec.p4.df$r0.med), 
                          c(NA,NA, eec0.1.p1.df$r0.med,   
                            NA,NA, eec0.1.p2.df$r0.med,   
                            NA,NA, eec0.1.p3.df$r0.med,   
                            NA,NA, eec0.1.p4.df$r0.med)),
             lower =cbind(c(NA,NA, eec.p1.df$r0.25,
                            NA,NA, eec.p2.df$r0.25,
                            NA,NA, eec.p3.df$r0.25,
                            NA,NA, eec.p4.df$r0.25),
                          c(NA,NA, eec0.1.p1.df$r0.25,
                            NA,NA, eec0.1.p2.df$r0.25,
                            NA,NA, eec0.1.p3.df$r0.25,
                            NA,NA, eec0.1.p4.df$r0.25)),
             upper =cbind(c(NA,NA, eec.p1.df$r0.75,
                            NA,NA, eec.p2.df$r0.75,
                            NA,NA, eec.p3.df$r0.75,
                            NA,NA, eec.p4.df$r0.75),
                          c(NA,NA, eec0.1.p1.df$r0.75,
                            NA,NA, eec0.1.p2.df$r0.75,
                            NA,NA, eec0.1.p3.df$r0.75,
                            NA,NA, eec0.1.p4.df$r0.75)),
             new_page = TRUE,
             is.summary=c(TRUE,TRUE,
                          rep(FALSE,nrow(eec.p1.df)+1),
                          TRUE,
                          rep(FALSE,nrow(eec.p2.df)+1),
                          TRUE,
                          rep(FALSE,nrow(eec.p3.df)+1),
                          TRUE,
                          rep(FALSE,nrow(eec.p4.df)+1)),
             hrzl_lines = list('2'= gpar(lwd = 1, col="black", columns = 1:4),
                               '3'= gpar(col = 'grey'),
                               '11'= gpar(col = 'grey'),
                               '42'= gpar(col = 'grey'),
                               '54'= gpar(col = 'grey')),
             txt_gp = fpTxtGp(xlab = gpar(cex = 1.2),
                              ticks = gpar(cex = 1.1)),
             vertices = TRUE,
             boxsize = 0.3,
             col = fpColors(box = c('black', 'blue'), lines = 'black'),
             clip=c(-Inf,Inf),
             xlab = expression(paste(R['0'])),
             zero = r0.He()[3],
             xticks = c(0, 1, 2, 3, 4, 5))
  dev.off()
  
#With relative change in r0 +/- SD at both EEC and 10% EEC values #########
dev.off()
tiff(paste('~/RemaisWork/Schisto/Agro_Review/Figures/EEC_forest/Pathways_forest_', date, '.tiff', sep = ''),
     width = 4200, height = 6037, res = 300)
forestplot(tabtext, 
           legend_args = fpLegend(pos = list(x = 0.25, y = 0.99)),
           legend = c('EEC', '10% EEC'),
           fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
           line.margin = 0.1,
           mean = cbind(c(NA,NA, eec.p1.df$relr0,   
                          NA,NA, eec.p2.df$relr0,   
                          NA,NA, eec.p3.df$relr0,   
                          NA,NA, eec.p4.df$relr0), 
                        c(NA,NA, eec0.1.p1.df$relr0,   
                          NA,NA, eec0.1.p2.df$relr0,   
                          NA,NA, eec0.1.p3.df$relr0,   
                          NA,NA, eec0.1.p4.df$relr0)),
           lower =cbind(c(NA,NA, eec.p1.df$relr0.lo,
                          NA,NA, eec.p2.df$relr0.lo,
                          NA,NA, eec.p3.df$relr0.lo,
                          NA,NA, eec.p4.df$relr0.lo),
                        c(NA,NA, eec0.1.p1.df$relr0.lo,
                          NA,NA, eec0.1.p2.df$relr0.lo,
                          NA,NA, eec0.1.p3.df$relr0.lo,
                          NA,NA, eec0.1.p4.df$relr0.lo)),
           upper =cbind(c(NA,NA, eec.p1.df$relr0.up,
                          NA,NA, eec.p2.df$relr0.up,
                          NA,NA, eec.p3.df$relr0.up,
                          NA,NA, eec.p4.df$relr0.up),
                        c(NA,NA, eec0.1.p1.df$relr0.up,
                          NA,NA, eec0.1.p2.df$relr0.up,
                          NA,NA, eec0.1.p3.df$relr0.up,
                          NA,NA, eec0.1.p4.df$relr0.up)),
           new_page = TRUE,
           is.summary=c(TRUE,TRUE,
                        rep(FALSE,nrow(eec.p1.df)+1),
                        TRUE,
                        rep(FALSE,nrow(eec.p2.df)+1),
                        TRUE,
                        rep(FALSE,nrow(eec.p3.df)+1),
                        TRUE,
                        rep(FALSE,nrow(eec.p4.df)+1)),
           hrzl_lines = list('2'= gpar(lwd = 1, col="black", columns = 1:4),
                             '3'= gpar(col = 'grey'),
                             '9'= gpar(col = 'grey'),
                             '40'= gpar(col = 'grey'),
                             '54'= gpar(col = 'grey')),
           txt_gp = fpTxtGp(xlab = gpar(cex = 1.2),
                            ticks = gpar(cex = 1.1)),
           vertices = TRUE,
           boxsize = 0.3,
           col = fpColors(box = c('black', 'blue'), lines = 'black'),
           clip=c(-Inf,Inf),
           xlab = expression(paste(Delta, R['0'], ' (%)')),
           xticks = c(-100, -50,0,50,100,200,300,400))
dev.off()

#With relative change in r0 +/- 95% CI at both EEC and 10% EEC values #########
dev.off()
tiff(paste('~/RemaisWork/Schisto/Agro_Review/Figures/EEC_forest/Pathways_forest_95CIs_', date, '.tiff', sep = ''),
     width = 4200, height = 6037, res = 300)
forestplot(tabtext, 
           legend_args = fpLegend(pos = list(x = 0.25, y = 0.99)),
           legend = c('EEC', '10% EEC'),
           fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
           line.margin = 0.1,
           mean = cbind(c(NA,NA, eec.p1.df$relr0,   
                          NA,NA, eec.p2.df$relr0,   
                          NA,NA, eec.p3.df$relr0,   
                          NA,NA, eec.p4.df$relr0), 
                        c(NA,NA, eec0.1.p1.df$relr0,   
                          NA,NA, eec0.1.p2.df$relr0,   
                          NA,NA, eec0.1.p3.df$relr0,   
                          NA,NA, eec0.1.p4.df$relr0)),
           lower =cbind(c(NA,NA, eec.p1.df$relr0.025,
                          NA,NA, eec.p2.df$relr0.025,
                          NA,NA, eec.p3.df$relr0.025,
                          NA,NA, eec.p4.df$relr0.025),
                        c(NA,NA, eec0.1.p1.df$relr0.025,
                          NA,NA, eec0.1.p2.df$relr0.025,
                          NA,NA, eec0.1.p3.df$relr0.025,
                          NA,NA, eec0.1.p4.df$relr0.025)),
           upper =cbind(c(NA,NA, eec.p1.df$relr0.975,
                          NA,NA, eec.p2.df$relr0.975,
                          NA,NA, eec.p3.df$relr0.975,
                          NA,NA, eec.p4.df$relr0.975),
                        c(NA,NA, eec0.1.p1.df$relr0.975,
                          NA,NA, eec0.1.p2.df$relr0.975,
                          NA,NA, eec0.1.p3.df$relr0.975,
                          NA,NA, eec0.1.p4.df$relr0.975)),
           new_page = TRUE,
           is.summary=c(TRUE,TRUE,
                        rep(FALSE,nrow(eec.p1.df)+1),
                        TRUE,
                        rep(FALSE,nrow(eec.p2.df)+1),
                        TRUE,
                        rep(FALSE,nrow(eec.p3.df)+1),
                        TRUE,
                        rep(FALSE,nrow(eec.p4.df)+1)),
           hrzl_lines = list('2'= gpar(lwd = 1, col="black", columns = 1:4),
                             '3'= gpar(col = 'grey'),
                             '9'= gpar(col = 'grey'),
                             '40'= gpar(col = 'grey'),
                             '54'= gpar(col = 'grey')),
           txt_gp = fpTxtGp(xlab = gpar(cex = 1.2),
                            ticks = gpar(cex = 1.1)),
           vertices = TRUE,
           boxsize = 0.3,
           col = fpColors(box = c('black', 'blue'), lines = 'black'),
           clip=c(-Inf,Inf),
           xlab = expression(paste(Delta, R['0'], ' (%)')),
           xticks = c(-100, -50,0,50,100,200,300,400))

dev.off()


#With raw change in r0 at EEC values ######
windows(width = 36, height = 28)
forestplot(tabtext, 
           vals.raw,
           new_page = TRUE,
           is.summary=c(TRUE,TRUE,
                        rep(FALSE,nrow(eec.p1.df)+1),
                        TRUE,
                        rep(FALSE,nrow(eec.p2.df)+1),
                        TRUE,
                        rep(FALSE,nrow(eec.p3.df)+1),
                        TRUE,
                        rep(FALSE,nrow(eec.p4.df)+1)),
           hrzl_lines = list('2'= gpar(lwd = 1, col="black", columns = 1:3),
                             '3'= gpar(col = 'grey'),
                             '9'= gpar(col = 'grey'),
                             '40'= gpar(col = 'grey'),
                             '54'= gpar(col = 'grey')),
           txt_gp = fpTxtGp(xlab = gpar(cex = 1.2),
                            ticks = gpar(cex = 1.1)),
           vertices = TRUE,
           boxsize = 0.25,
           col = fpColors(box = c('black', 'blue'), lines = 'black'),
           clip=c(-Inf,Inf),
           xlab = expression(paste(Delta, R['0'])))

#With raw change in r0 at both EEC and 50% EEC values #########
windows(width = 24, height = 15)
forestplot(tabtext, 
           fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
           line.margin = 0.1,
           mean = cbind(c(NA,NA, eec.p1.df$deltar0,   
                          NA,NA, eec.p2.df$deltar0,   
                          NA,NA, eec.p3.df$deltar0,   
                          NA,NA, eec.p4.df$deltar0), 
                        c(NA,NA, atr.0.5eec.df$deltar0,   
                          NA,NA, ch.0.5eec.df$deltar0,   
                          NA,NA, gly.0.5eec.df$deltar0,   
                          NA,NA, mal.0.5eec.df$deltar0)),
           lower =cbind(c(NA,NA, eec.p1.df$deltar0.lo,
                          NA,NA, eec.p2.df$deltar0.lo,
                          NA,NA, eec.p3.df$deltar0.lo,
                          NA,NA, eec.p4.df$deltar0.lo),
                        c(NA,NA, atr.0.5eec.df$deltar0.lo,
                          NA,NA, ch.0.5eec.df$deltar0.lo,
                          NA,NA, gly.0.5eec.df$deltar0.lo,
                          NA,NA, mal.0.5eec.df$deltar0.lo)),
           upper =cbind(c(NA,NA, eec.p1.df$deltar0.up,
                          NA,NA, eec.p2.df$deltar0.up,
                          NA,NA, eec.p3.df$deltar0.up,
                          NA,NA, eec.p4.df$deltar0.up),
                        c(NA,NA, atr.0.5eec.df$deltar0.up,
                          NA,NA, ch.0.5eec.df$deltar0.up,
                          NA,NA, gly.0.5eec.df$deltar0.up,
                          NA,NA, mal.0.5eec.df$deltar0.up)),
           new_page = TRUE,
           is.summary=c(TRUE,TRUE,
                        rep(FALSE,nrow(eec.p1.df)+1),
                        TRUE,
                        rep(FALSE,nrow(eec.p2.df)+1),
                        TRUE,
                        rep(FALSE,nrow(eec.p3.df)+1),
                        TRUE,
                        rep(FALSE,nrow(eec.p4.df)+1)),
           hrzl_lines = list('2'= gpar(lwd = 1, col="black", columns = 1:3),
                             '3'= gpar(col = 'grey'),
                             '12'= gpar(col = 'grey'),
                             '24'= gpar(col = 'grey'),
                             '34'= gpar(col = 'grey')),
           txt_gp = fpTxtGp(xlab = gpar(cex = 1.2),
                            ticks = gpar(cex = 1.1)),
           vertices = TRUE,
           boxsize = 0.25,
           col = fpColors(box = c('black', 'darkred'), lines = 'black'),
           clip=c(-Inf,Inf),
           xlab = expression(paste(Delta, R['0'])))

