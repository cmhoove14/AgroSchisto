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
source('Review_models/r0_of_q.R')

#Load the eec data frames ##########
load('Review_models/r0_EECs/atr.eec.df.RData')
  atr.eec.df = subset(atr.eec.df, study != 'piC meta')
  atr.eec.df = atr.eec.df[c(order(as.character(atr.eec.df$Parameter)[-nrow(atr.eec.df)]),nrow(atr.eec.df)),]
  #Get raw change in R0
    atr.eec.df$deltar0 = atr.eec.df$r0 - r0.fix()[3]
    atr.eec.df$deltar0.up = atr.eec.df$r0.up - r0.fix()[3]
    atr.eec.df$deltar0.lo = atr.eec.df$r0.lo - r0.fix()[3]
  #Get relevant change in R0
    atr.eec.df$relr0 = atr.eec.df$r0 / r0.fix()[3]
    atr.eec.df$relr0.up = atr.eec.df$r0.up / r0.fix()[3]
    atr.eec.df$relr0.lo = atr.eec.df$r0.lo / r0.fix()[3]
  
load('Review_models/r0_EECs/ch.eec.df.RData')
  ch.eec.df = ch.eec.df[c(order(as.character(ch.eec.df$Parameter)[-nrow(ch.eec.df)]),nrow(ch.eec.df)),]
  #Get raw change in R0
    ch.eec.df$deltar0 = ch.eec.df$r0 - r0.fix()[3]
    ch.eec.df$deltar0.up = ch.eec.df$r0.up - r0.fix()[3]
    ch.eec.df$deltar0.lo = ch.eec.df$r0.lo - r0.fix()[3]
  #Get relevant change in R0
    ch.eec.df$relr0 = ch.eec.df$r0 / r0.fix()[3]
    ch.eec.df$relr0.up = ch.eec.df$r0.up / r0.fix()[3]
    ch.eec.df$relr0.lo = ch.eec.df$r0.lo / r0.fix()[3]
  
load('Review_models/r0_EECs/gly.eec.df.RData')
  gly.eec.df = subset(gly.eec.df, study != 'Combined2')
  gly.eec.df = gly.eec.df[c(order(as.character(gly.eec.df$Parameter)[-nrow(gly.eec.df)]),nrow(gly.eec.df)),]
  #Get raw change in R0
    gly.eec.df$deltar0 = gly.eec.df$r0 - r0.fix()[3]
    gly.eec.df$deltar0.up = gly.eec.df$r0.up - r0.fix()[3]
    gly.eec.df$deltar0.lo = gly.eec.df$r0.lo - r0.fix()[3]
  #Get relative change in R0
    gly.eec.df$relr0 = gly.eec.df$r0 / r0.fix()[3]
    gly.eec.df$relr0.up = gly.eec.df$r0.up / r0.fix()[3]
    gly.eec.df$relr0.lo = gly.eec.df$r0.lo / r0.fix()[3]
  
load('Review_models/r0_EECs/mal.eec.df.RData')
  mal.eec.df = mal.eec.df[c(order(as.character(mal.eec.df$Parameter)[-nrow(mal.eec.df)]),nrow(mal.eec.df)),]
  #Get raw change in r0
    mal.eec.df$deltar0 = mal.eec.df$r0 - r0.fix()[3]
    mal.eec.df$deltar0.up = mal.eec.df$r0.up - r0.fix()[3]
    mal.eec.df$deltar0.lo = mal.eec.df$r0.lo - r0.fix()[3]
  #Get relative change in r0  
    mal.eec.df$relr0 = mal.eec.df$r0 / r0.fix()[3]
    mal.eec.df$relr0.up = mal.eec.df$r0.up / r0.fix()[3]
    mal.eec.df$relr0.lo = mal.eec.df$r0.lo / r0.fix()[3]

#Data frames of values for forest plot ##############  
#Df for raw values plot     
vals.raw = data.frame(coef = c(NA,NA, atr.eec.df$deltar0,   NA, NA,ch.eec.df$deltar0,   NA, NA, gly.eec.df$deltar0,   NA, NA,mal.eec.df$deltar0),
                      low =  c(NA,NA, atr.eec.df$deltar0.lo,NA, NA,ch.eec.df$deltar0.lo,NA, NA, gly.eec.df$deltar0.lo,NA, NA,mal.eec.df$deltar0.lo),
                      high = c(NA,NA, atr.eec.df$deltar0.up,NA, NA,ch.eec.df$deltar0.up,NA, NA, gly.eec.df$deltar0.up,NA, NA,mal.eec.df$deltar0.up))

#Df for relative values
vals.rel = data.frame(coef = c(NA,NA, atr.eec.df$relr0,   NA, NA,ch.eec.df$relr0,   NA, NA, gly.eec.df$relr0,   NA, NA,mal.eec.df$relr0),
                      low =  c(NA,NA, atr.eec.df$relr0.lo,NA, NA,ch.eec.df$relr0.lo,NA, NA, gly.eec.df$relr0.lo,NA, NA,mal.eec.df$relr0.lo),
                      high = c(NA,NA, atr.eec.df$relr0.up,NA, NA,ch.eec.df$relr0.up,NA, NA, gly.eec.df$relr0.up,NA, NA,mal.eec.df$relr0.up))

#Text list for forest plot ########
tabtext =list(list('Study', 
                   
                   'Atrazine',
                   '    Bakry et al, 2012',
                   expression(bold('    Omran & Salama, 2013')),
                   expression(bold('    Rohr et al, 2012')),
                   '    Griggs & Belden, 2008',
                   '    Koprivnikar et al, 2007',
                   expression(bold('    Rohr et al, 2008')),
                   expression(bold('  Combined')), ' ',
                   
                   'Chlorpyrifos',
                                      '    Hasheesh & Mohamed, 2011',
                   expression(bold('    Ibrahim et al, 1992')),
                   expression(bold('    Hasheesh & Mohamed, 2011')),
                   '    Ibrahim et al, 1992', 
                   expression(bold('    Halstead et al, 2015')),
                   '    Satapornvanit et al, 2009',
                   expression(bold('    Hasheesh & Mohamed, 2011')),
                   expression(bold('    Hasheesh & Mohamed, 2011')),
                   '    Satapornvanit et al, 2009',
                   expression(bold('  Combined')), ' ',
                   
                   'Glyphosate',
                   expression(bold('    Abdel-Ghaffar et al, 2016')),
                   '    Bakry et al, 2012',
                   expression(bold('    Abdel-Ghaffar et al, 2016')),
                   '    Bakry et al, 2012',
                   '    Omran & Salama 2013',
                   expression(bold('    Abdel-Ghaffar et al, 2016')),
                   expression(bold('    Abdel-Ghaffar et al, 2016')), 
                   expression(bold('  Combined')), ' ',
                   
                   'Malathion',
                   '    Bakry et al, 2011',
                   expression(bold('    Tchounwou et al, 1991b')),
                   '    Bakry et al, 2011',
                   expression(bold('    Tchounwou et al, 1991b')),
                   expression(bold('    Halstead et al, 2015')), 
                   expression(bold('    Tchounwou et al, 1992')),
                   expression(bold('    Tchounwou et al, 1991a')),
                   expression(bold('  Combined'))),
               
              list('Species', ' ', #atrazine
                   expression(italic('Biomphalaria alexandrina')),expression(italic('Biomphalaria alexandrina')),expression(italic('Physella spp.')),
                   expression(italic('Echinistoma trivolvis')), expression(italic('Echinistoma trivolvis')),expression(italic('Echinistoma trivolvis')),
                  ' ', ' ', ' ',  #Chlorpyrifos
                  expression(italic('Bulinus truncatus')), expression(italic('Biomphalaria alexandrina')), expression(italic('Bulinus truncatus')),
                  expression(italic('Biomphalaria alexandrina')), expression(italic('Procambarus clarkii')),expression(italic('Macrobrachium rosenbergii')),
                  expression(italic('Schistosoma haematobium')), expression(italic('Schistosoma haematobium')),expression(italic('Macrobrachium rosenbergii')),
                 ' ', ' ', ' ',   #Glyphosate
                 expression(italic('Biomphalaria alexandrina')),expression(italic('Biomphalaria alexandrina')),expression(italic('Biomphalaria alexandrina')),
                 expression(italic('Biomphalaria alexandrina')),expression(italic('Biomphalaria alexandrina')),
                 expression(italic('Schistosoma mansoni')),expression(italic('Schistosoma mansoni')),
                 ' ', ' ', ' ',   #Malathion
                 expression(italic('Helisoma duryi')),expression(italic('Bulinus havenensis')),expression(italic('Helisoma duryi')),
                 expression(italic('Bulinus havenensis')),expression(italic('Procambarus clarkii')),expression(italic('Schistosoma mansoni')),
                 expression(italic('Schistosoma mansoni')),' '),
               
              list('Parameter', ' ', #Atrazine
                   expression(paste(mu[N])), expression(paste(mu[N])),expression(paste(Phi[N])),
                   expression(paste(pi[C])), expression(paste(pi[C])),expression(paste(pi[C])),
                   ' ', ' ', ' ',    #Chlorpyrifos
                   expression(paste('f'[N])),expression(paste('f'[N])),expression(paste(mu[N])), 
                   expression(paste(mu[N])),expression(paste(mu[P])),expression(paste(mu[P])),
                   expression(paste(pi[C])),expression(paste(pi[M])),expression(paste(psi[P])),
                   ' ', ' ', ' ',     #Glyphosate
                   expression(paste('f'[N])),expression(paste('f'[N])),expression(paste(mu[N])),
                    expression(paste(mu[N])),expression(paste(mu[N])),expression(paste(pi[C])),expression(paste(pi[M])), 
                   ' ',' ', ' ',      #Malathion
                   expression(paste('f'[N])),expression(paste('f'[N])),expression(paste(mu[N])),expression(paste(mu[N])),
                   expression(paste(mu[P])),expression(paste(pi[C])),expression(paste(pi[M])),' ')
               
)

#Forest plots #########
#With raw change in r0
windows(width = 24, height = 15)
forestplot(tabtext, 
           vals.raw,
           new_page = TRUE,
           is.summary=c(TRUE,TRUE,
                        rep(FALSE,nrow(atr.eec.df)+1),
                        TRUE,
                        rep(FALSE,nrow(ch.eec.df)+1),
                        TRUE,
                        rep(FALSE,nrow(gly.eec.df)+1),
                        TRUE,
                        rep(FALSE,nrow(mal.eec.df)+1)),
           hrzl_lines = list('2'= gpar(lwd = 1, col="black", columns = 1:3),
                             '3'= gpar(col = 'grey'),
                             '12'= gpar(col = 'grey'),
                             '24'= gpar(col = 'grey'),
                             '34'= gpar(col = 'grey')),
           txt_gp = fpTxtGp(xlab = gpar(cex = 1.2),
                            ticks = gpar(cex = 1.1)),
           vertices = TRUE,
           boxsize = 0.25,
           col = fpColors(box = c('black', 'blue'), lines = 'black'),
           clip=c(-Inf,Inf),
           xlab = expression(paste(Delta, R['0'])))

#With relative change in r0
windows(width = 24, height = 15)
forestplot(tabtext, 
           vals.rel,
           new_page = TRUE,
           is.summary=c(TRUE,TRUE,
                        rep(FALSE,nrow(atr.eec.df)+1),
                        TRUE,
                        rep(FALSE,nrow(ch.eec.df)+1),
                        TRUE,
                        rep(FALSE,nrow(gly.eec.df)+1),
                        TRUE,
                        rep(FALSE,nrow(mal.eec.df)+1)),
           hrzl_lines = list('2'= gpar(lwd = 1, col="black", columns = 1:3),
                             '3'= gpar(col = 'grey'),
                             '12'= gpar(col = 'grey'),
                             '24'= gpar(col = 'grey'),
                             '34'= gpar(col = 'grey')),
           txt_gp = fpTxtGp(xlab = gpar(cex = 1.2),
                            ticks = gpar(cex = 1.1)),
           vertices = TRUE,
           boxsize = 0.25,
           col = fpColors(box = 'black', lines = 'black'),
           clip=c(0,2),
           xlab = expression(paste('Relative ', italic(R[0]))))

#Deprecated #########

tabtext2 = cbind(c('Study', ' Atrazine','  Bakry et al, 2012','  Griggs & Belden, 2008',
                  '  Koprivnikar et al, 2007','  Omran & Salama, 2013','  Rohr et al, 2008',
                  '  Rohr et al, 2012','  Combined',' Chlorpyrifos','  Halstead et al, 2015',
                  '  Hasheesh & Mohamed, 2011','  Hasheesh & Mohamed, 2011',
                  '  Hasheesh & Mohamed, 2011','  Hasheesh & Mohamed, 2011','  Ibrahim et al, 1992', 
                  '  Ibrahim et al, 1992','  Satapornvanit et al, 2009','  Satapornvanit et al, 2009','  Combined',
                  ' Glyphosate','  Abdel-Ghaffar et al, 2016','  Abdel-Ghaffar et al, 2016',
                  '  Abdel-Ghaffar et al, 2016','  Abdel-Ghaffar et al, 2016', '  Bakry et al, 2012',
                  '  Bakry et al, 2012','  Omran & Salama 2013','  Combined',' Malathion','  Bakry et al, 2011',
                  '  Bakry et al, 2011','  Halstead et al, 2015', '  Tchounwou et al, 1991a',
                  '  Tchounwou et al, 1991b','  Tchounwou et al, 1991b','  Tchounwou et al, 1992',
                  '  Combined'),
                
                c('Species', ' ','Biomphalaria alexandrina','Echinistoma trivolvis',
                  'Echinistoma trivolvis', 'Biomphalaria alexandrina', 'Echinistoma trivolvis',
                  'Physella spp.',' ', ' ',  'Procambarus clarkii', 'Schistosoma haematobium',
                  'Schistosoma haematobium', 'Bulinus truncatus', 'Bulinus truncatus',
                  'Biomphalaria alexandrina', 'Biomphalaria alexandrina',
                  'Macrobrachium rosenbergii', 'Macrobrachium rosenbergii',' ',' ',
                  'Schistosoma mansoni', 'Schistosoma mansoni', 'Biomphalaria alexandrina',
                  'Biomphalaria alexandrina', 'Biomphalaria alexandrina', 'Biomphalaria alexandrina',
                  'Biomphalaria alexandrina',' ', ' ', 'Helisoma duryi', 'Helisoma duryi',
                  'Procambarus clarkii', 'Schistosoma mansoni', 'Bulinus havenensis',
                  'Bulinus havenensis', 'Schistosoma mansoni',' '),
                
                c('Parameter', ' ', as.character(atr.eec.df$Parameter),  ' ', as.character(ch.eec.df$Parameter), ' ', 
                  as.character(gly.eec.df$Parameter),' ', as.character(mal.eec.df$Parameter))
                
)



tabtext = list(c('Study', 'Species', 'Parameter'),
               
               c(' Atrazine', ' ', ' '),
               
               c('  Bakry et al, 2012',  expression(italic('Biomphalaria alexandrina')), expression(paste(mu[N]))),
               c('  Griggs & Belden, 2008',  expression(italic('Echinistoma trivolvis')), expression(paste(pi[C]))),
               c('  Koprivnikar et al, 2007',  expression(italic('Echinistoma trivolvis')), expression(paste(pi[C]))),
               c(expression(bold('  Omran & Salama, 2013')),  expression(italic('Biomphalaria alexandrina')), expression(paste(mu[N]))),
               c(expression(bold('  Rohr et al, 2008')),  expression(italic('Echinistoma trivolvis')), expression(paste(pi[C]))),
               c(expression(bold('  Rohr et al, 2012')),  expression(italic('Physella spp.')), expression(paste(Phi[N]))),
               c(expression(bold('  Combined')),  ' ', ' '),
               
               c(' Chlorpyrifos', ' ', ' '),
               c(expression(bold('  Halstead et al, 2015')),  expression(italic('Procambarus clarkii')), expression(paste(mu[P]))),
               c(expression(bold('  Hasheesh & Mohamed, 2011')),  expression(italic('Schistosoma haematobium')), expression(paste(pi[C]))),
               c(expression(bold('  Hasheesh & Mohamed, 2011')),  expression(italic('Schistosoma haematobium')), expression(paste(pi[M]))),
               c(expression(bold('  Hasheesh & Mohamed, 2011')),  expression(italic('Bulinus truncatus')), expression(paste(mu[N]))),
               c('  Hasheesh & Mohamed, 2011',  expression(italic('Bulinus truncatus')), expression(paste('f'[N]))),
               c('  Ibrahim et al, 1992',  expression(italic('Biomphalaria alexandrina')), expression(paste(mu[N]))),
               c(expression(bold('  Ibrahim et al, 1992')),  expression(italic('Biomphalaria alexandrina')), expression(paste('f'[N]))),
               c('  Satapornvanit et al, 2009',  expression(italic('Macrobrachium rosenbergii')), expression(paste(mu[P]))),
               c('  Satapornvanit et al, 2009',  expression(italic('Macrobrachium rosenbergii')), expression(paste(psi[P]))),
               c(expression(bold('  Combined')), ' ', ' '),
               
               c(' Glyphosate', ' ', ' '),
               c(expression(bold('  Abdel-Ghaffar et al, 2016')),  expression(italic('Schistosoma mansoni')), expression(paste(pi[M]))),
               c(expression(bold('  Abdel-Ghaffar et al, 2016')),  expression(italic('Schistosoma mansoni')), expression(paste(pi[C]))),
               c(expression(bold('  Abdel-Ghaffar et al, 2016')),  expression(italic('Biomphalaria alexandrina')), expression(paste(mu[N]))),
               c(expression(bold('  Abdel-Ghaffar et al, 2016')),  expression(italic('Biomphalaria alexandrina')), expression(paste('f'[N]))),
               c('  Bakry et al, 2012',  expression(italic('Biomphalaria alexandrina')), expression(paste(mu[N]))),
               c('  Bakry et al, 2012',  expression(italic('Biomphalaria alexandrina')), expression(paste('f'[N]))),
               c('  Omran & Salama 2013',  expression(italic('Biomphalaria alexandrina')), expression(paste(mu[N]))),
               c(expression(bold('  Combined')), ' ', ' '),
               
               c(' Malathion', ' ', ' '),
               c('  Bakry et al, 2011',  expression(italic('Helisoma duryi')), expression(paste('f'[N]))),
               c('  Bakry et al, 2011',  expression(italic('Helisoma duryi')), expression(paste(mu[N]))),
               c(expression(bold('  Halstead et al, 2015')),  expression(italic('Procambarus clarkii')), expression(paste(mu[P]))),
               c(expression(bold('  Tchounwou et al, 1991a')),  expression(italic('Schistosoma mansoni')), expression(paste(pi[M]))),
               c(expression(bold('  Tchounwou et al, 1991b')),  expression(italic('Bulinus havenensis')), expression(paste(mu[N]))),
               c(expression(bold('  Tchounwou et al, 1991b')),  expression(italic('Bulinus havenensis')), expression(paste('f'[N]))),
               c(expression(bold('  Tchounwou et al, 1992')),  expression(italic('Schistosoma mansoni')), expression(paste(pi[C]))),
               c(expression(bold('  Combined')), ' ', ' ')
               )
